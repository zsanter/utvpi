#include "mine.h"

int main(int argc, char * argv[]){
  clock_t beginning = clock();
  if( argc != 2 && argc != 3 ){
    fprintf( stderr, "Proper use is %s [input file] {output file}.\nIf no output file is specified, output is to stdout.\n", argv[0] );
    exit(1);
  }
  FILE * input = fopen( argv[1], "r");
  if( input == NULL ){
    fprintf( stderr, "File \"%s\" could not be found.\n", argv[1] );
    exit(1);
  }
  FILE * output;
  if( argc == 2 ){
    output = stdout;
  }
  else{
    output = fopen( argv[2], "w");
    if( output == NULL ){
      fprintf( stderr, "File \"%s\" could not be created.\n", argv[2] );
      exit(1);
    }
  }
  System system;
  bool parseSuccessful = parseFile(input, &system, initializeSystem, addConstraint);
  fclose(input);
  if( !parseSuccessful ){
    fprintf( stderr, "Parsing of \"%s\" failed.\n", argv[1] );
    exit(1);
  }
  time_t afterSystemCreation = clock();
  Edge * negativeCycle = bellmanFord(&system);
  if( negativeCycle != NULL ){
    fputs("The following negative cost cycle was detected:\n", output);
    Edge * edge = negativeCycle;
    while( edge->backtrackSeen == true ){
      fputEdge( edge, output );
      edge->backtrackSeen = false;
      edge = edge->tail->L;
    }
  }
  else{
    fputs("Linear solution:\n", output);
    for(int i = 0; i < system.n; i++){
      half_int solution = intDivBy2ToHalfInt( system.graph[POSITIVE][i].D - system.graph[NEGATIVE][i].D );
      system.graph[POSITIVE][i].x = solution;
      system.graph[NEGATIVE][i].x = solution;
      fprintf(output, "x%i = %.1f\n", i + 1, halfIntToDouble( solution ) );
    }
  }
  fclose(output);
  time_t afterSolutionOutput = clock();
  cleanup(&system);
  time_t afterCleanup = clock();
  printf("%f,", ((double)(afterSystemCreation - beginning))/CLOCKS_PER_SEC );
  printf("%f,,", ((double)(afterSolutionOutput - afterSystemCreation))/CLOCKS_PER_SEC );
  printf("%f,", ((double)(afterCleanup - afterSolutionOutput))/CLOCKS_PER_SEC );
  printf("%f,", ((double)(afterCleanup - beginning))/CLOCKS_PER_SEC );
  return 0;
}

void fputEdge(Edge * edge, FILE * output){
  if( edge->tail->index == edge->head->index ){
    char sign;
    switch( edge->tail->sign ){
    case POSITIVE:
      sign = '-';
      break;
    case NEGATIVE:
      sign = '+';
    }
    fprintf(output, "%cx%d <= %d\n", sign, edge->tail->index, (edge->weight)/2);
  }
  else{
    char sign[2];
    if( edge->tail->sign == edge->head->sign ){
      switch( edge->tail->sign ){
      case POSITIVE:
        sign[0] = '-';
        sign[1] = '+';
        break;
      case NEGATIVE:
        sign[0] = '+';
        sign[1] = '-';
      }
    }
    else{
      switch( edge->tail->sign ){
      case POSITIVE:
        sign[0] = '-';
        sign[1] = '-';
        break;
      case NEGATIVE:
        sign[0] = '+';
        sign[1] = '+';
      }
    }
    fprintf(output, "%cx%d %cx%d <= %d\n", sign[0], edge->tail->index, sign[1], edge->head->index, edge->weight);
  }
}

void initializeSystem(void * object, int n, Parser * parser){
  System * system = (System *) object;
  system->n = n;
  system->graph[POSITIVE] = (Vertex *) malloc( sizeof(Vertex) * n );
  system->graph[NEGATIVE] = (Vertex *) malloc( sizeof(Vertex) * n );
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < system->n; j++){
      system->graph[i][j].index = j + 1;
      system->graph[i][j].sign = i;
      system->graph[i][j].L = NULL;
      system->graph[i][j].D = 0;
      system->graph[i][j].first = NULL;
      system->graph[i][j].x = 0;
    }
  }
}

void addConstraint(void * object, Constraint * constraint, Parser * parser){
  if( constraint->index[0] == constraint->index[1] ){
    parseError(parser, "Both indices in a constraint can not be the same." );
  }
  else{
    System * system = (System *) object;
    if( constraint->sign[1] == NONE ){
      Edge * newEdge = (Edge *) malloc( sizeof(Edge) );
      newEdge->weight = 2 * constraint->weight;
      newEdge->backtrackSeen = false;
      VertexSign tailSign, headSign;
      if( constraint->sign[0] == PLUS ){
        tailSign = NEGATIVE;
        headSign = POSITIVE;
      }
      else {
        tailSign = POSITIVE;
        headSign = NEGATIVE;
      }
      newEdge->tail = &system->graph[tailSign][ constraint->index[0] - 1 ];
      newEdge->head = &system->graph[headSign][ constraint->index[0] - 1 ];
      newEdge->next = newEdge->tail->first;
      newEdge->tail->first = newEdge;
    }
    else{
      Edge * newEdges[2];
      newEdges[0] = (Edge *) malloc( sizeof(Edge) );
      newEdges[0]->weight = constraint->weight;
      newEdges[0]->backtrackSeen = false;
      newEdges[1] = (Edge *) malloc( sizeof(Edge) );
      newEdges[1]->weight = constraint->weight;
      newEdges[1]->backtrackSeen = false;
      if( constraint->sign[0] != constraint->sign[1] ){
        int positiveIndex, negativeIndex;
        if( constraint->sign[0] == PLUS ){
          positiveIndex = 0;
          negativeIndex = 1;
        }
        else{
          positiveIndex = 1;
          negativeIndex = 0;
        }
        newEdges[0]->tail = &system->graph[POSITIVE][ constraint->index[ negativeIndex ] - 1 ];
        newEdges[0]->head = &system->graph[POSITIVE][ constraint->index[ positiveIndex ] - 1 ];
        newEdges[1]->tail = &system->graph[NEGATIVE][ constraint->index[ positiveIndex ] - 1 ];
        newEdges[1]->head = &system->graph[NEGATIVE][ constraint->index[ negativeIndex ] - 1 ];
      }
      else if( constraint->sign[0] == PLUS ){
        newEdges[0]->tail = &system->graph[NEGATIVE][ constraint->index[1] - 1 ];
        newEdges[0]->head = &system->graph[POSITIVE][ constraint->index[0] - 1 ];
        newEdges[1]->tail = &system->graph[NEGATIVE][ constraint->index[0] - 1 ];
        newEdges[1]->head = &system->graph[POSITIVE][ constraint->index[1] - 1 ];
      }
      else {
        newEdges[0]->tail = &system->graph[POSITIVE][ constraint->index[1] - 1 ];
        newEdges[0]->head = &system->graph[NEGATIVE][ constraint->index[0] - 1 ];
        newEdges[1]->tail = &system->graph[POSITIVE][ constraint->index[0] - 1 ];
        newEdges[1]->head = &system->graph[NEGATIVE][ constraint->index[1] - 1 ];
      }
      newEdges[0]->next = newEdges[0]->tail->first;
      newEdges[0]->tail->first = newEdges[0];
      newEdges[1]->next = newEdges[1]->tail->first;
      newEdges[1]->tail->first = newEdges[1];
    }
  }
}

Edge * bellmanFord(System * system){
  for(int i = 1; i <= (2 * system->n - 1); i++){
    for(VertexSign j = POSITIVE; j <= NEGATIVE; j++){
      for(int k = 0; k < system->n; k++){
        Edge * edge = system->graph[j][k].first;
        while(edge != NULL){
          relax(edge);
          edge = edge->next;
        }
      }
    }
  }
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < system->n; j++){
      Edge * edge = system->graph[i][j].first;
      while(edge != NULL){
        if( edge->head->D > edge->tail->D + edge->weight ){
          return backtrack(edge);
        }
        edge = edge->next;
      }
    }
  }
  return NULL;
}

void relax(Edge * edge){
  if( edge->head->D > edge->tail->D + edge->weight ){
    edge->head->D = edge->tail->D + edge->weight;
    edge->head->L = edge;
  }
}

Edge * backtrack(Edge * edge){
  while( edge->backtrackSeen == false ){
    edge->backtrackSeen = true;
    edge = edge->tail->L;
  }
  return edge;
}

void cleanup(System * system){
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < system->n; j++){
      Edge * edge = system->graph[i][j].first;
      while(edge != NULL){
        Edge * oldEdge = edge;
        edge = edge->next;
        free( oldEdge );
      }
    }
    free( system->graph[i] );
  }
}