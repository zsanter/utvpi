#include "subWojInt.h"
#include <time.h>

int main(int argc, char * argv[]){
  clock_t start = clock();
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
  finishSystemCreation(&system);
  clock_t beforeLinear = clock();
  int f;
  bool linearlyFeasible = relaxNetwork(&system);
  clock_t beforeIntegral = clock();
  if( !linearlyFeasible ){
    f = 0;
    fputs("The following negative gray cycle was detected:\n", output);
    EdgeRefListNode * R = system.infeasibilityProof->first;
    while( R != NULL ){
      fputEdge(R->edge, output);
      R = R->next;
    }
  }
  else{
    fputs("Linear solution:\n", output);
    for(int i = 1; i < system.vertexCount; i++){
      system.graph[i].a = intDivBy2ToHalfInt( system.graph[i].D[WHITE] - system.graph[i].D[BLACK] );
      fprintf(output, "x%i = %.1f\n", i, halfIntToDouble( system.graph[i].a ) );
    }
    bool integrallyFeasible = produceIntegerSolution(&system);
    if( !integrallyFeasible ){
      f = 1;
      fputs("\nSystem is not integrally feasible.\nProof:\n", output);
      EdgeRefListNode * O = system.infeasibilityProof->first;
      while( O != NULL ){
        fputEdge(O->edge, output);
        O = O->next;
      }
    }
    else {
      f = 2;
      fputs("\nIntegral solution:\n", output);
      for(int i = 1; i < system.vertexCount; i++){
        fprintf(output, "x%i = %i\n", i, system.graph[i].Z[FINAL] );
      }
    }
  }
  clock_t beforeCleanup = clock();
  freeSystem(&system);
  fclose(output);
  clock_t end = clock();
  printf("%i,", f);
  printf("%f,", ((double)(beforeLinear - start))/CLOCKS_PER_SEC);
  printf("%f,", ((double)(beforeIntegral - beforeLinear))/CLOCKS_PER_SEC);
  printf("%f,", ((double)(beforeCleanup - beforeIntegral))/CLOCKS_PER_SEC);
  printf("%f,", ((double)(end - beforeCleanup))/CLOCKS_PER_SEC);
  printf("%f,", ((double)(end - start))/CLOCKS_PER_SEC);
  return 0;
}

EdgeType reverseEdgeType(EdgeType input){ 
  EdgeType output;
  switch(input){
  case GRAY_FORWARD:
    output = GRAY_REVERSE;
    break;
  case GRAY_REVERSE:
    output = GRAY_FORWARD;
    break;
  default:
    output = input;
  }
  return output;
}

void fputEdge(Edge * edge, FILE * output){
  char sign[2];
  switch(edge->type){
  case WHITE:
    sign[0] = '+';
    sign[1] = '+';
    break;
  case BLACK:
    sign[0] = '-';
    sign[1] = '-';
    break;
  case GRAY_FORWARD:
    sign[0] = '-';
    sign[1] = '+';
    break;
  case GRAY_REVERSE:
    sign[0] = '+';
    sign[1] = '-';
    break;
  }
  if( edge->tail->index != 0 ){
    fprintf(output, "%cx%i ", sign[0], edge->tail->index);
  }
  if( edge->head->index != 0 ){
    fprintf(output, "%cx%i ", sign[1], edge->head->index);
  }
  fprintf(output, "<= %i\n", edge->weight);
}

void initializeSystem(void * object, int n, Parser * parser){
  System * system = (System *) object;
  system->vertexCount = n + 1;
  system->graph = (Vertex *) malloc( sizeof(Vertex) * system->vertexCount );
  system->n = n;
  system->C = 0;
  system->allEdgeFirst = NULL;
  system->infeasibilityProof = NULL;
  system->T[DOWN] = NULL;
  system->T[UP] = NULL;
  system->graph[0].index = 0;
  system->graph[0].E[NEG_ONE] = NULL;
  system->graph[0].E[POS_ONE] = NULL;
  system->graph[0].a = 0;
  system->graph[0].Z[FINAL] = 0;
  system->graph[0].Z[TEMP] = 0;
  system->graph[0].integerTreeVertex[DOWN] = NULL;
  system->graph[0].integerTreeVertex[UP] = NULL;
  for(EdgeType i = WHITE; i <= GRAY_REVERSE; i++){
    system->graph[0].L[i] = NULL;
    system->graph[0].D[i] = 0;
    system->graph[0].first[i] = NULL;
    system->graph[0].edgeCount[i] = 0;
  }
  for(int i = 1; i < system->vertexCount; i++){
    system->graph[i].index = i;
    system->graph[i].E[NEG_ONE] = NULL;
    system->graph[i].E[POS_ONE] = NULL;
    system->graph[i].a = 0;
    system->graph[i].Z[FINAL] = INT_MAX;
    system->graph[i].Z[TEMP] = INT_MAX;
    system->graph[i].integerTreeVertex[DOWN] = NULL;
    system->graph[i].integerTreeVertex[UP] = NULL;
    for(EdgeType j = WHITE; j <= GRAY_REVERSE; j++){
      system->graph[i].L[j] = NULL;
      system->graph[i].D[j] = INT_MAX;
      system->graph[i].first[j] = NULL;
      system->graph[i].edgeCount[j] = 0;
    }
  }
}

void addConstraint(void * object, Constraint * constraint, Parser * parser){
  if( constraint->index[0] == constraint->index[1] ){
    parseError(parser, "Both indices in a constraint can not be the same." );
  }
  else{
    System * system = (System *) object;
    if( abs( constraint->weight ) > system->C ){
      system->C = abs( constraint->weight );
    }
    if( constraint->sign[1] == NONE ){
      constraint->index[1] = 0;
      constraint->sign[1] = PLUS;
      addEdge( system, constraint );
      constraint->sign[1] = MINUS;
      addEdge( system, constraint );
      //Calls to addEdge place the new edges at first[edgeType] for both vertices
      if( constraint->sign[0] == PLUS ){
        if( constraint->weight < system->graph[ constraint->index[0] ].D[WHITE] ){
          system->graph[ constraint->index[0] ].D[WHITE] = constraint->weight;
          system->graph[ constraint->index[0] ].L[WHITE] = system->graph[0].first[WHITE];
          system->graph[ constraint->index[0] ].D[GRAY_FORWARD] = constraint->weight;
          system->graph[ constraint->index[0] ].L[GRAY_FORWARD] = system->graph[0].first[GRAY_FORWARD];
        }
      }
      else{
        if( constraint->weight < system->graph[ constraint->index[0] ].D[BLACK] ){
          system->graph[ constraint->index[0] ].D[BLACK] = constraint->weight;
          system->graph[ constraint->index[0] ].L[BLACK] = system->graph[0].first[BLACK];
          system->graph[ constraint->index[0] ].D[GRAY_REVERSE] = constraint->weight;
          system->graph[ constraint->index[0] ].L[GRAY_REVERSE] = system->graph[0].first[GRAY_REVERSE]; //.first[GRAY_FORWARD]; ? No.
        }
      }
    }
    else{
      addEdge( system, constraint );
    }
  }
}

void addEdge(System * system, Constraint * constraint){
  Edge * newEdges[2];
  newEdges[0] = (Edge *) malloc( sizeof(Edge) );
  newEdges[1] = (Edge *) malloc( sizeof(Edge) );
  newEdges[0]->weight = constraint->weight;
  newEdges[0]->reverse = newEdges[1];
  newEdges[0]->next = NULL;
  newEdges[1]->weight = constraint->weight;
  newEdges[1]->reverse = newEdges[0];
  newEdges[1]->next = NULL;
  if( constraint->sign[0] == constraint->sign[1] ){
    EdgeType edgeType;
    if( constraint->sign[0] == PLUS ){
      edgeType = WHITE;
    }
    else{
      edgeType = BLACK;
    }
    newEdges[0]->type = edgeType;
    newEdges[0]->tail = &system->graph[ constraint->index[0] ];
    newEdges[0]->head = &system->graph[ constraint->index[1] ];
    newEdges[1]->type = edgeType;
    newEdges[1]->tail = &system->graph[ constraint->index[1] ];
    newEdges[1]->head = &system->graph[ constraint->index[0] ];
  }
  else{
    int negativeIndex, positiveIndex;
    if( constraint->sign[0] == PLUS ){
      positiveIndex = 0;
      negativeIndex = 1;
    }
    else{
      positiveIndex = 1;
      negativeIndex = 0;
    }
    newEdges[0]->type = GRAY_FORWARD;
    newEdges[0]->tail = &system->graph[ constraint->index[ negativeIndex ] ];
    newEdges[0]->head = &system->graph[ constraint->index[ positiveIndex ] ];
    newEdges[1]->type = GRAY_REVERSE;
    newEdges[1]->tail = &system->graph[ constraint->index[ positiveIndex ] ];
    newEdges[1]->head = &system->graph[ constraint->index[ negativeIndex ] ];
  }
  newEdges[0]->next = newEdges[0]->tail->first[ newEdges[0]->type ];
  newEdges[0]->tail->first[ newEdges[0]->type ] = newEdges[0];
  newEdges[0]->tail->edgeCount[ newEdges[0]->type ]++;
  newEdges[0]->allNext = system->allEdgeFirst;
  newEdges[0]->allPrev = NULL;
  newEdges[0]->inAllEdgeList = true;
  if( system->allEdgeFirst != NULL ){
    system->allEdgeFirst->allPrev = newEdges[0];
  }
  system->allEdgeFirst = newEdges[0];
  newEdges[1]->next = newEdges[1]->tail->first[ newEdges[1]->type ];
  newEdges[1]->tail->first[ newEdges[1]->type ] = newEdges[1];
  newEdges[1]->tail->edgeCount[ newEdges[1]->type ]++;
  newEdges[1]->allNext = NULL;
  newEdges[1]->allPrev = NULL;
  newEdges[1]->inAllEdgeList = false;
}

void finishSystemCreation(System * system){
  int sourceEdgeWeights = (2 * system->n + 1) * system->C;
  for(int i = 0; i < system->vertexCount; i++){
    for(EdgeType j = WHITE; j <= GRAY_REVERSE; j++){
      if( system->graph[i].D[j] == INT_MAX ){
        system->graph[i].D[j] = sourceEdgeWeights;
      }
      if( system->graph[i].edgeCount[j] > 1 ){
        int edgeSortArrayLength = system->graph[i].edgeCount[j];
        Edge * edgeSortArray[ edgeSortArrayLength ];
        Edge * edge = system->graph[i].first[j];
        for(int k = 0; k < edgeSortArrayLength; k++){
          edgeSortArray[k] = edge;
          edge = edge->next;
        }
        qsort(edgeSortArray, edgeSortArrayLength, sizeof(Edge *), edgeCompare);
        system->graph[i].first[j] = edgeSortArray[0];
        Edge * beforePrior = NULL;
        Edge * prior = edgeSortArray[0];
        for(int k = 1; k < edgeSortArrayLength; k++){
          if( prior->head == edgeSortArray[k]->head ){
            if( prior->weight <= edgeSortArray[k]->weight ){
              removeFromAllEdgeList( system, edgeSortArray[k] );
              free( edgeSortArray[k] );
            }
            else{
              if( system->graph[i].first[j] == prior ){
                system->graph[i].first[j] = edgeSortArray[k];
              }
              else {
                beforePrior->next = edgeSortArray[k];
              }
              removeFromAllEdgeList( system, prior );
              free( prior );
              prior = edgeSortArray[k];
            }
            system->graph[i].edgeCount[j]--;
          }
          else{
            prior->next = edgeSortArray[k];
            beforePrior = prior;
            prior = edgeSortArray[k];
          }
        }
        prior->next = NULL;
      }
    }
  }
}

int edgeCompare(const void * edge1, const void * edge2){
  return (*(Edge **)edge1)->head->index - (*(Edge **)edge2)->head->index;
}

void removeFromAllEdgeList(System * system, Edge * edge){
  if( edge->inAllEdgeList == true ){
    if( edge->allNext != NULL ){
      edge->allNext->allPrev = edge->allPrev;
    }
    if( edge->allPrev != NULL ){
      edge->allPrev->allNext = edge->allNext;
    }
    else{
      system->allEdgeFirst = edge->allNext;
    }
    edge->allNext = NULL;
    edge->allPrev = NULL;
    edge->inAllEdgeList = false;
  }
}

bool relaxNetwork(System * system){
  //Lines 3-6 of algorithm implemented in finishSystemCreation().
  for(int r = 1; r <= 2 * system->n; r++){
    Edge * e = system->allEdgeFirst;
    while( e != NULL ){
      relaxEdge(e);
      e = e->allNext;
    }
  }
  for(int i = 0; i < system->vertexCount; i++){
    Edge * e = system->graph[i].first[WHITE];
    while( e != NULL ){
      if( e->head->D[GRAY_REVERSE] + e->weight < e->tail->D[WHITE] ){
        return backtrack(system, e->head, GRAY_REVERSE, e->reverse);
      }
      if( e->head->D[BLACK] + e->weight < e->tail->D[GRAY_FORWARD] ){
        return backtrack(system, e->head, BLACK, e->reverse);
      }
      if( e->tail->D[GRAY_REVERSE] + e->weight < e->head->D[WHITE] ){
        return backtrack(system, e->tail, GRAY_REVERSE, e);
      }
      if( e->tail->D[BLACK] + e->weight < e->head->D[GRAY_FORWARD] ){
        return backtrack(system, e->tail, BLACK, e);
      }
      e = e->next;
    }
  }
  for(int i = 0; i < system->vertexCount; i++){
    Edge * e = system->graph[i].first[BLACK];
    while( e != NULL ){
      if( e->head->D[GRAY_FORWARD] + e->weight < e->tail->D[BLACK] ){
        return backtrack(system, e->head, GRAY_FORWARD, e->reverse);
      }
      if( e->head->D[WHITE] + e->weight < e->tail->D[GRAY_REVERSE] ){
        return backtrack(system, e->head, WHITE, e->reverse);
      }
      if( e->tail->D[GRAY_FORWARD] + e->weight < e->head->D[BLACK] ){
        return backtrack(system, e->tail, GRAY_FORWARD, e);
      }
      if( e->tail->D[WHITE] + e->weight < e->head->D[GRAY_REVERSE] ){
        return backtrack(system, e->tail, WHITE, e);
      }
      e = e->next;
    }
  }
  for(int i = 0; i < system->vertexCount; i++){
    Edge * e = system->graph[i].first[GRAY_FORWARD];
    while( e != NULL ){
      if( e->head->D[GRAY_REVERSE] + e->weight < e->tail->D[GRAY_REVERSE] ){
        return backtrack(system, e->head, GRAY_REVERSE, e->reverse);
      }
      if( e->head->D[BLACK] + e->weight < e->tail->D[BLACK] ){
        return backtrack(system, e->head, BLACK, e->reverse);
      }
      if( e->tail->D[GRAY_FORWARD] + e->weight < e->head->D[GRAY_FORWARD] ){
        return backtrack(system, e->tail, GRAY_FORWARD, e);
      }
      if( e->tail->D[WHITE] + e->weight < e->head->D[WHITE] ){
        return backtrack(system, e->tail, WHITE, e);
      }
      e = e->next;
    }
  }
  for(int i = 0; i < system->vertexCount; i++){
    Edge * e = system->graph[i].first[GRAY_REVERSE];
    while( e != NULL ){
      if( e->head->D[GRAY_FORWARD] + e->weight < e->tail->D[GRAY_FORWARD] ){
        return backtrack(system, e->head, GRAY_FORWARD, e->reverse);
      }
      if( e->head->D[WHITE] + e->weight < e->tail->D[WHITE] ){
        return backtrack(system, e->head, WHITE, e->reverse);
      }
      if( e->tail->D[GRAY_REVERSE] + e->weight < e->head->D[GRAY_REVERSE] ){
        return backtrack(system, e->tail, GRAY_REVERSE, e);
      }
      if( e->tail->D[BLACK] + e->weight < e->head->D[BLACK] ){
        return backtrack(system, e->tail, BLACK, e);
      }
      e = e->next;
    }
  }
  return true;
}

void relaxEdge(Edge * e){
  switch(e->type){
  case WHITE:
    if( e->head->D[GRAY_REVERSE] + e->weight < e->tail->D[WHITE] ){
      e->tail->D[WHITE] = e->head->D[GRAY_REVERSE] + e->weight;
      e->tail->L[WHITE] = e->reverse;
    }
    if( e->head->D[BLACK] + e->weight < e->tail->D[GRAY_FORWARD] ){
      e->tail->D[GRAY_FORWARD] = e->head->D[BLACK] + e->weight;
      e->tail->L[GRAY_FORWARD] = e->reverse;
    }
    if( e->tail->D[GRAY_REVERSE] + e->weight < e->head->D[WHITE] ){
      e->head->D[WHITE] = e->tail->D[GRAY_REVERSE] + e->weight;
      e->head->L[WHITE] = e;
    }
    if( e->tail->D[BLACK] + e->weight < e->head->D[GRAY_FORWARD] ){
      e->head->D[GRAY_FORWARD] = e->tail->D[BLACK] + e->weight;
      e->head->L[GRAY_FORWARD] = e;
    }
    break;
  case BLACK:
    if( e->head->D[GRAY_FORWARD] + e->weight < e->tail->D[BLACK] ){
      e->tail->D[BLACK] = e->head->D[GRAY_FORWARD] + e->weight;
      e->tail->L[BLACK] = e->reverse;
    }
    if( e->head->D[WHITE] + e->weight < e->tail->D[GRAY_REVERSE] ){
      e->tail->D[GRAY_REVERSE] = e->head->D[WHITE] + e->weight;
      e->tail->L[GRAY_REVERSE] = e->reverse;
    }
    if( e->tail->D[GRAY_FORWARD] + e->weight < e->head->D[BLACK] ){
      e->head->D[BLACK] = e->tail->D[GRAY_FORWARD] + e->weight;
      e->head->L[BLACK] = e;
    }
    if( e->tail->D[WHITE] + e->weight < e->head->D[GRAY_REVERSE] ){
      e->head->D[GRAY_REVERSE] = e->tail->D[WHITE] + e->weight;
      e->head->L[GRAY_REVERSE] = e;
    }
    break;
  case GRAY_FORWARD:
    if( e->head->D[GRAY_REVERSE] + e->weight < e->tail->D[GRAY_REVERSE] ){
      e->tail->D[GRAY_REVERSE] = e->head->D[GRAY_REVERSE] + e->weight;
      e->tail->L[GRAY_REVERSE] = e->reverse;
    }
    if( e->head->D[BLACK] + e->weight < e->tail->D[BLACK] ){
      e->tail->D[BLACK] = e->head->D[BLACK] + e->weight;
      e->tail->L[BLACK] = e->reverse;
    }
    if( e->tail->D[GRAY_FORWARD] + e->weight < e->head->D[GRAY_FORWARD] ){
      e->head->D[GRAY_FORWARD] = e->tail->D[GRAY_FORWARD] + e->weight;
      e->head->L[GRAY_FORWARD] = e;
    }
    if( e->tail->D[WHITE] + e->weight < e->head->D[WHITE] ){
      e->head->D[WHITE] = e->tail->D[WHITE] + e->weight;
      e->head->L[WHITE] = e;
    }
    break;
  case GRAY_REVERSE:
    if( e->head->D[GRAY_FORWARD] + e->weight < e->tail->D[GRAY_FORWARD] ){
      e->tail->D[GRAY_FORWARD] = e->head->D[GRAY_FORWARD] + e->weight;
      e->tail->L[GRAY_FORWARD] = e->reverse;
    }
    if( e->head->D[WHITE] + e->weight < e->tail->D[WHITE] ){
      e->tail->D[WHITE] = e->head->D[WHITE] + e->weight;
      e->tail->L[WHITE] = e->reverse;
    }
    if( e->tail->D[GRAY_REVERSE] + e->weight < e->head->D[GRAY_REVERSE] ){
      e->head->D[GRAY_REVERSE] = e->tail->D[GRAY_REVERSE] + e->weight;
      e->head->L[GRAY_REVERSE] = e;
    }
    if( e->tail->D[BLACK] + e->weight < e->head->D[BLACK] ){
      e->head->D[BLACK] = e->tail->D[BLACK] + e->weight;
      e->head->L[BLACK] = e;
    }
  }
}

bool backtrack(System * system, Vertex * x_i, EdgeType t, Edge * e){
  Vertex * x_c = x_i;
  Edge * e_c = e;
  EdgeType t_c = t;
  BacktrackingIndex a_c;
  switch(t_c){
  case WHITE:
  case GRAY_FORWARD:
    a_c = POS_ONE;
    break;
  case BLACK:
  case GRAY_REVERSE:
    a_c = NEG_ONE;
    break;
  }
  while( x_c->E[a_c] == NULL ){
    x_c->E[a_c] = e_c;
    e_c = x_c->L[t_c];
    switch(t_c){
    case WHITE:
      switch(e_c->type){
      case WHITE:
        t_c = GRAY_REVERSE;
        a_c = NEG_ONE;
        break;
      case GRAY_FORWARD:
        t_c = WHITE;
        a_c = POS_ONE;
        break;
      }
      break;
    case BLACK:
      switch(e_c->type){
      case BLACK:
        t_c = GRAY_FORWARD;
        a_c = POS_ONE;
        break;
      case GRAY_REVERSE:
        t_c = BLACK;
        a_c = NEG_ONE;
        break;
      }
      break;
    case GRAY_FORWARD:
      switch(e_c->type){
      case GRAY_FORWARD:
        t_c = GRAY_FORWARD;
        a_c = POS_ONE;
        break;
      case WHITE:
        t_c = BLACK;
        a_c = NEG_ONE;
        break;
      }
      break;
    case GRAY_REVERSE:
      switch(e_c->type){
      case GRAY_REVERSE:
        t_c = GRAY_REVERSE;
        a_c = NEG_ONE;
        break;
      case BLACK:
        t_c = WHITE;
        a_c = POS_ONE;
        break;
      }
      break;
    }
    x_c = e_c->tail;
  }
  Vertex * x_f = x_c;
  BacktrackingIndex a_f = a_c;
  system->infeasibilityProof = generateEdgeRefList();
  addEdgeToEdgeRefListEnd( system->infeasibilityProof, e_c );
  switch(e_c->type){
  case WHITE:
  case GRAY_FORWARD:
    a_c = POS_ONE;
    break;
  case BLACK:
  case GRAY_REVERSE:
    a_c = NEG_ONE;
    break;
  }
  x_c = e_c->head;
  while(a_c != a_f || x_c != x_f){
    addEdgeToEdgeRefListEnd( system->infeasibilityProof, x_c->E[a_c] );
    Vertex * x_c_next = x_c->E[a_c]->head;
    switch( x_c->E[a_c]->type ){
    case WHITE:
    case GRAY_FORWARD:
      a_c = POS_ONE;
      break;
    case BLACK:
    case GRAY_REVERSE:
      a_c = NEG_ONE;
      break;
    }
    x_c = x_c_next;
  }
  return false;
}

bool produceIntegerSolution(System * system){

  bool integrallyFeasible;
  system->infeasibilityProof = generateEdgeRefList();
  system->T[DOWN] = generateIntegerTree(system, DOWN);

  for(int i = 1; i < system->vertexCount; i++){
    if( halfIntIsIntegral( system->graph[i].a ) ){
      system->graph[i].Z[FINAL] = halfIntToInt( system->graph[i].a );
    }
    else{
      integrallyFeasible = forcedRounding( system, &system->graph[i] );
      if( !integrallyFeasible ){
        return false;
      }
    }
  } 

  Vertex * vertex = pollIntegerTreeQueue( system->T[DOWN] );
  while( vertex != NULL ){
    checkDependencies(system, system->T[DOWN], vertex, FINAL);
    vertex = pollIntegerTreeQueue( system->T[DOWN] );
  }

  integrallyFeasible = checkAllConstraints(system, &system->graph[0], FINAL, DOWN);
  if( !integrallyFeasible ){
    return false;
  }

  freeEdgeRefList( system->infeasibilityProof );
  freeIntegerTree( system->T[DOWN] );

  for(int i = 0; i < system->vertexCount; i++){
    for(EdgeType j = WHITE; j <= GRAY_REVERSE; j++){
      Edge * prior = NULL;
      Edge * edge = system->graph[i].first[j];
      while( edge != NULL ){
        if( edge->tail->Z[FINAL] == INT_MAX && edge->head->Z[FINAL] == INT_MAX ){
          if( prior == NULL ){
            system->graph[i].first[j] = edge;
          }
          else {
            prior->next = edge;
          }
          prior = edge;
        }
        else {
          removeFromAllEdgeList(system, edge);
          free(edge);
        }
	edge = edge->next;
      }
    }
  }

  integrallyFeasible = optionalRoundings(system);
  
  return integrallyFeasible;
}

bool forcedRounding(System * system, Vertex * x_i){
  bool forcedDown = false;
  Edge * whiteEdge = x_i->first[WHITE];
  Edge * grayReverseEdge = x_i->first[GRAY_REVERSE];
  while( whiteEdge != NULL && grayReverseEdge != NULL && !forcedDown ){
    if( grayReverseEdge->head->index == whiteEdge->head->index ){
      if( (intToHalfInt( whiteEdge->weight ) == x_i->a + whiteEdge->head->a) 
          && (intToHalfInt( grayReverseEdge->weight ) == x_i->a - whiteEdge->head->a) ){
      
        x_i->Z[FINAL] = halfIntToInt( halfIntFloor( x_i->a ) );
        Edge * newEdge = generateAbsoluteConstraint(system, system->T[DOWN], x_i, x_i->Z[FINAL], WHITE);
        expandIntegerTree( system, system->T[DOWN], x_i, &system->graph[0], whiteEdge, grayReverseEdge, newEdge );
        forcedDown = true;
      }
      else {
        whiteEdge = whiteEdge->next;
        grayReverseEdge = grayReverseEdge->next;
      }
    }
    if( whiteEdge != NULL ){
      while( grayReverseEdge != NULL && grayReverseEdge->head->index < whiteEdge->head->index ){
        grayReverseEdge = grayReverseEdge->next;
      }
    }
    if( grayReverseEdge != NULL ){
      while( whiteEdge != NULL && whiteEdge->head->index < grayReverseEdge->head->index ){
        whiteEdge = whiteEdge->next;
      }
    }
  }
  bool forcedUp = false;
  Edge * blackEdge = x_i->first[BLACK];
  Edge * grayForwardEdge = x_i->first[GRAY_FORWARD];
  while( blackEdge != NULL && grayForwardEdge != NULL && !forcedUp ){
    if( grayForwardEdge->head->index == blackEdge->head->index ){ 
      if( (intToHalfInt( blackEdge->weight ) == -x_i->a - blackEdge->head->a) 
          && (intToHalfInt( grayForwardEdge->weight ) == -x_i->a + blackEdge->head->a) ){
      
        x_i->Z[FINAL] = halfIntToInt( halfIntCeil( x_i->a ) );
        Edge * newEdge = generateAbsoluteConstraint(system, system->T[DOWN], x_i, -x_i->Z[FINAL], BLACK);
        expandIntegerTree( system, system->T[DOWN], x_i, &system->graph[0], blackEdge, grayForwardEdge, newEdge );
        forcedUp = true;
      }
      else {
        blackEdge = blackEdge->next;
        grayForwardEdge = grayForwardEdge->next;
      }
    }
    if( blackEdge != NULL ){
      while( grayForwardEdge != NULL && grayForwardEdge->head->index < blackEdge->head->index ){
        grayForwardEdge = grayForwardEdge->next;
      }
    }
    if( grayForwardEdge != NULL ){
      while( blackEdge != NULL && blackEdge->head->index < grayForwardEdge->head->index ){
        blackEdge = blackEdge->next;
      }
    }
  }
  if( forcedDown && forcedUp ){
    copyTreeEdgesToList(system->infeasibilityProof, x_i->integerTreeVertex[DOWN]);
    return false;
  }
  return true;
}

bool optionalRoundings(System * system){  
  for(int i = 1; i < system->vertexCount; i++){
    if( system->graph[i].Z[FINAL] == INT_MAX ){
      system->T[DOWN] = generateIntegerTree(system, DOWN);
      system->T[UP] = generateIntegerTree(system, UP);
      system->infeasibilityProof = generateEdgeRefList();
      Edge * newEdge;
      Vertex * vertex;
      
      for(int j = 1; j < system->vertexCount; j++){
        system->graph[j].Z[TEMP] = INT_MAX;
      }
      
      system->graph[i].Z[TEMP] = halfIntToInt( halfIntFloor( system->graph[i].a ) );
      newEdge = generateAbsoluteConstraint(system, system->T[DOWN], &system->graph[i], system->graph[i].Z[TEMP], WHITE);
      expandIntegerTree(system, system->T[DOWN], &system->graph[i], &system->graph[0], newEdge, NULL, NULL);

      vertex = pollIntegerTreeQueue( system->T[DOWN] );
      while( vertex != NULL ){
	checkDependencies(system, system->T[DOWN], vertex, TEMP);
        vertex = pollIntegerTreeQueue( system->T[DOWN] );
      }

      bool floorFeasible = checkAllConstraints(system, &system->graph[i], TEMP, DOWN);
      
      if( floorFeasible ){
        for(int j = 1; j < system->vertexCount; j++){
          if( system->graph[j].Z[TEMP] != INT_MAX ){
            system->graph[j].Z[FINAL] = system->graph[j].Z[TEMP];
          }
        }  
      }
      else {

        for(int j = 1; j < system->vertexCount; j++){
          system->graph[j].Z[TEMP] = INT_MAX;
        }
      
        system->graph[i].Z[TEMP] = halfIntToInt( halfIntCeil( system->graph[i].a ) );
        newEdge = generateAbsoluteConstraint( system, system->T[UP], &system->graph[i], -system->graph[i].Z[TEMP], BLACK );
        expandIntegerTree( system, system->T[UP], &system->graph[i], &system->graph[0], newEdge, NULL, NULL);

        vertex = pollIntegerTreeQueue( system->T[UP] );
        while( vertex != NULL ){
          checkDependencies( system, system->T[UP], vertex, TEMP );
          vertex = pollIntegerTreeQueue( system->T[UP] );
        }

        bool ceilFeasible = checkAllConstraints(system, &system->graph[i], TEMP, UP );
        
        if( ceilFeasible ){
          for(int j = 1; j < system->vertexCount; j++){
            if( system->graph[j].Z[TEMP] != INT_MAX ){
              system->graph[j].Z[FINAL] = system->graph[j].Z[TEMP];
            }
          }
        }
        else {
          return false;
        }
      
      }
      
      freeIntegerTree(system->T[DOWN]);
      freeIntegerTree(system->T[UP]);
      freeEdgeRefList(system->infeasibilityProof);
    }
  }
  system->T[DOWN] = NULL;
  system->T[UP] = NULL;
  system->infeasibilityProof = NULL;
  return true;
}

//x_i - must never represent a Vertex where x_i->a is integral. The first for loop in produceIntegerSolution
//ensures this
void checkDependencies(System * system, IntegerTree * T, Vertex * x_i, IntegerType integerType){
  if( x_i->Z[integerType] == halfIntToInt( halfIntFloor( x_i->a ) ) ){
    Edge * grayForwardEdge = x_i->first[GRAY_FORWARD];
    while( grayForwardEdge != NULL ){
      if( intToHalfInt( grayForwardEdge->weight ) == -x_i->a + grayForwardEdge->head->a 
          && grayForwardEdge->head->Z[integerType] == INT_MAX ){
        grayForwardEdge->head->Z[integerType] = halfIntToInt( halfIntFloor( grayForwardEdge->head->a ) );
        expandIntegerTree(system, T, grayForwardEdge->head, x_i, grayForwardEdge, NULL, NULL);
      }
      grayForwardEdge = grayForwardEdge->next;
    }
    Edge * blackEdge = x_i->first[BLACK];
    while( blackEdge != NULL ){
      if( intToHalfInt( blackEdge->weight ) == -x_i->a - blackEdge->head->a 
          && blackEdge->head->Z[integerType] == INT_MAX ){
        blackEdge->head->Z[integerType] = halfIntToInt( halfIntCeil( blackEdge->head->a ) );
        expandIntegerTree(system, T, blackEdge->head, x_i, blackEdge, NULL, NULL);
      }
      blackEdge = blackEdge->next;
    }
  }
  else {
    Edge * whiteEdge = x_i->first[WHITE];
    while( whiteEdge != NULL ){
      if( intToHalfInt( whiteEdge->weight ) == x_i->a + whiteEdge->head->a 
          && whiteEdge->head->Z[integerType] == INT_MAX ){
        whiteEdge->head->Z[integerType] = halfIntToInt( halfIntFloor( whiteEdge->head->a ) );
        expandIntegerTree(system, T, whiteEdge->head, x_i, whiteEdge, NULL, NULL);
      }
      whiteEdge = whiteEdge->next;
    }
    Edge * grayReverseEdge = x_i->first[GRAY_REVERSE];
    while( grayReverseEdge != NULL ){
      if( intToHalfInt( grayReverseEdge->weight ) == x_i->a - grayReverseEdge->head->a 
          && grayReverseEdge->head->Z[integerType] == INT_MAX ){
        grayReverseEdge->head->Z[integerType] = halfIntToInt( halfIntCeil( grayReverseEdge->head->a ) );
        expandIntegerTree(system, T, grayReverseEdge->head, x_i, grayReverseEdge, NULL, NULL);
      }
      grayReverseEdge = grayReverseEdge->next;
    }
  }
}

bool checkAllConstraints(System * system, Vertex * toVertex, IntegerType integerType, IntegerTreeType integerTreeType){
  Edge * edge = system->allEdgeFirst;
  while( edge != NULL ){
    if( edge->tail->Z[integerType] != INT_MAX && edge->head->Z[integerType] != INT_MAX ){
      bool feasible = true;
      switch( edge->type ){
      case WHITE:
        if( edge->tail->Z[integerType] + edge->head->Z[integerType] > edge->weight ){
          feasible = false;
        }
        break;
      case BLACK:
        if( -edge->tail->Z[integerType] - edge->head->Z[integerType] > edge->weight ){
          feasible = false;
        }
        break;
      case GRAY_FORWARD:
        if( -edge->tail->Z[integerType] + edge->head->Z[integerType] > edge->weight ){
          feasible = false;
        }
        break;
      case GRAY_REVERSE:
        if( edge->tail->Z[integerType] - edge->head->Z[integerType] > edge->weight ){
          feasible = false;
        }
      }
      if( !feasible ){
        addEdgeToEdgeRefListBeginning(system->infeasibilityProof, edge);
	integerTreeBacktrack(system->infeasibilityProof, edge->tail->integerTreeVertex[integerTreeType], toVertex->integerTreeVertex[integerTreeType], false);
        integerTreeBacktrack(system->infeasibilityProof, edge->head->integerTreeVertex[integerTreeType], toVertex->integerTreeVertex[integerTreeType], true);
        return false;
      }
    }
    edge = edge->allNext;
  }
  return true;
}

//type : WHITE or GRAY_REVERSE for positive variable coefficient
//       BLACK or GRAY_FORWARD for negative variable coefficient
Edge * generateAbsoluteConstraint(System * system, IntegerTree * T, Vertex * x_i, int weight, EdgeType type){
  Edge * newEdge = (Edge *) malloc( sizeof(Edge) );
  newEdge->weight = weight;
  newEdge->type = type;
  newEdge->tail = x_i;
  newEdge->head = &system->graph[0];
  newEdge->reverse = NULL;
  newEdge->next = T->additionsFirst;
  T->additionsFirst = newEdge;
  newEdge->allNext = NULL;
  newEdge->allPrev = NULL;
  newEdge->inAllEdgeList = false;
  return newEdge;
}

IntegerTree * generateIntegerTree(System * system, IntegerTreeType type){
  
  IntegerTreeVertex * x0Root = (IntegerTreeVertex *) malloc( sizeof(IntegerTreeVertex) );
  x0Root->parent = NULL;
  x0Root->nextSibling = NULL;
  x0Root->firstChild = NULL;
  x0Root->queueNewer = NULL;
  x0Root->graphEdges = generateEdgeRefList();
  x0Root->graphVertex = &system->graph[0];
  system->graph[0].integerTreeVertex[type] = x0Root;
  
  IntegerTree * T = (IntegerTree *) malloc( sizeof(IntegerTree) );
  T->treeRoot = x0Root;
  T->queueNewest = NULL;
  T->queueOldest = NULL;
  T->type = type;
  T->additionsFirst = NULL;
  
  return T;
}

Vertex * pollIntegerTreeQueue(IntegerTree * tree){
  Vertex * output = NULL;
  if( tree->queueOldest != NULL ){
    output = tree->queueOldest->graphVertex;
    tree->queueOldest = tree->queueOldest->queueNewer;
    if( tree->queueOldest == NULL ){
      tree->queueNewest = NULL;
    }
  }
  return output;
}

void expandIntegerTree(System * system, IntegerTree * T, Vertex * active, Vertex * parent, Edge * edge0, Edge * edge1, Edge * edge2){
  
  IntegerTreeVertex * itv = active->integerTreeVertex[T->type];
  if( itv == NULL ){
    itv = (IntegerTreeVertex *) malloc( sizeof(IntegerTreeVertex) );
    itv->parent = parent->integerTreeVertex[T->type];
    itv->nextSibling = itv->parent->firstChild;
    itv->parent->firstChild = itv;
    itv->firstChild = NULL;
    itv->queueNewer = NULL;
    if( T->queueNewest != NULL ){
      T->queueNewest->queueNewer = itv;
    }
    else {
      T->queueOldest = itv;
    }
    T->queueNewest = itv;
    itv->graphEdges = generateEdgeRefList();
    itv->graphVertex = active;
    active->integerTreeVertex[T->type] = itv;
  }
  
  if( edge0 != NULL ){
    addEdgeToEdgeRefListEnd( itv->graphEdges, edge0 );
    if( edge1 != NULL ){
      addEdgeToEdgeRefListEnd( itv->graphEdges, edge1 );
      if( edge2 != NULL ){
        addEdgeToEdgeRefListEnd( itv->graphEdges, edge2 );
      }
    }
  }

}

void integerTreeBacktrack(EdgeRefList * list, IntegerTreeVertex * fromVertex, IntegerTreeVertex * toVertex, bool includeToVertex){
  if( fromVertex != NULL ){ //== NULL probably impossible
    while( fromVertex != toVertex ){
      copyTreeEdgesToList(list, fromVertex);
      fromVertex = fromVertex->parent;
    }
    if( includeToVertex ){
      copyTreeEdgesToList(list, fromVertex);
    }
  }
}

void copyTreeEdgesToList(EdgeRefList * list, IntegerTreeVertex * itv){
  EdgeRefListNode * treeEdge = itv->graphEdges->first;
  while( treeEdge != NULL ){
    addEdgeToEdgeRefListBeginning( list, treeEdge->edge );
    treeEdge = treeEdge->next;
  }
}

void freeIntegerTree(IntegerTree * tree){
  if( tree != NULL ){
    IntegerTreeVertex * itv = tree->treeRoot;
    while( itv != NULL ){
      IntegerTreeVertex * parent = NULL; //initialization unnecessary, but compiler might require it
      while( itv != NULL ){
        while( itv->firstChild != NULL ){
          itv = itv->firstChild;
        }
        parent = itv->parent;
        freeEdgeRefList( itv->graphEdges );
        itv->graphVertex->integerTreeVertex[ tree->type ] = NULL;
        IntegerTreeVertex * priorSibling = itv;
        itv = itv->nextSibling;
        free(priorSibling);
      }
      itv = parent;    
      if( itv != NULL ){
        itv->firstChild = NULL; //So we don't try to go back down the tree over nodes we've already freed
      }
    }
    Edge * edge = tree->additionsFirst;
    while( edge != NULL ){
      Edge * oldEdge = edge;
      edge = edge->next;
      free(oldEdge);
    }
    free(tree);
  }
}

EdgeRefList * generateEdgeRefList(){
  EdgeRefList * newERL = (EdgeRefList *) malloc( sizeof(EdgeRefList) );
  newERL->first = NULL;
  newERL->last = NULL;
  return newERL;
}

void addEdgeToEdgeRefListEnd(EdgeRefList * erl, Edge * edge){
  EdgeRefListNode * newERLN = (EdgeRefListNode *) malloc( sizeof(EdgeRefListNode) );
  newERLN->edge = edge;
  newERLN->next = NULL;
  if( erl->first == NULL ){
    erl->first = newERLN;
  }
  else {
    erl->last->next = newERLN;
  }
  erl->last = newERLN;
}

void addEdgeToEdgeRefListBeginning(EdgeRefList * erl, Edge * edge){
  EdgeRefListNode * newERLN = (EdgeRefListNode *) malloc( sizeof(EdgeRefListNode) );
  newERLN->edge = edge;
  newERLN->next = erl->first;
  erl->first = newERLN;
  if( erl->last == NULL ){
    erl->last = newERLN;
  }
}

void freeEdgeRefList(EdgeRefList * erl){
  if( erl != NULL ){
    EdgeRefListNode * erln = erl->first;
    while( erln != NULL ){
      EdgeRefListNode * oldERLN = erln;
      erln = erln->next;
      free(oldERLN);
    }
    free(erl);
  }
}

void freeSystem(System * system){
  for(EdgeType i = WHITE; i <= GRAY_REVERSE; i++){
    for(int j = 0; j < system->vertexCount; j++){
      Edge * e = system->graph[j].first[i];
      while(e != NULL){
        Edge * oldE = e;
        e = e->next;
        free( oldE );
      }
    }
  }
  freeEdgeRefList( system->infeasibilityProof );
  freeIntegerTree( system->T[DOWN] );
  freeIntegerTree( system->T[UP] );
  free( system->graph );
}
