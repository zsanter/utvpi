#include "lahiriMine.h"
#include <string.h>

void testFibHeap();

int main(int argc, char * argv[]){
  //printf("NEGATIVE = %d, POSITIVE = %d, !NEGATIVE = %d , !POSITIVE = %d\n", NEGATIVE, POSITIVE, !NEGATIVE, !POSITIVE );
  //printf("NULL = %p\n", NULL );
  clock_t beginning = clock();
  if( argc != 2 && argc != 3 ){
    fprintf( stderr, "Proper use is %s [input file] {output file}.\nIf no output file is specified, output is to stdout.\n", argv[0] );
    exit(1);
  }
  if( strcmp( argv[1], "testFibHeap") == 0 ){
    testFibHeap();
    return 0;
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
      double solution = ((double)( system.graph[POSITIVE][i].D - system.graph[NEGATIVE][i].D )) / 2.0;
      fprintf( output, "x%i = %.1f\n", i + 1, solution );
    }
    int infeasibleVertexIndex = lahiri(&system);
    if( infeasibleVertexIndex >= 0 ){
      fputs("\nSystem is not integrally feasible.\nProof:\n", output);
      int k = system.graph[NEGATIVE][infeasibleVertexIndex].D - system.graph[POSITIVE][infeasibleVertexIndex].D;
      fprintf(output, "-x%d <= %d\n", infeasibleVertexIndex + 1, (k-1)/2 );
      fprintf(output, "+x%d <= %d\n", infeasibleVertexIndex + 1, (-k-1)/2 );
    }
    else {
      fputs("\nIntegral solution:\n", output);
      for(int i = 0; i < system.n; i++){
        fprintf(output, "x%i = %i\n", i + 1, system.graph[POSITIVE][i].rho );
      }
    }
  }
  fclose(output);
  time_t afterSolutionOutput = clock();
  freeSystem(&system);
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
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    system->graph[i] = (Vertex *) malloc( sizeof(Vertex) * n );
    for(int j = 0; j < system->n; j++){
      system->graph[i][j].index = j + 1;
      system->graph[i][j].sign = i;
      system->graph[i][j].L = NULL;
      system->graph[i][j].D = 0;
      system->graph[i][j].first = NULL;
      system->graph[i][j].rho = INT_MAX;
      system->graph[i][j].dfsColor = WHITE;
      system->graph[i][j].discoveryTime = 0;
      system->graph[i][j].finishingTime = 0;
      system->graph[i][j].sccNumber = 0;
      system->graph[i][j].h = 0;
      system->graph[i][j].dijkstraFinalized = false;
      system->graph[i][j].fibHeapNode = NULL;
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

int lahiri(System * Gphi){
  System GphiPrime;
  //puts("1");
  onlySlacklessEdges(Gphi, &GphiPrime);
  //puts("2");
  stronglyConnectedComponents(&GphiPrime);
  //puts("3");
  int output = -1;
  for(int i = 0; i < Gphi->n && output < 0; i++){
    if( GphiPrime.graph[POSITIVE][i].sccNumber == GphiPrime.graph[NEGATIVE][i].sccNumber
        && ( Gphi->graph[POSITIVE][i].D - Gphi->graph[NEGATIVE][i].D ) % 2 != 0 ){
      output = i;
    }
  }
  //puts("4");
  freeSystem( &GphiPrime );
  if( output >= 0 ){
    return output;
  }
  //puts("5");
  johnsonAllPairs( Gphi );
  //puts("6");
  noHeadIndicesHigherThanTailIndeces( Gphi );
  //puts("7");
  for(int j = 0; j < Gphi->n; j++){
    int upperBound = INT_MAX;
    int lowerBound = INT_MIN;
    for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
      Edge * edge = Gphi->graph[i][j].first;
      while( edge != NULL ){
        int potentialNewUpperBound = INT_MAX;
        int potentialNewLowerBound = INT_MIN;
        if( edge->tail->index == edge->head->index ){
          //char sign;
          switch( edge->tail->sign ){
          case POSITIVE:
            potentialNewLowerBound = -edge->weight/2;
            //sign = '-';
            break;
          case NEGATIVE:
            potentialNewUpperBound = edge->weight/2;
            //sign = '+';
          }
          //fprintf(output, "%cx%d <= %d\n", sign, edge->tail->index, (edge->weight)/2);
        }
        else{
          //char sign[2];
          if( edge->tail->sign == edge->head->sign ){
            switch( edge->tail->sign ){
            case POSITIVE:
              potentialNewLowerBound = -edge->weight + edge->head->rho;
              //sign[0] = '-';
              //sign[1] = '+';
              break;
            case NEGATIVE:
              potentialNewUpperBound = edge->weight + edge->head->rho;
              //sign[0] = '+';
              //sign[1] = '-';
            }
          }
          else{
            switch( edge->tail->sign ){
            case POSITIVE:
              potentialNewLowerBound = -edge->weight - edge->head->rho;
              //sign[0] = '-';
              //sign[1] = '-';
              break;
            case NEGATIVE:
              potentialNewUpperBound = edge->weight - edge->head->rho;
              //sign[0] = '+';
              //sign[1] = '+';
            }
          }
          //fprintf(output, "%cx%d %cx%d <= %d\n", sign[0], edge->tail->index, sign[1], edge->head->index, edge->weight);
        }
        if( potentialNewUpperBound < upperBound ){
          upperBound = potentialNewUpperBound;
        }
        if( potentialNewLowerBound > lowerBound ){
          lowerBound = potentialNewLowerBound;
        }
        
        edge = edge->next;
      }  
    }
    int solution;
    if( upperBound == INT_MAX && lowerBound == INT_MIN ){
      solution = 0;
    }
    else if( upperBound == INT_MAX ){
      solution = lowerBound;
    }
    else if( lowerBound == INT_MIN ){
      solution = upperBound;
    }
    else {
      solution = upperBound/2 + lowerBound/2;
    }
    Gphi->graph[POSITIVE][j].rho = solution;
    Gphi->graph[NEGATIVE][j].rho = solution;
  }
  
  return output;
}

/*
 * Doesn't copy Vertex elements L, D, rho, dfsColor, discoveryTime, finishingTime, sccNumber, h, 
 * dijkstraFinalized, or fibHeapNode, because this is unnecessary for our purposes
 */
void onlySlacklessEdges(System * original, System * subgraph){
  initializeSystem( subgraph, original->n, NULL );
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < original->n; j++){
      Edge * originalEdge = original->graph[i][j].first;
      while( originalEdge != NULL ){
        if( originalEdge->head->D - originalEdge->tail->D == originalEdge->weight ){
          Edge * subgraphEdge = (Edge *) malloc( sizeof(Edge) );
          subgraphEdge->weight = originalEdge->weight;
          subgraphEdge->head = &subgraph->graph[ originalEdge->head->sign ][ originalEdge->head->index - 1 ];
          subgraphEdge->tail = &subgraph->graph[ originalEdge->tail->sign ][ originalEdge->tail->index - 1 ];
          subgraphEdge->next = subgraphEdge->tail->first;
          subgraphEdge->tail->first = subgraphEdge;
          subgraphEdge->backtrackSeen = false;
        }
        originalEdge = originalEdge->next;
      }
    }
  }
}

void stronglyConnectedComponents(System * system){
  int time = 0;
  int sccNumber = 0;
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < system->n; j++){
      if( system->graph[i][j].dfsColor == WHITE ){
        dfsVisit( &system->graph[i][j], &time, sccNumber );
      }
    }
  }
  System transpose;
  transposeSystem( system, &transpose );
  int vertexSortArrayLength = system->n * VERTEX_SIGN_COUNT;
  Vertex * vertexSortArray[ vertexSortArrayLength ];
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < system->n; j++){
      transpose.graph[i][j].finishingTime = system->graph[i][j].finishingTime;
      vertexSortArray[ i * system->n + j ] = &transpose.graph[i][j];
    }
  }
  qsort( vertexSortArray, vertexSortArrayLength, sizeof(Vertex *), vertexCompareFinishingTimes);
  time = 0;
  for(int i = vertexSortArrayLength - 1; i >= 0; i--){
    if( vertexSortArray[i]->dfsColor == WHITE ){
      dfsVisit( vertexSortArray[i], &time, sccNumber );
      sccNumber++;
    }
  }
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < system->n; j++){
      system->graph[i][j].sccNumber = transpose.graph[i][j].sccNumber;
    }
  }
  freeSystem( &transpose );
}

void dfsVisit(Vertex * vertex, int * time, int sccNumber){
  (*time)++;
  vertex->discoveryTime = *time;
  vertex->dfsColor = GRAY;
  vertex->sccNumber = sccNumber;
  Edge * edge = vertex->first;
  while( edge != NULL ){
    if( edge->head->dfsColor == WHITE ){
      dfsVisit( edge->head, time, sccNumber );
    }
    edge = edge->next;
  }
  vertex->dfsColor = BLACK;
  (*time)++;
  vertex->finishingTime = *time;
}

/*
 * Doesn't copy Vertex elements L, D, rho, dfsColor, discoveryTime, finishingTime, sccNumber, h, 
 * dijkstraFinalized, or fibHeapNode, because this is unnecessary for our purposes
 */
void transposeSystem(System * original, System * transpose){
  initializeSystem( transpose, original->n, NULL );
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < original->n; j++){
      Edge * originalEdge = original->graph[i][j].first;
      while( originalEdge != NULL ){
        Edge * transposeEdge = (Edge *) malloc( sizeof(Edge) );
        transposeEdge->weight = originalEdge->weight;
        transposeEdge->head = &transpose->graph[ originalEdge->tail->sign ][ originalEdge->tail->index - 1 ];
        transposeEdge->tail = &transpose->graph[ originalEdge->head->sign ][ originalEdge->head->index - 1 ];
        transposeEdge->next = transposeEdge->tail->first;
        transposeEdge->tail->first = transposeEdge;
        transposeEdge->backtrackSeen = false;
        
        originalEdge = originalEdge->next;
      }
    }
  }
}

int vertexCompareFinishingTimes(const void * vertex1, const void * vertex2){
  return (*(Vertex **)vertex1)->finishingTime - (*(Vertex **)vertex2)->finishingTime;
}

/*
 * Adds required absolute constraints to system, instead of creating
 * an n x n array whose entries will mostly be ignored.
 */
void johnsonAllPairs(System * system){
  System copy;
  copySystem(system, &copy);
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < copy.n; j++){
      copy.graph[i][j].h = system->graph[i][j].D;
    }
  }
  //puts("5.1");
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < copy.n; j++){
      Edge * edge = copy.graph[i][j].first;
      while( edge != NULL ){
        edge->weight += edge->tail->h - edge->head->h;
	if( edge->weight < 0 ){
	  fputEdge( edge, stdout );
	}
        edge = edge->next;
      }
    }
  }
  //puts("5.2");
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < system->n; j++){
      Vertex * copyTail = &copy.graph[i][j];
      Vertex * copyHead = &copy.graph[!i][j];
      dijkstra( &copy, copyTail );
      int weight = copyHead->D;
      if( weight != INT_MAX ){
        weight += copyHead->h - copyTail->h;
        if( weight % 2 != 0 ){
          weight--;
        }
        Vertex * systemTail = &system->graph[i][j];
        Vertex * systemHead = &system->graph[!i][j];
        Edge * newEdge = (Edge *) malloc( sizeof(Edge) );
        newEdge->weight = weight;
        newEdge->head = systemHead;
        newEdge->tail = systemTail;
        newEdge->next = systemTail->first;
        systemTail->first = newEdge;
        newEdge->backtrackSeen = false;
      }
    }
  }
  freeSystem( &copy );
}

/*
 * Doesn't copy Vertex elements L, D, rho, dfsColor, discoveryTime, finishingTime, sccNumber, h, 
 * dijkstraFinalized, or fibHeapNode, because this is unnecessary for our purposes
 */
void copySystem(System * original, System * copy){
  initializeSystem( copy, original->n, NULL );
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < original->n; j++){
      Edge * originalEdge = original->graph[i][j].first;
      while( originalEdge != NULL ){
        Edge * copyEdge = (Edge *) malloc( sizeof(Edge) );
        copyEdge->weight = originalEdge->weight;
        copyEdge->head = &copy->graph[ originalEdge->head->sign ][ originalEdge->head->index - 1 ];
        copyEdge->tail = &copy->graph[ originalEdge->tail->sign ][ originalEdge->tail->index - 1 ];
        copyEdge->next = copyEdge->tail->first;
        copyEdge->tail->first = copyEdge;
        copyEdge->backtrackSeen = false;
        
        originalEdge = originalEdge->next;
      }
    }
  }
}

void dijkstra(System * system, Vertex * vertex){
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < system->n; j++){
      system->graph[i][j].D = INT_MAX; //Initialized to 0 to avoid having an actual source node for Bellman-Ford.
      system->graph[i][j].dijkstraFinalized = false;
    }
  }
  vertex->D = 0;
  
  FibHeap pQueue;
  pQueue.min = NULL;
  pQueue.n = 0;
  //puts("5.2.1");
  fibHeapInsert( &pQueue, vertex );
  //puts("5.2.2");
  int i3 = 0;
  while( pQueue.n > 0 ){
    printf("i3 = %d, pQueue.n = %d\n", i3, pQueue.n );
    i3++;
    //puts("5.2.3");
    Vertex * current = fibHeapExtractMin( &pQueue );
    if( current == NULL ){
      puts("current == NULL");
    }
    //puts("5.2.4");
    current->dijkstraFinalized = true;
    Edge * outbound = current->first;
    while( outbound != NULL ){
      if( outbound->head->dijkstraFinalized == false 
          && outbound->head->D > outbound->tail->D + outbound->weight ){
        outbound->head->D = outbound->tail->D + outbound->weight;
        if( outbound->head->fibHeapNode == NULL ){
	  //puts("5.2.5");
	  fibHeapInsert( &pQueue, outbound->head );
	  //puts("5.2.6");
        }
        else {
	  //puts("5.2.7");
          fibHeapDecreaseKey( &pQueue, outbound->head );
	  //puts("5.2.8");
	}
      }
      outbound = outbound->next;
    }
  }
  
}

void testFibHeap(){
  Vertex vertices[100];

  FibHeap fibHeap;
  fibHeap.min = NULL;
  fibHeap.n = 0;

  char command[100];
  bool exit = false;

  puts("new [index] [D]\ndecrease [index] [new D]\nextract\nexit\nindeces 0-99 only");
  
  while( !exit ){
    fgets(command, 100, stdin);
    char * word = strtok(command, " \t\n\r");
    if( strcmp( word, "new") == 0 ){
      int index = atoi( strtok( NULL, " \t\n\r") );
      int D = atoi( strtok( NULL, " \t\n\r"  ) );
      vertices[index].D = D;
      vertices[index].index = index;
      vertices[index].sign = POSITIVE;
      fibHeapInsert(&fibHeap, &vertices[index]);
    }
    else if( strcmp( word, "decrease" ) == 0 ){
      int index = atoi( strtok( NULL, " \t\n\r") );
      int D = atoi( strtok( NULL, " \t\n\r"  ) );
      vertices[index].D = D;
      fibHeapDecreaseKey(&fibHeap, &vertices[index]);
    }
    else if( strcmp( word, "extract" ) == 0 ){
      Vertex * vertex = fibHeapExtractMin( &fibHeap );
      printf("Extracted x%d, D = %d\n", vertex->index, vertex->D);
    }
    else if( strcmp( word, "exit") == 0){
      exit = true;
    }
    if( fibHeap.min > 0 ){
      printf("min = x%d, D = %d\n", fibHeap.min->vertex->index, fibHeap.min->vertex->D);
    }
    else {
      puts("fibHeap shows empty");
    }
  }
  puts("Sure hope this accomplished something.");
}

void fibHeapInsert(FibHeap * fibHeap, Vertex * vertex){
  printf("fibHeapInsert() inserting %cx%d - %p\n", /*vertex->sign*/ (vertex->sign == POSITIVE) ? '+' : '-' , vertex->index, vertex );
  FibHeapNode * newFHN = (FibHeapNode *) malloc( sizeof(FibHeapNode) );
  newFHN->degree = 0;
  newFHN->parent = NULL;
  newFHN->child = NULL;
  newFHN->mark = false;
  newFHN->rootListTraverseSeen = false;
  newFHN->vertex = vertex;
  vertex->fibHeapNode = newFHN;
  if( fibHeap->min == NULL ){
    newFHN->right = newFHN;
    newFHN->left = newFHN;
    fibHeap->min = newFHN;
  }
  else {
    fibHeap->min->right->left = newFHN;
    newFHN->right = fibHeap->min->right;
    fibHeap->min->right = newFHN;
    newFHN->left = fibHeap->min;
    if( vertex->D < fibHeap->min->vertex->D ){
      fibHeap->min = newFHN;
    }
  }
  fibHeap->n++;
}

Vertex * fibHeapExtractMin(FibHeap * fibHeap){
  Vertex * output;
  int i1;
  if( fibHeap->min != NULL ){
    output = fibHeap->min->vertex;
    FibHeapNode * oldFHN = fibHeap->min;
    printf("fibHeapExtractMin() extracting %cx%d - %p\n", /*oldFHN->vertex->sign*/(oldFHN->vertex->sign == POSITIVE) ? '+' : '-' , oldFHN->vertex->index , output );
    FibHeapNode * fhn = fibHeap->min;
    puts("original main heap list:");
    i1 = 0;
    do {
      printf("%cx%d - %p > ", (fhn->vertex->sign == POSITIVE) ? '+' : '-' , fhn->vertex->index , fhn );
      fhn = fhn->left;
      i1++;
    } while( fhn != fibHeap->min && i1 < 20 );
    putchar('\n');
    if( fibHeap->min->child != NULL ){
      puts("if( fibHeap->min->child != NULL )");
      FibHeapNode * fhn = fibHeap->min->child;
      puts("original child list:");
      puts("5.2.3.1");
      while( fhn->parent != NULL ){
	printf("%cx%d - %p > ", (fhn->vertex->sign == POSITIVE) ? '+' : '-' , fhn->vertex->index , fhn );
        fhn->parent == NULL;
        fhn = fhn->left;
      }
      putchar('\n');
      if( i1 >= 20 ){
	exit(1);
      }
      puts("5.2.3.2");
      if( fibHeap->min->right != fibHeap->min ){
	puts("if( fibHeap->min->right != fibHeap->min )");
        fibHeap->min->right->left = fibHeap->min->child;
        fibHeap->min->child->right->left = fibHeap->min->left;
        fibHeap->min->left->right = fibHeap->min->child->right;
        fibHeap->min->child->right = fibHeap->min->right;
      }
      fibHeap->min = fibHeap->min->child;
    }
    else if( fibHeap->min->right != fibHeap->min ){
      /*puts("else if( fibHeap->min->right != fibHeap->min )");
      printf("fibHeap->min = %cx%d\n", (fibHeap->min->vertex->sign == POSITIVE) ? '+' : '-' , fibHeap->min->vertex->index);
      printf("fibHeap->min->right = %cx%d\n", (fibHeap->min->right->vertex->sign == POSITIVE) ? '+' : '-' , fibHeap->min->right->vertex->index);
      printf("fibHeap->min->left = %cx%d\n", (fibHeap->min->left->vertex->sign == POSITIVE) ? '+' : '-' , fibHeap->min->left->vertex->index);
      printf("fibHeap->min->right->left = %cx%d\n", (fibHeap->min->right->left->vertex->sign == POSITIVE) ? '+' : '-' , fibHeap->min->right->left->vertex->index);
      printf("fibHeap->min->left->right = %cx%d\n", (fibHeap->min->left->right->vertex->sign == POSITIVE) ? '+' : '-' , fibHeap->min->left->right->vertex->index);*/
      fibHeap->min->right->left = fibHeap->min->left;
      fibHeap->min->left->right = fibHeap->min->right;
      fibHeap->min = fibHeap->min->right;
      /*printf("fibHeap->min = %cx%d\n", (fibHeap->min->vertex->sign == POSITIVE) ? '+' : '-' , fibHeap->min->vertex->index);
      printf("fibHeap->min->right = %cx%d\n", (fibHeap->min->right->vertex->sign == POSITIVE) ? '+' : '-' , fibHeap->min->right->vertex->index);
      printf("fibHeap->min->left = %cx%d\n", (fibHeap->min->left->vertex->sign == POSITIVE) ? '+' : '-' , fibHeap->min->left->vertex->index);*/
    }
    else {
      puts("else");
      fibHeap->min = NULL;
    }
    output->fibHeapNode = NULL;
    if( fibHeap->min != NULL ){
      i1 = 0;
      puts("before consolidate main heap list:");
      fhn = fibHeap->min;
      do {
	printf("%cx%d - %p > ", (fhn->vertex->sign == POSITIVE) ? '+' : '-' , fhn->vertex->index , fhn );
	fhn = fhn->left;
	i1++;
      } while( fhn != fibHeap->min && i1 < 20 );
      putchar('\n');
      if( i1 >= 20 ){
	exit(1);
      }
      fibHeapConsolidate( fibHeap );
      i1 = 0;
      puts("after consolidate main heap list:");
      fhn = fibHeap->min;
      do {
	printf("%cx%d - %p > ", (fhn->vertex->sign == POSITIVE) ? '+' : '-' , fhn->vertex->index , fhn );
	fhn = fhn->left;
	i1++;
      } while( fhn != fibHeap->min && i1 < 20 );
      putchar('\n');
      if( i1 >= 20 ){
	exit(1);
      }
    }
    free( oldFHN );
    fibHeap->n--;
  }
  else {
    puts( "fibHeapExtractMin() returns NULL");
    output = NULL;
  }
  return output;
}

void fibHeapConsolidate(FibHeap * fibHeap){
  double phi = ( 1.0 + sqrt( 5.0 ) ) / 2.0;
  int Alength = ((int)( log( fibHeap->n ) / log( phi ) )) + 1;
  FibHeapNode * A[ Alength ];
  for(int i = 0; i < Alength; i++){
    A[i] = NULL;
  }
  FibHeapNode * w = fibHeap->min;
  puts("5.2.3.2.1");
  int i2 = 0;
  do {
    printf("%cx%d - %p > ", (w->vertex->sign == POSITIVE) ? '+' : '-' , w->vertex->index , w );
    w->rootListTraverseSeen = false;
    w = w->left;
    i2++;
  } while( w != fibHeap->min && i2 < 20 );
  if( i2 >= 20 ){
    exit(1);
  }
  putchar('\n');
  puts("5.2.3.2.2");
  while( w->rootListTraverseSeen == false ){
    w->rootListTraverseSeen = true;
    FibHeapNode * nextW = w->left;
    FibHeapNode * x = w;
    int d = x->degree;
    while( A[d] != NULL ){
      FibHeapNode * y = A[d];
      if( x->vertex->D > y->vertex->D ){
        FibHeapNode * temp = y;
        y = x;
        x = temp;
      }
      fibHeapLink(fibHeap, y, x);
      A[d] = NULL;
      d++;
    }
    A[d] = x;
    w = nextW;
  }
  puts("5.2.3.2.3");
  fibHeap->min = NULL;
  for(int i = 0; i < Alength; i++){
    if( A[i] != NULL ){
      if( fibHeap->min == NULL ){
        A[i]->parent = NULL;
        A[i]->left = A[i];
        A[i]->right = A[i];
        fibHeap->min = A[i];
      }
      else {
        fibHeap->min->right->left = A[i];
        A[i]->right = fibHeap->min->right;
        fibHeap->min->right = A[i];
        A[i]->left = fibHeap->min;
        if( A[i]->vertex->D < fibHeap->min->vertex->D ){
          fibHeap->min = A[i];
        }
      }
    }
  }
  if( fibHeap->min == NULL ){
    puts("fibHeap->min == NULL at end of consolidate()");
    exit(9);
  }
  puts("5.2.3.2.4");
}

void fibHeapLink(FibHeap * fibHeap, FibHeapNode * y, FibHeapNode * x){
  y->left->right = y->right;
  y->right->left = y->left;
  //printf("x = %p, x->child = %p\n", x, x->child);
  if( x->child != NULL ){
    y->right = x->child->right;
    x->child->right->left = y;
    y->left = x->child;
    x->child->right = y;
  }
  else {
    y->right = y;
    y->left = y;
    x->child = y;
  }
  x->degree++;
  y->mark = false;
}

/*
 * Will silently fail if vertex->D is actually set higher than it had been.
 * This scenario can cause the heap to violate the min-heap property, making it useless. 
 */
void fibHeapDecreaseKey(FibHeap * fibHeap, Vertex * vertex){
  FibHeapNode * x = vertex->fibHeapNode;
  FibHeapNode * y = x->parent;
  if( y != NULL && x->vertex->D < y->vertex->D ){
    fibHeapCut(fibHeap, x, y);
    fibHeapCascadingCut(fibHeap, y);
  }
  if( x->vertex->D < fibHeap->min->vertex->D ){
    fibHeap->min = x;
  }
}

void fibHeapCut(FibHeap * fibHeap, FibHeapNode * x, FibHeapNode * y){
  if( y->child == y->child->right ){
    y->child = NULL;
  }
  else if( y->child == x ){
    y->child = y->child->right;
  }
  x->left->right = x->right;
  x->right->left = x->left;
  x->right = fibHeap->min->right;
  fibHeap->min->right->left = x;
  x->left = fibHeap->min;
  fibHeap->min->right = x;
  y->degree--;
  x->parent = NULL;
  x->mark = false;
}

void fibHeapCascadingCut(FibHeap * fibHeap, FibHeapNode * y){
  FibHeapNode * z = y->parent;
  if( z != NULL ){
    if( y->mark == false ){
      y->mark = true;
    }
    else {
      fibHeapCut(fibHeap, y, z);
      fibHeapCascadingCut(fibHeap, z);
    }
  }
}

//Cut out all edges where head->index > tail->index
void noHeadIndicesHigherThanTailIndeces( System * system ){
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < system->n; j++){
      Edge * prior = NULL;
      Edge * edge = system->graph[i][j].first;
      while( edge != NULL ){
        Edge * next = edge->next;
        if( edge->tail->index >= edge->head->index ){
          if( prior == NULL ){
            system->graph[i][j].first = edge;
          }
          else {
            prior->next = edge;
          }
          prior = edge;
        }
        else {
          free(edge);
        }
        edge = next;
      }
      if( prior == NULL ){
        system->graph[i][j].first = NULL;
      }
      else {
        prior->next = NULL;
      }
    }
  }
}

void freeSystem(System * system){
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
