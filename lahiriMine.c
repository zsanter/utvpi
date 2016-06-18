#include "lahiriMine.h"

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
    int infeasibleVertexIndex = lahiri(&system);
    if( infeasibleVertexIndex >= 0 ){
      fputs("\nSystem is not integrally feasible.\nProof:\n", output);
      int k = system->graph[NEGATIVE][i].D - system->graph[POSITIVE][i].D;
      fprintf(output, "-x%d <= %d\n", infeasibleVertexIndex + 1, (k-1)/2 );
      fprintf(output, "+x%d <= %d\n", infeasibleVertexIndex + 1, (-k-1)/2 );
    }
    else {
      //system integrally feasible
      //respond accordingly
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
      system->graph[i][j].dfsColor = WHITE;
      system->graph[i][j].discoveryTime = 0;
      system->graph[i][j].finishingTime = 0;
      system->graph[i][j].sccNumber = 0;
      system->graph[i][j].h = 0;
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
  inducedSubgraph(Gphi, &GphiPrime);
  stronglyConnectedComponents(&GphiPrime);
  for(int i = 0; i < Gphi->n; i++){
    if( GphiPrime.graph[POSITIVE][i].sccNumber == GphiPrime.graph[NEGATIVE][i].sccNumber
        && ( Gphi->graph[POSITIVE][i].D - Gphi->graph[NEGATIVE][i].D ) % 2 != 0 ){
      return i;
    }
  }
  freeSystem( &GphiPrime );
  johnsonAllPairs( Gphi );
  
  return -1;
  
  
}

//Doesn't copy Vertex elements L, D, x, dfsColor, discoveryTime, finishingTime, sccNumber, or h, 
//because this is unnecessary for our purposes
void inducedSubgraph(System * original, System * subgraph){
  initializeSystem( subgraph, original->n, NULL );
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < original->n; j++){
      Edge * originalEdge = original->graph[i][j].first;
      while( originalEdge != NULL ){
        if( edge->head->D - edge->tail->D == edge->weight ){
          Edge * subgraphEdge = (Edge *) malloc( sizeof(Edge) );
          subgraphEdge->weight = originalEdge->weight;
          subgraphEdge->head = &subgraph->graph[ originalEdge->head->vertexSign ][ originalEdge->head->index - 1 ];
          subgraphEdge->tail = &subgraph->graph[ originalEdge->tail->vertexSign ][ originalEdge->tail->index - 1 ];
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
  generateTranspose( system, &transpose );
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
      dfsVisit( edge->head, time );
    }
    edge = edge->next;
  }
  vertex->dfsColor = BLACK;
  (*time)++;
  vertex->finishingTime = *time;
}

//Doesn't copy Vertex elements L, D, x, dfsColor, discoveryTime, finishingTime, sccNumber, or h, 
//because this is unnecessary for our purposes
void generateTranspose(System * original, System * transpose){
  initializeSystem( transpose, original->n, NULL );
  for(VertexSign i = POSITVE; i <= NEGATIVE; i++){
    for(int j = 0; j < original->n; j++){
      Edge * originalEdge = original->graph[i][j].first;
      while( originalEdge != NULL ){
        Edge * transposeEdge = (Edge *) malloc( sizeof(Edge) );
        transposeEdge->weight = originalEdge->weight;
        transposeEdge->head = &transpose->graph[ originalEdge->tail->vertexSign ][ originalEdge->tail->index - 1 ];
        transposeEdge->tail = &transpose->graph[ originalEdge->head->vertexSign ][ originalEdge->head->index - 1 ];
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

void johnsonAllPairs(System * system){
  System copy;
  copySystem(system, &copy);
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < copy->n; j++){
      copy->graph[i][j].h = system->graph[i][j].D;
    }
  }
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < copy->n; j++){
      Edge * edge = copy->graph[i][j].first;
      while( edge != NULL ){
        edge->weight = edge->weight + edge->tail->h - edge->head->h;
      }
    }
  }
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < system->n; j++){
      dijkstra( &copy, &copy->graph[i][j] );
      Vertex * copyHead = &copy->graph[!i][j];
      int weight = copyHead->D;
      if( weight != INT_MAX ){
        if( weight % 2 != 0 ){
          weight--;
        }
        Vertex * systemHead = &system->graph[!i][j];
        Vertex * systemTail = &system->graph[i][j];
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
  
}

//Doesn't copy Vertex elements L, D, x, dfsColor, discoveryTime, finishingTime, sccNumber, or h, 
//because this is unnecessary for our purposes
void copySystem(System * original, System * copy){
  initializeSystem( copy, original->n, NULL );
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < original->n; j++){
      Edge * originalEdge = original->graph[i][j].first;
      while( originalEdge != NULL ){
        Edge * copyEdge = (Edge *) malloc( sizeof(Edge) );
        copyEdge->weight = originalEdge->weight;
        copyEdge->head = &copy->graph[ originalEdge->head->vertexSign ][ originalEdge->head->index - 1 ];
        copyEdge->tail = &copy->graph[ originalEdge->tail->vertexSign ][ originalEdge->tail->index - 1 ];
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
      system->graph[i][j].D = INT_MAX;
      system->graph[i][j].dijkstraFinalized = false;
    }
  }
  vertex->D = 0;
  
  FibHeap pQueue;
  pQueue.min = NULL;
  pQueue.n = 0;
  
  fibHeapInsert( &pQueue, vertex );
  while( pQueue.n > 0 ){
    Vertex * current = fibHeapExtractMin( &pQueue );
    current->finalized = true;
    Edge * outbound = current->first;
    while( outbound != NULL ){
      if( outbound->head->finalized == false 
          && outbound->head->D > outbound->tail->D + outbound->weight ){
        outbound->head->D = ooutbound->tail->D + outbound->weight;
        if( outbound->head->fibHeapNode == NULL ){
          fibHeapInsert( &pQueue, outbound->head );
        }
        else {
          fibHeapDecreaseKey( &pQueue, outbound->head->fibHeapNode );
        }
      }
      outbound = outbound->next;
    }
  }
  
}

void fibHeapInsert(FibHeap * fibHeap, Vertex * vertex){
  FibHeapNode * newFHN = (FibHeapNode *) malloc( sizeof(FibHeapNode) );
  newFHN->degree = 0;
  newFHN->parent = NULL;
  newFHN->child = NULL;
  newFHN->mark = false;
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
  if( fibHeap->min != NULL ){
    output = fibHeap->min->vertex;
    FibHeapNode * oldFHN = fibHeap->min;
    if( fibHeap->min->child != NULL ){
      FibHeapNode * fhn = fibHeap->min->child;
      while( fhn->parent != NULL ){
        fhn->parent == NULL;
        fhn = fhn->left;
      }
      if( fibHeap->min->right != fibHeap->min ){
        fibHeap->min->right->left = fibHeap->min->child;
        fibHeap->min->child->right->left = fibHeap->min->left;
        fibHeap->min->left->right = fibHeap->min->child->right;
        fibHeap->min->child->right = fibheap->min->right;
      }
      fibHeap->min = fibHeap->min->child;
    }
    else if( fibHeap->min->right != fibHeap->min ){
      fibHeap->min == fibHeap->min->right;
    }
    else {
      fibHeap->min = NULL;
    }
    free( oldFHN );
    fibHeap->n--;
    output->fibHeapNode = NULL;
    if( fibHeap->min != NULL ){
      fibHeapConsolidate( fibHeap );
    }
  }
  else {
    output = NULL;
  }
  return output;
}

void fibHeapConsolidate(FibHeap * fibHeap){
  double phi = (1.0 + sqrt( 5.0 ) )/2.0;
  int Alength = ((int)( log(n)/log(phi) )) + 1;
  FibHeapNode * A[ Alength ];
  for(int i = 0; i < Alength; i++){
    A[i] = NULL;
  }
  FibHeapNode * w = fibHeap->min;
  do {
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
    w = w->left;
  } while( w != fibHeap->min );
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
}

void fibHeapLink(FibHeap * fibHeap, FibHeapNode * y, FibHeapNode * x){
  y->left->right = y->right;
  y->right->left = y->left;
  y->right = x->child->right;
  x->child->right->left = y;
  y->left = x->child;
  x->child->right = y;
  x->degree++;
  y->mark = false;
}

void fibHeapDecreaseKey(FibHeap * fibHeap, FibHeapNode * x){
  FibHeapNode * y = x->parent;
  if( y != NULL && x->vertex->D < y->vertex->D ){
    fibHeapCut(fibHeap, x, y);
    fibHeapCascadingCut(fibHeap, y);
  }
  if( x->vertex->D < fibHeap->min->vretex->D ){
    fibHeap->min = x;
  }
}

void fibHeapCut(FibHeap * fibHeap, FibHeapNode * x, FibHeapNode * y){
  if( y->child == y->child->right ){
    y->child == NULL;
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
      fibHeapCut(fibHeap. y, z);
      fibHeapCascadingCut(fibHeap, z);
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