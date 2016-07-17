#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include "utvpiInterpreter.h"

#define VERTEX_SIGN_COUNT 2
typedef enum VertexSign {
  POSITIVE,
  NEGATIVE,
} VertexSign;

typedef enum DFScolor {
  WHITE,
  GRAY,
  BLACK,
} DFScolor;

typedef struct System System;
typedef struct Vertex Vertex;
typedef struct Edge Edge;
typedef struct FibHeap FibHeap;
typedef struct FibHeapNode FibHeapNode;

struct System {
  Vertex * graph[VERTEX_SIGN_COUNT];
  int n;
  int falsePositives;
  int mainLoopIterations;
  int negativeCycleEdgeCount;
};

struct Vertex {
  int index;
  VertexSign sign;
  Edge * L;
  int D;
  Edge * cycleOriginator;
  Edge * first;
  union {
    struct {
      DFScolor dfsColor;
      int finishingTime;
      int sccNumber;
    };
    struct {
      int h;
      int rho;
      FibHeapNode * fibHeapNode;
    };
  };
};

struct Edge {
  int weight;
  Vertex * head;
  Vertex * tail;
  Edge * next;
  bool backtrackSeen;
};

struct FibHeap {
  FibHeapNode * min;
  int n;
};

struct FibHeapNode {
  FibHeapNode * parent;
  FibHeapNode * left;
  FibHeapNode * right;
  FibHeapNode * child;
  Vertex * vertex;
  int degree;
  bool mark;
  bool rootListTraverseSeen;
};

int main(int argc, char * argv[]);
static void fputEdge(Edge * edge, FILE * output);
static void initializeSystem(void * object, int n, Parser * parser);
static void setSystemForJohnson(System * system);
static void addConstraint(void * object, Constraint * constraint, Parser * parser);
static Edge * bellmanFord(System * system);
static Edge * backtrack(Edge * edge);
static int lahiri(System * Gphi);
static void onlySlacklessEdges(System * original, System * subgraph);
static void stronglyConnectedComponents(System * system);
static void dfsVisit(Vertex * vertex, int * time, int sccNumber);
static void transposeSystem(System * original, System * transpose);
static int vertexCompareFinishingTimes(const void * vertex1, const void * vertex2);
static void johnsonAllPairs(System * Gphi);
static void dijkstra(System * system, Vertex * vertex);
static void fibHeapInsert(FibHeap * fibHeap, Vertex * vertex);
static Vertex * fibHeapExtractMin(FibHeap * fibHeap);
static void fibHeapConsolidate(FibHeap * fibHeap);
static void fibHeapLink(FibHeap * fibHeap, FibHeapNode * y, FibHeapNode * x);
static void fibHeapDecreaseKey(FibHeap * fibHeap, Vertex * vertex);
static void fibHeapCut(FibHeap * fibHeap, FibHeapNode * x, FibHeapNode * y);
static void fibHeapCascadingCut(FibHeap * fibHeap, FibHeapNode * y);
static void freeSystem(System * system);

#ifdef __HPC__
  static void diff(struct timespec * start, struct timespec * end, struct timespec * difference);
#endif

int main(int argc, char * argv[]){
  #ifdef __HPC__
    struct timespec start;
    clock_gettime(CLOCK_MONOTONIC_RAW, &start);
  #else
    clock_t start = clock();
  #endif
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
  #ifdef __HPC__
    struct timespec beforeLinear;
    clock_gettime(CLOCK_MONOTONIC_RAW, &beforeLinear);
  #else
    clock_t beforeLinear = clock();
  #endif
  Edge * negativeCycle = bellmanFord(&system);
  #ifdef __HPC__
    struct timespec beforeIntegral;
    clock_gettime(CLOCK_MONOTONIC_RAW, &beforeIntegral);
  #else
    clock_t beforeIntegral = clock();
  #endif
  int f;
  if( negativeCycle != NULL ){
    f = 0;
    int negativeCycleEdgeCount = 0;
    fputs("The following negative cost cycle was detected:\n", output);
    Edge * edge = negativeCycle;
    while( edge->backtrackSeen == true ){
      fputEdge( edge, output );
      edge->backtrackSeen = false;
      edge = edge->tail->L;
      system.negativeCycleEdgeCount++;
    }
    fprintf( output, "\n%d false positives\n", system.falsePositives );
    fprintf( output, "%d main loop iterations\n", system.mainLoopIterations );
    fprintf( output, "%d negative cycle edges\n", system.negativeCycleEdgeCount );
  }
  else{
    fputs("Linear solution:\n", output);
    for(int i = 0; i < system.n; i++){
      double solution = ((double)( system.graph[POSITIVE][i].D - system.graph[NEGATIVE][i].D )) / 2.0;
      fprintf( output, "x%i = %.1f\n", i + 1, solution );
    }
    fprintf( output, "\n%d false positives\n", system.falsePositives );
    fprintf( output, "%d main loop iterations\n", system.mainLoopIterations );
    int infeasibleVertexIndex = lahiri(&system);
    if( infeasibleVertexIndex >= 0 ){
      f = 1;
      fputs("\nSystem is not integrally feasible.\nProof:\n", output);
      int k = system.graph[NEGATIVE][infeasibleVertexIndex].D - system.graph[POSITIVE][infeasibleVertexIndex].D;
      fprintf(output, "-x%d <= %d\n", infeasibleVertexIndex + 1, (k-1)/2 );
      fprintf(output, "+x%d <= %d\n", infeasibleVertexIndex + 1, (-k-1)/2 );
    }
    else {
      f = 2;
      fputs("\nIntegral solution:\n", output);
      for(int i = 0; i < system.n; i++){
        fprintf(output, "x%i = %i\n", i + 1, system.graph[POSITIVE][i].rho );
      }
    }
  }
  #ifdef __HPC__
    struct timespec beforeCleanup;
    clock_gettime(CLOCK_MONOTONIC_RAW, &beforeCleanup);
  #else
    clock_t beforeCleanup = clock();
  #endif
  fclose(output);
  freeSystem(&system);
  #ifdef __HPC__
    struct timespec end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
  #else
    clock_t end = clock();
  #endif
  printf("%d,", f);
  printf("%d,", system.falsePositives);
  printf("%d,", system.mainLoopIterations);
  printf("%d,", system.negativeCycleEdgeCount);
  #ifdef __HPC__
    struct timespec setup;
    diff(&start, &beforeLinear, &setup);
    printf("%d.%09d,", (int)setup.tv_sec, (int)setup.tv_nsec);
    struct timespec linear;
    diff(&beforeLinear, &beforeIntegral, &linear);
    printf("%d.%09d,", (int)linear.tv_sec, (int)linear.tv_nsec);
    struct timespec integral;
    diff(&beforeIntegral, &beforeCleanup, &integral);
    printf("%d.%09d,", (int)integral.tv_sec, (int)integral.tv_nsec);
    struct timespec cleanup;
    diff(&beforeCleanup, &end, &cleanup);
    printf("%d.%09d,", (int)cleanup.tv_sec, (int)cleanup.tv_nsec);
    struct timespec total;
    diff(&start, &end, &total);
    printf("%d.%09d,", (int)total.tv_sec, (int)total.tv_nsec);
  #else
    printf("%f,", ((double)(beforeLinear - start))/CLOCKS_PER_SEC);
    printf("%f,", ((double)(beforeIntegral - beforeLinear))/CLOCKS_PER_SEC);
    printf("%f,", ((double)(beforeCleanup - beforeIntegral))/CLOCKS_PER_SEC);
    printf("%f,", ((double)(end - beforeCleanup))/CLOCKS_PER_SEC);
    printf("%f,", ((double)(end - start))/CLOCKS_PER_SEC);
  #endif
  return 0;
}

#ifdef __HPC__
  /*
   * Copied from https://www.guyrutenberg.com/2007/09/22/profiling-code-using-clock_gettime/
   * and modified.
   */
  static void diff(struct timespec * start, struct timespec * end, struct timespec * difference){
    if ( ( end->tv_nsec - start->tv_nsec ) < 0 ) {
      difference->tv_sec = end->tv_sec - start->tv_sec - 1;
      difference->tv_nsec = 1000000000 + end->tv_nsec - start->tv_nsec;
    }
    else {
      difference->tv_sec = end->tv_sec - start->tv_sec;
      difference->tv_nsec = end->tv_nsec - start->tv_nsec;
    }
  }
#endif

static void fputEdge(Edge * edge, FILE * output){
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

static void initializeSystem(void * object, int n, Parser * parser){
  System * system = (System *) object;
  system->n = n;
  system->falsePositives = 0;
  system->mainLoopIterations = INT_MAX;
  system->negativeCycleEdgeCount = 0;
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    system->graph[i] = (Vertex *) malloc( sizeof(Vertex) * n );
    for(int j = 0; j < system->n; j++){
      system->graph[i][j].index = j + 1;
      system->graph[i][j].sign = i;
      system->graph[i][j].L = NULL;
      system->graph[i][j].cycleOriginator = NULL;
      system->graph[i][j].D = 0;
      system->graph[i][j].first = NULL;
      system->graph[i][j].dfsColor = WHITE;
      system->graph[i][j].finishingTime = 0;
      system->graph[i][j].sccNumber = 0;
    }
  }
}

static void setSystemForJohnson(System * system){
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < system->n; j++){
      system->graph[i][j].h = 0;
      system->graph[i][j].fibHeapNode = NULL;
      system->graph[i][j].rho = INT_MAX;
    }
  }
}

static void addConstraint(void * object, Constraint * constraint, Parser * parser){
  System * system = (System *) object;
  if( constraint->sign[1] == CONSTRAINT_NONE ){
    Edge * newEdge = (Edge *) malloc( sizeof(Edge) );
    newEdge->weight = 2 * constraint->weight;
    newEdge->backtrackSeen = false;
    VertexSign tailSign, headSign;
    if( constraint->sign[0] == CONSTRAINT_PLUS ){
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
      if( constraint->sign[0] == CONSTRAINT_PLUS ){
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
    else if( constraint->sign[0] == CONSTRAINT_PLUS ){
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

static Edge * bellmanFord(System * system){
  bool anyChange = true;
  for(int i = 1; i <= (2 * system->n - 1) && anyChange; i++){
    anyChange = false;
    for(VertexSign j = POSITIVE; j <= NEGATIVE; j++){
      for(int k = 0; k < system->n; k++){
        Edge * edge = system->graph[j][k].first;
        while(edge != NULL){
          if( edge->head->D > edge->tail->D + edge->weight ){
            Edge * priorL = edge->head->L;
            edge->head->D = edge->tail->D + edge->weight;
            edge->head->L = edge;
            if( priorL != edge ){
              edge->head->cycleOriginator = edge;
            }
            else {
              edge->head->cycleOriginator = edge->tail->cycleOriginator;
              if( edge == edge->head->cycleOriginator ){
                Edge * negativeCycle = backtrack( edge );
                if( negativeCycle == NULL ){
                  system->falsePositives++;
                }
                else {
                  system->mainLoopIterations = i;
                  return negativeCycle;
                }
              }
            }
            anyChange = true;
          }
          edge = edge->next;
        }
      }
    }
    if( !anyChange ){
      system->mainLoopIterations = i;
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

static Edge * backtrack(Edge * edge){
  Edge * input = edge;
  while( edge != NULL && edge->backtrackSeen == false ){
    edge->backtrackSeen = true;
    edge = edge->tail->L;
  }
  if( edge == NULL ){
    edge = input;
    while( edge != NULL && edge->backtrackSeen == true ){
      edge->backtrackSeen = false;
      edge = edge->tail->L;
    }
  }
  return edge;
}

static int lahiri(System * Gphi){
  System GphiPrime;
  onlySlacklessEdges(Gphi, &GphiPrime);
  stronglyConnectedComponents(&GphiPrime);
  int output = -1;
  for(int i = 0; i < Gphi->n && output < 0; i++){
    if( GphiPrime.graph[POSITIVE][i].sccNumber == GphiPrime.graph[NEGATIVE][i].sccNumber
        && ( Gphi->graph[POSITIVE][i].D - Gphi->graph[NEGATIVE][i].D ) % 2 != 0 ){
      output = i;
    }
  }
  freeSystem( &GphiPrime );
  if( output >= 0 ){
    return output;
  }
  setSystemForJohnson( Gphi );
  johnsonAllPairs( Gphi );
  return output;
}

/*
 * Doesn't copy Vertex elements L, D, rho, dfsColor, discoveryTime, finishingTime, sccNumber, h, 
 * dijkstraFinalized, or fibHeapNode, because this is unnecessary for our purposes
 */
static void onlySlacklessEdges(System * original, System * subgraph){
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

static void stronglyConnectedComponents(System * system){
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

static void dfsVisit(Vertex * vertex, int * time, int sccNumber){
  (*time)++;
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
static void transposeSystem(System * original, System * transpose){
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

static int vertexCompareFinishingTimes(const void * vertex1, const void * vertex2){
  return (*(Vertex **)vertex1)->finishingTime - (*(Vertex **)vertex2)->finishingTime;
}

/*
 * Lahiri's Model Generation algorithm is integrated into Johnson All-Pairs, such that this function leaves Gphi with an integral 
 * solution in Vertex.rho, without generating a graphical representation of the transitive and tight closure C*, thus avoiding the
 * O(n^2) space complexity term.
 * Leaves Gphi with modified edge weights, which is fine, because Gphi's edges are no longer used after this function call.
 */
static void johnsonAllPairs(System * Gphi){
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < Gphi->n; j++){
      Gphi->graph[i][j].h = Gphi->graph[i][j].D;
    }
  }
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < Gphi->n; j++){
      Edge * edge = Gphi->graph[i][j].first;
      while( edge != NULL ){
        edge->weight += edge->tail->h - edge->head->h;
        edge = edge->next;
      }
    }
  }
  for(int j = 0; j < Gphi->n; j++){
    int upperBound = INT_MAX;
    int lowerBound = INT_MIN;
    for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
      dijkstra( Gphi, &Gphi->graph[i][j] );
      for(int m = 0; m < j; m++){
        for(VertexSign k = POSITIVE; k <= NEGATIVE; k++){
          int weight = Gphi->graph[k][m].D;
          if( weight != INT_MAX ){
            weight += Gphi->graph[k][m].h - Gphi->graph[i][j].h;
            if( i == k ){
              switch( i ){
              case POSITIVE:
                ;
                int potentialNewLowerBound = -weight + Gphi->graph[k][m].rho;
                if( potentialNewLowerBound > lowerBound ){
                  lowerBound = potentialNewLowerBound;
                }
                break;
              case NEGATIVE:
                ;
                int potentialNewUpperBound = weight + Gphi->graph[k][m].rho;
                if( potentialNewUpperBound < upperBound ){
                  upperBound = potentialNewUpperBound;
                }
              }
            }
            else {
              switch( i ){
              case POSITIVE:
                ;
                int potentialNewLowerBound = -weight - Gphi->graph[k][m].rho;
                if( potentialNewLowerBound > lowerBound ){
                  lowerBound = potentialNewLowerBound;
                }
                break;
              case NEGATIVE:
                ;
                int potentialNewUpperBound = weight - Gphi->graph[k][m].rho;
                if( potentialNewUpperBound < upperBound ){
                  upperBound = potentialNewUpperBound;
                }
              }
            }
          }
        }
      }
      int weight = Gphi->graph[!i][j].D;
      if( weight != INT_MAX ){
        weight += Gphi->graph[!i][j].h - Gphi->graph[i][j].h;
        if( weight % 2 != 0 ){
          weight--;
        }
        switch( i ){
        case POSITIVE:
          ;
          int potentialNewLowerBound = -weight/2;
          if( potentialNewLowerBound > lowerBound ){
            lowerBound = potentialNewLowerBound;
          }
          break;
        case NEGATIVE:
          ;
          int potentialNewUpperBound = weight/2;
          if( potentialNewUpperBound < upperBound ){
            upperBound = potentialNewUpperBound;
          }
        }
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
      solution = (upperBound + lowerBound)/2;
    }
    Gphi->graph[POSITIVE][j].rho = solution;
    Gphi->graph[NEGATIVE][j].rho = solution;
  }
}

static void dijkstra(System * system, Vertex * source){
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < system->n; j++){
      system->graph[i][j].D = INT_MAX;
    }
  }
  source->D = 0;
  
  FibHeap pQueue;
  pQueue.min = NULL;
  pQueue.n = 0;

  fibHeapInsert( &pQueue, source );
  while( pQueue.n > 0 ){
    Vertex * current = fibHeapExtractMin( &pQueue );
    Edge * edge = current->first;
    while( edge != NULL ){
      if( edge->head->D > edge->tail->D + edge->weight ){
        edge->head->D = edge->tail->D + edge->weight;
        if( edge->head->fibHeapNode == NULL ){
          fibHeapInsert( &pQueue, edge->head );
        }
        else {
          fibHeapDecreaseKey( &pQueue, edge->head );
        }
      }
      edge = edge->next;
    }
  }
  
}

static void fibHeapInsert(FibHeap * fibHeap, Vertex * vertex){
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

static Vertex * fibHeapExtractMin(FibHeap * fibHeap){
  Vertex * output;
  if( fibHeap->min != NULL ){
    output = fibHeap->min->vertex;
    FibHeapNode * oldFHN = fibHeap->min;
    if( fibHeap->min->child != NULL ){
      FibHeapNode * fhn = fibHeap->min->child;
      while( fhn->parent != NULL ){
        fhn->parent = NULL;
        fhn = fhn->left;
      }
      if( fibHeap->min->right != fibHeap->min ){
        fibHeap->min->right->left = fibHeap->min->child;
        fibHeap->min->child->right->left = fibHeap->min->left;
        fibHeap->min->left->right = fibHeap->min->child->right;
        fibHeap->min->child->right = fibHeap->min->right;
      }
      fibHeap->min = fibHeap->min->child;
    }
    else if( fibHeap->min->right != fibHeap->min ){
      fibHeap->min->right->left = fibHeap->min->left;
      fibHeap->min->left->right = fibHeap->min->right;
      fibHeap->min = fibHeap->min->right;
    }
    else {
      fibHeap->min = NULL;
    }
    output->fibHeapNode = NULL;
    if( fibHeap->min != NULL ){
      fibHeapConsolidate( fibHeap );
    }
    free( oldFHN );
    fibHeap->n--;
  }
  else {
    output = NULL;
  }
  return output;
}

static void fibHeapConsolidate(FibHeap * fibHeap){
  double phi = ( 1.0 + sqrt( 5.0 ) ) / 2.0;
  int Alength = ((int)( log( fibHeap->n ) / log( phi ) )) + 1;
  FibHeapNode * A[ Alength ];
  for(int i = 0; i < Alength; i++){
    A[i] = NULL;
  }
  FibHeapNode * w = fibHeap->min;
  do {
    w->rootListTraverseSeen = false;
    w = w->left;
  } while( w != fibHeap->min );
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

static void fibHeapLink(FibHeap * fibHeap, FibHeapNode * y, FibHeapNode * x){
  y->left->right = y->right;
  y->right->left = y->left;
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
  y->parent = x;
  x->degree++;
  y->mark = false;
}

/*
 * Will silently fail if vertex->D is actually set higher than it had been.
 * This scenario can cause the heap to violate the min-heap property, making it useless. 
 */
static void fibHeapDecreaseKey(FibHeap * fibHeap, Vertex * vertex){
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

static void fibHeapCut(FibHeap * fibHeap, FibHeapNode * x, FibHeapNode * y){
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

static void fibHeapCascadingCut(FibHeap * fibHeap, FibHeapNode * y){
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

static void freeSystem(System * system){
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
