#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "utvpiInterpreter.h"
#include "halfint.h"

#define EDGE_TYPE_COUNT 4
typedef enum EdgeType {
  WHITE,
  BLACK,
  GRAY_FORWARD,
  GRAY_REVERSE,
} EdgeType;

#define BACKTRACKING_INDEX_COUNT 2
typedef enum BacktrackingIndex {
  NEG_ONE,
  POS_ONE,
} BacktrackingIndex;

#define INTEGER_TYPE_COUNT 2
typedef enum IntegerType {
  FINAL,
  TEMP,
} IntegerType;

typedef struct System System;
typedef struct Vertex Vertex;
typedef struct Edge Edge;
typedef struct EdgeRefList EdgeRefList;
typedef struct EdgeRefListNode EdgeRefListNode;
typedef struct IntegerTree IntegerTree;
typedef struct IntegerTreeVertex IntegerTreeVertex;

struct System {
  Vertex * graph;
  int vertexCount;
  int n;
  int C;
  Edge * allEdgeFirst;
  EdgeRefList * infeasibilityProof;
  IntegerTree * T;
};

struct Vertex {
  int index;
  Edge * L[EDGE_TYPE_COUNT];
  int D[EDGE_TYPE_COUNT];
  Edge * E[BACKTRACKING_INDEX_COUNT];
  half_int a;
  int Z[INTEGER_TYPE_COUNT];
  Edge * first[EDGE_TYPE_COUNT];
  int edgeCount[EDGE_TYPE_COUNT];
  IntegerTreeVertex * integerTreeVertex;
};

struct Edge {
  int weight;
  EdgeType type;
  Vertex * tail;
  Vertex * head;
  Edge * reverse;
  Edge * next;
  Edge * allNext;
  Edge * allPrev;
  bool inAllEdgeList;
};

struct EdgeRefList {
  EdgeRefListNode * first;
  EdgeRefListNode * last;
};

struct EdgeRefListNode {
  Edge * edge;
  EdgeRefListNode * next;
};

struct IntegerTree {
  IntegerTreeVertex * treeRoot;
  IntegerTreeVertex * queueNewest;
  IntegerTreeVertex * queueOldest;
  Edge * additionsFirst;
};

struct IntegerTreeVertex {
  IntegerTreeVertex * parent;
  IntegerTreeVertex * queueNewer;
  EdgeRefList * graphEdges;
  Vertex * graphVertex;
};

int main(int argc, char * argv[]);
static EdgeType reverseEdgeType(EdgeType input);
static void fputEdge(Edge * edge, FILE * output);
static void initializeSystem(void * object, int n, Parser * parser);
static void addConstraint(void * object, Constraint * constraint, Parser * parser);
static void addEdge(System * system, Constraint * constraint);
static void finishSystemCreation(System * system);
static int edgeCompare(const void * edge1, const void * edge2);
static void removeFromAllEdgeList(System * system, Edge * edge);
static bool relaxNetwork(System * system);
static void relaxEdge(Edge * e);
static bool backtrack(System * system, Vertex * x_i, EdgeType t, Edge * e);
static bool produceIntegerSolution(System * system);
static bool forcedRounding(System * system, Vertex * x_i);
static bool optionalRoundings(System * system);
static void checkDependencies(IntegerTree * T, Vertex * x_i, IntegerType integerType);
static bool checkAllConstraints(System * system, Vertex * toVertex, IntegerType integerType);
static void systemSubset(System * system);
static Edge * generateAbsoluteConstraint(System * system, Vertex * x_i, int weight, EdgeType type);
static IntegerTree * generateIntegerTree(System * system);
static Vertex * pollIntegerTreeQueue(IntegerTree * tree);
static void expandIntegerTree(IntegerTree * T, Vertex * active, Vertex * parent, Edge * edge0, Edge * edge1, Edge * edge2);
static void integerTreeBacktrack(EdgeRefList * list, IntegerTreeVertex * fromVertex, IntegerTreeVertex * toVertex, bool includeToVertex);
static void copyTreeEdgesToList(EdgeRefList * list, IntegerTreeVertex * itv);
static void freeIntegerTree(IntegerTree * tree);
static EdgeRefList * generateEdgeRefList();
static void edgeRefListAppend(EdgeRefList * erl, Edge * edge);
static void edgeRefListPrepend(EdgeRefList * erl, Edge * edge);
static void freeEdgeRefList(EdgeRefList * erl);
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
  finishSystemCreation(&system);
  #ifdef __HPC__
    struct timespec beforeLinear;
    clock_gettime(CLOCK_MONOTONIC_RAW, &beforeLinear);
  #else
    clock_t beforeLinear = clock();
  #endif
  int f;
  bool linearlyFeasible = relaxNetwork(&system);
  #ifdef __HPC__
    struct timespec beforeIntegral;
    clock_gettime(CLOCK_MONOTONIC_RAW, &beforeIntegral);
  #else
    clock_t beforeIntegral = clock();
  #endif
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
  #ifdef __HPC__
    struct timespec beforeCleanup;
    clock_gettime(CLOCK_MONOTONIC_RAW, &beforeCleanup);
  #else
    clock_t beforeCleanup = clock();
  #endif
  freeSystem(&system);
  fclose(output);
  #ifdef __HPC__
    struct timespec end;
    clock_gettime(CLOCK_MONOTONIC_RAW, &end);
  #else
    clock_t end = clock();
  #endif
  printf("%d,", f);
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

static EdgeType reverseEdgeType(EdgeType input){ 
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

static void fputEdge(Edge * edge, FILE * output){
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

static void initializeSystem(void * object, int n, Parser * parser){
  System * system = (System *) object;
  system->vertexCount = n + 1;
  system->graph = (Vertex *) malloc( sizeof(Vertex) * system->vertexCount );
  system->n = n;
  system->C = 0;
  system->allEdgeFirst = NULL;
  system->infeasibilityProof = NULL;
  system->T = NULL;
  system->graph[0].index = 0;
  system->graph[0].E[NEG_ONE] = NULL;
  system->graph[0].E[POS_ONE] = NULL;
  system->graph[0].a = 0;
  system->graph[0].Z[FINAL] = 0;
  system->graph[0].Z[TEMP] = 0;
  system->graph[0].integerTreeVertex = NULL;
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
    system->graph[i].integerTreeVertex = NULL;
    for(EdgeType j = WHITE; j <= GRAY_REVERSE; j++){
      system->graph[i].L[j] = NULL;
      system->graph[i].D[j] = INT_MAX;
      system->graph[i].first[j] = NULL;
      system->graph[i].edgeCount[j] = 0;
    }
   }
}

static void addConstraint(void * object, Constraint * constraint, Parser * parser){
  System * system = (System *) object;
  if( abs( constraint->weight ) > system->C ){
    system->C = abs( constraint->weight );
  }
  if( constraint->sign[1] == CONSTRAINT_NONE ){
    constraint->index[1] = 0;
    constraint->sign[1] = CONSTRAINT_PLUS;
    addEdge( system, constraint );
    constraint->sign[1] = CONSTRAINT_MINUS;
    addEdge( system, constraint );
    constraint->sign[1] = CONSTRAINT_NONE;
    //Calls to addEdge place the new edges at first[edgeType] for both vertices
    if( constraint->sign[0] == CONSTRAINT_PLUS ){
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

static void addEdge(System * system, Constraint * constraint){
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
    if( constraint->sign[0] == CONSTRAINT_PLUS ){
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
    if( constraint->sign[0] == CONSTRAINT_PLUS ){
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

static void finishSystemCreation(System * system){
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

static int edgeCompare(const void * edge1, const void * edge2){
  return (*(Edge **)edge1)->head->index - (*(Edge **)edge2)->head->index;
}

static void removeFromAllEdgeList(System * system, Edge * edge){
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

static bool relaxNetwork(System * system){
  //Lines 3-6 of algorithm implemented in finishSystemCreation().
  for(int r = 1; r <= 2 * system->n; r++){
    Edge * e = system->allEdgeFirst;
    while( e != NULL ){
      relaxEdge(e);
      e = e->allNext;
    }
  }
  Edge * e = system->allEdgeFirst;
  while( e != NULL ){
    switch( e->type ){
    case WHITE:
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
      break;
    case BLACK:
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
      break;
    case GRAY_FORWARD:
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
      break;
    case GRAY_REVERSE:
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
    }
    e = e->allNext;
  }
  return true;
}

static void relaxEdge(Edge * e){
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

static bool backtrack(System * system, Vertex * x_i, EdgeType t, Edge * e){
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
  edgeRefListAppend( system->infeasibilityProof, e_c );
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
    edgeRefListAppend( system->infeasibilityProof, x_c->E[a_c] );
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

static bool produceIntegerSolution(System * system){

  bool integrallyFeasible;
  system->infeasibilityProof = generateEdgeRefList();
  system->T = generateIntegerTree(system);

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

  Vertex * vertex = pollIntegerTreeQueue( system->T );
  while( vertex != NULL ){
    checkDependencies(system->T, vertex, FINAL);
    vertex = pollIntegerTreeQueue( system->T );
  }

  integrallyFeasible = checkAllConstraints(system, &system->graph[0], FINAL);
  if( !integrallyFeasible ){
    return false;
  }

  freeEdgeRefList( system->infeasibilityProof );
  freeIntegerTree( system->T );

  systemSubset( system );

  integrallyFeasible = optionalRoundings(system);
  
  return integrallyFeasible;
}

static bool forcedRounding(System * system, Vertex * x_i){
  bool forcedDown = false;
  Edge * whiteEdge = x_i->first[WHITE];
  Edge * grayReverseEdge = x_i->first[GRAY_REVERSE];
  while( whiteEdge != NULL && grayReverseEdge != NULL && !forcedDown ){
    if( grayReverseEdge->head->index == whiteEdge->head->index ){
      if( (intToHalfInt( whiteEdge->weight ) == x_i->a + whiteEdge->head->a) 
          && (intToHalfInt( grayReverseEdge->weight ) == x_i->a - whiteEdge->head->a) ){
      
        x_i->Z[FINAL] = halfIntToInt( halfIntFloor( x_i->a ) );
        Edge * newEdge = generateAbsoluteConstraint(system, x_i, x_i->Z[FINAL], WHITE);
        expandIntegerTree( system->T, x_i, &system->graph[0], whiteEdge, grayReverseEdge, newEdge );
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
        Edge * newEdge = generateAbsoluteConstraint(system, x_i, -x_i->Z[FINAL], BLACK);
        expandIntegerTree( system->T, x_i, &system->graph[0], blackEdge, grayForwardEdge, newEdge );
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
    copyTreeEdgesToList(system->infeasibilityProof, x_i->integerTreeVertex);
    return false;
  }
  return true;
}

static bool optionalRoundings(System * system){  
  for(int i = 1; i < system->vertexCount; i++){
    if( system->graph[i].Z[FINAL] == INT_MAX ){
      system->T = generateIntegerTree(system);
      system->infeasibilityProof = generateEdgeRefList();
      Edge * newEdge;
      Vertex * vertex;
      
      for(int j = 1; j < system->vertexCount; j++){
        system->graph[j].Z[TEMP] = INT_MAX;
      }
      
      system->graph[i].Z[TEMP] = halfIntToInt( halfIntFloor( system->graph[i].a ) );
      newEdge = generateAbsoluteConstraint(system, &system->graph[i], system->graph[i].Z[TEMP], WHITE);
      expandIntegerTree( system->T, &system->graph[i], &system->graph[0], newEdge, NULL, NULL);

      vertex = pollIntegerTreeQueue( system->T );
      while( vertex != NULL ){
        checkDependencies(system->T, vertex, TEMP);
        vertex = pollIntegerTreeQueue( system->T );
      }

      bool floorFeasible = checkAllConstraints(system, &system->graph[i], TEMP);
      
      if( floorFeasible ){
        for(int j = 1; j < system->vertexCount; j++){
          if( system->graph[j].Z[TEMP] != INT_MAX ){
            system->graph[j].Z[FINAL] = system->graph[j].Z[TEMP];
          }
        }  
      }
      else {
        
        freeIntegerTree( system->T );
        system->T = generateIntegerTree( system );
        
        for(int j = 1; j < system->vertexCount; j++){
          system->graph[j].Z[TEMP] = INT_MAX;
        }
        
        system->graph[i].Z[TEMP] = halfIntToInt( halfIntCeil( system->graph[i].a ) );
        newEdge = generateAbsoluteConstraint( system, &system->graph[i], -system->graph[i].Z[TEMP], BLACK );
        expandIntegerTree( system->T, &system->graph[i], &system->graph[0], newEdge, NULL, NULL);

        vertex = pollIntegerTreeQueue( system->T );
        while( vertex != NULL ){
          checkDependencies( system->T, vertex, TEMP );
          vertex = pollIntegerTreeQueue( system->T );
        }

        bool ceilFeasible = checkAllConstraints( system, &system->graph[i], TEMP );
        
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
      
      freeIntegerTree(system->T);
      freeEdgeRefList(system->infeasibilityProof);
    }
  }
  system->T = NULL;
  system->infeasibilityProof = NULL;
  return true;
}

//x_i - must never represent a Vertex where x_i->a is integral. The first for loop in produceIntegerSolution
//ensures this can never occur
static void checkDependencies(IntegerTree * T, Vertex * x_i, IntegerType integerType){
  if( x_i->Z[integerType] == halfIntToInt( halfIntFloor( x_i->a ) ) ){
    Edge * grayForwardEdge = x_i->first[GRAY_FORWARD];
    while( grayForwardEdge != NULL ){
      if( intToHalfInt( grayForwardEdge->weight ) == -x_i->a + grayForwardEdge->head->a 
          && grayForwardEdge->head->Z[integerType] == INT_MAX ){
        grayForwardEdge->head->Z[integerType] = halfIntToInt( halfIntFloor( grayForwardEdge->head->a ) );
        expandIntegerTree(T, grayForwardEdge->head, x_i, grayForwardEdge, NULL, NULL);
      }
      grayForwardEdge = grayForwardEdge->next;
    }
    Edge * blackEdge = x_i->first[BLACK];
    while( blackEdge != NULL ){
      if( intToHalfInt( blackEdge->weight ) == -x_i->a - blackEdge->head->a 
          && blackEdge->head->Z[integerType] == INT_MAX ){
        blackEdge->head->Z[integerType] = halfIntToInt( halfIntCeil( blackEdge->head->a ) );
        expandIntegerTree(T, blackEdge->head, x_i, blackEdge, NULL, NULL);
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
        expandIntegerTree(T, whiteEdge->head, x_i, whiteEdge, NULL, NULL);
      }
      whiteEdge = whiteEdge->next;
    }
    Edge * grayReverseEdge = x_i->first[GRAY_REVERSE];
    while( grayReverseEdge != NULL ){
      if( intToHalfInt( grayReverseEdge->weight ) == x_i->a - grayReverseEdge->head->a 
          && grayReverseEdge->head->Z[integerType] == INT_MAX ){
        grayReverseEdge->head->Z[integerType] = halfIntToInt( halfIntCeil( grayReverseEdge->head->a ) );
        expandIntegerTree(T, grayReverseEdge->head, x_i, grayReverseEdge, NULL, NULL);
      }
      grayReverseEdge = grayReverseEdge->next;
    }
  }
}

static bool checkAllConstraints(System * system, Vertex * toVertex, IntegerType integerType){
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
        edgeRefListPrepend(system->infeasibilityProof, edge);
        integerTreeBacktrack(system->infeasibilityProof, edge->tail->integerTreeVertex, toVertex->integerTreeVertex, false);
        integerTreeBacktrack(system->infeasibilityProof, edge->head->integerTreeVertex, toVertex->integerTreeVertex, true);
        return false;
      }
    }
    edge = edge->allNext;
  }
  return true;
}

static void systemSubset(System * system){
  for(int i = 0; i < system->vertexCount; i++){
    for(EdgeType j = WHITE; j <= GRAY_REVERSE; j++){
      Edge * prior = NULL;
      Edge * edge = system->graph[i].first[j];
      while( edge != NULL ){
        Edge * next = edge->next;
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
        edge = next;
      }
      if( prior == NULL ){
        system->graph[i].first[j] = NULL;
      }
      else {
        prior->next = NULL;
      }
    }
  }
}

//type : WHITE or GRAY_REVERSE for positive variable coefficient
//       BLACK or GRAY_FORWARD for negative variable coefficient
static Edge * generateAbsoluteConstraint(System * system, Vertex * x_i, int weight, EdgeType type){
  Edge * newEdge = (Edge *) malloc( sizeof(Edge) );
  newEdge->weight = weight;
  newEdge->type = type;
  newEdge->tail = x_i;
  newEdge->head = &system->graph[0];
  newEdge->reverse = NULL;
  newEdge->next = system->T->additionsFirst;
  system->T->additionsFirst = newEdge;
  newEdge->allNext = NULL;
  newEdge->allPrev = NULL;
  newEdge->inAllEdgeList = false;
  return newEdge;
}

static IntegerTree * generateIntegerTree(System * system){
  
  IntegerTreeVertex * x0Root = (IntegerTreeVertex *) malloc( sizeof(IntegerTreeVertex) );
  x0Root->parent = NULL;
  x0Root->queueNewer = NULL;
  x0Root->graphEdges = generateEdgeRefList();
  x0Root->graphVertex = &system->graph[0];
  system->graph[0].integerTreeVertex = x0Root;
  
  IntegerTree * T = (IntegerTree *) malloc( sizeof(IntegerTree) );
  T->treeRoot = x0Root;
  T->queueNewest = x0Root;
  T->queueOldest = NULL;
  T->additionsFirst = NULL;
  
  return T;
}

static Vertex * pollIntegerTreeQueue(IntegerTree * tree){
  Vertex * output;
  if( tree->queueOldest != NULL ){
    output = tree->queueOldest->graphVertex;
    tree->queueOldest = tree->queueOldest->queueNewer;
  }
  else {
    output = NULL;
  }
  return output;
}

static void expandIntegerTree(IntegerTree * T, Vertex * active, Vertex * parent, Edge * edge0, Edge * edge1, Edge * edge2){
  
  IntegerTreeVertex * itv = active->integerTreeVertex;
  if( itv == NULL ){
    itv = (IntegerTreeVertex *) malloc( sizeof(IntegerTreeVertex) );
    itv->parent = parent->integerTreeVertex;
    itv->queueNewer = NULL;
    T->queueNewest->queueNewer = itv;
    T->queueNewest = itv;
    if( T->queueOldest == NULL ){
      T->queueOldest = itv;
    }
    itv->graphEdges = generateEdgeRefList();
    itv->graphVertex = active;
    active->integerTreeVertex = itv;
  }
  
  if( edge0 != NULL ){
    edgeRefListAppend( itv->graphEdges, edge0 );
    if( edge1 != NULL ){
      edgeRefListAppend( itv->graphEdges, edge1 );
      if( edge2 != NULL ){
        edgeRefListAppend( itv->graphEdges, edge2 );
      }
    }
  }

}

static void integerTreeBacktrack(EdgeRefList * list, IntegerTreeVertex * fromVertex, IntegerTreeVertex * toVertex, bool includeToVertex){
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

static void copyTreeEdgesToList(EdgeRefList * list, IntegerTreeVertex * itv){
  EdgeRefListNode * treeEdge = itv->graphEdges->first;
  while( treeEdge != NULL ){
    edgeRefListPrepend( list, treeEdge->edge );
    treeEdge = treeEdge->next;
  }
}

static void freeIntegerTree(IntegerTree * tree){
  if( tree != NULL ){
    IntegerTreeVertex * itv = tree->treeRoot;
    while( itv != NULL ){
      freeEdgeRefList( itv->graphEdges );
      itv->graphVertex->integerTreeVertex = NULL;
      IntegerTreeVertex * oldITV = itv;
      itv = itv->queueNewer;
      free( oldITV );
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

static EdgeRefList * generateEdgeRefList(){
  EdgeRefList * newERL = (EdgeRefList *) malloc( sizeof(EdgeRefList) );
  newERL->first = NULL;
  newERL->last = NULL;
  return newERL;
}

static void edgeRefListAppend(EdgeRefList * erl, Edge * edge){
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

static void edgeRefListPrepend(EdgeRefList * erl, Edge * edge){
  EdgeRefListNode * newERLN = (EdgeRefListNode *) malloc( sizeof(EdgeRefListNode) );
  newERLN->edge = edge;
  newERLN->next = erl->first;
  erl->first = newERLN;
  if( erl->last == NULL ){
    erl->last = newERLN;
  }
}

static void freeEdgeRefList(EdgeRefList * erl){
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

static void freeSystem(System * system){
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
  freeIntegerTree( system->T );
  free( system->graph );
}
