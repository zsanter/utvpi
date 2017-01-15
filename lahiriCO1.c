/*
 * lahiriCO1.c
 * The Lahiri-Musuvathi integral UTVPI system solver, with first cycle-originator arrangement
 * An Efficient Decision Procedure for UTVPI Constraints - Lahiri, Musuvathi
 * Empirical results given in the VMCAI paper
 *
 * Cycle-originators are very reliable in this version, at least when it comes to the large, linearly infeasible (f0) input
 * systems. These systems consist of a large number of negative cost edges, leading to a situation where many edges could be part
 * of a number of different negative cost cycles. As edges are relaxed, cycle-originators are passed through a constantly-changing
 * predecessor structure. If they still somehow find their way back to where they started, chances are good that there will be a
 * negative cycle somewhere back through the predecessor structure, even if the current edge isn't within one. This cycle-
 * originator setup likely can not be guaranteed to work in any given situation.
 *
 * Call with [executable] [input file] {output file}
 * [input file] must be properly formatted to be read by utvpiInterpreter.h
 * {output file} will contain a linear solution, if one exists, followed by an integral solution, if one exists. If the system is
 *   not linearly feasible, a proof of linear infeasibility - a negative cost cycle - will be output. If the system is linearly
 *   feasible, but not integrally feasible, a proof of integral infeasibility will be output. If {output file} is not specified,
 *   output will be to stdout.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include "delay.h"

/*
 * This algorithm represents each variable with two vertices, one representing the positive occurrence of the variable, and the
 * other representing its negative occurrence. VertexSign specifies the sign of each Vertex.
 */
#define VERTEX_SIGN_COUNT 2
typedef enum VertexSign {
  POSITIVE,
  NEGATIVE,
} VertexSign;

/*
 * DFScolor is used to specify the color of a Vertex during the depth-first search process used when detecting strongly-connected
 * components.
 */
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

/*
 * The System struct contains the graph representation used by the system solver, along with other information about the system.
 *
 * graph - pointers to POSITIVE and NEGATIVE arrays of Vertex structs, defining the structure of the graph representation
 * n - the number of variables represented within the system, equal to the number of Vertex structs stored within each graph array
 * falsePositives - number of false positives thrown by the cycle-originator negative cycle detection mechanism
 * mainLoopIterations - number of iterations of the main relaxation loop within bellmanFord()
 * negativeCycleEdgeCount - number of edges within a detected negative cost cycle
 */
struct System {
  Vertex * graph[VERTEX_SIGN_COUNT];
  int n;
  #ifdef __HPC__
    struct timespec beforeLinear;
  #else
    clock_t beforeLinear;
  #endif
  int falsePositives;
  int mainLoopIterations;
  int negativeCycleEdgeCount;
};

/*
 * The Vertex struct contains all information about one signed occurrence of a variable, represented by a vertex within the graph.
 *
 * index - index of the variable
 * sign - the sign of the variable occurrence this Vertex represents
 * L - pointer to the predecessor Edge between this Vertex and the Vertex's predecessor Vertex
 * D - distance label
 * cycleOriginator - pointer to an edge serving as a cycle-originator. Cycle-originators are used to detect negative cost cycles
 *   before the total number of relaxation loop iterations specified by the algorithm run
 * first - pointer to the first Edge whose tail Vertex is this Vertex. The remainder of such edges are connected together in a
 *   singly-linked list
 * dfsColor - the color of a Vertex during the depth-first search process used when detecting strongly-connected components
 * finishingTime - depth-first search Vertex finishing time
 * sccNumber - vertices within each strongly-connected component are assigned the same sccNumber as each other
 * h - Johnson All-Pairs reweighting value
 * rho - integral-solution variable value
 * fibHeapNode - pointer to the FibHeapNode corresponding to this Vertex within the FibHeap priority queue data structure
 *
 * Sets of data members that do not need to exist at the same time as each other are organized into anonymous structs within an
 * anonymous union. This necessitates initializing each set of data members immediately before it is used.
 */
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

/*
 * The Edge struct contains all information about a specific constraint, represented by an edge within the graph.
 *
 * weight - the weight of the edge, corresponding to the defining constant of the constraint the edge represents
 * head - pointer to the Edge's head Vertex
 * tail - pointer to the Edge's tail Vertex
 * next - pointer to the next Edge in the singly-linked list of Edges with the same tail Vertex as one another
 * backtrackSeen - boolean flag indicating whether or not this Edge has been encountered by the backtrack() function as it
 *   attempts to detect a negative cost cycle
 */
struct Edge {
  int weight;
  Vertex * head;
  Vertex * tail;
  Edge * next;
  bool backtrackSeen;
};

/*
 * FibHeap is an implementation of the Fibonacci Heap priority queue structure, as defined in Introduction To Algorithms, Third
 * Edition
 *
 * min - pointer to the FibHeapNode representing the Vertex with the lowest distance label D
 * n - the number of FibHeapNodes within the overall FibHeap data structure
 */
struct FibHeap {
  FibHeapNode * min;
  int n;
};

/*
 * FibHeapNode is a single node within the FibHeap, as defined in Introduction to Algorithms, Third Edition
 *
 * parent - pointer to this FibHeapNode's parent FibHeapNode
 * left - pointer to the FibHeapNode to the left of this FibHeapNode, in the circular, doubly-linked list of sibling FibHeapNodes
 * right - pointer to the FibHeapNode to the right of this FibHeapNode, in the circular, doubly-linked list of sibling
 *   FibHeapNodes
 * child - pointer to one of this FibHeapNode's child FibHeapNodes
 * vertex - pointer to the Vertex struct that this FibHeapNode represents within the FibHeap
 * degree - the number of children of this FibHeapNode
 * mark - used to indicate if a FibHeapNode has lost two of its child FibHeapNodes in a row
 * rootListTraverseSeen - boolean flag indicating whether or not this FibHeapNode has been encountered by fibHeapConsolidate() as
 *   it traverses the root list
 */
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
static void edgeToConstraint(Edge * edge, Constraint * constraint);
static void initializeSystem(void * object, int n, Parser * parser);
static void setSystemForJohnson(System * system);
static void addConstraint(void * object, Constraint * constraint, Parser * parser);
static Edge * bellmanFord(System * system);
static Edge * backtrack(Edge * edge);
static ConstraintRefList * lahiri(System * Gphi);
static ConstraintRefList * generateProof(System * Gphi, System * GphiPrime, int infeasibleVertexIndex);
static void onlySlacklessEdges(System * original, System * subgraph);
static void stronglyConnectedComponents(System * system);
static bool dfsVisit(Vertex * vertex, int * time, int sccNumber, Vertex * destination);
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

/*
 * main()
 * - handles file input and output
 * - calls utvpiInterpreter's parseFile() function, which calls initializeSystem() and addConstraint() to build the graph
 *   representation corresponding to the constraint system in the input file
 * - calls bellmanFord()
 * - calls lahiri()
 * - prints profiling information formatted for csv-file input to stdout
 */
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
  bool parseSuccessful = parseFileDelayedSysGen(input, &system, initializeSystem, addConstraint);
  fclose(input);
  if( !parseSuccessful ){
    fprintf( stderr, "Parsing of \"%s\" failed.\n", argv[1] );
    exit(1);
  }
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
    ConstraintRefList * infeasibilityProof = lahiri(&system);
    if( infeasibilityProof != NULL ){
      f = 1;
      fputs("\nSystem is not integrally feasible.\nProof:\n", output);
      Constraint * constraint = constraintRefListNext( infeasibilityProof );
      while( constraint != NULL ){
        fputConstraint(constraint, output);
        constraint = constraintRefListNext( infeasibilityProof );;
      }
      freeConstraintRefList( infeasibilityProof );
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
    diff(&start, &system.beforeLinear, &setup);
    printf("%d.%09d,", (int)setup.tv_sec, (int)setup.tv_nsec);
    struct timespec linear;
    diff(&system.beforeLinear, &beforeIntegral, &linear);
    printf("%d.%09d,", (int)linear.tv_sec, (int)linear.tv_nsec);
    struct timespec integral;
    diff(&beforeIntegral, &beforeCleanup, &integral);
    printf("%d.%09d,", (int)integral.tv_sec, (int)integral.tv_nsec);
    struct timespec cleanup;
    diff(&beforeCleanup, &end, &cleanup);
    printf("%d.%09d,", (int)cleanup.tv_sec, (int)cleanup.tv_nsec);
    struct timespec total;
    diff(&system.beforeLinear, &end, &total);
    printf("%d.%09d,", (int)total.tv_sec, (int)total.tv_nsec);
  #else
    printf("%f,", ((double)(system.beforeLinear - start))/CLOCKS_PER_SEC);
    printf("%f,", ((double)(beforeIntegral - system.beforeLinear))/CLOCKS_PER_SEC);
    printf("%f,", ((double)(beforeCleanup - beforeIntegral))/CLOCKS_PER_SEC);
    printf("%f,", ((double)(end - beforeCleanup))/CLOCKS_PER_SEC);
    printf("%f,", ((double)(end - system.beforeLinear))/CLOCKS_PER_SEC);
  #endif
  return 0;
}

#ifdef __HPC__
  /*
   * diff() takes the difference of two struct timespecs
   * Copied from https://www.guyrutenberg.com/2007/09/22/profiling-code-using-clock_gettime/ and modified.
   *
   * start - pointer to a struct timespec from the beginning of a time period
   * end - pointer to a struct timespec from the end of a time period
   * difference - pointer to an already-declared struct timespec to be filled with the time difference between start and end
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

/*
 * fputEdge() prints the constraint equation corresponding to edge to output
 *
 * edge - pointer to an Edge struct to convert to a constraint equation
 * output - FILE pointer to print the constraint equation to
 */
static void fputEdge(Edge * edge, FILE * output){
  Constraint constraint;
  edgeToConstraint(edge, &constraint);
  fputConstraint(&constraint, output);
}

/*
 * edgeToConstraint() fills in the referenced Constraint struct with the constraint represented by the referenced Edge.
 *
 * edge - pointer to an Edge to be interpreted as a constraint
 * constraint - pointer to a Constraint to be filled with the constraint represented by edge
 */
static void edgeToConstraint(Edge * edge, Constraint * constraint){
  if( edge->tail->index == edge->head->index ){
    switch( edge->tail->sign ){
    case POSITIVE:
      constraint->sign[0] = CONSTRAINT_MINUS;
      break;
    case NEGATIVE:
      constraint->sign[0] = CONSTRAINT_PLUS;
    }
    constraint->index[0] = edge->tail->index;
    constraint->sign[1] = CONSTRAINT_NONE;
    constraint->index[1] = 0;
    constraint->weight = (edge->weight)/2;
  }
  else{
    if( edge->tail->sign == edge->head->sign ){
      switch( edge->tail->sign ){
      case POSITIVE:
        constraint->sign[0] = CONSTRAINT_MINUS;
        constraint->sign[1] = CONSTRAINT_PLUS;
        break;
      case NEGATIVE:
        constraint->sign[0] = CONSTRAINT_PLUS;
        constraint->sign[1] = CONSTRAINT_MINUS;
      }
    }
    else{
      switch( edge->tail->sign ){
      case POSITIVE:
        constraint->sign[0] = CONSTRAINT_MINUS;
        constraint->sign[1] = CONSTRAINT_MINUS;
        break;
      case NEGATIVE:
        constraint->sign[0] = CONSTRAINT_PLUS;
        constraint->sign[1] = CONSTRAINT_PLUS;
      }
    }
    constraint->index[0] = edge->tail->index;
    constraint->index[1] = edge->head->index;
    constraint->weight = edge->weight;
  }
}

/*
 * initializeSystem() initializes an already-declared System struct, given the number of variables that the system must represent
 *
 * object - a void pointer pointing to an already-declared System struct
 * n - the number of variables that the graph representation held by the System struct must represent
 * parser - pointer to the Parser struct that utvpiInterpreter uses during the input file parsing process, so that parseError()
 *   can be called, if need be
 */
static void initializeSystem(void * object, int n, Parser * parser){
  System * system = (System *) object;
  #ifdef __HPC__
    clock_gettime(CLOCK_MONOTONIC_RAW, &system->beforeLinear);
  #else
    system->beforeLinear = clock();
  #endif
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

/*
 * setSystemForJohnson() initializes Vertex struct anonymous union data members only used during the process of generating an
 * integral solution, for all Vertices.
 *
 * system - pointer to the overall System struct containing the graph representation
 */
static void setSystemForJohnson(System * system){
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < system->n; j++){
      system->graph[i][j].h = 0;
      system->graph[i][j].fibHeapNode = NULL;
      system->graph[i][j].rho = INT_MAX;
    }
  }
}

/*
 * addConstraint() adds a constraint to the graph representation held by the System struct.
 *
 * object - a void pointer pointing to an already-initialized System struct
 * constraint - pointer to a Constraint struct describing a constraint
 * parser - pointer to the Parser struct that utvpiInterpreter uses during the input file parsing process, so that parseError()
 *   can be called, if need be
 */
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

/*
 * bellmanFord() implements BELLMAN-FORD() and RELAX(), as defined in Introduction to Algorithms, Third Edition, with the
 * modification of adding cycle-originators that are initialized during the first round of relaxations, then passed through the
 * predecessor structure as distance labels continue to decrease. For linearly infeasible systems, a pointer to one Edge within
 * the detected negative cost cycle is returned. For linearly feasible systems, NULL is returned.
 *
 * system - pointer to the overall System struct containing the graph representation
 */
static Edge * bellmanFord(System * system){
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < system->n; j++){
      Edge * edge = system->graph[i][j].first;
      while(edge != NULL){
        if( edge->head->D > edge->tail->D + edge->weight ){
          edge->head->D = edge->tail->D + edge->weight;
          edge->head->L = edge;
          edge->head->cycleOriginator = edge;
        }
        edge = edge->next;
      }
    }
  }
  bool anyChange = true;
  for(int i = 1; i <= (2 * system->n - 2 /*1*/) && anyChange; i++){
    anyChange = false;
    for(VertexSign j = POSITIVE; j <= NEGATIVE; j++){
      for(int k = 0; k < system->n; k++){
        Edge * edge = system->graph[j][k].first;
        while(edge != NULL){
          if( edge->head->D > edge->tail->D + edge->weight ){
            edge->head->D = edge->tail->D + edge->weight;
            edge->head->L = edge;
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

/*
 * backtrack() implements a single backtrack through the predecessor structure, as defined in Network Flows - Ahuja, Magnanti,
 * Orlin. If a negative cost cycle is detected, returns a pointer to an Edge within that cycle. Otherwise, cleans up after itself
 * and returns NULL.
 *
 * edge - pointer to an Edge to backtrack from
 */
static Edge * backtrack(Edge * edge){
  Edge * input = edge;
  while( edge != NULL && edge->backtrackSeen == false ){
    edge->backtrackSeen = true;
    edge = edge->tail->L;
  }
  if( edge == NULL ){
    edge = input;
    while( edge != NULL ){
      edge->backtrackSeen = false;
      edge = edge->tail->L;
    }
  }
  return edge;
}

/*
 * lahiri()
 * - implements the NCD-UTVPI() algorithm from An Efficient Decision Procedure for UTVPI Constraints - Lahiri, Musuvathi
 * - calls johnsonAllPairs(), into which the Model Generation algorithm from this paper was integrated, in order to significantly
 *   reduce the space complexity of this implementation
 * - if the system is found to be integrally infeasible, returns a pointer to a ConstraintRefList storing a proof of this
 * - otherwise, returns NULL
 *
 * Gphi - pointer to the overall System struct containing the graph representation
 */
static ConstraintRefList * lahiri(System * Gphi){
  System GphiPrime;
  onlySlacklessEdges(Gphi, &GphiPrime);
  stronglyConnectedComponents(&GphiPrime);
  int infeasibleVertexIndex = -1;
  for(int i = 0; i < Gphi->n && infeasibleVertexIndex < 0; i++){
    if( GphiPrime.graph[POSITIVE][i].sccNumber == GphiPrime.graph[NEGATIVE][i].sccNumber
        && ( Gphi->graph[POSITIVE][i].D - Gphi->graph[NEGATIVE][i].D ) % 2 != 0 ){
      infeasibleVertexIndex = i;
    }
  }
  if( infeasibleVertexIndex >= 0 ){
    ConstraintRefList * infeasibilityProof = generateProof(Gphi, &GphiPrime, infeasibleVertexIndex);
    freeSystem( &GphiPrime );
    return infeasibilityProof;
  }
  freeSystem( &GphiPrime );
  setSystemForJohnson( Gphi );
  johnsonAllPairs( Gphi );
  return NULL;
}

/*
 * generateProof() generates a proof of integral infeasibility, as defined in An Efficient Decision Procedure for UTVPI
 * Constraints - Lahiri, Musuvathi. Returns a pointer to a ConstraintRefList storing this proof.
 *
 * Gphi - pointer to the System struct containing the original graph representation
 * GphiPrime - pointer to the system struct containing only the slackless edges from Gphi
 * infeasibleVertexIndex - integer representing the index of the variable causing the contradiction
 */
static ConstraintRefList * generateProof(System * Gphi, System * GphiPrime, int infeasibleVertexIndex){
  ConstraintRefList * proof = generateConstraintRefList();
  int k = Gphi->graph[NEGATIVE][infeasibleVertexIndex].D - Gphi->graph[POSITIVE][infeasibleVertexIndex].D;

  Constraint * constraint;
  Edge * edge;

  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < GphiPrime->n; j++){
      GphiPrime->graph[i][j].dfsColor = WHITE;
    }
  }
  int time = 0;
  GphiPrime->graph[NEGATIVE][infeasibleVertexIndex].L = NULL;
  dfsVisit( &GphiPrime->graph[NEGATIVE][infeasibleVertexIndex], &time, 0, &GphiPrime->graph[POSITIVE][infeasibleVertexIndex] );
  edge = GphiPrime->graph[POSITIVE][infeasibleVertexIndex].L;
  while( edge != NULL ){
    constraint = (Constraint *) malloc( sizeof(Constraint) );
    edgeToConstraint(edge, constraint);
    constraintRefListPrepend(proof, constraint);
    edge = edge->tail->L;
  }

  constraint = (Constraint *) malloc( sizeof(Constraint) );
  constraint->sign[0] = CONSTRAINT_PLUS;
  constraint->index[0] = infeasibleVertexIndex + 1;
  constraint->sign[1] = CONSTRAINT_NONE;
  constraint->index[1] = 0;
  constraint->weight = (-k-1)/2;
  constraintRefListPrepend(proof, constraint);

  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < GphiPrime->n; j++){
      GphiPrime->graph[i][j].dfsColor = WHITE;
    }
  }
  time = 0;
  GphiPrime->graph[POSITIVE][infeasibleVertexIndex].L = NULL;
  dfsVisit( &GphiPrime->graph[POSITIVE][infeasibleVertexIndex], &time, 0, &GphiPrime->graph[NEGATIVE][infeasibleVertexIndex] );
  edge = GphiPrime->graph[NEGATIVE][infeasibleVertexIndex].L;
  while( edge != NULL ){
    constraint = (Constraint *) malloc( sizeof(Constraint) );
    edgeToConstraint(edge, constraint);
    constraintRefListPrepend(proof, constraint);
    edge = edge->tail->L;
  }

  constraint = (Constraint *) malloc( sizeof(Constraint) );
  constraint->sign[0] = CONSTRAINT_MINUS;
  constraint->index[0] = infeasibleVertexIndex + 1;
  constraint->sign[1] = CONSTRAINT_NONE;
  constraint->index[1] = 0;
  constraint->weight = (k-1)/2;
  constraintRefListPrepend(proof, constraint);

  return proof;
}

/*
 * onlySlacklessEdges() copies the subgraph of original consisting only of edges (u,v) where D(v) - D(u) = w(u,v) into subgraph.
 * Does not copy Vertex data members.
 *
 * original - pointer to the System struct containing the original graph representation
 * subgraph - pointer to a declared yet uninitialized System struct, to be filled with the subgraph of original consisting only of
 *   slackless edges
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

/*
 * stronglyConnectedComponents() finds strongly-connected components within system, and gives all Vertices within each strongly-
 * connected component the same Vertex.sccNumber value, which will differ from that within all Vertices not within that strongly-
 * connected component. Implements STRONGLY-CONNECTED-COMPONENTS() and DFS(), as defined in Introduction to Algorithms, Third
 * Edition.
 *
 * system - pointer to the System struct in which strongly-connected components are to be found
 */
static void stronglyConnectedComponents(System * system){
  int time = 0;
  int sccNumber = 0;
  for(VertexSign i = POSITIVE; i <= NEGATIVE; i++){
    for(int j = 0; j < system->n; j++){
      if( system->graph[i][j].dfsColor == WHITE ){
        dfsVisit( &system->graph[i][j], &time, sccNumber, NULL );
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
      dfsVisit( vertexSortArray[i], &time, sccNumber, NULL );
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

/*
 * dfsVisit() implements DFS-VISIT(), as defined in Introduction to Algorithms, Third Edition. Also sets Vertex.sccNumber for all
 * Vertices within a set of recurrences with the input sccNumber. This value is not useful unless this function is called on a
 * transposed System being examined in order of decreasing Vertex.finishingTime, as found when searching the original System, and
 * the input sccNumber is set to be different for each set of recurrences. Additionally, the set of recurrences of this function
 * will stop and return as soon as the destination Vertex is found. If this destination Vertex is found, the function returns
 * true. Otherwise, returns false.
 *
 * vertex - pointer to the Vertex to visit
 * time - pointer to an integer value for time, a value defined as global in Introduction to Algorithms
 * sccNumber - an integer value to label each Vertex within each strongly-connected component with
 * destination - pointer to a Vertex to be found. Set NULL if all Vertices must be found.
 */
static bool dfsVisit(Vertex * vertex, int * time, int sccNumber, Vertex * destination){
  if( vertex == destination ){
    return true;
  }
  (*time)++;
  vertex->dfsColor = GRAY;
  vertex->sccNumber = sccNumber;
  Edge * edge = vertex->first;
  while( edge != NULL ){
    if( edge->head->dfsColor == WHITE ){
      edge->head->L = edge;
      bool destinationFound = dfsVisit( edge->head, time, sccNumber, destination );
      if( destinationFound ){
        return true;
      }
    }
    edge = edge->next;
  }
  vertex->dfsColor = BLACK;
  (*time)++;
  vertex->finishingTime = *time;
  return false;
}

/*
 * transposeSystem() fills transpose with the transpose of original. Does not copy Vertex data members.
 *
 * original - pointer to the System struct containing the original graph representation
 * transpose - pointer to a declared but uninitialized System struct, to be filled with the transpose of original
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

/*
 * vertexCompareFinishingTimes() is the comparison function used by qsort() in stronglyConnectedComponents().
 */
static int vertexCompareFinishingTimes(const void * vertex1, const void * vertex2){
  return (*(Vertex **)vertex1)->finishingTime - (*(Vertex **)vertex2)->finishingTime;
}

/*
 * johnsonAllPairs() is an implementation of JOHNSON(), as given in Introduction to Algorithms, Third Edition, with the Model
 * Generation algorithm from An Efficient Decision Procedure for UTVPI Constraints - Lahiri, Musuvathi integrated directly into
 * it. Rather than generating and representing the transitive and tight closure C* as a graph stored in memory, C is processed by
 * johnsonAllPairs() in the same order as C* is processed by the Model Generation algorithm. Each time an edge within C* is
 * discovered, it is immediately used to modify bounds, then forgotten. This reduces the asymptotic space complexity of this
 * implementation from O(n^2) to O(m*n), reducing the practical memory requirements significantly.
 *
 * The call to BELLMAN-FORD() within JOHNSON() effectively occurs before this function is called, and is handled by bellmanFord().
 * This function leaves Gphi with an integral solution in Vertex.rho; rho for both POSITIVE and NEGATIVE Vertices are set to the
 * integral value for the variable whose positive or negative occurrence is represented by that Vertex. The function leaves Gphi
 * with modified edge weights; Gphi's edges are no longer used after this function call.
 *
 * Gphi - pointer to the overall System struct containing the graph representation of C
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

/*
 * dijkstra() implements DIJKSTRA(), as given in Introduction to Algorithms, Third Edition. This implementation differs from the
 * algorithm presented in the book in the fact that Vertices are only added to the priority queue when their distance labels are
 * first updated, as opposed to all Vertices being added to the priority queue before any distance labels are updated. Any vertex
 * unreachable from source will be left with a D distance label of INT_MAX.
 *
 * system - pointer to the System struct storing the graph representation to be processed
 * source - pointer to the source Vertex
 */
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

/*
 * fibHeapInsert() implements FIB-HEAP-INSERT(), as given in Introduction to Algorithms, Third Edition. No attempt is made to
 * detect if an input Vertex is already present within the FibHeap. This function must only be called on a given Vertex once
 * during the life of the FibHeap. This function sets Vertex.fibHeapNode for the input Vertex. This value should be initialized
 * to NULL before the FibHeap is used. As such, it can be determined whether or not a given Vertex is already present within the
 * FibHeap by testing whether or not Vertex.fibHeapNode != NULL.
 *
 * fibHeap - pointer to the FibHeap to which a Vertex is to be added
 * vertex - pointer to the Vertex to add to the FibHeap
 */
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

/*
 * fibHeapExtractMin() implements FIB-HEAP-EXTRACT-MIN(), as given in Introduction to Algorithms, Third Edition. Returns a pointer
 * to the Vertex within the FibHeap with the lowest distnace label D, removing this Vertex from the FibHeap, and setting its
 * Vertex.fibHeapNode member to NULL. If the FibHeap is empty, returns NULL.
 *
 * fibHeap - pointer to the FibHeap from which a Vertex with minimal distance label D is to be removed
 */
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

/*
 * fibHeapConsolidate() implements CONSOLIDATE(), as given in Introduction To Algorithms, Third Edition. Not to be called outside
 * of FibHeap implementation.
 *
 * fibHeap
 */
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

/*
 * fibHeapLink() implements FIB-HEAP-LINK(), as given in Introduction to Algorithms, Third Edition. Not to be called outside of
 * FibHeap implementation.
 *
 * fibHeap
 * y
 * x
 */
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
 * fibHeapDecreaseKey() implements FIB-HEAP-DECREASE-KEY(), as given in Introduction To Algorithms, Third Edition. The function
 * requires that the Vertex whose distance label D is decreased already exists within the FibHeap and that the distance label D is
 * decreased prior to this function being called. The function makes no attempt to check either. It can be determined whether or
 * not a given Vertex is already present within the FibHeap by testing whether or not Vertex.fibHeapNode != NULL. If the Vertex
 * does not already exist within the FibHeap, will attempt to access NULL. If the distance label D is actually increased prior to
 * calling this function, the FibHeap may be caused to violate the min-heap property, making it no longer useful as a priority
 * queue.
 *
 * fibHeap - pointer to the FibHeap that must be modified to correctly place vertex
 * vertex - pointer to the Vertex whose distance label has decreased and whose placement within the FibHeap therefore must change
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

/*
 * fibHeapCut() implements CUT(), as defined in Introduction to Algorithms, Third Edition. Not to be called outside of FibHeap
 * implementation.
 *
 * fibHeap
 * x
 * y
 */
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

/*
 * fibHeapCascadingCut() implements CASCADING-CUT(), as defined in Introduction to Algorithms, Third Edition. Not to be called
 * outside of FibHeap implementation.
 *
 * fibHeap
 * y
 */
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

/*
 * freeSystem() frees the graph representation stored within a System struct
 *
 * system - pointer to a System whose graph representation is to be freed
 */
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
