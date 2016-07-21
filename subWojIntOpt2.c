/*
 * subWojIntOpt2.c
 * The Subramani-Wojciechowski integral UTVPI system solver, with proof-backed optimization
 * Not yet referenced in a paper.
 *
 * Call with [executable] [input file] {output file}
 * [input file] must be properly formatted to be read by utvpiInterpreter.h
 * {output file} will contain a linear solution, if one exists, followed by an integral solution, if one exists. If the system is
 *   not linearly feasible, a proof of linear infeasibility - a negative cost cycle - will be output. If the system is linearly 
 *   feasible, but not integrally feasible, a proof of integral infeasibility will be output. If {output file} is not specified,
 *   output will be to stdout.
 */

#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "utvpiInterpreter.h"
#include "halfint.h"

/*
 * EdgeType is used to specify the color of an edge. The algorithm defines white, black, and gray edges. White and black edges are
 * non-directional, and gray edges are directional. All edges can be traversed in either direction, so each edge is represented by
 * a pair of Edge structs, with alternating head and tail pointers. 
 *
 * GRAY_FORWARD corresponds to an Edge struct representing a gray edge, where the tail pointer points to the vertex on the black 
 * side of the black and white box shown in the algorithm, and the head pointer points to the vertex on the white side of this 
 * box. GRAY_REVERSE corresponds to the opposite case. "Gray-Out" in the algorithm corresponds to GRAY_FORWARD, and "Gray-In"
 * corresponds to GRAY_REVERSE.
 *
 * For each white edge that enters the system, two WHITE Edge structs are created. For each black edge that enters the system, two 
 * BLACK Edge structs are created. For each gray edge that enters the system, one GRAY_FORWARD and one GRAY_REVERSE Edge struct
 * are created.
 *
 * These enums are also used to index into arrays.
 */
#define EDGE_TYPE_COUNT 4
typedef enum EdgeType {
  WHITE,
  BLACK,
  GRAY_FORWARD,
  GRAY_REVERSE,
} EdgeType;

/*
 * In the Backtrack algorithm, -1 and 1 are used as indeces. The members of the BacktrackingIndex enum take the place of these 
 * numbers for this purpose.
 */
#define BACKTRACKING_INDEX_COUNT 2
typedef enum BacktrackingIndex {
  NEG_ONE,
  POS_ONE,
} BacktrackingIndex;

/*
 * The integral algorithms make use of separate Z and Z^T arrays. Here, these values are stored in two-element Z arrays, contained
 * within each Vertex struct. Z[FINAL] refers to a variable's Z value, and Z[TEMP] refers to a variable's Z^T value.
 */
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

/*
 * The System struct contains references to all dynamically-allocated data structures used within the system solver, along with
 * other information about the system.
 *
 * graph - pointer to an array of Vertex structs, defining the structure of the graph representation
 * vertexCount - the number of Vertex structs within graph
 * n - the number of variables represented within the system - always one less than vertexCount, but used for the sake of clarity
 * C - the greatest absolute value of any edge weight which entered the system as part of a constraint
 * allEdgeFirst - pointer to the first Edge in a doubly-linked list where only one Edge struct representing each edge is present
 * infeasibilityProof - pointer to an EdgeRefList storing the current proof of infeasibility
 * T - pointer to the tree structure used in the integral algorithm
 * falsePositives - number of false positives thrown by the cycle-originator negative cycle detection mechanism
 * mainLoopIterations - number of iterations of the main relaxation loop within relaxNetwork()
 * negativeCycleEdgeCount - number of edges within a detected negative weight cycle
 */
struct System {
  Vertex * graph;
  int vertexCount;
  int n;
  int C;
  Edge * allEdgeFirst;
  EdgeRefList * infeasibilityProof;
  IntegerTree * T;
  int falsePositives;
  int mainLoopIterations;
  int negativeCycleEdgeCount;
};

/*
 * The Vertex struct contains all information about a specific variable, represented by a vertex within the graph.
 *
 * index - index of the variable
 * L - pointer to an edge between this vertex and the vertex's predecessor vertex, for each EdgeType-color path between the source
 *   vertex and this vertex
 * D - distance label for each EdgeType-color path between the source vertex and this vertex
 * cycleOriginator - pointer to an edge serving as a cycle-originator. Cycle-originators are used to detect negative cost cycles 
 *   before the total number of relaxation loop iterations specified by the algorithm run
 * E - pointer to the last Edge traversed on the path from x_i to x, in backtrack()
 * a - half-integral linear-solution variable value
 * Z - TEMP and FINAL integral-solution variable values
 * first - pointer to the first Edge of each EdgeType whose tail Vertex is this Vertex. The remainder of such edges are connected
 *   together in a singly-linked list
 * edgeCount - the number of each EdgeType edge whose tail Vertex is this Vertex
 * integerTreeVertex - pointer to an IntegerTreeVertex representing the same variable as this Vertex, within the IntegerTree data 
 *   structure
 */
struct Vertex {
  int index;
  Edge * L[EDGE_TYPE_COUNT];
  int D[EDGE_TYPE_COUNT];
  Edge * cycleOriginator[EDGE_TYPE_COUNT];
  Edge * E[BACKTRACKING_INDEX_COUNT];
  half_int a;
  int Z[INTEGER_TYPE_COUNT];
  Edge * first[EDGE_TYPE_COUNT];
  int edgeCount[EDGE_TYPE_COUNT];
  IntegerTreeVertex * integerTreeVertex;
};

/*
 * The Edge struct contains all information about a specific constraint, represented by an edge within the graph.
 *
 * weight - the weight of the edge, corresponding to the defining constant of the constraint the edge represents
 * type - the color type of the edge
 * tail - pointer to the Edge's tail Vertex
 * head - pointer to the Edge's head Vertex
 * reverse - pointer to the other Edge representing the same constraint, with reversed tail and head Vertices
 * next - pointer to the next Edge in the singly-linked list of Edges with the same type and tail Vertex as one another
 * allNext - pointer to the next Edge in the doubly-linked list pointed to by System.allEdgeFirst
 * allPrev - pointer to the previous Edge in the doubly-linked list pointed to by System.allEdgeFirst
 * inAllEdgeList - boolean flag indicating whether or not this Edge occurs within the doubly-linked list pointed to by 
 *   System.allEdgeFirst
 */
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

/*
 * The EdgeRefList is a way to reference a list of Edges without added pointers within each Edge struct
 *
 * first - pointer to the first EdgeRefListNode within the EdgeRefList
 * last - pointer to the last EdgeRefListNode within the EdgeRefList
 */
struct EdgeRefList {
  EdgeRefListNode * first;
  EdgeRefListNode * last;
};

/*
 * The EdgeRefListNode is one node within the EdgeRefList
 *
 * edge - pointer to an Edge
 * next - pointer to the next EdgeRefListNode within the EdgeRefList
 */
struct EdgeRefListNode {
  Edge * edge;
  EdgeRefListNode * next;
};

/*
 * The IntegerTree data structure is a combination of the T tree specified in the integral portion of the algorithm and a simple 
 * FIFO queue used to specify the order in which Vertices are processed.
 *
 * treeRoot - pointer to the IntegerTreeVertex that is the root of the overall IntegerTree
 * queueNewest - pointer to the most-recently-added member of the IntegerTree FIFO queue
 * queueOldest - pointer to the earliest-added member of the IntegerTree FIFO queue
 * additionsFirst - whenever an absolute constraint is generated by the integer algorithms, it is represented by a new Edge
 *   struct. Edge structs generated in this way are linked together through their next pointers and the beginning of this singly-
 *   linked list is pointed to by additionsFirst, so that these Edges may be properly freed when the IntegerTree is freed.
 */
struct IntegerTree {
  IntegerTreeVertex * treeRoot;
  IntegerTreeVertex * queueNewest;
  IntegerTreeVertex * queueOldest;
  Edge * additionsFirst;
};

/*
 * The IntegerTreeVertex is a single vertex within the IntegerTree, corresponding to a single Vertex within the graph.
 *
 * parent - pointer to this IntegerTreeVertex's parent vertex
 * queueNewer - pointer to the IntegerTreeVertex added to the IntegerTree FIFO queue immediately after this one
 * graphEdges - pointer to the EdgeRefList storing graph Edges listed under this IntegerTreeVertex
 * graphVertex - pointer to the graph Vertex corresponding to this IntegerTreeVertex
 */
struct IntegerTreeVertex {
  IntegerTreeVertex * parent;
  IntegerTreeVertex * queueNewer;
  EdgeRefList * graphEdges;
  Vertex * graphVertex;
};

int main(int argc, char * argv[]);
static void fputEdge(Edge * edge, FILE * output);
static void initializeSystem(void * object, int n, Parser * parser);
static void addConstraint(void * object, Constraint * constraint, Parser * parser);
static void addEdge(System * system, Constraint * constraint);
static void finishSystemCreation(System * system);
static int edgeCompare(const void * edge1, const void * edge2);
static void removeFromAllEdgeList(System * system, Edge * edge);
static bool relaxNetwork(System * system);
static bool relaxEdge(System * system, Edge * e, bool * anyChange);
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

/*
 * main()
 * - handles file input and output
 * - calls utvpiInterpreter's parseFile() function, which calls initializeSystem() and addConstraint() to build the graph 
 *   representation corresponding to the constraint system in the input file
 * - calls finishSystemCreation()
 * - implements UTVPI-LINEAR-FEAS()
 * - calls produceIntegerSolution()
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
      system.negativeCycleEdgeCount++;
    }
    fprintf( output, "\n%d false positives\n", system.falsePositives );
    fprintf( output, "%d main loop iterations\n", system.mainLoopIterations );
    fprintf( output, "%d negative cycle edges\n", system.negativeCycleEdgeCount );
  }
  else{
    fputs("Linear solution:\n", output);
    for(int i = 1; i < system.vertexCount; i++){
      system.graph[i].a = intDivBy2ToHalfInt( system.graph[i].D[WHITE] - system.graph[i].D[BLACK] );
      fprintf(output, "x%i = %.1f\n", i, halfIntToDouble( system.graph[i].a ) );
    }
    fprintf( output, "\n%d false positives\n", system.falsePositives );
    fprintf( output, "%d main loop iterations\n", system.mainLoopIterations );
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
   * diff() takes the difference of two struct timespecs
   * Copied from https://www.guyrutenberg.com/2007/09/22/profiling-code-using-clock_gettime/
   * and modified.
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
 * edge - pointer to an Edge struct to convert to a constraint equation
 * output - FILE pointer to print the constraint equation to
 */
static void fputEdge(Edge * edge, FILE * output){
  if( edge == NULL ){
    fputs("(nil)\n", output);
  }
  else {
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
}

/*
 * initializeSystem() initializes an already-declared System struct, given the number of variables that the system must represent
 * object - a void pointer pointing to an already-declared System struct
 * n - the number of variables that the graph representation held by the System struct must represent
 * parser - pointer to the Parser struct that utvpiInterpreter uses during the input file parsing process, so that parseError() 
 *   can be called, if need be
 */
static void initializeSystem(void * object, int n, Parser * parser){
  System * system = (System *) object;
  system->vertexCount = n + 1;
  system->graph = (Vertex *) malloc( sizeof(Vertex) * system->vertexCount );
  system->n = n;
  system->C = 0;
  system->allEdgeFirst = NULL;
  system->infeasibilityProof = NULL;
  system->T = NULL;
  system->falsePositives = 0;
  system->mainLoopIterations = INT_MAX;
  system->negativeCycleEdgeCount = 0;
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
    system->graph[0].cycleOriginator[i] = NULL;
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
      system->graph[i].cycleOriginator[j] = NULL;
      system->graph[i].first[j] = NULL;
      system->graph[i].edgeCount[j] = 0;
    }
  }
}

/*
 * addConstraint() adds a constraint to the graph representation held by the System struct. Also sets non-source Vertex distance 
 *   and predecessor labels associated with absolute constraints.
 * object - a void pointer pointing to an already-initialized System struct
 * constraint - pointer to a Constraint struct describing a constraint
 * parser - pointer to the Parser struct that utvpiInterpreter uses during the input file parsing process, so that parseError() 
 *   can be called, if need be
 */
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
        system->graph[ constraint->index[0] ].cycleOriginator[WHITE] = system->graph[0].first[WHITE];
        system->graph[ constraint->index[0] ].D[GRAY_FORWARD] = constraint->weight;
        system->graph[ constraint->index[0] ].L[GRAY_FORWARD] = system->graph[0].first[GRAY_FORWARD];
        system->graph[ constraint->index[0] ].cycleOriginator[GRAY_FORWARD] = system->graph[0].first[GRAY_FORWARD];
      }
    }
    else{
      if( constraint->weight < system->graph[ constraint->index[0] ].D[BLACK] ){
        system->graph[ constraint->index[0] ].D[BLACK] = constraint->weight;
        system->graph[ constraint->index[0] ].L[BLACK] = system->graph[0].first[BLACK];
        system->graph[ constraint->index[0] ].cycleOriginator[BLACK] = system->graph[0].first[BLACK];
        system->graph[ constraint->index[0] ].D[GRAY_REVERSE] = constraint->weight;
        system->graph[ constraint->index[0] ].L[GRAY_REVERSE] = system->graph[0].first[GRAY_REVERSE];
        system->graph[ constraint->index[0] ].cycleOriginator[GRAY_REVERSE] = system->graph[0].first[GRAY_REVERSE];
      }
    }
  }
  else{
    addEdge( system, constraint );
  }
}

/*
 * addEdge() adds a pair of Edge structs representing the input constraint to the graph representation held by the System struct.
 *   Each Edge struct will be the reverse of the other. This function can not handle absolute constraints. In the case of an 
 *   absolute constraint, addConstraint() modifies the input Constraint struct twice, passing the modified Constraint struct to 
 *   this function each time.
 * system - pointer to the overall System struct
 * constraint - pointer to a constraint struct describing a constraint, which can not be an absolute constraint
 */
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

/*
 * finishSystemCreation()
 * - sets all distance labels not set by absolute constraints
 * - sorts all edge lists (not including the allEdge doubly-linked list) by increasing head vertex index
 * - removes duplicate edges with equal or greater weight
 * system - pointer to the overall System struct containing the graph representation
 */
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

/*
 * edgeCompare() is the comparison function used by qsort() in finishSystemCreation().
 */
static int edgeCompare(const void * edge1, const void * edge2){
  return (*(Edge **)edge1)->head->index - (*(Edge **)edge2)->head->index;
}

/*
 * removeFromAllEdgeList() removes edge from the allEdge doubly-linked list stored in system.
 * system - pointer to the System struct storing the overall graph representation
 * edge - pointer to the edge to be removed from the allEdge list
 */
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

/*
 * relaxNetwork() implements RELAX-NETWORK(), with proof-backed optimizations. This runs the relaxation loop until 
 *   - no distance-label change occurs during an entire iteration,
 *   - a negative cost cycle is correctly identified, or
 *   - as many iterations as could possibly be necessary for a feasible system are completed
 *   If a distance label changes after the end of the main relaxation loop, a negative cost cycle is detected. The function 
 *   returns true for a linearly feasible system. False otherwise. If at any stage a negative cost cycle is detected, it is placed
 *   in system->infeasibilityProof.
 * system - pointer to the System struct storing the overall graph representation
 */
static bool relaxNetwork(System * system){
  //Lines 3-6 of RELAX-NETWORK() implemented in finishSystemCreation().
  bool anyChange = true;
  for(int r = 1; r <= 2 * system->n && anyChange; r++){
    //printf("Round %d\n", r);
    anyChange = false;
    Edge * e = system->allEdgeFirst;
    while( e != NULL ){
      bool edgeFeasible = relaxEdge(system, e, &anyChange);
      if( !edgeFeasible ){
        system->mainLoopIterations = r;
        return false;
      }
      e = e->allNext;
    }
    if( !anyChange ){
      system->mainLoopIterations = r;
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

/*
 * relaxEdge() implements RELAX-EDGE(), with proof-backed optimizations. If a negative cost cycle is detected, system->
 *   infeasibilityProof is filled with the negative cost cycle and false is returned. Otherwise, true is returned.
 * system - pointer to the System struct storing the overall graph representation
 * e - pointer to the Edge to be relaxed
 * anyChange - pointer to a boolean which will be set to true if any distance label is changed during this call to relaxEdge()
 */
static bool relaxEdge(System * system, Edge * e, bool * anyChange){
  /*fputs("Relaxing ", stdout); fputEdge(e, stdout);
  puts("Edge state before:");
  printf("Head: x%d\n", e->head->index);
  printf("D[WHITE] = %d\n", e->head->D[WHITE]);
  fputs("L[WHITE] = ", stdout); fputEdge(e->head->L[WHITE], stdout);
  fputs("cycleOriginator[WHITE] = ", stdout); fputEdge(e->head->cycleOriginator[WHITE], stdout);
  printf("D[BLACK] = %d\n", e->head->D[BLACK]);
  fputs("L[BLACK] = ", stdout); fputEdge(e->head->L[BLACK], stdout);
  fputs("cycleOriginator[BLACK] = ", stdout); fputEdge(e->head->cycleOriginator[BLACK], stdout);
  printf("D[GRAY_FORWARD] = %d\n", e->head->D[GRAY_FORWARD]);
  fputs("L[GRAY_FORWARD] = ", stdout); fputEdge(e->head->L[GRAY_FORWARD], stdout);
  fputs("cycleOriginator[GRAY_FORWARD] = ", stdout); fputEdge(e->head->cycleOriginator[GRAY_FORWARD], stdout);
  printf("D[GRAY_REVERSE] = %d\n", e->head->D[GRAY_REVERSE]);
  fputs("L[GRAY_REVERSE] = ", stdout); fputEdge(e->head->L[GRAY_REVERSE], stdout);
  fputs("cycleOriginator[GRAY_REVERSE] = ", stdout); fputEdge(e->head->cycleOriginator[GRAY_REVERSE], stdout);
  printf("Tail: x%d\n", e->tail->index);
  printf("D[WHITE] = %d\n", e->tail->D[WHITE]);
  fputs("L[WHITE] = ", stdout); fputEdge(e->tail->L[WHITE], stdout);
  fputs("cycleOriginator[WHITE] = ", stdout); fputEdge(e->tail->cycleOriginator[WHITE], stdout);
  printf("D[BLACK] = %d\n", e->tail->D[BLACK]);
  fputs("L[BLACK] = ", stdout); fputEdge(e->tail->L[BLACK], stdout);
  fputs("cycleOriginator[BLACK] = ", stdout); fputEdge(e->tail->cycleOriginator[BLACK], stdout);
  printf("D[GRAY_FORWARD] = %d\n", e->tail->D[GRAY_FORWARD]);
  fputs("L[GRAY_FORWARD] = ", stdout); fputEdge(e->tail->L[GRAY_FORWARD], stdout);
  fputs("cycleOriginator[GRAY_FORWARD] = ", stdout); fputEdge(e->tail->cycleOriginator[GRAY_FORWARD], stdout);
  printf("D[GRAY_REVERSE] = %d\n", e->tail->D[GRAY_REVERSE]);
  fputs("L[GRAY_REVERSE] = ", stdout); fputEdge(e->tail->L[GRAY_REVERSE], stdout);
  fputs("cycleOriginator[GRAY_REVERSE] = ", stdout); fputEdge(e->tail->cycleOriginator[GRAY_REVERSE], stdout);*/
  switch( e->type ){
  case WHITE:
    if( e->head->D[GRAY_REVERSE] + e->weight < e->tail->D[WHITE] ){
      Edge * priorL = e->tail->L[WHITE];
      e->tail->D[WHITE] = e->head->D[GRAY_REVERSE] + e->weight;
      e->tail->L[WHITE] = e->reverse;
      if( e->tail->index == 0 ){
        //puts("1");
        bool falsePositive = backtrack(system, e->head, GRAY_REVERSE, e->reverse);
        if( !falsePositive ){
          return false;
        }
      }
      else if( priorL != e->reverse ){
        e->tail->cycleOriginator[WHITE] = e->reverse;
      }
      else {
        e->tail->cycleOriginator[WHITE] = e->head->cycleOriginator[GRAY_REVERSE];
        if( e->reverse == e->tail->cycleOriginator[WHITE] ){
          bool falsePositive = backtrack(system, e->head, GRAY_REVERSE, e->reverse);
          if( !falsePositive ){
            return false;
          }
        }
      }
      *anyChange = true;
    }
    if( e->head->D[BLACK] + e->weight < e->tail->D[GRAY_FORWARD] ){
      Edge * priorL = e->tail->L[GRAY_FORWARD];
      e->tail->D[GRAY_FORWARD] = e->head->D[BLACK] + e->weight;
      e->tail->L[GRAY_FORWARD] = e->reverse;
      if( e->tail->index == 0 ){
        //puts("2");
        bool falsePositive = backtrack(system, e->head, BLACK, e->reverse);
        if( !falsePositive ){
          return false;
        }
      }
      else if( priorL != e->reverse ){
        e->tail->cycleOriginator[GRAY_FORWARD] = e->reverse;
      }
      else {
        e->tail->cycleOriginator[GRAY_FORWARD] = e->head->cycleOriginator[BLACK];
        if( e->reverse == e->tail->cycleOriginator[GRAY_FORWARD] ){
          bool falsePositive = backtrack(system, e->head, BLACK, e->reverse);
          if( !falsePositive ){
            return false;
          }
        }
      }
      *anyChange = true;
    }
    if( e->tail->D[GRAY_REVERSE] + e->weight < e->head->D[WHITE] ){
      Edge * priorL = e->head->L[WHITE];
      e->head->D[WHITE] = e->tail->D[GRAY_REVERSE] + e->weight;
      e->head->L[WHITE] = e;
      if( e->head->index == 0 ){
        //puts("3");
        bool falsePositive = backtrack(system, e->tail, GRAY_REVERSE, e);
        if( !falsePositive ){
          return false;
        }
      }
      else if( priorL != e ){
        e->head->cycleOriginator[WHITE] = e;
      }
      else {
        e->head->cycleOriginator[WHITE] = e->tail->cycleOriginator[GRAY_REVERSE];
        if( e == e->head->cycleOriginator[WHITE] ){
          bool falsePositive = backtrack(system, e->tail, GRAY_REVERSE, e);
          if( !falsePositive ){
            return false;
          }
        }
      }
      *anyChange = true;
    }
    if( e->tail->D[BLACK] + e->weight < e->head->D[GRAY_FORWARD] ){
      Edge * priorL = e->head->L[GRAY_FORWARD];
      e->head->D[GRAY_FORWARD] = e->tail->D[BLACK] + e->weight;
      e->head->L[GRAY_FORWARD] = e;
      if( e->head->index == 0 ){
        //puts("4");
        bool falsePositive = backtrack(system, e->tail, BLACK, e);
        if( !falsePositive ){
          return false;
        }
      }
      else if( priorL != e ){
        e->head->cycleOriginator[GRAY_FORWARD] = e;
      }
      else {
        e->head->cycleOriginator[GRAY_FORWARD] = e->tail->cycleOriginator[BLACK];
        if( e == e->head->cycleOriginator[GRAY_FORWARD] ){
          bool falsePositive = backtrack(system, e->tail, BLACK, e);
          if( !falsePositive ){
            return false;
          }
        }
      }
      *anyChange = true;
    }
    break;
  case BLACK:
    if( e->head->D[GRAY_FORWARD] + e->weight < e->tail->D[BLACK] ){
      Edge * priorL = e->tail->L[BLACK];
      e->tail->D[BLACK] = e->head->D[GRAY_FORWARD] + e->weight;
      e->tail->L[BLACK] = e->reverse;
      if( e->tail->index == 0 ){
        //puts("5");
        bool falsePositive = backtrack(system, e->head, GRAY_FORWARD, e->reverse);
        if( !falsePositive ){
          return false;
        }
      }
      else if( priorL != e->reverse ){
        e->tail->cycleOriginator[BLACK] = e->reverse;
      }
      else {
        e->tail->cycleOriginator[BLACK] = e->head->cycleOriginator[GRAY_FORWARD];
        if( e->reverse == e->tail->cycleOriginator[BLACK] ){
          bool falsePositive = backtrack(system, e->head, GRAY_FORWARD, e->reverse);
          if( !falsePositive ){
            return false;
          }
        }
      }
      *anyChange = true;
    }
    if( e->head->D[WHITE] + e->weight < e->tail->D[GRAY_REVERSE] ){
      Edge * priorL = e->tail->L[GRAY_REVERSE];
      e->tail->D[GRAY_REVERSE] = e->head->D[WHITE] + e->weight;
      e->tail->L[GRAY_REVERSE] = e->reverse;
      if( e->tail->index == 0 ){
        //puts("6");
        bool falsePositive = backtrack(system, e->head, WHITE, e->reverse);
        if( !falsePositive ){
          return false;
        }
      }
      else if( priorL != e->reverse ){
        e->tail->cycleOriginator[GRAY_REVERSE] = e->reverse;
      }
      else {
        e->tail->cycleOriginator[GRAY_REVERSE] = e->head->cycleOriginator[WHITE];
        if( e->reverse == e->tail->cycleOriginator[GRAY_REVERSE] ){
          bool falsePositive = backtrack(system, e->head, WHITE, e->reverse);
          if( !falsePositive ){
            return false;
          }
        }
      }
      *anyChange = true;
    }
    if( e->tail->D[GRAY_FORWARD] + e->weight < e->head->D[BLACK] ){
      Edge * priorL = e->head->L[BLACK];
      e->head->D[BLACK] = e->tail->D[GRAY_FORWARD] + e->weight;
      e->head->L[BLACK] = e;
      if( e->head->index == 0 ){
        //puts("7");
        bool falsePositive = backtrack(system, e->tail, GRAY_FORWARD, e);
        if( !falsePositive ){
          return false;
        }
      }
      else if( priorL != e ){
        e->head->cycleOriginator[BLACK] = e;
      }
      else {
        e->head->cycleOriginator[BLACK] = e->tail->cycleOriginator[GRAY_FORWARD];
        if( e == e->head->cycleOriginator[BLACK] ){
          bool falsePositive = backtrack(system, e->tail, GRAY_FORWARD, e);
          if( !falsePositive ){
            return false;
          }
        }
      }
      *anyChange = true;
    }
    if( e->tail->D[WHITE] + e->weight < e->head->D[GRAY_REVERSE] ){
      Edge * priorL = e->head->L[GRAY_REVERSE];
      e->head->D[GRAY_REVERSE] = e->tail->D[WHITE] + e->weight;
      e->head->L[GRAY_REVERSE] = e;
      if( e->head->index == 0 ){
        //puts("8");
        bool falsePositive = backtrack(system, e->tail, WHITE, e);
        if( !falsePositive ){
          return false;
        }
      }
      else if( priorL != e ){
        e->head->cycleOriginator[GRAY_REVERSE] = e;
      }
      else {
        e->head->cycleOriginator[GRAY_REVERSE] = e->tail->cycleOriginator[WHITE];
        if( e == e->head->cycleOriginator[GRAY_REVERSE] ){
          bool falsePositive = backtrack(system, e->tail, WHITE, e);
          if( !falsePositive ){
            return false;
          }
        }
      }
      *anyChange = true;
    }
    break;
  case GRAY_FORWARD:
    if( e->head->D[GRAY_REVERSE] + e->weight < e->tail->D[GRAY_REVERSE] ){
      Edge * priorL = e->tail->L[GRAY_REVERSE];
      e->tail->D[GRAY_REVERSE] = e->head->D[GRAY_REVERSE] + e->weight;
      e->tail->L[GRAY_REVERSE] = e->reverse;
      if( e->tail->index == 0 ){
        //puts("9");
        bool falsePositive = backtrack(system, e->head, GRAY_REVERSE, e->reverse);
        if( !falsePositive ){
          return false;
        }
      }
      else if( priorL != e->reverse ){
        e->tail->cycleOriginator[GRAY_REVERSE] = e->reverse;
      }
      else {
        e->tail->cycleOriginator[GRAY_REVERSE] = e->head->cycleOriginator[GRAY_REVERSE];
        if( e->reverse == e->tail->cycleOriginator[GRAY_REVERSE] ){
          bool falsePositive = backtrack(system, e->head, GRAY_REVERSE, e->reverse);
          if( !falsePositive ){
            return false;
          }
        }
      }
      *anyChange = true;
    }
    if( e->head->D[BLACK] + e->weight < e->tail->D[BLACK] ){
      Edge * priorL = e->tail->L[BLACK];
      e->tail->D[BLACK] = e->head->D[BLACK] + e->weight;
      e->tail->L[BLACK] = e->reverse;
      if( e->tail->index == 0 ){
        //puts("10");
        bool falsePositive = backtrack(system, e->head, BLACK, e->reverse);
        if( !falsePositive ){
          return false;
        }
      }
      else if( priorL != e->reverse ){
        e->tail->cycleOriginator[BLACK] = e->reverse;
      }
      else {
        e->tail->cycleOriginator[BLACK] = e->head->cycleOriginator[BLACK];
        if( e->reverse == e->tail->cycleOriginator[BLACK] ){
          bool falsePositive = backtrack(system, e->head, BLACK, e->reverse);
          if( !falsePositive ){
            return false;
          }
        }
      }
      *anyChange = true;
    }
    if( e->tail->D[GRAY_FORWARD] + e->weight < e->head->D[GRAY_FORWARD] ){
      Edge * priorL = e->head->L[GRAY_FORWARD];
      e->head->D[GRAY_FORWARD] = e->tail->D[GRAY_FORWARD] + e->weight;
      e->head->L[GRAY_FORWARD] = e;
      if( e->head->index == 0 ){
        //puts("11");
        bool falsePositive = backtrack(system, e->tail, GRAY_FORWARD, e);
        if( !falsePositive ){
          return false;
        }
      }
      else if( priorL != e ){
        e->head->cycleOriginator[GRAY_FORWARD] = e;
      }
      else {
        e->head->cycleOriginator[GRAY_FORWARD] = e->tail->cycleOriginator[GRAY_FORWARD];
        if( e == e->head->cycleOriginator[GRAY_FORWARD] ){
          bool falsePositive = backtrack(system, e->tail, GRAY_FORWARD, e);
          if( !falsePositive ){
            return false;
          }
        }
      }
      *anyChange = true;
    }
    if( e->tail->D[WHITE] + e->weight < e->head->D[WHITE] ){
      Edge * priorL = e->head->L[WHITE];
      e->head->D[WHITE] = e->tail->D[WHITE] + e->weight;
      e->head->L[WHITE] = e;
      if( e->head->index == 0 ){
        //puts("12");
        bool falsePositive = backtrack(system, e->tail, WHITE, e);
        if( !falsePositive ){
          return false;
        }
      }
      else if( priorL != e ){
        e->head->cycleOriginator[WHITE] = e;
      }
      else {
        e->head->cycleOriginator[WHITE] = e->tail->cycleOriginator[WHITE];
        if( e == e->head->cycleOriginator[WHITE] ){
          bool falsePositive = backtrack(system, e->tail, WHITE, e);
          if( !falsePositive ){
            return false;
          }
        }
      }
      *anyChange = true;
    }
    break;
  case GRAY_REVERSE:
    if( e->head->D[GRAY_FORWARD] + e->weight < e->tail->D[GRAY_FORWARD] ){
      Edge * priorL = e->tail->L[GRAY_FORWARD];
      e->tail->D[GRAY_FORWARD] = e->head->D[GRAY_FORWARD] + e->weight;
      e->tail->L[GRAY_FORWARD] = e->reverse;
      if( e->tail->index == 0 ){
        //puts("13");
        bool falsePositive = backtrack(system, e->head, GRAY_FORWARD, e->reverse);
        if( !falsePositive ){
          return false;
        }
      }
      else if( priorL != e->reverse ){
        e->tail->cycleOriginator[GRAY_FORWARD] = e->reverse;
      }
      else {
        e->tail->cycleOriginator[GRAY_FORWARD] = e->head->cycleOriginator[GRAY_FORWARD];
        if( e->reverse == e->tail->cycleOriginator[GRAY_FORWARD] ){
          bool falsePositive = backtrack(system, e->head, GRAY_FORWARD, e->reverse);
          if( !falsePositive ){
            return false;
          }
        }
      }
      *anyChange = true;
    }
    if( e->head->D[WHITE] + e->weight < e->tail->D[WHITE] ){
      Edge * priorL = e->tail->L[WHITE];
      e->tail->D[WHITE] = e->head->D[WHITE] + e->weight;
      e->tail->L[WHITE] = e->reverse;
      if( e->tail->index == 0 ){
        //puts("14");
        bool falsePositive = backtrack(system, e->head, WHITE, e->reverse);
        if( !falsePositive ){
          return false;
        }
      }
      else if( priorL != e->reverse ){
        e->tail->cycleOriginator[WHITE] = e->reverse;
      }
      else {
        e->tail->cycleOriginator[WHITE] = e->head->cycleOriginator[WHITE];
        if( e->reverse == e->tail->cycleOriginator[WHITE] ){
          bool falsePositive = backtrack(system, e->head, WHITE, e->reverse);
          if( !falsePositive ){
            return false;
          }
        }
      }
      *anyChange = true;
    }
    if( e->tail->D[GRAY_REVERSE] + e->weight < e->head->D[GRAY_REVERSE] ){
      Edge * priorL = e->head->L[GRAY_REVERSE];
      e->head->D[GRAY_REVERSE] = e->tail->D[GRAY_REVERSE] + e->weight;
      e->head->L[GRAY_REVERSE] = e;
      if( e->head->index == 0 ){
        //puts("15");
        bool falsePositive = backtrack(system, e->tail, GRAY_REVERSE, e);
        if( !falsePositive ){
          return false;
        }
      }
      else if( priorL != e ){
        e->head->cycleOriginator[GRAY_REVERSE] = e;
      }
      else {
        e->head->cycleOriginator[GRAY_REVERSE] = e->tail->cycleOriginator[GRAY_REVERSE];
        if( e == e->head->cycleOriginator[GRAY_REVERSE] ){
          bool falsePositive = backtrack(system, e->tail, GRAY_REVERSE, e);
          if( !falsePositive ){
            return false;
          }
        }
      }
      *anyChange = true;
    }
    if( e->tail->D[BLACK] + e->weight < e->head->D[BLACK] ){
      Edge * priorL = e->head->L[BLACK];
      e->head->D[BLACK] = e->tail->D[BLACK] + e->weight;
      e->head->L[BLACK] = e;
      if( e->head->index == 0 ){
        //puts("16");
        bool falsePositive = backtrack(system, e->tail, BLACK, e);
        if( !falsePositive ){
          return false;
        }
      }
      else if( priorL != e ){
        e->head->cycleOriginator[BLACK] = e;
      }
      else {
        e->head->cycleOriginator[BLACK] = e->tail->cycleOriginator[BLACK];
        if( e == e->head->cycleOriginator[BLACK] ){
          bool falsePositive = backtrack(system, e->tail, BLACK, e);
          if( !falsePositive ){
            return false;
          }
        }
      }
      *anyChange = true;
    }
  }
  /*puts("Edge state after:");
  printf("Head: x%d\n", e->head->index);
  printf("D[WHITE] = %d\n", e->head->D[WHITE]);
  fputs("L[WHITE] = ", stdout); fputEdge(e->head->L[WHITE], stdout);
  fputs("cycleOriginator[WHITE] = ", stdout); fputEdge(e->head->cycleOriginator[WHITE], stdout);
  printf("D[BLACK] = %d\n", e->head->D[BLACK]);
  fputs("L[BLACK] = ", stdout); fputEdge(e->head->L[BLACK], stdout);
  fputs("cycleOriginator[BLACK] = ", stdout); fputEdge(e->head->cycleOriginator[BLACK], stdout);
  printf("D[GRAY_FORWARD] = %d\n", e->head->D[GRAY_FORWARD]);
  fputs("L[GRAY_FORWARD] = ", stdout); fputEdge(e->head->L[GRAY_FORWARD], stdout);
  fputs("cycleOriginator[GRAY_FORWARD] = ", stdout); fputEdge(e->head->cycleOriginator[GRAY_FORWARD], stdout);
  printf("D[GRAY_REVERSE] = %d\n", e->head->D[GRAY_REVERSE]);
  fputs("L[GRAY_REVERSE] = ", stdout); fputEdge(e->head->L[GRAY_REVERSE], stdout);
  fputs("cycleOriginator[GRAY_REVERSE] = ", stdout); fputEdge(e->head->cycleOriginator[GRAY_REVERSE], stdout);
  printf("Tail: x%d\n", e->tail->index);
  printf("D[WHITE] = %d\n", e->tail->D[WHITE]);
  fputs("L[WHITE] = ", stdout); fputEdge(e->tail->L[WHITE], stdout);
  fputs("cycleOriginator[WHITE] = ", stdout); fputEdge(e->tail->cycleOriginator[WHITE], stdout);
  printf("D[BLACK] = %d\n", e->tail->D[BLACK]);
  fputs("L[BLACK] = ", stdout); fputEdge(e->tail->L[BLACK], stdout);
  fputs("cycleOriginator[BLACK] = ", stdout); fputEdge(e->tail->cycleOriginator[BLACK], stdout);
  printf("D[GRAY_FORWARD] = %d\n", e->tail->D[GRAY_FORWARD]);
  fputs("L[GRAY_FORWARD] = ", stdout); fputEdge(e->tail->L[GRAY_FORWARD], stdout);
  fputs("cycleOriginator[GRAY_FORWARD] = ", stdout); fputEdge(e->tail->cycleOriginator[GRAY_FORWARD], stdout);
  printf("D[GRAY_REVERSE] = %d\n", e->tail->D[GRAY_REVERSE]);
  fputs("L[GRAY_REVERSE] = ", stdout); fputEdge(e->tail->L[GRAY_REVERSE], stdout);
  fputs("cycleOriginator[GRAY_REVERSE] = ", stdout); fputEdge(e->tail->cycleOriginator[GRAY_REVERSE], stdout);*/
  return true;
}

/*
 * backtrack() implements BACKTRACK(), with modifications to handle false positives generated by cycle-originators. If the 
 *   predecessor structure is backtracked through, up to the source node, indicated with a NULL predecessor, the E arrays are
 *   reinitialized, and true is returned. Otherwise, system->infeasibilityProof is filled with a negative cost cycle and false is
 *   returned.
 * system - pointer to the System struct storing the overall graph representation
 * x_i - pointer to the initial Vertex
 * t - initial path type
 * e - pointer to the initial edge
 */
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
  while( x_c->E[a_c] == NULL && x_c->L[t_c] != NULL ){
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
  
  if( x_c->L[t_c] == NULL ){ //due to cycleOriginator false positive
    for(int i = 0; i < system->vertexCount; i++){
      system->graph[i].E[NEG_ONE] = NULL;
      system->graph[i].E[POS_ONE] = NULL;
    }
    system->falsePositives++;
    return true;
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

/*
 * produceIntegerSolution() implements PRODUCE-INTEGER-SOLUTION(). If the system is not integrally feasible, system->
 *   infeasibilityProof is filled with a proof of integral infeasibility and false is returned. Otherwise, true is returned.
 * system - pointer to the System struct storing the overall graph representation
 */
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

/*
 * forcedRounding() implements FORCED-ROUNDING(). If a contradiction is detected within this function, the proof of integral 
 *   infeasibility arising from it is placed in system->infeasibilityProof, and false is returned. Otherwise, true is returned.
 * system - pointer to the System struct storing the overall graph representation
 * x_i - pointer to the current variable
 */
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

/*
 * optionalRoundings() implements OPTIONAL-ROUNDINGS(). If contradictions are detected when rounding a variable down and when 
 *   rounding the same variable up, the proof of integral infeasibility arising from this is placed in system->infeasibilityProof,
 *   and false is returned. Otherwise, true is returned.
 * system - pointer to the System struct storing the overall graph representation
 */
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

/*
 * checkDependencies() implements CHECK-DEPENDENCIES().
 * T - pointer to the IntegerTree to be expanded with dependencies
 * x_i - pointer to the variable whose dependencies are checked. Must never represent a variable where x_i->a is integral. The 
 *   first for loop in produceIntegerSolution() ensures this can never occur.
 * integerType - type of Z value being filled in during this call
 */
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

/*
 * checkAllConstraints() checks all constraints to determine if they are satisfied by integral assignments. If all are, true is 
 *   returned. Otherwise, system->infeasibilityProof is filled with the set of edges leading to a contradiction, and false is 
 *   returned. The edges added to system->infeasibilityProof during one call to this function make up no more than one helf of a 
 *   full proof of integral infeasibility.
 * system - pointer to the System struct storing the overall graph representation
 * toVertex - pointer to the vertex where integerTreeBacktrack() should stop when backtracking up the IntegerTree
 * integerType - type of Z value being checked against constraints during this call
 */
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

/*
 * systemSubset() removes all edges from the system where tail->Z[FINAL] and head->Z[FINAL] don't both equal INT_MAX, in 
 *   accordance with line 23 of PRODUCE-INTEGER-SOLUTION().
 * system - pointer to the System struct storing the overall graph representation
 */
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

/*
 * generateAbsoluteConstraint() allocates and initializes an Edge struct representing an absolute constraint. This Edge struct can
 *   only be used for the purpose of representing a constraint derived during the production of an integral solution. The Edge 
 *   struct is not added into the overall graph representation, and is instead added into a list of similarly-generated absolute
 *   constraint Edges tied together through their next pointers and referenced by the additionsFirst pointer in the IntegerTree
 *   struct. Returns a pointer to the newly allocated and initialized Edge struct.
 * system - pointer to the System struct storing the overall graph representation
 * x_i - pointer to the vertex involved in the absolute constraint
 * weight - weight of the absolute constraint
 * type - type of absolute constraint edge. WHITE or GRAY_REVERSE for positive variable coefficient. BLACK or GRAY_FORWARD for 
 *   negative variable coefficient.
 */
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

/*
 * generateIntegerTree() allocates and initializes an IntegerTree rooted at x0, returning a pointer to the new IntegerTree.
 * system - pointer to the System struct storing the overall graph representation
 */
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

/*
 * pollIntegerTreeQueue() polls the FIFO queue stored within the input IntegerTree, returning a pointer to the next Vertex to 
 *   check for dependencies. Returns NULL if the FIFO queue is empty. (tree->queueNewest is never set to NULL, so that all 
 *   IntegerTreeVertices may be tied together through their queueNewer pointers, for freeing later.)
 * tree - pointer to the IntegerTree whose FIFO queue is to be polled
 */
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

/*
 * expandIntegerTree() adds a new IntegerTreeVertex corresponding to active to the input IntegerTree as the child of the 
 *   IntegerTreeVertex corresponding to parent, if there is not already an IntegerTreeVertex corresponding to active. The new 
 *   IntegerTreeVertex is also added to the FIFO queue running through the overall IntegerTree. It then adds references to edges 
 *   edge0, edge1, and edge2 under the IntegerTreeVertex corresponding to active. The first of edge0, edge1, and edge2 set to NULL
 *   will cause all further edge inputs to be ignored. If only one edge is to be added under the IntegerTreeVertex corresponding 
 *   to active, edge0 should be set to it, and edge1 and edge2 should be set to NULL.
 * T - pointer to the IntegerTree to be expanded
 * active - pointer to the graph Vertex which should have a corresponding IntegerTreeVertex added for it, if one is not already 
 *   present
 * parent - pointer to the graph Vertex whose corresponding IntegerTreeVertex should serve as the parent IntegerTreeVertex for the
 *   active graph Vertex
 * edge0 - pointer to the first Edge to be added under the IntegerTreeVertex corresponding to the active graph Vertex. If NULL, no
 *   edges are added.
 * edge1 - pointer to the second Edge to be added under the IntegerTreeVertex corresponding to the active graph Vertex. If NULL,
 *   only edge0 is added.
 * edge2 - pointer to the third Edge to be added under the IntegerTreeVertex corresponding to the active graph Vertex. If NULL, 
 *   only edge0 and edge1 are added.
 */
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

/*
 * integerTreeBacktrack() copies Edge pointers from the IntegerTree data structure to the input list, backtracking from fromVertex
 *   through parent vertices, to toVertex.
 * list - pointer to an EdgeRefList to which Edge pointers are copied from the IntegerTree data structure
 * fromVertex - pointer to the IntegerTreeVertex where backtracking begins
 * toVertex - pointer to the IntegerTreeVertex where backtracking ends
 * includeToVertex - true if Edge pointers from the toVertex should be added to list. false otherwise.
 */
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

/*
 * copyTreeEdgesToList() copies any Edge pointers listed under itv to list.
 * list - pointer to an EdgeRefList to copy Edge pointers to
 * itv - pointer to an IntegerTreeVertex to copy Edge pointers from
 */
static void copyTreeEdgesToList(EdgeRefList * list, IntegerTreeVertex * itv){
  EdgeRefListNode * treeEdge = itv->graphEdges->first;
  while( treeEdge != NULL ){
    edgeRefListPrepend( list, treeEdge->edge );
    treeEdge = treeEdge->next;
  }
}

/*
 * freeIntegerTree() frees the input IntegerTree, including all Edges representing absolute constraints derived within this 
 *   IntegerTree. Also sets integerTreeVertex in each associated graph Vertex to NULL.
 * tree - pointer to the IntegerTree to be freed
 */
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

/*
 * generateEdgeRefList() allocates and initializes an EdgeRefList, returning a pointer to it.
 */
static EdgeRefList * generateEdgeRefList(){
  EdgeRefList * newERL = (EdgeRefList *) malloc( sizeof(EdgeRefList) );
  newERL->first = NULL;
  newERL->last = NULL;
  return newERL;
}

/*
 * edgeRefListAppend() adds an Edge pointer to the end of an EdgeRefList.
 * erl - pointer to an EdgeRefList to append an Edge pointer to
 * edge - pointer to an Edge, which will be added to the end of erl
 */
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

/*
 * edgeRefListPrepend() adds an Edge pointer to the beginning of an EdgeRefList.
 * erl - pointer to an EdgeRefList to prepend with an Edge pointer
 * edge - pointer to an Edge, which will be added to the beginning of erl
 */
static void edgeRefListPrepend(EdgeRefList * erl, Edge * edge){
  EdgeRefListNode * newERLN = (EdgeRefListNode *) malloc( sizeof(EdgeRefListNode) );
  newERLN->edge = edge;
  newERLN->next = erl->first;
  erl->first = newERLN;
  if( erl->last == NULL ){
    erl->last = newERLN;
  }
}

/*
 * freeEdgeRefList() frees the input EdgeRefList.
 * erl - pointer to the EdgeRefList to be freed
 */
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

/*
 * freeSystem() frees the graph representation stored within system, along with system->infeasibilityProof and system->T, if 
 *   those have yet to be freed.
 * system - pointer to the System struct, whose contents are to be freed
 */
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
