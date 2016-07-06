/*
 * subWojLin.c
 * The Subramani-Wojciechowski linear UTVPI system solver
 *
 * Call with [executable] [input file] {output file}
 * [input file] must be properly formatted to be read by utvpiInterpreter.h
 * {output file} will contain a linear solution, if one exists. Otherwise, a proof of linear infeasibility - a negative cost cycle
 *   - will be output. If {output file} is not specified, output will be to stdout.
 *
 * Modifications have been made to the linear portion of the algorithm implementation in subWojInt and subWojIntOpt in the time 
 * since this implementation was completed.
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
 * GRAY_FORWARD corresponds to an Edge struct representing a gray edge, where the tail pointer points to the vertex at the tail of
 * the arrow shown in the algorithm, and the head pointer points to the vertex at the head of the arrow shown in the algorithm. 
 * GRAY_REVERSE corresponds to the opposite case. A right-facing arrow in the algorithm corresponds to GRAY_FORWARD, and a left-
 * facing arrow corresponds to GRAY_REVERSE.
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

typedef struct System System;
typedef struct Vertex Vertex;
typedef struct Edge Edge;
typedef struct EdgeRefList EdgeRefList;

/*
 * The System struct contains references to all dynamically-allocated data structures used within the system solver, along with
 * other information about the system.
 *
 * graph - pointer to an array of Vertex structs, defining the structure of the graph representation
 * vertexCount - the number of Vertex structs within graph
 * n - the number of variables represented within the system - always one less than vertexCount, but used for the sake of clarity
 * C - the greatest absolute value of any edge weight which entered the system as part of a constraint
 * allEdgeFirst - pointer to the first Edge in a singly-linked list where only one Edge struct representing each edge is present
 */
struct System {
  Vertex * graph;
  int vertexCount;
  int n;
  int C;
  Edge * allEdgeFirst;
};

/*
 * The Vertex struct contains all information about a specific variable, represented by a vertex within the graph.
 *
 * index - index of the variable
 * L - pointer to an edge between this vertex and the vertex's predecessor vertex, for each EdgeType-color path between the source
 *   vertex and this vertex
 * D - distance label for each EdgeType-color path between the source vertex and this vertex
 * E - pointer to the last Edge traversed on the path from x_i to x, in backtrack()
 * x - half-integral linear-solution variable value
 * first - pointer to the first Edge of each EdgeType whose tail Vertex is this Vertex. The remainder of such edges are connected
 *   together in a singly-linked list
 */
struct Vertex {
  int index;
  Edge * L[EDGE_TYPE_COUNT];
  int D[EDGE_TYPE_COUNT];
  Edge * E[BACKTRACKING_INDEX_COUNT];
  half_int x;
  Edge * first[EDGE_TYPE_COUNT];
};

/*
 * The Edge struct contains all information about a specific constraint, represented by an edge within the graph.
 *
 * weight - the weight of the edge
 * type - the color type of the edge
 * tail - pointer to the Edge's tail Vertex
 * head - pointer to the Edge's head Vertex
 * reverse - pointer to the other Edge representing the same constraint, with reversed tail and head Vertices
 * next - pointer to the next Edge in the singly-linked list of Edges with the same type and tail Vertex as one another
 * allNext - pointer to the next Edge in the singly-linked list pointed to by System.allEdgeFirst
 */
struct Edge {
  int weight;
  EdgeType type;
  Vertex * tail;
  Vertex * head;
  Edge * reverse;
  Edge * next;
  Edge * allNext;
};

/*
 * The EdgeRefList is a way to reference a list of Edges without added pointers within each Edge struct
 *
 * edge - pointer to an Edge
 * next - pointer to the next EdgeRefList struct within the overall EdgeRefList
 */
struct EdgeRefList {
  Edge * edge;
  EdgeRefList * next;
};

int main(int argc, char * argv[]);
static void fputEdge(Edge * edge, FILE * output);
//static EdgeType reverseEdgeType(EdgeType input);
static void initializeSystem(void * object, int n, Parser * parser);
static void addConstraint(void * object, Constraint * constraint, Parser * parser);
static void addEdge(System * system, Constraint * constraint);
static void finishSystemCreation(System * system);
static EdgeRefList * relaxNetwork(System * system);
static void relaxEdge(Edge * e);
static EdgeRefList * backtrack(Vertex * x_i, EdgeType t, Edge * e);
static void cleanupSystem(System * system);

/*
 * main()
 * - handles file input and output
 * - calls utvpiInterpreter's parseFile() function, which calls initializeSystem() and addConstraint() to build the graph 
 *   representation corresponding to the constraint system in the input file
 * - calls finishSystemCreation()
 * - implements UTVPI-LINEAR-FEAS()
 * - prints profiling information formatted for csv-file input to stdout
 */
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
  finishSystemCreation(&system);
  time_t afterSystemCreation = clock();
  EdgeRefList * R = relaxNetwork(&system);
  if( R != NULL ){
    fputs("The following negative gray cycle was detected:\n", output);
    while( R != NULL ){
      fputEdge(R->edge, output);
      EdgeRefList * oldR = R;
      R = R->next;
      free(oldR);
    }
  }
  else{
    fputs("Linear solution:\n", output);
    for(int i = 1; i < system.vertexCount; i++){
      system.graph[i].x = intDivBy2ToHalfInt( system.graph[i].D[WHITE] - system.graph[i].D[BLACK] );
      fprintf(output, "x%i = %.1f\n", i, halfIntToDouble( system.graph[i].x ) );
    }
  }
  fclose(output);
  time_t afterSolutionOutput = clock();
  cleanupSystem(&system);
  time_t afterCleanup = clock();
  printf("%f,", ((double)(afterSystemCreation - beginning))/CLOCKS_PER_SEC );
  printf("%f,,", ((double)(afterSolutionOutput - afterSystemCreation))/CLOCKS_PER_SEC );
  printf("%f,", ((double)(afterCleanup - afterSolutionOutput))/CLOCKS_PER_SEC );
  printf("%f,", ((double)(afterCleanup - beginning))/CLOCKS_PER_SEC );
  return 0;
}

/*
 * fputEdge() prints the constraint equation corresponding to edge to output
 * edge - pointer to an Edge struct to convert to a constraint equation
 * output - FILE pointer to print the constraint equation to
 */
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

/*static EdgeType reverseEdgeType(EdgeType input){ /////////
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
}*/

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
  system->graph[0].index = 0;
  system->graph[0].E[NEG_ONE] = NULL;
  system->graph[0].E[POS_ONE] = NULL;
  system->graph[0].x = 0;
  for(EdgeType i = WHITE; i <= GRAY_REVERSE; i++){
    system->graph[0].L[i] = NULL;
    system->graph[0].D[i] = 0;
    system->graph[0].first[i] = NULL;
  }
  for(int i = 1; i < system->vertexCount; i++){
    system->graph[i].index = i;
    system->graph[i].E[NEG_ONE] = NULL;
    system->graph[i].E[POS_ONE] = NULL;
    system->graph[i].x = 0;
    for(EdgeType j = WHITE; j <= GRAY_REVERSE; j++){
      system->graph[i].L[j] = NULL;
      system->graph[i].D[j] = INT_MAX;
      system->graph[i].first[j] = NULL;
    }
  }
}

/*
 * addConstraint() adds a constraint to the graph representation held by the System struct. Also sets distance and predecessor 
 *   labels associated with absolute constraints.
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
    addEdge(system, constraint);
    constraint->sign[1] = CONSTRAINT_MINUS;
    addEdge(system, constraint);
    constraint->sign[1] = CONSTRAINT_NONE;
    //Calls to addEdge place the new edges at first[edgeType] for both vertices
    if( constraint->sign[0] == CONSTRAINT_PLUS ){
      if( constraint->weight < system->graph[ constraint->index[0] ].D[WHITE] ){
        system->graph[ constraint->index[0] ].D[WHITE] = constraint->weight;
        system->graph[ constraint->index[0] ].L[WHITE] = system->graph[0].first[WHITE];
        system->graph[ constraint->index[0] ].D[GRAY_FORWARD] = constraint->weight;
        system->graph[ constraint->index[0] ].L[GRAY_FORWARD] = system->graph[0].first[GRAY_FORWARD]; //.first[GRAY_REVERSE]; ?
      }
    }
    else{
      if( constraint->weight < system->graph[ constraint->index[0] ].D[BLACK] ){
        system->graph[ constraint->index[0] ].D[BLACK] = constraint->weight;
        system->graph[ constraint->index[0] ].L[BLACK] = system->graph[0].first[BLACK];
        system->graph[ constraint->index[0] ].D[GRAY_REVERSE] = constraint->weight;
        system->graph[ constraint->index[0] ].L[GRAY_REVERSE] = system->graph[0].first[GRAY_REVERSE]; //.first[GRAY_FORWARD]; ?
      }
    }
  }
  else{
    addEdge(system, constraint);
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
  newEdges[0]->allNext = system->allEdgeFirst;
  system->allEdgeFirst = newEdges[0];
  newEdges[1]->next = newEdges[1]->tail->first[ newEdges[1]->type ];
  newEdges[1]->tail->first[ newEdges[1]->type ] = newEdges[1];
  newEdges[1]->allNext = NULL;
}

/*
 * finishSystemCreation() sets all distance labels not set by an absolute constraint
 * system - pointer to the overall System struct containing the graph representation
 */
static void finishSystemCreation(System * system){
  int sourceEdgeWeights = (2 * system->n + 1) * system->C;
  for(EdgeType i = WHITE; i <= GRAY_REVERSE; i++){
    for(int j = 1; j < system->vertexCount; j++){
      if( system->graph[j].D[i] == INT_MAX ){
        system->graph[j].D[i] = sourceEdgeWeights;
      }
    }
  }
}

/*
 * relaxNetwork() implements RELAX-NETWORK(). The function returns NULL for a linearly feasible system. Otherwise, a pointer to an
 *   EdgeRefList containing a negative cost cycle is returned.
 * system - pointer to the System struct storing the overall graph representation
 */
static EdgeRefList * relaxNetwork(System * system){
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
        return backtrack(e->head, GRAY_REVERSE, e->reverse);
      }
      if( e->head->D[BLACK] + e->weight < e->tail->D[GRAY_FORWARD] ){
        return backtrack(e->head, BLACK, e->reverse);
      }
      if( e->tail->D[GRAY_REVERSE] + e->weight < e->head->D[WHITE] ){
        return backtrack(e->tail, GRAY_REVERSE, e);
      }
      if( e->tail->D[BLACK] + e->weight < e->head->D[GRAY_FORWARD] ){
        return backtrack(e->tail, BLACK, e);
      }
      e = e->next;
    }
  }
  for(int i = 0; i < system->vertexCount; i++){
    Edge * e = system->graph[i].first[BLACK];
    while( e != NULL ){
      if( e->head->D[GRAY_FORWARD] + e->weight < e->tail->D[BLACK] ){
        return backtrack(e->head, GRAY_FORWARD, e->reverse);
      }
      if( e->head->D[WHITE] + e->weight < e->tail->D[GRAY_REVERSE] ){
        return backtrack(e->head, WHITE, e->reverse);
      }
      if( e->tail->D[GRAY_FORWARD] + e->weight < e->head->D[BLACK] ){
        return backtrack(e->tail, GRAY_FORWARD, e);
      }
      if( e->tail->D[WHITE] + e->weight < e->head->D[GRAY_REVERSE] ){
        return backtrack(e->tail, WHITE, e);
      }
      e = e->next;
    }
  }
  for(int i = 0; i < system->vertexCount; i++){
    Edge * e = system->graph[i].first[GRAY_FORWARD];
    while( e != NULL ){
      if( e->head->D[GRAY_REVERSE] + e->weight < e->tail->D[GRAY_REVERSE] ){
        return backtrack(e->head, GRAY_REVERSE, e->reverse);
      }
      if( e->head->D[BLACK] + e->weight < e->tail->D[BLACK] ){
        return backtrack(e->head, BLACK, e->reverse);
      }
      if( e->tail->D[GRAY_FORWARD] + e->weight < e->head->D[GRAY_FORWARD] ){
        return backtrack(e->tail, GRAY_FORWARD, e);
      }
      if( e->tail->D[WHITE] + e->weight < e->head->D[WHITE] ){
        return backtrack(e->tail, WHITE, e);
      }
      e = e->next;
    }
  }
  for(int i = 0; i < system->vertexCount; i++){
    Edge * e = system->graph[i].first[GRAY_REVERSE];
    while( e != NULL ){
      if( e->head->D[GRAY_FORWARD] + e->weight < e->tail->D[GRAY_FORWARD] ){
        return backtrack(e->head, GRAY_FORWARD, e->reverse);
      }
      if( e->head->D[WHITE] + e->weight < e->tail->D[WHITE] ){
        return backtrack(e->head, WHITE, e->reverse);
      }
      if( e->tail->D[GRAY_REVERSE] + e->weight < e->head->D[GRAY_REVERSE] ){
        return backtrack(e->tail, GRAY_REVERSE, e);
      }
      if( e->tail->D[BLACK] + e->weight < e->head->D[BLACK] ){
        return backtrack(e->tail, BLACK, e);
      }
      e = e->next;
    }
  }
  return NULL;
}

/*
 * relaxEdge() implements RELAX-EDGE().
 * e - pointer to the Edge to be relaxed
 */
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

/*
 * backtrack() implements BACKTRACK(). This function returns a pointer to an EdgeRefList storing a negative cost cycle.
 * x_i - pointer to the initial Vertex
 * t - initial path type
 * e - pointer to the initial edge
 */
static EdgeRefList * backtrack(Vertex * x_i, EdgeType t, Edge * e){
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
  EdgeRefList * R = (EdgeRefList *) malloc( sizeof(EdgeRefList) );
  R->edge = e_c;
  R->next = NULL;
  EdgeRefList * R_last = R;
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
    EdgeRefList * newEdgeRef = (EdgeRefList *) malloc( sizeof(EdgeRefList) );
    newEdgeRef->edge = x_c->E[a_c];
    newEdgeRef->next = NULL;
    R_last->next = newEdgeRef;
    R_last = newEdgeRef;
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
  return R;
}

/*
 * cleanupSystem() frees the graph representation stored within system.
 * system - pointer to the System struct, whose contents are to be freed
 */
static void cleanupSystem(System * system){
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
  free( system->graph );
}
