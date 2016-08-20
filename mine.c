/*
 * mine.c
 * The Miné linear UTVPI system solver, with simple Bellman-Ford implementation
 * The Octagon Abstract Domain - Antoine Miné
 * This serves as an implementation of only the linear portion of the lahiri* implementations.
 *
 * Call with [executable] [input file] {output file}
 * [input file] must be properly formatted to be read by utvpiInterpreter.h
 * {output file} will contain a linear solution, if one exists, followed by an integral solution, if one exists. If the system is
 *   not linearly feasible, a proof of linear infeasibility - a negative cost cycle - will be output. If the system is linearly 
 *   feasible, but not integrally feasible, a proof of integral infeasibility will be output. If {output file} is not specified,
 *   output will be to stdout.
 *
 * Modifications have been made to the linear portion of the algorithm implementation in lahiri and onward in the time since this
 * implementation was completed.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "utvpiInterpreter.h"
#include "halfint.h"

/*
 * This algorithm represents each variable with two vertices, one representing the positive occurrence of the variable, and the 
 * other representing its negative occurrence. VertexSign specifies the sign of each Vertex.
 */
#define VERTEX_SIGN_COUNT 2
typedef enum VertexSign {
  POSITIVE,
  NEGATIVE,
} VertexSign;

typedef struct System System;
typedef struct Vertex Vertex;
typedef struct Edge Edge;

/*
 * The System struct contains the graph representation used by the system solver, along with other information about the system.
 *
 * graph - pointers to POSITIVE and NEGATIVE arrays of Vertex structs, defining the structure of the graph representation
 * n - the number of variables represented within the system, equal to the number of Vertex structs stored within each graph array
 */
struct System {
  Vertex * graph[VERTEX_SIGN_COUNT];
  int n;
};

/*
 * The Vertex struct contains all information about one signed occurrence of a variable, represented by a vertex within the graph.
 *
 * index - index of the variable
 * sign - the sign of the variable occurrence this Vertex represents
 * L - pointer to the predecessor Edge between this Vertex and the Vertex's predecessor Vertex
 * D - distance label
 * first - pointer to the first Edge whose tail Vertex is this Vertex. The remainder of such edges are connected together in a 
 *   singly-linked list
 * x - half-integral linear-solution variable value, not particular to a signed occurrence
 */
struct Vertex {
  int index;
  VertexSign sign;
  Edge * L;
  int D;
  Edge * first;
  half_int x;
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

int main(int argc, char * argv[]);
static void fputEdge(Edge * edge, FILE * output);
static void initializeSystem(void * object, int n, Parser * parser);
static void addConstraint(void * object, Constraint * constraint, Parser * parser);
static Edge * bellmanFord(System * system);
static void relax(Edge * edge);
static Edge * backtrack(Edge * edge);
static void cleanup(System * system);

/*
 * main()
 * - handles file input and output
 * - calls utvpiInterpreter's parseFile() function, which calls initializeSystem() and addConstraint() to build the graph 
 *   representation corresponding to the constraint system in the input file
 * - calls bellmanFord()
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

/*
 * fputEdge() prints the constraint equation corresponding to edge to output
 *
 * edge - pointer to an Edge struct to convert to a constraint equation
 * output - FILE pointer to print the constraint equation to
 */
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
 * bellmanFord() implements BELLMAN-FORD(), as defined in Introduction to Algorithms, Third Edition. For linearly infeasible 
 * systems, a pointer to one Edge within the detected negative cost cycle is returned. For linearly feasible systems, NULL is 
 * returned.
 *
 * system - pointer to the overall System struct containing the graph representation
 */
static Edge * bellmanFord(System * system){
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

/*
 * relax() implements RELAX(), as defined in Introduction to Algorithms, Third Edition. Relaxes a single edge.
 *
 * edge - pointer to the edge to be relaxed
 */
static void relax(Edge * edge){
  if( edge->head->D > edge->tail->D + edge->weight ){
    edge->head->D = edge->tail->D + edge->weight;
    edge->head->L = edge;
  }
}

/*
 * backtrack() implements a single backtrack through the predecessor structure, as defined in Network Flows - Ahuja, Magnanti, 
 * Orlin. This function should only be called in the situation that it is guaranteed that the function will detect a negative cost
 * cycle. Returns a pointer to an Edge within the detected negative cost cycle. 
 * 
 * edge - pointer to an Edge to backtrack from
 */
static Edge * backtrack(Edge * edge){
  while( edge->backtrackSeen == false ){
    edge->backtrackSeen = true;
    edge = edge->tail->L;
  }
  return edge;
}

/*
 * cleanup() frees the graph representation stored within a System struct
 *
 * system - pointer to a System whose graph representation is to be freed
 */
static void cleanup(System * system){
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