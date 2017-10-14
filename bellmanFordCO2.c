/*
 * bellmanFordCO2.c
 * The Bellman-Ford difference constraint system solver, with second cycle-originator arrangement
 *
 * Cycle-originators are much less reliable in this version. In systems where many edges could be part of a number of different 
 * negative cost cycles, as edges are relaxed, rather than passing cycle-originators through fully-formed and unchanging negative-
 * cost cycles, predecessors for each edge continue to change, making this cycle-originator setup ineffective. It is only 
 * guaranteed to work if there is only one negative cost cycle within a graph, likely not a particularly common occurrence in 
 * real-world applications.
 *
 * Call with [executable] [input file] {output file}
 * [input file] must be properly formatted to be read by utvpiInterpreter.h
 * {output file} will contain a solution, if one exists. If the system is not feasible, a proof of infeasibility - a negative cost
 *   cycle - will be output. If {output file} is not specified, output will be to stdout.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include "utvpiInterpreter.h"

typedef struct System System;
typedef struct Vertex Vertex;
typedef struct Edge Edge;

/*
 * The System struct contains the graph representation used by the system solver, along with other information about the system.
 *
 * graph - pointer to an array of Vertex structs, defining the structure of the graph representation
 * n - the number of variables represented within the system, equal to the number of Vertex structs stored within the graph array
 * falsePositives - number of false positives thrown by the cycle-originator negative cycle detection mechanism
 * mainLoopIterations - number of iterations of the main relaxation loop within bellmanFord()
 * negativeCycleEdgeCount - number of edges within a detected negative cost cycle
 */
struct System {
  Vertex * graph;
  int n;
  int falsePositives;
  int mainLoopIterations;
  int negativeCycleEdgeCount;
};

/*
 * The Vertex struct contains all information about a variable, represented by a vertex within the graph.
 *
 * index - index of the variable
 * D - distance label
 * L - pointer to the predecessor Edge between this Vertex and the Vertex's predecessor Vertex
 * cycleOriginator - pointer to an edge serving as a cycle-originator. Cycle-originators are used to detect negative cost cycles 
 *   before the total number of relaxation loop iterations specified by the algorithm run
 * first - pointer to the first Edge whose tail Vertex is this Vertex. The remainder of such edges are connected together in a 
 *   singly-linked list
 */
struct Vertex {
  int index;
  int D;
  Edge * L;
  Edge * cycleOriginator;
  Edge * first;
};

/*
 * The Edge struct contains all information about a specific constraint, represented by an edge within the graph.
 *
 * weight - the weight of the edge, corresponding to the defining constant of the constraint the edge represents
 * tail - pointer to the Edge's tail Vertex
 * head - pointer to the Edge's head Vertex
 * next - pointer to the next Edge in the singly-linked list of Edges with the same tail Vertex as one another
 * backtrackSeen - boolean flag indicating whether or not this Edge has been encountered by the backtrack() function as it 
 *   attempts to detect a negative cost cycle
 */
struct Edge {
  int weight;
  Vertex * tail;
  Vertex * head;
  Edge * next;
  bool backtrackSeen;
};

int main(int argc, char * argv[]);
static void fputEdge(Edge * edge, FILE * output);
static void initializeSystem(void * object, int n, Parser * parser);
static void addEdge(void * object, Constraint * constraint, Parser * parser);
static Edge * bellmanFord(System * system);
static Edge * backtrack(Edge * edge);
static void freeSystem(System * system);

/*
 * main()
 * - handles file input and output
 * - calls utvpiInterpreter's parseFile() function, which calls initializeSystem() and addConstraint() to build the graph 
 *   representation corresponding to the constraint system in the input file
 * - calls bellmanFord()
 */
int main(int argc, char * argv[]){
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
  bool parseSuccessful = parseFile(input, &system, initializeSystem, addEdge);
  fclose(input);
  if( !parseSuccessful ){
    fprintf( stderr, "Parsing of \"%s\" failed.\n", argv[1] );
    exit(1);
  }
  Edge * negativeCycle = bellmanFord(&system);
  if( negativeCycle != NULL ){
    fputs("The following negative cost cycle was detected:\n", output);
    Edge * edge = negativeCycle;
    while( edge->backtrackSeen == true ){
      fputEdge( edge, output );
      edge->backtrackSeen = false;
      system.negativeCycleEdgeCount++;
      edge = edge->tail->L;
    }
    fprintf( output, "\n%d false positives\n", system.falsePositives );
    fprintf( output, "%d main loop iterations\n", system.mainLoopIterations );
    fprintf( output, "%d negative cycle edges\n", system.negativeCycleEdgeCount );
  }
  else {
    fputs("The following solution satisfies the given constraints:\n", output);
    for(int i = 0; i < system.n; i++){
      fprintf(output, "x%i = %i\n", i+1, system.graph[i].D); 
    }
    fprintf( output, "\n%d false positives\n", system.falsePositives );
    fprintf( output, "%d main loop iterations\n", system.mainLoopIterations );
  }
  fclose(output);
  freeSystem(&system);
  return 0;
}

/*
 * fputEdge() prints the constraint equation corresponding to edge to output
 *
 * edge - pointer to an Edge struct to convert to a constraint equation
 * output - FILE pointer to print the constraint equation to
 */
static void fputEdge(Edge * edge, FILE * output){
  fprintf(output, "+x%d -x%d <= %d\n", edge->head->index, edge->tail->index, edge->weight);
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
  system->graph = (Vertex *) malloc( sizeof(Vertex) * system->n );
  system->falsePositives = 0;
  system->mainLoopIterations = INT_MAX;
  system->negativeCycleEdgeCount = 0;
  for(int i = 0; i < system->n; i++){
    system->graph[i].index = i+1;
    system->graph[i].D = 0;
    system->graph[i].L = NULL;
    system->graph[i].cycleOriginator = NULL;
    system->graph[i].first = NULL;
  }
}

/*
 * addEdge() adds a constraint to the graph representation held by the System struct.
 *
 * object - a void pointer pointing to an already-initialized System struct
 * constraint - pointer to a Constraint struct describing a constraint
 * parser - pointer to the Parser struct that utvpiInterpreter uses during the input file parsing process, so that parseError() 
 *   can be called, if need be
 */
static void addEdge(void * object, Constraint * constraint, Parser * parser){
  System * system = (System *) object;
  if( constraint->sign[1] == CONSTRAINT_NONE || constraint->sign[0] == constraint->sign[1] ){
    parseError(parser, "The Bellman-Ford algorithm can only handle difference constraints." );
  }
  else{
    int tailIndex, headIndex;
    if( constraint->sign[0] == CONSTRAINT_MINUS ) {
      tailIndex = constraint->index[0] - 1;
      headIndex = constraint->index[1] - 1;
    }
    else{
      tailIndex = constraint->index[1] - 1;
      headIndex = constraint->index[0] - 1;
    }
    Edge * newEdge = (Edge *) malloc( sizeof(Edge) );
    newEdge->weight = constraint->weight;
    newEdge->tail = &system->graph[tailIndex];
    newEdge->head = &system->graph[headIndex];
    newEdge->next = system->graph[tailIndex].first;
    system->graph[tailIndex].first = newEdge;
    newEdge->backtrackSeen = false;
  }
}

/*
 * bellmanFord() implements BELLMAN-FORD() and RELAX(), as defined in Introduction to Algorithms, Third Edition, with the 
 * modification of adding cycle-originators that are reinitialized every time a Vertex's predecessor Edge changes and which are 
 * passed through the predecessor structure each time a distance label changes without the predecessor edge changing. For linearly
 * infeasible systems, a pointer to one Edge within the detected negative cost cycle is returned. For linearly feasible systems,
 * NULL is returned.
 *
 * system - pointer to the overall System struct containing the graph representation
 */
static Edge * bellmanFord(System * system){
  bool anyChange = true;
  for(int i = 1; ( i <= system->n - 1 ) && anyChange; i++){
    anyChange = false;
    for(int j = 0; j < system->n; j++){
      Edge * edge = system->graph[j].first;
      while( edge != NULL ){
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
    if( !anyChange ){
      system->mainLoopIterations = i;
    }
  }
  for(int i = 0; i < system->n; i++){
    Edge * edge = system->graph[i].first;
    while( edge != NULL ){
      if( edge->head->D > edge->tail->D + edge->weight ){
        return backtrack( edge );
      }
      edge = edge->next;
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
 * freeSystem() frees the graph representation stored within a System struct
 *
 * system - pointer to a System whose graph representation is to be freed
 */
static void freeSystem(System * system){
  for(int i = 0; i < system->n; i++){
    Edge * edge = system->graph[i].first;
    while( edge != NULL ){
      Edge * oldEdge = edge;
      edge = edge->next;
      free( oldEdge );
    }
  }
  free( system->graph );
}
