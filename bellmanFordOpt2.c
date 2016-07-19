#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include "utvpiInterpreter.h"

typedef struct System System;
typedef struct Vertex Vertex;
typedef struct Edge Edge;

struct System {
  Vertex * graph;
  int n;
  int falsePositives;
  int mainLoopIterations;
  int negativeCycleEdgeCount;
};

struct Vertex {
  int index;
  int D;
  Edge * L;
  Edge * cycleOriginator;
  Edge * first;
};

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

static void fputEdge(Edge * edge, FILE * output){
  fprintf(output, "+x%d -x%d <= %d\n", edge->head->index, edge->tail->index, edge->weight);
}

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
    system->graph[i].first = NULL;
  }
}

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
    Edge * current = system->graph[i].first;
    while( current != NULL ){
      if( current->head->D > current->tail->D + current->weight ){
        return backtrack( current );
      }
      current = current->next;
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
    while( edge != NULL ){
      edge->backtrackSeen = false;
      edge = edge->tail->L;
    }
  }
  return edge;
}

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
