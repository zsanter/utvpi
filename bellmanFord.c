#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "utvpiInterpreter.h"

typedef struct System System;
typedef struct Vertex Vertex;
typedef struct Edge Edge;

struct System {
  Vertex * graph;
  int n;
};

struct Vertex {
  int index;
  int D;
  Edge * L;
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
static void relax(Edge * edge);
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
  }
  else {
    fputs("The following solution satisfies the given constraints:\n", output);
    for(int i = 0; i < system.vertexCount; i++){
      fprintf(output, "x%i = %i\n", i+1, system.graph[i].D); 
    }
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
  system->graph = (Vertex *) malloc( sizeof(Vertex) * system->vertexCount );
  for(int i = 0; i < system->vertexCount; i++){
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
  for(int i = 0; i < system->n - 1; i++){
    for(int j = 0; j < system->n; j++){
      Edge * current = system->graph[j].first;
      while( current != NULL ){
        relax( current );
        current = current->next;
      }
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

static void relax(Edge * edge){
  if( edge->head->D > edge->tail->D + edge->weight ){
    edge->head->D = edge->tail->D + edge->weight;
    edge->head->L = edge;
  }
}

static Edge * backtrack(Edge * edge){
  while( edge->backtrackSeen == false ){
    edge->backtrackSeen = true;
    edge = edge->tail->L;
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
