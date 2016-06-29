#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "utvpiInterpreter.h"

typedef struct System System;
typedef struct Vertex Vertex;
typedef struct Edge Edge;

struct System {
  Vertex * graph;
  int vertexCount;
};

struct Vertex {
  int D;
  Edge * first;
};

struct Edge {
  int weight;
  Vertex * tail;
  Vertex * head;
  Edge * next;
};

int main(int argc, char * argv[]);
static void initializeSystem(void * object, int n, Parser * parser);
static void addEdge(void * object, Constraint * constraint, Parser * parser);
static bool bellmanFord(System * system);
static void relax(Edge * edge);
static void cleanup(System * system);

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
  bool solutionExists = bellmanFord(&system);
  if( solutionExists ){
    fputs("The following solution satisfies the given constraints:\n", output);
    for(int i = 0; i < system.vertexCount; i++){
      fprintf(output, "x%i = %i\n", i+1, system.graph[i].D); 
    }
  }
  else{
    fputs("Negative cost cycle detected.\n", output);
  }
  fclose(output);
  cleanup(&system);
  return 0;
}

static void initializeSystem(void * object, int n, Parser * parser){
  System * system = (System *) object;
  system->vertexCount = n;
  system->graph = (Vertex *) malloc( sizeof(Vertex) * system->vertexCount );
  for(int i = 0; i < system->vertexCount; i++){
    system->graph[i].D = 0;
    system->graph[i].first = NULL;
  }
}

static void addEdge(void * object, Constraint * constraint, Parser * parser){
  System * system = (System *) object;
  if( constraint->sign[1] == CONSTRAINT_NONE || constraint->sign[0] == constraint->sign[1] ){
    parseError(parser, "The Bellman-Ford algorithm cannot handle this type of constraint." );
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
  }
}

static bool bellmanFord(System * system){
  for(int i = 0; i < system->vertexCount - 1; i++){
    for(int j = 0; j < system->vertexCount; j++){
      Edge * current = system->graph[j].first;
      while( current != NULL ){
        relax( current );
        current = current->next;
      }
    }
  }
  for(int i = 0; i < system->vertexCount; i++){
    Edge * current = system->graph[i].first;
    while( current != NULL ){
      if( current->head->D > current->tail->D + current->weight ){
        return false;
      }
      current = current->next;
    }
  }
  return true;
}

static void relax(Edge * edge){
  if( edge->head->D > edge->tail->D + edge->weight ){
    edge->head->D = edge->tail->D + edge->weight;
  }
}

static void cleanup(System * system){
  for(int i = 0; i < system->vertexCount; i++){
    Edge * edge = system->graph[i].first;
    while( edge != NULL ){
      Edge * oldEdge = edge;
      edge = edge->next;
      free( oldEdge );
    }
  }
  free( system->graph );
}
