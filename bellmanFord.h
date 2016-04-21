#ifndef _BELLMAN_FORD_H
#define _BELLMAN_FORD_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "UTVPIinterpreter.h"

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
void initializeSystem(void * object, int n, Parser * parser);
void addEdge(void * object, Constraint * constraint, Parser * parser);
bool bellmanFord(System * system);
void relax(Edge * edge);
void cleanup(System * system);

#endif
