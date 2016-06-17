#ifndef _LAHIRI_MINE_H
#define _LAHIRI_MINE_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <time.h>
#include "UTVPIinterpreter.h"
#include "halfint.h"

#define VERTEX_SIGN_COUNT 2
typedef enum VertexSign {
  POSITIVE,
  NEGATIVE,
} VertexSign;

typedef enum DFScolor {
  WHITE,
  GRAY,
  BLACK,
} DFScolor;

typedef struct System System;
typedef struct Vertex Vertex;
typedef struct Edge Edge;

struct System {
  Vertex * graph[VERTEX_SIGN_COUNT];
  int n;
};

struct Vertex {
  int index;
  VertexSign sign;
  Edge * L;
  int D;
  Edge * first;
  half_int x;
  //int rho;
  DFScolor dfsColor;
  int discoveryTime;
  int finishingTime;
  int sccNumber;
  int h;
};

struct Edge {
  int weight;
  Vertex * head;
  Vertex * tail;
  Edge * next;
  bool backtrackSeen;
};

int main(int argc, char * argv[]);
void fputEdge(Edge * edge, FILE * output);
void initializeSystem(void * object, int n, Parser * parser);
void addConstraint(void * object, Constraint * constraint, Parser * parser);
Edge * bellmanFord(System * system);
void relax(Edge * edge);
Edge * backtrack(Edge * edge);
void cleanup(System * system);

#endif