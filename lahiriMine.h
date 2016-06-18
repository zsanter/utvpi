#ifndef _LAHIRI_MINE_H
#define _LAHIRI_MINE_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <math.h>
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
typedef struct FibHeap FibHeap;
typedef struct FibHeapNode FibHeapNode;

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
  bool dijkstraFinalized;
  FibHeapNode * fibHeapNode;
};

struct Edge {
  int weight;
  Vertex * head;
  Vertex * tail;
  Edge * next;
  bool backtrackSeen;
};

struct FibHeap {
  FibHeapNode * min;
  int n;
};

struct FibHeapNode {
  FibHeapNode * parent;
  FibHeapNode * left;
  FibHeapNode * right;
  FibHeapNode * child;
  Vertex * vertex;
  int degree;
  bool mark;
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