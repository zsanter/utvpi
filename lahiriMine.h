#ifndef _LAHIRI_MINE_H
#define _LAHIRI_MINE_H

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include "UTVPIinterpreter.h"

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
  int rho;
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
  bool rootListTraverseSeen;
};

int main(int argc, char * argv[]);
void fputEdge(Edge * edge, FILE * output);
void initializeSystem(void * object, int n, Parser * parser);
void addConstraint(void * object, Constraint * constraint, Parser * parser);
Edge * bellmanFord(System * system);
void relax(Edge * edge);
Edge * backtrack(Edge * edge);
int lahiri(System * Gphi);
void onlySlacklessEdges(System * original, System * subgraph);
void stronglyConnectedComponents(System * system);
void dfsVisit(Vertex * vertex, int * time, int sccNumber);
void transposeSystem(System * original, System * transpose);
int vertexCompareFinishingTimes(const void * vertex1, const void * vertex2);
void johnsonAllPairs(System * Gphi, System * Cstar);
//void copySystem(System * original, System * copy);
void dijkstra(System * system, Vertex * vertex);
void fibHeapInsert(FibHeap * fibHeap, Vertex * vertex);
Vertex * fibHeapExtractMin(FibHeap * fibHeap);
void fibHeapConsolidate(FibHeap * fibHeap);
void fibHeapLink(FibHeap * fibHeap, FibHeapNode * y, FibHeapNode * x);
void fibHeapDecreaseKey(FibHeap * fibHeap, Vertex * vertex);
void fibHeapCut(FibHeap * fibHeap, FibHeapNode * x, FibHeapNode * y);
void fibHeapCascadingCut(FibHeap * fibHeap, FibHeapNode * y);
//void noHeadIndicesHigherThanTailIndeces( System * system );
void freeSystem(System * system);

#endif
