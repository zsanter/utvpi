#ifndef _SUB_WOJ_LIN_H
#define _SUB_WOJ_LIN_H

#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdbool.h>
#include <time.h>
#include "UTVPIinterpreter.h"
#include "halfint.h"

#define EDGE_TYPE_COUNT 4
typedef enum EdgeType {
  WHITE,
  BLACK,
  GRAY_FORWARD,
  GRAY_REVERSE,
} EdgeType;

#define BACKTRACKING_COEFFICIENT_COUNT 2
typedef enum BacktrackingCoefficient {
  NEG_ONE,
  POS_ONE,
} BacktrackingCoefficient;

typedef struct System System;
typedef struct Vertex Vertex;
typedef struct Edge Edge;
typedef struct EdgeRefList EdgeRefList;

struct System {
  Vertex * graph;
  int vertexCount;
  int n;
  int C; //largestWeight
  Edge * allEdgeFirst;
};

struct Vertex {
  int index;
  Edge * L[EDGE_TYPE_COUNT];
  int D[EDGE_TYPE_COUNT];
  Edge * E[BACKTRACKING_COEFFICIENT_COUNT];
  half_int x;
  Edge * first[EDGE_TYPE_COUNT];
};

struct Edge {
  int weight;
  EdgeType type;
  Vertex * tail;
  Vertex * head;
  Edge * reverse;
  Edge * next;
  Edge * allNext;
};

struct EdgeRefList {
  Edge * edge;
  EdgeRefList * next;
};

int main(int argc, char * argv[]);
void fputEdge(Edge * edge, FILE * output);
EdgeType reverseEdgeType(EdgeType input);
void initializeSystem(void * object, int n, Parser * parser);
void addConstraint(void * object, Constraint * constraint, Parser * parser);
void addEdge(System * system, Constraint * constraint);
void finishSystemCreation(System * system);
EdgeRefList * relaxNetwork(System * system);
void relaxEdge(Edge * e);
EdgeRefList * backtrack(Vertex * x_i, EdgeType t, Edge * e);
void cleanupSystem(System * system);

#endif
