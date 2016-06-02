#ifndef _SUB_WOJ_INT_H
#define _SUB_WOJ_INT_H

#include <stdio.h>
#include <limits.h>
#include <stdlib.h>
#include <stdbool.h>
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
typedef struct IntegerTreeVertex IntegerTreeVertex;
typedef struct IntegerTreeEdge IntegerTreeEdge;

struct System {
  Vertex * graph;
  int vertexCount;
  int n;
  int C; //largestWeight
  Edge * allEdgeFirst;
  IntegerTreeEdge * x0FirstChild;
  IntegerTreeVertex * integerQueueFirst;
  IntegerTreeVertex * integerQueueLast;
};

struct Vertex {
  int index;
  Edge * L[EDGE_TYPE_COUNT];
  int D[EDGE_TYPE_COUNT];
  Edge * E[BACKTRACKING_COEFFICIENT_COUNT];
  //half_int x; //May go unused
  half_int a; //This seems to stand in for x in the integer feasibility code.
  int Z;
  //int Z_T;
  Edge * first[EDGE_TYPE_COUNT];
  int edgeCount[EDGE_TYPE_COUNT];
};

struct Edge {
  int weight;
  EdgeType type;
  Vertex * tail;
  Vertex * head;
  Edge * reverse;
  Edge * next;
  Edge * allNext;
  Edge * allPrev;
  bool inAllEdgeList;
};

struct EdgeRefList {
  Edge * edge;
  EdgeRefList * next;
};

struct IntegerTreeVertex {
  IntegerTreeEdge * parent;
  IntegerTreeEdge * firstChild;
  IntegerTreeVertex * queueNext;
  Vertex * graphVertex;
};

struct IntegerTreeEdge {
  IntegerTreeVertex * parent;
  IntegerTreeVertex * child;
  IntegerTreeEdge * nextSibling;
  Edge * graphEdge;
}

int main(int argc, char * argv[]);
void fputEdge(Edge * edge, FILE * output);
EdgeType reverseEdgeType(EdgeType input);
void initializeSystem(void * object, int n, Parser * parser);
void addConstraint(void * object, Constraint * constraint, Parser * parser);
void addEdge(Systen * system, Constraint * constraint);
void finishSystemCreation(System * system);
int edgeCompare(const void * edge1, const void * edge2);
void removeFromAllEdgeList(System * system, Edge * edge);
EdgeRefList * relaxNetwork(System * system);
void relaxEdge(Edge * e);
EdgeRefList * backtrack(Vertex * x_i, EdgeType t, Edge * e);
void produceIntegerSolution(System * system);
void forcedRounding(Vertex * x_i);

#endif
