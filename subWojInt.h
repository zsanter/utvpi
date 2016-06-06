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

#define BACKTRACKING_INDEX_COUNT 2
typedef enum BacktrackingIndex {
  NEG_ONE,
  POS_ONE,
} BacktrackingIndex;

#define INTEGER_TREE_TYPE_COUNT 2
typedef enum IntegerTreeType {
  DOWN,
  UP,
} IntegerTreeType;

#define INTEGER_TYPE_COUNT 2
typedef enum IntegerType {
  FINAL,
  TEMP,
} IntegerType;

typedef struct System System;
typedef struct Vertex Vertex;
typedef struct Edge Edge;
typedef struct EdgeRefList EdgeRefList;
typedef struct IntegerTree IntegerTree;
typedef struct IntegerTreeVertex IntegerTreeVertex;
//typedef struct IntegerTreeEdge IntegerTreeEdge;

struct System {
  Vertex * graph;
  int vertexCount;
  int n;
  int C; //largestWeight
  Edge * allEdgeFirst;
  //Edge * integerTreeAdditionsFirst;
};

struct Vertex {
  int index;
  Edge * L[EDGE_TYPE_COUNT];
  int D[EDGE_TYPE_COUNT];
  Edge * E[BACKTRACKING_INDEX_COUNT];
  //half_int x; //May go unused
  half_int a; //This seems to stand in for x in the integer feasibility code.
  int Z[INTEGER_TYPE_COUNT];
  Edge * first[EDGE_TYPE_COUNT];
  int edgeCount[EDGE_TYPE_COUNT];
  IntegerTreeVertex * integerTreeVertex[INTEGER_TREE_TYPE_COUNT];
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

struct IntegerTree {
  IntegerTreeVertex * treeRoot;
  IntegerTreeVertex * queueFirst;
  IntegerTreeVertex * queueLast;
  IntegerTreeType type;
  Edge * additionsFirst;
}

struct IntegerTreeVertex {
  IntegerTreeVertex * parent;
  IntegerTreeVertex * nextSibling;
  IntegerTreeVertex * firstChild;
  IntegerTreeVertex * queueNext;
  EdgeRefList * firstEdge;
  Vertex * graphVertex;
};
/*
struct IntegerTreeEdge {
  IntegerTreeVertex * parent;
  IntegerTreeVertex * child;
  IntegerTreeEdge * nextSibling;
  Edge * graphEdge;
}
*/
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
EdgeRefList * produceIntegerSolution(System * system);
bool constraintIsFeasible(Edge * edge);
EdgeRefList * forcedRounding(System * system, IntegerTree * T, Vertex * x_i);
Edge * generateAbsoluteConstraint(System * system, IntegerTree * T, Vertex * x_i, int weight, EdgeType type);
IntegerTree * generateIntegerTree(System * system, IntegerTreeType type);
EdgeRefList * optionalRoundings(System * system, IntegerTree * T);
void checkDependencies(System * system, IntegerTree * T, Vertex * x_i);
void expandIntegerTree(System * system, IntegerTree * T, Vertex * active, Vertex * parent, Edge * newEdge0, Edge * newEdge1, Edge * newEdge2);
EdgeRefList * integerTreeBacktrack(EdgeRefList * list, IntegerTreeVertex * fromVertex, IntegerTreeVertex * toVertex);
EdgeRefList * copyTreeEdgesToList(EdgeRefList * list, IntegerTreeVertex * itv);
void freeIntegerTree(IntegerTree * integerTree);
void freeSystem(System * system);

#endif
