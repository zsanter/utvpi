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
typedef struct EdgeRefListNode EdgeRefListNode;
typedef struct IntegerTree IntegerTree;
typedef struct IntegerTreeVertex IntegerTreeVertex;

struct System {
  Vertex * graph;
  int vertexCount;
  int n;
  int C;
  Edge * allEdgeFirst;
  EdgeRefList * infeasibilityProof;
  IntegerTree * T[INTEGER_TREE_TYPE_COUNT];
};

struct Vertex {
  int index;
  Edge * L[EDGE_TYPE_COUNT];
  int D[EDGE_TYPE_COUNT];
  Edge * E[BACKTRACKING_INDEX_COUNT];
  half_int a;
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
  EdgeRefListNode * first;
  EdgeRefListNode * last;
};

struct EdgeRefListNode {
  Edge * edge;
  EdgeRefListNode * next;
};

struct IntegerTree {
  IntegerTreeVertex * treeRoot;
  IntegerTreeVertex * queueNewest;
  IntegerTreeVertex * queueOldest;
  IntegerTreeType type;
  Edge * additionsFirst;
};

struct IntegerTreeVertex {
  IntegerTreeVertex * parent;
  IntegerTreeVertex * nextSibling;
  IntegerTreeVertex * firstChild;
  IntegerTreeVertex * queueNewer;
  EdgeRefList * graphEdges;
  Vertex * graphVertex;
};

int main(int argc, char * argv[]);
EdgeType reverseEdgeType(EdgeType input);
void fputEdge(Edge * edge, FILE * output);
void initializeSystem(void * object, int n, Parser * parser);
void addConstraint(void * object, Constraint * constraint, Parser * parser);
void addEdge(System * system, Constraint * constraint);
void finishSystemCreation(System * system);
int edgeCompare(const void * edge1, const void * edge2);
void removeFromAllEdgeList(System * system, Edge * edge);
bool relaxNetwork(System * system);
void relaxEdge(Edge * e);
bool backtrack(System * system, Vertex * x_i, EdgeType t, Edge * e);
bool produceIntegerSolution(System * system);
bool forcedRounding(System * system, Vertex * x_i);
Edge * generateAbsoluteConstraint(System * system, IntegerTree * T, Vertex * x_i, int weight, EdgeType type);
bool optionalRoundings(System * system);
bool checkAllConstraints(System * system, Vertex * toVertex, IntegerType integerType, IntegerTreeType integerTreeType);
void checkDependencies(System * system, IntegerTree * T, Vertex * x_i, IntegerType integerType);
IntegerTree * generateIntegerTree(System * system, IntegerTreeType type);
Vertex * pollIntegerTreeQueue(IntegerTree * tree);
void expandIntegerTree(System * system, IntegerTree * T, Vertex * active, Vertex * parent, Edge * edge0, Edge * edge1, Edge * edge2);
void integerTreeBacktrack(EdgeRefList * list, IntegerTreeVertex * fromVertex, IntegerTreeVertex * toVertex, bool includeToVertex);
void copyTreeEdgesToList(EdgeRefList * list, IntegerTreeVertex * itv);
void freeIntegerTree(IntegerTree * tree);
EdgeRefList * generateEdgeRefList();
void addEdgeToEdgeRefList(EdgeRefList * erl, Edge * edge);
void freeEdgeRefList(EdgeRefList * erl);
void freeSystem(System * system);


#endif
