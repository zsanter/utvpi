#include "subWojInt.h"

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
  bool parseSuccessful = parseFile(input, &system, initializeSystem, addConstraint);
  fclose(input);
  if( !parseSuccessful ){
    fprintf( stderr, "Parsing of \"%s\" failed.\n", argv[1] );
    exit(1);
  }
  finishSystemCreation(&system);
  EdgeRefList * R = relaxNetwork(&system);
  if( R != NULL ){
    fputs("The following negative gray cycle was detected:\n", output);
    while( R != NULL ){
      fputEdge(R->edge, output);
      EdgeRefList * oldR = R;
      R = R->next;
      free(oldR);
    }
  }
  else{
    fputs("Linear solution:\n", output);
    for(int i = 1; i < system.vertexCount; i++){
      system.graph[i].a = intDivBy2ToHalfInt( system.graph[i].D[WHITE] - system.graph[i].D[BLACK] );
      fprintf(output, "x%i = %.1f\n", i, halfIntToDouble( system.graph[i].a ) );
    }
    EdgeRefList * O = produceIntegerSolution(&system);
    if(O != NULL){
      fputs("\nSystem is not integrally feasible.\nProof:\n", output);
      while( O != NULL ){
        fputEdge(O->edge, output);
        EdgeRefList * oldO = O;
        O = O->next;
        free(oldO);
      }
    }
    else {
      fputs("\nIntegral solution:\n", output);
      for(int i = 1; i < system.vertexCount; i++){
        fprintf(output, "x%i = %i\n", i, system.graph[i].Z );
      }
    }
  }
  
  fclose(output);
  return 0;
}

EdgeType reverseEdgeType(EdgeType input){ 
  EdgeType output;
  switch(input){
  case GRAY_FORWARD:
    output = GRAY_REVERSE;
    break;
  case GRAY_REVERSE:
    output = GRAY_FORWARD;
    break;
  default:
    output = input;
  }
  return output;
}

void fputEdge(Edge * edge, FILE * output){
  char sign[2];
  switch(edge->type){
  case WHITE:
    sign[0] = '+';
    sign[1] = '+';
    break;
  case BLACK:
    sign[0] = '-';
    sign[1] = '-';
    break;
  case GRAY_FORWARD:
    sign[0] = '-';
    sign[1] = '+';
    break;
  case GRAY_REVERSE:
    sign[0] = '+';
    sign[1] = '-';
    break;
  }
  if( edge->tail->index != 0 ){
    fprintf(output, "%cx%i ", sign[0], edge->tail->index);
  }
  if( edge->head->index != 0 ){
    fprintf(output, "%cx%i ", sign[1], edge->head->index);
  }
  fprintf(output, "<= %i\n", edge->weight);
}

void initializeSystem(void * object, int n, Parser * parser){
  System * system = (System *) object;
  system->vertexCount = n + 1;
  system->graph = (Vertex *) malloc( sizeof(Vertex) * system->vertexCount );
  system->n = n;
  system->C = 0;
  system->allEdgeFirst = NULL;
  system->graph[0].index = 0;
  system->graph[0].E[NEG_ONE] = NULL;
  system->graph[0].E[POS_ONE] = NULL;
  system->graph[0].a = 0;
  //system->graph[0].x = 0;
  system->graph[0].Z = 0;
  //system->graph[0].Z_T = 0;
  system->graph[0].integerTreeVertex = NULL;
  for(EdgeType i = WHITE; i <= GRAY_REVERSE; i++){
    system->graph[0].L[i] = NULL;
    system->graph[0].D[i] = 0;
    system->graph[0].first[i] = NULL
    system->graph[0].edgeCount[i] = 0;
  }
  for(int i = 1; i < system->vertexCount; i++){
    system->graph[i].index = i;
    system->graph[i].E[NEG_ONE] = NULL;
    system->graph[i].E[POS_ONE] = NULL;
    system->graph[i].a = 0;
    //system->graph[i].x = 0;
    system->graph[i].Z = INT_MAX;
    //system->graph[i].Z_T = 0;
    system->graph[i].integerTreeVertex = NULL;
    for(EdgeType j = WHITE; j <= GRAY_REVERSE; j++){
      system->graph[i].L[j] = NULL;
      system->graph[i].D[j] = INT_MAX;
      system->graph[i].first[j] = NULL;
      system->graph[i].edgeCount[j] = 0;
    }
  }
}

void addConstraint(void * object, Constraint * constraint, Parser * parser){
  if( constraint->index[0] == constraint->index[1] ){
    parseError(parser, "Both indices in a constraint can not be the same." );
  }
  else{
    System * system = (System *) object;
    if( abs( constraint->weight ) > system->C ){
      system->C = abs( constraint->weight );
    }
    if( constraint->sign[1] == NONE ){
      constraint->index[1] = 0;
      constraint->sign[1] = PLUS;
      addEdge( system, constraint );
      constraint->sign[1] = MINUS;
      addEdge( system, constraint );
      //Calls to addEdge place the new edges at first[edgeType] for both vertices
      if( constraint->sign[0] == PLUS ){
        if( constraint->weight < system->graph[ constraint->index[0] ].D[WHITE] ){
          system->graph[ constraint->index[0] ].D[WHITE] = constraint->weight;
          system->graph[ constraint->index[0] ].L[WHITE] = system->graph[0].first[WHITE];
          system->graph[ constraint->index[0] ].D[GRAY_FORWARD] = constraint->weight;
          system->graph[ constraint->index[0] ].L[GRAY_FORWARD] = system->graph[0].first[GRAY_FORWARD]; //.first[GRAY_REVERSE]; ? No.
        }
      }
      else{
        if( constraint->weight < system->graph[ constraint->index[0] ].D[BLACK] ){
          system->graph[ constraint->index[0] ].D[BLACK] = constraint->weight;
          system->graph[ constraint->index[0] ].L[BLACK] = system->graph[0].first[BLACK];
          system->graph[ constraint->index[0] ].D[GRAY_REVERSE] = constraint->weight;
          system->graph[ constraint->index[0] ].L[GRAY_REVERSE] = system->graph[0].first[GRAY_REVERSE]; //.first[GRAY_FORWARD]; ? No.
        }
      }
    }
    else{
      addEdge( system, constraint );
    }
  }
}

void addEdge(Systen * system, Constraint * constraint){
  Edge * newEdges[2];
  newEdges[0] = (Edge *) malloc( sizeof(Edge) );
  newEdges[1] = (Edge *) malloc( sizeof(Edge) );
  newEdges[0]->weight = constraint->weight;
  newEdges[0]->reverse = newEdges[1];
  newEdges[0]->next = NULL;
  newEdges[1]->weight = constraint->weight;
  newEdges[1]->reverse = newEdges[0];
  newEdges[1]->next = NULL;
  if( constraint->sign[0] == constraint->sign[1] ){
    EdgeType edgeType;
    if( constraint->sign[0] == PLUS ){
      edgeType = WHITE;
    }
    else{
      edgeType = BLACK;
    }
    newEdges[0]->type = edgeType;
    newEdges[0]->tail = &system->graph[ constraint->index[0] ];
    newEdges[0]->head = &system->graph[ constraint->index[1] ];
    newEdges[1]->type = edgeType;
    newEdges[1]->tail = &system->graph[ constraint->index[1] ];
    newEdges[1]->head = &system->graph[ constraint->index[0] ];
  }
  else{
    int negativeIndex, positiveIndex;
    if( constraint->sign[0] == PLUS ){
      positiveIndex = 0;
      negativeIndex = 1;
    }
    else{
      positiveIndex = 1;
      negativeIndex = 0;
    }
    newEdges[0]->type = GRAY_FORWARD;
    newEdges[0]->tail = &system->graph[ constraint->index[ negativeIndex ] ];
    newEdges[0]->head = &system->graph[ constraint->index[ positiveIndex ] ];
    newEdges[1]->type = GRAY_REVERSE;
    newEdges[1]->tail = &system->graph[ constraint->index[ positiveIndex ] ];
    newEdges[1]->head = &system->graph[ constraint->index[ negativeIndex ] ];
  }
  newEdges[0]->next = newEdges[0]->tail->first[ newEdges[0]->type ];
  newEdges[0]->tail->first[ newEdges[0]->type ] = newEdges[0];
  newEdges[0]->tail->edgeCount[ newEdges[0]->type ]++;
  newEdges[0]->allNext = system->allEdgeFirst;
  newEdges[0]->allPrev = NULL;
  newEdges[0]->inAllEdgeList = true;
  if( system->allEdgeFirst != NULL ){
    system->allEdgeFirst->allPrev = newEdges[0];
  }
  system->allEdgeFirst = newEdges[0];
  newEdges[1]->next = newEdges[1]->tail->first[ newEdges[1]->type ];
  newEdges[1]->tail->first[ newEdges[1]->type ] = newEdges[1];
  newEdges[1]->tail->edgeCount[ newEdges[1]->type ]++;
  newEdges[1]->allNext = NULL;
  newEdges[1]->allPrev = NULL;
  newEdges[1]->inAllEdgeList = false;
}

void finishSystemCreation(System * system){
  int sourceEdgeWeights = (2 * system->n + 1) * system->C;
  for(int i = 0; i < system->vertexCount; i++){
    for(EdgeType j = WHITE; j <= GRAY_REVERSE; j++){
      if( system->graph[i].D[j] == INT_MAX ){
        system->graph[i].D[j] = sourceEdgeWeights;
      }
      if( system->graph[i].edgeCount[j] > 1 ){
        int edgeSortArrayLength = system->graph[i].edgeCount[j];
        Edge * edgeSortArray[ edgeSortArrayLength ];
        Edge * edge = system->graph[i].first[j];
        for(int k = 0; k < edgeSortArrayLength; k++){
          edgeSortArray[k] = edge;
          edge = edge->next;
        }
        qsort(edgeSortArray, edgeSortArrayLength, sizeof(Edge *), edgeCompare);
        system->graph[i].first[j] = edgeSortArray[0];
        Edge * beforePrior = NULL;
        Edge * prior = edgeSortArray[0];
        for(int k = 1; k < edgeSortArrayLength; k++){
          if( prior->head == edgeSortArray[k]->head ){
            if( prior->weight <= edgeSortArray[k]->weight ){
              //puts("free 2 before");
              //fputEdge( edgeSortArray[k], stdout);
              removeFromAllEdgeList( system, edgeSortArray[k] );
              free( edgeSortArray[k] );
              //puts("free 2 after");
            }
            else{
              if( system->graph[i].first[j] == prior ){
                system->graph[i].first[j] = edgeSortArray[k];
              }
              else {
                beforePrior->next = edgeSortArray[k];
              }
              //puts("free 3 before");
              removeFromAllEdgeList( system, prior );
              free( prior );
              //puts("free 3 after");
              prior = edgeSortArray[k];
            }
            system->graph[i].edgeCount[j]--;
          }
          else{
            prior->next = edgeSortArray[k];
            beforePrior = prior;
            prior = edgeSortArray[k];
          }
        }
        prior->next = NULL;
      }
      /*
      Edge * edge = system->graph[i].first[j];
      while(edge != NULL){
        fputEdge(edge, stdout);
        edge = edge->next;
      }
      */
    }
  }
}

int edgeCompare(const void * edge1, const void * edge2){
  return (*(Edge **)edge1)->head->index - (*(Edge **)edge2)->head->index;
}

void removeFromAllEdgeList(System * system, Edge * edge){
  if( edge->inAllEdgeList == true ){
    if( edge->allNext != NULL ){
      edge->allNext->allPrev = edge->allPrev;
    }
    if( edge->allPrev != NULL ){
      edge->allPrev->allNext = edge->allNext;
    }
    else{
      system->allEdgeFirst = edge->allNext;
    }
  }
}

EdgeRefList * relaxNetwork(System * system){
  //Lines 3-6 of algorithm implemented in finishSystemCreation().
  for(int r = 1; r <= 2 * system->n; r++){
    Edge * e = system->allEdgeFirst;
    while( e != NULL ){
      relaxEdge(e);
      e = e->allNext;
    }
  }
  for(int i = 0; i < system->vertexCount; i++){
    Edge * e = system->graph[i].first[WHITE];
    while( e != NULL ){
      if( e->head->D[GRAY_REVERSE] + e->weight < e->tail->D[WHITE] ){
        return backtrack(e->head, GRAY_REVERSE, e->reverse);
      }
      if( e->head->D[BLACK] + e->weight < e->tail->D[GRAY_FORWARD] ){
        return backtrack(e->head, BLACK, e->reverse);
      }
      if( e->tail->D[GRAY_REVERSE] + e->weight < e->head->D[WHITE] ){
        return backtrack(e->tail, GRAY_REVERSE, e);
      }
      if( e->tail->D[BLACK] + e->weight < e->head->D[GRAY_FORWARD] ){
        return backtrack(e->tail, BLACK, e);
      }
      e = e->next;
    }
  }
  for(int i = 0; i < system->vertexCount; i++){
    Edge * e = system->graph[i].first[BLACK];
    while( e != NULL ){
      if( e->head->D[GRAY_FORWARD] + e->weight < e->tail->D[BLACK] ){
        return backtrack(e->head, GRAY_FORWARD, e->reverse);
      }
      if( e->head->D[WHITE] + e->weight < e->tail->D[GRAY_REVERSE] ){
        return backtrack(e->head, WHITE, e->reverse);
      }
      if( e->tail->D[GRAY_FORWARD] + e->weight < e->head->D[BLACK] ){
        return backtrack(e->tail, GRAY_FORWARD, e);
      }
      if( e->tail->D[WHITE] + e->weight < e->head->D[GRAY_REVERSE] ){
        return backtrack(e->tail, WHITE, e);
      }
      e = e->next;
    }
  }
  for(int i = 0; i < system->vertexCount; i++){
    Edge * e = system->graph[i].first[GRAY_FORWARD];
    while( e != NULL ){
      if( e->head->D[GRAY_REVERSE] + e->weight < e->tail->D[GRAY_REVERSE] ){
        return backtrack(e->head, GRAY_REVERSE, e->reverse);
      }
      if( e->head->D[BLACK] + e->weight < e->tail->D[BLACK] ){
        return backtrack(e->head, BLACK, e->reverse);
      }
      if( e->tail->D[GRAY_FORWARD] + e->weight < e->head->D[GRAY_FORWARD] ){
        return backtrack(e->tail, GRAY_FORWARD, e);
      }
      if( e->tail->D[WHITE] + e->weight < e->head->D[WHITE] ){
        return backtrack(e->tail, WHITE, e);
      }
      e = e->next;
    }
  }
  for(int i = 0; i < system->vertexCount; i++){
    Edge * e = system->graph[i].first[GRAY_REVERSE];
    while( e != NULL ){
      if( e->head->D[GRAY_FORWARD] + e->weight < e->tail->D[GRAY_FORWARD] ){
        return backtrack(e->head, GRAY_FORWARD, e->reverse);
      }
      if( e->head->D[WHITE] + e->weight < e->tail->D[WHITE] ){
        return backtrack(e->head, WHITE, e->reverse);
      }
      if( e->tail->D[GRAY_REVERSE] + e->weight < e->head->D[GRAY_REVERSE] ){
        return backtrack(e->tail, GRAY_REVERSE, e);
      }
      if( e->tail->D[BLACK] + e->weight < e->head->D[BLACK] ){
        return backtrack(e->tail, BLACK, e);
      }
      e = e->next;
    }
  }
  return NULL;
}

void relaxEdge(Edge * e){
  switch(e->type){
  case WHITE:
    if( e->head->D[GRAY_REVERSE] + e->weight < e->tail->D[WHITE] ){
      e->tail->D[WHITE] = e->head->D[GRAY_REVERSE] + e->weight;
      e->tail->L[WHITE] = e->reverse;
    }
    if( e->head->D[BLACK] + e->weight < e->tail->D[GRAY_FORWARD] ){
      e->tail->D[GRAY_FORWARD] = e->head->D[BLACK] + e->weight;
      e->tail->L[GRAY_FORWARD] = e->reverse;
    }
    if( e->tail->D[GRAY_REVERSE] + e->weight < e->head->D[WHITE] ){
      e->head->D[WHITE] = e->tail->D[GRAY_REVERSE] + e->weight;
      e->head->L[WHITE] = e;
    }
    if( e->tail->D[BLACK] + e->weight < e->head->D[GRAY_FORWARD] ){
      e->head->D[GRAY_FORWARD] = e->tail->D[BLACK] + e->weight;
      e->head->L[GRAY_FORWARD] = e;
    }
    break;
  case BLACK:
    if( e->head->D[GRAY_FORWARD] + e->weight < e->tail->D[BLACK] ){
      e->tail->D[BLACK] = e->head->D[GRAY_FORWARD] + e->weight;
      e->tail->L[BLACK] = e->reverse;
    }
    if( e->head->D[WHITE] + e->weight < e->tail->D[GRAY_REVERSE] ){
      e->tail->D[GRAY_REVERSE] = e->head->D[WHITE] + e->weight;
      e->tail->L[GRAY_REVERSE] = e->reverse;
    }
    if( e->tail->D[GRAY_FORWARD] + e->weight < e->head->D[BLACK] ){
      e->head->D[BLACK] = e->tail->D[GRAY_FORWARD] + e->weight;
      e->head->L[BLACK] = e;
    }
    if( e->tail->D[WHITE] + e->weight < e->head->D[GRAY_REVERSE] ){
      e->head->D[GRAY_REVERSE] = e->tail->D[WHITE] + e->weight;
      e->head->L[GRAY_REVERSE] = e;
    }
    break;
  case GRAY_FORWARD:
    if( e->head->D[GRAY_REVERSE] + e->weight < e->tail->D[GRAY_REVERSE] ){
      e->tail->D[GRAY_REVERSE] = e->head->D[GRAY_REVERSE] + e->weight;
      e->tail->L[GRAY_REVERSE] = e->reverse;
    }
    if( e->head->D[BLACK] + e->weight < e->tail->D[BLACK] ){
      e->tail->D[BLACK] = e->head->D[BLACK] + e->weight;
      e->tail->L[BLACK] = e->reverse;
    }
    if( e->tail->D[GRAY_FORWARD] + e->weight < e->head->D[GRAY_FORWARD] ){
      e->head->D[GRAY_FORWARD] = e->tail->D[GRAY_FORWARD] + e->weight;
      e->head->L[GRAY_FORWARD] = e;
    }
    if( e->tail->D[WHITE] + e->weight < e->head->D[WHITE] ){
      e->head->D[WHITE] = e->tail->D[WHITE] + e->weight;
      e->head->L[WHITE] = e;
    }
    break;
  case GRAY_REVERSE:
    if( e->head->D[GRAY_FORWARD] + e->weight < e->tail->D[GRAY_FORWARD] ){
      e->tail->D[GRAY_FORWARD] = e->head->D[GRAY_FORWARD] + e->weight;
      e->tail->L[GRAY_FORWARD] = e->reverse;
    }
    if( e->head->D[WHITE] + e->weight < e->tail->D[WHITE] ){
      e->tail->D[WHITE] = e->head->D[WHITE] + e->weight;
      e->tail->L[WHITE] = e->reverse;
    }
    if( e->tail->D[GRAY_REVERSE] + e->weight < e->head->D[GRAY_REVERSE] ){
      e->head->D[GRAY_REVERSE] = e->tail->D[GRAY_REVERSE] + e->weight;
      e->head->L[GRAY_REVERSE] = e;
    }
    if( e->tail->D[BLACK] + e->weight < e->head->D[BLACK] ){
      e->head->D[BLACK] = e->tail->D[BLACK] + e->weight;
      e->head->L[BLACK] = e;
    }
  }
}

EdgeRefList * backtrack(Vertex * x_i, EdgeType t, Edge * e){
  Vertex * x_c = x_i;
  Edge * e_c = e;
  EdgeType t_c = t;
  BacktrackingCoefficient a_c;
  switch(t_c){
  case WHITE:
  case GRAY_FORWARD:
    a_c = POS_ONE;
    break;
  case BLACK:
  case GRAY_REVERSE:
    a_c = NEG_ONE;
    break;
  }
  while( x_c->E[a_c] == NULL ){
    x_c->E[a_c] = e_c;
    e_c = x_c->L[t_c];
    switch(t_c){
    case WHITE:
      switch(e_c->type){
      case WHITE:
        t_c = GRAY_REVERSE;
        a_c = NEG_ONE;
        break;
      case GRAY_FORWARD:
        t_c = WHITE;
        a_c = POS_ONE;
        break;
      }
      break;
    case BLACK:
      switch(e_c->type){
      case BLACK:
        t_c = GRAY_FORWARD;
        a_c = POS_ONE;
        break;
      case GRAY_REVERSE:
        t_c = BLACK;
        a_c = NEG_ONE;
        break;
      }
      break;
    case GRAY_FORWARD:
      switch(e_c->type){
      case GRAY_FORWARD:
        t_c = GRAY_FORWARD;
        a_c = POS_ONE;
        break;
      case WHITE:
        t_c = BLACK;
        a_c = NEG_ONE;
        break;
      }
      break;
    case GRAY_REVERSE:
      switch(e_c->type){
      case GRAY_REVERSE:
        t_c = GRAY_REVERSE;
        a_c = NEG_ONE;
        break;
      case BLACK:
        t_c = WHITE;
        a_c = POS_ONE;
        break;
      }
      break;
    }
    x_c = e_c->tail;
  }
  Vertex * x_f = x_c;
  BacktrackingCoefficient a_f = a_c;
  EdgeRefList * R = (EdgeRefList *) malloc( sizeof(EdgeRefList) );
  R->edge = e_c;
  R->next = NULL;
  EdgeRefList * R_last = R;
  switch(e_c->type){
  case WHITE:
  case GRAY_FORWARD:
    a_c = POS_ONE;
    break;
  case BLACK:
  case GRAY_REVERSE:
    a_c = NEG_ONE;
    break;
  }
  x_c = e_c->head;
  while(a_c != a_f || x_c != x_f){
    EdgeRefList * newEdgeRef = (EdgeRefList *) malloc( sizeof(EdgeRefList) );
    newEdgeRef->edge = x_c->E[a_c];
    newEdgeRef->next = NULL;
    R_last->next = newEdgeRef;
    R_last = newEdgeRef;
    Vertex * x_c_next = x_c->E[a_c]->head;
    switch( x_c->E[a_c]->type ){
    case WHITE:
    case GRAY_FORWARD:
      a_c = POS_ONE;
      break;
    case BLACK:
    case GRAY_REVERSE:
      a_c = NEG_ONE;
      break;
    }
    x_c = x_c_next;
  }
  return R;
}

EdgeRefList * produceIntegerSolution(System * system){
  
  IntegerTreeVertex * x0Root = (IntegerTreeVertex *) malloc( sizeof(IntegerTreeVertex) );
  x0Root->parent = NULL;
  x0Root->nextSibling = NULL;
  x0Root->firstChild = NULL;
  x0Root->queueNext = NULL;
  x0Root->firstEdge = NULL;
  x0Root->graphVertex = &system->graph[0];
  system->graph[0].integerTreeVertex = x0Root;
  
  IntegerTree T;
  T.treeRoot = x0Root;
  T.queueFirst = NULL;
  T.queueLast = NULL;
  T.type = X0;
  
  for(int i = 1; i < system->vertexCount; i++){
    if( halfIntIsIntegral( system->graph[i].a ) ){
      system->graph[i].Z = halfIntToInt( system->graph[i].a );
    }
    else{
      EdgeRefList * infeasibilityPrrof = forcedRounding( system, &system->graph[i], &T );
      if( infeasibilityProof != NULL ){
        return infeasibilityProof;
      }
    }
  }
  
  while( T.queueLast != NULL ){
    checkDependencies(system, T.queueLast->graphVertex, &T);
    T.queueLast = T.queueLast->queueNext;
  }
  
  return NULL;
}

EdgeRefList * forcedRounding(System * system, Vertex * x_i, IntegerTree * T){
  bool forcedDown = false;
  Edge * whiteEdge = x_i->first[WHITE];
  Edge * grayReverseEdge = x_i->first[GRAY_REVERSE];
  while( whiteEdge != NULL && grayReverseEdge != NULL && !forcedDown ){
    if( (grayReverseEdge->head->index == whiteEdge->head->index) 
        && (intToHalfInt( whiteEdge->weight ) == x_i->a + whiteEdge->head->a) 
        && (intToHalfInt( grayReverseEdge->weight ) == x_i->a - whiteEdge->head->a) ){
      
      x_i->Z = halfIntToInt( halfIntFloor( x_i->a ) );
      
      Edge * newEdge = (Edge *) malloc( sizeof(Edge) );
      newEdge->weight = x_i->Z;
      newEdge->type = WHITE;
      newEdge->tail = x_i;
      newEdge->head = &system->graph[0];
      newEdge->reverse = NULL;
      newEdge->next = system->integerTreeAdditionsFirst;
      system->integerTreeAdditionsFirst = newEdge;
      newEdge->allNext = NULL;
      newEdge->allPrev = NULL;
      newEdge->inAllEdgeList = false;
      
      modifyIntegerTree( system, T, x_i, &system->graph[0], whiteEdge, grayReverseEdge, newEdge );
      
      forcedDown = true;
      
    }
    else {
      while( grayReverseEdge != NULL && grayReverseEdge->head->index < whiteEdge->head->index ){
        grayReverseEdge = grayReverseEdge->next;
      }
      if( grayReverseEdge != NULL ){
        while( whiteEdge != NULL && whiteEdge->head->index < grayReverseEdge->head->index ){
          whiteEdge = whiteEdge->next;
        }
      }
    }
  }
  bool forcedUp = false;
  Edge * blackEdge = x_i->first[BLACK];
  Edge * grayForwardEdge = x_i->first[GRAY_FORWARD];
  while( blackEdge != NULL && grayForwardEdge != NULL && !forcedUp ){
    if( (grayForwardEdge->head->index == blackEdge->head->index) 
        && (intToHalfInt( blackEdge->weight ) == -x_i->a - blackEdge->head->a) 
        && (intToHalfInt( grayForwardEdge->weight ) == -x_i->a + whiteEdge->head->a) ){
      
      x_i->Z = halfIntToInt( halfIntCeil( x_i->a ) );
      
      Edge * newEdge = (Edge *) malloc( sizeof(Edge) );
      newEdge->weight = -x_i->Z;
      newEdge->type = BLACK;
      newEdge->tail = x_i;
      newEdge->head = &system->graph[0];
      newEdge->reverse = NULL;
      newEdge->next = system->integerTreeAdditionsFirst;
      system->integerTreeAdditionsFirst = newEdge;
      newEdge->allNext = NULL;
      newEdge->allPrev = NULL;
      newEdge->inAllEdgeList = false;
      
      modifyIntegerTree( system, T, x_i, &system->graph[0], blackEdge, grayForwardEdge, newEdge );
      
      forcedUp = true;
      
    }
    else {
      while( grayForwardEdge != NULL && grayForwardEdge->head->index < blackEdge->head->index ){
        grayForwardEdge = grayForwardEdge->next;
      }
      if( grayForwardEdge != NULL ){
        while( blackEdge != NULL && blackEdge->head->index < grayForwardEdge->head->index ){
          blackEdge = blackEdge->next;
        }
      }
    }
  }
  if( forcedDown && forcedUp ){
    return copyTreeEdgesToList(NULL, x_i->integerTreeVertex);
  }
  return NULL;
}

void checkDependencies(System * system, Vertex * x_i, IntegerTree * T){
  
}

void modifyIntegerTree(System * system, IntegerTree * T, Vertex * active, Vertex * parent, Edge * edge0, Edge * edge1, Edge * edge2){
  
  IntegerTreeVertex * itv = x_i->integerTreeVertex[T->type];
  if( itv == NULL ){
    itv = (IntegerTreeVertex *) malloc( sizeof(IntegerTreeVertex) );
    itv->parent = parent->integerTreeVertex;
    itv->nextSibling = itv->parent->firstChild;
    itv->parent->firstChild = itv;
    itv->firstChild = NULL;
    itv->queueNext = NULL; //queueNext points at newer queue entries, as opposed to every other linked list used here.
    if( T->queueFirst != NULL ){
      T->queueFirst->queueNext = itv;
    }
    T->queueFirst = itv;
    x_i->integerTreeVertex[T->type] = itv;
  }
  
  if( edge0 != NULL ){
    EdgeRefList * newEdgeRefList = (EdgeRefList *) malloc( sizeof(EdgeRefList) );
    newEdgeRefList->edge = edge0;
    newEdgeRefList->next = itv->firstEdge;
    itv->firstEdge = newEdgeRefList;
    if( edge1 != NULL ){
      newEdgeRefList = (EdgeRefList *) malloc( sizeof(EdgeRefList) );
      newEdgeRefList->edge = edge1;
      newEdgeRefList->next = itv->firstEdge;
      itv->firstEdge = newEdgeRefList;
      if( edge2 != NULL ){
        newEdgeReflist = (EdgeRefList *) malloc( sizeof(EdgeReflist) );
        newEdgeRefList->edge = edge2;
        newEdgeRefList->next = itv->firstEdge;
        itv->firstEdge = newEdgeRefList;
      }
    }
  }

}

EdgeRefList * copyTreeEdgesToList(EdgeRefList * list, IntegerTreeVertex * itv){
  EdgeRefList * treeEdge = itv->firstEdge;
  while( treeEdge != NULL ){
    EdgeRefList * newEdgeRefList = (EdgeRefList *) malloc( sizeof(EdgeRefList) );
    newEdgeRefList->edge = treeEdge->edge;
    newEdgeRefList->next = list;
    list = newEdgeRefList;
    treeEdge = treeEdge->next;
  }
  return list;
}

void cleanupIntegerTree(IntegerTree * tree){
  //Fill in later.
}

void cleanupSystem(System * system){
  for(EdgeType i = WHITE; i <= GRAY_REVERSE; i++){
    for(int j = 0; j < system->vertexCount; j++){
      Edge * e = system->graph[j].first[i];
      while(e != NULL){
        Edge * oldE = e;
        e = e->next;
        free( oldE );
      }
    }
  }
  free( system->graph );
}
