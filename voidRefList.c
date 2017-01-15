#include "voidRefList.h"

typedef struct VoidRefListNode VoidRefListNode;

/*
 * The VoidRefList is a way to reference a list of void pointers
 *
 * first - pointer to the first VoidRefListNode within the VoidRefList
 * last - pointer to the last VoidRefListNode within the VoidRefList
 * iterator - pointer to the VoidReflistNode to be considered by voidRefListNext()
 */
struct VoidRefList {
  VoidRefListNode * first;
  VoidRefListNode * last;
  VoidRefListNode * iterator;
};

/*
 * The VoidRefListNode is one node within the VoidRefList
 *
 * pointer - void pointer
 * next - pointer to the next VoidRefListNode within the VoidRefList
 */
struct VoidRefListNode {
  void * pointer;
  VoidRefListNode * next;
};

/*
 * generateVoidRefList() allocates and initializes a VoidRefList, returning a pointer to it.
 */
VoidRefList * generateVoidRefList(){
  VoidRefList * newVRL = (VoidRefList *) malloc( sizeof(VoidRefList) );
  newVRL->first = NULL;
  newVRL->last = NULL;
  newVRL->iterator = NULL;
  return newVRL;
}

/*
 * voidRefListAppend() adds a void pointer to the end of a VoidRefList.
 *
 * vrl - pointer to a VoidRefList to append a void pointer to
 * pointer - void pointer, which will be added to the end of vrl
 */
void voidRefListAppend(VoidRefList * vrl, void * pointer){
  VoidRefListNode * newVRLN = (VoidRefListNode *) malloc( sizeof(VoidRefListNode) );
  newVRLN->pointer = pointer;
  newVRLN->next = NULL;
  if( vrl->first == NULL ){
    vrl->first = newVRLN;
    vrl->iterator = newVRLN;
  }
  else {
    vrl->last->next = newVRLN;
  }
  vrl->last = newVRLN;
}

/*
 * voidRefListPrepend() adds a void pointer to the beginning of a VoidRefList.
 *
 * vrl - pointer to a VoidRefList to prepend with a void pointer
 * pointer - void pointer, which will be added to the beginning of vrl
 */
void voidRefListPrepend(VoidRefList * vrl, void * pointer){
  VoidRefListNode * newVRLN = (VoidRefListNode *) malloc( sizeof(VoidRefListNode) );
  newVRLN->pointer = pointer;
  newVRLN->next = vrl->first;
  vrl->first = newVRLN;
  vrl->iterator = newVRLN;
  if( vrl->last == NULL ){
    vrl->last = newVRLN;
  }
}

/*
 * voidRefListNext() allows for the VoidRefList to be iterated through. Each call to this function returns another element within
 * the VoidRefList, until the end of the list is reached, when NULL is returned. New elements should not be added to the
 * VoidRefList while this process is ongoing. Reaching the end of the list resets the internal iterator which enables this
 * process. If iteration through the list ends before the list has been fully traversed, voidRefListIteratorReset() should be
 * called.
 *
 * vrl - the VoidRefLIst to be iterated through
 */
void * voidRefListNext(VoidRefList * vrl){
  if( vrl->iterator == NULL ){
    vrl->iterator = vrl->first;
    return NULL;
  }
  else{
    void * output = vrl->iterator->pointer;
    vrl->iterator = vrl->iterator->next;
    return output;
  }
}

/*
 * voidRefListIteratorReset() resets the VoidRefList's internal iterator after a partial traverse of the VoidRefLIst has been
 * conducted using voidRefListNext().
 *
 * vrl - the VoidRefList
 */
void voidRefListIteratorReset(VoidRefList * vrl){
  vrl->iterator = vrl->first;
}

/*
 * freeVoidRefList() frees the input VoidRefList, as well as referenced pointers, if so desired.
 *
 * vrl - pointer to the VoidRefList to be freed
 * freeMemberPointers - true if data stored in the list should be f reed, along with the list itself. False, otherwise.
 */
void freeConstraintRefList(ConstraintRefList * crl, bool freeMemberPointers){
  if( crl != NULL ){
    VoidRefListNode * vrln = vrl->first;
    while( vrln != NULL ){
      VoidRefListNode * oldVRLN = vrln;
      vrln = vrln->next;
      if( freeMemberPointers ){
        free( oldVRLN->pointer );
      }
      free( oldVRLN );
    }
    free( vrl );
  }
}
