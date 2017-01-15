/*
 * voidRefList.h
 *
 * A basic list structure, storing void pointers, such that particular use cases
 * can call these functions and cast to a particular kind of pointer.
 */

#ifndef _VOID_REF_LIST_H
#define _VOID_REF_LIST_H

#include <stdbool.h>

typedef struct VoidRefList VoidRefList;

/*
 * generateVoidRefList() allocates and initializes a VoidRefList, returning a pointer to it.
 */
VoidRefList * generateVoidRefList();

/*
 * voidRefListAppend() adds a void pointer to the end of a ConstraintRefList.
 *
 * vrl - pointer to a VoidRefList to append a void pointer to
 * pointer - void pointer, which will be added to the end of vrl
 */
void voidRefListAppend(VoidRefList * vrl, void * pointer);

/*
 * voidRefListPrepend() adds a void pointer to the beginning of a VoidRefList.
 *
 * vrl - pointer to a VoidRefList to prepend with a void pointer
 * pointer - void pointer, which will be added to the beginning of vrl
 */
void voidRefListPrepend(VoidRefList * vrl, void * pointer);

/*
 * voidRefListNext() allows for the VoidRefList to be iterated through. Each call to this function returns another element within
 * the VoidRefList, until the end of the list is reached, when NULL is returned. New elements should not be added to the
 * VoidRefList while this process is ongoing. Reaching the end of the list resets the internal iterator which enables this
 * process. If iteration through the list ends before the list has been fully traversed, voidRefListIteratorReset() should be
 * called.
 *
 * vrl - the VoidRefLIst to be iterated through
 */
void * voidRefListNext(VoidRefList * vrl);

/*
 * voidRefListIteratorReset() resets the VoidRefList's internal iterator after a partial traverse of the VoidRefLIst has been
 * conducted using voidRefListNext().
 *
 * vrl - the VoidRefList
 */
void voidRefListIteratorReset(VoidRefList * vrl);

/*
 * freeVoidRefList() frees the input VoidRefList, as well as referenced pointers, if so desired.
 *
 * vrl - pointer to the VoidRefList to be freed
 * freeMemberPointers - true if data stored in the list should be f reed, along with the list itself. False, otherwise.
 */
void freevoidRefList(VoidRefList * vrl, bool freeMemberPointers);
