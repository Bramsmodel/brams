#ifndef __spAllocate__
#define __spAllocate__
/*
 * PROTOTIPES
 */
int RecordAllocation(MatrixPtr, char *);
int InitializeElementBlocks(MatrixPtr, int, int);
int AllocateBlockOfAllocationList(MatrixPtr);
#endif
