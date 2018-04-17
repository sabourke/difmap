#include <stdlib.h>
#include <stdio.h>
#include "utils.h"

static void insert_in_heap(float arrin[], int indx[], int node, int num_node,
			   int new_el);

/*.......................................................................
  This function takes a data array, arrin[] with npts elements and
  returns, via the argument list, an index array *indx. The sort method
  is the heap-sort, an NlogN process.  IMPORTANT: The indx arrays is
  allocated in this function, so the calling routine MUST deallocate it
  after use. On error in this routine, -1 is returned and the indx
  array is automatically zapped.
*/
int indexx(int npts, float arrin[], int **indx)
{
        void insert_in_heap(float arrin[], int indx[], int start_node, int num_node, int new_el);
	int indxt,i;
/*
  Allocate memory for the index array.
*/
	if( (*indx = (int *) calloc(npts+1, sizeof(int))) == NULL) {
	  lprintf(stderr, "sort: Memory allocation of index array failed.\n");
	  return -1;
	};
/*
  Initialize the index array elements with their element numbers.
  The data array will be indexed through this array for the
  comparisons between data array elements.
*/
	for (i=0; i<npts; i++) (*indx)[i]=i;
/*
  The algorithm fails for npts=1. However no sorting is required in this case
  anyway so return with no error.
*/
	if(npts==1)
	  return 0;
/*
  Re-arrange the index array to index the data array as a heap.
  Start at the bottom node of the heap tree, i=npts/2-1 and
  work upwards, rearanging the values in the tree until all
  branches have lower values than there nodes.
*/
	for (i = npts/2 - 1 ; i >= 0; i--)
	  insert_in_heap(arrin, *indx, i, npts, (*indx)[i]);
/*
  Now parse the tree, from the top node downwards, removing the
  value of the top node, and recursively promoting values in sub-nodes to
  fill the gap. The removed value is placed at the top of the index array,
  in the element vacated by the last promotion. 'i' keeps a record
  of the number of nodes still to be removed.
*/
	for(i = npts-1; i > 0; i--) {
/*
  Remove the root value and copy to its resting place at the end
  of the un-treated portion of the array. This overwrites the last
  element, so keep a temporary record of it in indxt.
*/
	  indxt = (*indx)[i];
	  (*indx)[i] = (*indx)[0];
/*
  Now follow down the branches, promoting the highest valued branches.
*/
	  insert_in_heap(arrin, *indx, 0, i, indxt);
	};
	return 0;
}

/*.......................................................................
  This function is called by indexx. Given the element number, new_el of
  the data array, arrin[], search arrin[] as a heap, from the node at
  element number, node, to find the correct position for the new
*/
static void insert_in_heap(float arrin[], int indx[], int node, int num_node, int new_el)
{
        int branch, right;
	float temp_value;
/*
  node holds a record of the current node being treated. At the start
  of each iteration, branch points to the leftmost branch of that node.
*/
	branch = node + node + 1;
/*
  Keep a record of the value pointed to by new_el in the data array.
*/
	temp_value = arrin[new_el];
/*
  Follow the tree branches, promoting branch values where necessary,
  until the appropriate node is found for the index to the new value
  recorded in new_value.
*/
	while (branch < num_node) {
/*
  Make 'branch' point at the branch with the highest value in it.
  Initially it points at the left branch.
*/
	  right = branch+1;
	  if (right < num_node && arrin[indx[branch]] < arrin[indx[right]]) 
	    branch = right;
/*
  Stop looking when the root value exceeds the values of both branches.
*/
	  if(temp_value >= arrin[indx[branch]])
	    break;
/*
  Otherwise promote the highest branch value and continue the search
  down that branch.
*/
	  indx[node] = indx[branch];
	  node = branch;
	  branch += branch+1;
	};
/*
  Install the index of the temp_value at the newly found location in
  the index array.
*/
	indx[node]=new_el;
	return;
}

/*.......................................................................
  Given the axis numbers specified in axis[3] and the number of elements
  per axis contained in ndim[3], calculate the element
  increments required to step through an array of the size indicated by
  ndim[], in the order specified in axis[3]. To use the resulting values,
  returned in add[3], three nested for() loops must be used with the
  first add[] nested most deeply. ndim[3] should be in the order of the
  original array.
*/
void get_increments(int axis[3], int ndim[3], int add[3])
{
        int fly[3],i,j;
	for(j=0; j<3; j++) {
	  add[j] = 1;
/*
  Increment along axis j.
*/
	  for(i=0; i<axis[j]; i++)
	    add[j] *= ndim[i];
/*
  Determine the required flyback from the end of the last axis.
*/
	  if(j>0)
	    fly[j] = add[j-1] * ndim[axis[j-1]];
	};
/*
  Include the flyback contribution.
*/
	for(j=1; j<3; j++)
	  add[j] -= fly[j];
	return;
}
