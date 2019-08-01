// Code to produce a tree from an index using prufer codes
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include "timerc.h"
#include <omp.h>
#include <sys/resource.h>
#include <sys/time.h>

// we probably will never consider trees with more than 13 nodes so we can use this number as a fixed bound which helps optimize the code
#define MAX_N_NODES 13
#define delta (1e6)
#define topk 100

// this allows us to change the precision with which we compute our results
typedef float realnumber;
typedef unsigned long int longint;
typedef short int shortint;



// this allows us to catch errors produced by the GPU
#define gerror(ans) { gpuAssert((ans), __FILE__, __LINE__); }
inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort=true)
{
    if (code != cudaSuccess)
    {
        fprintf(stderr,"GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}



// maybe we should change this to a simple array later.
typedef struct {
    shortint first;
    shortint second;
} edge;


//
//
//
//
// ------------------  PRIOR ON TREE COST
//
//      ---------------------- THIS MIGHT BE USED TO EXCLUDE SOME TOPOLOGIES THAT WE KNOW WILL NOT HAPPEN
//				----- IF WE WANT TO EXCLUDE SOME TOPOLOGY, WE CAN SIMPLY RETURN A VALUE OF 10000 OR ABOVE
//

__device__ __host__ realnumber tree_prior(shortint num_tree_vertices, shortint *adjacency_mat, shortint * adjacency_list, shortint * degrees){

	return 0;

}




//
//
//
// ***************************** START OF SORTING CODE ***********************************************************************************************************
//
//
// we need to do several sorting operations so here we have efficient code to do sorting

//This sorts smaller to larger
__device__ __host__ void insertion_sort(realnumber *a, shortint n)
{
    for (shortint i = 1; i < n; i++)
    {
        realnumber tmp = a[i];
        shortint j = i;
        for (; j && tmp < a[j - 1]; --j)
            a[j] = a[j - 1];
        a[j] = tmp;
    }
}

#define CMP_SWAP(i, j) if (a[i] > a[j])  { realnumber tmp = a[i]; a[i] = a[j]; a[j] = tmp; }

// for the sorting networks the following websites are useful    http://jgamble.ripco.net/cgi-bin/nw.cgi?inputs=7&algorithm=best&output=svg       and      http://stackoverflow.com/questions/4770651/what-is-the-fastest-sorting-algorithm-for-a-small-number-of-integers

__device__ __host__ void sort14_network(realnumber *a){
    CMP_SWAP(0,1);CMP_SWAP(2,3);CMP_SWAP(4,5);CMP_SWAP(6,7);CMP_SWAP(8,9);CMP_SWAP(10,11);CMP_SWAP(12,13);
    CMP_SWAP(0,2);CMP_SWAP(4,6);CMP_SWAP(8,10);CMP_SWAP(1,3);CMP_SWAP(5,7);CMP_SWAP(9,11);
    CMP_SWAP(0,4);CMP_SWAP(8,12);CMP_SWAP(1,5);CMP_SWAP(9,13);CMP_SWAP(2,6);CMP_SWAP(3,7);
    CMP_SWAP(0,8);CMP_SWAP(1,9);CMP_SWAP(2,10);CMP_SWAP(3,11);CMP_SWAP(4,12);CMP_SWAP(5,13);
    CMP_SWAP(5,10);CMP_SWAP(6,9);CMP_SWAP(3,12);CMP_SWAP(7,11);CMP_SWAP(1,2);CMP_SWAP(4,8);
    CMP_SWAP(1,4);CMP_SWAP(7,13);CMP_SWAP(2,8);CMP_SWAP(5,6);CMP_SWAP(9,10);
    CMP_SWAP(2,4);CMP_SWAP(11,13);CMP_SWAP(3,8);CMP_SWAP(7,12);
    CMP_SWAP(6,8);CMP_SWAP(10,12);CMP_SWAP(3,5);CMP_SWAP(7,9);
    CMP_SWAP(3,4);CMP_SWAP(5,6);CMP_SWAP(7,8);CMP_SWAP(9,10);CMP_SWAP(11,12);
    CMP_SWAP(6,7);CMP_SWAP(8,9);
}


__device__ __host__ void sort13_network(realnumber *a){
    CMP_SWAP(1,7);CMP_SWAP(9,11);CMP_SWAP(3,4);CMP_SWAP(5,8);CMP_SWAP(0,12);CMP_SWAP(2,6);
    CMP_SWAP(0,1);CMP_SWAP(2,3);CMP_SWAP(4,6);CMP_SWAP(8,11);CMP_SWAP(7,12);CMP_SWAP(5,9);
    CMP_SWAP(0,2);CMP_SWAP(3,7);CMP_SWAP(10,11);CMP_SWAP(1,4);CMP_SWAP(6,12);
    CMP_SWAP(7,8);CMP_SWAP(11,12);CMP_SWAP(4,9);CMP_SWAP(6,10);
    CMP_SWAP(3,4);CMP_SWAP(5,6);CMP_SWAP(8,9);CMP_SWAP(10,11);CMP_SWAP(1,7);
    CMP_SWAP(2,6);CMP_SWAP(9,11);CMP_SWAP(1,3);CMP_SWAP(4,7);CMP_SWAP(8,10);CMP_SWAP(0,5);
    CMP_SWAP(2,5);CMP_SWAP(6,8);CMP_SWAP(9,10);
    CMP_SWAP(1,2);CMP_SWAP(3,5);CMP_SWAP(7,8);CMP_SWAP(4,6);
    CMP_SWAP(2,3);CMP_SWAP(4,5);CMP_SWAP(6,7);CMP_SWAP(8,9);
    CMP_SWAP(3,4);CMP_SWAP(5,6);
}

__device__ __host__ void sort12_network(realnumber *a){
    CMP_SWAP(0,1);CMP_SWAP(2,3);CMP_SWAP(4,5);CMP_SWAP(6,7);CMP_SWAP(8,9);CMP_SWAP(10,11);
    CMP_SWAP(1,3);CMP_SWAP(5,7);CMP_SWAP(9,11);CMP_SWAP(0,2);CMP_SWAP(4,6);CMP_SWAP(8,10);
    CMP_SWAP(1,2);CMP_SWAP(5,6);CMP_SWAP(9,10);CMP_SWAP(0,4);CMP_SWAP(7,11);
    CMP_SWAP(1,5);CMP_SWAP(6,10);CMP_SWAP(3,7);CMP_SWAP(4,8);
    CMP_SWAP(5,9);CMP_SWAP(2,6);CMP_SWAP(0,4);CMP_SWAP(7,11);CMP_SWAP(3,8);
    CMP_SWAP(1,5);CMP_SWAP(6,10);CMP_SWAP(2,3);CMP_SWAP(8,9);
    CMP_SWAP(1,4);CMP_SWAP(7,10);CMP_SWAP(3,5);CMP_SWAP(6,8);
    CMP_SWAP(2,4);CMP_SWAP(7,9);CMP_SWAP(5,6);
    CMP_SWAP(3,4);CMP_SWAP(7,8);
}


__device__ __host__ void sort11_network(realnumber *a){
    CMP_SWAP(0,1);CMP_SWAP(2,3);CMP_SWAP(4,5);CMP_SWAP(6,7);CMP_SWAP(8,9);
    CMP_SWAP(1,3);CMP_SWAP(5,7);CMP_SWAP(0,2);CMP_SWAP(4,6);CMP_SWAP(8,10);
    CMP_SWAP(1,2);CMP_SWAP(5,6);CMP_SWAP(9,10);CMP_SWAP(0,4);CMP_SWAP(3,7);
    CMP_SWAP(1,5);CMP_SWAP(6,10);CMP_SWAP(4,8);
    CMP_SWAP(5,9);CMP_SWAP(2,6);CMP_SWAP(0,4);CMP_SWAP(3,8);
    CMP_SWAP(1,5);CMP_SWAP(6,10);CMP_SWAP(2,3);CMP_SWAP(8,9);
    CMP_SWAP(1,4);CMP_SWAP(7,10);CMP_SWAP(3,5);CMP_SWAP(6,8);
    CMP_SWAP(2,4);CMP_SWAP(7,9);CMP_SWAP(5,6);
    CMP_SWAP(3,4);CMP_SWAP(7,8);
}


__device__ __host__ void sort10_network(realnumber *a){
    CMP_SWAP(4,9);CMP_SWAP(3,8);CMP_SWAP(2,7);CMP_SWAP(1,6);CMP_SWAP(0,5);
    CMP_SWAP(1,4);CMP_SWAP(6,9);CMP_SWAP(0,3);CMP_SWAP(5,8);
    CMP_SWAP(0,2);CMP_SWAP(3,6);CMP_SWAP(7,9);
    CMP_SWAP(0,1);CMP_SWAP(2,4);CMP_SWAP(5,7);CMP_SWAP(8,9);
    CMP_SWAP(1,2);CMP_SWAP(4,6);CMP_SWAP(7,8);CMP_SWAP(3,5);
    CMP_SWAP(2,5);CMP_SWAP(6,8);CMP_SWAP(1,3);CMP_SWAP(4,7);
    CMP_SWAP(2,3);CMP_SWAP(6,7);
    CMP_SWAP(3,4);CMP_SWAP(5,6);
    CMP_SWAP(4,5);
}

__device__ __host__ void sort9_network(realnumber *a){
    CMP_SWAP(0,1);CMP_SWAP(3,4);CMP_SWAP(6,7);
    CMP_SWAP(1,2);CMP_SWAP(4,5);CMP_SWAP(7,8);
    CMP_SWAP(0,1);CMP_SWAP(3,4);CMP_SWAP(6,7);CMP_SWAP(2,5);
    CMP_SWAP(0,3);CMP_SWAP(1,4);CMP_SWAP(5,8);
    CMP_SWAP(3,6);CMP_SWAP(4,7);CMP_SWAP(2,5);
    CMP_SWAP(0,3);CMP_SWAP(1,4);CMP_SWAP(5,7);CMP_SWAP(2,6);
    CMP_SWAP(1,3);CMP_SWAP(4,6);
    CMP_SWAP(2,4);CMP_SWAP(5,6);
    CMP_SWAP(2,3);
}

__device__ __host__ void sort8_network(realnumber *a){
    CMP_SWAP(0, 1); CMP_SWAP(2, 3); CMP_SWAP(4, 5); CMP_SWAP(6, 7);
    CMP_SWAP(0, 2); CMP_SWAP(1, 3); CMP_SWAP(4, 6); CMP_SWAP(5, 7);
    CMP_SWAP(1, 2); CMP_SWAP(5, 6); CMP_SWAP(0, 4); CMP_SWAP(3, 7);
    CMP_SWAP(1, 5); CMP_SWAP(2, 6);
    CMP_SWAP(1, 4); CMP_SWAP(3, 6);
    CMP_SWAP(2, 4); CMP_SWAP(3, 5);
    CMP_SWAP(3, 4);
}

__device__ __host__ void sort7_network(realnumber *a){
    CMP_SWAP(1,2);CMP_SWAP(3,4);CMP_SWAP(5,6);
    CMP_SWAP(0,2);CMP_SWAP(3,5);CMP_SWAP(4,6);
    CMP_SWAP(0,1);CMP_SWAP(4,5);CMP_SWAP(2,6);
    CMP_SWAP(0,4);CMP_SWAP(1,5);
    CMP_SWAP(0,3);CMP_SWAP(2,5);
    CMP_SWAP(1,3);CMP_SWAP(2,4);
    CMP_SWAP(2,3);
}

__device__ __host__ void sort6_network(realnumber *a){
    CMP_SWAP(1,2);CMP_SWAP(4,5);
    CMP_SWAP(0,2);CMP_SWAP(3,5);
    CMP_SWAP(0,1);CMP_SWAP(3,4);CMP_SWAP(2,5);
    CMP_SWAP(0,3);CMP_SWAP(1,4);
    CMP_SWAP(2,4);CMP_SWAP(1,3);
    CMP_SWAP(2,3);
}

__device__ __host__ void sort5_network(realnumber *a){
    CMP_SWAP(0,1);CMP_SWAP(3,4);
    CMP_SWAP(2,4);
    CMP_SWAP(2,3);CMP_SWAP(1,4);
    CMP_SWAP(0,3);
    CMP_SWAP(0,2);CMP_SWAP(1,3);
    CMP_SWAP(1,2);
}

__device__ __host__ void sort4_network(realnumber *a){
    CMP_SWAP(0,1);CMP_SWAP(2,3);
    CMP_SWAP(0,2);CMP_SWAP(1,3);
    CMP_SWAP(1,2);
}

__device__ __host__ void sort3_network(realnumber *a){
    CMP_SWAP(1,2);
    CMP_SWAP(0,2);
    CMP_SWAP(0,1);
}

__device__ __host__ void sort2_network(realnumber *a){
    CMP_SWAP(0,1);
}


//This sorts smaller to larger
__device__ __host__ void sort_using_networks(realnumber *a, shortint n){
    
    // we only process if the number of elements is  between 2 and 14
    if (n==1 || n > 14){
        return;
    }
    
    switch(n) {
            
        case 2  :
            sort2_network(a);
            break;
            
        case 3  :
            sort3_network(a);
            break;
            
        case 4  :
            sort4_network(a);
            break;
            
        case 5  :
            sort5_network(a);
            break;
            
        case 6  :
            sort6_network(a);
            break;
            
        case 7  :
            sort7_network(a);
            break;
            
        case 8  :
            sort8_network(a);
            break;
            
        case 9  :
            sort9_network(a);
            break;
            
        case 10  :
            sort10_network(a);
            break;
            
        case 11  :
            sort11_network(a);
            break;
            
        case 12  :
            sort12_network(a);
            break;
            
        case 13  :
            sort13_network(a);
            break;
            
        case 14  :
            sort14_network(a);
            break;
    }
    
}



//
//
//
// ***************************** END OF SORTING CODE ***********************************************************************************************************
//
//






//
//
//
//****************************** START OF BEI'S CODE WITH ARRAYS NO MALLOCS NO STRUCTS ******************************
//
//
//
//
//

/*
 LIST OF ARRAYS ABOUT TREE:
 realnumber value[], size = 2*MAX_N_NODES + 1
 realnumber z[], size = 2*MAX_N_NODES + 1
 realnumber status[], size = 2*MAX_N_NODES + 1
 status: (0 : free), (1 : branch), (2 : fixed), (3 : effective)
 int father[], size = 2*MAX_N_NODES + 1
 int num_kids[], size = MAX_N_NODES + 1 ---> Use final_degrees[] instead
 int kids[], size = (MAX_N_NODES + 1) * (MAX_N_NODES + 1)
 */
/*
 QUEUE OPERATIONS:
 push_back: back = (back + 1) % capacity;
 ++size;
 pop_front: front = (front + 1) % capacity;
 --size;
 */
/*
 STACK OPERATIONS:
 push_back: back = (back + 1) % capacity;
 ++size;
 pop_back: --size;
 back = (back - 1) % capacity;
 */

// Some path finding functions:

/*
 This function finds the path from the node at start_ind to the node at node_ind, and store the path in holder (a stack)
 */
__host__ __device__ void up_to_node_ind(shortint start_ind, shortint node_ind, shortint status[], shortint level[], shortint father[], shortint num_kids[], shortint kids[], shortint holder[], shortint *holder_front, shortint *holder_back, shortint *holder_size) {
    if (start_ind == node_ind) {
        return;
    }
    shortint current_ind = father[start_ind];
    if (current_ind < 0) {
        return;
    }
    *holder_front = 0;
    *holder_back = 0;
    *holder_size = 1;
    holder[*holder_back] = current_ind;
    while (current_ind != node_ind) {
        current_ind = father[current_ind];
        if (current_ind < 0) {
            return;
        }
        *holder_back = (*holder_back + 1) % (MAX_N_NODES + 1);
        holder[*holder_back] = current_ind;
        ++(*holder_size);
    }
}

/*
 This function finds the path from the node at start_ind to the nearest branch node above it, and store the path in holder (a stack)
 */
__host__ __device__ void up_to_branch_ind(shortint start_ind, shortint status[], shortint level[], shortint father[], shortint num_kids[], shortint kids[], shortint holder[], shortint *holder_front, shortint *holder_back, shortint *holder_size) {
    shortint current_ind = father[start_ind];
    if (current_ind < 0) {
        return;
    }
    *holder_front = 0;
    *holder_back = 0;
    *holder_size = 1;
    holder[*holder_back] = current_ind;
    while (status[current_ind] != 1) {
        current_ind = father[current_ind];
        if (current_ind < 0) {
            return;
        }
        *holder_back = (*holder_back + 1) % (MAX_N_NODES + 1);
        holder[*holder_back] = current_ind;
        ++(*holder_size);
    }
}

/*
 This function uses BFS on the entire tree, find all fixed nodes and store them in fixed_queue[], and levels of all nodes and stored in level[]
 */
__host__ __device__ void bfs_get_level(shortint status[], shortint level[], shortint num_kids[], shortint kids[], shortint fixed_queue[], shortint *fixed_queue_front, shortint *fixed_queue_back, shortint *fixed_queue_size) {
    shortint queue[MAX_N_NODES + 1];
    shortint queue_front = 0;
    queue[queue_front] = 0;
    shortint queue_back = 0;
    shortint queue_size = 1;
    level[0] = 0;
    while (queue_size != 0) {
        shortint current_ind = queue[queue_front];
        queue_front = (queue_front + 1) % (MAX_N_NODES + 1);
        --queue_size;
        if (status[current_ind] == 2) {
            *fixed_queue_back = (*fixed_queue_back + 1) % (MAX_N_NODES + 1);
            fixed_queue[*fixed_queue_back] = current_ind;
            ++(*fixed_queue_size);
        }
        for (shortint i = 0; i < num_kids[current_ind]; ++i) {
            // Note: this BFS will never look at effective nodes produced during compression
            queue_back = (queue_back + 1) % (MAX_N_NODES + 1);
            shortint kid_ind = kids[current_ind * (MAX_N_NODES + 1) + i];
            queue[queue_back] = kid_ind;
            ++queue_size;
            level[kid_ind] = level[current_ind] + 1;
        }
    }
}

/*
 This function uses BFS to find nodes up to fixed nodes (not inlcuding effective nodes), sore found branch nodes in branch_stack[], and store found fixed nodes in fixed_stack[]:
 */
__host__ __device__ void bfs_array_to_fixed(shortint start_ind, shortint status[], shortint num_kids[], shortint kids[], shortint branch_stack[], shortint *branch_stack_back, shortint *branch_stack_size, shortint fixed_stack[], shortint *fixed_stack_back, shortint *fixed_stack_size) {
    shortint queue[MAX_N_NODES + 1];
    shortint queue_front = 0;
    queue[queue_front] = start_ind;
    shortint queue_back = 0;
    shortint queue_size = 1;
    while (queue_size != 0) {
        shortint current_ind = queue[queue_front];
        queue_front = (queue_front + 1) % (MAX_N_NODES + 1);
        --queue_size;
        if (status[current_ind] == 1) {
            *branch_stack_back = (*branch_stack_back + 1) % (MAX_N_NODES + 1);
            branch_stack[*branch_stack_back] = current_ind;
            ++(*branch_stack_size);
        }
        for (shortint i = 0; i < num_kids[current_ind]; ++i) {
            // Note: this BFS will never look at effective nodes produced during compression
            shortint kid_ind = kids[current_ind * (MAX_N_NODES + 1) + i];
            if (status[kid_ind] != 2 && status[kid_ind] != 3) {
                queue_back = (queue_back + 1) % (MAX_N_NODES + 1);
                queue[queue_back] = kid_ind;
                ++queue_size;
            }
            else {
                *fixed_stack_back = (*fixed_stack_back + 1) % (MAX_N_NODES + 1);
                fixed_stack[*fixed_stack_back] = kid_ind;
                ++(*fixed_stack_size);
            }
        }
    }
}

/*
 This function uses BFS to find nodes up to fixed nodes (inlcuding effective nodes), sore found branch nodes in branch_stack[], and store found fixed nodes (and effective nodes) in fixed_stack[]:
 */
__host__ __device__ void bfs_array_to_fixed_effective(shortint start_ind, shortint status[], shortint num_kids[], shortint kids[], shortint branch_stack[], shortint *branch_stack_back, shortint *branch_stack_size, shortint fixed_stack[], shortint *fixed_stack_back, shortint *fixed_stack_size) {
    shortint queue[MAX_N_NODES + 1];
    shortint queue_front = 0;
    queue[queue_front] = start_ind;
    shortint queue_back = 0;
    shortint queue_size = 1;
    while (queue_size != 0) {
        shortint current_ind = queue[queue_front];
        queue_front = (queue_front + 1) % (MAX_N_NODES + 1);
        --queue_size;
        if (status[current_ind] == 1) {
            *branch_stack_back = (*branch_stack_back + 1) % (MAX_N_NODES + 1);
            branch_stack[*branch_stack_back] = current_ind;
            ++(*branch_stack_size);
        }
        if (num_kids[current_ind] > 0) {
            for (shortint i = 0; i <= num_kids[current_ind]; ++i) {
                // Note: this BFS will look at effective nodes produced during compression
                shortint kid_ind = kids[current_ind * (MAX_N_NODES + 1) + i];
                // If that last spot has no effective node in it:
                if (kid_ind == -2) {
                    break;
                }
                else if (status[kid_ind] != 2 && status[kid_ind] != 3) {
                    queue_back = (queue_back + 1) % (MAX_N_NODES + 1);
                    queue[queue_back] = kid_ind;
                    ++queue_size;
                }
                else {
                    // Keep record of ending fixed nodes (including effective nodes)
                    *fixed_stack_back = (*fixed_stack_back + 1) % (MAX_N_NODES + 1);
                    fixed_stack[*fixed_stack_back] = kid_ind;
                    ++(*fixed_stack_size);
                }
            }
        }
    }
}

/*
 This function uses BFS to find positions (stored in z[]) and speeds (stored in v[]) of unstretched nodes
 */
__host__ __device__ void bfs_array_unstretched(shortint status[], shortint level[], shortint father[], shortint num_kids[], shortint kids[], shortint seen[], realnumber z[], realnumber v[]) {
    shortint queue[MAX_N_NODES + 1];
    shortint queue_front = 0;
    queue[queue_front] = 0;
    shortint queue_back = 0;
    shortint queue_size = 1;
    while (queue_size != 0) {
        shortint current_ind = queue[queue_front];
        queue_front = (queue_front + 1) % (MAX_N_NODES + 1);
        --queue_size;
        if (seen[current_ind] == 0 && status[current_ind] != 2) {
            z[current_ind] = z[father[current_ind]];
            v[current_ind] = v[father[current_ind]];
            seen[current_ind] = 1;
        }
        for (shortint i = 0; i < num_kids[current_ind]; ++i) {
            // Note: this BFS will never look at effective nodes produced during compression
            queue_back = (queue_back + 1) % (MAX_N_NODES + 1);
            shortint kid_ind = kids[current_ind * (MAX_N_NODES + 1) + i];
            queue[queue_back] = kid_ind;
            ++queue_size;
        }
    }
}

/*
 This function check wether the node at node_ind is a branch node or not
 */
__host__ __device__ shortint is_branch_array(shortint node_ind, shortint status[], shortint father[], shortint num_kids[], shortint kids[], shortint dim) {
    if (status[node_ind] == 2) {
        return 0;
    }
    // Define a branch_stack;
    shortint branch_stack[MAX_N_NODES + 1];
    shortint branch_stack_back = MAX_N_NODES;
    shortint branch_stack_size = 0;
    // Define a fixed_stack;
    shortint fixed_stack[MAX_N_NODES + 1];
    shortint fixed_stack_back = MAX_N_NODES;
    shortint fixed_stack_size = 0;
    bfs_array_to_fixed(node_ind, status, num_kids, kids, branch_stack, &branch_stack_back, &branch_stack_size, fixed_stack, &fixed_stack_back, &fixed_stack_size);
    shortint count_fixed = 0;
    while (fixed_stack_size != 0) {
        shortint current_ind = fixed_stack[fixed_stack_back];
        fixed_stack_back = (fixed_stack_back - 1) % (MAX_N_NODES + 1);
        --(fixed_stack_size);
        if (current_ind > 0 && current_ind <= dim) {
            ++count_fixed;
        }
        if (count_fixed > 1) {
            return 1;
        }
    }
    return 0;
}

/*
 This function implement the compression procedure starting with a branching node at branch_ind, create and store the data of a new effective node (associated with this branch node) from compressing:
 */
__host__ __device__ void compress_array(shortint branch_ind, realnumber value[], realnumber z[], shortint status[], shortint level[], shortint father[], shortint num_kids[], shortint kids[], shortint seen_fixed[], shortint holder[], shortint *holder_front, shortint *holder_back, shortint *holder_size, shortint *eff_ind) {
    ++(*eff_ind);
    // Define a branch_stack;
    shortint branch_stack[MAX_N_NODES + 1];
    shortint branch_stack_back = MAX_N_NODES;
    shortint branch_stack_size = 0;
    // Define a fixed_stack;
    shortint fixed_stack[MAX_N_NODES + 1];
    shortint fixed_stack_back = MAX_N_NODES;
    shortint fixed_stack_size = 0;
    // Define a fixed_queue;
    shortint fixed_queue[MAX_N_NODES + 1];
    shortint fixed_queue_front = 0;
    shortint fixed_queue_back = MAX_N_NODES;
    shortint fixed_queue_size = 0;
    // Find other branch nodes and fixed nodes in effective tree:
    bfs_array_to_fixed_effective(branch_ind, status, num_kids, kids, branch_stack, &branch_stack_back, &branch_stack_size, fixed_stack, &fixed_stack_back, &fixed_stack_size);
    realnumber denominator = 0;
    realnumber numerator = 0;
    *holder_front = 0;
    *holder_back = MAX_N_NODES;
    *holder_size = 0;
    while (fixed_stack_size != 0) {
        shortint current_ind = fixed_stack[fixed_stack_back];
        fixed_stack_back = (fixed_stack_back - 1) % (MAX_N_NODES + 1);
        --(fixed_stack_size);
        fixed_queue_back = (fixed_queue_back + 1) % (MAX_N_NODES + 1);
        fixed_queue[fixed_queue_back] = current_ind;
        ++(fixed_queue_size);
        if (seen_fixed[current_ind] == 0 && status[current_ind] != 3) {
            denominator += 1.0 / (level[current_ind] - level[branch_ind]);
        }
        else if (seen_fixed[current_ind] == 0 && status[current_ind] == 3) {
            up_to_branch_ind(current_ind, status, level, father, num_kids, kids, holder, holder_front, holder_back, holder_size);
            realnumber diff = level[holder[*holder_back]] - level[branch_ind];
            realnumber dist = (1.0 / value[current_ind]) + diff;
            denominator += 1.0 / dist;
        }
    }
    *holder_front = 0;
    *holder_back = MAX_N_NODES;
    *holder_size = 0;
    while (fixed_queue_size != 0) {
        shortint current_ind = fixed_queue[fixed_queue_front];
        fixed_queue_front = (fixed_queue_front + 1) % (MAX_N_NODES + 1);
        --(fixed_queue_size);
        if (seen_fixed[current_ind] == 0 && status[current_ind] != 3) {
            numerator += z[current_ind] / (level[current_ind] - level[branch_ind]);
            seen_fixed[current_ind] = 1;
        }
        else if (seen_fixed[current_ind] == 0 && status[current_ind] == 3) {
            up_to_branch_ind(current_ind, status, level, father, num_kids, kids, holder, holder_front, holder_back, holder_size);
            realnumber diff = level[holder[*holder_back]] - level[branch_ind];
            realnumber dist = (1.0 / value[current_ind]) + diff;
            numerator += z[current_ind] / dist;
            seen_fixed[current_ind] = 1;
        }
    }
    // Create a new effective node from compressing:
    status[*eff_ind] = 3;
    value[*eff_ind] = denominator;
    z[*eff_ind] = numerator / denominator;
    father[*eff_ind] = branch_ind;
}

/*
 This function update positions z[] AND speeds v[] at a given t:
 */
__host__ __device__ void update_z_array(realnumber t, realnumber v[], shortint fixed_queue[], shortint *fixed_queue_front, shortint *fixed_queue_back, shortint *fixed_queue_size, realnumber value[], realnumber z[], shortint status[], shortint level[], shortint father[], shortint num_kids[], shortint kids[], shortint dim) {
    shortint seen[MAX_N_NODES + 1] = {0};
    // Find z for all fixed nodes:
    for (shortint i = 0; i < *fixed_queue_size; ++i) {
        shortint current_ind = fixed_queue[*fixed_queue_front];
        *fixed_queue_front = (*fixed_queue_front + 1) % (MAX_N_NODES + 1);
        --(*fixed_queue_size);
        if (current_ind == 0) {
            z[current_ind] = 0;
            v[current_ind] = 0;
        }
        else {
            z[current_ind] = t - value[current_ind];
            v[current_ind] = 1;
        }
        seen[current_ind] = 1;
        *fixed_queue_back = (*fixed_queue_back + 1) % (MAX_N_NODES + 1);
        fixed_queue[*fixed_queue_back] = current_ind;
        ++(*fixed_queue_size);
    }
    
    // Define a branch_stack;
    shortint branch_stack[MAX_N_NODES + 1];
    shortint branch_stack_back = MAX_N_NODES;
    shortint branch_stack_size = 0;
    // Define a fixed_stack;
    shortint fixed_stack[MAX_N_NODES + 1];
    shortint fixed_stack_back = MAX_N_NODES;
    shortint fixed_stack_size = 0;
    // Define a effective_stack;
    shortint effective_stack[MAX_N_NODES + 1];
    shortint effective_stack_back = MAX_N_NODES;
    shortint effective_stack_size = 0;
    // Define a branch_stack_reverse;
    shortint branch_stack_reverse[MAX_N_NODES + 1];
    shortint branch_stack_reverse_back = MAX_N_NODES;
    shortint branch_stack_reverse_size = 0;
    // Define a universal holder:
    shortint holder[MAX_N_NODES + 1];
    shortint holder_front = 0;
    shortint holder_back = MAX_N_NODES;
    shortint holder_size = 0;
    
    // Find z of all branch nodes and stretched free nodes:
    for (shortint i = 0; i < *fixed_queue_size; ++i) {
        shortint current_fixed_ind = fixed_queue[*fixed_queue_front];
        *fixed_queue_front = (*fixed_queue_front + 1) % (MAX_N_NODES + 1);
        --(*fixed_queue_size);
        shortint current_num_kids = num_kids[current_fixed_ind];
        for (shortint j = 0; j < current_num_kids; ++j) {
            shortint current_kid_ind = kids[current_fixed_ind * (MAX_N_NODES + 1) + j];
            if (status[current_kid_ind] != 2) {
                // Reset containers:
                branch_stack_back = MAX_N_NODES;
                branch_stack_size = 0;
                fixed_stack_back = MAX_N_NODES;
                fixed_stack_size = 0;
                effective_stack_back = MAX_N_NODES;
                effective_stack_size = 0;
                branch_stack_reverse_back = MAX_N_NODES;
                branch_stack_reverse_size = 0;
                holder_front = 0;
                holder_back = MAX_N_NODES;
                holder_size = 0;
                // An effective tree for each kids of this fixed node:
                bfs_array_to_fixed(current_kid_ind, status, num_kids, kids, branch_stack, &branch_stack_back, &branch_stack_size, fixed_stack, &fixed_stack_back, &fixed_stack_size);
                // If there is no branch node but there are fixed nodes in this effective tree (which will be a chain):
                if (branch_stack_size == 0 && fixed_stack_size != 0) {
                    // There must be only 1 ending fixed node
                    up_to_node_ind(fixed_stack[fixed_stack_back], current_fixed_ind, status, level, father, num_kids, kids, holder, &holder_front, &holder_back, &holder_size);
                    while (holder_size != 0) {
                        shortint current_node_ind = holder[holder_front];
                        holder_front = (holder_front + 1) % (MAX_N_NODES + 1);
                        --(holder_size);
                        z[current_node_ind] = z[fixed_stack[fixed_stack_back]] + (level[fixed_stack[fixed_stack_back]] - level[current_node_ind]) * (z[current_fixed_ind] - z[fixed_stack[fixed_stack_back]]) / (level[fixed_stack[fixed_stack_back]] - level[current_fixed_ind]);
                        if (current_fixed_ind == 0) {
                            v[current_node_ind] = z[current_node_ind] / z[fixed_stack[fixed_stack_back]];
                        }
                        else {
                            v[current_node_ind] = 1;
                        }
                        seen[current_node_ind] = 1;
                    }
                }
                // If there are branch nodes, do compression:
                else if (branch_stack_size != 0) {
                    // Compressing bottom up in this effective tree:
                    shortint seen_fixed[2*MAX_N_NODES + 1] = {0};
                    shortint eff_ind = dim;
                    while (branch_stack_size != 0) {
                        shortint branch_ind = branch_stack[branch_stack_back];
                        branch_stack_back = (branch_stack_back - 1) % (MAX_N_NODES + 1);
                        --(branch_stack_size);
                        compress_array(branch_ind, value, z, status, level, father, num_kids, kids, seen_fixed, holder, &holder_front, &holder_back, &holder_size, &eff_ind);
                        // Add effective node to kids of its branch node:
                        kids[branch_ind * (MAX_N_NODES + 1) + num_kids[branch_ind]] = eff_ind;
                        effective_stack_back = (effective_stack_back + 1) % (MAX_N_NODES + 1);
                        effective_stack[effective_stack_back] = eff_ind;
                        ++(effective_stack_size);
                        branch_stack_reverse_back = (branch_stack_reverse_back + 1) % (MAX_N_NODES + 1);
                        branch_stack_reverse[branch_stack_reverse_back] = branch_ind;
                        ++(branch_stack_reverse_size);
                    }
                    // Note: bijection between branch_stack and effective_stack
                    // Find z of the top most branch node:
                    shortint effective_ind = effective_stack[effective_stack_back];
                    effective_stack_back = (effective_stack_back - 1) % (MAX_N_NODES + 1);
                    --(effective_stack_size);
                    shortint branch_ind = branch_stack_reverse[branch_stack_reverse_back];
                    branch_stack_reverse_back = (branch_stack_reverse_back - 1) % (MAX_N_NODES + 1);
                    --(branch_stack_reverse_size);
                    realnumber dist = (level[branch_ind] - level[current_fixed_ind]) + 1.0 / value[effective_ind];
                    realnumber unit_len = (z[current_fixed_ind] - z[effective_ind]) / dist;
                    z[branch_ind] = z[effective_ind] + (1.0 / value[effective_ind]) * unit_len;
                    if (current_fixed_ind == 0) {
                        v[branch_ind] = z[branch_ind] / z[effective_ind];
                    }
                    else {
                        v[branch_ind] = 1;
                    }
                    // ! Remove effective node from tree !:
                    kids[branch_ind * (MAX_N_NODES + 1) + num_kids[branch_ind]] = -2;
                    seen[branch_ind] = 1;
                    // Find z of free nodes between top most branch node and root of current effective tree:
                    up_to_node_ind(branch_ind, current_fixed_ind, status, level, father, num_kids, kids, holder, &holder_front, &holder_back, &holder_size);
                    while (holder_size != 0) {
                        shortint current_node_ind = holder[holder_front];
                        holder_front = (holder_front + 1) % (MAX_N_NODES + 1);
                        --(holder_size);
                        z[current_node_ind] = z[branch_ind] + (level[branch_ind] - level[current_node_ind]) * unit_len;
                        if (current_fixed_ind == 0) {
                            v[current_node_ind] = z[current_node_ind] / z[branch_ind];
                        }
                        else {
                            v[current_node_ind] = 1;
                        }
                        seen[current_node_ind] = 1;
                    }
                    // Decompressing:
                    while (effective_stack_size != 0) {
                        shortint current_eff_ind = effective_stack[effective_stack_back];
                        effective_stack_back = (effective_stack_back - 1) % (MAX_N_NODES + 1);
                        --(effective_stack_size);
                        shortint current_branch_ind = branch_stack_reverse[branch_stack_reverse_back];
                        branch_stack_reverse_back = (branch_stack_reverse_back - 1) % (MAX_N_NODES + 1);
                        --(branch_stack_reverse_size);
                        up_to_branch_ind(current_branch_ind, status, level, father, num_kids, kids, holder, &holder_front, &holder_back, &holder_size);
                        dist = (level[current_branch_ind] - level[holder[holder_back]]) + 1.0 / value[current_eff_ind];
                        unit_len = (z[holder[holder_back]] - z[current_eff_ind]) / dist;
                        z[current_branch_ind] = z[current_eff_ind] + (1.0 / value[current_eff_ind]) * unit_len;
                        if (current_fixed_ind == 0) {
                            v[current_branch_ind] = z[current_branch_ind] / z[current_eff_ind];
                        }
                        else {
                            v[current_branch_ind] = 1;
                        }
                        seen[current_branch_ind] = 1;
                        // Find z of all the free nodes between current branch node and its closest up branch node:
                        holder_front = (holder_front + 1) % (MAX_N_NODES + 1);
                        --(holder_size);
                        while (holder_size != 0) {
                            shortint free_ind = holder[holder_front];
                            holder_front = (holder_front + 1) % (MAX_N_NODES + 1);
                            --(holder_size);
                            z[free_ind] = z[current_branch_ind] + (level[current_branch_ind] - level[free_ind]) * unit_len;
                            if (current_fixed_ind == 0) {
                                v[free_ind] = z[free_ind] / z[current_branch_ind];
                            }
                            else {
                                v[free_ind] = 1;
                            }
                            seen[free_ind] = 1;
                        }
                        // ! Remove effective node from tree !:
                        kids[current_branch_ind * (MAX_N_NODES + 1) + num_kids[current_branch_ind]] = -2;
                    }
                    // Find z of the rest stretched free nodes in current effective tree:
                    while (fixed_stack_size != 0) {
                        shortint a_fixed_ind = fixed_stack[fixed_stack_back];
                        fixed_stack_back = (fixed_stack_back - 1) % (MAX_N_NODES + 1);
                        --(fixed_stack_size);
                        if (a_fixed_ind != 0) {
                            up_to_branch_ind(a_fixed_ind, status, level, father, num_kids, kids, holder, &holder_front, &holder_back, &holder_size);
                            shortint holder_back_ind = holder[holder_back];
                            while (holder_size != 0) {
                                shortint a_free_ind = holder[holder_front];
                                holder_front = (holder_front + 1) % (MAX_N_NODES + 1);
                                --(holder_size);
                                z[a_free_ind] = z[a_fixed_ind] + (level[a_fixed_ind] - level[a_free_ind]) * (z[holder_back_ind] - z[a_fixed_ind]) / (level[a_fixed_ind] - level[holder_back_ind]);
                                if (current_fixed_ind == 0) {
                                    v[a_free_ind] = z[a_free_ind] / z[a_fixed_ind];
                                }
                                else {
                                    v[a_free_ind] = 1;
                                }
                                seen[a_free_ind] = 1;
                            }
                        }
                    }
                }
            }
        }
        // Add the current fixed node back to keep fixed_queue intact:
        *fixed_queue_back = (*fixed_queue_back + 1) % (MAX_N_NODES + 1);
        fixed_queue[*fixed_queue_back] = current_fixed_ind;
        ++(*fixed_queue_size);
    }
    // Now we have seen all those stretched nodes
    // Find z of all other (free) nodes (unstretched):
    bfs_array_unstretched(status, level, father, num_kids, kids, seen, z, v);
    v[0] = 0;
}

/*
 This function compute the next turing point, then the next node to be fixed and store its index in t_data[1]:
 */
__host__ __device__ void next_turn_array(realnumber t_data[], realnumber t_pre, realnumber z_pre[], realnumber v[], realnumber value[], shortint status[], shortint level[], shortint dim) {
    v[0] = 0;
    realnumber t_possilbe[MAX_N_NODES + 1];
    for (shortint i = 0; i <= dim; ++i) {
        t_possilbe[i] = -INFINITY;
    }
    // Find intersections:
    for (shortint i = 1; i <= dim; ++i) {
        if (status[i] != 2) {
            t_possilbe[i] = (z_pre[i] - v[i]*t_pre + value[i]) / (1 - v[i]);
        }
    }
    // Find the max t among all intersections:
    realnumber max = t_possilbe[0];
    shortint max_ind = 0;
    for (shortint i = 1; i <= dim; ++i) {
        if (t_possilbe[i] >= max) {
            max = t_possilbe[i];
        }
    }
    shortint top_level = dim + 1;
    for (shortint i = 1; i <= dim; ++i) {
        if (t_possilbe[i] == max && level[i] < top_level) {
            max_ind = i;
            top_level = level[i]; // Fix the node with smallest level first
        }
    }
    realnumber min_v = 1;
    for (shortint i = 1; i <= dim; ++i) {
        if (v[i] <= min_v) {
            min_v = v[i];
        }
    }
    // If there won't be any new fixed node:
    if (min_v == 1) {
        t_data[0] = -INFINITY;
        t_data[1] = 0;
    }
    // If there will be a next fixed node:
    else {
        t_data[0] = max;
        t_data[1] = max_ind;
    }
}

/*
 This function computes the value of df(t) at t:
 */
__host__ __device__ realnumber df_exact_array(shortint father[], shortint num_kids[], shortint kids[], realnumber z_pre[], realnumber v[], shortint dim) {
    v[0] = 0;
    realnumber df = 0;
    for (shortint i = 1; i <= dim; ++i) {
        realnumber holder = z_pre[i] - z_pre[father[i]];
        for (shortint j = 0; j < num_kids[i]; ++j) {
            holder -= z_pre[kids[i * (MAX_N_NODES + 1) + j]] - z_pre[i];
        }
        df += holder * v[i];
    }
    return df;
}

/*
 This function uses all the above functions to look for solutions to df(t) = -1:
 */
__host__ __device__ void best_tree_array(realnumber value[], realnumber z[], shortint status[], shortint father[], shortint num_kids[], shortint kids[], shortint dim) {
    shortint level[MAX_N_NODES + 1] = {0};
    realnumber v[MAX_N_NODES + 1] = {1};
    realnumber z_pre[2*MAX_N_NODES + 1]; // For computing v[]
    for (shortint i = 0; i <= dim; ++i) {
        z[i] = 0;
        level[i] = 0;
        v[i] = 1;
        z_pre[i] = z[i];
    }
    v[0] = 0;
    shortint fixed_queue[MAX_N_NODES + 1];
    shortint fixed_queue_front = 0;
    shortint fixed_queue_back = MAX_N_NODES;
    shortint fixed_queue_size = 0;
    bfs_get_level(status, level, num_kids, kids, fixed_queue, &fixed_queue_front, &fixed_queue_back, &fixed_queue_size);
    // First turning point:
    realnumber n_max = value[1];
    shortint ind = 0;
    for (shortint i = 1; i <= dim; ++i) {
        if (value[i] >= n_max) {
            n_max = value[i];
        }
    }
    shortint top_level = dim + 1;
    for (shortint i = 1; i <= dim; ++i) {
        if (value[i] == n_max && level[i] < top_level) {
            ind = i;
            top_level = level[i]; // Fix the node with smallest level first
        }
    }
    status[ind] = 2;
    fixed_queue_back = (fixed_queue_back + 1) % (MAX_N_NODES + 1);
    fixed_queue[fixed_queue_back] = ind;
    ++(fixed_queue_size);
    realnumber t_now = n_max;
    realnumber t_pre = n_max;
    // Update speed by looking ahead wiht CURRENT status:
    update_z_array(t_now-delta, v, fixed_queue, &fixed_queue_front, &fixed_queue_back, &fixed_queue_size, value, z, status, level, father, num_kids, kids, dim);
    realnumber df_now = 0;
    realnumber df_pre = 0;
    // Other turning points:
    shortint next_possible_turn_current_status = 0;
    while (df_now > -1 && t_now > -INFINITY) {
        t_pre = t_now;
        df_pre = df_now;
        realnumber t_data[2] = {0};
        next_turn_array(t_data, t_pre, z_pre, v, value, status, level, dim);
        t_now = t_data[0];
        if (t_now == -INFINITY) {
            break;
        }
        ind = t_data[1];
        next_possible_turn_current_status = status[ind];
        status[ind] = 2;
        fixed_queue_back = (fixed_queue_back + 1) % (MAX_N_NODES + 1);
        fixed_queue[fixed_queue_back] = ind;
        ++(fixed_queue_size);
        // Update branching status of nodes:
        for (shortint i = 1; i <= dim; ++i) {
            if (is_branch_array(i, status, father, num_kids, kids, dim) == 1) {
                status[i] = 1;
            }
        }
        // Update z:
        update_z_array(t_now, v, fixed_queue, &fixed_queue_front, &fixed_queue_back, &fixed_queue_size, value, z, status, level, father, num_kids, kids, dim);
        for (shortint i = 1; i <= dim; ++i) {
            z_pre[i] = z[i];
        }
        df_now = df_exact_array(father, num_kids, kids, z_pre, v, dim);
        // Update speed by looking ahead wiht CURRENT status:
        update_z_array(t_now-delta, v, fixed_queue, &fixed_queue_front, &fixed_queue_back, &fixed_queue_size, value, z, status, level, father, num_kids, kids, dim);
        /*
         for (shortint i = 0; i <= dim; ++i) {
         printf("(%d, %f, %f)\n", i, z_pre[i], v[i]);
         }
         printf("________\n");
         */
    }
    // Finally get what we really want:
    if (t_now != -INFINITY) {
        // Reverse the last fixed node to its previous status:
        status[ind] = next_possible_turn_current_status;
        fixed_queue_back = (fixed_queue_back - 1) % (MAX_N_NODES + 1);
        --(fixed_queue_size);
        realnumber k = (df_pre - df_now) / (t_pre - t_now);
        realnumber t = t_pre - (df_pre + 1)  / k;
        update_z_array(t, v, fixed_queue, &fixed_queue_front, &fixed_queue_back, &fixed_queue_size, value, z, status, level, father, num_kids, kids, dim);
        z[0] = t; // Covention: z[0] = t, and z[i] = tree[i]->z for i != 0
    }
    else {
        realnumber t_last = t_pre - delta;
        update_z_array(t_last, v, fixed_queue, &fixed_queue_front, &fixed_queue_back, &fixed_queue_size, value, z, status, level, father, num_kids, kids, dim);
        for (shortint i = 1; i <= dim; ++i) {
            z_pre[i] = z[i];
        }
        realnumber df_last = df_exact_array(father, num_kids, kids, z_pre, v, dim);
        realnumber k = (df_pre - df_last) / (t_pre - t_last);
        realnumber t = (k*t_pre - df_pre - 1) / k;
        update_z_array(t, v, fixed_queue, &fixed_queue_front, &fixed_queue_back, &fixed_queue_size, value, z, status, level, father, num_kids, kids, dim);
        z[0] = t; // Covention: z[0] = t, and z[i] = tree[i]->z for i != 0
    }
}

/*
 This function uses all the above functions to look for solutions to df(t) = -1, but with slightly different constraints:
 */
__host__ __device__ void best_tree_array_inner(realnumber value[], realnumber z[], shortint status[], shortint father[], shortint num_kids[], shortint kids[], shortint dim) {
    shortint level[MAX_N_NODES + 1] = {0};
    realnumber v[MAX_N_NODES + 1] = {1};
    realnumber z_pre[2*MAX_N_NODES + 1]; // For computing v[]
    for (shortint i = 0; i <= dim; ++i) {
        z[i] = 0;
        level[i] = 0;
        v[i] = 1;
        z_pre[i] = z[i];
    }
    v[0] = 0;
    shortint fixed_queue[MAX_N_NODES + 1];
    shortint fixed_queue_front = 0;
    shortint fixed_queue_back = MAX_N_NODES;
    shortint fixed_queue_size = 0;
    bfs_get_level(status, level, num_kids, kids, fixed_queue, &fixed_queue_front, &fixed_queue_back, &fixed_queue_size);
    // First turning point:
    realnumber n_max = value[1];
    shortint ind = 0;
    for (shortint i = 1; i <= dim; ++i) {
        if (value[i] >= n_max) {
            n_max = value[i];
        }
    }
    shortint top_level = dim + 1;
    for (shortint i = 1; i <= dim; ++i) {
        if (value[i] == n_max && level[i] < top_level) {
            ind = i;
            top_level = level[i]; // Fix the node with smallest level first
        }
    }
    status[ind] = 2;
    fixed_queue_back = (fixed_queue_back + 1) % (MAX_N_NODES + 1);
    fixed_queue[fixed_queue_back] = ind;
    ++(fixed_queue_size);
    realnumber t_now = n_max;
    realnumber t_pre = n_max;
    // Update speed by looking ahead wiht CURRENT status:
    update_z_array(t_now-delta, v, fixed_queue, &fixed_queue_front, &fixed_queue_back, &fixed_queue_size, value, z, status, level, father, num_kids, kids, dim);
    realnumber df_now = 0;
    realnumber df_pre = 0;
    // Other turning points:
    shortint next_possible_turn_current_status = 0;
    while (df_now > -1 && t_now > -INFINITY) {
        t_pre = t_now;
        df_pre = df_now;
        realnumber t_data[2] = {0};
        next_turn_array(t_data, t_pre, z_pre, v, value, status, level, dim);
        t_now = t_data[0];
        if (t_now == -INFINITY) {
            break;
        }
        ind = t_data[1];
        next_possible_turn_current_status = status[ind];
        status[ind] = 2;
        fixed_queue_back = (fixed_queue_back + 1) % (MAX_N_NODES + 1);
        fixed_queue[fixed_queue_back] = ind;
        ++(fixed_queue_size);
        // Update branching status of nodes:
        for (shortint i = 1; i <= dim; ++i) {
            if (is_branch_array(i, status, father, num_kids, kids, dim) == 1) {
                status[i] = 1;
            }
        }
        // Update z:
        update_z_array(t_now, v, fixed_queue, &fixed_queue_front, &fixed_queue_back, &fixed_queue_size, value, z, status, level, father, num_kids, kids, dim);
        for (shortint i = 1; i <= dim; ++i) {
            z_pre[i] = z[i];
        }
        df_now = df_exact_array(father, num_kids, kids, z_pre, v, dim);
        if (t_now < 0) {
            break;
        }
        // Update speed by looking ahead wiht CURRENT status:
        update_z_array(t_now-delta, v, fixed_queue, &fixed_queue_front, &fixed_queue_back, &fixed_queue_size, value, z, status, level, father, num_kids, kids, dim);
    }
    // Finally get what we really want:
    if (df_now > -1 && t_now != -INFINITY) {
        // Reverse the last fixed node to its previous status:
        status[ind] = next_possible_turn_current_status;
        fixed_queue_back = (fixed_queue_back - 1) % (MAX_N_NODES + 1);
        --(fixed_queue_size);
        update_z_array(0, v, fixed_queue, &fixed_queue_front, &fixed_queue_back, &fixed_queue_size, value, z, status, level, father, num_kids, kids, dim);
        z[0] = 0; // Covention: z[0] = t, and z[i] = tree[i]->z for i != 0
    }
    else if (t_now != -INFINITY) {
        // Reverse the last fixed node to its previous status:
        status[ind] = next_possible_turn_current_status;
        fixed_queue_back = (fixed_queue_back - 1) % (MAX_N_NODES + 1);
        --(fixed_queue_size);
        realnumber k = (df_pre - df_now) / (t_pre - t_now);
        realnumber t = t_pre - (df_pre + 1)  / k;
        if (t > 0) {
            update_z_array(t, v, fixed_queue, &fixed_queue_front, &fixed_queue_back, &fixed_queue_size, value, z, status, level, father, num_kids, kids, dim);
            z[0] = t; // Covention: z[0] = t, and z[i] = tree[i]->z for i != 0
        }
        else {
            update_z_array(0, v, fixed_queue, &fixed_queue_front, &fixed_queue_back, &fixed_queue_size, value, z, status, level, father, num_kids, kids, dim);
            z[0] = 0; // Covention: z[0] = t, and z[i] = tree[i]->z for i != 0
        }
    }
    else {
        realnumber t_last = t_pre - delta;
        update_z_array(t_last, v, fixed_queue, &fixed_queue_front, &fixed_queue_back, &fixed_queue_size, value, z, status, level, father, num_kids, kids, dim);
        for (shortint i = 1; i <= dim; ++i) {
            z_pre[i] = z[i];
        }
        realnumber df_last = df_exact_array(father, num_kids, kids, z_pre, v, dim);
        realnumber k = (df_pre - df_last) / (t_pre - t_last);
        realnumber t = (k*t_pre - df_pre - 1) / k;
        if (t > 0) {
            update_z_array(t, v, fixed_queue, &fixed_queue_front, &fixed_queue_back, &fixed_queue_size, value, z, status, level, father, num_kids, kids, dim);
            z[0] = t; // Covention: z[0] = t, and z[i] = tree[i]->z for i != 0
        }
        else {
            update_z_array(0, v, fixed_queue, &fixed_queue_front, &fixed_queue_back, &fixed_queue_size, value, z, status, level, father, num_kids, kids, dim);
            z[0] = 0; // Covention: z[0] = t, and z[i] = tree[i]->z for i != 0
        }
    }
}

//
//
//
//
//
//****************************** END OF BEI's CODE WITH ARRAYS ******************************
//
//
//







/*
 This function uses DFS to compute the father nodes of each node, and store the result in fathers_list[]:
 */
__host__ __device__ void dfs_tree_compute_fathers_non_recursive_array(shortint num_nodes,shortint *fathers_list, shortint *adj, shortint *deg, shortint root, shortint *stack, shortint *visited){
    fathers_list[0] = -1;
    shortint curr = root;
    shortint stack_depth = 0;
    stack[stack_depth] = curr;
    stack_depth = stack_depth + 1;
    fathers_list[root + 1] = 0;
    
    while(stack_depth > 0){
        stack_depth = stack_depth - 1;
        curr = stack[stack_depth];
        visited[curr] = 1 - visited[curr];
        for (shortint j = 0; j < deg[curr]; j++){
            shortint child_ix = adj[curr*num_nodes + j];
            if (visited[child_ix] == 1 - visited[root]){
                stack[stack_depth] = child_ix;
                stack_depth = stack_depth + 1;
                fathers_list[child_ix + 1] = curr + 1;
            }
        }
    }
}

/*
 This function uses DFS to compute ntilde from tree structure and input data:
 */
__host__ __device__ void dfs_tree_compute_ntilde_non_recursive_array(shortint num_nodes, shortint *adj, shortint *deg, shortint root, shortint *stack, shortint *visited, realnumber *ntilde, realnumber *data,  shortint T, shortint t){
    ntilde[0] = 0;
    shortint curr = root;
    shortint stack_depth = 0;
    stack[stack_depth] = curr;
    stack_depth = stack_depth + 1;
    ntilde[curr + 1] = data[curr*T + t];
    
    while(stack_depth > 0){
        stack_depth = stack_depth - 1;
        curr = stack[stack_depth];
        visited[curr] = 1 - visited[curr];
        for (shortint j = 0; j < deg[curr]; j++){
            shortint child_ix = adj[curr*num_nodes + j];
            if (visited[child_ix] == 1 - visited[root]){
                stack[stack_depth] = child_ix;
                stack_depth = stack_depth + 1;
                ntilde[child_ix + 1] = ntilde[curr + 1] + data[child_ix*T + t];
            }
        }
    }
}

/*
 This function converts tree structure data into the format needed for cost1 and cost2
 */
__host__ __device__ void convert_tree_data(shortint num_nodes, shortint *final_degrees, shortint *adj_list, shortint *father_list, shortint root_node, shortint num_kids[], shortint kids[]) {
    // Virtural node:
    num_kids[0] = 1;
    kids[0] = root_node + 1;
    // num_nodes:
    for (int i = 1; i <= num_nodes; ++i) {
        if (i == root_node + 1) {
            num_kids[i] = final_degrees[i - 1];
        }
        else {
            num_kids[i] = final_degrees[i - 1] - 1;
        }
    }
    // kids:
    for (shortint i = 1; i <= num_nodes; ++i) {
        if (i == root_node + 1) {
            for (shortint j = 0; j <= num_kids[i]; ++j) {
                // Add one extra spot to accomadate possible effective node
                // Default to -2 indicating no effective node present
                if (j == num_kids[i]) {
                    kids[i * (MAX_N_NODES + 1) + num_kids[i]] = -2;
                }
                else {
                    kids[i * (MAX_N_NODES + 1) + j] = adj_list[(i-1)*num_nodes + j] + 1;
                }
            }
        }
        else {
            shortint father_ind = father_list[i];
            shortint j = 0;
            while (j < num_kids[i] && adj_list[(i-1)*num_nodes + j] + 1 != father_ind) {
                kids[i * (MAX_N_NODES + 1) + j] = adj_list[(i-1)*num_nodes + j] + 1;
                ++j;
            }
            while (j < num_kids[i]) {
                kids[i * (MAX_N_NODES + 1) + j] = adj_list[(i-1)*num_nodes + j + 1] + 1;
                ++j;
            }
            // Add one extra spot to accomadate possible effective node
            // Default to -2 indicating no effective node present
            if (j == num_kids[i]) {
                kids[i * (MAX_N_NODES + 1) + num_kids[i]] = -2;
            }
        }
    }
}

/*
 This function uses DFS to compute the entire cost of a given tree at a single time t
 */
__host__ __device__ realnumber dfs_tree_cost_from_z_non_recursive(shortint num_nodes, shortint *adj_list, shortint *final_degrees, shortint root_node, shortint *stack, shortint *visited, realnumber *z, realnumber *data, shortint T, shortint t) {
    
    shortint curr = root_node;
    shortint stack_depth = 0;
    stack[stack_depth] = curr;
    stack_depth = stack_depth + 1;
    // Root:
    realnumber temp2 = (z[curr + 1] + data[curr*T + t]);
    realnumber temp = (temp2*temp2);
    // Other nodes:
    while(stack_depth > 0){
        stack_depth = stack_depth - 1;
        curr = stack[stack_depth];
        visited[curr] = 1 - visited[curr];
        for (shortint j = 0; j < final_degrees[curr]; j++){
            shortint child_ix = adj_list[curr*num_nodes + j];
            if (visited[child_ix] == 1 - visited[root_node]){
                stack[stack_depth] = child_ix;
                stack_depth = stack_depth + 1;
                temp2 = ((z[child_ix + 1] - z[curr + 1]) + data[child_ix*T + t]);
                temp +=  (temp2*temp2);
            }
        }
    }
    return temp;
}

/*
 This function uses DFS to compute the entire cost of a given tree over all time t in [0, T)
 */
__host__ __device__ realnumber tree_cost_bei_array(shortint num_nodes, shortint T, realnumber *data, shortint root_node, edge *tree, shortint *adjacency_mat, shortint *final_degrees, shortint *adj_list){
    
    realnumber z[2*MAX_N_NODES + 1];
    realnumber ntilde[2*MAX_N_NODES + 1];
    shortint stack[MAX_N_NODES];
    shortint visited[MAX_N_NODES];
    shortint father_list[2*MAX_N_NODES + 1];
    shortint num_kids[MAX_N_NODES + 1];
    shortint kids[(MAX_N_NODES + 1) * (MAX_N_NODES + 1)];
    
    // we need to set visited to be all zeros. We only need to do this once.
    // the rest does not need to be initialized
    for (shortint i = 0; i < num_nodes; i++){
        visited[i] = 0;
    }
    
    realnumber error_tree_model = 0.0;
    
    dfs_tree_compute_fathers_non_recursive_array(num_nodes, father_list, adj_list, final_degrees, root_node, stack, visited);
    
    convert_tree_data(num_nodes, final_degrees, adj_list, father_list, root_node, num_kids, kids); //maybe remove in the future
    
    // go over all the nodes and compute f_i - sum_children f_j
    for (shortint t = 0; t < T; ++t) {
        // Reset status of all nodes to be free for the next iteration:
        shortint status[2*MAX_N_NODES + 1] = {0};
        status[0] = 2;
        
        // we explore the tree using BFS to compute tiln = U'*F;
        dfs_tree_compute_ntilde_non_recursive_array(num_nodes, adj_list, final_degrees, root_node, stack, visited, ntilde, data, T, t);
        
        // Compute z_i for all i:
        best_tree_array(ntilde, z, status, father_list, num_kids, kids, num_nodes);
        
        // Compute square_norm((z_i - z_i_father) + F_i):
        realnumber tmp_cost = dfs_tree_cost_from_z_non_recursive(num_nodes, adj_list, final_degrees, root_node, stack, visited, z, data, T, t);
        
        error_tree_model += tmp_cost;
    }
    
    return error_tree_model;
}

__host__ __device__ realnumber tree_cost_bei_array_inner(shortint num_nodes, shortint T, realnumber *data, shortint root_node, edge *tree, shortint *adjacency_mat, shortint *final_degrees, shortint *adj_list){
    
    realnumber z[2*MAX_N_NODES + 1];
    realnumber ntilde[2*MAX_N_NODES + 1];
    shortint stack[MAX_N_NODES];
    shortint visited[MAX_N_NODES];
    shortint father_list[2*MAX_N_NODES + 1];
    shortint num_kids[MAX_N_NODES + 1];
    shortint kids[(MAX_N_NODES + 1) * (MAX_N_NODES + 1)];
    
    // we need to set visited to be all zeros. We only need to do this once.
    // the rest does not need to be initialized
    for (shortint i = 0; i < num_nodes; i++){
        visited[i] = 0;
    }
    
    realnumber error_tree_model = 0.0;
    
    dfs_tree_compute_fathers_non_recursive_array(num_nodes, father_list, adj_list, final_degrees, root_node, stack, visited);
    
    convert_tree_data(num_nodes, final_degrees, adj_list, father_list, root_node, num_kids, kids); //maybe remove in the future
    
    // go over all the nodes and compute f_i - sum_children f_j
    for (shortint t = 0; t < T; ++t) {
        // Reset status of all nodes to be free for the next iteration:
        shortint status[2*MAX_N_NODES + 1] = {0};
        status[0] = 2;
        
        // we explore the tree using BFS to compute tiln = U'*F;
        dfs_tree_compute_ntilde_non_recursive_array(num_nodes, adj_list, final_degrees, root_node, stack, visited, ntilde, data, T, t);
        
        // Compute z_i for all i:
        best_tree_array_inner(ntilde, z, status, father_list, num_kids, kids, num_nodes);
        
        // Compute square_norm((z_i - z_i_father) + F_i):
        realnumber tmp_cost = dfs_tree_cost_from_z_non_recursive(num_nodes, adj_list, final_degrees, root_node, stack, visited, z, data, T, t);
        
        error_tree_model += tmp_cost;
    }
    
    return error_tree_model;
}

/*
 Update and maintain a max heap to store the smallest k costs and keep an array of corresponding tree indices:
 */
__host__ __device__ void heapify(realnumber cost, longint tree_index, realnumber max_heap[], longint smallest_trees[], shortint k) {
    if (cost >= max_heap[0]) {
        return;
    }
    else {
        // Remove the original max:
        realnumber temp = max_heap[k - 1];
        max_heap[0] = temp;
        longint temp_ind = smallest_trees[k - 1];
        smallest_trees[0] = temp_ind;
        shortint current_ind = 0;
        // Bubble down:
        while (1) {
            shortint left_ind = 2*current_ind + 1;
            shortint right_ind = 2*current_ind + 2;
            // Make sure we are reaching within boundary:
            if (left_ind < k && right_ind < k) {
                if (temp < max_heap[left_ind] && temp < max_heap[right_ind]) {
                    if (max_heap[left_ind] >= max_heap[right_ind]) {
                        max_heap[current_ind] = max_heap[left_ind];
                        max_heap[left_ind] = temp;
                        smallest_trees[current_ind] = smallest_trees[left_ind];
                        smallest_trees[left_ind] = temp_ind;
                        current_ind = left_ind;
                    }
                    else {
                        max_heap[current_ind] = max_heap[right_ind];
                        max_heap[right_ind] = temp;
                        smallest_trees[current_ind] = smallest_trees[right_ind];
                        smallest_trees[right_ind] = temp_ind;
                        current_ind = right_ind;
                    }
                }
                else if (temp < max_heap[left_ind] && temp >= max_heap[right_ind]) {
                    max_heap[current_ind] = max_heap[left_ind];
                    max_heap[left_ind] = temp;
                    smallest_trees[current_ind] = smallest_trees[left_ind];
                    smallest_trees[left_ind] = temp_ind;
                    current_ind = left_ind;
                }
                else if (temp < max_heap[right_ind] && temp >= max_heap[left_ind]) {
                    max_heap[current_ind] = max_heap[right_ind];
                    max_heap[right_ind] = temp;
                    smallest_trees[current_ind] = smallest_trees[right_ind];
                    smallest_trees[right_ind] = temp_ind;
                    current_ind = right_ind;
                }
                else {
                    break;
                }
            }
            else if (left_ind < k && temp < max_heap[left_ind]) {
                max_heap[current_ind] = max_heap[left_ind];
                max_heap[left_ind] = temp;
                smallest_trees[current_ind] = smallest_trees[left_ind];
                smallest_trees[left_ind] = temp_ind;
                current_ind = left_ind;
            }
            else if (right_ind < k && temp < max_heap[right_ind]) {
                max_heap[current_ind] = max_heap[right_ind];
                max_heap[right_ind] = temp;
                smallest_trees[current_ind] = smallest_trees[right_ind];
                smallest_trees[right_ind] = temp_ind;
                current_ind = right_ind;
            }
            else {
                break;
            }
        }
        // Put cost at the end of array:
        max_heap[k - 1] = cost;
        smallest_trees[k - 1] = tree_index;
        current_ind = k - 1;
        // Bubble up:
        while (1) {
            shortint father_ind = (current_ind - 1) / 2; // Integer division acts like floor function
            if (cost > max_heap[father_ind]) {
                max_heap[current_ind] = max_heap[father_ind];
                max_heap[father_ind] = cost;
                smallest_trees[current_ind] = smallest_trees[father_ind];
                smallest_trees[father_ind] = tree_index;
                current_ind = father_ind;
            }
            else {
                break;
            }
        }
    }
}

//
//
//
//
//
//****************************** END OF BEI's CODE ******************************
//
//
//












longint pow_int(shortint n, shortint p){
    longint r = 1;
    for (shortint i = 1; i <= p; i++){
        r = r*n;
    }
    return r;
}

// maybe this will have to be implemented in a non recursive manner
// depth-first-search to form the vector U for a given time slice
// the notion of "visited" is relative so that we do not need to reinitialize the vector visited each time we call bfs_tree
// the bfs_tree will induce a direction on the tree starting at the root node
__host__ __device__ void dfs_tree_cost(shortint num_nodes, shortint * adj, shortint * deg, shortint curr, shortint * visited, realnumber * U, realnumber * data,  shortint T, shortint t){
    
    visited[curr] = 1 - visited[curr]; // this marks the node as visited
    
    U[curr] = data[curr*T + t];
    
    for (shortint j = 0; j < deg[curr]; j++){
        shortint child_ix = adj[curr*num_nodes + j];
        if ( visited[child_ix] == 1 - visited[curr] ){ // we compare the current status with the status of the node we just visited
            
            U[curr] = U[curr] - data[child_ix*T + t];
            
            dfs_tree_cost(num_nodes, adj, deg, child_ix, visited, U, data,T , t);
        }
        
    }
    
    
}

// this is a non-recursive version of the depth-first-search algorithm
// maybe the GPU code will work better with non recursive algorithms
__host__ __device__ void dfs_tree_cost_non_recursive(shortint num_nodes, shortint * adj, shortint * deg, shortint root, shortint * stack, shortint * visited, realnumber * U, realnumber * data,  shortint T, shortint t){
    
    shortint curr = root;
    
    shortint stack_depth = 0;
    
    stack[stack_depth] = curr;
    stack_depth = stack_depth + 1;
    
    //printf("LIST: ");
    while(stack_depth > 0){
        
        stack_depth = stack_depth - 1;
        curr = stack[stack_depth];
        visited[curr] = 1 - visited[curr];
        
        //printf("(%d, %d) ", curr,deg[curr]);
        
        U[curr] = data[curr*T + t];
        
        for (shortint j = 0; j < deg[curr]; j++){
            shortint child_ix = adj[curr*num_nodes + j];
            //printf(":%d:%d ",child_ix,visited[child_ix]);
            if (visited[child_ix] == 1 - visited[root]){
                stack[stack_depth] = child_ix;
                stack_depth = stack_depth + 1;
                
                U[curr] = U[curr] - data[child_ix*T + t];
            }
        }
    }
    //printf("\n");
}


// this projects an element into the simplex where vectors are non-negative and must sum to one
// the algorithm runs in N unsigned N time
// the algorithm requires temporary arrays nsorted and cumsum. These arrays do not need to be initialized
// the vectors n and x need to be allocated but only n needs to be initialized
// if nooutput == 1 then we ignore x and simply ouput the error
// we return the distance squared because we will have to sum this distance with other distances when we are dealing with multiple time instants
__host__ __device__ realnumber projection_onto_simplex(shortint nooutput, realnumber *x, realnumber *n, shortint numvars, realnumber *nsorted, realnumber *cumsum){
    
    realnumber distance_squared = 0;
    
    for (shortint i = 0; i < numvars; i++){
        nsorted[i] = -n[ i ];
    }
    //sort
    sort_using_networks(nsorted, numvars);
    
    cumsum[0] = -nsorted[0];
    
    for (shortint i = 1; i < numvars; i++){
        cumsum[i] = (-nsorted[i]) + cumsum[i-1];
    }
    //find switch ix
    shortint ix = 0;
    
    realnumber maxval = 100000000;
    for (shortint i = 0; i < numvars; i++){
        realnumber tmp = (-nsorted[i]) + (1.0/ (1.0 + (realnumber) i))*(1 - cumsum[i]);
        
        if (tmp > 0 && tmp < maxval){
            maxval = tmp;
            ix = i;
        }
    }
    //find lambda
    
    realnumber lambda =  (1.0/ (1.0 + (realnumber) ix))*(1 - cumsum[ix]);
    
    //create final x
    for (shortint i = 0; i < numvars; i++){
        realnumber tmp = 	n[i] + lambda;
        
        realnumber xi;
        if (tmp<0){
            xi = 0;
        }else{
            xi = tmp;
        }
        
        distance_squared = distance_squared + (xi - n[i])*(xi - n[i]);
        
        // we might not be interested in the output in which we case we skip this part
        if (nooutput == 0){
            x[i] = xi;
        }
        
    }
    
    return distance_squared;
    
}

// this does the projection according to Condats' algorithm which is substantially faster than ours
__device__ __host__ realnumber projection_onto_simplex_condat(shortint nooutput, realnumber* x,  realnumber* y, const unsigned int length) {

	realnumber distance_squared = 0;
	const realnumber a = 1;    
	realnumber*    aux = x;
    realnumber*  aux0=aux;
    int        auxlength=1;
    int        auxlengthold=-1;
    realnumber    tau=(*aux=*y)-a;
    int     i=1;
    for (; i<length; i++)
        if (y[i]>tau) {
            if ((tau+=((aux[auxlength]=y[i])-tau)/(auxlength-auxlengthold))
                <=y[i]-a) {
                tau=y[i]-a;
                auxlengthold=auxlength-1;
            }
            auxlength++;
        }
    if (auxlengthold>=0) {
        auxlength-=++auxlengthold;
        aux+=auxlengthold;
        while (--auxlengthold>=0)
            if (aux0[auxlengthold]>tau)
                tau+=((*(--aux)=aux0[auxlengthold])-tau)/(++auxlength);
    }
    do {
        auxlengthold=auxlength-1;
        for (i=auxlength=0; i<=auxlengthold; i++)
            if (aux[i]>tau)
                aux[auxlength++]=aux[i];
            else
                tau+=(tau-aux[i])/(auxlengthold-i+auxlength);
    } while (auxlength<=auxlengthold);
    
	for (i=0; i<length; i++){
		realnumber xi = (y[i]>tau ? y[i]-tau : 0.0);

	 	distance_squared = distance_squared + (xi - y[i])*(xi - y[i]);

		if (nooutput == 0){
			x[i]=xi;		
		}
	    
	}
	
	return distance_squared;

}


// this is just like the function above by the restriction sum(x) == 1 is now replaced by sum(x) <= 1
__host__ __device__ realnumber projection_onto_inner_simplex(shortint nooutput, realnumber *x, realnumber *n, shortint numvars, realnumber *nsorted, realnumber *cumsum){
    
    realnumber distance_squared = 0;
    
    for (shortint i = 0; i < numvars; i++){
        nsorted[i] = -n[ i ];
    }
    //sort
    sort_using_networks(nsorted, numvars);
    
    cumsum[0] = -nsorted[0];
    
    for (shortint i = 1; i < numvars; i++){
        cumsum[i] = (-nsorted[i]) + cumsum[i-1];
    }
    //find switch ix
    shortint ix = 0;
    
    realnumber maxval = 100000000;
    for (shortint i = 0; i < numvars; i++){
        realnumber tmp = (-nsorted[i]) + (1.0/ (1.0 + (realnumber) i))*(1 - cumsum[i]);
        
        if (tmp > 0 && tmp < maxval){
            maxval = tmp;
            ix = i;
        }
    }
    //find lambda
    realnumber lambda =  (1.0/ (1.0 + (realnumber) ix))*(1 - cumsum[ix]);
    
    // threshold because we now have sum(x) <=1 instead of sum(x) == 1
    if (lambda > 0){
        lambda = 0;
    }
    
    //create final x
    for (shortint i = 0; i < numvars; i++){
        realnumber tmp = 	n[i] + lambda;
        
        realnumber xi;
        if (tmp<0){
            xi = 0;
        }else{
            xi = tmp;
        }
        
        distance_squared = distance_squared + (xi - n[i])*(xi - n[i]);
        
        // we might not be interested in the output in which we case we skip this part
        if (nooutput == 0){
            x[i] = xi;
        }
        
    }
    
    return distance_squared;
    
}





__device__ __host__ realnumber projection_onto_inner_simplex_condat(shortint nooutput, realnumber* x,  realnumber* y, const unsigned int length) {

	realnumber distance_squared = 0;
	const realnumber a = 1;    
	realnumber*    aux = x;
    realnumber*  aux0=aux;
    int        auxlength=1;
    int        auxlengthold=-1;
    realnumber    tau=(*aux=*y)-a;
    int     i=1;
    for (; i<length; i++)
        if (y[i]>tau) {
            if ((tau+=((aux[auxlength]=y[i])-tau)/(auxlength-auxlengthold))
                <=y[i]-a) {
                tau=y[i]-a;
                auxlengthold=auxlength-1;
            }
            auxlength++;
        }
    if (auxlengthold>=0) {
        auxlength-=++auxlengthold;
        aux+=auxlengthold;
        while (--auxlengthold>=0)
            if (aux0[auxlengthold]>tau)
                tau+=((*(--aux)=aux0[auxlengthold])-tau)/(++auxlength);
    }
    do {
        auxlengthold=auxlength-1;
        for (i=auxlength=0; i<=auxlengthold; i++)
            if (aux[i]>tau)
                aux[auxlength++]=aux[i];
            else
                tau+=(tau-aux[i])/(auxlengthold-i+auxlength);
    } while (auxlength<=auxlengthold);
    
	if (tau < 0){
		tau = 0;	
	}

	for (i=0; i<length; i++){
		realnumber xi = (y[i]>tau ? y[i]-tau : 0.0);

	 	distance_squared = distance_squared + (xi - y[i])*(xi - y[i]);

		if (nooutput == 0){
			x[i]=xi;		
		}
	    
	}
	
	return distance_squared;

}



// here we will write several functions to evaluate the quality of the trees being generated
// all the functrions receive a liset of edges, an adjacency matrix, a list of degrees and an adjacancy list
// however, they might not make use of all of these
// the functions also receive an array that contains, per row, the evolution of a mutate position aunsigned time
// the function might make use of a scrape memory that needs to be initialized before hand, appropriately
// data is a matrix with num_nodes rows and T columns
__host__ __device__ realnumber tree_cost(shortint num_nodes, shortint T, realnumber * data, shortint root_node, edge *tree, shortint *adjacency_mat , shortint * final_degrees, shortint *adj_list){
    
    // extract from the scrap memory the different components we want to make use of
    //realnumber *nsorted = (realnumber *) scrapmem;
    //realnumber *cumsum = &(nsorted[num_nodes]);
    //realnumber *U_vec = &(cumsum[num_nodes]);
    //shortint *visited = (shortint *) &(U_vec[num_nodes]);
    //shortint *stack = (shortint *) &(visited[num_nodes]);
    
    realnumber nsorted[MAX_N_NODES];
    //realnumber cumsum[MAX_N_NODES];
    realnumber U_vec[MAX_N_NODES];
    shortint visited[MAX_N_NODES];
    shortint stack[MAX_N_NODES];
    
    // we need to set visited to be all zeros. We only need to do this once.
    // the rest does not need to be initialized
    for (int i = 0; i < num_nodes; i++){
        visited[i] = 0;
    }
    
    realnumber error_tree_model = 0;
    
    // go over all the nodes and compute f_i - sum_children f_j
    for (shortint t = 0; t < T; t++){
        
        // we explore the tree using BFS
        //dfs_tree_cost(num_nodes, adj_list, final_degrees, root_node, visited, U_vec, data,  T, t);
        dfs_tree_cost_non_recursive(num_nodes,adj_list, final_degrees,  root_node, stack, visited, U_vec, data,   T,  t);
        //printf(" U VECT\n");
        
        //for (int i = 0; i < num_nodes ; i++){
        //    printf("%f ",U_vec[i]);
        //}
        //printf("\n");
        
        
        // here the variable U has, for every point in time, a measure of the error in the tree model
        // note that it is important for the projection funtion to return the distance squared so that we can add these up
        //error_tree_model = error_tree_model + projection_onto_simplex(1, NULL, U_vec, num_nodes, nsorted, cumsum);
        error_tree_model = error_tree_model + projection_onto_simplex_condat(1, nsorted, U_vec, num_nodes);        

    }
    
    return error_tree_model;
}

__host__ __device__ realnumber tree_cost_inner(shortint num_nodes, shortint T, realnumber * data, shortint root_node, edge *tree, shortint *adjacency_mat , shortint * final_degrees, shortint *adj_list){
    
    // extract from the scrap memory the different components we want to make use of
    //realnumber *nsorted = (realnumber *) scrapmem;
    //realnumber *cumsum = &(nsorted[num_nodes]);
    //realnumber *U_vec = &(cumsum[num_nodes]);
    //shortint *visited = (shortint *) &(U_vec[num_nodes]);
    //shortint *stack = (shortint *) &(visited[num_nodes]);
    
    realnumber nsorted[MAX_N_NODES];
    //realnumber cumsum[MAX_N_NODES];
    realnumber U_vec[MAX_N_NODES];
    shortint visited[MAX_N_NODES];
    shortint stack[MAX_N_NODES];
    
    // we need to set visited to be all zeros. We only need to do this once.
    // the rest does not need to be initialized
    for (int i = 0; i < num_nodes; i++){
        visited[i] = 0;
    }
    
    realnumber error_tree_model = 0;
    
    // go over all the nodes and compute f_i - sum_children f_j
    for (shortint t = 0; t < T; t++){
        
        // we explore the tree using BFS
        //dfs_tree_cost(num_nodes, adj_list, final_degrees, root_node, visited, U_vec, data,  T, t);
        dfs_tree_cost_non_recursive(num_nodes,adj_list, final_degrees,  root_node, stack, visited, U_vec, data,   T,  t);
        //printf(" U VECT\n");
        
        //for (int i = 0; i < num_nodes ; i++){
        //    printf("%f ",U_vec[i]);
        //}
        //printf("\n");
        
        
        // here the variable U has, for every point in time, a measure of the error in the tree model
        // note that it is important for the projection funtion to return the distance squared so that we can add these up
        //error_tree_model = error_tree_model + projection_onto_inner_simplex(1, NULL, U_vec, num_nodes, nsorted, cumsum);
		error_tree_model = error_tree_model + projection_onto_inner_simplex_condat(1, nsorted, U_vec, num_nodes);
        
    }
    
    return error_tree_model;
}


// The perm has size n-2
// Each element in perm has numbers between 0 and num_tree_vertices - 1
__host__ __device__ void permutation(shortint num_tree_vertices, longint perm_index, shortint *perm)
{
    for (shortint k = num_tree_vertices - 2; k >= 1; --k )
    {
        perm[k - 1] = perm_index % num_tree_vertices;
        perm_index = perm_index / num_tree_vertices;
    }
}

// Degrees is a scrap memory position. In the end it will not contain anything useful
// The degrees arraw does not need to be initialized
// We assume that the smallest element of "code" starts at 0
// note that none of the arrays needs to be pre-initialized to any values. only the memory needs to be available

__host__ __device__ void prufer_tree(edge *tree, shortint *degrees, shortint *code, shortint num_tree_vertices)
{
    shortint i, j;
    
    shortint last_edge_ix = 0; // this is next free index to use for edges
    
    shortint code_len = num_tree_vertices - 2;
    
    
    // Start with 1 in all of them
    for (i = 0; i < num_tree_vertices; i++) {
        degrees[i] = 1;
    }
    
    // Add to the degrees the number of occurrences of each node index in the code
    for (i = 0; i < code_len; i++) {
        degrees[code[i]] = degrees[code[i]] + 1;
    }
    
    // Add edges to nodes in the code
    for (i = 0; i < code_len; i++) {
        // Find the lowest-numbered node with degree 1
        // note that the following codes requires j to be global withing this function
        for (j = 0; degrees[j] != 1; j++);
        
        // Add the edge
        tree[last_edge_ix].first = j; // Note that these indices start at zero
        tree[last_edge_ix].second = code[i]; // Note that these indices start at zero
        
        degrees[j] = degrees[j] - 1;
        degrees[code[i]] = degrees[code[i]] - 1;
        last_edge_ix = last_edge_ix + 1;
    }
    /* Find the last 2 degree-1 nodes */
    for (i = 0; degrees[i] != 1; i++);
    for (j = i + 1; degrees[j] != 1; j++);
    
    /* Add the last edge */
    tree[last_edge_ix].first = i;
    tree[last_edge_ix].second = j;
}

// the memory "integer_array" will never change
// degrees is a scrap memory position. in the end it will not contain anything useful
// regarding this function the permutation_array will also be scrap memory
// note that none of the arrays needs to be pre-initialized to any values. only the memory needs to be available
__host__ __device__ void generate_tree_from_index(shortint num_tree_vertices, longint tree_index, shortint *permutation_array,  shortint *degrees, edge *tree, shortint *adjacency_mat , shortint * final_degrees, shortint *adj_list)
{
    // General purpose counters
    shortint j = 0;
    
    shortint source; //source node for an edge
    shortint destination;
    
    // Getting permutation from function
    permutation( num_tree_vertices, tree_index, permutation_array );
    
    // Function to provide Prufer_code tree structure
    prufer_tree( tree, degrees, permutation_array, num_tree_vertices );
    
    // Filling up the adjacency matrix
    // Note that a tree with   num_tree_vertices    vertices has (num_tree_vertices - 1) edges
    
    for (j = 0; j < num_tree_vertices ; j++){
        final_degrees[ j ] = 0;
    }
    for (j = 0; j < num_tree_vertices*num_tree_vertices ; j++){
        adjacency_mat[ j ] = 0;
    }
    
    for ( j = 0; j < num_tree_vertices - 1; j++ )
    {
        source = tree[j].first;
        destination = tree[j].second;
        // use 1D representation for the adjacency matrix
        adjacency_mat[   num_tree_vertices*source   +   destination] = 1;
        adjacency_mat[   num_tree_vertices*destination   +   source] = 1;
        
        // use a 1D representation for hte adjacency list that is over initialized since we have enough space
        adj_list[source*num_tree_vertices + final_degrees[source]] = destination;
        adj_list[destination*num_tree_vertices + final_degrees[destination]] = source;
        final_degrees[source] = final_degrees[source] + 1;
        final_degrees[destination] = final_degrees[destination] + 1;
        
    }
    
}

// CUDA kernels with top k, one for each cost fucntion:
__global__ void kernel_to_compute_optimal_tree_only_global_mem_1k(int chunck_per_cycle, shortint num_devices, shortint device_index, longint total_num_trees, shortint num_tree_vertices, shortint  T, shortint root_node, realnumber *data, void *scrape_memory, shortint shift_size, realnumber *max_heap_device, longint *smallest_trees_device, int k){
    
    int baseix = blockIdx.x*blockDim.x + threadIdx.x;

    char * shifted_pointer = (char *) scrape_memory;
    shifted_pointer = shifted_pointer + shift_size*baseix;
    
    // Allocate top k stuff:
    realnumber max_heap[topk];
    longint smallest_trees[topk];
    for (int i = 0; i < topk; ++i) {
        max_heap[i] = topk * MAX_N_NODES + i;
        smallest_trees[i] = 0;
    }
    
    // we only use global memory here
    shortint * adjacency_mat = (shortint *) shifted_pointer;
    shortint * adjacency_list = (shortint *) &(adjacency_mat[num_tree_vertices*num_tree_vertices]);
    edge *tree = (edge *) &(adjacency_list[num_tree_vertices*num_tree_vertices]);
    shortint * final_degrees = (shortint *) &(tree[num_tree_vertices-1]);
    shortint * perm_scrape = (shortint*) &(final_degrees[num_tree_vertices-2]);
    shortint * deg_scrape = (shortint*) &(perm_scrape[num_tree_vertices]);
    
    for (long int ix = device_index  +  baseix*num_devices  ; ix < total_num_trees ; ix = ix + num_devices*chunck_per_cycle){
        
        generate_tree_from_index( num_tree_vertices, ix, perm_scrape, deg_scrape, tree, adjacency_mat , final_degrees, adjacency_list );
        
        realnumber val  = tree_prior(num_tree_vertices, adjacency_mat, adjacency_list, final_degrees); //here we include a prior cost on the tree  
		// if the prior says that the tree will topology is not allowed, we do not even compute the rest 			
		if (val < 10000){
			val +=  tree_cost_bei_array(num_tree_vertices, T, data, root_node, tree, adjacency_mat, final_degrees, adjacency_list);			
		}

        // Push into max_heap:
        if (val == val && val < max_heap[0]) { // Avoid NaN
            heapify(val, ix, max_heap, smallest_trees, k);
        }
        
    }
    
    // Write top k stuff to global memory:
    for (shortint i = 0; i < k; ++i) {
        max_heap_device[baseix + i * chunck_per_cycle] = max_heap[i];
        smallest_trees_device[baseix + i * chunck_per_cycle] = smallest_trees[i];
    }
}

__global__ void kernel_to_compute_optimal_tree_only_global_mem_2k(int chunck_per_cycle, shortint num_devices, shortint device_index, longint total_num_trees, shortint num_tree_vertices, shortint  T, shortint root_node, realnumber *data, void *scrape_memory , shortint shift_size, realnumber *max_heap_device, longint *smallest_trees_device, int k){
    
    int baseix = blockIdx.x*blockDim.x + threadIdx.x;
    
    char * shifted_pointer = (char *) scrape_memory;
    shifted_pointer = shifted_pointer + shift_size*baseix;
    
    // Allocate top k stuff:
    realnumber max_heap[topk];
    longint smallest_trees[topk];
    for (int i = 0; i < topk; ++i) {
        max_heap[i] = topk * MAX_N_NODES + i;
        smallest_trees[i] = 0;
    }
    
    // we only use global memory here
    shortint * adjacency_mat = (shortint *) shifted_pointer;
    shortint * adjacency_list = (shortint *) &(adjacency_mat[num_tree_vertices*num_tree_vertices]);
    edge *tree = (edge *) &(adjacency_list[num_tree_vertices*num_tree_vertices]);
    shortint * final_degrees = (shortint *) &(tree[num_tree_vertices-1]);
    shortint * perm_scrape = (shortint*) &(final_degrees[num_tree_vertices-2]);
    shortint * deg_scrape = (shortint*) &(perm_scrape[num_tree_vertices]);
    
    for (long int ix = device_index  +  baseix*num_devices  ; ix < total_num_trees ; ix = ix + num_devices*chunck_per_cycle){
            
        generate_tree_from_index( num_tree_vertices, ix, perm_scrape, deg_scrape, tree, adjacency_mat , final_degrees, adjacency_list );
        
		realnumber val  = tree_prior(num_tree_vertices, adjacency_mat, adjacency_list, final_degrees); //here we include a prior cost on the tree  
		// if the prior says that the tree will topology is not allowed, we do not even compute the rest 			
		if (val < 10000){
			val +=  tree_cost_bei_array_inner(num_tree_vertices, T, data, root_node, tree, adjacency_mat, final_degrees, adjacency_list);		
		}

        // Push into max_heap:
        if (val == val && val < max_heap[0]) { // Avoid NaN
            heapify(val, ix, max_heap, smallest_trees, k);
        }

    }
    
    // Write top k stuff to global memory:
    for (shortint i = 0; i < k; ++i) {
        max_heap_device[baseix + i * chunck_per_cycle] = max_heap[i];
        smallest_trees_device[baseix + i * chunck_per_cycle] = smallest_trees[i];
    }
}

__global__ void kernel_to_compute_optimal_tree_only_global_mem_3k(int chunck_per_cycle,shortint num_devices, shortint device_index,  longint total_num_trees, shortint num_tree_vertices, shortint  T, shortint root_node, realnumber *data, void *scrape_memory , shortint shift_size, realnumber *max_heap_device, longint *smallest_trees_device, int k){
    
    int baseix = blockIdx.x*blockDim.x + threadIdx.x;
    
    char * shifted_pointer = (char *) scrape_memory;
    shifted_pointer = shifted_pointer + shift_size*baseix;
    
    // Allocate top k stuff:
    realnumber max_heap[topk];
    longint smallest_trees[topk];
    for (int i = 0; i < topk; ++i) {
        max_heap[i] = topk * MAX_N_NODES + i;
        smallest_trees[i] = 0;
    }
    
    // we only use global memory here
    shortint * adjacency_mat = (shortint *) shifted_pointer;
    shortint * adjacency_list = (shortint *) &(adjacency_mat[num_tree_vertices*num_tree_vertices]);
    edge *tree = (edge *) &(adjacency_list[num_tree_vertices*num_tree_vertices]);
    shortint * final_degrees = (shortint *) &(tree[num_tree_vertices-1]);
    shortint * perm_scrape = (shortint*) &(final_degrees[num_tree_vertices-2]);
    shortint * deg_scrape = (shortint*) &(perm_scrape[num_tree_vertices]);
    
    for (long int ix = device_index  +  baseix*num_devices  ; ix < total_num_trees ; ix = ix + num_devices*chunck_per_cycle){
            
        generate_tree_from_index( num_tree_vertices, ix, perm_scrape, deg_scrape, tree, adjacency_mat , final_degrees, adjacency_list );
        
		realnumber val  = tree_prior(num_tree_vertices, adjacency_mat, adjacency_list, final_degrees); //here we include a prior cost on the tree  
		// if the prior says that the tree will topology is not allowed, we do not even compute the rest 			
		if (val < 10000){
			val += tree_cost(num_tree_vertices, T, data, root_node, tree, adjacency_mat, final_degrees, adjacency_list);	
		}

        // Push into max_heap:
        if (val == val && val < max_heap[0]) { // Avoid NaN
            heapify(val, ix, max_heap, smallest_trees, k);
        }

    }
    
    // Write top k stuff to global memory:
    for (shortint i = 0; i < k; ++i) {
        max_heap_device[baseix + i * chunck_per_cycle] = max_heap[i];
        smallest_trees_device[baseix + i * chunck_per_cycle] = smallest_trees[i];
    }
}

__global__ void kernel_to_compute_optimal_tree_only_global_mem_4k(int chunck_per_cycle, shortint num_devices, shortint device_index, longint total_num_trees, shortint num_tree_vertices, shortint  T, shortint root_node, realnumber *data, void *scrape_memory , shortint shift_size, realnumber *max_heap_device, longint *smallest_trees_device, int k){
    
    int baseix = blockIdx.x*blockDim.x + threadIdx.x;
    
    char * shifted_pointer = (char *) scrape_memory;
    shifted_pointer = shifted_pointer + shift_size*baseix;
    
    // Allocate top k stuff:
    realnumber max_heap[topk];
    longint smallest_trees[topk];
    for (int i = 0; i < topk; ++i) {
        max_heap[i] = topk * MAX_N_NODES + i;
        smallest_trees[i] = 0;
    }
    
    // we only use global memory here
    shortint * adjacency_mat = (shortint *) shifted_pointer;
    shortint * adjacency_list = (shortint *) &(adjacency_mat[num_tree_vertices*num_tree_vertices]);
    edge *tree = (edge *) &(adjacency_list[num_tree_vertices*num_tree_vertices]);
    shortint * final_degrees = (shortint *) &(tree[num_tree_vertices-1]);
    shortint * perm_scrape = (shortint*) &(final_degrees[num_tree_vertices-2]);
    shortint * deg_scrape = (shortint*) &(perm_scrape[num_tree_vertices]);
    
    for (long int ix = device_index  +  baseix*num_devices  ; ix < total_num_trees ; ix = ix + num_devices*chunck_per_cycle){
            
        generate_tree_from_index( num_tree_vertices, ix, perm_scrape, deg_scrape, tree, adjacency_mat , final_degrees, adjacency_list );
        
		realnumber val  = tree_prior(num_tree_vertices, adjacency_mat, adjacency_list, final_degrees); //here we include a prior cost on the tree  
		// if the prior says that the tree will topology is not allowed, we do not even compute the rest 			
		if (val < 10000){
			val += tree_cost_inner(num_tree_vertices, T, data, root_node, tree, adjacency_mat, final_degrees, adjacency_list);	
		}

        // Push into max_heap:
        if (val == val && val < max_heap[0]) { // Avoid NaN
            heapify(val, ix, max_heap, smallest_trees, k);
        }

    }
    
    // Write top k stuff to global memory:
    for (shortint i = 0; i < k; ++i) {
        max_heap_device[baseix + i * chunck_per_cycle] = max_heap[i];
        smallest_trees_device[baseix + i * chunck_per_cycle] = smallest_trees[i];
    }
}


//*********************************   TESTS


void test_projection(){


	realnumber input[] = {0.01, 0.01266, 0.0035,0.02,-1.3,-1.3,0.165,0.4,0.1,0.01,0.1,0.005};
	int num_nodes = 12;

	realnumber output[num_nodes];
	realnumber output_condat[num_nodes];
	realnumber nsorted[num_nodes];
	realnumber cumsum[num_nodes];

	int i = 0;	
	int j = 0;
	
	realnumber error_proj_simplex_our = 0;
	realnumber error_proj_simplex_condat = 0;

	//call functions
	float time;
	cstart();
	for (int i = 1; i < 100000; i++)
		error_proj_simplex_our = projection_onto_inner_simplex( 1, NULL, input, num_nodes, nsorted, cumsum );
	cend(&time);

	printf( "error_proj_simplex our_version = %f, totaltime %f \n", error_proj_simplex_our, time );
	for ( i = 0; i < num_nodes; i++ )
	{
		
		printf( "%f ", output[i] );
	}
 	printf( "\n" );

	cstart();
	for (int i = 1; i < 100000; i++)
		error_proj_simplex_condat = projection_onto_inner_simplex_condat( 1, nsorted, input, num_nodes );
	cend(&time);


	printf( "error_proj_simplex condat_version = %f, totaltime %f\n", error_proj_simplex_condat, time );
	for ( i = 0; i < num_nodes; i++ )
	{
		
		printf( "%f ", output_condat[i] );
	}
 	printf( "\n" );

	// compare
	int err = 0;
	for ( j = 0; j < num_nodes; j++ )
	{
		
		if ( output[j] != output_condat[j] )
		{
			err = 1;
		}
		
	}
	printf( "Error flag = %d\n", err );

}


//******************************************




int main(int argc, const char * argv[]) {


    // argv[1] = cpu, cpu_multithread, gpu
    // argv[2] = cost1 (bei), cost2 (bei_inner), cost3 (ray), cost4 (ray_inner)
    // argv[3] = num_nodes
    // argv[4] = T
    // argv[5] = input_file
    // argv[6] = output_best_k
    // argv[7] = k
    // argv[8] = user_GPU_choice
    // argv[9] = number of CPU cores to use (might not be respected, the number input will be fixed to be a divisor of the number of trees)
    // argv[10] = number of devices that we will be using
    // argv[11] = which particular subset of trees the current device will work on
    // argv[12] = number of thread per block in CUDA
    // argv[13] = number of blocks
    

    printf( "Number of arguments = %d\n", argc );
    // check the number of arguments
    //if (argc != 12 || argc != 14) {
    //   printf("Invalid number of inputs!\n");
    //    return 0;
    //}
	if ( argc != 12 ) {
		if ( argc != 14 ) {
			printf("Invalid number of inputs!\n");
			return 0;
		}
	}
    
    // if number of arguments is OK, continue
    srand (time(NULL));
    
    // Get user input of device choice:
    char device_choice = 'C';
    if (strcmp(argv[1], "cpu") == 0) {
        device_choice = 'C';
    }
    else if (strcmp(argv[1], "cpu_multithread") == 0) {
        device_choice = 'M';
    }
    else if (strcmp(argv[1], "gpu") == 0) {
        device_choice = 'G';
    }
    else {
        printf("Invalid choice of device!\n");
        return 0;
    }
    
    // Get user input of cost function choice:
    char cost_choice = '1';
    if (strcmp(argv[2], "cost1") == 0) {
        cost_choice = '1';
    }
    else if (strcmp(argv[2], "cost2") == 0) {
        cost_choice = '2';
    }
    else if (strcmp(argv[2], "cost3") == 0) {
        cost_choice = '3';
    }
    else if (strcmp(argv[2], "cost4") == 0) {
        cost_choice = '4';
    }
    else {
        printf("Invalid choice of cost!\n");
        return 0;
    }
    
    FILE *input_file = fopen(argv[5], "r");
    FILE *output_best_k = fopen(argv[6], "w+");
    //FILE *run_time = fopen("run_time.txt", "w+");
    
	// read the input only until the correct size
	shortint num_tree_vertices = atoi(argv[3]);
    shortint T = atoi(argv[4]);

    // Read input data. Only read what we will use
    realnumber *input_data = (realnumber*) malloc(T*num_tree_vertices*sizeof(realnumber));
    for (int i = 0; i < T*num_tree_vertices; i++){
    	fscanf(input_file, "%f", &input_data[i]);
    	realnumber x = ((realnumber)rand() / (realnumber)(RAND_MAX/1e-4)) - 5e-5; // Random float between -5e-5 and 5e-5
    	input_data[i] += x; // Add small noise to data
    }

    
    float cpu_time;
    float gpu_time;
    

    shortint root_node = 0;  // in the end of the day we are looking for rooted trees so we need to know which node is the root. Normally, a virtual node that represents the background and wildtype mutant
    longint total_num_trees = pow_int(num_tree_vertices, num_tree_vertices - 2);
    printf("Total number of trees: %lu\n", total_num_trees);
    
    //realnumber smallest_cost = 100000;
    
    // Create an array to act as max heap
    const int k = atoi(argv[7]);
    if (k > total_num_trees) {
        printf("k is too large!\n");
        return 0;
    }
    
    realnumber max_heap[k];
    // Create an array to store trees with smallest costs:
    longint smallest_trees[k];
    for (shortint i = 0; i < k; ++i) {
        max_heap[i] = 100 * MAX_N_NODES - i;
        smallest_trees[i] = 0;
    }
    
    shortint num_threads_for_openmp = 1;    //it is important that this is a divisor of the total number of trees which is n^(n-2)
    shortint num_cpu_core = atoi(argv[9]);

    // the number of threads can be large or smaller than the number of trees and does not have to be a multiple of the total number of trees
    // we are using an interleaved distribution of tree indices per device and per thread.
    num_threads_for_openmp = num_cpu_core;
    
    
    printf("Number of OpenMP threads: %d\n", num_threads_for_openmp);
    omp_set_num_threads(num_threads_for_openmp);
    
    // we are going to allow the current program being lauched to just process a subset of the total number of trees so that we can exploit multiple devices,  GPUs, CPUs, etc.
    shortint num_devices = atoi(argv[10]);
    shortint device_index = atoi(argv[11]) - 1;


    longint tree_index = 0;
    
    // Create space for the Prufer code and also for the scrape degree array
    shortint code_length = num_tree_vertices - 2;
    shortint *permutation_array = (shortint *) malloc( code_length * sizeof(shortint) );
    shortint *degrees = (shortint *) malloc( num_tree_vertices * sizeof(shortint) );
    
    // Create space for the list of edges in the treee
    edge *tree = (edge*) malloc((num_tree_vertices - 1) * sizeof(edge));
    
    // Create space for the adjacency matrix
    shortint* adjacency_mat = (shortint *) malloc(num_tree_vertices*num_tree_vertices*sizeof(shortint));
    
    // Create space for the adjacency list
    shortint *final_degrees = (shortint *) malloc( num_tree_vertices * sizeof(shortint) );
    shortint* adjacency_list = (shortint *) malloc(num_tree_vertices*num_tree_vertices*sizeof(shortint)); // we are over initializing
    
    // Allocate arrays for top k trees for all OpenMP thread:
    realnumber *max_heap_mp = (realnumber *) malloc(k * num_threads_for_openmp * sizeof(realnumber));
    longint *smallest_trees_mp = (longint *) malloc(k * num_threads_for_openmp * sizeof(longint));
    for (shortint i = 0; i < k * num_threads_for_openmp; ++i) {
        max_heap_mp[i] = 100 * k * num_threads_for_openmp - i;
        smallest_trees_mp[i] = 0;
    }
    
    if (device_choice == 'C' || device_choice == 'M') {
        if (device_choice == 'C') {
            // ************************* START OF SINGLE CORE CODE **********************************

            cstart();
            for (tree_index = device_index; tree_index < total_num_trees ; tree_index = tree_index + num_devices){
                // Master function contains methods to generate permutations from index and generate prufer trees
                generate_tree_from_index(num_tree_vertices, tree_index, permutation_array, degrees, tree, adjacency_mat, final_degrees, adjacency_list);
                
            	//printf("%d, ",(int) tree_index);


                realnumber cost = 0;
                cost  = tree_prior(num_tree_vertices, adjacency_mat, adjacency_list, final_degrees); //here we include a prior cost on the tree  
                
                if (cost < 10000){
                    switch (cost_choice) {
                        case '1':
                            cost = tree_cost_bei_array(num_tree_vertices, T, input_data, root_node, tree, adjacency_mat, final_degrees, adjacency_list);
                            break;
                        case '2':
                            cost = tree_cost_bei_array_inner(num_tree_vertices, T, input_data, root_node, tree, adjacency_mat, final_degrees, adjacency_list);
                            break;
                        case '3':
                            cost = tree_cost(num_tree_vertices, T, input_data, root_node, tree, adjacency_mat, final_degrees, adjacency_list);
                            break;
                        case '4':
                            cost = tree_cost_inner(num_tree_vertices, T, input_data, root_node, tree, adjacency_mat, final_degrees, adjacency_list);
                            break;
                        default:
                            printf("Invalid cost choice!\n");
                            break;
                    }
                }
                
                if (cost == cost && cost < max_heap[0]) { // Avoid NaN
                    heapify(cost, tree_index, max_heap, smallest_trees, k);
                }
            }
            cend(&cpu_time);
        }
        else if (device_choice == 'M') {

            // ************************* START OF MULTI CORE CODE **********************************

            realnumber start = omp_get_wtime();
            #pragma omp parallel for //schedule(static)
            for (int part_num = 0; part_num < num_threads_for_openmp; part_num++){
                // we define severel local objects that need to be created for each of the parallel threads
                shortint local_permutation_array[MAX_N_NODES];
                shortint local_degrees[MAX_N_NODES];
                edge local_tree[MAX_N_NODES];
                shortint local_adjacency_mat[MAX_N_NODES*MAX_N_NODES];
                shortint local_final_degrees[MAX_N_NODES];
                shortint local_adjacency_list[MAX_N_NODES*MAX_N_NODES];
                
                // Top k stuff for each thread:
                realnumber max_heap_curr[k];
                longint smallest_trees_curr[k];
                for (shortint i = 0; i < k; ++i) {
                    max_heap_curr[i] = 10 * MAX_N_NODES * k + i;
                    smallest_trees_curr[i] = 0;
                }
                
                // this inner loop is not parallelized
                for (longint local_tree_index = device_index + part_num*num_devices; local_tree_index < total_num_trees; local_tree_index = local_tree_index + num_threads_for_openmp*num_devices){
                    

                    generate_tree_from_index( num_tree_vertices, local_tree_index, local_permutation_array, local_degrees,local_tree, local_adjacency_mat , local_final_degrees, local_adjacency_list );
                    
                    realnumber cost = 0;
                    cost  = tree_prior(num_tree_vertices, local_adjacency_mat, local_adjacency_list, local_final_degrees); //here we include a prior cost on the tree  

                    if (cost < 10000){
                        switch (cost_choice) {
                            case '1':
                                cost += tree_cost_bei_array(num_tree_vertices, T, input_data, root_node, tree, local_adjacency_mat, local_final_degrees, local_adjacency_list);
                                break;
                            case '2':
                                cost += tree_cost_bei_array_inner(num_tree_vertices, T, input_data, root_node, tree, local_adjacency_mat, local_final_degrees, local_adjacency_list);
                                break;
                            case '3':
                                cost += tree_cost(num_tree_vertices, T, input_data, root_node, tree, local_adjacency_mat, local_final_degrees, local_adjacency_list);
                                break;
                            case '4':
                                cost += tree_cost_inner(num_tree_vertices, T, input_data, root_node, tree, local_adjacency_mat, local_final_degrees, local_adjacency_list);
                                break;
                            default:
                                printf("Invalid cost choice!\n");
                                break;
                        }
                    }

                    
                    if (cost == cost && cost < max_heap_curr[0]) { // Avoid NaN
                        heapify(cost, local_tree_index, max_heap_curr, smallest_trees_curr, k);
                    }
                }
                
                // Write top k stuff to global memory:
                longint thread_id = omp_get_thread_num();
                for (shortint i = 0; i < k; ++i) {
                    max_heap_mp[thread_id + i * num_threads_for_openmp] = max_heap_curr[i]; // part_num = thread index
                    smallest_trees_mp[thread_id + i * num_threads_for_openmp] = smallest_trees_curr[i];
                    //max_heap_mp[thread_id * k + i] = max_heap_curr[i]; // part_num = thread index
                    //smallest_trees_mp[thread_id * k + i] = smallest_trees_curr[i];
                }
            }
            realnumber end = omp_get_wtime();
            cpu_time = (float) 1000*(end - start);
        }
        
        printf("%f = CPU time\n", cpu_time);
        //fprintf(run_time, "%f\n", cpu_time);
        
        shortint stack[MAX_N_NODES];
        shortint visited[MAX_N_NODES];
        shortint father_list[2*num_tree_vertices + 1];
        for (shortint i = 0; i < num_tree_vertices; i++){
            visited[i] = 0;
        }
        
        if (device_choice == 'C'){
            printf("%d trees with smallest costs (tree index, cost): \n", k);
            for (shortint i = 0; i < k; ++i) {
                printf("(%lu, %f)\n", smallest_trees[i], max_heap[i]);
                generate_tree_from_index(num_tree_vertices, smallest_trees[i], permutation_array, degrees, tree, adjacency_mat, final_degrees, adjacency_list);
                dfs_tree_compute_fathers_non_recursive_array(num_tree_vertices, father_list, adjacency_list, final_degrees, root_node, stack, visited);
                fprintf(output_best_k, "%f\n", max_heap[i]);
                for (shortint j = 0 ; j < num_tree_vertices ; ++j){
                    if (j != root_node) {
                        printf("(%d, %d)\n", j + 1, father_list[j + 1]);
                        fprintf(output_best_k, "(%d, %d)\n", j + 1, father_list[j + 1]);
                    }
                }
                printf("\n");
                fprintf(output_best_k, "\n");
            }
            printf("\n");
        }
        else {
            // Merge all local top k from all threads together:
            for (int i = 0; i < k * num_threads_for_openmp; ++i) {
                heapify(max_heap_mp[i], smallest_trees_mp[i], max_heap, smallest_trees, k);
            }
            printf("%d trees with smallest costs (tree index, cost): \n", k);
            for (shortint i = 0; i < k; ++i) {
                printf("(%lu, %f)\n", smallest_trees[i], max_heap[i]);
                generate_tree_from_index(num_tree_vertices, smallest_trees[i], permutation_array, degrees, tree, adjacency_mat, final_degrees, adjacency_list);
                dfs_tree_compute_fathers_non_recursive_array(num_tree_vertices, father_list, adjacency_list, final_degrees, root_node, stack, visited);
                fprintf(output_best_k, "%f\n", max_heap[i]);
                for (shortint j = 0 ; j < num_tree_vertices ; ++j){
                    if (j != root_node) {
                        printf("(%d, %d)\n", j + 1, father_list[j + 1]);
                        fprintf(output_best_k, "(%d, %d)\n", j + 1, father_list[j + 1]);
                    }
                }
                printf("\n");
                fprintf(output_best_k, "\n");
            }
            printf("\n");
        }
    }
    //
    //
    // ***************************** START OF GPU CODE *****************************************************************************
    //
    //
    //
    else if (device_choice == 'G') {
        cudaSetDevice(atoi(argv[8]));
        cudaDeviceReset();
        
        int numthreadsperblock = 32;
        int numblocks =  128;

        if (argc == 14){
            numthreadsperblock = atoi(argv[12]);
            numblocks = atoi(argv[13]);
        }

        printf("Cuda is using %d threads per block and %d blocks\n",numthreadsperblock,numblocks);


        int chunck_per_cycle = numblocks*numthreadsperblock;
        //int num_trees_per_thread = ((total_num_trees + chunck_per_cycle - 1)/ chunck_per_cycle );
        
        // Allocate arrays for local top k trees for all CUDA threads:
        realnumber *max_heap_host = (realnumber *) malloc(topk * chunck_per_cycle * sizeof(realnumber));
        longint *smallest_trees_host = (longint *) malloc(topk * chunck_per_cycle * sizeof(longint));
        for (int i = 0; i < topk * chunck_per_cycle; ++i) {
            max_heap_host[i] = 100 * topk * chunck_per_cycle - i;
            smallest_trees_host[i] = 0;
        }
        
        // create space in the GPU global memory to store the best/top k values found
        realnumber *max_heap_device;
        longint *smallest_trees_device;
        cudaMalloc((void **)&max_heap_device, topk * chunck_per_cycle * sizeof(realnumber));
        cudaMalloc((void **)&smallest_trees_device, topk * chunck_per_cycle * sizeof(longint));
        
        // Initialize GPU global max_heap stuff to have heap property:
        cudaMemcpy((void*) max_heap_device, (void*) max_heap_host, topk * chunck_per_cycle * sizeof(realnumber), cudaMemcpyHostToDevice);
        cudaMemcpy((void*) smallest_trees_device, (void*) smallest_trees_host, topk * chunck_per_cycle * sizeof(longint), cudaMemcpyHostToDevice);
        
        // create space in the GPU to store the mutation frequency data
        realnumber * data_device;
        cudaMalloc((void **)&data_device, T*num_tree_vertices*sizeof(realnumber) );
        cudaMemcpy( (void*) data_device , (void*) input_data , T*num_tree_vertices*sizeof(realnumber) , cudaMemcpyHostToDevice );
        
        // now we need to create space in the GPU to accomodate all the scrap memory we need
        
        int amount_of_mem_per_thread = sizeof(edge)*(num_tree_vertices-1) + sizeof(shortint)*(num_tree_vertices*num_tree_vertices + num_tree_vertices + (num_tree_vertices)*(num_tree_vertices) + num_tree_vertices-2 + num_tree_vertices); // + sizeof_treecost_scrapmem
        
        void * device_scrape_memory;
        cudaMalloc((void **)&device_scrape_memory,amount_of_mem_per_thread*numblocks*numthreadsperblock ); // here we allocate all the scrape memory that we need
        
        gstart();
        
        switch (cost_choice) {
            case '1':
                kernel_to_compute_optimal_tree_only_global_mem_1k<<<numblocks,numthreadsperblock>>>(chunck_per_cycle, num_devices, device_index, total_num_trees,  num_tree_vertices, T, root_node, data_device, device_scrape_memory, amount_of_mem_per_thread, max_heap_device, smallest_trees_device, k);
                break;
            case '2':
                kernel_to_compute_optimal_tree_only_global_mem_2k<<<numblocks,numthreadsperblock>>>(chunck_per_cycle, num_devices, device_index, total_num_trees,  num_tree_vertices, T, root_node, data_device, device_scrape_memory, amount_of_mem_per_thread, max_heap_device, smallest_trees_device, k);
                break;
            case '3':
                kernel_to_compute_optimal_tree_only_global_mem_3k<<<numblocks,numthreadsperblock>>>(chunck_per_cycle, num_devices, device_index, total_num_trees,  num_tree_vertices, T, root_node, data_device, device_scrape_memory, amount_of_mem_per_thread, max_heap_device, smallest_trees_device, k);
                break;
            case '4':
                kernel_to_compute_optimal_tree_only_global_mem_4k<<<numblocks,numthreadsperblock>>>(chunck_per_cycle, num_devices, device_index, total_num_trees,  num_tree_vertices, T, root_node, data_device, device_scrape_memory, amount_of_mem_per_thread, max_heap_device, smallest_trees_device, k);
                break;
            default:
                printf("Invalid cost choice!\n");
                break;
        }
        
        gend(&gpu_time);
        printf("%f = GPU time\n",gpu_time);
        //fprintf(run_time, "%f\n", gpu_time);
        
        // Copy all local top k stuff back to CPU memory:
        cudaMemcpy((void*) max_heap_host, (void*) max_heap_device, topk * chunck_per_cycle * sizeof(realnumber), cudaMemcpyDeviceToHost);
        cudaMemcpy((void*) smallest_trees_host, (void*) smallest_trees_device, topk * chunck_per_cycle * sizeof(longint), cudaMemcpyDeviceToHost);
        
        shortint stack[MAX_N_NODES];
        shortint visited[MAX_N_NODES];
        shortint father_list[2*num_tree_vertices + 1];
        for (shortint i = 0; i < num_tree_vertices; i++){
            visited[i] = 0;
        }
        
        if (k <= 100 && k <= total_num_trees) {
            // Merge all local top k from all threads together:
            for (int i = 0; i < topk * chunck_per_cycle; ++i) {
                heapify(max_heap_host[i], smallest_trees_host[i], max_heap, smallest_trees, k);
            }
            printf("%d trees with smallest costs (tree index, cost): \n", k);
            for (shortint i = 0; i < k; ++i) {
                printf("(%lu, %f)\n", smallest_trees[i], max_heap[i]);
                generate_tree_from_index(num_tree_vertices, smallest_trees[i], permutation_array, degrees,tree, adjacency_mat, final_degrees, adjacency_list);
                dfs_tree_compute_fathers_non_recursive_array(num_tree_vertices, father_list, adjacency_list, final_degrees, root_node, stack, visited);
                fprintf(output_best_k, "%f\n", max_heap[i]);
                for (shortint j = 0 ; j < num_tree_vertices ; ++j){
                    if (j != root_node) {
                        printf("(%d, %d)\n", j + 1, father_list[j + 1]);
                        fprintf(output_best_k, "(%d, %d)\n", j + 1, father_list[j + 1]);
                    }
                }
                printf("\n");
                fprintf(output_best_k, "\n");
            }
            printf("\n");
        }
        else {
            printf("For GPU, k <= min(100, total number of trees)!\n");
            fprintf(output_best_k, "For GPU, k <= min(100, total number of trees)!\n");
        }
        
        gerror( cudaPeekAtLastError() );
        cudaDeviceSynchronize();
        
        cudaFree(data_device);
        cudaFree(max_heap_device);
        cudaFree(smallest_trees_device);
        free(max_heap_host);
        free(smallest_trees_host);
    }
    //
    //
    // ***************************** END OF GPU CODE *****************************************************************************
    //
    //
    else {
        printf("Invalid choice of device!\n");
    }
    
    //free all the dynamically allocated memory in the CPU
    free(tree);
    free(degrees);
    free(permutation_array);
    //free(treecost_scrapmem);
    free(adjacency_list);
    free(adjacency_mat);
    free(final_degrees);
    free(max_heap_mp);
    free(smallest_trees_mp);
    free(input_data);

    // Close files:
    fclose(input_file);
    fclose(output_best_k);
    //fclose(run_time);
    
    
    return 0;
}
