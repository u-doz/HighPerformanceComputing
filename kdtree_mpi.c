#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <mpi.h>
#include <math.h>
#include <omp.h>
#include<time.h>
#define COUNT 10
#define MAX_DIM 2


struct kdnode{
    double x[MAX_DIM];
    int d;
    struct kdnode *left, *right;
};

void print_kdnode(struct kdnode *n, int size, int dim) {
    printf("Array of size %d : \n", size);
    int i,j;
    for(i = 0; i< size; i++){
         for(j = 0; j<dim; j++){
                printf(" : %f", n->x[j]);
         }
         printf(" \n");
         n = n+1;
    }
    printf(" \n");

}
    
void print2DUtil(struct kdnode *root, int space)
{
    // Base case
    if (root == NULL)
        return;
    int dim = root->d;
    int i;
    // Increase distance between levels
    space += COUNT;

    // Process right child first
    print2DUtil(root->right, space);

    // Print current node after space
    printf("\n");
    for (i = COUNT; i < space; i++)
        printf(" ");

    printf("%f (dim: %d ) \n", root->x[dim], dim);

    // Process left child
    print2DUtil(root->left, space);
}

/* prints the kd-tree from the root to the leaves, level by level.
 * The output can be viewed in horizontal from left to right */
void print2D(struct kdnode *root)
{
   print2DUtil(root, 0);
}
double drand ( double low, double high )
{
    return ( (double)rand() * ( high - low ) ) / (double)RAND_MAX + low;
}

/* populate the values of a list of kd-nodes (of lenght size) with random numbers between 0 and 1000 */
void getInput(struct kdnode wp[], int size, int dim){
   int i, j;
   for(i=0; i<size; i++) {
      for(j=0;j<dim;j++) {
         wp->x[j] = drand(0, 1000);
      }
      wp = wp+1;
   }
}

/* create a copy (by value) of the list of kd-nodes n on t_right */   
void getRightTree(struct kdnode t_right[], int size, int dim, struct kdnode *n){
   int i;
   for(i=0; i<size; i++) {
      double n0 = (n+i)->x[0];
      double n1 = (n+i)->x[1];

      t_right->x[0] = n0;
      t_right->x[1] = n1;
      t_right = t_right+1;
   }
}

void swap(struct kdnode *xp, struct kdnode *yp)
{
    struct kdnode temp = *xp;
    *xp = *yp;
    *yp = temp;
}

/* sort a list of kd-nodes with the bubble sort algorithm */
void bubbleSort(struct kdnode *t, int n, int dim)
{
    {
   int i, j;
   for (i = 0; i < n-1; i++) {
       for (j = 0; j < n-i-1; j++) {
           struct kdnode *current, *next;
           current = t+j;
           next = t+j+1;
           if (current->x[dim] > next->x[dim]){
              swap(t+j, t+j+1);
           }
       }
   }
   }
}

/* find the median kd-node of a list of kd-nodes */
struct kdnode* find_median2(struct kdnode *start, struct kdnode *end)
{
    if (end <= start) return NULL;
    if (end == start + 1){
        return start;
    }

    struct kdnode  *md = start + (end - start) / 2;
    return md;
}

    
    
    /* merge two sorted arrays v1, v2 of lengths n1, n2, respectively */
void  merge(struct kdnode * v1, int n1, struct kdnode * v2, int n2, int dim, struct kdnode * t, int id)
{
  int i = 0;
  int j = 0;
  int k, h;

  for (k = 0; k < n1 + n2; k++) {
    struct kdnode *temp;
    if (i >= n1) {
      temp = v2+j;
      j++;
    }
    else if (j >= n2) {
      temp = v1+i;
      i++;
    }
    else if ((v1+i)->x[dim] < (v2+j)->x[dim]) {
      temp = v1+i;
      i++;
    }
    else { // v2[j] <= v1[i]
      temp = v2+j;
      j++;
    }
    for(h=0;h<2;h++) {
         t->x[h] = temp->x[h];
      }
   t++;
   }
}

/* create the kd-tree structure by recursion */
struct kdnode * make_tree(struct kdnode *t, int len, int i, int dim, int p, int id)
{

    i = (i + 1) % dim;
    int step, c, s, o;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Status status;
    MPI_Bcast(&len, 1, MPI_INT, 0, MPI_COMM_WORLD);
    printf("I am processor %d and I START the tree  \n", id);
    if(len > 1){
       // compute the chunk size for all the processors (approximated upwards)
       c = (len %p == 0) ?(len/p) :(len /(p - 1));

       // compute the chunk size specific for this processor
       s = (len >= c * (id+1)) ? c : len - c * id;

       //scatter the chunks to the processors
       struct kdnode *chunk = (struct kdnode*)malloc(c*sizeof(struct kdnode));
       MPI_Scatter(t, c*sizeof(struct kdnode), MPI_BYTE, chunk, c*sizeof(struct kdnode), MPI_BYTE, 0, MPI_COMM_WORLD);

       // sort the chunks
       bubbleSort(chunk, s, i);
       MPI_Barrier(MPI_COMM_WORLD);


       // up to log_2 p merge steps
       for (step = 1; step < p; step = 2*step) {
           if (id % (2*step)!=0) {
              // id is no multiple of 2*step: send chunk to id-step and exit loop
              if(s>0){
                   MPI_Send(chunk, s*sizeof(struct kdnode), MPI_BYTE, (id-step), 0, MPI_COMM_WORLD);
              }
              break;
            }
            // id is multiple of 2*step: merge in chunk from id+step (if it exists)
            if (id+step < p) {
               // compute size of chunk to be received
               o = (len >= c * (id+2*step)) ? c * step : len - c * (id+step);
               // receive other chunk
               if(o>0){
                   struct kdnode* other = (struct kdnode*)malloc(o*sizeof(struct kdnode));
                   MPI_Recv(other, o*sizeof(struct kdnode), MPI_BYTE, (id+step), 0, MPI_COMM_WORLD, &status);

                   // merge the chunk of the processor with the chunk received
                   struct kdnode *t = (struct kdnode*)malloc((s+o)*sizeof(struct kdnode));
                   merge(chunk, s, other, o, i, t, id);

                   free(other);
                   chunk = t;
                   s = s + o;
                   printf("I am processor %d and I MERGED. Now my array has size %d \n", id, s);
                }
             }
          }
    }
    
    //now the list of kd-nodes is sorted and it's possible to find the median and start the recursion 
    struct kdnode *n;
    if(id == 0){
        n = find_median2(t, t + len);
    }
    else{
       t = NULL;
       n = NULL;
    }
    if (!len){
       printf("I am processor %d and I RETURN NULL \n", id);
       n= NULL;
    return 0;
    }
    else if(len == 1){
       printf("I am processor %d and I RETURN a leaf \n", id);
       return n;
    }
    else{
       int size_left, size_right;
       size_left = (int)(n-t);
       size_right = (int)(t + len - (n+1));

       struct kdnode t_right[size_right];
       if(id == 0){
          getRightTree(t_right, size_right,2,n+1);
       }

    printf("I am processor %d and I should start the left recursion \n", id);
    if(size_left > 0)
        n->left = make_tree(t, size_left, i, dim, p, id);
    printf("I am processor %d and I should start the right recursion \n", id);
    if(id == 0){
        printf("KDNODE FOR RIGHT RECURSION \n");
        print_kdnode(t_right, size_right, 2);
    }
    if(size_right > 0)
        n->right = make_tree(t_right, size_right, i, dim, p, id);

    //associate the dimension to the kd-node  
    n->d = i;
    }

    return n;
}


int main (int argc, char *argv[]) {

    MPI_Init(&argc, &argv);
    int size, id, p;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    if (id == 0) {
        printf("Insert the number of nodes \n");
        scanf("%d", &size);
   }

    int dim = 2;
    struct kdnode wp[size];
    //only processor 0 will have the input and so the kd-tree will be computed only once
    if (id == 0) {
       getInput(wp, size, dim);
       print_kdnode(wp, size, 2);
       printf("I am building the tree with P = %d \n", p);
    }
    double tstart = MPI_Wtime();
    struct kdnode *tree;

    tree =  make_tree(wp, size, 1, dim, p, id);

    if (id == 0) {
       double tend = MPI_Wtime();
       printf("Building tree DONE! \n");
       if(size < 20)
          print2D(tree);
       printf("Time: %g of wall-clock time \n", tend-tstart);
       printf("N PROC = %d \n", p);
    }
    MPI_Finalize();
    return 0;
}

