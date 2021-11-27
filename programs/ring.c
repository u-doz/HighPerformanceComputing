#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char** argv) {

  int rank, size;
  
  MPI_Init(&argc, &argv);
  //what is the ID of this processor ?
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  //how many processors are associated with a communicator?
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  // Receive from the lower process and send to the higher process. To prevent deadlock,
  // the first process to send the message is process 0 and after that, it will be
  // able to receive, and the other processes will be able to receive and then to send
  int msgleft, msgright, tagleft, tagright;
  int np = 0;
  MPI_Status status;
  if (world_rank != 0) {
    MPI_Recv(&msgright, 1, MPI_INT, (rank - 1), MPI_ANY_TAG, MPI_COMM_WORLD,
             &status);
    tagright = status.MPI_TAG;
    MPI_Recv(&msgleft, 1, MPI_INT, (rank + 1), MPI_ANY_TAG, MPI_COMM_WORLD,
             &status);
    tagleft = status.MPI_TAG;
    printf("Process %d received token %d from process %d with tag %d \n", rank, msgright,
           rank - 1, tagright);
    printf("Process %d received token %d from process %d with tag %d \n", rank, msgleft,
           rank + 1, tagleft);
    np += 2;  
  } 
  
  MPI_Send(rank, 1, MPI_INT, (rank - 1), rank*10, MPI_COMM_WORLD);
  MPI_Send(-rank, 1, MPI_INT, (rank + 1), rank*10, MPI_COMM_WORLD);
  // Now process 0 can receive from the last process. This makes sure that at
  // least one MPI_Send is initialized before all MPI_Recvs (again, to prevent
  // deadlock)
  if (world_rank == 0) {
    MPI_Recv(&msgright, 1, MPI_INT, (rank - 1), MPI_ANY_TAG, MPI_COMM_WORLD,
             &status);
    tagright = status.MPI_TAG;
    MPI_Recv(&msgleft, 1, MPI_INT, (rank + 1), MPI_ANY_TAG, MPI_COMM_WORLD,
             &status);
    tagleft = status.MPI_TAG;
    printf("Process %d received token %d from process %d with tag %d \n", rank, msgright,
           rank - 1, tagright);
    printf("Process %d received token %d from process %d with tag %d \n", rank, msgleft,
           rank + 1, tagleft);
    np += 2;  
  }
  
  printf("I am process %d and I have received %d messages. My final messages have tag %d, %d and value %d, %d \n", rank, np, tagleft, tagright, msgleft, msgright);
  MPI_Finalize();
}
