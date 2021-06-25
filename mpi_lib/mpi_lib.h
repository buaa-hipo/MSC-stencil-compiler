
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

// global status and request for overlap of communication and computation
MPI_Status status_s, status_r;
MPI_Request handler_s, handler_r;
// buffer
#define BUFFER(Type, size)\
Type * buffer_send_##Type [size];\
Type * buffer_recv_##Type [size];

BUFFER(float, 8)
BUFFER(double, 8)
BUFFER(int, 8)


// 1D

#define PRINT_GRID_1D_BARRIER(Type, flag)   \
void print_grid_1D_##Type##_barrier(Type * input, int d1, int halo_size, const int pid){\
  /*barrier*/       \
  MPI_Barrier(MPI_COMM_WORLD);      \
  /* delay*/      \
  for(int j=0;j<pid*pid*pid*10000;j++){             \
    double ans=(pid*4524523523454+100)/452352545425*pid+j;    \
  }                                         \
  /* print_grid*/                                         \
  printf("<------------------------->\n");    \
  printf("from process : %d\n", pid);       \
  Type * grid = input;                    \
  for(int j=0;j<d1+2*halo_size;j++){                    \
    printf( #flag , grid[j]);                  \
    printf( " " );                                    \
  }                                                     \
  printf("\n");                                       \
}                                      

PRINT_GRID_1D_BARRIER(float, %.f)
PRINT_GRID_1D_BARRIER(double, %.f)
PRINT_GRID_1D_BARRIER(int, %d)

int within_bound_1D(const int ind1, const int D1){
  if( 0<=ind1 && ind1<D1 ) return 1;
  return 0;
}

#define EXCHANGE_HALO_1D(Type, MPI_TYPE)\
void exchange_halo_1D_##Type(Type * input, int dt, int d1, int halo_size, int pid, int D1, int mpi_size){\
  Type * grid = input;\
  /* start exchange halos*/     \
  int dest_pid;\
  int tag=0;\
  /*left*/   \
  if( within_bound_1D(pid-1, D1) ){  \
    /*isend && irecv*/  \
    dest_pid = pid-1; \
    int cnt = halo_size; \
    MPI_Isend (grid+halo_size, cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (grid          , cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_r);  \
    /*wait*/ \
    MPI_Wait(&handler_s, &status_s); \
    MPI_Wait(&handler_r, &status_r); \
  }\
  /*right*/ \
  if( within_bound_1D(pid+1, D1) ){ \
    /*isend && irecv*/ \
    dest_pid = pid+1; \
    int cnt = halo_size; \
    MPI_Isend (grid+halo_size+d1-halo_size, cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_s); \
    MPI_Irecv (grid+halo_size+d1          , cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_r); \
    /*wait*/ \
    MPI_Wait(&handler_s, &status_s); \
    MPI_Wait(&handler_r, &status_r); \
  } \
}

EXCHANGE_HALO_1D(float, MPI_FLOAT)
EXCHANGE_HALO_1D(double, MPI_DOUBLE)
EXCHANGE_HALO_1D(int, MPI_INT)

#define EXCHANGE_HALO_1D_START(Type, MPI_TYPE)\
void exchange_halo_1D_##Type##_start(Type * input, int dt, int d1, int halo_size, int pid, int D1, int mpi_size){\
  Type * grid = input;\
  print_grid_1D_##Type##_barrier(input, d1, halo_size, pid);\
  /* start exchange halos*/     \
  int dest_pid;\
  int tag=0;\
  /*left*/   \
  if( within_bound_1D(pid-1, D1) ){  \
    /*isend && irecv*/  \
    dest_pid = pid-1; \
    int cnt = halo_size; \
    MPI_Isend (grid+halo_size, cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (grid          , cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_r);  \
  }\
  /*right*/ \
  if( within_bound_1D(pid+1, D1) ){ \
    /*isend && irecv*/ \
    dest_pid = pid+1; \
    int cnt = halo_size; \
    MPI_Isend (grid+halo_size+d1-halo_size, cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_s); \
    MPI_Irecv (grid+halo_size+d1          , cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_r); \
  } \
}

EXCHANGE_HALO_1D_START(float, MPI_FLOAT)
EXCHANGE_HALO_1D_START(double, MPI_DOUBLE)
EXCHANGE_HALO_1D_START(int, MPI_INT)

#define EXCHANGE_HALO_1D_END(Type, MPI_TYPE)\
void exchange_halo_1D_##Type##_end(Type * input, int dt, int d1, int halo_size, int pid, int D1, int mpi_size){\
  Type * grid = input;\
  /* start exchange halos*/     \
  int dest_pid;\
  int tag=0;\
  /*left*/   \
  if( within_bound_1D(pid-1, D1) ){  \
    /*wait*/ \
    MPI_Wait(&handler_s, &status_s); \
    MPI_Wait(&handler_r, &status_r); \
  }\
  /*right*/ \
  if( within_bound_1D(pid+1, D1) ){ \
    /*wait*/ \
    MPI_Wait(&handler_s, &status_s); \
    MPI_Wait(&handler_r, &status_r); \
  } \
  print_grid_1D_##Type##_barrier(input, d1, halo_size, pid); \
}

EXCHANGE_HALO_1D_END(float, MPI_FLOAT)
EXCHANGE_HALO_1D_END(double, MPI_DOUBLE)
EXCHANGE_HALO_1D_END(int, MPI_INT)

// 2D

void pid2indice_2D(const int pid, const int D2, int *ind1, int *ind2){
  *ind1 = pid / D2;
  *ind2 = pid % D2;
}

void indice2pid_2D(int *pid, const int D2, const int ind1, const int ind2){
  *pid = ind1 * D2 + ind2;
}

void get_neighbors_2D(const int ind1, const int ind2, int D1, int D2, int * nei_ind1, int * nei_ind2, int * num_neighbors){
  *num_neighbors = 0;
  for(int i=ind1-1; i<= ind1+1; i++){
    for(int j=ind2-1; j<= ind2+1; j++){
      if(i==ind1 && j== ind2)continue;
      if( 0<=i && i<D1 && 0<=j && j<D2){
        nei_ind1[ *num_neighbors ] = i;
        nei_ind2[ *num_neighbors ] = j;
        (*num_neighbors)+=1;
      }
    }
  }
}

int within_bound_2D(const int ind1, const int ind2, const int D1, const int D2){
  if( 0<=ind1 && ind1<D1 && 0<=ind2 && ind2<D2 ) return 1;
  return 0;
}

int is_neighbor_2D(const int ind1, const int ind2, int * nei_ind1, int * nei_ind2, int num_neighbors){
  int flag=0;
  for(int i=0; i< num_neighbors; i++){
    if( nei_ind1[i]==ind1 && nei_ind2[i]==ind2)flag=1;
  }
  return flag;
}


#define PRINT_GRID_2D_BARRIER(Type, flag)   \
void print_grid_2D_##Type##_barrier(Type * input, int d1, int d2, int halo_size, const int pid){\
  /*barrier*/       \
  MPI_Barrier(MPI_COMM_WORLD);      \
  /* delay*/      \
  for(int j=0;j<pid*pid*pid*10000;j++){             \
    double ans=(pid*4524523523454+100)/452352545425*pid+j;    \
  }                                         \
  /* print_grid*/                                         \
  printf("<------------------------->\n");    \
  printf("from process : %d\n", pid);       \
  Type (*grid) [d2+2*halo_size] = (Type (*)[d2+2*halo_size]) input;\
  for(int j=0;j<d1+2*halo_size;j++){                    \
    for(int k=0;k<d2+2*halo_size;k++){                  \
      printf( #flag , grid[j][k]);                  \
      printf( " " );                                    \
    }                                                   \
    printf("\n");                                       \
  }                                                     \
}                                      

PRINT_GRID_2D_BARRIER(float, %.f)
PRINT_GRID_2D_BARRIER(double, %.f)
PRINT_GRID_2D_BARRIER(int, %d)

#define PRINT_GRID_2D(Type, flag)\
void print_grid_2D_##Type(Type * input, int d1, int d2, int halo_size, const int pid, int * nei_ind1, int * nei_ind2, int num_neighbors){\
  Type (*grid) [d2+2*halo_size] = (Type (*)[d2+2*halo_size]) input;\
  printf("<------------------------->\n");\
  printf("from process : %d\n", pid);\
  for(int i=0;i<num_neighbors;i++){\
    printf("%d, %d\n",nei_ind1[i], nei_ind2[i]);\
  }\
  for(int j=0;j<d1+2*halo_size;j++){\
    for(int k=0;k<d2+2*halo_size;k++){\
      printf( #flag, grid[j][k]);\
      printf( " " );\
    }\
    printf("\n");\
  }\
}

PRINT_GRID_2D(float, %.f)
PRINT_GRID_2D(double, %.f)
PRINT_GRID_2D(int, %d)

#define EXCHANGE_HALO_2D(Type, MPI_TYPE)\
void exchange_halo_2D_##Type(Type * input, int dt, int d1, int d2, int halo_size, int pid, int D1, int D2, int mpi_size){\
  Type (*grid) [d2+2*halo_size] = (Type (*)[d2+2*halo_size]) input;   \
  int nei_ind1[9];   \
  int nei_ind2[9];   \
  int num_neighbors;   \
  int ind1;   \
  int ind2;   \
  pid2indice_2D(pid, D2, &ind1, &ind2);   \
  get_neighbors_2D(ind1, ind2, D1, D2, nei_ind1, nei_ind2, &num_neighbors);   \
  /* start exchange halos*/   \
  const int sym_tags = 8;   \
  int dest_pid;   \
  int tag=0;   \
  MPI_Status status_s, status_r;   \
  MPI_Request handler_s, handler_r;   \
  /* upper : isend tag 0 , irecv 8 = sym_tags-tag*/   \
  if( within_bound_2D(ind1-1, ind2, D1, D2) ){   \
    /*buffer*/   \
    Type * buffer_send = (Type *) malloc ( d2*halo_size*sizeof(Type) );   \
    Type * buffer_recv = (Type *) malloc ( d2*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;   \
    for(int i=halo_size; i<2*halo_size; i++){   \
      for(int j=halo_size; j<halo_size+d2; j++){   \
        buffer_send[cnt++] = grid[i][j];   \
      }   \
    }   \
    /*isend && irecv*/   \
    tag=0;   \
    indice2pid_2D( &dest_pid, D2, ind1-1, ind2);   \
    MPI_Isend (buffer_send, cnt, MPI_TYPE, dest_pid, tag,          MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv, cnt, MPI_TYPE, dest_pid, sym_tags-tag, MPI_COMM_WORLD, &handler_r);   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=0; i<halo_size; i++){   \
      for(int j=halo_size; j<halo_size+d2; j++){   \
        grid[i][j] = buffer_recv[cnt++];   \
      }   \
    }   \
    free(buffer_send);   \
    free(buffer_recv);   \
  }   \
  /* bottom : isend tag 8 , irecv 0 = sym_tags-tag*/   \
  if( within_bound_2D(ind1+1, ind2, D1, D2) ){   \
    /*buffer*/   \
    Type * buffer_send = (Type *) malloc ( d2*halo_size*sizeof(Type) );   \
    Type * buffer_recv = (Type *) malloc ( d2*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;   \
    for(int i=halo_size+d1-halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size; j<halo_size+d2; j++){   \
        buffer_send[cnt++] = grid[i][j];   \
      }   \
    }   \
    /*isend && irecv*/   \
    tag=8;   \
    indice2pid_2D( &dest_pid, D2, ind1+1, ind2);   \
    MPI_Isend (buffer_send, cnt, MPI_TYPE, dest_pid, tag,          MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv, cnt, MPI_TYPE, dest_pid, sym_tags-tag, MPI_COMM_WORLD, &handler_r);   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=halo_size+d1; i<2*halo_size+d1; i++){   \
      for(int j=halo_size; j<halo_size+d2; j++){   \
        grid[i][j] = buffer_recv[cnt++];   \
      }   \
    }   \
    /**/   \
    free(buffer_send);   \
    free(buffer_recv);   \
  }   \
  /* left : isend tag 1 , irecv 7 = sym_tags-tag*/   \
  if( within_bound_2D(ind1, ind2-1, D1, D2) ){   \
    /*buffer*/   \
    Type * buffer_send = (Type *) malloc ( d1*halo_size*sizeof(Type) );   \
    Type * buffer_recv = (Type *) malloc ( d1*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;   \
    for(int i=halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size; j<2*halo_size; j++){   \
        buffer_send[cnt++] = grid[i][j];   \
      }   \
    }   \
    /*isend && irecv*/   \
    tag=1;   \
    indice2pid_2D( &dest_pid, D2, ind1, ind2-1);   \
    MPI_Isend (buffer_send, cnt, MPI_TYPE, dest_pid, tag,          MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv, cnt, MPI_TYPE, dest_pid, sym_tags-tag, MPI_COMM_WORLD, &handler_r);   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=halo_size; i<halo_size+d1; i++){   \
      for(int j=0; j<halo_size; j++){   \
        grid[i][j] = buffer_recv[cnt++];   \
      }   \
    }   \
    free(buffer_send);   \
    free(buffer_recv);   \
  }   \
  /* right : isend tag 7 , irecv 1 = sym_tags-tag*/   \
  if( within_bound_2D(ind1, ind2+1, D1, D2) ){   \
    /*buffer*/   \
    Type * buffer_send = (Type *) malloc ( d1*halo_size*sizeof(Type) );   \
    Type * buffer_recv = (Type *) malloc ( d1*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;   \
    for(int i=halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size+d2-halo_size; j<halo_size+d2; j++){   \
        buffer_send[cnt++] = grid[i][j];   \
      }   \
    }   \
    /*isend && irecv*/   \
    tag=7;   \
    indice2pid_2D( &dest_pid, D2, ind1, ind2+1);   \
    MPI_Isend (buffer_send, cnt, MPI_TYPE, dest_pid, tag,          MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv, cnt, MPI_TYPE, dest_pid, sym_tags-tag, MPI_COMM_WORLD, &handler_r);   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size+d2; j<2*halo_size+d2; j++){   \
        grid[i][j] = buffer_recv[cnt++];   \
      }   \
    }   \
    /**/   \
    free(buffer_send);   \
    free(buffer_recv);   \
  }   \
  /* upper left : isend tag 2 , irecv 6 = sym_tags-tag*/   \
  if( within_bound_2D(ind1-1, ind2-1, D1, D2) ){   \
    /*buffer*/   \
    Type * buffer_send = (Type *) malloc ( halo_size*halo_size*sizeof(Type) );   \
    Type * buffer_recv = (Type *) malloc ( halo_size*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;   \
    for(int i=halo_size; i<2*halo_size; i++){   \
      for(int j=halo_size; j<2*halo_size; j++){   \
        buffer_send[cnt++] = grid[i][j];   \
      }   \
    }   \
    /*isend && irecv*/   \
    tag=2;   \
    indice2pid_2D( &dest_pid, D2, ind1-1, ind2-1);   \
    MPI_Isend (buffer_send, cnt, MPI_TYPE, dest_pid, tag,          MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv, cnt, MPI_TYPE, dest_pid, sym_tags-tag, MPI_COMM_WORLD, &handler_r);   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=0; i<halo_size; i++){   \
      for(int j=0; j<halo_size; j++){   \
        grid[i][j] = buffer_recv[cnt++];   \
      }   \
    }   \
    free(buffer_send);   \
    free(buffer_recv);   \
  }   \
  /* bottom right : isend tag 6 , irecv 2 = sym_tags-tag*/   \
  if( within_bound_2D(ind1+1, ind2+1, D1, D2) ){   \
    /*buffer*/   \
    Type * buffer_send = (Type *) malloc ( halo_size*halo_size*sizeof(Type) );   \
    Type * buffer_recv = (Type *) malloc ( halo_size*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;   \
    for(int i=halo_size+d1-halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size+d2-halo_size; j<halo_size+d2; j++){   \
        buffer_send[cnt++] = grid[i][j];   \
      }   \
    }   \
    /*isend && irecv*/   \
    tag=6;   \
    indice2pid_2D( &dest_pid, D2, ind1+1, ind2+1);   \
    MPI_Isend (buffer_send, cnt, MPI_TYPE, dest_pid, tag,          MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv, cnt, MPI_TYPE, dest_pid, sym_tags-tag, MPI_COMM_WORLD, &handler_r);   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=halo_size+d1; i<2*halo_size+d1; i++){   \
      for(int j=halo_size+d2; j<2*halo_size+d2; j++){   \
        grid[i][j] = buffer_recv[cnt++];   \
      }   \
    }   \
    /**/   \
    free(buffer_send);   \
    free(buffer_recv);   \
  }   \
  /* bottom left : isend tag 3 , irecv 5 = sym_tags-tag*/   \
  if( within_bound_2D(ind1+1, ind2-1, D1, D2) ){   \
    /*buffer*/   \
    Type * buffer_send = (Type *) malloc ( halo_size*halo_size*sizeof(Type) );   \
    Type * buffer_recv = (Type *) malloc ( halo_size*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;   \
    for(int i=halo_size+d1-halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size; j<2*halo_size; j++){   \
        buffer_send[cnt++] = grid[i][j];   \
      }   \
    }   \
    /*isend && irecv*/   \
    tag=3;   \
    indice2pid_2D( &dest_pid, D2, ind1+1, ind2-1);   \
    MPI_Isend (buffer_send, cnt, MPI_TYPE, dest_pid, tag,          MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv, cnt, MPI_TYPE, dest_pid, sym_tags-tag, MPI_COMM_WORLD, &handler_r);   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=halo_size+d1; i<2*halo_size+d1; i++){   \
      for(int j=0; j<halo_size; j++){   \
        grid[i][j] = buffer_recv[cnt++];   \
      }   \
    }   \
    /**/   \
    free(buffer_send);   \
    free(buffer_recv);   \
  }   \
  /* upper right : isend tag 5 , irecv 3 = sym_tags-tag*/   \
  if( within_bound_2D(ind1-1, ind2+1, D1, D2) ){   \
    /*buffer*/   \
    Type * buffer_send = (Type *) malloc ( halo_size*halo_size*sizeof(Type) );   \
    Type * buffer_recv = (Type *) malloc ( halo_size*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;   \
    for(int i=halo_size; i<2*halo_size; i++){   \
      for(int j=halo_size+d2-halo_size; j<halo_size+d2; j++){   \
        buffer_send[cnt++] = grid[i][j];   \
      }   \
    }   \
    /*isend && irecv*/   \
    tag=5;   \
    indice2pid_2D( &dest_pid, D2, ind1-1, ind2+1);   \
    MPI_Isend (buffer_send, cnt, MPI_TYPE, dest_pid, tag,          MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv, cnt, MPI_TYPE, dest_pid, sym_tags-tag, MPI_COMM_WORLD, &handler_r);   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=0; i<halo_size; i++){   \
      for(int j=halo_size+d2; j<2*halo_size+d2; j++){   \
        grid[i][j] = buffer_recv[cnt++];   \
      }   \
    }   \
    /**/   \
    free(buffer_send);   \
    free(buffer_recv);   \
  }   \
}

EXCHANGE_HALO_2D(float, MPI_FLOAT)
EXCHANGE_HALO_2D(double, MPI_DOUBLE)
EXCHANGE_HALO_2D(int, MPI_INT)

//blow is asyn
  
#define EXCHANGE_HALO_2D_START(Type, MPI_TYPE)\
void exchange_halo_2D_##Type##_start(Type * input, int dt, int d1, int d2, int halo_size, int pid, int D1, int D2, int mpi_size){\
  Type (*grid) [d2+2*halo_size] = (Type (*)[d2+2*halo_size]) input;   \
  int nei_ind1[9];   \
  int nei_ind2[9];   \
  int num_neighbors;   \
  int ind1;   \
  int ind2;   \
  pid2indice_2D(pid, D2, &ind1, &ind2);   \
  get_neighbors_2D(ind1, ind2, D1, D2, nei_ind1, nei_ind2, &num_neighbors);   \
  /* barrier*/ \
  MPI_Barrier(MPI_COMM_WORLD);      \
  /* delay*/ \
  for(int j=0;j<pid*pid*pid*10000;j++){      \
    double ans=(pid*4524523523454+100)/452352545425*pid+j;   \
  }   \
  /* print_grid*/   \
  print_grid_2D_##Type(input, d1, d2, halo_size, pid, nei_ind1, nei_ind2, num_neighbors);   \
  /* barrier*/ \
  MPI_Barrier(MPI_COMM_WORLD);      \
  /* start exchange halos*/   \
  const int sym_tags = 8;   \
  int dest_pid;   \
  int tag=0;   \
  /* upper : isend tag 0 , irecv 8 = sym_tags-tag*/   \
  if( within_bound_2D(ind1-1, ind2, D1, D2) ){   \
    /*buffer*/   \
    buffer_send_##Type[0] = (Type *) malloc ( d2*halo_size*sizeof(Type) );   \
    buffer_recv_##Type[0] = (Type *) malloc ( d2*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;   \
    for(int i=halo_size; i<2*halo_size; i++){   \
      for(int j=halo_size; j<halo_size+d2; j++){   \
        buffer_send_##Type[0][cnt++] = grid[i][j];   \
      }   \
    }   \
    /*isend && irecv*/   \
    tag=0;   \
    indice2pid_2D( &dest_pid, D2, ind1-1, ind2);   \
    MPI_Isend (buffer_send_##Type[0], cnt, MPI_TYPE, dest_pid, tag,          MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv_##Type[0], cnt, MPI_TYPE, dest_pid, sym_tags-tag, MPI_COMM_WORLD, &handler_r);   \
  }   \
  /* bottom : isend tag 8 , irecv 0 = sym_tags-tag*/   \
  if( within_bound_2D(ind1+1, ind2, D1, D2) ){   \
    /*buffer*/   \
    buffer_send_##Type[1] = (Type *) malloc ( d2*halo_size*sizeof(Type) );   \
    buffer_recv_##Type[1] = (Type *) malloc ( d2*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;   \
    for(int i=halo_size+d1-halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size; j<halo_size+d2; j++){   \
        buffer_send_##Type[1][cnt++] = grid[i][j];   \
      }   \
    }   \
    /*isend && irecv*/   \
    tag=8;   \
    indice2pid_2D( &dest_pid, D2, ind1+1, ind2);   \
    MPI_Isend (buffer_send_##Type[1], cnt, MPI_TYPE, dest_pid, tag,          MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv_##Type[1], cnt, MPI_TYPE, dest_pid, sym_tags-tag, MPI_COMM_WORLD, &handler_r);   \
  }   \
  /* left : isend tag 1 , irecv 7 = sym_tags-tag*/   \
  if( within_bound_2D(ind1, ind2-1, D1, D2) ){   \
    /*buffer*/   \
    buffer_send_##Type[2] = (Type *) malloc ( d1*halo_size*sizeof(Type) );   \
    buffer_recv_##Type[2] = (Type *) malloc ( d1*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;   \
    for(int i=halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size; j<2*halo_size; j++){   \
        buffer_send_##Type[2][cnt++] = grid[i][j];   \
      }   \
    }   \
    /*isend && irecv*/   \
    tag=1;   \
    indice2pid_2D( &dest_pid, D2, ind1, ind2-1);   \
    MPI_Isend (buffer_send_##Type[2], cnt, MPI_TYPE, dest_pid, tag,          MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv_##Type[2], cnt, MPI_TYPE, dest_pid, sym_tags-tag, MPI_COMM_WORLD, &handler_r);   \
  }   \
  /* right : isend tag 7 , irecv 1 = sym_tags-tag*/   \
  if( within_bound_2D(ind1, ind2+1, D1, D2) ){   \
    /*buffer*/   \
    buffer_send_##Type[3] = (Type *) malloc ( d1*halo_size*sizeof(Type) );   \
    buffer_recv_##Type[3] = (Type *) malloc ( d1*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;   \
    for(int i=halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size+d2-halo_size; j<halo_size+d2; j++){   \
        buffer_send_##Type[3][cnt++] = grid[i][j];   \
      }   \
    }   \
    /*isend && irecv*/   \
    tag=7;   \
    indice2pid_2D( &dest_pid, D2, ind1, ind2+1);   \
    MPI_Isend (buffer_send_##Type[3], cnt, MPI_TYPE, dest_pid, tag,          MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv_##Type[3], cnt, MPI_TYPE, dest_pid, sym_tags-tag, MPI_COMM_WORLD, &handler_r);   \
  }   \
  /* upper left : isend tag 2 , irecv 6 = sym_tags-tag*/   \
  if( within_bound_2D(ind1-1, ind2-1, D1, D2) ){   \
    /*buffer*/   \
    buffer_send_##Type[4] = (Type *) malloc ( halo_size*halo_size*sizeof(Type) );   \
    buffer_recv_##Type[4] = (Type *) malloc ( halo_size*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;   \
    for(int i=halo_size; i<2*halo_size; i++){   \
      for(int j=halo_size; j<2*halo_size; j++){   \
        buffer_send_##Type[4][cnt++] = grid[i][j];   \
      }   \
    }   \
    /*isend && irecv*/   \
    tag=2;   \
    indice2pid_2D( &dest_pid, D2, ind1-1, ind2-1);   \
    MPI_Isend (buffer_send_##Type[4], cnt, MPI_TYPE, dest_pid, tag,          MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv_##Type[4], cnt, MPI_TYPE, dest_pid, sym_tags-tag, MPI_COMM_WORLD, &handler_r);   \
  }   \
  /* bottom right : isend tag 6 , irecv 2 = sym_tags-tag*/   \
  if( within_bound_2D(ind1+1, ind2+1, D1, D2) ){   \
    /*buffer*/   \
    buffer_send_##Type[5] = (Type *) malloc ( halo_size*halo_size*sizeof(Type) );   \
    buffer_recv_##Type[5] = (Type *) malloc ( halo_size*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;   \
    for(int i=halo_size+d1-halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size+d2-halo_size; j<halo_size+d2; j++){   \
        buffer_send_##Type[5][cnt++] = grid[i][j];   \
      }   \
    }   \
    /*isend && irecv*/   \
    tag=6;   \
    indice2pid_2D( &dest_pid, D2, ind1+1, ind2+1);   \
    MPI_Isend (buffer_send_##Type[5], cnt, MPI_TYPE, dest_pid, tag,          MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv_##Type[5], cnt, MPI_TYPE, dest_pid, sym_tags-tag, MPI_COMM_WORLD, &handler_r);   \
  }   \
  /* bottom left : isend tag 3 , irecv 5 = sym_tags-tag*/   \
  if( within_bound_2D(ind1+1, ind2-1, D1, D2) ){   \
    /*buffer*/   \
    buffer_send_##Type[6] = (Type *) malloc ( halo_size*halo_size*sizeof(Type) );   \
    buffer_recv_##Type[6] = (Type *) malloc ( halo_size*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;   \
    for(int i=halo_size+d1-halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size; j<2*halo_size; j++){   \
        buffer_send_##Type[6][cnt++] = grid[i][j];   \
      }   \
    }   \
    /*isend && irecv*/   \
    tag=3;   \
    indice2pid_2D( &dest_pid, D2, ind1+1, ind2-1);   \
    MPI_Isend (buffer_send_##Type[6], cnt, MPI_TYPE, dest_pid, tag,          MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv_##Type[6], cnt, MPI_TYPE, dest_pid, sym_tags-tag, MPI_COMM_WORLD, &handler_r);   \
  }   \
  /* upper right : isend tag 5 , irecv 3 = sym_tags-tag*/   \
  if( within_bound_2D(ind1-1, ind2+1, D1, D2) ){   \
    /*buffer*/   \
    buffer_send_##Type[7] = (Type *) malloc ( halo_size*halo_size*sizeof(Type) );   \
    buffer_recv_##Type[7] = (Type *) malloc ( halo_size*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;   \
    for(int i=halo_size; i<2*halo_size; i++){   \
      for(int j=halo_size+d2-halo_size; j<halo_size+d2; j++){   \
        buffer_send_##Type[7][cnt++] = grid[i][j];   \
      }   \
    }   \
    /*isend && irecv*/   \
    tag=5;   \
    indice2pid_2D( &dest_pid, D2, ind1-1, ind2+1);   \
    MPI_Isend (buffer_send_##Type[7], cnt, MPI_TYPE, dest_pid, tag,          MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv_##Type[7], cnt, MPI_TYPE, dest_pid, sym_tags-tag, MPI_COMM_WORLD, &handler_r);   \
  }   \
}

EXCHANGE_HALO_2D_START(float, MPI_FLOAT)
EXCHANGE_HALO_2D_START(double, MPI_DOUBLE)
EXCHANGE_HALO_2D_START(int, MPI_INT)

 
#define EXCHANGE_HALO_2D_END(Type, MPI_TYPE)\
void exchange_halo_2D_##Type##_end(Type * input, int dt, int d1, int d2, int halo_size, int pid, int D1, int D2, int mpi_size){\
  Type (*grid) [d2+2*halo_size] = (Type (*)[d2+2*halo_size]) input;   \
  int nei_ind1[9];   \
  int nei_ind2[9];   \
  int num_neighbors;   \
  int ind1;   \
  int ind2;   \
  pid2indice_2D(pid, D2, &ind1, &ind2);   \
  get_neighbors_2D(ind1, ind2, D1, D2, nei_ind1, nei_ind2, &num_neighbors);   \
  /* start exchange halos*/   \
  const int sym_tags = 8;   \
  int dest_pid;   \
  int tag=0;   \
  int cnt=0;   \
  /* upper : isend tag 0 , irecv 8 = sym_tags-tag*/   \
  if( within_bound_2D(ind1-1, ind2, D1, D2) ){   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=0; i<halo_size; i++){   \
      for(int j=halo_size; j<halo_size+d2; j++){   \
        grid[i][j] = buffer_recv_##Type[0][cnt++];   \
      }   \
    }   \
    free(buffer_send_##Type[0]);   \
    free(buffer_recv_##Type[0]);   \
  }   \
  /* bottom : isend tag 8 , irecv 0 = sym_tags-tag*/   \
  if( within_bound_2D(ind1+1, ind2, D1, D2) ){   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=halo_size+d1; i<2*halo_size+d1; i++){   \
      for(int j=halo_size; j<halo_size+d2; j++){   \
        grid[i][j] = buffer_recv_##Type[1][cnt++];   \
      }   \
    }   \
    /**/   \
    free(buffer_send_##Type[1]);   \
    free(buffer_recv_##Type[1]);   \
  }   \
  /* left : isend tag 1 , irecv 7 = sym_tags-tag*/   \
  if( within_bound_2D(ind1, ind2-1, D1, D2) ){   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=halo_size; i<halo_size+d1; i++){   \
      for(int j=0; j<halo_size; j++){   \
        grid[i][j] = buffer_recv_##Type[2][cnt++];   \
      }   \
    }   \
    free(buffer_send_##Type[2]);   \
    free(buffer_recv_##Type[2]);   \
  }   \
  /* right : isend tag 7 , irecv 1 = sym_tags-tag*/   \
  if( within_bound_2D(ind1, ind2+1, D1, D2) ){   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size+d2; j<2*halo_size+d2; j++){   \
        grid[i][j] = buffer_recv_##Type[3][cnt++];   \
      }   \
    }   \
    /**/   \
    free(buffer_send_##Type[3]);   \
    free(buffer_recv_##Type[3]);   \
  }   \
  /* upper left : isend tag 2 , irecv 6 = sym_tags-tag*/   \
  if( within_bound_2D(ind1-1, ind2-1, D1, D2) ){   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=0; i<halo_size; i++){   \
      for(int j=0; j<halo_size; j++){   \
        grid[i][j] = buffer_recv_##Type[4][cnt++];   \
      }   \
    }   \
    free(buffer_send_##Type[4]);   \
    free(buffer_recv_##Type[4]);   \
  }   \
  /* bottom right : isend tag 6 , irecv 2 = sym_tags-tag*/   \
  if( within_bound_2D(ind1+1, ind2+1, D1, D2) ){   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=halo_size+d1; i<2*halo_size+d1; i++){   \
      for(int j=halo_size+d2; j<2*halo_size+d2; j++){   \
        grid[i][j] = buffer_recv_##Type[5][cnt++];   \
      }   \
    }   \
    /**/   \
    free(buffer_send_##Type[5]);   \
    free(buffer_recv_##Type[5]);   \
  }   \
  /* bottom left : isend tag 3 , irecv 5 = sym_tags-tag*/   \
  if( within_bound_2D(ind1+1, ind2-1, D1, D2) ){   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=halo_size+d1; i<2*halo_size+d1; i++){   \
      for(int j=0; j<halo_size; j++){   \
        grid[i][j] = buffer_recv_##Type[6][cnt++];   \
      }   \
    }   \
    /**/   \
    free(buffer_send_##Type[6]);   \
    free(buffer_recv_##Type[6]);   \
  }   \
  /* upper right : isend tag 5 , irecv 3 = sym_tags-tag*/   \
  if( within_bound_2D(ind1-1, ind2+1, D1, D2) ){   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=0; i<halo_size; i++){   \
      for(int j=halo_size+d2; j<2*halo_size+d2; j++){   \
        grid[i][j] = buffer_recv_##Type[7][cnt++];   \
      }   \
    }   \
    /**/   \
    free(buffer_send_##Type[7]);   \
    free(buffer_recv_##Type[7]);   \
  }   \
  /* barrier*/   \
  MPI_Barrier(MPI_COMM_WORLD);   \
  /* delay*/   \
  for(int j=0;j<pid*pid*pid*10000;j++){   \
    double ans=(pid*4524523523454+100)/452352545425*pid+j;   \
  }   \
  /* print_grid*/   \
  print_grid_2D_##Type(input, d1, d2, halo_size, pid, nei_ind1, nei_ind2, num_neighbors);   \
}

EXCHANGE_HALO_2D_END(float, MPI_FLOAT)
EXCHANGE_HALO_2D_END(double, MPI_DOUBLE)
EXCHANGE_HALO_2D_END(int, MPI_INT)


// 3D
 
void pid2indice_3D(const int pid, const int D2, const int D3, int *ind1, int *ind2, int *ind3){
  *ind1 = pid / (D2*D3);
  *ind2 = (pid % (D2*D3) ) / D3;
  *ind3 = pid % D3;
}

void indice2pid_3D(int *pid, const int D2, const int D3, const int ind1, const int ind2, const int ind3){
  *pid = ind1 * D2 * D3 + ind2 * D3 + ind3;
}

int within_bound_3D(const int ind1, const int ind2, const int ind3, const int D1, const int D2, const int D3){
  if( 0<=ind1 && ind1<D1 && 0<=ind2 && ind2<D2 && 0<=ind3 && ind3<D3 ) return 1;
  return 0;
}



#define PRINT_GRID_3D_BARRIER(Type, flag)   \
void print_grid_3D_##Type##_barrier(Type * input, int d1, int d2, int d3, int halo_size, const int pid){\
  /*barrier*/       \
  MPI_Barrier(MPI_COMM_WORLD);      \
  /* delay*/      \
  for(int j=0;j<pid*pid*pid*10000;j++){             \
    double ans=(pid*4524523523454+100)/452352545425*pid+j;    \
  }                                         \
  /* print_grid*/                                         \
  printf("<------------------------->\n");    \
  printf("from process : %d\n", pid);       \
  Type (*grid) [d2+2*halo_size][d3+2*halo_size] = (Type (*)[d2+2*halo_size][d3+2*halo_size]) input;\
  printf("{");                                       \
  for(int j=0;j<d1+2*halo_size;j++){                    \
    printf("{");                                       \
    for(int k=0;k<d2+2*halo_size;k++){                  \
      printf("{");                                       \
      for(int p=0;p<d3+2*halo_size;p++){\
        printf( #flag , grid[j][k][p]);                  \
        printf( " " );                                    \
      }\
      printf("}\n");                                       \
    }                                                   \
    printf("}\n");                                       \
  }                                                     \
  printf("}\n");                                       \
}                                      

PRINT_GRID_3D_BARRIER(float, %.f)
PRINT_GRID_3D_BARRIER(double, %.f)
PRINT_GRID_3D_BARRIER(int, %d)


#define EXCHANGE_HALO_3D(Type, MPI_TYPE)\
void exchange_halo_3D_##Type(Type * input, int dt, int d1, int d2, int d3, int halo_size, int pid, int D1, int D2, int D3, int mpi_size){\
  Type (*grid) [d2+2*halo_size][d3+2*halo_size] = (Type (*)[d2+2*halo_size][d3+2*halo_size]) input;   \
  int ind1;   \
  int ind2;   \
  int ind3;   \
  pid2indice_3D(pid, D2, D3, &ind1, &ind2, &ind3);   \
  /* start exchange halos*/   \
  int dest_pid;   \
  int tag=0;   \
  /* front*/   \
  if( within_bound_3D(ind1, ind2-1, ind3, D1, D2, D3) ){   \
    /*buffer*/      \
    buffer_send_##Type[0] = (Type *) malloc ( d1*d3*halo_size*sizeof(Type) );\
    buffer_recv_##Type[0] = (Type *) malloc ( d1*d3*halo_size*sizeof(Type) );\
    /*pack*/\
    int cnt=0;\
    for(int i=halo_size; i<halo_size+d1; i++){\
      for(int j=halo_size; j<2*halo_size; j++){\
        for(int k=halo_size; k<halo_size+d3; k++){   \
          buffer_send_##Type[0][cnt++] = grid[i][j][k];   \
        }\
      }   \
    }   \
    /*isend && irecv*/   \
    indice2pid_3D( &dest_pid, D2, D3, ind1, ind2-1, ind3);   \
    MPI_Isend (buffer_send_##Type[0], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv_##Type[0], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_r);   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=halo_size; i<halo_size+d1; i++){   \
      for(int j=0; j<halo_size; j++){   \
        for(int k=halo_size; k<halo_size+d3; k++){   \
           grid[i][j][k] = buffer_recv_##Type[0][cnt++];   \
        }\
      }   \
    }   \
    free(buffer_send_##Type[0]);   \
    free(buffer_recv_##Type[0]);  \
  }\
  /* back*/\
  if( within_bound_3D(ind1, ind2+1, ind3, D1, D2, D3) ){\
    /*buffer*/   \
    buffer_send_##Type[1] = (Type *) malloc ( d1*d3*halo_size*sizeof(Type) );   \
    buffer_recv_##Type[1] = (Type *) malloc ( d1*d3*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;\
    for(int i=halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size+d2-halo_size; j<halo_size+d2; j++){   \
        for(int k=halo_size; k<halo_size+d3; k++){   \
          buffer_send_##Type[1][cnt++] = grid[i][j][k];   \
        }\
      }  \
    }   \
    /*isend && irecv*/   \
    indice2pid_3D( &dest_pid, D2, D3, ind1, ind2+1, ind3);   \
    MPI_Isend (buffer_send_##Type[1], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv_##Type[1], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_r);   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size+d2; j<2*halo_size+d2; j++){   \
        for(int k=halo_size; k<halo_size+d3; k++){   \
           grid[i][j][k] = buffer_recv_##Type[1][cnt++];   \
        }\
      }   \
    }   \
    free(buffer_send_##Type[1]);   \
    free(buffer_recv_##Type[1]);  \
  }\
  /* left*/\
  if( within_bound_3D(ind1, ind2, ind3-1, D1, D2, D3) ){\
    /*buffer*/   \
    buffer_send_##Type[2] = (Type *) malloc ( d1*d2*halo_size*sizeof(Type) );   \
    buffer_recv_##Type[2] = (Type *) malloc ( d1*d2*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;   \
    for(int i=halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size; j<halo_size+d2; j++){   \
        for(int k=halo_size; k<2*halo_size; k++){   \
          buffer_send_##Type[2][cnt++] = grid[i][j][k];   \
        }\
      }   \
    }   \
    /*isend && irecv*/   \
    indice2pid_3D( &dest_pid, D2, D3, ind1, ind2, ind3-1);   \
    MPI_Isend (buffer_send_##Type[2], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv_##Type[2], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_r);   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size; j<halo_size+d2; j++){  \
        for(int k=0; k<halo_size; k++){   \
           grid[i][j][k] = buffer_recv_##Type[2][cnt++];   \
        }\
      }   \
    }   \
    free(buffer_send_##Type[2]);   \
    free(buffer_recv_##Type[2]);  \
  }\
  /* right*/    \
  if( within_bound_3D(ind1, ind2, ind3+1, D1, D2, D3) ){    \
    /*buffer*/       \
    buffer_send_##Type[3] = (Type *) malloc ( d1*d2*halo_size*sizeof(Type) );       \
    buffer_recv_##Type[3] = (Type *) malloc ( d1*d2*halo_size*sizeof(Type) );       \
    /*pack*/       \
    int cnt=0;       \
    for(int i=halo_size; i<halo_size+d1; i++){       \
      for(int j=halo_size; j<halo_size+d2; j++){       \
        for(int k=halo_size+d3-halo_size; k<halo_size+d3; k++){       \
          buffer_send_##Type[3][cnt++] = grid[i][j][k];       \
        }    \
      }       \
    }       \
    /*isend && irecv*/       \
    indice2pid_3D( &dest_pid, D2, D3, ind1, ind2, ind3+1);       \
    MPI_Isend (buffer_send_##Type[3], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_s);       \
    MPI_Irecv (buffer_recv_##Type[3], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_r);       \
    /*wait*/       \
    MPI_Wait(&handler_s, &status_s);       \
    MPI_Wait(&handler_r, &status_r);       \
    /*unpack*/       \
    cnt=0;       \
    for(int i=halo_size; i<halo_size+d1; i++){       \
      for(int j=halo_size; j<halo_size+d2; j++){       \
        for(int k=halo_size+d3; k<2*halo_size+d3; k++){       \
           grid[i][j][k] = buffer_recv_##Type[3][cnt++];       \
        }    \
      }       \
    }       \
    free(buffer_send_##Type[3]);       \
    free(buffer_recv_##Type[3]);      \
  }    \
  /* up*/    \
  if( within_bound_3D(ind1+1, ind2, ind3, D1, D2, D3) ){    \
    /*buffer*/       \
    buffer_send_##Type[4] = (Type *) malloc ( d2*d3*halo_size*sizeof(Type) );       \
    buffer_recv_##Type[4] = (Type *) malloc ( d2*d3*halo_size*sizeof(Type) );       \
    /*pack*/       \
    int cnt=0;       \
    for(int i=halo_size+d1-halo_size; i<halo_size+d1; i++){       \
      for(int j=halo_size; j<halo_size+d2; j++){       \
        for(int k=halo_size; k<halo_size+d3; k++){       \
          buffer_send_##Type[4][cnt++] = grid[i][j][k];       \
        }    \
      }       \
    }       \
    /*isend && irecv*/       \
    indice2pid_3D( &dest_pid, D2, D3, ind1+1, ind2, ind3);       \
    MPI_Isend (buffer_send_##Type[4], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_s);       \
    MPI_Irecv (buffer_recv_##Type[4], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_r);       \
    /*wait*/       \
    MPI_Wait(&handler_s, &status_s);       \
    MPI_Wait(&handler_r, &status_r);       \
    /*unpack*/       \
    cnt=0;       \
    for(int i=halo_size+d1; i<2*halo_size+d1; i++){       \
      for(int j=halo_size; j<halo_size+d2; j++){       \
        for(int k=halo_size; k<halo_size+d3; k++){       \
           grid[i][j][k] = buffer_recv_##Type[4][cnt++];       \
        }    \
      }       \
    }       \
    free(buffer_send_##Type[4]);       \
    free(buffer_recv_##Type[4]);      \
  }    \
  /* down*/    \
  if( within_bound_3D(ind1-1, ind2, ind3, D1, D2, D3) ){    \
    /*buffer*/       \
    buffer_send_##Type[5] = (Type *) malloc ( d2*d3*halo_size*sizeof(Type) );       \
    buffer_recv_##Type[5] = (Type *) malloc ( d2*d3*halo_size*sizeof(Type) );       \
    /*pack*/       \
    int cnt=0;       \
    for(int i=halo_size; i<2*halo_size; i++){       \
      for(int j=halo_size; j<halo_size+d2; j++){       \
        for(int k=halo_size; k<halo_size+d3; k++){       \
          buffer_send_##Type[5][cnt++] = grid[i][j][k];       \
        }    \
      }       \
    }       \
    /*isend && irecv*/       \
    indice2pid_3D( &dest_pid, D2, D3, ind1-1, ind2, ind3);       \
    MPI_Isend (buffer_send_##Type[5], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_s);       \
    MPI_Irecv (buffer_recv_##Type[5], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_r);       \
    /*wait*/       \
    MPI_Wait(&handler_s, &status_s);       \
    MPI_Wait(&handler_r, &status_r);       \
    /*unpack*/       \
    cnt=0;       \
    for(int i=0; i<halo_size; i++){       \
      for(int j=halo_size; j<halo_size+d2; j++){       \
        for(int k=halo_size; k<halo_size+d3; k++){       \
           grid[i][j][k] = buffer_recv_##Type[5][cnt++];       \
        }    \
      }       \
    }       \
    free(buffer_send_##Type[5]);       \
    free(buffer_recv_##Type[5]);      \
  }    \
}

EXCHANGE_HALO_3D(float, MPI_FLOAT)
EXCHANGE_HALO_3D(double, MPI_DOUBLE)
EXCHANGE_HALO_3D(int, MPI_INT)


#define EXCHANGE_HALO_3D_START(Type, MPI_TYPE)\
void exchange_halo_3D_##Type##_start(Type * input, int dt, int d1, int d2, int d3, int halo_size, int pid, int D1, int D2, int D3, int mpi_size){\
  Type (*grid) [d2+2*halo_size][d3+2*halo_size] = (Type (*)[d2+2*halo_size][d3+2*halo_size]) input;   \
  int ind1;   \
  int ind2;   \
  int ind3;   \
  pid2indice_3D(pid, D2, D3, &ind1, &ind2, &ind3);   \
  /* print_grid*/   \
  print_grid_3D_##Type##_barrier(input, d1, d2, d3, halo_size, pid);   \
  /* start exchange halos*/   \
  int dest_pid;   \
  int tag=0;   \
  /* front*/   \
  if( within_bound_3D(ind1, ind2-1, ind3, D1, D2, D3) ){   \
    /*buffer*/      \
    buffer_send_##Type[0] = (Type *) malloc ( d1*d3*halo_size*sizeof(Type) );\
    buffer_recv_##Type[0] = (Type *) malloc ( d1*d3*halo_size*sizeof(Type) );\
    /*pack*/\
    int cnt=0;\
    for(int i=halo_size; i<halo_size+d1; i++){\
      for(int j=halo_size; j<2*halo_size; j++){\
        for(int k=halo_size; k<halo_size+d3; k++){   \
          buffer_send_##Type[0][cnt++] = grid[i][j][k];   \
        }\
      }   \
    }   \
    /*isend && irecv*/   \
    indice2pid_3D( &dest_pid, D2, D3, ind1, ind2-1, ind3);   \
    MPI_Isend (buffer_send_##Type[0], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv_##Type[0], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_r);   \
  }\
  /* back*/\
  if( within_bound_3D(ind1, ind2+1, ind3, D1, D2, D3) ){\
    /*buffer*/   \
    buffer_send_##Type[1] = (Type *) malloc ( d1*d3*halo_size*sizeof(Type) );   \
    buffer_recv_##Type[1] = (Type *) malloc ( d1*d3*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;\
    for(int i=halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size+d2-halo_size; j<halo_size+d2; j++){   \
        for(int k=halo_size; k<halo_size+d3; k++){   \
          buffer_send_##Type[1][cnt++] = grid[i][j][k];   \
        }\
      }  \
    }   \
    /*isend && irecv*/   \
    indice2pid_3D( &dest_pid, D2, D3, ind1, ind2+1, ind3);   \
    MPI_Isend (buffer_send_##Type[1], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv_##Type[1], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_r);   \
  }\
  /* left*/\
  if( within_bound_3D(ind1, ind2, ind3-1, D1, D2, D3) ){\
    /*buffer*/   \
    buffer_send_##Type[2] = (Type *) malloc ( d1*d2*halo_size*sizeof(Type) );   \
    buffer_recv_##Type[2] = (Type *) malloc ( d1*d2*halo_size*sizeof(Type) );   \
    /*pack*/   \
    int cnt=0;   \
    for(int i=halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size; j<halo_size+d2; j++){   \
        for(int k=halo_size; k<2*halo_size; k++){   \
          buffer_send_##Type[2][cnt++] = grid[i][j][k];   \
        }\
      }   \
    }   \
    /*isend && irecv*/   \
    indice2pid_3D( &dest_pid, D2, D3, ind1, ind2, ind3-1);   \
    MPI_Isend (buffer_send_##Type[2], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_s);   \
    MPI_Irecv (buffer_recv_##Type[2], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_r);   \
  }\
  /* right*/    \
  if( within_bound_3D(ind1, ind2, ind3+1, D1, D2, D3) ){    \
    /*buffer*/       \
    buffer_send_##Type[3] = (Type *) malloc ( d1*d2*halo_size*sizeof(Type) );       \
    buffer_recv_##Type[3] = (Type *) malloc ( d1*d2*halo_size*sizeof(Type) );       \
    /*pack*/       \
    int cnt=0;       \
    for(int i=halo_size; i<halo_size+d1; i++){       \
      for(int j=halo_size; j<halo_size+d2; j++){       \
        for(int k=halo_size+d3-halo_size; k<halo_size+d3; k++){       \
          buffer_send_##Type[3][cnt++] = grid[i][j][k];       \
        }    \
      }       \
    }       \
    /*isend && irecv*/       \
    indice2pid_3D( &dest_pid, D2, D3, ind1, ind2, ind3+1);       \
    MPI_Isend (buffer_send_##Type[3], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_s);       \
    MPI_Irecv (buffer_recv_##Type[3], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_r);       \
  }    \
  /* up*/    \
  if( within_bound_3D(ind1+1, ind2, ind3, D1, D2, D3) ){    \
    /*buffer*/       \
    buffer_send_##Type[4] = (Type *) malloc ( d2*d3*halo_size*sizeof(Type) );       \
    buffer_recv_##Type[4] = (Type *) malloc ( d2*d3*halo_size*sizeof(Type) );       \
    /*pack*/       \
    int cnt=0;       \
    for(int i=halo_size+d1-halo_size; i<halo_size+d1; i++){       \
      for(int j=halo_size; j<halo_size+d2; j++){       \
        for(int k=halo_size; k<halo_size+d3; k++){       \
          buffer_send_##Type[4][cnt++] = grid[i][j][k];       \
        }    \
      }       \
    }       \
    /*isend && irecv*/       \
    indice2pid_3D( &dest_pid, D2, D3, ind1+1, ind2, ind3);       \
    MPI_Isend (buffer_send_##Type[4], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_s);       \
    MPI_Irecv (buffer_recv_##Type[4], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_r);       \
  }    \
  /* down*/    \
  if( within_bound_3D(ind1-1, ind2, ind3, D1, D2, D3) ){    \
    /*buffer*/       \
    buffer_send_##Type[5] = (Type *) malloc ( d2*d3*halo_size*sizeof(Type) );       \
    buffer_recv_##Type[5] = (Type *) malloc ( d2*d3*halo_size*sizeof(Type) );       \
    /*pack*/       \
    int cnt=0;       \
    for(int i=halo_size; i<2*halo_size; i++){       \
      for(int j=halo_size; j<halo_size+d2; j++){       \
        for(int k=halo_size; k<halo_size+d3; k++){       \
          buffer_send_##Type[5][cnt++] = grid[i][j][k];       \
        }    \
      }       \
    }       \
    /*isend && irecv*/       \
    indice2pid_3D( &dest_pid, D2, D3, ind1-1, ind2, ind3);       \
    MPI_Isend (buffer_send_##Type[5], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_s);       \
    MPI_Irecv (buffer_recv_##Type[5], cnt, MPI_TYPE, dest_pid, tag, MPI_COMM_WORLD, &handler_r);       \
  }    \
}


EXCHANGE_HALO_3D_START(float, MPI_FLOAT)
EXCHANGE_HALO_3D_START(double, MPI_DOUBLE)
EXCHANGE_HALO_3D_START(int, MPI_INT)

#define EXCHANGE_HALO_3D_END(Type, MPI_TYPE)\
void exchange_halo_3D_##Type##_end(Type * input, int dt, int d1, int d2, int d3, int halo_size, int pid, int D1, int D2, int D3, int mpi_size){\
  Type (*grid) [d2+2*halo_size][d3+2*halo_size] = (Type (*)[d2+2*halo_size][d3+2*halo_size]) input;   \
  int ind1;   \
  int ind2;   \
  int ind3;   \
  pid2indice_3D(pid, D2, D3, &ind1, &ind2, &ind3);   \
  /* start exchange halos*/   \
  int dest_pid;   \
  int tag=0;   \
  int cnt=0;   \
  /* front*/   \
  if( within_bound_3D(ind1, ind2-1, ind3, D1, D2, D3) ){   \
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=halo_size; i<halo_size+d1; i++){   \
      for(int j=0; j<halo_size; j++){   \
        for(int k=halo_size; k<halo_size+d3; k++){   \
           grid[i][j][k] = buffer_recv_##Type[0][cnt++];   \
        }\
      }   \
    }   \
    free(buffer_send_##Type[0]);   \
    free(buffer_recv_##Type[0]);  \
  }\
  /* back*/\
  if( within_bound_3D(ind1, ind2+1, ind3, D1, D2, D3) ){\
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size+d2; j<2*halo_size+d2; j++){   \
        for(int k=halo_size; k<halo_size+d3; k++){   \
           grid[i][j][k] = buffer_recv_##Type[1][cnt++];   \
        }\
      }   \
    }   \
    free(buffer_send_##Type[1]);   \
    free(buffer_recv_##Type[1]);  \
  }\
  /* left*/\
  if( within_bound_3D(ind1, ind2, ind3-1, D1, D2, D3) ){\
    /*wait*/   \
    MPI_Wait(&handler_s, &status_s);   \
    MPI_Wait(&handler_r, &status_r);   \
    /*unpack*/   \
    cnt=0;   \
    for(int i=halo_size; i<halo_size+d1; i++){   \
      for(int j=halo_size; j<halo_size+d2; j++){  \
        for(int k=0; k<halo_size; k++){   \
           grid[i][j][k] = buffer_recv_##Type[2][cnt++];   \
        }\
      }   \
    }   \
    free(buffer_send_##Type[2]);   \
    free(buffer_recv_##Type[2]);  \
  }\
  /* right*/    \
  if( within_bound_3D(ind1, ind2, ind3+1, D1, D2, D3) ){    \
    /*wait*/       \
    MPI_Wait(&handler_s, &status_s);       \
    MPI_Wait(&handler_r, &status_r);       \
    /*unpack*/       \
    cnt=0;       \
    for(int i=halo_size; i<halo_size+d1; i++){       \
      for(int j=halo_size; j<halo_size+d2; j++){       \
        for(int k=halo_size+d3; k<2*halo_size+d3; k++){       \
           grid[i][j][k] = buffer_recv_##Type[3][cnt++];       \
        }    \
      }       \
    }       \
    free(buffer_send_##Type[3]);       \
    free(buffer_recv_##Type[3]);      \
  }    \
  /* up*/    \
  if( within_bound_3D(ind1+1, ind2, ind3, D1, D2, D3) ){    \
    /*wait*/       \
    MPI_Wait(&handler_s, &status_s);       \
    MPI_Wait(&handler_r, &status_r);       \
    /*unpack*/       \
    cnt=0;       \
    for(int i=halo_size+d1; i<2*halo_size+d1; i++){       \
      for(int j=halo_size; j<halo_size+d2; j++){       \
        for(int k=halo_size; k<halo_size+d3; k++){       \
           grid[i][j][k] = buffer_recv_##Type[4][cnt++];       \
        }    \
      }       \
    }       \
    free(buffer_send_##Type[4]);       \
    free(buffer_recv_##Type[4]);      \
  }    \
  /* down*/    \
  if( within_bound_3D(ind1-1, ind2, ind3, D1, D2, D3) ){    \
    /*wait*/       \
    MPI_Wait(&handler_s, &status_s);       \
    MPI_Wait(&handler_r, &status_r);       \
    /*unpack*/       \
    cnt=0;       \
    for(int i=0; i<halo_size; i++){       \
      for(int j=halo_size; j<halo_size+d2; j++){       \
        for(int k=halo_size; k<halo_size+d3; k++){       \
           grid[i][j][k] = buffer_recv_##Type[5][cnt++];       \
        }    \
      }       \
    }       \
    free(buffer_send_##Type[5]);       \
    free(buffer_recv_##Type[5]);      \
  }    \
  /* print_grid*/    \
  print_grid_3D_##Type##_barrier(input, d1, d2, d3, halo_size, pid);    \
}

EXCHANGE_HALO_3D_END(float, MPI_FLOAT)
EXCHANGE_HALO_3D_END(double, MPI_DOUBLE)
EXCHANGE_HALO_3D_END(int, MPI_INT)



