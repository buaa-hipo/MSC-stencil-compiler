#include "mpi_lib.h"

void test(int mpi_rank, int mpi_size){
  const int halo_size=1;
  const int d1=3, d2=4;
  const int temporal_size=2;
  const int d1h=d1+2*halo_size;
  const int d2h=d2+2*halo_size;
  float input[2][5][6]=
  {
    {
      {0,0,0,0,0,0},
      {0,1,1,1,1,0},
      {0,1,2,2,1,0},
      {0,1,1,1,1,0},
      {0,0,0,0,0,0}
    },
   
    {
      {0,0,0,0,0,0},
      {0,5,1,1,8,0},
      {0,3,8,8,4,0},
      {0,7,2,2,6,0},
      {0,0,0,0,0,0}
    }
  };
  int pid=mpi_rank;
  int D1=2;
  int D2=mpi_size/D1;
  exchange_halo_2D_float(&input[1][0][0], temporal_size, d1, d2 , halo_size, pid, D1, D2, mpi_size);
}

void test2(int mpi_rank, int mpi_size){
  const int halo_size=2;
  const int d1=5, d2=6;
  const int temporal_size=1;
  const int d1h=d1+2*halo_size;
  const int d2h=d2+2*halo_size;
  float input[1][9][10]=
  {
    {
      {0,0, 0, 0, 0, 0, 0, 0,0,0},
      {0,0, 0, 0, 0, 0, 0, 0,0,0},
      {0,0, 1, 2,17,18,13,14,0,0},
      {0,0, 3, 4,19,20,15,16,0,0},
      {0,0,25,26,29,30,27,28,0,0},
      {0,0, 9,10,21,22, 5, 6,0,0},
      {0,0,11,12,23,24, 7, 8,0,0},
      {0,0, 0, 0, 0, 0, 0, 0,0,0},
      {0,0, 0, 0, 0, 0, 0, 0,0,0}
    }
  };
  int pid=mpi_rank;
  int D1=2;
  int D2=mpi_size/D1;
  exchange_halo_2D_float(&input[0][0][0], temporal_size, d1, d2 , halo_size, pid, D1, D2, mpi_size);
  for(int i=0;i<temporal_size;i++){
    for(int j=0;j<d1h;j++){
      for(int k=0;k<d2h;k++){
        if(halo_size<=j&&j<halo_size+d1 &&halo_size<=k&&k<halo_size+d2)
        input[i][j][k]+=1;
      }
    }
  }
  exchange_halo_2D_float(&input[0][0][0], temporal_size, d1, d2 , halo_size, pid, D1, D2, mpi_size);
}

void test2_double(int mpi_rank, int mpi_size){
  const int halo_size=2;
  const int d1=5, d2=6;
  const int temporal_size=1;
  const int d1h=d1+2*halo_size;
  const int d2h=d2+2*halo_size;
  double input[1][9][10]=
  {
    {
      {0,0, 0, 0, 0, 0, 0, 0,0,0},
      {0,0, 0, 0, 0, 0, 0, 0,0,0},
      {0,0, 1, 2,17,18,13,14,0,0},
      {0,0, 3, 4,19,20,15,16,0,0},
      {0,0,25,26,29,30,27,28,0,0},
      {0,0, 9,10,21,22, 5, 6,0,0},
      {0,0,11,12,23,24, 7, 8,0,0},
      {0,0, 0, 0, 0, 0, 0, 0,0,0},
      {0,0, 0, 0, 0, 0, 0, 0,0,0}
    }
  };
  int pid=mpi_rank;
  int D1=2;
  int D2=mpi_size/D1;
  exchange_halo_2D_double(&input[0][0][0], temporal_size, d1, d2 , halo_size, pid, D1, D2, mpi_size);
  for(int i=0;i<temporal_size;i++){
    for(int j=0;j<d1h;j++){
      for(int k=0;k<d2h;k++){
        if(halo_size<=j&&j<halo_size+d1 &&halo_size<=k&&k<halo_size+d2)
        input[i][j][k]+=1;
      }
    }
  }
  exchange_halo_2D_double(&input[0][0][0], temporal_size, d1, d2 , halo_size, pid, D1, D2, mpi_size);
}

void test2_int(int mpi_rank, int mpi_size){
  const int halo_size=2;
  const int d1=5, d2=6;
  const int temporal_size=1;
  const int d1h=d1+2*halo_size;
  const int d2h=d2+2*halo_size;
  int input[1][9][10]=
  {
    {
      {0,0, 0, 0, 0, 0, 0, 0,0,0},
      {0,0, 0, 0, 0, 0, 0, 0,0,0},
      {0,0, 1, 2,17,18,13,14,0,0},
      {0,0, 3, 4,19,20,15,16,0,0},
      {0,0,25,26,29,30,27,28,0,0},
      {0,0, 9,10,21,22, 5, 6,0,0},
      {0,0,11,12,23,24, 7, 8,0,0},
      {0,0, 0, 0, 0, 0, 0, 0,0,0},
      {0,0, 0, 0, 0, 0, 0, 0,0,0}
    }
  };
  int pid=mpi_rank;
  int D1=2;
  int D2=mpi_size/D1;
  exchange_halo_2D_int(&input[0][0][0], temporal_size, d1, d2 , halo_size, pid, D1, D2, mpi_size);
  for(int i=0;i<temporal_size;i++){
    for(int j=0;j<d1h;j++){
      for(int k=0;k<d2h;k++){
        if(halo_size<=j&&j<halo_size+d1 &&halo_size<=k&&k<halo_size+d2)
        input[i][j][k]+=1;
      }
    }
  }
  exchange_halo_2D_int(&input[0][0][0], temporal_size, d1, d2 , halo_size, pid, D1, D2, mpi_size);
}

void test_int_asyn(int mpi_rank, int mpi_size){
  const int halo_size=2;
  const int d1=5, d2=6;
  const int temporal_size=1;
  const int d1h=d1+2*halo_size;
  const int d2h=d2+2*halo_size;
  int input[1][9][10]=
  {
    {
      {0,0, 0, 0, 0, 0, 0, 0,0,0},
      {0,0, 0, 0, 0, 0, 0, 0,0,0},
      {0,0, 1, 2,17,18,13,14,0,0},
      {0,0, 3, 4,19,20,15,16,0,0},
      {0,0,25,26,29,30,27,28,0,0},
      {0,0, 9,10,21,22, 5, 6,0,0},
      {0,0,11,12,23,24, 7, 8,0,0},
      {0,0, 0, 0, 0, 0, 0, 0,0,0},
      {0,0, 0, 0, 0, 0, 0, 0,0,0}
    }
  };
  int pid=mpi_rank;
  int D1=2;
  int D2=mpi_size/D1;
  exchange_halo_2D_int_start(&input[0][0][0], temporal_size, d1, d2 , halo_size, pid, D1, D2, mpi_size);
  for(int i=0;i<temporal_size;i++){
    for(int j=0;j<d1h;j++){
      for(int k=0;k<d2h;k++){
        if(2*halo_size<=j&&j<halo_size+d1-halo_size &&2*halo_size<=k&&k<halo_size+d2-halo_size)
        input[i][j][k]+=1;
      }
    }
  }
  exchange_halo_2D_int_end(&input[0][0][0], temporal_size, d1, d2 , halo_size, pid, D1, D2, mpi_size);
  for(int i=0;i<temporal_size;i++){
    for(int j=0;j<d1h;j++){
      for(int k=0;k<d2h;k++){
        if(halo_size<=j&&j<halo_size+d1 &&halo_size<=k&&k<halo_size+d2){
          if(2*halo_size<=j&&j<halo_size+d1-halo_size &&2*halo_size<=k&&k<halo_size+d2-halo_size){
            continue;
          }
          else
          {
            input[i][j][k]-=1;
          }
        }
      }
    }
  }
  print_grid_2D_int_barrier(&input[0][0][0], d1, d2, halo_size, mpi_rank);
}


void test_float_1D(int mpi_rank, int mpi_size){
  const int halo_size=2;
  const int d1=6;
  const int temporal_size=1;
  const int d1h=d1+2*halo_size;
  float input[1][10]=
  {
    {0,0,1,2,3,4,5,6,0,0}
  };
  int pid=mpi_rank;
  int D1=mpi_size;
  exchange_halo_1D_float(&input[0][0], temporal_size, d1, halo_size, pid, D1, mpi_size);
  for(int i=0;i<temporal_size;i++){
    for(int j=0;j<d1h;j++){
      if(2*halo_size<=j&&j<halo_size+d1-halo_size)
        input[i][j]+=1;
    }
  }
  for(int i=0;i<temporal_size;i++){
    for(int j=0;j<d1h;j++){
        if(halo_size<=j&&j<halo_size+d1){
          if(2*halo_size<=j&&j<halo_size+d1-halo_size){
            continue;
          }
          else
          {
            input[i][j]-=1;
          }
        }
    }
  }
  exchange_halo_1D_float(&input[0][0], temporal_size, d1, halo_size, pid, D1, mpi_size);
}

void test_double_1D(int mpi_rank, int mpi_size){
  const int halo_size=2;
  const int d1=6;
  const int temporal_size=1;
  const int d1h=d1+2*halo_size;
  double input[1][10]=
  {
    {0,0,1,2,3,4,5,6,0,0}
  };
  int pid=mpi_rank;
  int D1=mpi_size;
  exchange_halo_1D_double(&input[0][0], temporal_size, d1, halo_size, pid, D1, mpi_size);
  for(int i=0;i<temporal_size;i++){
    for(int j=0;j<d1h;j++){
      if(2*halo_size<=j&&j<halo_size+d1-halo_size)
        input[i][j]+=1;
    }
  }
  for(int i=0;i<temporal_size;i++){
    for(int j=0;j<d1h;j++){
        if(halo_size<=j&&j<halo_size+d1){
          if(2*halo_size<=j&&j<halo_size+d1-halo_size){
            continue;
          }
          else
          {
            input[i][j]-=1;
          }
        }
    }
  }
  exchange_halo_1D_double(&input[0][0], temporal_size, d1, halo_size, pid, D1, mpi_size);
}

void test_int_1D(int mpi_rank, int mpi_size){
  const int halo_size=2;
  const int d1=6;
  const int temporal_size=1;
  const int d1h=d1+2*halo_size;
  int input[1][10]=
  {
    {0,0,1,2,3,4,5,6,0,0}
  };
  int pid=mpi_rank;
  int D1=mpi_size;
  exchange_halo_1D_int(&input[0][0], temporal_size, d1, halo_size, pid, D1, mpi_size);
  for(int i=0;i<temporal_size;i++){
    for(int j=0;j<d1h;j++){
      if(2*halo_size<=j&&j<halo_size+d1-halo_size)
        input[i][j]+=1;
    }
  }
  for(int i=0;i<temporal_size;i++){
    for(int j=0;j<d1h;j++){
        if(halo_size<=j&&j<halo_size+d1){
          if(2*halo_size<=j&&j<halo_size+d1-halo_size){
            continue;
          }
          else
          {
            input[i][j]-=1;
          }
        }
    }
  }
  exchange_halo_1D_int(&input[0][0], temporal_size, d1, halo_size, pid, D1, mpi_size);
}

void test_float_1D_asyn(int mpi_rank, int mpi_size){
  const int halo_size=2;
  const int d1=6;
  const int temporal_size=1;
  const int d1h=d1+2*halo_size;
  float input[1][10]=
  {
    {0,0,1,2,3,4,5,6,0,0}
  };
  int pid=mpi_rank;
  int D1=mpi_size;
  exchange_halo_1D_float_start(&input[0][0], temporal_size, d1, halo_size, pid, D1, mpi_size);
  for(int i=0;i<temporal_size;i++){
    for(int j=0;j<d1h;j++){
      if(2*halo_size<=j&&j<halo_size+d1-halo_size)
        input[i][j]+=1;
    }
  }
  exchange_halo_1D_float_end(&input[0][0], temporal_size, d1, halo_size, pid, D1, mpi_size);
  for(int i=0;i<temporal_size;i++){
    for(int j=0;j<d1h;j++){
        if(halo_size<=j&&j<halo_size+d1){
          if(2*halo_size<=j&&j<halo_size+d1-halo_size){
            continue;
          }
          else
          {
            input[i][j]-=1;
          }
        }
    }
  }
  print_grid_1D_float_barrier(&input[0][0], d1, halo_size, mpi_rank);
}

void test_float_3D(int mpi_rank, int mpi_size){
  const int halo_size=1;
  const int d1=3, d2=3, d3=4;
  const int temporal_size=1;
  const int d1h=d1+2*halo_size;
  const int d2h=d2+2*halo_size;
  const int d3h=d3+2*halo_size;
  float input[1][5][5][6]=
  {
    {
      {
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0}
      },
      {
        {0, 0, 0, 0, 0, 0},
        {0, 1, 2, 3, 4, 0},
        {0, 5, 6, 7, 8, 0},
        {0, 9,10,11,12, 0},
        {0, 0, 0, 0, 0, 0}
      },
      {
        {0, 0, 0, 0, 0, 0},
        {0,13,14,15,16, 0},
        {0,17,18,19,20, 0},
        {0,21,22,23,24, 0},
        {0, 0, 0, 0, 0, 0}
      },
      {
        {0, 0, 0, 0, 0, 0},
        {0,25,26,27,28, 0},
        {0,29,30,31,32, 0},
        {0,33,34,35,36, 0},
        {0, 0, 0, 0, 0, 0}
      },
      {
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0}
      }
    }
  };
  int pid=mpi_rank;
  int D1=2;
  int D2=2;
  int D3=mpi_size/(D1*D2);
  exchange_halo_3D_float(&input[0][0][0][0], temporal_size, d1, d2 , d3, halo_size, pid, D1, D2, D3, mpi_size);
  for(int i=0;i<temporal_size;i++){
    for(int j=0;j<d1h;j++){
      for(int k=0;k<d2h;k++){
        for(int p=0;p<d3h;p++){
          if(halo_size<=j&&j<halo_size+d1 &&halo_size<=k&&k<halo_size+d2 &&halo_size<=p&&p<halo_size+d3)
            input[i][j][k][p]+=1;
        }
      }
    }
  }
  exchange_halo_3D_float(&input[0][0][0][0], temporal_size, d1, d2, d3, halo_size, pid, D1, D2, D3, mpi_size);
}

void test_double_3D(int mpi_rank, int mpi_size){
  const int halo_size=1;
  const int d1=3, d2=3, d3=4;
  const int temporal_size=1;
  const int d1h=d1+2*halo_size;
  const int d2h=d2+2*halo_size;
  const int d3h=d3+2*halo_size;
  double input[1][5][5][6]=
  {
    {
      {
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0}
      },
      {
        {0, 0, 0, 0, 0, 0},
        {0, 1, 2, 3, 4, 0},
        {0, 5, 6, 7, 8, 0},
        {0, 9,10,11,12, 0},
        {0, 0, 0, 0, 0, 0}
      },
      {
        {0, 0, 0, 0, 0, 0},
        {0,13,14,15,16, 0},
        {0,17,18,19,20, 0},
        {0,21,22,23,24, 0},
        {0, 0, 0, 0, 0, 0}
      },
      {
        {0, 0, 0, 0, 0, 0},
        {0,25,26,27,28, 0},
        {0,29,30,31,32, 0},
        {0,33,34,35,36, 0},
        {0, 0, 0, 0, 0, 0}
      },
      {
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0}
      }
    }
  };
  int pid=mpi_rank;
  int D1=2;
  int D2=2;
  int D3=mpi_size/(D1*D2);
  exchange_halo_3D_double(&input[0][0][0][0], temporal_size, d1, d2 , d3, halo_size, pid, D1, D2, D3, mpi_size);
  for(int i=0;i<temporal_size;i++){
    for(int j=0;j<d1h;j++){
      for(int k=0;k<d2h;k++){
        for(int p=0;p<d3h;p++){
          if(halo_size<=j&&j<halo_size+d1 &&halo_size<=k&&k<halo_size+d2 &&halo_size<=p&&p<halo_size+d3)
            input[i][j][k][p]+=1;
        }
      }
    }
  }
  exchange_halo_3D_double(&input[0][0][0][0], temporal_size, d1, d2, d3, halo_size, pid, D1, D2, D3, mpi_size);
}

void test_double_3D_asyn(int mpi_rank, int mpi_size){
  const int halo_size=1;
  const int d1=3, d2=3, d3=4;
  const int temporal_size=1;
  const int d1h=d1+2*halo_size;
  const int d2h=d2+2*halo_size;
  const int d3h=d3+2*halo_size;
  double input[1][5][5][6]=
  {
    {
      {
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0}
      },
      {
        {0, 0, 0, 0, 0, 0},
        {0, 1, 2, 3, 4, 0},
        {0, 5, 6, 7, 8, 0},
        {0, 9,10,11,12, 0},
        {0, 0, 0, 0, 0, 0}
      },
      {
        {0, 0, 0, 0, 0, 0},
        {0,13,14,15,16, 0},
        {0,17,18,19,20, 0},
        {0,21,22,23,24, 0},
        {0, 0, 0, 0, 0, 0}
      },
      {
        {0, 0, 0, 0, 0, 0},
        {0,25,26,27,28, 0},
        {0,29,30,31,32, 0},
        {0,33,34,35,36, 0},
        {0, 0, 0, 0, 0, 0}
      },
      {
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0},
        {0, 0, 0, 0, 0, 0}
      }
    }
  };
  int pid=mpi_rank;
  int D1=2;
  int D2=2;
  int D3=mpi_size/(D1*D2);
  exchange_halo_3D_double_start(&input[0][0][0][0], temporal_size, d1, d2 , d3, halo_size, pid, D1, D2, D3, mpi_size);
  for(int i=0;i<temporal_size;i++){
    for(int j=0;j<d1h;j++){
      for(int k=0;k<d2h;k++){
        for(int p=0;p<d3h;p++){
          if(2*halo_size<=j&&j<halo_size+d1-halo_size &&2*halo_size<=k&&k<halo_size+d2-halo_size &&2*halo_size<=p&&p<halo_size+d3-halo_size)
            input[i][j][k][p]+=1;
        }
      }
    }
  }
  exchange_halo_3D_double_end(&input[0][0][0][0], temporal_size, d1, d2, d3, halo_size, pid, D1, D2, D3, mpi_size);
  for(int i=0;i<temporal_size;i++){
    for(int j=0;j<d1h;j++){
      for(int k=0;k<d2h;k++){
        for(int p=0;p<d3h;p++){
          if(2*halo_size<=j&&j<halo_size+d1-halo_size &&2*halo_size<=k&&k<halo_size+d2-halo_size &&2*halo_size<=p&&p<halo_size+d3-halo_size){
            continue;
          }
          else if(halo_size<=j&&j<halo_size+d1 &&halo_size<=k&&k<halo_size+d2 &&halo_size<=p&&p<halo_size+d3){
            input[i][j][k][p]-=1;
          }
        }
      }
    }
  }
  /* print_grid*/   
  print_grid_3D_double_barrier(&input[0][0][0][0], d1, d2, d3, halo_size, pid);  

}





int main(){
  int rank;   //进程标识
  int size;   //进程总数

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

  MPI_Status status_s, status_r;  //status object (Status)
  MPI_Request handle_s, handle_r; //MPI request (handle)

  //test(rank, size);
  //test2(rank, size);
  //test2_double(rank, size);
  //test2_int(rank, size);
  //test_float_1D(rank, size);
  //test_double_1D(rank, size);
  //test_int_1D(rank, size);
  //test_float_1D_asyn(rank, size);
  //test_float_3D(rank, size);
  test_double_3D_asyn(rank, size);

  //
  MPI_Finalize(); /*MPI的结束函数*/
  return 0;
}
