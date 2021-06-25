#include <iostream>
#include <vector>
#include "../ir.h"
#include <cmath>
#include <string.h>

//config
const std::string input_file = "../data/rand.data";
//set datatype in here 
#define DT f32
//

void bench_2d9pt_star(std::string target, int num_threads, int D1, int D2,
        int M, int N, int tile_size_x,  int tile_size_y, int use_schedule=1, int is_assigned=0){
  init();
  const int halo_size = 2;
  const int time_window_size = 2;
  DefTensor2D_TimeWin(B, time_window_size, halo_size, DT, M, N)  ;
  DefVar(k,i32);
  DefVar(j,i32);
  DefVar_Value(c0, DT, 0.1);
  DefVar_Value(c1, DT, 0.2);
  DefVar_Value(c2, DT, 0.3);
  DefVar_Value(c3, DT, 0.4);
  DefVar_Value(c4, DT, 0.5);
  DefVar_Value(c5, DT, 0.4);
  DefVar_Value(c6, DT, 0.3);
  DefVar_Value(c7, DT, 0.2);
  DefVar_Value(c8, DT, 0.1);
  //
  Kernel A( (k,j) , c0* B[k,j] + c1* B[k,j-1] + c2* B[k,j-2] + c3* B[k,j+1] + c4* B[k,j+2] + c5* B[k-1,j] + c6* B[k-2,j] + c7* B[k+1,j] + c8* B[k+2,j] , schedule );
  //
  if (use_schedule){
    Axis xo, yo, xi, yi;
    A.tile(tile_size_x, tile_size_y, xo, xi, yo, yi);
    A.reorder(xo, yo, xi, yi);
    if (strcmp(target.c_str(), "sunway") == 0){
      CacheRead  buffer_read;
      CacheWrite buffer_write;
      A.cache_read  (B, buffer_read,  "global");
      A.cache_write (buffer_write, "global");
      A.compute_at(buffer_read,  yo);
      A.compute_at(buffer_write, yo);
      assert(num_threads==64);
    }
    A.parallel(xo, num_threads);
    A.build( target );
  }
  //
  Result Res( (k,j) , B[k,j] );
  //
  DefShapeMPI2D(shape_mpi, D1, D2);
  auto t = Stencil::t;
  Stencil stencil( shape_mpi, B, (k,j) , Res[t] << A[t-1], 0, is_assigned );

  stencil.input( shape_mpi, B, input_file);
  stencil.run(1,100);
  const int time_step10 = 100;
  //stencil.output(B, time_step10, "../data/out.data");
  //PrintIR();
  const bool use_timing = true;
  char func_name[100];
  sprintf(func_name, "2d9pt_star_%s_%d_%d_%d_%d_%d_%d_%d_%dschedule", target.c_str(), num_threads, D1,D2, M,N, tile_size_x, tile_size_y, use_schedule);
  stencil.compile_to_source_code_mpi(func_name, target, use_timing);
}


void bench_2d9pt_box(std::string target, int num_threads, int D1, int D2,
  int M, int N, int tile_size_x,  int tile_size_y, int use_schedule=1, int is_assigned=0){
  init();
  const int halo_size = 1;
  const int time_window_size = 2;
  DefTensor2D_TimeWin(B, time_window_size, halo_size, DT, M, N)  ;
  DefVar(k,i32);
  DefVar(j,i32);
  DefVar_Value(c0, DT, 0.1);
  DefVar_Value(c1, DT, 0.2);
  DefVar_Value(c2, DT, 0.3);
  DefVar_Value(c3, DT, 0.4);
  DefVar_Value(c4, DT, 0.5);
  DefVar_Value(c5, DT, 0.4);
  DefVar_Value(c6, DT, 0.3);
  DefVar_Value(c7, DT, 0.2);
  DefVar_Value(c8, DT, 0.1);
  //
  Kernel A( (k,j) , c0* B[k,j] + c1* B[k,j-1] + c2* B[k-1,j] + c3* B[k,j+1] + c4* B[k+1,j] + c5* B[k-1,j+1] + c6* B[k+1,j-1] + c7* B[k-1,j-1] + c8* B[k+1,j+1] , schedule );
  if (use_schedule){
    Axis xo, yo, xi, yi;
    A.tile(tile_size_x, tile_size_y, xo, xi, yo, yi);
    A.reorder(xo, yo, xi, yi);
    if (strcmp(target.c_str(), "sunway") == 0){
      CacheRead  buffer_read;
      CacheWrite buffer_write;
      A.cache_read  (B, buffer_read,  "global");
      A.cache_write (buffer_write, "global");
      A.compute_at(buffer_read,  yo);
      A.compute_at(buffer_write, yo);
      assert(num_threads==64);
    }
    A.parallel(xo, num_threads);
    A.build( target );
  }
  //
  Result Res( (k,j) , B[k,j] );
  //
  DefShapeMPI2D(shape_mpi, D1, D2);
  auto t = Stencil::t;
  Stencil stencil( shape_mpi, B, (k,j) , Res[t] << A[t-1], 0, is_assigned );

  stencil.input( shape_mpi, B, input_file);
  stencil.run(1,100);
  const int time_step10 = 100;
  //stencil.output(B, time_step10, "../data/out.data");
  //PrintIR();
  const bool use_timing = true;
  char func_name[100];
  sprintf(func_name, "2d9pt_box_%s_%d_%d_%d_%d_%d_%d_%d_%dschedule", target.c_str(), num_threads, D1,D2, M,N, tile_size_x, tile_size_y, use_schedule);
  stencil.compile_to_source_code_mpi(func_name, target, use_timing);
}


void bench_3d7pt_star(std::string target, int num_threads, int D1, int D2, int D3,
  int M, int N, int P, int tile_size_x,  int tile_size_y, int tile_size_z, int use_schedule=1, int is_assigned=0){
  init();
  const int halo_size = 1;
  const int time_window_size = 2;
  DefTensor3D_TimeWin(B, time_window_size, halo_size, DT, M, N, P)  ;
  DefVar(k,i32);
  DefVar(j,i32);
  DefVar(i,i32);
  DefVar_Value(c0, DT, 0.1);
  DefVar_Value(c1, DT, 0.2);
  DefVar_Value(c2, DT, 0.3);
  DefVar_Value(c3, DT, 0.4);
  DefVar_Value(c4, DT, 0.3);
  DefVar_Value(c5, DT, 0.2);
  DefVar_Value(c6, DT, 0.1);
  //
  Kernel A( (k,j,i) , c0* B[k,j,i] + c1* B[k,j,i-1] + c2* B[k,j,i+1] + c3* B[k-1,j,i] + c4* B[k+1,j,i] + c5* B[k,j-1,i] + c6* B[k,j+1,i] , schedule );
  if (use_schedule){
    Axis xo, yo, xi, yi, zo, zi;
    A.tile(tile_size_x, tile_size_y, tile_size_z, xo, xi, yo, yi, zo, zi);
    A.reorder(xo, yo, zo, xi, yi, zi);
    if (strcmp(target.c_str(), "sunway") == 0){
      CacheRead  buffer_read;
      CacheWrite buffer_write;
      A.cache_read  (B, buffer_read,  "global");
      A.cache_write (buffer_write, "global");
      A.compute_at(buffer_read,  zo);
      A.compute_at(buffer_write, zo);
      assert(num_threads==64);
    }
    A.parallel(xo, num_threads);
    A.build( target );
  }
  //
  Result Res( (k,j,i) , B[k,j,i] );
  //
  DefShapeMPI3D(shape_mpi, D1, D2, D3);
  auto t = Stencil::t;
  Stencil stencil( shape_mpi, B, (k,j,i) , Res[t] << A[t-1], 0, is_assigned );

  stencil.input( shape_mpi, B, input_file);
  stencil.run(1,100);
  const int time_step10 = 100;
  //stencil.output(B, time_step10, "../data/out.data");
  //PrintIR();
stencil.print_ir();
  const bool use_timing = true;
  char func_name[100];
  sprintf(func_name, "3d7pt_star_%s_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%dschedule", target.c_str(), num_threads, D1,D2,D3, M,N,P, tile_size_x, tile_size_y, tile_size_z, use_schedule);
  stencil.compile_to_source_code_mpi(func_name, target, use_timing);
}


void bench_3d13pt_star(std::string target, int num_threads, int D1, int D2, int D3, 
            int M, int N, int P, int tile_size_x,  int tile_size_y, int tile_size_z, int use_schedule=1, int is_assigned=0){
  init();
  const int halo_size = 2;
  const int time_window_size = 2;
  DefTensor3D_TimeWin(B, time_window_size, halo_size, DT, M, N, P)  ;
  DefVar(k,i32);
  DefVar(j,i32);
  DefVar(i,i32);
  DefVar_Value(a, DT, 0.1);
  DefVar_Value(d, DT, 0.2);
  DefVar_Value(e, DT, 0.3);
  DefVar_Value(f, DT, 0.4);
  //
  Kernel A( (k,j,i) , a* B[k,j,i] + d* (B[k-2,j,i] + B[k,j-2,i] + B[k,j,i-2] + B[k,j,i+2] + B[k,j+2,i] + B[k+2,j,i] ) + e* (B[k-1,j,i] + B[k,j-1,i] + B[k,j,i-1] + B[k,j,i+1] + B[k,j+1,i] + B[k+1,j,i] ) - f * B[k,j,i] , schedule );
  //
  if (use_schedule){
    Axis xo, yo, xi, yi, zo, zi;
    A.tile(tile_size_x, tile_size_y, tile_size_z, xo, xi, yo, yi, zo, zi);
    A.reorder(xo, yo, zo, xi, yi, zi);

    if (strcmp(target.c_str(), "sunway") == 0){
      CacheRead  buffer_read;
      CacheWrite buffer_write;
      A.cache_read  (B, buffer_read,  "global");
      A.cache_write (buffer_write, "global");
      A.compute_at(buffer_read,  zo);
      A.compute_at(buffer_write, zo);
      assert(num_threads==64);
    }
    A.parallel(xo, num_threads);
    A.build( target );
  }
  //
  Result Res( (k,j,i) , B[k,j,i] );
  //
  DefShapeMPI3D(shape_mpi, D1, D2, D3);
  auto t = Stencil::t;
  Stencil stencil( shape_mpi, B, (k,j,i) , Res[t] << A[t-1], 0, is_assigned );

  stencil.input( shape_mpi, B, input_file);
  stencil.run(1,100);
  const int time_step10 = 100;
  //stencil.output(B, time_step10, "../data/out.data");
  //PrintIR();
  const bool use_timing = true;
  char func_name[100];
  sprintf(func_name, "3d13pt_star_%s_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%dschedule", target.c_str(), num_threads, D1,D2,D3, M,N,P, tile_size_x, tile_size_y, tile_size_z, use_schedule);
  stencil.compile_to_source_code_mpi(func_name, target, use_timing);
}

void bench_2d121pt_box(std::string target, int num_threads, int D1, int D2,
            int M, int N, int tile_size_x,  int tile_size_y, int use_schedule=1, int is_assigned=0)
{
    init();
    const int halo_size = 5;
    const int time_window_size = 2;
    //const int M = 256; // x dim
    //const int N = 256; // y dim

    DefTensor2D_TimeWin(B, time_window_size, halo_size, DT, M, N);
    DefVar(x, i32);
    DefVar(y, i32);
    
    // 0
    DefVar_Value(c0_0, DT, 0); 
    DefVar_Value(c0_1, DT, 1);
    DefVar_Value(c0_2, DT, 2); 
    DefVar_Value(c0_3, DT, 3);
    DefVar_Value(c0_4, DT, 4); 
    DefVar_Value(c0_5, DT, 5);
    DefVar_Value(c0_6, DT, 6); 
    DefVar_Value(c0_7, DT, 7);
    DefVar_Value(c0_8, DT, 8); 
    DefVar_Value(c0_9, DT, 9);
    DefVar_Value(c0_10, DT, 10);
    // 1
    DefVar_Value(c1_0, DT, 0); 
    DefVar_Value(c1_1, DT, 1);
    DefVar_Value(c1_2, DT, 2); 
    DefVar_Value(c1_3, DT, 3);
    DefVar_Value(c1_4, DT, 4); 
    DefVar_Value(c1_5, DT, 5);
    DefVar_Value(c1_6, DT, 6); 
    DefVar_Value(c1_7, DT, 7);
    DefVar_Value(c1_8, DT, 8); 
    DefVar_Value(c1_9, DT, 9);
    DefVar_Value(c1_10, DT, 10);
    // 2
    DefVar_Value(c2_0, DT, 0); 
    DefVar_Value(c2_1, DT, 1);
    DefVar_Value(c2_2, DT, 2); 
    DefVar_Value(c2_3, DT, 3);
    DefVar_Value(c2_4, DT, 4); 
    DefVar_Value(c2_5, DT, 5);
    DefVar_Value(c2_6, DT, 6); 
    DefVar_Value(c2_7, DT, 7);
    DefVar_Value(c2_8, DT, 8); 
    DefVar_Value(c2_9, DT, 9);
    DefVar_Value(c2_10, DT, 10);
    // 3
    DefVar_Value(c3_0, DT, 0); 
    DefVar_Value(c3_1, DT, 1);
    DefVar_Value(c3_2, DT, 2); 
    DefVar_Value(c3_3, DT, 3);
    DefVar_Value(c3_4, DT, 4); 
    DefVar_Value(c3_5, DT, 5);
    DefVar_Value(c3_6, DT, 6); 
    DefVar_Value(c3_7, DT, 7);
    DefVar_Value(c3_8, DT, 8); 
    DefVar_Value(c3_9, DT, 9);
    DefVar_Value(c3_10, DT, 10);
    // 4
    DefVar_Value(c4_0, DT, 0); 
    DefVar_Value(c4_1, DT, 1);
    DefVar_Value(c4_2, DT, 2); 
    DefVar_Value(c4_3, DT, 3);
    DefVar_Value(c4_4, DT, 4); 
    DefVar_Value(c4_5, DT, 5);
    DefVar_Value(c4_6, DT, 6); 
    DefVar_Value(c4_7, DT, 7);
    DefVar_Value(c4_8, DT, 8); 
    DefVar_Value(c4_9, DT, 9);
    DefVar_Value(c4_10, DT, 10);
    // 5
    DefVar_Value(c5_0, DT, 0); 
    DefVar_Value(c5_1, DT, 1);
    DefVar_Value(c5_2, DT, 2); 
    DefVar_Value(c5_3, DT, 3);
    DefVar_Value(c5_4, DT, 4); 
    DefVar_Value(c5_5, DT, 5);
    DefVar_Value(c5_6, DT, 6); 
    DefVar_Value(c5_7, DT, 7);
    DefVar_Value(c5_8, DT, 8); 
    DefVar_Value(c5_9, DT, 9);
    DefVar_Value(c5_10, DT, 10);
    // 6
    DefVar_Value(c6_0, DT, 0); 
    DefVar_Value(c6_1, DT, 1);
    DefVar_Value(c6_2, DT, 2); 
    DefVar_Value(c6_3, DT, 3);
    DefVar_Value(c6_4, DT, 4); 
    DefVar_Value(c6_5, DT, 5);
    DefVar_Value(c6_6, DT, 6); 
    DefVar_Value(c6_7, DT, 7);
    DefVar_Value(c6_8, DT, 8); 
    DefVar_Value(c6_9, DT, 9);
    DefVar_Value(c6_10, DT, 10);
    // 7
    DefVar_Value(c7_0, DT, 0); 
    DefVar_Value(c7_1, DT, 1);
    DefVar_Value(c7_2, DT, 2); 
    DefVar_Value(c7_3, DT, 3);
    DefVar_Value(c7_4, DT, 4); 
    DefVar_Value(c7_5, DT, 5);
    DefVar_Value(c7_6, DT, 6); 
    DefVar_Value(c7_7, DT, 7);
    DefVar_Value(c7_8, DT, 8); 
    DefVar_Value(c7_9, DT, 9);
    DefVar_Value(c7_10, DT, 10);
    // 8
    DefVar_Value(c8_0, DT, 0); 
    DefVar_Value(c8_1, DT, 1);
    DefVar_Value(c8_2, DT, 2); 
    DefVar_Value(c8_3, DT, 3);
    DefVar_Value(c8_4, DT, 4); 
    DefVar_Value(c8_5, DT, 5);
    DefVar_Value(c8_6, DT, 6); 
    DefVar_Value(c8_7, DT, 7);
    DefVar_Value(c8_8, DT, 8); 
    DefVar_Value(c8_9, DT, 9);
    DefVar_Value(c8_10, DT, 10);
    // 9
    DefVar_Value(c9_0, DT, 0); 
    DefVar_Value(c9_1, DT, 1);
    DefVar_Value(c9_2, DT, 2); 
    DefVar_Value(c9_3, DT, 3);
    DefVar_Value(c9_4, DT, 4); 
    DefVar_Value(c9_5, DT, 5);
    DefVar_Value(c9_6, DT, 6); 
    DefVar_Value(c9_7, DT, 7);
    DefVar_Value(c9_8, DT, 8); 
    DefVar_Value(c9_9, DT, 9);
    DefVar_Value(c9_10, DT, 10);
    // 10
    DefVar_Value(c10_0, DT, 0); 
    DefVar_Value(c10_1, DT, 1);
    DefVar_Value(c10_2, DT, 2); 
    DefVar_Value(c10_3, DT, 3);
    DefVar_Value(c10_4, DT, 4); 
    DefVar_Value(c10_5, DT, 5);
    DefVar_Value(c10_6, DT, 6); 
    DefVar_Value(c10_7, DT, 7);
    DefVar_Value(c10_8, DT, 8); 
    DefVar_Value(c10_9, DT, 9);
    DefVar_Value(c10_10, DT, 10);

    Kernel A((x, y),
            0.1*c0_0*B[x-5, y-5] + c0_1*B[x-4, y-5] +c0_2*B[x-3, y-5] + c0_3*B[x-2, y-5] - c0_4*B[x-1, y-5] + c0_5*B[x, y-5] - c0_6*B[x+1, y-5] + c0_7*B[x+2, y-5] - c0_8*B[x+3, y-5] + c0_9*B[x+4, y-5] + c0_10*B[x+5, y-5]
            -   c1_0*B[x-5, y-4] + c1_1*B[x-4, y-4] +c1_2*B[x-3, y-4] + c1_3*B[x-2, y-4] - c1_4*B[x-1, y-4] + c1_5*B[x, y-4] - c1_6*B[x+1, y-4] + c1_7*B[x+2, y-4] - c1_8*B[x+3, y-4] + c1_9*B[x+4, y-4] + c1_10*B[x+5, y-4]
            -   c2_0*B[x-5, y-3] + c2_1*B[x-4, y-3] +c2_2*B[x-3, y-3] + c2_3*B[x-2, y-3] - c2_4*B[x-1, y-3] + c2_5*B[x, y-3] - c2_6*B[x+1, y-3] + c2_7*B[x+2, y-3] - c2_8*B[x+3, y-3] + c2_9*B[x+4, y-3] + c2_10*B[x+5, y-3]
            -   c3_0*B[x-5, y-2] + c3_1*B[x-4, y-2] +c3_2*B[x-3, y-2] + c3_3*B[x-2, y-2] - c3_4*B[x-1, y-2] + c3_5*B[x, y-2] - c3_6*B[x+1, y-2] + c3_7*B[x+2, y-2] - c3_8*B[x+3, y-2] + c3_9*B[x+4, y-2] + c3_10*B[x+5, y-2]
            -   c4_0*B[x-5, y-1] + c4_1*B[x-4, y-1] +c4_2*B[x-3, y-1] + c4_3*B[x-2, y-1] - c4_4*B[x-1, y-1] + c4_5*B[x, y-1] - c4_6*B[x+1, y-1] + c4_7*B[x+2, y-1] - c4_8*B[x+3, y-1] + c4_9*B[x+4, y-1] + c4_10*B[x+5, y-1]
            -   c5_0*B[x-5, y] + c5_1*B[x-4, y] +c5_2*B[x-3, y] + c5_3*B[x-2, y] - c5_4*B[x-1, y] + c5_5*B[x, y] - c5_6*B[x+1, y] + c5_7*B[x+2, y] - c5_8*B[x+3, y] + c5_9*B[x+4, y] + c5_10*B[x+5, y]
            -   c6_0*B[x-5, y+1] + c6_1*B[x-4, y+1] +c6_2*B[x-3, y+1] + c6_3*B[x-2, y+1] - c6_4*B[x-1, y+1] + c6_5*B[x, y+1] - c6_6*B[x+1, y+1] + c6_7*B[x+2, y+1] - c6_8*B[x+3, y+1] + c6_9*B[x+4, y+1] + c6_10*B[x+5, y+1]
            -   c7_0*B[x-5, y+2] + c7_1*B[x-4, y+2] +c7_2*B[x-3, y+2] + c7_3*B[x-2, y+2] - c7_4*B[x-1, y+2] + c7_5*B[x, y+2] - c7_6*B[x+1, y+2] + c7_7*B[x+2, y+2] - c7_8*B[x+3, y+2] + c7_9*B[x+4, y+2] + c7_10*B[x+5, y+2]
            -   c8_0*B[x-5, y+3] + c8_1*B[x-4, y+3] +c8_2*B[x-3, y+3] + c8_3*B[x-2, y+3] - c8_4*B[x-1, y+3] + c8_5*B[x, y+3] - c8_6*B[x+1, y+3] + c8_7*B[x+2, y+3] - c8_8*B[x+3, y+3] + c8_9*B[x+4, y+3] + c8_10*B[x+5, y+3]
            -   c9_0*B[x-5, y+4] + c9_1*B[x-4, y+4] +c9_2*B[x-3, y+4] + c9_3*B[x-2, y+4] - c9_4*B[x-1, y+4] + c9_5*B[x, y+4] - c9_6*B[x+1, y+4] + c9_7*B[x+2, y+4] - c9_8*B[x+3, y+4] + c9_9*B[x+4, y+4] + c9_10*B[x+5, y+4]
            -   c10_0*B[x-5, y+5] + c10_1*B[x-4, y+5] +c10_2*B[x-3, y+5] + c10_3*B[x-2, y+5] - c10_4*B[x-1, y+5] + c10_5*B[x, y+5] - c10_6*B[x+1, y+5] + c10_7*B[x+2, y+5] - c10_8*B[x+3, y+5] + c10_9*B[x+4, y+5] + c10_10*B[x+5, y+5], schedule);

    //const int tile_size_x = 2;
    //const int tile_size_y = 8;

    if (use_schedule){
      Axis xo, yo, xi, yi;
      A.tile(tile_size_x, tile_size_y, xo, xi, yo, yi);
      //A.reorder(xo, xi, yo, yi);  // 外层
       A.reorder(xo, yo, xi, yi);  // 内层

      if (strcmp(target.c_str(), "sunway") == 0){
        CacheRead buffer_read;
        CacheWrite buffer_write;

        A.cache_read(B, buffer_read, "global");
        A.cache_write(buffer_write, "global");
        //A.compute_at(buffer_read, xo);
        //A.compute_at(buffer_write, xo);  // 外层
        A.compute_at(buffer_read, yo);
        A.compute_at(buffer_write, yo);  // 内层

        assert(num_threads==64);
      }
      A.parallel(xo, num_threads);
      A.build(target);
    }

    Result Res((x, y), B[x, y]);
    DefShapeMPI2D(shape_mpi, D1, D2);
    auto t = Stencil::t;
    Stencil stencil(shape_mpi, B, (x, y), Res[t] << A[t-1], 0, is_assigned);

    stencil.input(shape_mpi, B, input_file);
    stencil.run(1, 100);

    const bool use_timing = true;
    char func_name[100];
    sprintf(func_name, "2d121pt_box_%s_%d_%d_%d_%d_%d_%d_%d_%dschedule", target.c_str(), num_threads, D1,D2, M,N, tile_size_x, tile_size_y, use_schedule);
    stencil.compile_to_source_code_mpi(func_name, target, use_timing);
}

void bench_2d169pt_box(std::string target, int num_threads, int D1, int D2,
            int M, int N, int tile_size_x,  int tile_size_y, int use_schedule=1, int is_assigned=0)
{
    init();
    const int halo_size = 6;
    const int time_window_size = 2;
    //const int M = 1024; // x dim
    //const int N = 1024; // y dim

    DefTensor2D_TimeWin(B, time_window_size, halo_size, DT, M, N);
    DefVar(x, i32);
    DefVar(y, i32);
    
    // 0
    DefVar_Value(c0_0, DT, 0); 
    DefVar_Value(c0_1, DT, 1);
    DefVar_Value(c0_2, DT, 2); 
    DefVar_Value(c0_3, DT, 3);
    DefVar_Value(c0_4, DT, 4); 
    DefVar_Value(c0_5, DT, 5);
    DefVar_Value(c0_6, DT, 6); 
    DefVar_Value(c0_7, DT, 7);
    DefVar_Value(c0_8, DT, 8); 
    DefVar_Value(c0_9, DT, 9);
    DefVar_Value(c0_10, DT, 10);
    DefVar_Value(c0_11, DT, 11);
    DefVar_Value(c0_12, DT, 12);
    // 1
    DefVar_Value(c1_0, DT, 0); 
    DefVar_Value(c1_1, DT, 1);
    DefVar_Value(c1_2, DT, 2); 
    DefVar_Value(c1_3, DT, 3);
    DefVar_Value(c1_4, DT, 4); 
    DefVar_Value(c1_5, DT, 5);
    DefVar_Value(c1_6, DT, 6); 
    DefVar_Value(c1_7, DT, 7);
    DefVar_Value(c1_8, DT, 8); 
    DefVar_Value(c1_9, DT, 9);
    DefVar_Value(c1_10, DT, 10);
    DefVar_Value(c1_11, DT, 11);
    DefVar_Value(c1_12, DT, 12);
    // 2
    DefVar_Value(c2_0, DT, 0); 
    DefVar_Value(c2_1, DT, 1);
    DefVar_Value(c2_2, DT, 2); 
    DefVar_Value(c2_3, DT, 3);
    DefVar_Value(c2_4, DT, 4); 
    DefVar_Value(c2_5, DT, 5);
    DefVar_Value(c2_6, DT, 6); 
    DefVar_Value(c2_7, DT, 7);
    DefVar_Value(c2_8, DT, 8); 
    DefVar_Value(c2_9, DT, 9);
    DefVar_Value(c2_10, DT, 10);    
    DefVar_Value(c2_11, DT, 11);
    DefVar_Value(c2_12, DT, 12);

    // 3
    DefVar_Value(c3_0, DT, 0); 
    DefVar_Value(c3_1, DT, 1);
    DefVar_Value(c3_2, DT, 2); 
    DefVar_Value(c3_3, DT, 3);
    DefVar_Value(c3_4, DT, 4); 
    DefVar_Value(c3_5, DT, 5);
    DefVar_Value(c3_6, DT, 6); 
    DefVar_Value(c3_7, DT, 7);
    DefVar_Value(c3_8, DT, 8); 
    DefVar_Value(c3_9, DT, 9);
    DefVar_Value(c3_10, DT, 10);
    DefVar_Value(c3_11, DT, 11);
    DefVar_Value(c3_12, DT, 12);
    // 4
    DefVar_Value(c4_0, DT, 0); 
    DefVar_Value(c4_1, DT, 1);
    DefVar_Value(c4_2, DT, 2); 
    DefVar_Value(c4_3, DT, 3);
    DefVar_Value(c4_4, DT, 4); 
    DefVar_Value(c4_5, DT, 5);
    DefVar_Value(c4_6, DT, 6); 
    DefVar_Value(c4_7, DT, 7);
    DefVar_Value(c4_8, DT, 8); 
    DefVar_Value(c4_9, DT, 9);
    DefVar_Value(c4_10, DT, 10);
    DefVar_Value(c4_11, DT, 11);
    DefVar_Value(c4_12, DT, 12);
    // 5
    DefVar_Value(c5_0, DT, 0); 
    DefVar_Value(c5_1, DT, 1);
    DefVar_Value(c5_2, DT, 2); 
    DefVar_Value(c5_3, DT, 3);
    DefVar_Value(c5_4, DT, 4); 
    DefVar_Value(c5_5, DT, 5);
    DefVar_Value(c5_6, DT, 6); 
    DefVar_Value(c5_7, DT, 7);
    DefVar_Value(c5_8, DT, 8); 
    DefVar_Value(c5_9, DT, 9);
    DefVar_Value(c5_10, DT, 10);
    DefVar_Value(c5_11, DT, 11);
    DefVar_Value(c5_12, DT, 12);
    // 6
    DefVar_Value(c6_0, DT, 0); 
    DefVar_Value(c6_1, DT, 1);
    DefVar_Value(c6_2, DT, 2); 
    DefVar_Value(c6_3, DT, 3);
    DefVar_Value(c6_4, DT, 4); 
    DefVar_Value(c6_5, DT, 5);
    DefVar_Value(c6_6, DT, 6); 
    DefVar_Value(c6_7, DT, 7);
    DefVar_Value(c6_8, DT, 8); 
    DefVar_Value(c6_9, DT, 9);
    DefVar_Value(c6_10, DT, 10);
    DefVar_Value(c6_11, DT, 11);
    DefVar_Value(c6_12, DT, 12);
    // 7
    DefVar_Value(c7_0, DT, 0); 
    DefVar_Value(c7_1, DT, 1);
    DefVar_Value(c7_2, DT, 2); 
    DefVar_Value(c7_3, DT, 3);
    DefVar_Value(c7_4, DT, 4); 
    DefVar_Value(c7_5, DT, 5);
    DefVar_Value(c7_6, DT, 6); 
    DefVar_Value(c7_7, DT, 7);
    DefVar_Value(c7_8, DT, 8); 
    DefVar_Value(c7_9, DT, 9);
    DefVar_Value(c7_10, DT, 10);
    DefVar_Value(c7_11, DT, 11);
    DefVar_Value(c7_12, DT, 12);
    // 8
    DefVar_Value(c8_0, DT, 0); 
    DefVar_Value(c8_1, DT, 1);
    DefVar_Value(c8_2, DT, 2); 
    DefVar_Value(c8_3, DT, 3);
    DefVar_Value(c8_4, DT, 4); 
    DefVar_Value(c8_5, DT, 5);
    DefVar_Value(c8_6, DT, 6); 
    DefVar_Value(c8_7, DT, 7);
    DefVar_Value(c8_8, DT, 8); 
    DefVar_Value(c8_9, DT, 9);
    DefVar_Value(c8_10, DT, 10);
    DefVar_Value(c8_11, DT, 11);
    DefVar_Value(c8_12, DT, 12);
    // 9
    DefVar_Value(c9_0, DT, 0); 
    DefVar_Value(c9_1, DT, 1);
    DefVar_Value(c9_2, DT, 2); 
    DefVar_Value(c9_3, DT, 3);
    DefVar_Value(c9_4, DT, 4); 
    DefVar_Value(c9_5, DT, 5);
    DefVar_Value(c9_6, DT, 6); 
    DefVar_Value(c9_7, DT, 7);
    DefVar_Value(c9_8, DT, 8); 
    DefVar_Value(c9_9, DT, 9);
    DefVar_Value(c9_10, DT, 10);
    DefVar_Value(c9_11, DT, 11);
    DefVar_Value(c9_12, DT, 12);
    // 10
    DefVar_Value(c10_0, DT, 0); 
    DefVar_Value(c10_1, DT, 1);
    DefVar_Value(c10_2, DT, 2); 
    DefVar_Value(c10_3, DT, 3);
    DefVar_Value(c10_4, DT, 4); 
    DefVar_Value(c10_5, DT, 5);
    DefVar_Value(c10_6, DT, 6); 
    DefVar_Value(c10_7, DT, 7);
    DefVar_Value(c10_8, DT, 8); 
    DefVar_Value(c10_9, DT, 9);
    DefVar_Value(c10_10, DT, 10);
    DefVar_Value(c10_11, DT, 11);
    DefVar_Value(c10_12, DT, 12);
    // 11
    DefVar_Value(c11_0, DT, 0); 
    DefVar_Value(c11_1, DT, 1);
    DefVar_Value(c11_2, DT, 2); 
    DefVar_Value(c11_3, DT, 3);
    DefVar_Value(c11_4, DT, 4); 
    DefVar_Value(c11_5, DT, 5);
    DefVar_Value(c11_6, DT, 6); 
    DefVar_Value(c11_7, DT, 7);
    DefVar_Value(c11_8, DT, 8); 
    DefVar_Value(c11_9, DT, 9);
    DefVar_Value(c11_10, DT, 10);
    DefVar_Value(c11_11, DT, 11);
    DefVar_Value(c11_12, DT, 12);
    // 12
    DefVar_Value(c12_0, DT, 0); 
    DefVar_Value(c12_1, DT, 1);
    DefVar_Value(c12_2, DT, 2); 
    DefVar_Value(c12_3, DT, 3);
    DefVar_Value(c12_4, DT, 4); 
    DefVar_Value(c12_5, DT, 5);
    DefVar_Value(c12_6, DT, 6); 
    DefVar_Value(c12_7, DT, 7);
    DefVar_Value(c12_8, DT, 8); 
    DefVar_Value(c12_9, DT, 9);
    DefVar_Value(c12_10, DT, 10);
    DefVar_Value(c12_11, DT, 11);
    DefVar_Value(c12_12, DT, 12);

    Kernel A((x, y),
            0.1*c0_0*B[x-6, y-6] + c0_1*B[x-5, y-6] +c0_2*B[x-4, y-6] + c0_3*B[x-3, y-6] - c0_4*B[x-2, y-6] + c0_5*B[x-1, y-6] - c0_6*B[x, y-6] + c0_7*B[x+1, y-6] - c0_8*B[x+2, y-6] + c0_9*B[x+3, y-6] + c0_10*B[x+4, y-6] + c0_11*B[x+5, y-6] + c0_12*B[x+6, y-6]
            -   c1_0*B[x-6, y-5] + c1_1*B[x-5, y-5] +c1_2*B[x-4, y-5] + c1_3*B[x-3, y-5] - c1_4*B[x-2, y-5] + c1_5*B[x-1, y-5] - c1_6*B[x, y-5] + c1_7*B[x+1, y-5] - c1_8*B[x+2, y-5] + c1_9*B[x+3, y-5] + c1_10*B[x+4, y-5] + c1_11*B[x+5, y-5] + c1_12*B[x+6, y-5]
            -   c2_0*B[x-6, y-4] + c2_1*B[x-5, y-4] +c2_2*B[x-4, y-4] + c2_3*B[x-3, y-4] - c2_4*B[x-2, y-4] + c2_5*B[x-1, y-4] - c2_6*B[x, y-4] + c2_7*B[x+1, y-4] - c2_8*B[x+2, y-4] + c2_9*B[x+3, y-4] + c2_10*B[x+4, y-4] + c2_11*B[x+5, y-4] + c2_12*B[x+6, y-4]
            -   c3_0*B[x-6, y-3] + c3_1*B[x-5, y-3] +c3_2*B[x-4, y-3] + c3_3*B[x-3, y-3] - c3_4*B[x-2, y-3] + c3_5*B[x-1, y-3] - c3_6*B[x, y-3] + c3_7*B[x+1, y-3] - c3_8*B[x+2, y-3] + c3_9*B[x+3, y-3] + c3_10*B[x+4, y-3] + c3_11*B[x+5, y-3] + c3_12*B[x+6, y-3]
            -   c4_0*B[x-6, y-2] + c4_1*B[x-5, y-2] +c4_2*B[x-4, y-2] + c4_3*B[x-3, y-2] - c4_4*B[x-2, y-2] + c4_5*B[x-1, y-2] - c4_6*B[x, y-2] + c4_7*B[x+1, y-2] - c4_8*B[x+2, y-2] + c4_9*B[x+3, y-2] + c4_10*B[x+4, y-2] + c4_11*B[x+5, y-2] + c4_12*B[x+6, y-2]
            -   c5_0*B[x-6, y-1] + c5_1*B[x-5, y-1] +c5_2*B[x-4, y-1] + c5_3*B[x-3, y-1] - c5_4*B[x-2, y-1] + c5_5*B[x-1, y-1] - c5_6*B[x, y-1] + c5_7*B[x+1, y-1] - c5_8*B[x+2, y-1] + c5_9*B[x+3, y-1] + c5_10*B[x+4, y-1] + c5_11*B[x+5, y-1] + c5_12*B[x+6, y-1]
            -   c6_0*B[x-6, y] + c6_1*B[x-5, y] +c6_2*B[x-4, y] + c6_3*B[x-3, y] - c6_4*B[x-2, y] + c6_5*B[x-1, y] - c6_6*B[x, y] + c6_7*B[x+1, y] - c6_8*B[x+2, y] + c6_9*B[x+3, y] + c6_10*B[x+4, y] + c6_11*B[x+5, y] + c6_12*B[x+6, y]
            -   c7_0*B[x-6, y+1] + c7_1*B[x-5, y+1] +c7_2*B[x-4, y+1] + c7_3*B[x-3, y+1] - c7_4*B[x-2, y+1] + c7_5*B[x-1, y+1] - c7_6*B[x, y+1] + c7_7*B[x+1, y+1] - c7_8*B[x+2, y+1] + c7_9*B[x+3, y+1] + c7_10*B[x+4, y+1] + c7_11*B[x+5, y+1] + c7_12*B[x+6, y+1]
            -   c8_0*B[x-6, y+2] + c8_1*B[x-5, y+2] +c8_2*B[x-4, y+2] + c8_3*B[x-3, y+2] - c8_4*B[x-2, y+2] + c8_5*B[x-1, y+2] - c8_6*B[x, y+2] + c8_7*B[x+1, y+2] - c8_8*B[x+2, y+2] + c8_9*B[x+3, y+2] + c8_10*B[x+4, y+2] + c8_11*B[x+5, y+2] + c8_12*B[x+6, y+2]
            -   c9_0*B[x-6, y+3] + c9_1*B[x-5, y+3] +c9_2*B[x-4, y+3] + c9_3*B[x-3, y+3] - c9_4*B[x-2, y+3] + c9_5*B[x-1, y+3] - c9_6*B[x, y+3] + c9_7*B[x+1, y+3] - c9_8*B[x+2, y+3] + c9_9*B[x+3, y+3] + c9_10*B[x+4, y+3] + c9_11*B[x+5, y+3] + c9_12*B[x+6, y+3]
            -   c10_0*B[x-6, y+4] + c10_1*B[x-5, y+4] +c10_2*B[x-4, y+4] + c10_3*B[x-3, y+4] - c10_4*B[x-2, y+4] + c10_5*B[x-1, y+4] - c10_6*B[x, y+4] + c10_7*B[x+1, y+4] - c10_8*B[x+2, y+4] + c10_9*B[x+3, y+4] + c10_10*B[x+4, y+4] + c10_11*B[x+5, y+4] + c10_12*B[x+6, y+4]
            -   c11_0*B[x-6, y-6] + c11_1*B[x-5, y+5] +c11_2*B[x-4, y+5] + c11_3*B[x-3, y+5] - c11_4*B[x-2, y+5] + c11_5*B[x-1, y+5] - c11_6*B[x, y+5] + c11_7*B[x+1, y+5] - c11_8*B[x+2, y+5] + c11_9*B[x+3, y+5] + c11_10*B[x+4, y+5] + c11_11*B[x+5, y+5] + c11_12*B[x+6, y+5]
            -   c12_0*B[x-6, y+6] + c12_1*B[x-5, y+6] +c12_2*B[x-4, y+6] + c12_3*B[x-3, y+6] - c12_4*B[x-2, y+6] + c12_5*B[x-1, y+6] - c12_6*B[x, y+6] + c12_7*B[x+1, y+6] - c12_8*B[x+2, y+6] + c12_9*B[x+3, y+6] + c12_10*B[x+4, y+6] + c12_11*B[x+5, y+6] + c12_12*B[x+6, y+6], schedule);

    //const int tile_size_x = 32;
    //const int tile_size_y = 32;

    if (use_schedule){
      Axis xo, yo, xi, yi;
      A.tile(tile_size_x, tile_size_y, xo, xi, yo, yi);
      A.reorder(xo, yo, xi, yi);

      if (strcmp(target.c_str(), "sunway") == 0){
        CacheRead buffer_read;
        CacheWrite buffer_write;

        A.cache_read(B, buffer_read, "global");
        A.cache_write(buffer_write, "global");
        A.compute_at(buffer_read, yo);
        A.compute_at(buffer_write, yo);

        assert(num_threads == 64);
      }
      A.parallel(xo, num_threads);
      A.build(target);
    }

    Result Res((x, y), B[x, y]);
    DefShapeMPI2D(shape_mpi, D1, D2);
    auto t = Stencil::t;
    Stencil stencil(shape_mpi, B, (x, y), Res[t] << A[t-1], 0, is_assigned);

    stencil.input(shape_mpi, B, input_file);
    stencil.run(1, 100);

    const bool use_timing = true;
    char func_name[100];
    sprintf(func_name, "2d169pt_box_%s_%d_%d_%d_%d_%d_%d_%d_%dschedule", target.c_str(), num_threads, D1,D2, M,N, tile_size_x, tile_size_y, use_schedule);
    stencil.compile_to_source_code_mpi(func_name, target, use_timing);
}

void bench_3d25pt_star(std::string target, int num_threads, int D1, int D2, int D3, 
            int M, int N, int P, int tile_size_x,  int tile_size_y, int tile_size_z, int use_schedule=1, int is_assigned=0)
{
    init();

    const int halo_size = 4;
    const int time_window_size = 2;
    //const int M = 1024;
    //const int N = 16;
    //const int P = 16;

    DefTensor3D_TimeWin(B, time_window_size, halo_size, DT, M, N, P);
    DefVar(x, i32);
    DefVar(y, i32);
    DefVar(z, i32);
    DefVar_Value(dxinv1, DT, 0.1);
    DefVar_Value(dxinv2, DT, 0.2);
    DefVar_Value(dxinv3, DT, 0.3);
    DefVar_Value(a, DT, -1);
    DefVar_Value(b, DT, 0.8);
    DefVar_Value(c, DT, 0.2);
    DefVar_Value(d, DT, 0.038);
    DefVar_Value(e, DT, 0.0035);

    Kernel A((x, y, z), 
                a*(b*(B[x+1, y, z] - B[x-1, y, z])
                -c*(B[x+2, y, z] - B[x-2, y, z])
                +d*(B[x+3, y, z] - B[x-3, y, z])
                -e*(B[x+4, y, z] - B[x-4, y, z]))*dxinv1

                -(b*(B[x, y+1, z] - B[x, y-1, z])
                -c*(B[x, y+2, z] - B[x, y-2, z])
                +d*(B[x, y+3, z] - B[x, y-3, z])
                -e*(B[x, y+4, z] - B[x, y-4, z]))*dxinv2

                -(b*(B[x, y, z+1] - B[x, y, z-1])
                -c*(B[x, y, z+2] - B[x, y, z-2])
                +d*(B[x, y, z+3] - B[x, y, z-3])
                -e*(B[x, y, z+4] - B[x, y, z-4]))*dxinv3, schedule);

    //const int tile_size_x = 2;
    //const int tile_size_y = 4;
    //const int tile_size_z = 8;
    if (use_schedule){
      Axis xo, xi, yo, yi, zo, zi;
      A.tile(tile_size_x, tile_size_y, tile_size_z, xo, xi, yo, yi, zo, zi);
      // A.reorder(xo, xi, yo, zo, yi, zi); // 最外层
      // A.reorder(xo, yo, zo, xi, yi, zi);  // 中层
      A.reorder(xo, yo, zo, xi, yi, zi);  // 内层

      if (strcmp(target.c_str(), "sunway") == 0){
        CacheRead buffer_read;
        CacheWrite buffer_write;
        A.cache_read(B, buffer_read, "global");
        A.cache_write(buffer_write, "global");
        // A.compute_at(buffer_read, xo);
        // A.compute_at(buffer_write, xo);  // 最外层
        // A.compute_at(buffer_read, yo);
        // A.compute_at(buffer_write, yo);     // 中层
        A.compute_at(buffer_read, zo);
        A.compute_at(buffer_write, zo);     // 内层

        assert(num_threads == 64);
      }
      A.parallel(xo, num_threads);
      A.build(target);
    }

    Result Res((x, y, z), B[x, y, z]);
    DefShapeMPI3D(shape_mpi, D1, D2, D3);
    auto t = Stencil::t;
    Stencil stencil(shape_mpi, B, (x, y, z), Res[t] << A[t-1], 0, is_assigned);

    stencil.input(shape_mpi, B, input_file);
    stencil.run(1, 100);
    const bool use_timing = true;
    char func_name[100];
    sprintf(func_name, "3d25pt_star_%s_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%dschedule", target.c_str(), num_threads, D1,D2,D3, M,N,P, tile_size_x, tile_size_y, tile_size_z, use_schedule);
    stencil.compile_to_source_code_mpi(func_name, target, use_timing);
}

void bench_3d31pt_star(std::string target, int num_threads, int D1, int D2, int D3, 
            int M, int N, int P, int tile_size_x,  int tile_size_y, int tile_size_z, int use_schedule=1, int is_assigned=0)
{
    init();

    const int halo_size = 5;
    const int time_window_size = 2;
    //const int M = 1024;
    //const int N = 32;
    //const int P = 32;

    DefTensor3D_TimeWin(B, time_window_size, halo_size, DT, M, N, P);
    DefVar(x, i32);
    DefVar(y, i32);
    DefVar(z, i32);
    DefVar_Value(dxinv1, DT, 0.1);
    DefVar_Value(dxinv2, DT, 0.2);
    DefVar_Value(dxinv3, DT, 0.3);
    DefVar_Value(a, DT, -1);
    DefVar_Value(b, DT, 0.8);
    DefVar_Value(c, DT, 0.2);
    DefVar_Value(d, DT, 0.038);
    DefVar_Value(e, DT, 0.0035);
    DefVar_Value(f, DT, 0.00028);

    Kernel A((x, y, z), 
                a*(b*(B[x+1, y, z] - B[x-1, y, z])
                -c*(B[x+2, y, z] - B[x-2, y, z])
                +d*(B[x+3, y, z] - B[x-3, y, z])
                -e*(B[x+4, y, z] - B[x-4, y, z])
                -f*(B[x+5, y, z] - B[x-5, y, z]))*dxinv1

                -(b*(B[x, y+1, z] - B[x, y-1, z])
                -c*(B[x, y+2, z] - B[x, y-2, z])
                +d*(B[x, y+3, z] - B[x, y-3, z])
                -e*(B[x, y+4, z] - B[x, y-4, z])
                -f*(B[x, y+5, z] - B[x, y-5, z]))*dxinv2

                -(b*(B[x, y, z+1] - B[x, y, z-1])
                -c*(B[x, y, z+2] - B[x, y, z-2])
                +d*(B[x, y, z+3] - B[x, y, z-3])
                -e*(B[x, y, z+4] - B[x, y, z-4])
                -f*(B[x, y, z+5] - B[x, y, z-5]))*dxinv3, schedule);

    //const int tile_size_x = 1;
    //const int tile_size_y = 8;
    //const int tile_size_z = 8;
    if (use_schedule){
      Axis xo, xi, yo, yi, zo, zi;
      A.tile(tile_size_x, tile_size_y, tile_size_z, xo, xi, yo, yi, zo, zi);
      A.reorder(xo, yo, zo, xi, yi, zi);

      if (strcmp(target.c_str(), "sunway") == 0){
        CacheRead buffer_read;
        CacheWrite buffer_write;
        A.cache_read(B, buffer_read, "global");
        A.cache_write(buffer_write, "global");
        A.compute_at(buffer_read, zo);
        A.compute_at(buffer_write, zo);

        assert(num_threads == 64);
      }
      A.parallel(xo, num_threads);
    }
    A.build(target);

    Result Res((x, y, z), B[x, y, z]);
    DefShapeMPI3D(shape_mpi, D1, D2, D3);
    auto t = Stencil::t;
    Stencil stencil(shape_mpi, B, (x, y, z), Res[t] << A[t-1], 0, is_assigned);

    stencil.input(shape_mpi, B, input_file);
    stencil.run(1, 100);
    const bool use_timing = true;
    char func_name[100];
    sprintf(func_name, "3d31pt_star_%s_%d_%d_%d_%d_%d_%d_%d_%d_%d_%d_%dschedule", target.c_str(), num_threads, D1,D2,D3, M,N,P, tile_size_x, tile_size_y, tile_size_z, use_schedule);
    stencil.compile_to_source_code_mpi(func_name, target, use_timing);
}

int main(int argc, char * argv[]){

  if (argc == 2+8+2){
    void (*func)(std::string target, int num_threads, int D1, int D2,
        int M, int N, int tile_size_x,  int tile_size_y, int use_schedule, int is_assigned);
    if (strcmp(argv[1], "2d9pt_star") == 0)              func = bench_2d9pt_star;
    else if (strcmp(argv[1], "2d9pt_box") == 0)          func = bench_2d9pt_box;
    else if (strcmp(argv[1], "2d121pt_box") == 0)        func = bench_2d121pt_box;
    else if (strcmp(argv[1], "2d169pt_box") == 0)        func = bench_2d169pt_box;
    else assert(0 && "## 2d argv[1] name error");

    (*func)(std::string(argv[2]), 
            atoi(argv[3]),
            atoi(argv[4]),
            atoi(argv[5]),
            atoi(argv[6]),
            atoi(argv[7]),
            atoi(argv[8]),
            atoi(argv[9]),
            atoi(argv[10]),
            atoi(argv[11])
           );
  }

  else if (argc == 2+11+2){
   void (*func)(std::string target, int num_threads, int D1, int D2, int D3,
       int M, int N, int P, int tile_size_x,  int tile_size_y, int tile_size_z, int use_schedule, int is_assigned);
   if (strcmp(argv[1], "3d7pt_star") == 0)             func = bench_3d7pt_star;
   else if (strcmp(argv[1], "3d13pt_star") == 0)       func = bench_3d13pt_star;
   else if (strcmp(argv[1], "3d25pt_star") == 0)       func = bench_3d25pt_star;
   else if (strcmp(argv[1], "3d31pt_star") == 0)       func = bench_3d31pt_star;
   else assert(0 && "## 3d stencil name error");

   (*func)(std::string(argv[2]), 
           atoi(argv[3]),
           atoi(argv[4]),
           atoi(argv[5]),
           atoi(argv[6]),
           atoi(argv[7]),
           atoi(argv[8]),
           atoi(argv[9]),
           atoi(argv[10]),
           atoi(argv[11]),
           atoi(argv[12]),
           atoi(argv[13]),
           atoi(argv[14])
          );
  } 

  else if (argc == 1){
    bench_2d9pt_star("sunway", 64, 1, 1, 4096, 4096, 32, 64);
    //bench_2d9pt_box("sunway", 64, 1, 1, 4096, 4096, 32, 64);
    //bench_2d121pt_box("sunway", 64, 1, 1, 4096, 4096, 16, 32);
    //bench_2d169pt_box("sunway", 64, 1, 1, 4096, 4096, 16, 32);

    //bench_3d7pt_star("sunway", 64, 1, 1, 1, 256, 256, 256, 2, 8, 64);
    //bench_3d13pt_star("sunway", 64, 1, 1, 1, 256, 256, 256, 2, 8, 64);
    //bench_3d25pt_star("sunway", 64, 1, 1, 1, 256, 256, 256, 2, 4, 32);
    //bench_3d31pt_star("sunway", 64, 1, 1, 1, 256, 256, 256, 2, 4, 32);

    //bench_3d7pt_star("sunway", 64, 1, 1, 1, 256, 256, 256, 1, 8, 8);
    //bench_3d7pt_star("sunway", 64, 8, 4, 4, 1024, 32, 32, 1, 8, 8);
    //bench_3d7pt_star("sunway", 64, 8, 8, 8, 4096, 64, 32, 128, 8, 8);
  }
  else assert(0 && "## args error");

  return 0;
}

