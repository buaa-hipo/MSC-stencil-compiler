
const std::string headers_cpu =
    "#include <stdio.h>\n"
    "#include <stdlib.h>\n"
    "#include <math.h>\n"
    "#include <time.h>\n"
    "#include <sys/time.h>\n"
    "#include <string.h>\n"
    "#include <stdint.h>\n"
    "#include \"../IO.h\"\n"
    "\n"
    "double get_time() {\n"
    "  struct timeval tv;\n"
    "  gettimeofday(&tv, NULL);\n"
    "  return tv.tv_sec + 1e-6 * tv.tv_usec;\n"
    "}\n"
    "double time_start;\n"
    "double time_end;\n"
    "\n";

const std::string timing_start_mpi =
    "\n"
    " MPI_Barrier(MPI_COMM_WORLD);\n"
    " if(mpi_rank==0){\n"
    "   time_start = MPI_Wtime();\n"
    " }\n"
    ;

const std::string timing_end_mpi =
    " MPI_Barrier(MPI_COMM_WORLD);\n"
    " if(mpi_rank==0){\n"
    "   time_end = MPI_Wtime();\n"
    "   printf(\"\\n mpi rank: %d\\n total millisecond: %.2f\\n timestep: %d\\n avg millesecond: %.2f\\n\",mpi_rank,(time_end-time_start)*1000,timesteps,(time_end-time_start)/(timesteps)*1000);\n"
    " }\n"
    "\n";

const std::string timing_start =
    "\n"
    " time_start = get_time();\n"
    ;

const std::string timing_end =
    " time_end = get_time();\n"
    " printf(\"\\n total millisecond: %.2f\\n timestep: %d\\n avg millesecond: %.2f\\n\",(time_end-time_start)*1000,timesteps,(time_end-time_start)/(timesteps)*1000);\n"
    "\n"
    "\n";


const std::string headers_omp = 
    "#include <omp.h>\n"
    "\n";


const std::string headers_mpi = 
    "#include \"../mpi_lib/mpi_lib.h\"\n"
    "int mpi_rank ;\n"
    "int mpi_size ;\n"
    "\n";

const std::string headers_athread = 
    "#include <athread.h>\n"
    "\n";

const std::string headers_slave = 
    "#include \"slave.h\"\n"
    "#include <stdio.h>\n"
    "#include <stdlib.h>\n"
    "#include <math.h>\n"
    "#include <string.h>\n"
    "#include <stdint.h>\n"
    "\n"
    "__thread_local volatile unsigned long get_reply, put_reply;\n"
    "__thread_local volatile int my_id;\n"
    "\n";

const std::string headers_argument_host_slave = 
    "typedef struct Argfloat {\n"
    "  float * TeNode;\n"
    "  float * SpNodeTemp;\n"
    "} Arg_float;\n"
    "typedef struct Argdouble {\n"
    "  float * TeNode;\n"
    "  float * SpNodeTemp;\n"
    "} Arg_double;\n"
    "typedef struct Argint {\n"
    "  float * TeNode;\n"
    "  float * SpNodeTemp;\n"
    "} Arg_int;\n"
    "\n";

const std::string headers_dma_get_put_slave = 
    "#define Get(x, y, cnt) { \\\n"
    "  volatile unsigned long get_reply = 0 ; \\\n"
    "  athread_get(PE_MODE, x, y, cnt, &get_reply, 0, 0, 0); \\\n"
    "  while(get_reply != 1) ; \\\n"
    "}  \n"
    "#define Put(x, y, cnt) { \\\n"
    "  volatile unsigned long put_reply = 0 ; \\\n"
    "  athread_put(PE_MODE, y, x, cnt, &put_reply, 0, 0); \\\n"
    "  while(put_reply != 1) ; \\\n"
    "}  \n"
    "#define Get_Stride(x, y, z_dim_size, cnt, stride, bsize) { \\\n"
    "    volatile unsigned long get_reply = 0 ; \\\n"
    "    int z_iter; \\\n"
    "    for (z_iter = 0; z_iter < z_dim_size; z_iter++) \\\n"
    "       athread_get(PE_MODE, x, y, cnt, &get_reply, 0, stride, bsize); \\\n"
    "    while(get_reply != z_dim_size) ; \\\n"
    "}  \n"  
    "#define Put_Stride(x, y, z_dim_size, cnt, stride, bsize) { \\\n"
    "    volatile unsigned long put_reply = 0 ; \\\n"
    "    int z_iter; \\\n"
    "    for (z_iter = 0; z_iter < z_dim_size; z_iter++) \\\n"
    "       athread_put(PE_MODE, y, x, cnt, &put_reply, stride, bsize); \\\n"
    "    while(put_reply != z_dim_size) ; \\\n"
    "} \n"
    "\n";




