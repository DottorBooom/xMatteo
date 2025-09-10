// Remove the comment below for core affinity check
//#define _GNU_SOURCE
#include "stencil_parallel.h"

int main(int argc, char **argv)
{
  MPI_Comm myCOMM_WORLD;
  int Rank, Ntasks;
  int neighbours[4];

  int Niterations;
  int periodic;
  vec2_t S, N; // S : global size of the plate, N : grid decomposition of MPI tasks

  int Nsources;
  int Nsources_local;
  vec2_t *Sources_local;
  double energy_per_source;

  plane_t planes[2];
  buffers_t buffers[2];
  
  int output_energy_stat_perstep = 0;

  double comm_time = 0.0, comp_time = 0.0, total_time, init_time, wait_time;
	double start_time_comm, start_time_comp;
  
  /* initialize MPI envrionment */
  {
    int level_obtained;
    
    // NOTE: change MPI_FUNNELED if appropriate
    MPI_Init_thread( &argc, &argv, MPI_THREAD_FUNNELED, &level_obtained );
    if ( level_obtained < MPI_THREAD_FUNNELED ) 
    {
      printf("MPI_thread level obtained is %d instead of %d\n",
	    level_obtained, MPI_THREAD_FUNNELED );
      MPI_Finalize();
      exit(1); 
    }
    
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &Ntasks);
    MPI_Comm_dup (MPI_COMM_WORLD, &myCOMM_WORLD);
  }

  init_time = MPI_Wtime();
  /* argument checking and setting */
  int ret = initialize(&myCOMM_WORLD, Rank, Ntasks, argc, argv, 
                       &S, &N, &periodic, &output_energy_stat_perstep,
			                 neighbours, &Niterations, &Nsources, &Nsources_local, 
                       &Sources_local, &energy_per_source, &planes[0], &buffers[0]);

  /* After calling initialize():
   *
   * myCOMM_WORLD         : Duplicated MPI communicator for this job
   * Rank                 : Rank of this MPI process (0 ... Ntasks-1)
   * Ntasks               : Total number of MPI processes
   * neighbours[4]        : Ranks of neighboring processes (NORTH, EAST, SOUTH, WEST), or MPI_PROC_NULL if absent
   * Niterations          : Number of simulation iterations
   * periodic             : 1 if periodic boundaries are enabled, 0 otherwise
   * S                    : Global grid size [X, Y]
   * N                    : Process grid decomposition [Nx, Ny]
   * Nsources             : Total number of heat sources
   * Nsources_local       : Number of sources assigned to this process
   * Sources_local        : Array of local source coordinates (size Nsources_local)
   * energy_per_source    : Energy injected per source
   * planes[2]            : Two data planes (OLD/NEW), each allocated with size (mysize[0]+2) x (mysize[1]+2)
   * buffers[2]           : Communication buffers for neighbors (SEND/RECV for N, S, E, W)
   * output_energy_stat_perstep : If >0, print energy at every step
   */

  if ( ret )
    {
      printf("task %d is opting out with termination code %d\n",
	     Rank, ret );
      
      MPI_Finalize();
      return 0;
    }
  
  
  int current = OLD;

  uint ysize = planes[current].size[_y_];
	uint xsize = planes[current].size[_x_];

  init_time = MPI_Wtime() - init_time;

  /*
  // --- BEGIN CORE AFFINITY CHECK BLOCK ---
  // Print the core ID for each OpenMP thread in each MPI task
  if (Rank == 0) {
      printf("\n--- Core Affinity Check ---\n");
      printf("MPI Rank | OMP Thread | Core ID\n");
      printf("---------------------------\n");
      fflush(stdout);
  }

  // Ordered loop to avoid jumbled output
  for (int i = 0; i < Ntasks; i++) {
      MPI_Barrier(myCOMM_WORLD);
      if (i == Rank) {
          #pragma omp parallel
          {
              int thread_id = omp_get_thread_num();
              int core_id = sched_getcpu();
                // Use a buffer to create the string and print it atomically
              char buffer[100];
              sprintf(buffer, "   %d\t |    %d\t  |   %d\n", Rank, thread_id, core_id);
              printf("%s", buffer);
              fflush(stdout);
          }
      }
  }
  MPI_Barrier(myCOMM_WORLD);
  if (Rank == 0) {
      printf("--- End Core Affinity Check ---\n\n");
      fflush(stdout);
  }
  // --- END OF CORE AFFINITY CHECK BLOCK ---
  */

    total_time = MPI_Wtime();

  #define IDX(i, j) ((j) * (xsize + 2) + (i))
  for (int iter = 0; iter < Niterations; ++iter)
  {
    MPI_Request reqs[8];
    int nreqs = 0;
    for (int i = 0; i < 8; ++i)
      reqs[i] = MPI_REQUEST_NULL;
    
    /* new energy from sources */
    inject_energy(periodic, Nsources_local, Sources_local, energy_per_source, &planes[current], N);

    // [A] fill the buffers, and/or make the buffers' pointers pointing to the correct position
    if (neighbours[WEST] != MPI_PROC_NULL && buffers[SEND][WEST] != NULL) 
        for (uint j = 0; j < ysize; ++j) 
            buffers[SEND][WEST][j] = planes[current].data[IDX(1, j+1)];

    if (neighbours[EAST] != MPI_PROC_NULL && buffers[SEND][EAST] != NULL) 
        for (uint j = 0; j < ysize; ++j) 
            buffers[SEND][EAST][j] = planes[current].data[IDX(xsize, j+1)];

    // [B] perfoem the halo communications
    //     (1) use Send / Recv
    //     (2) use Isend / Irecv
    //         --> can you overlap communication and compution in this way?

    start_time_comm = MPI_Wtime();

    // For NORTH and SOUTH, we use direct pointers to contiguous data (no separate allocation needed)
		if (neighbours[NORTH] != MPI_PROC_NULL) {
			buffers[SEND][NORTH] = &(planes[current].data[xsize + 3]); 		// the first effective row
			buffers[RECV][NORTH] = &(planes[current].data[1]);
		}
		if (neighbours[SOUTH] != MPI_PROC_NULL) {
			buffers[SEND][SOUTH] = &(planes[current].data[ysize * (xsize + 2) + 1]); // the last effective row
			buffers[RECV][SOUTH] = &(planes[current].data[(ysize + 1) * (xsize + 2) + 1]);
		}

    if (neighbours[EAST] != MPI_PROC_NULL) {
			// optimization: if the neighbor is the same rank, we can just copy the data
			if (neighbours[EAST] == Rank) {
				for (uint i = 0; i < ysize; i++) {
					buffers[RECV][EAST][i] = buffers[SEND][EAST][i];
				}
			} else {
				MPI_Isend(buffers[SEND][EAST], (int)ysize, MPI_DOUBLE, neighbours[EAST], TAG_E, myCOMM_WORLD, &reqs[nreqs++]);
				MPI_Irecv(buffers[RECV][EAST], (int)ysize, MPI_DOUBLE, neighbours[EAST], TAG_W, myCOMM_WORLD, &reqs[nreqs++]);
			}
		}
    if (neighbours[WEST] != MPI_PROC_NULL) {
			if (neighbours[WEST] == Rank) {
				for (uint i = 0; i < ysize; i++) {
					buffers[RECV][WEST][i] = buffers[SEND][WEST][i];
				}
			} else {
				MPI_Isend(buffers[SEND][WEST], (int)ysize, MPI_DOUBLE, neighbours[WEST], TAG_W, myCOMM_WORLD, &reqs[nreqs++]);
				MPI_Irecv(buffers[RECV][WEST], (int)ysize, MPI_DOUBLE, neighbours[WEST], TAG_E, myCOMM_WORLD, &reqs[nreqs++]);
			}
		}
    if (neighbours[NORTH] != MPI_PROC_NULL) {
			if (neighbours[NORTH] == Rank) {
				for (uint i = 0; i < xsize; i++) {
					buffers[RECV][NORTH][i] = buffers[SEND][NORTH][i];
				}
			} else {
				MPI_Isend(buffers[SEND][NORTH], (int)xsize, MPI_DOUBLE, neighbours[NORTH], TAG_N, myCOMM_WORLD, &reqs[nreqs++]);
				MPI_Irecv(buffers[RECV][NORTH], (int)xsize, MPI_DOUBLE, neighbours[NORTH], TAG_S, myCOMM_WORLD, &reqs[nreqs++]);
			}
		}
		if (neighbours[SOUTH] != MPI_PROC_NULL) {
			if (neighbours[SOUTH] == Rank) {
				for (uint i = 0; i < xsize; i++) {
					buffers[RECV][SOUTH][i] = buffers[SEND][SOUTH][i];
				}
			} else {
				MPI_Isend(buffers[SEND][SOUTH], (int)xsize, MPI_DOUBLE, neighbours[SOUTH], TAG_S, myCOMM_WORLD, &reqs[nreqs++]);
				MPI_Irecv(buffers[RECV][SOUTH], (int)xsize, MPI_DOUBLE, neighbours[SOUTH], TAG_N, myCOMM_WORLD, &reqs[nreqs++]);
			}
		}
    comm_time += MPI_Wtime() - start_time_comm;

    wait_time += MPI_Wtime();
		MPI_Waitall(nreqs, reqs, MPI_STATUSES_IGNORE);
		wait_time = MPI_Wtime() - wait_time;

    // [C] copy the haloes data

    start_time_comp = MPI_Wtime();
    if (neighbours[WEST] != MPI_PROC_NULL)
    for (uint j = 0; j < planes[current].size[_y_]; ++j)
        planes[current].data[IDX(0, j+1)] = buffers[RECV][WEST][j];

    if (neighbours[EAST] != MPI_PROC_NULL)
        for (uint j = 0; j < planes[current].size[_y_]; ++j)
            planes[current].data[IDX(xsize+1, j+1)] = buffers[RECV][EAST][j];

    /* --------------------------------------  */
    /* update grid points */

    update_plane(periodic, N, &planes[current], &planes[!current]);
    comp_time += MPI_Wtime() - start_time_comp;

    /* output if needed */
    if ( output_energy_stat_perstep )
      output_energy_stat (iter, &planes[!current], (iter+1) * Nsources*energy_per_source, Rank, &myCOMM_WORLD);

    /* swap plane indexes for the new iteration */
    current = !current;  
  }
  total_time = MPI_Wtime() - total_time;
  #undef IDX

  output_energy_stat (-1, &planes[!current], Niterations * Nsources*energy_per_source, Rank, &myCOMM_WORLD );
  memory_release(planes, buffers);
  
  double max_comp_time, max_comm_time, min_comp_time, sum_comp_time, max_wait_time, max_total_time;
    MPI_Reduce(&comp_time, &max_comp_time, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD);
    MPI_Reduce(&comm_time, &max_comm_time, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD);
    MPI_Reduce(&comp_time, &min_comp_time, 1, MPI_DOUBLE, MPI_MIN, 0, myCOMM_WORLD);
    MPI_Reduce(&comp_time, &sum_comp_time, 1, MPI_DOUBLE, MPI_SUM, 0, myCOMM_WORLD);
    MPI_Reduce(&wait_time, &max_wait_time, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD);
    MPI_Reduce(&total_time, &max_total_time, 1, MPI_DOUBLE, MPI_MAX, 0, myCOMM_WORLD);
    
    // Clean up the duplicated communicator
    MPI_Comm_free(&myCOMM_WORLD);

  if (Rank == 0) 
  {
        // Calcolo del load imbalance
        double avg_comp_time = sum_comp_time / Ntasks;
        double load_imbalance = 0.0;
        if (avg_comp_time > 0) {
            load_imbalance = (max_comp_time - min_comp_time) / avg_comp_time;
        }

        printf("Total time: %f\n", total_time);
        printf("Initialization time: %f\n", init_time);
        printf("Computation time: %f\n", comp_time);
        printf("Communication time: %f\n", comm_time);
        printf("Communication/Total ratio: %.2f%%\n", (comm_time/total_time)*100.0);
        printf("Computation/Total ratio: %.2f%%\n", (comp_time/total_time)*100.0);
        printf("Wait time for communication: %f\n", wait_time);
        printf("Other time (overhead): %f (%.2f%%)\n", 
            total_time - comp_time - comm_time - wait_time, 
            ((total_time - comp_time - comm_time - wait_time)/total_time)*100.0);
        printf("Max total time: %f\n", max_total_time);
        printf("Max computation time: %f\n", max_comp_time);
        printf("Max communication time: %f\n", max_comm_time);
        printf("Max wait time for communication: %f\n", max_wait_time);
        printf("Load imbalance: %f\n", load_imbalance);
        printf("Load balance efficiency: %f\n", (avg_comp_time/max_comp_time));
        printf("Communication efficiency: %f\n", (max_comp_time/max_total_time));
  }

  MPI_Finalize();
  return 0;
}

/* ==========================================================================
   =                                                                        =
   =   routines called within the integration loop                          =
   ========================================================================== */

/* ==========================================================================
   =                                                                        =
   =   initialization                                                       =
   ========================================================================== */		      

int initialize (MPI_Comm *Comm,
                int      Me,                  // the rank of the calling process
                int      Ntasks,              // the total number of MPI ranks
                int      argc,                // the argc from command line
                char   **argv,                // the argv from command line
                vec2_t  *S,                   // the size of the plane
                vec2_t  *N,                   // two-uint array defining the MPI tasks' grid
                int     *periodic,            // periodic-boundary tag
                int     *output_energy_stat,  // output energy statistics flag
                int     *neighbours,          // four-int array that gives back the neighbours of the calling task
                int     *Niterations,         // how many iterations
                int     *Nsources,            // how many heat sources
                int     *Nsources_local,      // how many heat sources are local to each MPI task
                vec2_t **Sources_local,       // the heat sources
                double  *energy_per_source,   // how much heat per source
                plane_t *planes,              // the two planes, old and new  
                buffers_t *buffers)           // comunication buffers
{
  int halt = 0;
  int ret = 0;
  int verbose = 0;

  // set default values

  (*S)[_x_]         = 1000;
  (*S)[_y_]         = 1000;
  *periodic         = 0;
  *Nsources         = 5;
  *Nsources_local   = 0;
  *Sources_local    = NULL;
  *Niterations      = 100;
  *energy_per_source = 1.0;

  if ( planes == NULL ) {
    printf("Error: invalid pointer passed to memory_allocate\n");
    return 1;
  }

  planes[OLD].size[0] = planes[OLD].size[1] = 0;
  planes[NEW].size[0] = planes[NEW].size[1] = 0;

  for ( int i = 0; i < 4; i++ )
    neighbours[i] = MPI_PROC_NULL;

  for ( int b = 0; b < 2; b++ )
    for ( int d = 0; d < 4; d++ )
      buffers[b][d] = NULL;

  // process the command line
  while ( 1 )
  {
    int opt;
    while((opt = getopt(argc, argv, ":h:x:y:e:E:n:o:p:v:")) != -1)
    {
	    switch( opt )
      {
      case 'x': (*S)[_x_] = (uint)atoi(optarg);
        break;

      case 'y': (*S)[_y_] = (uint)atoi(optarg);
        break;

      case 'e': *Nsources = atoi(optarg);
        break;

      case 'E': *energy_per_source = atof(optarg);
        break;

      case 'n': *Niterations = atoi(optarg);
        break;

      case 'o': *output_energy_stat = (atoi(optarg) > 0);
        break;

      case 'p': *periodic = (atoi(optarg) > 0);
        break;

      case 'v': verbose = atoi(optarg);
        break;

      case 'h': {
        if ( Me == 0 )
          printf( "\nvalid options are ( values btw [] are the default values ):\n"
            "-x    x size of the plate [10000]\n"
            "-y    y size of the plate [10000]\n"
            "-e    how many energy sources on the plate [4]\n"
            "-E    how many energy sources on the plate [1.0]\n"
            "-n    how many iterations [1000]\n"
            "-p    whether periodic boundaries applies  [0 = false]\n\n"
            "-o    whether to print the energy budgets at every step [0 = false]\n"
            "-v    verbosity level [0]\n"
            );
        halt = 1; }
        break;
      
      case ':': printf( "option -%c requires an argument\n", optopt);
        break;
        
      case '?': printf(" -------- help unavailable ----------\n");
        break;
      }
    }

    if ( opt == -1 )
      break;
  }

  if ( halt )
    return 1;

  /*
  * find a suitable domain decomposition
  * very simple algorithm, you may want to
  * substitute it with a better one
  *
  * the plane Sx x Sy will be solved with a grid
  * of Nx x Ny MPI tasks
  */

  vec2_t Grid;
  double formfactor = ((*S)[_x_] >= (*S)[_y_] ? (double)(*S)[_x_]/(*S)[_y_] : (double)(*S)[_y_]/(*S)[_x_] );
  int    dimensions = 2 - (Ntasks <= ((int)formfactor+1) );

  if ( dimensions == 1 )
  {
    if ( (*S)[_x_] >= (*S)[_y_] )
	    Grid[_x_] = Ntasks, Grid[_y_] = 1;
    else
	    Grid[_x_] = 1, Grid[_y_] = Ntasks;
  }
  else
  {
    int   Nf;
    uint *factors;
    uint  first = 1;
    ret = simple_factorization( Ntasks, &Nf, &factors );
     
    if (ret != 0) {
			printf("Error: factorization failed\n");
			return ret;
		}

    for ( int i = 0; (i < Nf) && ((Ntasks/first)/first > formfactor); i++ )
	    first *= factors[i];

		uint factor1 = first;
		uint factor2 = (uint)Ntasks/first;
		uint who_max = (factor1 > factor2);

		if ( (*S)[_x_] >= (*S)[_y_] ) {
			// wide data, make grid wide
			Grid[_x_] = who_max ? factor1 : factor2;
			Grid[_y_] = who_max ? factor2 : factor1;
		} else {
			// tall or square data, make grid tall
			Grid[_y_] = who_max ? factor1 : factor2;
			Grid[_x_] = who_max ? factor2 : factor1;
		}
		
		// Free the factors array allocated by simple_factorization
		if (factors != NULL) {
			free(factors);
		}
  }

  (*N)[_x_] = Grid[_x_];
  (*N)[_y_] = Grid[_y_];

  // my coordinates in the grid of processors
  uint X = Me % Grid[_x_];
  uint Y = Me / Grid[_x_];

  // find my neighbours
	if ( *periodic ) {
		// Horizontal neighbours
		if (Grid[_x_] > 1 || *periodic) {
			neighbours[EAST]  = (int)(((uint)Y)*Grid[_x_] + (uint)(X + 1 ) % Grid[_x_]);
			neighbours[WEST]  = (int)(((uint)Y)*Grid[_x_] + (uint)(X - 1 + (int)Grid[_x_]) % Grid[_x_]);
		}
		// Vertical neighbours
		if (Grid[_y_] > 1 || *periodic) {
			neighbours[NORTH] = (int)(((uint)(Y - 1 + (int)Grid[_y_]) % Grid[_y_]) * Grid[_x_] + (uint)X);
			neighbours[SOUTH] = (int)(((uint)(Y + 1) % Grid[_y_]) * Grid[_x_] + (uint)X);
		}
	} else {
		// Horizontal neighbours
		if ( Grid[_x_] > 1 ) {  
			neighbours[EAST]  = ( X < (int)Grid[_x_]-1 ? Me+1 : MPI_PROC_NULL );
			neighbours[WEST]  = ( X > 0 ? Me-1 : MPI_PROC_NULL ); 
		}
		// Vertical neighbours
		if ( Grid[_y_] > 1 ) {
			neighbours[NORTH] = ( Y > 0 ? Me - (int)Grid[_x_]: MPI_PROC_NULL );
			neighbours[SOUTH] = ( Y < (int)Grid[_y_]-1 ? Me + (int)Grid[_x_] : MPI_PROC_NULL );
		}
	}

  // the size of my patch

  /*
  * every MPI task determines the size sx x sy of its own domain
  * REMIND: the computational domain will be embedded into a frame
  *         that is (sx+2) x (sy+2)
  *         the outern frame will be used for halo communication or
  */
  
  vec2_t mysize;
  uint s = (*S)[_x_] / Grid[_x_];
  uint r = (*S)[_x_] % Grid[_x_];
  mysize[_x_] = s + (X < r);
  s = (*S)[_y_] / Grid[_y_];
  r = (*S)[_y_] % Grid[_y_];
  mysize[_y_] = s + (Y < r);

  planes[OLD].size[0] = mysize[0];
  planes[OLD].size[1] = mysize[1];
  planes[NEW].size[0] = mysize[0];
  planes[NEW].size[1] = mysize[1];
  

  if ( verbose > 0 )
  {
    if ( Me == 0 ) {
		  printf("Tasks are decomposed in a grid %d x %d\n\n",
			Grid[_x_], Grid[_y_] );
		  fflush(stdout);
	  }
    MPI_Barrier(*Comm);
    
    if (Me == 0) {
      printf("Neighbours:\n\n");
      printf("   Task   N     E     S     W\n");
      printf("  ============================\n");
      fflush(stdout);
    }

    MPI_Barrier(*Comm);
    for (int t = 0; t < Ntasks; t++) {
      if (t == Me) {
        printf("%5d %5d %5d %5d %5d\n",
          Me,
          neighbours[NORTH],
          neighbours[EAST],
          neighbours[SOUTH],
          neighbours[WEST]
        );
        fflush(stdout);
      }
      MPI_Barrier(*Comm);
    }
    if (Me == 0) {
      printf("\n");
      fflush(stdout);
    }
    MPI_Barrier(*Comm);
  }

  // allocate the needed memory
  ret = memory_allocate(neighbours, buffers, planes);
	if (ret != 0) {
		printf("Error: failed to allocate memory for the buffers\n");
		return ret;
	}

  // allocate the heat sources
  ret = initialize_sources( Me, Ntasks, Comm, mysize, *Nsources, Nsources_local, Sources_local );
	if (ret != 0) {
		printf("Error: failed to initialize the sources\n");
		return ret;
	}

  return ret;
}

uint simple_factorization( uint A, int *Nfactors, uint **factors )
{
    /*
     * Algoritmo di fattorizzazione corretto e ottimizzato.
     * Esegue la scomposizione in un unico passaggio e gestisce tutti i casi.
     */
    uint temp_A = A;
    uint f = 2;
    int N = 0;
    // Usa un array temporaneo per memorizzare i fattori. 64 è un limite sicuro
    // dato che 2^64 è un numero enorme.
    uint temp_factors[64];

    // Scompone il numero e memorizza i fattori nell'array temporaneo
    while (f * f <= temp_A) {
        while (temp_A % f == 0) {
            temp_factors[N++] = f;
            temp_A /= f;
        }
        f++;
    }
    // Se è rimasto un numero > 1, è un fattore primo
    if (temp_A > 1) {
        temp_factors[N++] = temp_A;
    }

    *Nfactors = N;
    if (N == 0) {
        *factors = NULL;
        return 0; // Successo
    }

    // Alloca la memoria della dimensione esatta necessaria
    uint *_factors_ = (uint*)malloc((size_t)N * sizeof(uint));
    if (_factors_ == NULL) {
        perror("simple_factorization: malloc failed");
        return 1; // Errore
    }

    // Copia i fattori dall'array temporaneo a quello finale
    for (int i = 0; i < N; i++) {
        _factors_[i] = temp_factors[i];
    }

    *factors = _factors_;
    return 0; // Successo
}

int memory_allocate (const int *neighbours,
                     buffers_t *buffers_ptr,
                     plane_t *planes_ptr)
{
  /*
    here you allocate the memory buffers that you need to
    (i)  hold the results of your computation
    (ii) communicate with your neighbours

    The memory layout that I propose to you is as follows:

    (i) --- calculations
    you need 2 memory regions: the "OLD" one that contains
    results for the step (i-1)th, and the "NEW" one that will contain
    the updated results from the step ith.

    Then, the "NEW" will be treated as "OLD" and viceversa.

    These two memory regions are indexed by *plate_ptr:

    planew_ptr[0] ==> the "OLD" region
    plamew_ptr[1] ==> the "NEW" region


    (ii) --- communications

    you may need two buffers (one for sending and one for receiving)
    for each one of your neighnours, that are at most 4:
    north, south, east amd west.      

    To them you need to communicate at most mysizex or mysizey
    daouble data.

    These buffers are indexed by the buffer_ptr pointer so
    that

    (*buffers_ptr)[SEND][ {NORTH,...,WEST} ] = .. some memory regions
    (*buffers_ptr)[RECV][ {NORTH,...,WEST} ] = .. some memory regions
    
    --->> Of course you can change this layout as you prefer
  */

  if (planes_ptr == NULL ){
    printf("Error: invalid pointer passed to memory_allocate\n");
    return 1;
  }

  if (buffers_ptr == NULL ){
    printf("Error: invalid pointer passed to memory_allocate\n");
    return 1;
  }

  // ··················································
  // allocate memory for data
  // we allocate the space needed for the plane plus a contour frame
  // that will contains data form neighbouring MPI tasks

  /*
  
  [ X  R  R  R  X ] 
  [ S  D  D  D  R ]
  [ S  D  D  D  R ]
  [ S  D  D  D  R ]
  [ x  R  R  R  X ]

  D = Data; R = Halo for reciving; x = Contour frame, never used
  */
 
  uint sx = planes_ptr[OLD].size[_x_];
  uint sy = planes_ptr[OLD].size[_y_];

  unsigned int frame_size = ((sx + 2) * (sy + 2)); // -4 to exclude the corners

  planes_ptr[OLD].data = (double*)aligned_alloc(64, frame_size * sizeof(double));
  if ( planes_ptr[OLD].data == NULL ){
    printf("Error: aligned_alloc failed for OLD plane data\n");
    return 2;
  }
  //  memset ( planes_ptr[OLD].data, 0, frame_size * sizeof(double) );
	#pragma omp parallel for schedule(static)
	for (uint i = 0; i < frame_size; i++) {
		planes_ptr[OLD].data[i] = 0.0;
	}

  planes_ptr[NEW].data = (double*)aligned_alloc(64, frame_size * sizeof(double));
  if ( planes_ptr[NEW].data == NULL ){
    printf("Error: aligned_alloc failed for NEW plane data\n");
    return 2;
  }
  //memset ( planes_ptr[NEW].data, 0, frame_size * sizeof(double) );
  #pragma omp parallel for schedule(static)
	for (unsigned int i = 0; i < frame_size; i++) {
		planes_ptr[NEW].data[i] = 0.0;
	}

  planes_ptr[OLD].size[_x_] = sx;
	planes_ptr[OLD].size[_y_] = sy;
	planes_ptr[NEW].size[_x_] = sx;
	planes_ptr[NEW].size[_y_] = sy;

  // ··················································
  // buffers for north and south communication 
  // are not really needed
  //
  // in fact, they are already contiguous, just the
  // first and last line of every rank's plane
  //
  // you may just make some pointers pointing to the
  // correct positions

  // or, if you prefer, just go on and allocate buffers
  // also for north and south communications

  // ··················································
  // allocate buffers
  //

	if (neighbours[EAST] != MPI_PROC_NULL) {
		buffers_ptr[SEND][EAST] = (double*)aligned_alloc(64, sy * sizeof(double));
		buffers_ptr[RECV][EAST] = (double*)aligned_alloc(64, sy * sizeof(double));
		if (buffers_ptr[SEND][EAST] == NULL || buffers_ptr[RECV][EAST] == NULL) {
			printf("Error: failed to allocate memory for the EAST buffers\n");
			return 1;
		}
	}
	if (neighbours[WEST] != MPI_PROC_NULL) {
		buffers_ptr[SEND][WEST] = (double*)aligned_alloc(64, sy * sizeof(double));
		buffers_ptr[RECV][WEST] = (double*)aligned_alloc(64, sy * sizeof(double));
		if (buffers_ptr[SEND][WEST] == NULL || buffers_ptr[RECV][WEST] == NULL) {
			printf("Error: failed to allocate memory for the WEST buffers\n");
			return 1;
		}
	}

  return 0;
}

int initialize_sources(int Me,
                       int Ntasks,
                       MPI_Comm *Comm,
                       vec2_t mysize,
                       int Nsources,
                       int *Nsources_local,
                       vec2_t **Sources)
{
  srand48(0 + Me);
  int *tasks_with_sources = (int*)malloc( Nsources * sizeof(int) );
  
  if ( Me == 0 )
  {
    for ( int i = 0; i < Nsources; i++ )
  	  tasks_with_sources[i] = (int)lrand48() % Ntasks;
  }
  
  MPI_Bcast( tasks_with_sources, Nsources, MPI_INT, 0, *Comm );

  int nlocal = 0;
  for ( int i = 0; i < Nsources; i++ )
    nlocal += (tasks_with_sources[i] == Me);
  *Nsources_local = nlocal;
  
  if ( nlocal > 0 )
  {
    vec2_t * restrict helper = (vec2_t*)malloc( nlocal * sizeof(vec2_t) );      
    for ( int s = 0; s < nlocal; s++ )
	  {
	    helper[s][_x_] = 1 + lrand48() % mysize[_x_];
	    helper[s][_y_] = 1 + lrand48() % mysize[_y_];
	  }
    *Sources = helper;
  }
  
  free(tasks_with_sources);
  return 0;
}

int memory_release (plane_t *planes,
		                buffers_t *buffers)
{
  if ( planes != NULL )
  {
    if ( planes[OLD].data != NULL )
	    free (planes[OLD].data);
      
    if ( planes[NEW].data != NULL )
	    free (planes[NEW].data);
  }

  if ( buffers != NULL )
  {
    for ( int b = 0; b < 2; b++ )
      for ( int d = 2; d < 4; d++ )
      {
        if ( buffers[b][d] != NULL )
          free( buffers[b][d] );
      }
  }
      
  return 0;
}

int output_energy_stat ( int step, plane_t *plane, double budget, int Me, MPI_Comm *Comm )
{

  double system_energy = 0;
  double tot_system_energy = 0;
  get_total_energy ( plane, &system_energy );
  
  MPI_Reduce ( &system_energy, &tot_system_energy, 1, MPI_DOUBLE, MPI_SUM, 0, *Comm );
  
  if ( Me == 0 )
  {
    if ( step >= 0 ){
	    printf(" [ step %4d ] ", step ); 
      fflush(stdout);
    }

    printf( "total injected energy is %g, "
            "system energy is %g "
            "( in avg %g per grid point)\n",
            budget,
            tot_system_energy,
            tot_system_energy / (plane->size[_x_]*plane->size[_y_]) );
  }
  
  return 0;
}