/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/* See COPYRIGHT in top-level directory. */

#define _XOPEN_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <float.h>
#include <math.h>

#define NORTH 0
#define SOUTH 1
#define EAST  2
#define WEST  3

#define SEND 0
#define RECV 1

#define OLD 0
#define NEW 1

#define _x_ 0
#define _y_ 1

/* ==========================================================================
   =   Function Prototypes                                                  =
   ========================================================================== */
   
int initialize (int,
                char **,
                unsigned int *,
                int *,
                int *,
                int *,
                int **,
                double *,
                double **,
                int *,
                int *);

int memory_release (double *, int *);


extern int inject_energy (  int,
                            int,
                            int *,
                            double,
                            unsigned int [2],
                            double *); 

extern int update_plane (   int,
                            unsigned int [2],
                            double *,
		                    double *,
                            double *); //Timestamp

extern int get_total_energy( unsigned int [2],
                             double *,
                             double *);


int dump (  const double *, 
            unsigned int [2], 
            const char *, 
            double *, 
            double * );

int memory_allocate ( unsigned int [2],
		                  double **);

int initialize_sources( unsigned int [2],
			                  int ,
			                  int **);

void export_and_plot(double *, double *, double *,int);

/* ==========================================================================
   =   Main Function and Simulation Loop                                    =
   ========================================================================== */

inline int inject_energy (  int periodic,
                            int Nsources,
			                      int *Sources,
			                      double energy,
			                      unsigned int mysize[2],
                            double *plane)
{

    #define IDX( i, j ) ( (j)*(mysize[_x_]+2) + (i) )
    
    // Start timing
    // clock_t begin = clock();

    for (int s = 0; s < Nsources; s++) {
        
        unsigned x = Sources[2*s];
        unsigned y = Sources[2*s+1];
        plane[IDX(x, y)] += energy;

        if ( periodic )
            {
                if ( x == 1 )
                    plane[IDX(mysize[_x_]+1, y)] += energy;
                if ( x == mysize[_x_] )
                    plane[IDX(0, y)] += energy;
                if ( y == 1 )
                    plane[IDX(x, mysize[_y_]+1)] += energy;
                if ( y == mysize[_y_] )
                    plane[IDX(x, 0)] += energy;
            }
    }
    // End timing
    // clock_t end = clock();
    // *time_spent = (double)(end - begin) / CLOCKS_PER_SEC;

    #undef IDX
    return 0;
}

inline int update_plane (   int periodic,
                            unsigned int size[2],
			                double *old,
                            double *new,
                            double *time_spent)

/*
 * calculate the new energy values
 * the old plane contains the current data, the new plane
 * will store the updated data
 *
 * NOTE: in parallel, every MPI task will perform the
 *       calculation for its patch
 *
 */
{
    register const unsigned fxsize = size[_x_]+2;
    //register const unsigned fysize = size[_y_]+2;
    register const unsigned xsize = size[_x_];
    register const unsigned ysize = size[_y_];

    #define IDX( i, j ) ( (j)*fxsize + (i) )

    // HINT: you may attempt to
    //       (i)  manually unroll the loop
    //       (ii) ask the compiler to do it
    // for instance
    // #pragma GCC unroll 4
    //
    // HINT: in any case, this loop is a good candidate
    //       for openmp parallelization

    // Start timing
    clock_t begin = clock();

    for (unsigned int j = 1; j <= ysize; j++)
        for (unsigned int i = 1; i <= xsize; i++)
        {
            // five-points stencil formula
            // simpler stencil with no explicit diffusivity
            // always conserve the smoohed quantity
            // alpha here mimics how much "easily" the heat
            // travels
            
            double alpha = 0.6;
            double result = old[ IDX(i,j) ] *alpha;
            double sum_i  = (old[IDX(i-1, j)] + old[IDX(i+1, j)]) / 4.0 * (1-alpha);
            double sum_j  = (old[IDX(i, j-1)] + old[IDX(i, j+1)]) / 4.0 * (1-alpha);
            result += (sum_i + sum_j );
            
            // implentation from the derivation of
            // 3-points 2nd order derivatives
            // however, that should depends on an adaptive
            // time-stepping so that given a diffusivity
            // coefficient the amount of energy diffused is
            // "small"
            // however the implicit methods are not stable

            /*
            #define alpha_guess 0.5     // mimic the heat diffusivity

            double alpha = alpha_guess;
            double sum = old[IDX(i,j)];
            double result = old[ IDX(i,j) ] *alpha;

            int   done = 0;
            do
                {                
                    double sum_i = alpha * (old[IDX(i-1, j)] + old[IDX(i+1, j)] - 2*sum);
                    double sum_j = alpha * (old[IDX(i, j-1)] + old[IDX(i, j+1)] - 2*sum);
                    result = sum + ( sum_i + sum_j);
                    double ratio = fabs((result-sum)/(sum!=0? sum : 1.0));
                    done = ( (ratio < 2.0) && (result >= 0) );    // not too fast diffusion and
                                                                    // not so fast that the (i,j)
                                                                    // goes below zero energy
                    alpha /= 2;
                }
            while ( !done );
            */

            new[ IDX(i,j) ] = result;
        }
    
    // End timing
    clock_t end = clock();
    *time_spent = (double)(end - begin) / CLOCKS_PER_SEC;
        
    if ( periodic )
    /*
    * propagate boundaries if they are periodic
    *
    * NOTE: when is that needed in distributed memory, if any?
    */
    {
        // Vertical boundaries top <-> bottom
        for (unsigned int i = 1; i <= xsize; i++ )
            {
                new[ i ] = new[ IDX(i, ysize) ];
                new[ IDX(i, ysize+1) ] = new[ i ];
            }
        // Horizontal boundaries left <-> right
        for (unsigned int j = 1; j <= ysize; j++ )
            {
                new[ IDX( 0, j) ] = new[ IDX(xsize, j) ];
                new[ IDX( xsize+1, j) ] = new[ IDX(1, j) ];
            }
    }

    #undef IDX
    return 0;
}

inline int get_total_energy(unsigned int size[2],
                            double *plane,
                            double *energy)
/*
 * NOTE: this routine a good candiadate for openmp
 *       parallelization
 */
{
    register const int xsize = size[_x_];
    
    #define IDX( i, j ) ( (j)*(xsize+2) + (i) )

    #if defined(LONG_ACCURACY)    
        long double totenergy = 0;
    #else
        double totenergy = 0;    
    #endif

    // HINT: you may attempt to
    //       (i)  manually unroll the loop
    //       (ii) ask the compiler to do it
    // for instance
    // #pragma GCC unroll 4

    for (unsigned int j = 1; j <= size[_y_]; j++ )
        for (unsigned int i = 1; i <= size[_x_]; i++ )
            totenergy += plane[ IDX(i, j) ];
    
    #undef IDX

    *energy = (double)totenergy;
    return 0;
}
                            
/* ==========================================================================
   =   Initialization of Variables and Memory Allocation                    =
   ========================================================================== */

int initialize (int argc,                 // the argc from command line
                char **argv,              // the argv from command line
                unsigned int *S,          // two-unsigned-int array defining the x,y dimensions of the grid
                int *periodic,            // periodic-boundary tag
                int *Niterations,         // how many iterations
                int *Nsources,            // how many heat sources
                int **Sources,
                double *energy_per_source,// how much heat per source
                double **planes,
                int *output_energy_at_steps,
                int *injection_frequency)
{
  int ret;
  
  // set default values

  S[_x_]            = 1000;
  S[_y_]            = 1000;
  *periodic         = 0;
  *Nsources         = 1;
  *Niterations      = 99;
  *output_energy_at_steps = 0;
  *energy_per_source = 1.0;
  *injection_frequency = *Niterations;

  double freq = 0;

  // process the command line

  while ( 1 )
  {
    int opt;
    while((opt = getopt(argc, argv, ":x:y:e:E:f:n:p:o:")) != -1)
    {
      switch( opt )
      {
        case 'x': S[_x_] = (unsigned int)atoi(optarg);
          break;

        case 'y': S[_y_] = (unsigned int)atoi(optarg);
          break;

        case 'e': *Nsources = atoi(optarg);
          break;

        case 'E': *energy_per_source = atof(optarg);
          break;

        case 'n': *Niterations = atoi(optarg);
          break;

        case 'p': *periodic = (atoi(optarg) > 0);
          break;

        case 'o': *output_energy_at_steps = (atoi(optarg) > 0);
          break;

        case 'f': freq = atof(optarg);
          break;
          
        case 'h': printf( "valid options are ( values btw [] are the default values ):\n"
              "-x    x size of the plate [1000]\n"
              "-y    y size of the plate [1000]\n"
              "-e    how many energy sources on the plate [1]\n"
              "-E    how much energy per sources [1.0]\n"
              "-f    the frequency of energy injection [0.0]\n"
              "-n    how many iterations [100]\n"
              "-p    whether periodic boundaries applies  [0 = false]\n"
              "-o    whether to print the energy budgest at every step [0 = false]\n"
              );
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

  if ( freq == 0 )
    *injection_frequency = 1;
  else
  {
    freq = (freq > 1.0 ? 1.0 : freq );
    *injection_frequency = freq * *Niterations;
  }

  // allocate the needed memory
  ret = memory_allocate( S, planes ); 
  if ( ret != 0 )
    return ret;

  // allocate the heat sources
  ret = initialize_sources( S, *Nsources, Sources );
  if ( ret != 0 )
    return ret;

  return 0;
}


int memory_allocate ( unsigned int size[2],
		                  double **planes_ptr )
/*
 * allocate the memory for the planes
 * we need 2 planes: the first contains the
 * current data, the second the updated data
 *
 * in the integration loop then the roles are
 * swapped at every iteration
 *
 */
{
  if (planes_ptr == NULL ){
      // an invalid pointer has been passed
      // manage the situation
      printf("Error: invalid pointer passed to memory_allocate\n");
      return 1;
    }
  unsigned int bytes = (size[_x_]+2)*(size[_y_]+2);

  planes_ptr[OLD] = (double*)malloc( 2*bytes*sizeof(double) );
  memset ( planes_ptr[OLD], 0, 2*bytes*sizeof(double) );
  planes_ptr[NEW] = planes_ptr[OLD] + bytes;
      
  return 0;
}


int initialize_sources( unsigned int size[2],
			                  int Nsources,
			                  int **Sources)
/*
 * randomly spread heat suources
 *
 */
{
  *Sources = (int*)malloc( Nsources * 2 *sizeof(unsigned int) );
  for ( int s = 0; s < Nsources; s++ )
  {
    (*Sources)[s*2] = 1+ lrand48() % size[_x_];
    (*Sources)[s*2+1] = 1+ lrand48() % size[_y_];
  }
  return 0;
}

/* ==========================================================================
   =   Memory Release                                                       =
   ========================================================================== */

int memory_release ( double *data, int *sources )  
{
  if( data != NULL )
    free( data );

  if( sources != NULL )
    free( sources );

  return 0;
}

int dump ( const double *data, unsigned int size[2], const char *filename, double *min, double *max )
{
  if ( (filename != NULL) && (filename[0] != '\0') )
  {
    FILE *outfile = fopen( filename, "w" );
    if ( outfile == NULL )
	    return 2;
      
    float *array = (float*)malloc( size[0] * sizeof(float) );
      
    double _min_ = DBL_MAX;
    double _max_ = 0;

    for (unsigned int j = 0; j < size[1]; j++ )
	  {
      /*
      float y = (float)j / size[1];
      fwrite ( &y, sizeof(float), 1, outfile );
      */
      
      const double * restrict line = data + j*size[0];
      for (unsigned int i = 0; i < size[0]; i++ )
      {
        array[i] = (float)line[i];
        _min_ = ( line[i] < _min_? line[i] : _min_ );
        _max_ = ( line[i] > _max_? line[i] : _max_ ); 
      }
      
      fwrite( array, sizeof(float), size[0], outfile );
	  }
      
    free( array );
      
    fclose( outfile );

    if ( min != NULL )
	    *min = _min_;
    if ( max != NULL )
	    *max = _max_;
  }

  else return 1;

  return 0;
}

/* ==========================================================================
   =   Plot                                                                 =
   ========================================================================== */

void export_and_plot(double *injecting, double *updating, double *simulation, int Niterations) {
    FILE *fp = fopen("timing_data.csv", "w");
    fprintf(fp, "iteration,injecting,updating,total\n");

    double accI = 0.0, accU = 0.0, accS = 0.0;
    for (int i = 0; i < Niterations; i++) {
        accI += injecting[i];
        accU += updating[i];
        accS += simulation[i];
        fprintf(fp, "%d,%.8f,%.8f,%.8f\n", i, accI, accU, accS);
    }
    fclose(fp);

    
    FILE *gp = fopen("plot.gp", "w");
    fprintf(gp,
        "set datafile separator \",\"\n"
        "set title \"Cumulative Execution Time\"\n"
        "set xlabel \"Iteration\"\n"
        "set ylabel \"Cumulative Time (s)\"\n"
        "set grid\n"
        "set key outside\n"
        "plot \"timing_data.csv\" using 1:2 with lines lw 2 title \"Injecting\", \\\n"
        "     \"timing_data.csv\" using 1:3 with lines lw 2 title \"Updating\", \\\n"
        "     \"timing_data.csv\" using 1:4 with lines lw 2 title \"Total\"\n"
        "pause -1 \"Press Enter to close\"\n");
    fclose(gp);

    system("gnuplot plot.gp");
}