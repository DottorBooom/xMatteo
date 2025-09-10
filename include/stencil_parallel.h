/* -*- Mode: C; c-basic-offset:4 ; indent-tabs-mode:nil ; -*- */
/* See COPYRIGHT in top-level directory. */

#define _XOPEN_SOURCE
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <getopt.h>
#include <time.h>
#include <math.h>
#include <ctype.h>
#include <errno.h>
#include <sched.h>

#include <omp.h>
#include <mpi.h>

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

// tags for MPI communication
#define TAG_N 0
#define TAG_S 1
#define TAG_E 2
#define TAG_W 3

typedef unsigned int uint;

typedef uint    vec2_t[2];
typedef double *restrict buffers_t[4];

typedef struct {
    double *restrict data;
    vec2_t size;
} plane_t;

extern int inject_energy (const int,
                          const int,
			              const vec2_t *,
			              const double,
                          plane_t *,
                          const vec2_t);

extern int update_plane (const int,
                         const vec2_t,
                         const plane_t *,
                         plane_t *);

extern int get_total_energy(plane_t *,
                            double *);

int initialize (MPI_Comm *,
                int,
                int,
                int,
                char **,
                vec2_t *,
                vec2_t *,
                int *,
                int *,
                int *,
                int *,
                int *,
                int *,
                vec2_t **,
                double *,
                plane_t *,
                buffers_t * );

int memory_release (plane_t *planes,
		            buffers_t *buffers);

int output_energy_stat ( int,
                         plane_t *,
                         double,
                         int,
                         MPI_Comm *);

uint simple_factorization(uint, int *, uint **);

int initialize_sources(int,
                       int,
                       MPI_Comm *,
                       uint [2],
                       int,
                       int *,
                       vec2_t  **);

int memory_allocate (const int *,
		            buffers_t *,
		            plane_t *);

inline int inject_energy (const int periodic,
                          const int Nsources,
			              const vec2_t *Sources,
			              const double energy,
                          plane_t *plane,
                          const vec2_t N)
{
    register const int xsize = plane->size[_x_];
    register const int ysize = plane->size[_y_];
    register uint fxsize = plane->size[_x_]+2;

    double * restrict data = plane->data;

    #define IDX(i, j) ((j) * (fxsize) + (i))
    for (int s = 0; s < Nsources; s++)
    {
        int x = Sources[s][_x_];
        int y = Sources[s][_y_];
        
        data[IDX(x,y)] += energy;
        
        if (periodic)
        {
            if ((N[_x_] == 1))
            {
                data[IDX(0, y)] += data[IDX(xsize + 1, y)]; // West from East
                data[IDX(xsize + 1, y)] += data[IDX(1, y)]; // East from West
            }

            if ((N[_y_] == 1))
            {
                data[IDX(x, 0)] += data[IDX(x, ysize + 1)]; // North from South
                data[IDX(x, ysize + 1)] += data[IDX(x, 1)]; // South from North
            }
        }                
    }
    #undef IDX
    return 0;
}

inline int update_plane ( const int periodic, 
                          const vec2_t N,         // the grid of MPI tasks
                          const plane_t *oldplane,
                          plane_t *newplane)
{
    register uint fxsize = oldplane->size[_x_]+2;
    
    register uint xsize = oldplane->size[_x_];
    register uint ysize = oldplane->size[_y_];

    #define IDX( i, j ) ( (j)*fxsize + (i))
    
    // HINT: you may attempt to
    //       (i)  manually unroll the loop
    //       (ii) ask the compiler to do it
    // for instance
    // #pragma GCC unroll 4
    //
    // HINT: in any case, this loop is a good candidate
    //       for openmp parallelization

    double * restrict old = oldplane->data;
    double * restrict new = newplane->data;

    #pragma omp for collapse(2) schedule(static) nowait
    for (uint j = 1; j <= ysize; j++) {
        for (uint i = 1; i <= xsize; i++) 
        {
            double alpha = 0.6;
            double result = old[IDX(i,j)] * alpha;
            double sum_i = (old[IDX(i-1, j)] + old[IDX(i+1, j)]) / 4.0 * (1 - alpha);
            double sum_j = (old[IDX(i, j-1)] + old[IDX(i, j+1)]) / 4.0 * (1 - alpha);
            result += (sum_i + sum_j);
            new[IDX(i,j)] = result;
        }
    }

    if (periodic) {
        if (N[_x_] == 1) {
            for (uint i = 1; i <= xsize; i++) {
                new[i] = new[IDX(i, ysize)];
                new[IDX(i, ysize+1)] = new[i];
            }
        }
        if (N[_y_] == 1) {
            for (uint j = 1; j <= ysize; j++) {
                new[IDX(0, j)] = new[IDX(xsize, j)];
                new[IDX(xsize+1, j)] = new[IDX(1, j)];
            }
        }
    }
    #undef IDX
    return 0;
}

inline int get_total_energy( plane_t *plane,
                             double  *energy )
/*
 * NOTE: this routine a good candiadate for openmp
 *       parallelization
 */
{

    register const int xsize = plane->size[_x_];
    register const int ysize = plane->size[_y_];
    register const int fsize = xsize+2;

    double * restrict data = plane->data;
    
   #define IDX( i, j ) ( (j)*fsize + (i))

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
    #pragma omp parallel for collapse(2) reduction(+:totenergy) schedule(static)
    for ( int j = 1; j <= ysize; j++ )
        for ( int i = 1; i <= xsize; i++ )
            totenergy += data[ IDX(i, j) ];

    
   #undef IDX

    *energy = (double)totenergy;
    return 0;
}