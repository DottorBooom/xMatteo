#include "stencil_serial_nomp.h"

int main(int argc, char **argv)
{

  int  Niterations; // how many iterations
  int  periodic; // whether periodic boundaries apply
  unsigned int S[2]; // size of the plate

  int     Nsources; // how many heat sources
  int    *Sources; // the heat sources
  double  energy_per_source; // how much energy per source

  double *planes[2]; // the two planes, old and new
  
  double injected_heat = 0; // how much energy has been injected in the system

  int injection_frequency; // how often to inject energy
  int output_energy_at_steps = 0; // whether to print the energy budget at every step
   
  int ret; 

  /* argument checking and setting */
  ret = initialize ( argc, argv, &S[0], &periodic, &Niterations,
	       &Nsources, &Sources, &energy_per_source, &planes[0],
	       &output_energy_at_steps, &injection_frequency );

  if (ret == 1)
  {
    printf("Error during memory allocation\n");
    return 1;
  }
  int current = OLD;
  double time_spent_updating [Niterations];
  double time_spent_simulation [Niterations];

  if ( injection_frequency > 1 )
    inject_energy( periodic, Nsources, Sources, energy_per_source, S, planes[current]);

  for (int iter = 0; iter < Niterations; iter++)
  {

    // Start timing
    clock_t start = clock();

    /* new energy from sources */
    if ( iter % injection_frequency == 0 )
    {
      inject_energy( periodic, Nsources, Sources, energy_per_source, S, planes[current]);
      injected_heat += Nsources*energy_per_source;
    }
                  
    /* update grid points */
    update_plane(periodic, S, planes[current], planes[!current], &time_spent_updating[iter] );

    if ( output_energy_at_steps )
    {
      double system_heat;
      get_total_energy( S, planes[!current], &system_heat);
              
      printf("step %d :: injected energy is %g, updated system energy is %g\n", 
              iter, injected_heat, system_heat );

      char filename[100];
      sprintf( filename, "plane_%05d.bin", iter );
      //dump( planes[!current], S, filename, NULL, NULL );        
    }

    /* swap planes for the new iteration */
      current = !current;

    // End timing
    time_spent_simulation[iter] = (double)(clock() - start) / CLOCKS_PER_SEC;
  }

  double sumU = 0;
  double sumS = 0;
  for(int i = 0; i < Niterations; i++)
  {
    sumU += time_spent_updating[i];
    sumS += time_spent_simulation[i];
  }
  printf("Average time spent updating: %f ms\n", (sumU/Niterations)*1000);
  printf("Total time: %f s\n", sumS);

  /* get final heat in the system */
  double system_heat;
  get_total_energy( S, planes[current], &system_heat);

  printf("injected energy is %g, system energy is %g\n", injected_heat, system_heat );

  memory_release( planes[OLD], Sources );

  //export_and_plot(time_spent_injecting, time_spent_updating, time_spent_simulation, Niterations);

  return 0;
}