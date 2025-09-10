#include "stencil_serial.h"

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

  double init_time, total_time;

  init_time = omp_get_wtime();
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

  int nthreads = omp_get_num_threads();
  double thread_time [nthreads];

  if ( injection_frequency > 1 )
    inject_energy( periodic, Nsources, Sources, energy_per_source, S, planes[current]);

  init_time = omp_get_wtime() - init_time;
  total_time = omp_get_wtime();

  for (int iter = 0; iter < Niterations; iter++)
  {

    /* new energy from sources */
    if ( iter % injection_frequency == 0 )
    {
      inject_energy( periodic, Nsources, Sources, energy_per_source, S, planes[current]);
      injected_heat += Nsources*energy_per_source;
    }
                  
    /* update grid points */
    update_plane(periodic, S, planes[current], planes[!current], &time_spent_updating[iter], thread_time);

    if ( output_energy_at_steps )
    {
      double system_heat;
      get_total_energy( S, planes[!current], &system_heat);
              
      printf("step %d :: injected energy is %g, updated system energy is %g\n", 
              iter, injected_heat, system_heat );

      //char filename[100];
      //sprintf( filename, "plane_%05d.bin", iter );
      //dump( planes[!current], S, filename, NULL, NULL );        
    }

    /* swap planes for the new iteration */
      current = !current;

  }
  total_time = omp_get_wtime() - total_time;

  /* get final heat in the system */
  double system_heat;
  get_total_energy( S, planes[current], &system_heat);
  printf("injected energy is %g, system energy is %g\n", injected_heat, system_heat );
  memory_release( planes[OLD], Sources );

  double sumU = 0;

  #pragma omp parallel for reduction(+: sumU) schedule(static)
  for(int i = 0; i < Niterations; i++)
  {
    sumU += time_spent_updating[i];
  }

  // Calcolo load imbalance per l'ultima iterazione 
  double max_time = thread_time[0], min_time = thread_time[0], avg_time = 0.0;
  for (int i = 0; i < nthreads; i++) {
      if (thread_time[i] > max_time) max_time = thread_time[i];
      if (thread_time[i] < min_time) min_time = thread_time[i];
      avg_time += thread_time[i];
  }
  avg_time /= nthreads;
  double load_imbalance = (max_time - min_time) / avg_time;
  
  printf("Total time: %f\n", total_time);
  printf("Initialization time: %f\n", init_time);
  printf("Computation time: %f\n", sumU);
  printf("Computation/Total ratio: %.2f%%\n", (sumU/total_time)*100.0);
  printf("Other time (overhead): %f (%.2f%%)\n", 
    total_time - sumU, 
    ((total_time - sumU)/total_time)*100.0);
  printf("Load imbalance: %.2f\n", load_imbalance);

  //export_and_plot(time_spent_injecting, time_spent_updating, time_spent_simulation, Niterations);

  return 0;
}