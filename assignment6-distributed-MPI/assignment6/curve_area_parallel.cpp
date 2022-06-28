#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <mpi.h>


#define sqr(x) ((x) * (x))
#define DEFAULT_NUMBER_OF_POINTS "1000000000"
#define DEFAULT_A "2"
#define DEFAULT_B "1"
#define DEFAULT_RANDOM_SEED "1"

uint c_const = (uint)RAND_MAX + (uint)1;
inline double get_random_coordinate(uint *random_seed) {
  return ((double)rand_r(random_seed)) / c_const;  // thread-safe random number generator
}

unsigned long get_points_in_curve(unsigned long n, uint random_seed, float a, float b) {
  unsigned long local_curve_count = 0;
  double x_coord, y_coord;
  for (unsigned long i = 0; i < n; i++) {
    x_coord = ((2.0 * get_random_coordinate(&random_seed)) - 1.0);
    y_coord = ((2.0 * get_random_coordinate(&random_seed)) - 1.0);
    if ((a*sqr(x_coord) + b*sqr(sqr(y_coord))) <= 1.0)
      local_curve_count++;
  }
  return local_curve_count;
}

void curve_area_calculation_parallel(unsigned long n, float a, float b, uint r_seed) {
  timer serial_timer;
  double time_taken = 0.0;
  uint random_seed = r_seed;

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  r_seed+=world_rank;
  uint points_to_be_generated;
  uint min_points_per_process = n / world_size;
  uint excess_points = n % world_size;
  if (world_rank < excess_points)
    points_to_be_generated = min_points_per_process + 1;
  else
    points_to_be_generated = min_points_per_process;   

  serial_timer.start();
 
  unsigned long local_curve_points = get_points_in_curve(points_to_be_generated, r_seed, a, b);

  unsigned long global_curve_points;
  MPI_Reduce(&local_curve_points, &global_curve_points, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
  double area_value;
  if(world_rank==0){
    area_value = 4.0 * (double)global_curve_points / (double)n;
  }
  time_taken = serial_timer.stop();
  double * time_takens = (double*)malloc(sizeof(double) * world_size);
  uint * points_generated = (uint*)malloc(sizeof(uint)* world_size);
  unsigned long * curve_points_arr = (unsigned long*)malloc(sizeof(unsigned long) * world_size);
  //*------------------------------------------------------------------------
  if(world_rank==0){

    MPI_Gather(&time_taken, 1, MPI_DOUBLE, time_takens, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&points_to_be_generated, 1, MPI_UNSIGNED, points_generated, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Gather(&local_curve_points, 1, MPI_UNSIGNED_LONG, curve_points_arr, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  } else{
    MPI_Gather(&time_taken, 1, MPI_DOUBLE, NULL, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Gather(&points_to_be_generated, 1, MPI_UNSIGNED, NULL, 1, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
    MPI_Gather(&local_curve_points, 1, MPI_UNSIGNED_LONG, NULL, 1, MPI_UNSIGNED_LONG, 0, MPI_COMM_WORLD);
  }
  
  if (world_rank==0){
    std::cout << "rank, points_generated, curve_points, time_taken\n";
    for( int i = 0; i < world_size; i++){
      std::cout << i << ", " << points_generated[i] << ", "
                  << curve_points_arr[i] << ", " << std::setprecision(TIME_PRECISION)
                  << time_takens[i] << "\n";
    }

    std::cout << "Total points generated : " << n << "\n";
    std::cout << "Total points in curve : " << global_curve_points << "\n";
    std::cout << "Area : " << std::setprecision(VAL_PRECISION) << area_value
              << "\n";
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
              << time_taken << "\n";
  }
  free(time_takens);
  free(points_generated);
  free(curve_points_arr);
}

int main(int argc, char *argv[]) {
  // Initialize command line arguments
  cxxopts::Options options("Curve_area_calculation",
                           "Calculate area inside curve a x^2 + b y ^4 = 1 using serial and parallel execution");
  options.add_options(
      "custom",
      {
          {"nPoints", "Number of points",         
           cxxopts::value<unsigned long>()->default_value(DEFAULT_NUMBER_OF_POINTS)},
	        {"coeffA", "Coefficient a",
	         cxxopts::value<float>()->default_value(DEFAULT_A)},
          {"coeffB", "Coefficient b",
           cxxopts::value<float>()->default_value(DEFAULT_B)},
          {"rSeed", "Random Seed",
           cxxopts::value<uint>()->default_value(DEFAULT_RANDOM_SEED)}
      });
  auto cl_options = options.parse(argc, argv);
  unsigned long n_points = cl_options["nPoints"].as<unsigned long>();
  float a = cl_options["coeffA"].as<float>();
  float b = cl_options["coeffB"].as<float>();
  uint r_seed = cl_options["rSeed"].as<uint>();
  MPI_Init(NULL, NULL);
    // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  if(world_rank==0){
    std::cout << "Number of points : " << n_points << "\n";;
    std::cout << "A : " << a << "\n" << "B : " << b << "\n";
    std::cout << "Random Seed : " << r_seed << "\n";
  }

  curve_area_calculation_parallel(n_points, a, b, r_seed);
  MPI_Finalize();
  return 0;
}
