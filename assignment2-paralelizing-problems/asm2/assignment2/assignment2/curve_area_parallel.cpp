#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>

#define sqr(x) ((x) * (x))
#define DEFAULT_NUMBER_OF_THREADS "1"
#define DEFAULT_NUMBER_OF_POINTS "1000000000"
#define DEFAULT_A "2"
#define DEFAULT_B "1"
#define DEFAULT_RANDOM_SEED "1"

uint c_const = (uint)RAND_MAX + (uint)1;
inline double get_random_coordinate(uint *random_seed) {
  return ((double)rand_r(random_seed)) / c_const;  // thread-safe random number generator
}

void get_points_in_curve(unsigned long n, uint random_seed, float a, float b, double* time_taken, unsigned long * curve_count) {
  timer serial_timer;
  *curve_count = 0;
  unsigned long curve_point_count = 0;
  double x_coord, y_coord;
  serial_timer.start();
  for (unsigned long i = 0; i < n; i++) {
    x_coord = ((2.0 * get_random_coordinate(&random_seed)) - 1.0);
    y_coord = ((2.0 * get_random_coordinate(&random_seed)) - 1.0);
    if ((a*sqr(x_coord) + b*sqr(sqr(y_coord))) <= 1.0)
      curve_point_count++;
  }
  double timer_taken = serial_timer.stop();
  //std::cout<<timer_taken<<"\n";
  *curve_count = curve_point_count;
  *time_taken = timer_taken;

}

void curve_area_calculation_parallel(uint n_threads,unsigned long n, float a, float b, uint r_seed) {
  timer serial_timer;
  double time_taken = 0.0;
  uint random_seed = r_seed;
  unsigned long point_per_thread = n / n_threads;
  double *time_taken_of_thread = new double[n_threads];
  unsigned long *curve_points = new unsigned long [n_threads];
  std::thread *threads = new std::thread[n_threads];
  unsigned long total_curve_point = 0;

  serial_timer.start();
 
  //unsigned long curve_points = get_points_in_curve(n, r_seed, a, b);
  threads[0] = std::thread(get_points_in_curve, point_per_thread + (n % n_threads), r_seed, a, b, &time_taken_of_thread[0], &curve_points[0]);
  for (int i = 1; i < n_threads; i++){
  	threads[i] = std::thread(get_points_in_curve, point_per_thread, r_seed, a, b, &time_taken_of_thread[i], &curve_points[i]);
  }
  for (int i = 0; i < n_threads; i++){
  	threads[i].join();
  	total_curve_point+=curve_points[i];
  }
  time_taken = serial_timer.stop();

  
  double area_value = 4.0 * (double)total_curve_point / (double)n;


  std::cout << "thread_id, points_generated, curve_points, time_taken\n";
	  std::cout << 0 << ", "<< point_per_thread+ (n % n_threads) << ", "
                << 0  << curve_points[0] << ", " << std::setprecision(TIME_PRECISION)
              << time_taken_of_thread[0] << "\n";
  for (int i = 1; i < n_threads; i++){
  	  std::cout << i << ", "<< point_per_thread << ", "
   		<< i  << curve_points[i] << ", " << std::setprecision(TIME_PRECISION)
              << time_taken_of_thread[i] << "\n";
  }

  //*------------------------------------------------------------------------

  std::cout << "Total points generated : " <<n << "\n";
  std::cout << "Total points in curve : " << total_curve_point << "\n";
  std::cout << "Area : " << std::setprecision(VAL_PRECISION) << area_value
            << "\n";
  std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
            << time_taken << "\n";
}

int main(int argc, char *argv[]) {
  // Initialize command line arguments
  cxxopts::Options options("Curve_area_calculation",
                           "Calculate area inside curve a x^2 + b y ^4 = 1 using serial and parallel execution");
  options.add_options(
      "custom",
      {   {"nThreads", "Number of threads",
  		   cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_THREADS)},
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
  uint n_threads = cl_options["nThreads"].as<uint>();
  unsigned long n_points = cl_options["nPoints"].as<unsigned long>();
  float a = cl_options["coeffA"].as<float>();
  float b = cl_options["coeffB"].as<float>();
  uint r_seed = cl_options["rSeed"].as<uint>();
  std::cout << "Number of points : " << n_points << "\n";
  std::cout << "Number of threads : " << n_threads << "\n";
  std::cout << "A : " << a << "\n" << "B : " << b << "\n";
  std::cout << "Random Seed : " << r_seed << "\n";

  curve_area_calculation_parallel(n_threads, n_points, a, b, r_seed);
  return 0;
}
