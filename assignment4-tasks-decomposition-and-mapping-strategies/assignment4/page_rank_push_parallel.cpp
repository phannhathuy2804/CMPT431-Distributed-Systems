#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>
#include <mutex>

#ifdef USE_INT
#define INIT_PAGE_RANK 100000
#define EPSILON 1000
#define PAGE_RANK(x) (15000 + (5 * x) / 6)
#define CHANGE_IN_PAGE_RANK(x, y) std::abs(x - y)
typedef int64_t PageRankType;
#else
#define INIT_PAGE_RANK 1.0
#define EPSILON 0.01
#define DAMPING 0.85
#define PAGE_RANK(x) (1 - DAMPING + DAMPING * x)
#define CHANGE_IN_PAGE_RANK(x, y) std::fabs(x - y)
typedef double PageRankType;
#endif

//std::mutex mtx; 


typedef struct {
  PageRankType value;
  std::mutex lock;
} PR_NEXT;

typedef struct{
  Graph *g;
  int max_iters;
  int start;
  int end;
  double * time_taken;
  PageRankType *pr_curr;
  PR_NEXT * pr_next;
  CustomBarrier * barrier;
} DATA;




inline void pageRankCalculation(DATA data ){
  timer t1;
  t1.start();

  for (int iter = 0; iter < data.max_iters; iter++) {

  //std::unique_lock<std::mutex> u_lock(mtx);
    // for each vertex 'u', process all its outNeighbors 'v'
    for (uintV u = data.start; u <= data.end; u++) {
      uintE out_degree = data.g->vertices_[u].getOutDegree();
      for (uintE i = 0; i < out_degree; i++) {
        uintV v = data.g->vertices_[u].getOutNeighbor(i);
        bool ret= true;
        do{
          ret = data.pr_next[v].lock.try_lock();
        } while (!ret);
        data.pr_next[v].value += (data.pr_curr[u] / (PageRankType) out_degree);
        data.pr_next[v].lock.unlock();
      }
    }
    data.barrier->wait();
    for (uintV v = data.start; v <= data.end; v++) {
      PageRankType temp = data.pr_next[v].value;
      data.pr_next[v].value = PAGE_RANK(temp);

      // reset pr_curr for the next iteration
      data.pr_curr[v] = data.pr_next[v].value;
      data.pr_next[v].value = 0.0;
    }
    data.barrier->wait();
  }
  *data.time_taken = t1.stop();

}

void pageRankParallel(Graph *g, int max_iters, uint n_threads) {
  uintV n = g->n_;

  PageRankType *pr_curr = new PageRankType[n];
  PageRankType *pr_next = new PageRankType[n];
  std::vector<uint> start(n_threads);
  std::vector<uint> end(n_threads);
  uint min_vertices_for_each_thread = n /   n_threads;
  uint excess_vertices = n % n_threads;
  uint curr_vertex = 0;
  double *time_taken_of_thread = new double[n_threads];
  std::mutex *mtx = new std::mutex[n];
  PR_NEXT * pr_next1 = new PR_NEXT[n];
  DATA * data = new DATA[n_threads];





  for (uintV i = 0; i < n; i++) {
    pr_curr[i] = INIT_PAGE_RANK;
    pr_next[i] = 0.0;
    pr_next1[i].value = 0.0;
  }

  // Push based pagerank
  timer t1;
  double time_taken = 0.0;

  CustomBarrier * barrier = new CustomBarrier(n_threads);

  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  t1.start();

  for (uint i = 0; i < n_threads; i++) {
    start[i] = curr_vertex;
    if (excess_vertices > 0) {
      end[i] = curr_vertex + min_vertices_for_each_thread;
      excess_vertices--;
      } 
    else {
        end[i] = curr_vertex + min_vertices_for_each_thread - 1;
      }
    curr_vertex = end[i]+1;
  } 
  for (int i=0 ; i < n_threads; i++){
    data[i].g = g;
    data[i].max_iters = max_iters;
    data[i].start = start[i];
    data[i].end = end[i];
    data[i].time_taken = &time_taken_of_thread[i];
    data[i].pr_curr = pr_curr;
    data[i].pr_next = pr_next1;
    data[i].barrier = barrier;
  }
  //std::thread *threads = new std::thread[n_threads];
  //std::cout<<"size of prcurr"<<sizeof(pr_curr)<<"\n";
  //std::cout<<"size of data"<<sizeof(data[1])<<"\n";

  std::vector<std::thread> threads (n_threads);
  for (int i = 0; i < n_threads; i++){
      threads[i] = std::thread(pageRankCalculation, data[i]);
  }
  
  for (uint i = 0; i < n_threads; i++){
  	threads[i].join();
  }
  

  
  time_taken = t1.stop();
  // -------------------------------------------------------------------
  std::cout << "thread_id, time_taken" << std::endl;
  //std::cout << "0, " << time_taken << std::endl;
  // Print the above statistics for each thread
  // Example output for 2 threads:
  // thread_id, time_taken
  // 0, 0.12
  // 1, 0.12
  for (uint i = 0; i < n_threads; i++){
    std::cout << i << ", " << std::setprecision(TIME_PRECISION) << time_taken_of_thread[i] <<"\n";
  }

  PageRankType sum_of_page_ranks = 0;
  for (uintV u = 0; u < n; u++) {
    sum_of_page_ranks += pr_curr[u];
  }
  std::cout << "Sum of page ranks : " << sum_of_page_ranks << "\n";
  std::cout << "Time taken (in seconds) : " << time_taken << "\n";
  delete[] pr_curr;
  delete[] pr_next;
}

int main(int argc, char *argv[]) {
  cxxopts::Options options(
      "page_rank_push",
      "Calculate page_rank using serial and parallel execution");
  options.add_options(
      "",
      {
          {"nThreads", "Number of Threads",
           cxxopts::value<uint>()->default_value(DEFAULT_NUMBER_OF_THREADS)},
          {"nIterations", "Maximum number of iterations",
           cxxopts::value<uint>()->default_value(DEFAULT_MAX_ITER)},
          {"inputFile", "Input graph file path",
           cxxopts::value<std::string>()->default_value(
               "/scratch/input_graphs/roadNet-CA")},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_threads = cl_options["nThreads"].as<uint>();
  uint max_iterations = cl_options["nIterations"].as<uint>();
  std::string input_file_path = cl_options["inputFile"].as<std::string>();

#ifdef USE_INT
  std::cout << "Using INT" << std::endl;
#else
  std::cout << "Using DOUBLE" << std::endl;
#endif
  std::cout << std::fixed;
  std::cout << "Number of Threads : " << n_threads << std::endl;
  std::cout << "Number of Iterations: " << max_iterations << std::endl;

  Graph g;
  std::cout << "Reading graph\n";
  g.readGraphFromBinary<int>(input_file_path);
  std::cout << "Created graph\n";
  pageRankParallel(&g, max_iterations, n_threads);

  return 0;
}
