#include "core/graph.h"
#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>


#define DEFAULT_STRATEGY "1"
#define DEFAULT_GRANULARITY "1"
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

typedef struct{
  double time_taken_of_thread;
  double time_waited_at_barrier1;
  double time_waited_at_barrier2;
  double time_getNextVetex;
  uintE edges_processed;
  uintV vertices_processed;
} OUTPUT_DATA;

//-------------------------------------------------------------------
struct DynamicCustomBarrier {
  int num_of_threads_;
  int current_waiting_;
  int barrier_call_;
  std::mutex my_mutex_;
  std::condition_variable my_cv_;

  DynamicCustomBarrier(int t_num_of_threads)
      : num_of_threads_(t_num_of_threads), current_waiting_(0),
        barrier_call_(0) {}

  void wait(std::atomic <uintV> * current_static_v) {
    std::unique_lock<std::mutex> u_lock(my_mutex_);
    int c = barrier_call_;
    current_waiting_++;
    if (current_waiting_ == num_of_threads_) {
      current_waiting_ = 0;
      *current_static_v = 0;
      
      // unlock and send signal to wake up
      barrier_call_++;
      u_lock.unlock();
      my_cv_.notify_all();
      return;
    }
    my_cv_.wait(u_lock, [&] { return (c != barrier_call_); });
    //  Condition has been reached. return
  }
};
//-------------------------------------------------------------------


// uintV getNextVertexToBeProcessed(std::atomic<uintV> * current_static_v, uint n){
//   uintV old;
//   uintV desired;
//   do
//   {
//     old = current_static_v->load();
//     desired = old + 1;
//     //std::cout<<"desired = "<< desired <<<<" tid= "<< tid <<"\n";
//     //std::cout<<"n = "<< n <<"\n";
//     if (desired > n)
//       return -1;
//   }
//   while(!current_static_v->compare_exchange_weak(old, desired));
//   return old;
// }

uintV getNextVertexToBeProcessed(std::atomic<uintV> * current_static_v, uint n, uint granularity){
  uintV val = current_static_v->fetch_add(granularity, std::memory_order_relaxed);
  if (val >= n)
    return -1;
  return val;
}

inline void DynamicPageRankCalculation(Graph* g, int max_iters, PageRankType *pr_curr, PageRankType *pr_next, std::vector<uintV> vertices, uint tid, OUTPUT_DATA* data, DynamicCustomBarrier* barrier, std::atomic <uintV> * current_static_v, uintV n, uint k ){
  timer t1;
  timer t_barrier1;
  timer t_barrier2;
  timer t_getNextVertex;
  double local_getNextVerteex_time = 0;
  double local_barrier1_time =0;
  double local_barrier2_time =0;
  uint local_processed_edges = 0;
  uint local_processed_vertices = 0;
  t1.start();
  
  for (int iter = 0; iter < max_iters; iter++) {
    // for each vertex 'v', process all its inNeighbors 'u'
    while (true) {
      t_getNextVertex.start();
      uintV v = getNextVertexToBeProcessed(current_static_v, n, k);
      local_getNextVerteex_time += t_getNextVertex.stop();
      if (v==-1)
        break;
      for(int j = 0 ; j < k; j++ ){
        uintE in_degree = g->vertices_[v].getInDegree();
        local_processed_edges += in_degree;
        for (uintE i = 0; i < in_degree; i++) {
          uintV u = g->vertices_[v].getInNeighbor(i);
          uintE u_out_degree = g->vertices_[u].getOutDegree();
          if (u_out_degree > 0)
              pr_next[v] += (pr_curr[u] / (PageRankType) u_out_degree);
        }
        v++;
        if(v >= n) break;
      }
    }
    t_barrier1.start();
    barrier->wait(current_static_v);
    local_barrier1_time += t_barrier1.stop();


    while(true) {
      t_getNextVertex.start();
      uintV v = getNextVertexToBeProcessed(current_static_v, n, k);
      local_getNextVerteex_time += t_getNextVertex.stop();
      if (v==-1)
        break;
      for(int j = 0 ; j < k; j++ ){
        local_processed_vertices+=1;
        pr_next[v] = PAGE_RANK(pr_next[v]);
        // reset pr_curr for the next iteration
        pr_curr[v] = pr_next[v];
        pr_next[v] = 0.0;
        v++;
        if(v >= n) break;
      }
    }
    t_barrier2.start();
    barrier->wait(current_static_v);
    local_barrier2_time += t_barrier2.stop();
  }
  data->time_getNextVetex = local_getNextVerteex_time;
  data->time_waited_at_barrier1 = local_barrier1_time;
  data->time_waited_at_barrier2 = local_barrier2_time;
  data->edges_processed = local_processed_edges;
  data->vertices_processed = local_processed_vertices;
  data->time_taken_of_thread = t1.stop();
  //*time_taken = t1.stop();
}


inline void StaticPageRankCalculation(Graph* g, int max_iters, PageRankType *pr_curr, PageRankType *pr_next, std::vector<uintV> vertices, uint tid, OUTPUT_DATA *data, CustomBarrier* barrier){
  timer t1;
  timer t_barrier1;
  timer t_barrier2;
  double local_barrier1_time =0;
  double local_barrier2_time =0;
  uint local_processed_edges = 0;
  uint local_processed_vertices = 0;
  t1.start();
  
  for (int iter = 0; iter < max_iters; iter++) {
    // for each vertex 'v', process all its inNeighbors 'u'
    for (uint v : vertices) {
      uintE in_degree = g->vertices_[v].getInDegree();
      local_processed_edges += in_degree;
      for (uintE i = 0; i < in_degree; i++) {
        uintV u = g->vertices_[v].getInNeighbor(i);
        uintE u_out_degree = g->vertices_[u].getOutDegree();
        if (u_out_degree > 0)
            pr_next[v] += (pr_curr[u] / (PageRankType) u_out_degree);

      }
    }
    t_barrier1.start();
    barrier->wait();
    local_barrier1_time += t_barrier1.stop();

    for (uint v : vertices) {
      local_processed_vertices+=1;
      pr_next[v] = PAGE_RANK(pr_next[v]);
      // reset pr_curr for the next iteration
      pr_curr[v] = pr_next[v];
      pr_next[v] = 0.0;
    }
    t_barrier2.start();
    barrier->wait();
    local_barrier2_time += t_barrier2.stop();


  }
  data->time_waited_at_barrier1 = local_barrier1_time;
  data->time_waited_at_barrier2 = local_barrier2_time;
  data->edges_processed = local_processed_edges;
  data->vertices_processed = local_processed_vertices;
  data->time_taken_of_thread = t1.stop();
  //*time_taken = t1.stop();
}

void pageRankParallel(Graph *g, int max_iters, uint n_threads , uint strategy, uint granularity) {
  uintV n = g->n_;
  uintE m = g->m_;

  PageRankType *pr_curr = new PageRankType[n];
  PageRankType *pr_next = new PageRankType[n];
  std::vector<uint> start(n_threads);
  std::vector<uint> end(n_threads);
  uint excess_vertices = n % n_threads;
  uint curr_vertex = 0;
  //ouble *time_taken_of_thread = new double[n_threads];
  OUTPUT_DATA * data = new OUTPUT_DATA[n_threads];

  std::vector <std::vector<uintV>>  work_balancer (n_threads);

  for (uintV i = 0; i < n; i++) {
    pr_curr[i] = INIT_PAGE_RANK;
    pr_next[i] = 0.0;
  }

  for (uint i = 0; i < n_threads; i++){
    data[i].time_taken_of_thread = 0;
    data[i].time_waited_at_barrier1 = 0;
    data[i].time_waited_at_barrier2 = 0;
    data[i].time_getNextVetex = 0;
    data[i].edges_processed = 0;
    data[i].vertices_processed = 0;
  }

  // Pull based pagerank
  timer t1;
  double time_taken = 0.0;

  CustomBarrier * barrier = new CustomBarrier(n_threads);
  DynamicCustomBarrier * dynamic_barrier = new DynamicCustomBarrier(n_threads);

  // Create threads and distribute the work across T threads
  // -------------------------------------------------------------------
  t1.start();
  std::thread *threads = new std::thread[n_threads];

  switch(strategy){
    case 1:
    {
      uint min_vertices_for_each_thread = n / n_threads;
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
        for (uintV j = start[i]; j <= end[i]; j++){
          work_balancer[i].push_back(j);
        }
      }
      for (int i = 0; i < n_threads; i++){
        threads[i] = std::thread(StaticPageRankCalculation, g, max_iters, pr_curr, pr_next, work_balancer[i], i, &data[i], barrier);
      }
      break; 
    }
    case 2:
    {
      uint total_assigned_edges = 0;
      uintV current_v = 0;
      for (uint i = 0; i < n_threads; i++){
        while (current_v < n){
          if (total_assigned_edges < (i + 1) * m / n_threads){
            total_assigned_edges +=  g->vertices_[current_v].getInDegree();
            work_balancer[i].push_back(current_v);
            current_v++;
          }
          else
            break;
        }
      }
      for (int i = 0; i < n_threads; i++){
        threads[i] = std::thread(StaticPageRankCalculation, g, max_iters, pr_curr, pr_next, work_balancer[i], i, &data[i], barrier);
      }
      break;
    }
    case 3: 
    {
      std::atomic <uintV> current_static_v{0};
      for (int i = 0; i < n_threads; i++){
        threads[i] = std::thread(DynamicPageRankCalculation, g, max_iters, std::ref(pr_curr), std::ref(pr_next), work_balancer[i], i, &data[i], dynamic_barrier, &current_static_v, n, 1);
      }
      break;
    }
    case 4:
    {
      std::atomic <uintV> current_static_v{0};
      for (int i = 0; i < n_threads; i++){
        threads[i] = std::thread(DynamicPageRankCalculation, g, max_iters, std::ref(pr_curr), std::ref(pr_next), work_balancer[i], i, &data[i], dynamic_barrier, &current_static_v, n, granularity);
      }
      break;
    }
  }


  for (uint i = 0; i < n_threads; i++){
  	threads[i].join();
  }

  
  time_taken = t1.stop();
  // -------------------------------------------------------------------
  std::cout << "thread_id, num_vertices, num_edges, barrier1_time, barrier2_time, getNextVertex_time, total_time" << std::endl;
  // Print the above statistics for each thread
  // Example output for 2 threads:
  // thread_id, time_taken
  // 0, 0.12
  // 1, 0.12
  for (uint i = 0; i < n_threads; i++){
    std::cout<< i <<", "<<data[i].vertices_processed<<", "<<data[i].edges_processed<<", "<<data[i].time_waited_at_barrier1<<", "<<data[i].time_waited_at_barrier2 << ", "<<data[i].time_getNextVetex  <<", "<< data[i].time_taken_of_thread<<"\n";
  }

  PageRankType sum_of_page_ranks = 0;
  for (uintV u = 0; u < n; u++) {
      //std::cout << pr_curr[u] << "\n";

    sum_of_page_ranks += pr_curr[u];
  }
  std::cout << "Sum of page ranks : " << sum_of_page_ranks << "\n";
  std::cout << "Time taken (in seconds) : " << time_taken << "\n";
  delete[] pr_curr;
  delete[] pr_next;
}

int main(int argc, char *argv[]) {
  cxxopts::Options options(
      "page_rank_pull",
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
          {"strategy", "Strategy id", 
           cxxopts::value<uint>()->default_value(DEFAULT_STRATEGY)},
           {"granularity", "Granularity", 
           cxxopts::value<uint>()->default_value(DEFAULT_GRANULARITY)},
      });

  auto cl_options = options.parse(argc, argv);
  uint n_threads = cl_options["nThreads"].as<uint>();
  uint max_iterations = cl_options["nIterations"].as<uint>();
  std::string input_file_path = cl_options["inputFile"].as<std::string>();
  uint strategy = cl_options["strategy"].as<uint>();
  if (strategy < 0 || strategy > 4){
    std::cout << "Invalid strategy" << std::endl;
    return 0;
  }
  uint granularity = cl_options["granularity"].as<uint>();
  if (granularity <= 0){
    std::cout << "Invalid granularity" << std::endl;
    return 0;
  }

#ifdef USE_INT
  std::cout << "Using INT\n";
#else
  std::cout << "Using DOUBLE\n";
#endif
  std::cout << std::fixed;
  std::cout << "Number of Threads : " << n_threads << std::endl;
  std::cout << "Strategy : " << strategy << std::endl;
  std::cout << "Granularity : " << granularity << std::endl;
  std::cout << "Iterations : " << max_iterations << std::endl;

  Graph g;
  std::cout << "Reading graph\n";
  g.readGraphFromBinary<int>(input_file_path);
  std::cout << "Created graph\n";
  pageRankParallel(&g, max_iterations, n_threads, strategy, granularity);

  return 0;
}
