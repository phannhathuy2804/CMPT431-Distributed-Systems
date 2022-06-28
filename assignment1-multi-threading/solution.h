#include <pthread.h>
#include <stdlib.h>
#include <iomanip>
#include "core/circular_queue.h"
#include "core/utils.h"


class Producer {
    public:
    int No_producers;
    int producer_id;
    long No_items;
    long No_items_of_each_type[3] = {0};
    long value_items_of_each_type[3] = {0};
    CircularQueue * production_buffer;
    pthread_mutex_t * mutex;
    pthread_cond_t * empty_cond;
    pthread_cond_t * full_cond;
    int * active_producer_count;
    int * active_consumer_count;
    double time_taken;
    int remainder;


    //Producer(CircularQueue * production_buffer, pthread_mutex_t * mutex, pthread_cond_t * empty_cond, pthread_cond_t * full_cond);
};

class Consumer {
    public:
    int consumer_id;
    long No_items = 0;
    long No_items_of_each_type[3] = {0};
    long value_items_of_each_type[3] = {0};
    CircularQueue * production_buffer;
    pthread_mutex_t * mutex;
    pthread_cond_t * empty_cond;
    pthread_cond_t * full_cond;
    int * active_producer_count;
    int * active_consumer_count;
    double time_taken;


};

class ProducerConsumerProblem {
    long n_items;
    int n_producers;
    int n_consumers;
    CircularQueue production_buffer;

    // Dynamic array of thread identifiers for producer and consumer threads.
    // Use these identifiers while creating the threads and joining the threads.
    pthread_t *producer_threads;
    pthread_t *consumer_threads;

    Producer *producers;
    Consumer *consumers;

    int active_producer_count;
    int active_consumer_count;

    // define any other members, mutexes, condition variables here
    pthread_mutex_t mutex;
    pthread_cond_t empty_cond;
    pthread_cond_t full_cond;

   public:
    // The following 6 methods should be defined in the implementation file (solution.cpp)
    ProducerConsumerProblem(long _n_items, int _n_producers, int _n_consumers,
                            long _queue_size);
    ~ProducerConsumerProblem();
    void startProducers();
    void startConsumers();
    void joinProducers();
    void joinConsumers();
    void printStats();

};
