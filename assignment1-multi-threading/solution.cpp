#include "solution.h"
#include <unistd.h>

void *producerFunction(void *_arg) {
    // Parse the _arg passed to the function.
    // Enqueue `n` items into the `production_buffer`. The items produced should
    // be 0, 1, 2,..., (n-1).
    // Keep track of the number of items produced, their type
    // and the value produced by the thread
    // The producer that was last active should ensure that all the consumers have
    // finished. NOTE: Each thread will enqueue `n` items.
    // Use mutex variables and conditional variables as necessary.

    // Each producer enqueues `np` items where np=n/nProducers except for producer 0
    timer t1;
    t1.start();
    Producer * producer = (Producer*) _arg;
    long np = producer->No_items;

    long item = producer->producer_id;
    long items_produced = 0;
    long value_produced = 0;
    int i = 0;
    int type = 0;
    CircularQueue * production_buffer = producer->production_buffer;
    long No_items_of_each_type[3];

    for (int j = 0; j < 3; j++){
        No_items_of_each_type[j] = producer->No_items_of_each_type[j];
    }
    
    while (items_produced < np) {
        pthread_mutex_lock(producer->mutex);
        item = producer->producer_id + i * producer->No_producers;
        bool ret = production_buffer->enqueue(item, type, producer->producer_id);  //item includes value and type
        if (ret == true) {
            if (production_buffer->itemCount() == 1) {
                // The queue is no longer empty
                // Signal all consumers indicating queue is not empty
                pthread_cond_broadcast(producer->empty_cond);
            }
            pthread_mutex_unlock(producer->mutex);
            //value_produced += item;
            items_produced++;
            item = producer->producer_id + i * producer->No_producers;
            producer->value_items_of_each_type[type]+=item;
            i++;
            No_items_of_each_type[type]-=1;
            if (No_items_of_each_type[type]==0){
                value_produced = 0;
                type+=1;
            }
            // also update item count and value for item type produced
        } else {
            // production_buffer is full, so block on conditional variable waiting for consumer to signal.
            if(production_buffer->isFull()){
                pthread_cond_wait(producer->full_cond, producer->mutex);
                pthread_mutex_unlock(producer->mutex);
            }
        }
    }
    // After production is completed:
    // Update the number of producers that are currently active.
    //--active_producer_count;
    pthread_mutex_lock(producer->mutex);
    *producer->active_producer_count-=1;
    pthread_mutex_unlock(producer->mutex);
    // The producer that was last active (can be determined using `active_producer_count`) will keep signalling the consumers until all consumers have finished (can be determined using `active_consumer_count`).
    if(*producer->active_producer_count==0){
        while (*producer->active_consumer_count){
            pthread_cond_broadcast(producer->empty_cond);

        }
    }
    producer->time_taken  = t1.stop();
    pthread_exit(producer);
}

void *consumerFunction(void *_arg) {
    // Parse the _arg passed to the function.
    // The consumer thread will consume items by dequeueing the items from the
    // `production_buffer`.
    // Keep track of the number of items consumed and their value and type
    // Once the productions is complete and the queue is also empty, the thread
    // will exit. NOTE: The number of items consumed by each thread need not be the same
    // Use mutex variables and conditional variables as necessary.
    // Each consumer dequeues items from the `production_buffer`
    timer t1;
    t1.start();
    Consumer * consumer = (Consumer*) _arg;

    long np = consumer->No_items;

    long item;
    long items_consumed[3] = {0};
    long value_consumed[3] = {0};
    int type;
    int source;
    CircularQueue * production_buffer = consumer->production_buffer;


    while (true) {
        int ret_mutex = pthread_mutex_lock(consumer->mutex);
        bool ret = production_buffer->dequeue(&item, &type, &source);  // item includes value and type
        if (ret == true) {
            if (production_buffer->itemCount() == production_buffer->getCapacity() - 1) {
                // The queue is no longer full
                // Signal all producers indicating queue is not full
                pthread_cond_broadcast(consumer->full_cond);
            }
            pthread_mutex_unlock(consumer->mutex);
            consumer->value_items_of_each_type[type] += item;
            consumer->No_items_of_each_type[type]+=1;
            // also update value and count for consumed item type
        } else {
            // production_buffer is empty, so block on conditional variable waiting for producer to signal.
            pthread_cond_wait(consumer->empty_cond, consumer->mutex);
            pthread_mutex_unlock(consumer->mutex);
            // The thread can wake up because of 2 scenarios:
            // Scenario 1: There are no more active producers (i.e., production is complete) and the queue is empty. This is the exit condition for consumers, and at this point consumers should decrement `active_consumer_count`.
            if (!(*consumer->active_producer_count) && production_buffer->isEmpty()){
                pthread_mutex_lock(consumer->mutex);
                *consumer->active_consumer_count-=1;
                pthread_mutex_unlock(consumer->mutex);
                break;
            }

            // Scenario 2: The queue is not empty and/or the producers are active. Continue consuming.
        }
    }
    consumer->time_taken = t1.stop();
    pthread_exit(consumer);
}

ProducerConsumerProblem::ProducerConsumerProblem(long _n_items,
                                                 int _n_producers,
                                                 int _n_consumers,
                                                 long _queue_size)
    : n_items(_n_items), n_producers(_n_producers), n_consumers(_n_consumers), production_buffer(_queue_size) {
    std::cout << "Constructor\n";
    std::cout << "Number of producers: " << n_producers << "\n";
    std::cout << "Number of consumers: " << n_consumers << "\n";

    if (n_consumers) { 
        consumers = new Consumer[n_consumers];
        consumer_threads = new pthread_t[n_consumers];
    }
    if (n_producers) {
        producers = new Producer[n_producers];
        producer_threads = new pthread_t[n_producers];
    }

    // Initialize all mutex and conditional variables here.

    pthread_mutex_init(&mutex, NULL);
    pthread_cond_init(&full_cond, NULL);
    pthread_cond_init(&empty_cond, NULL);
}

ProducerConsumerProblem::~ProducerConsumerProblem() {
    std::cout << "Destructor\n";
    if (n_producers) {
        delete[] producers;
        delete[] producer_threads;
    }
    if (n_consumers) {
        delete[] consumers;
        delete[] consumer_threads;
    }
    // Destroy all mutex and conditional variables here.
}

void ProducerConsumerProblem::startProducers() {
    std::cout << "Starting Producers\n";
    active_producer_count = n_producers;
    // Compute number of items for each thread, and number of items per type per thread
    if (producers){
        producers[0].producer_id = 0;
        producers[0].No_items =  (n_items / n_producers) + (n_items % n_producers);
        producers[0].No_items_of_each_type[0] = producers[0].No_items/2;
        producers[0].No_items_of_each_type[1] = producers[0].No_items/3;
        producers[0].No_items_of_each_type[2] = producers[0].No_items - producers[0].No_items_of_each_type[0] - producers[0].No_items_of_each_type[1];
        producers[0].production_buffer = &production_buffer;
        producers[0].mutex = &this->mutex;
        producers[0].full_cond = &this->full_cond;
        producers[0].empty_cond = &this->empty_cond;
        producers[0].active_producer_count = &active_producer_count;
        producers[0].active_consumer_count = &active_consumer_count;
        producers[0].No_producers = n_producers;

        for(int i = 1; i < n_producers; i++){
            producers[i].producer_id = i;
            producers[i].No_items =  (n_items / n_producers);
            producers[i].No_items_of_each_type[0] = producers[i].No_items/2;
            producers[i].No_items_of_each_type[1] = producers[i].No_items/3;
            producers[i].No_items_of_each_type[2] = producers[i].No_items - producers[i].No_items_of_each_type[0] - producers[i].No_items_of_each_type[1];
            producers[i].production_buffer = &production_buffer;
            producers[i].mutex = &this->mutex;
            producers[i].full_cond = &this->full_cond;
            producers[i].empty_cond = &this->empty_cond;
            producers[i].active_producer_count = &active_producer_count;
            producers[i].active_consumer_count = &active_consumer_count;
            producers[i].No_producers = n_producers;

       }
    }
    // Create producer threads P1, P2, P3,.. using pthread_create.
    for (int i = 0; i < n_producers; i++){
        pthread_create(&producer_threads[i], NULL, producerFunction, &producers[i]);
    }
}

void ProducerConsumerProblem::startConsumers() {
    std::cout << "Starting Consumers\n";
    active_consumer_count = n_consumers;
    // Create consumer threads C1, C2, C3,.. using pthread_create.
    for (int i = 0; i < n_consumers; i++){
        consumers[i].consumer_id = i;
        consumers[i].production_buffer = &production_buffer;
        consumers[i].mutex = &this->mutex;
        consumers[i].full_cond = &this->full_cond;
        consumers[i].empty_cond = &this->empty_cond;
        consumers[i].active_producer_count = &active_producer_count;
        consumers[i].active_consumer_count = &active_consumer_count;

        pthread_create(&consumer_threads[i], NULL, consumerFunction, &consumers[i]);
    }
}

void ProducerConsumerProblem::joinProducers() {
    std::cout << "Joining Producers\n";
    // Join the producer threads with the main thread using pthread_join
    void * ret;
    for (int i = 0; i < n_producers; i++){
        int ret = pthread_join(producer_threads[i], NULL);
    }
}

void ProducerConsumerProblem::joinConsumers() {
    std::cout << "Joining Consumers\n";
    // Join the consumer threads with the main thread using pthread_join
    for (int i = 0; i < n_consumers; i++){
        int ret = pthread_join(consumer_threads[i], NULL);
    }
}

void ProducerConsumerProblem::printStats() {
    std::cout << "Producer stats\n";
    std::cout << "producer_id, items_produced_type0:value_type0, items_produced_type1:value_type1, items_produced_type2:value_type2, total_value_produced, time_taken\n";
    
    // Make sure you print the producer stats in the following manner
    //  0, 125000:31249750000, 83333:55555111112, 41667:38194638888, 124999500000, 0.973188
    //  1, 125000:31249875000, 83333:55555194445, 41667:38194680555, 124999750000, 1.0039
    //  2, 125000:31250000000, 83333:55555277778, 41667:38194722222, 125000000000, 1.02925
    //  3, 125000:31250125000, 83333:55555361111, 41667:38194763889, 125000250000, 0.999188

    long total_produced[3] = {0};        // total produced per type
    long total_value_produced[3] = {0};  // total value produced per type

    for (int i = 0; i < n_producers; i++){
        long total_value_of_thread = producers[i].value_items_of_each_type[0] + producers[i].value_items_of_each_type[1] + producers[i].value_items_of_each_type[2];
        std::cout << producers[i].producer_id << ", " <<producers[i].No_items_of_each_type[0] << ":" << producers[i].value_items_of_each_type[0]<<", " <<producers[i].No_items_of_each_type[1] << ":" << producers[i].value_items_of_each_type[1]<<", " <<producers[i].No_items_of_each_type[2] << ":" << producers[i].value_items_of_each_type[2]<< ", " <<total_value_of_thread<<", "<<producers->time_taken<<std::endl;
        
        total_value_produced[0] += producers[i].value_items_of_each_type[0];
        total_value_produced[1] += producers[i].value_items_of_each_type[1];
        total_value_produced[2] += producers[i].value_items_of_each_type[2];

        total_produced[0] += producers[i].No_items_of_each_type[0];
        total_produced[1] += producers[i].No_items_of_each_type[1];
        total_produced[2] += producers[i].No_items_of_each_type[2];

    }

    std::cout << "Total produced = " << total_produced[0] + total_produced[1] + total_produced[2] << "\n";
    std::cout << "Total value produced = " << total_value_produced[0] + total_value_produced[1] + total_value_produced[2] << "\n";
    std::cout << "Consumer stats\n";
    std::cout << "consumer_id, items_consumed_type0:value_type0, items_consumed_type1:value_type1, items_consumed_type2:value_type2, time_taken\n";

    // Make sure you print the consumer stats in the following manner
    // 0, 256488:63656791749, 163534:109699063438, 87398:79885550318, 1.02899
    // 1, 243512:61342958251, 169798:112521881008, 79270:72893255236, 1.02891

 long total_consumed[3] = {0};        // total consumed per type
    long total_value_consumed[3] = {0};  // total value consumed per type

    for (int i = 0; i < n_consumers; i++){
        long total_value_of_thread = consumers[i].value_items_of_each_type[0] + consumers[i].value_items_of_each_type[1] + consumers[i].value_items_of_each_type[2];
        std::cout << consumers[i].consumer_id << ", " <<consumers[i].No_items_of_each_type[0] << ":" << consumers[i].value_items_of_each_type[0]<<", " <<consumers[i].No_items_of_each_type[1] << ":" << consumers[i].value_items_of_each_type[1]<<", " <<consumers[i].No_items_of_each_type[2] << ":" << consumers[i].value_items_of_each_type[2]<< ", " <<consumers[i].time_taken<<std::endl;
        
        total_value_consumed[0] += consumers[i].value_items_of_each_type[0];
        total_value_consumed[1] += consumers[i].value_items_of_each_type[1];
        total_value_consumed[2] += consumers[i].value_items_of_each_type[2];

        total_consumed[0] += consumers[i].No_items_of_each_type[0];
        total_consumed[1] += consumers[i].No_items_of_each_type[1];
        total_consumed[2] += consumers[i].No_items_of_each_type[2];

    }

    std::cout << "Total consumed = " << total_consumed[0] + total_consumed[1] + total_consumed[2] << "\n";
    std::cout << "Total value consumed = " << total_value_consumed[0] + total_value_consumed[1] + total_value_consumed[2] << "\n";
    std::cout <<(long) production_buffer.itemCount();
}
