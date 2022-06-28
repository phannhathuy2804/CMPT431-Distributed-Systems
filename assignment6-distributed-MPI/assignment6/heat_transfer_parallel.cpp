#include "core/utils.h"
#include <iomanip>
#include <iostream>
#include <stdlib.h>
#include <thread>
#include <mpi.h>

#define DEFAULT_GRID_SIZE "1000"
#define DEFAULT_CX "1"
#define DEFAULT_CY "1"
#define DEFAULT_TIME_STEPS "1000"
#define DEFAULT_MIDDLE_TEMP "600"

class TemperatureArray {
private:
	uint size;
	uint step;
	double Cx;
	double Cy;
	double *CurrArray;
	double *PrevArray;
  void assign(double *A, uint x, uint y, double newvalue) {
    A[x*size+y] = newvalue;
  };
  double read(double *A, uint x, uint y) {
    return A[x*size+y];
  }; 
public:	
	TemperatureArray(uint input_size, double iCx, double iCy, double init_temp) { // create array of dimension sizexsize
		size = input_size;
 		Cx = iCx;
    		Cy = iCy;
		step = 0;
		CurrArray = (double *)malloc(size*size*sizeof(double));
		PrevArray = (double *)malloc(size*size*sizeof(double));
		for (uint i = 0; i < size; i++)
			for (uint j = 0; j < size; j++) {
				if ((i > size/3) && (i < 2*size/3) && (j > size/3) && (j < 2*size/3)) {
					assign(PrevArray, i, j, init_temp); assign(CurrArray, i, j, init_temp);
				}
				else {
					assign(PrevArray, i, j, 0); assign (CurrArray, i, j, 0);
				}	
			}
	};
 
	~TemperatureArray() {
		free (PrevArray);   free (CurrArray);
	};

	void IncrementStepCount() { step ++; };

	uint ReadStepCount() { return(step); };

	void ComputeNewTemp(uint x, uint y) {
		if ((x > 0) && (x < size-1) && (y > 0) && (y < size-1))
			assign(CurrArray, x, y , read(PrevArray,x,y)	
				+ Cx * (read(PrevArray, x-1, y) + read(PrevArray, x+1, y) - 2*read(PrevArray, x, y)) 
				+ Cy * (read(PrevArray, x, y-1) + read(PrevArray, x, y+1) - 2*read(PrevArray, x, y)));
	};

	void SwapArrays() {
		double *temp = PrevArray;
		PrevArray = CurrArray;
		CurrArray = temp;
	};	
 
	double temp(uint x, uint y) {
		return read(CurrArray, x, y);
	};

	void write(uint x, uint y, double newvalue) {
		assign(CurrArray, x, y, newvalue);
	};
};

inline void heat_transfer_calculation(uint size, uint start, uint end, double *time_taken, TemperatureArray* T, uint steps) {
  timer t1;
    // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  t1.start();
  uint stepcount;
  for (stepcount = 1; stepcount <= steps; stepcount ++) {
	  for (uint x = start; x <= end; x++) {
		  for (uint y = 0; y < size; y++) {
			  T->ComputeNewTemp(x, y);
		  }
	  }

    double * boundaryColumn = new double[size];

    if(world_rank % 2 ==0){
      if(world_rank < world_size -1){
        for(int i =0; i < size; i++){
          boundaryColumn[i] = T->temp(end,i);
        }
        MPI_Send(boundaryColumn, size, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD);
        MPI_Recv(boundaryColumn, size, MPI_DOUBLE, world_rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for(int i =0; i < size; i++){
          T->write(end+1, i, boundaryColumn[i]);
        }
      }
      
      if(world_rank > 0){
        MPI_Recv(boundaryColumn, size, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for(int i =0; i < size; i++){
          T->write(start-1, i, boundaryColumn[i]);
        }
        for(int i =0; i < size; i++){
          boundaryColumn[i] = T->temp(start,i);
        }
        MPI_Send(boundaryColumn, size, MPI_DOUBLE, world_rank-1, 1, MPI_COMM_WORLD);

      }
    } else {
      if(world_rank > 0){
        MPI_Recv(boundaryColumn, size, MPI_DOUBLE, world_rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for(int i =0; i < size; i++){
          T->write(start-1, i, boundaryColumn[i]);
        }
        for(int i =0; i < size; i++){
          boundaryColumn[i] = T->temp(start,i);
        }
        MPI_Send(boundaryColumn, size, MPI_DOUBLE, world_rank-1, 1, MPI_COMM_WORLD);

      }

      if(world_rank < world_size -1){
        for(int i =0; i < size; i++){
          boundaryColumn[i] = T->temp(end,i);
        }
        MPI_Send(boundaryColumn, size, MPI_DOUBLE, world_rank+1, 0, MPI_COMM_WORLD);
        MPI_Recv(boundaryColumn, size, MPI_DOUBLE, world_rank+1, 1, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        for(int i =0; i < size; i++){
          T->write(end+1, i, boundaryColumn[i]);
        }
      }
    }

	  T->SwapArrays();
	  T->IncrementStepCount();
  }  // end of current step
  *time_taken = t1.stop();
}

void heat_transfer_calculation_parallel(uint size, TemperatureArray* T, uint steps) {
  timer serial_timer;
  double time_taken = 0.0;
  uint startx = 0;
  uint endx = size - 1;

  // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

  uint min_columns = size / world_size;
  uint excess_columns = size % world_size;
  if (world_rank < excess_columns) {
    startx = world_rank * (min_columns + 1);
    endx = startx + min_columns;
  }
  else {
    startx = (excess_columns * (min_columns + 1)) + ((world_rank-excess_columns) * min_columns);
    endx = startx + min_columns - 1;
  }

  serial_timer.start();
  //*------------------------------------------------------------------------

  heat_transfer_calculation (size, startx, endx, &time_taken, T, steps);

  // Print these statistics for each thread 
  uint a=0;
  if (world_rank==0){
    std::cout << "rank, start_column, end_column, time_taken\n";
  }

  if(world_rank != 0 ){
    MPI_Recv(&a, 1, MPI_UNSIGNED, world_rank -1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }
  
  std::cout<< world_rank << ", "<<startx<<", " << endx << ", " << std::setprecision(TIME_PRECISION) << time_taken << "\n";
  if(world_rank==0){    
    MPI_Send(&a, 1, MPI_UNSIGNED, world_rank +1, 0, MPI_COMM_WORLD);
    MPI_Recv(&a, 1, MPI_UNSIGNED, world_size -1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);   

  }

  if (world_rank== world_size-1){
    MPI_Send(&a, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);

  }
  if(world_rank != 0 && world_rank != world_size-1){
    MPI_Send(&a, 1, MPI_UNSIGNED, world_rank +1, 0, MPI_COMM_WORLD);
  }
//===========================
  if(world_rank != 0 ){
    MPI_Recv(&a, 1, MPI_UNSIGNED, world_rank -1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
  }


  uint step = size/6;
  uint position = 0;
  for(uint i = 0; i < 6; i++) {
      if(startx <= position && position <= endx){
        std::cout<< "Temp[" << position << "," << position << "]=" << T->temp(position,position) << "\n";
      }
      position += step;
  } 
  std::cout<< "Temp[" << endx << "," << endx << "]=" << T->temp(endx,endx) << "\n";

  if (world_rank == world_size-1){
    MPI_Send(&a, 1, MPI_UNSIGNED, 0, 0, MPI_COMM_WORLD);

  }
  if(world_rank != 0 && world_rank != world_size-1){
    MPI_Send(&a, 1, MPI_UNSIGNED, world_rank +1, 0, MPI_COMM_WORLD);
  }
  if(world_rank ==0 ){
    MPI_Send(&a, 1, MPI_UNSIGNED, world_rank +1, 0, MPI_COMM_WORLD);
    MPI_Recv(&a, 1, MPI_UNSIGNED, world_size -1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
  
  

  // Print temparature at select boundary points;  

  //*------------------------------------------------------------------------
  time_taken = serial_timer.stop();
    std::cout << "Time taken (in seconds) : " << std::setprecision(TIME_PRECISION)
              << time_taken << "\n";
  }
  
}

int main(int argc, char *argv[]) {
  // Initialize command line arguments
  cxxopts::Options options("Heat_transfer_calculation",
                           "Model heat transfer in a grid using serial and parallel execution");
  options.add_options(
      "custom",
      {
          {"gSize", "Grid Size",         
           cxxopts::value<uint>()->default_value(DEFAULT_GRID_SIZE)},
          {"mTemp", "Temperature in middle of array",         
           cxxopts::value<double>()->default_value(DEFAULT_MIDDLE_TEMP)},
	        {"iCX", "Coefficient of horizontal heat transfer",
           cxxopts::value<double>()->default_value(DEFAULT_CX)},
          {"iCY", "Coefficient of vertical heat transfer",
           cxxopts::value<double>()->default_value(DEFAULT_CY)},
          {"tSteps", "Time Steps",
           cxxopts::value<uint>()->default_value(DEFAULT_TIME_STEPS)}
      });
  MPI_Init(NULL, NULL);
    // Get the number of processes
  int world_size;
  MPI_Comm_size(MPI_COMM_WORLD, &world_size);

  // Get the rank of the process
  int world_rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
  auto cl_options = options.parse(argc, argv);
  uint grid_size = cl_options["gSize"].as<uint>();
  double init_temp = cl_options["mTemp"].as<double>();
  double Cx = cl_options["iCX"].as<double>();
  double Cy = cl_options["iCY"].as<double>();
  uint steps = cl_options["tSteps"].as<uint>();
  if(world_rank ==0){
    std::cout << "Number of processes : " << world_size<< "\n";
    std::cout << "Grid Size : " << grid_size << "x" << grid_size << "\n";
    std::cout << "Cx : " << Cx << "\n" << "Cy : " << Cy << "\n";
    std::cout << "Temperature in the middle of grid : " << init_temp << "\n";
    std::cout << "Time Steps : " << steps << "\n";
    
    std::cout << "Initializing Temperature Array..." << "\n";
  }

  TemperatureArray *T = new TemperatureArray(grid_size, Cx, Cy, init_temp);
  if (!T) {
      std::cout << "Cannot Initialize Temperature Array...Terminating" << "\n";
      return 2;
  }
  heat_transfer_calculation_parallel (grid_size, T, steps);
  MPI_Finalize();
  delete T;
  return 0;
}
