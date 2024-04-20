# MPI-Learning
Explora ejercicios prácticos de MPI (Message Passing Interface) para mejorar tus habilidades en programación paralela. Desde conceptos básicos hasta desafíos avanzados, este repositorio ofrece oportunidades de aprendizaje progresivo en computación de alto rendimiento.

# Exercise 1: MPI Vector Addition

## Overview
This program demonstrates vector addition using MPI for parallel computing. MPI allows distributing tasks across multiple processes, enabling faster computation for large datasets.

## Explanation of MPI Functions

### MPI_Init
- **Functionality**: Initializes the MPI environment.
- **Parameters**: 
  - `argc`: Number of command line arguments.
  - `argv`: Command line arguments.
- **Usage**: `MPI_Init(&argc, &argv);`

### MPI_Comm_size
- **Functionality**: Retrieves the total number of processes in the MPI communicator.
- **Parameters**: 
  - `comm`: MPI communicator.
  - `size`: Output parameter to store the number of processes.
- **Usage**: `MPI_Comm_size(MPI_COMM_WORLD, &size);`

### MPI_Comm_rank
- **Functionality**: Retrieves the rank (identifier) of the calling process within the communicator.
- **Parameters**: 
  - `comm`: MPI communicator.
  - `rank`: Output parameter to store the process rank.
- **Usage**: `MPI_Comm_rank(MPI_COMM_WORLD, &rank);`

### MPI_Finalize
- **Functionality**: Finalizes the MPI environment.
- **Usage**: `MPI_Finalize();`

## Code Functionality
1. **Initialization**: 
   - MPI environment is initialized, and process information is retrieved.

2. **Work Division**: 
   - The vector addition task is divided among processes, with each assigned a portion of the vectors to add.

3. **Vector Addition**: 
   - Each process computes the addition of its assigned portion of the vectors.

4. **Printing Results**: 
   - Each process prints the portion of the vectors it computed.

5. **Printing Process Information**: 
   - Each process prints its ID, the number of elements it processed, the start and end indices of its portion of the vectors, and the time taken for computation.

6. **Finalization**: 
   - MPI environment is finalized.

## Running the Code
To compile and run the program:
```sh
mpicxx Exercise1_MPI.cpp -o Exercise1_MPI
mpiexec -n 4 ./Exercise1_MPI
```
Replace `4` with the desired number of processes.

# Exercise 2: MPI Sum of Elements of a Vector

## Overview
This program calculates the sum of elements of a vector using MPI for parallel computing. It distributes the workload among multiple processes to compute the sum more efficiently.

## Explanation of New MPI Functions

### MPI_Send
- **Functionality**: Sends data from one process to another.
- **Parameters**:
  - `buf`: Pointer to the data to be sent.
  - `count`: Number of elements in the data.
  - `datatype`: Data type of the elements.
  - `dest`: Rank of the destination process.
  - `tag`: Message tag.
  - `comm`: MPI communicator.
- **Usage**: `MPI_Send(&suma, 1, MPI_DOUBLE_PRECISION, 0, 1000, MPI_COMM_WORLD);`

### MPI_Recv
- **Functionality**: Receives data from another process.
- **Parameters**:
  - `buf`: Pointer to the buffer where received data will be stored.
  - `count`: Number of elements in the buffer.
  - `datatype`: Data type of the elements.
  - `source`: Rank of the source process.
  - `tag`: Message tag.
  - `comm`: MPI communicator.
  - `status`: Status object.
- **Usage**: `MPI_Recv(&sumaLocal, 1, MPI_DOUBLE_PRECISION, i, 1000, MPI_COMM_WORLD, MPI_STATUS_IGNORE);`

## Code Functionality
1. **Initialization**: 
   - MPI environment is initialized, and process information is retrieved.

2. **Work Division**: 
   - The vector sum task is divided among processes, with each assigned a portion of the vector elements.

3. **Sum Calculation**: 
   - Each process calculates the sum of its assigned portion of vector elements.

4. **Sending and Receiving Sums**: 
   - If there are multiple processes, each process sends its local sum to the root process (process 0), which collects and computes the global sum.

5. **Printing Results**: 
   - Each process prints its ID, the number of elements it processed, the start and end indices of its portion of the vector, its local sum, the global sum (if applicable), and the time taken for computation.

6. **Finalization**: 
   - MPI environment is finalized.

## Running the Code
To compile and run the program:
```sh
mpicxx Exercise2_MPI.cpp -o Exercise2_MPI
mpiexec -n 4 ./Exercise2_MPI
```

# Exercise 3: MPI Communication Methods

## Overview
This program demonstrates various MPI communication methods for parallel computing. It calculates the sum of elements in a vector using different MPI communication functions, such as send and recv, reduce, all reduce, gather, and all gather.

## Code Functionality
1. **Initialization**: 
   - MPI environment is initialized, and process information is retrieved.

2. **Input Interval**: 
   - Process 0 reads the number of intervals (n) from the user.

3. **Broadcast Interval**: 
   - Process 0 broadcasts the value of n to all other processes.

4. **Interval Calculation**: 
   - Each process calculates its portion of the interval and computes the sum of the function `f(x)`.

5. **MPI Communication Methods**:
   - Different MPI communication functions are used to exchange data or perform operations across processes.

6. **Print Results**: 
   - Process 0 prints the calculated value of π and the error compared to the known value of π (π = 3.141592653589793238462643).

7. **Print CPU Time**: 
   - Each process prints its CPU time for computation.

8. **MPI Finalization**: 
   - MPI environment is finalized.

## Explanation of MPI [4] - Communication Functions
### Send and Recv
- **Functionality**: Sends data from one process to another process and receives data.
- **Usage**: 
  ```cpp
  if (numtasks > 1){
      if (taskid != 0) {
            MPI_Send(&suma,1,MPI_DOUBLE_PRECISION,0,1000,MPI_COMM_WORLD);
        } else {
            printf(" \n");
            sumaglob = suma; 
            for (int i = 1; i < numtasks; i++){
                MPI_Recv(&sumaglobito,1,MPI_DOUBLE_PRECISION,i,1000,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
                sumaglob += sumaglobito;
            }
        }
   } else {
       sumaglob = suma;
   }
  ```
### Reduce
- **Functionality**: Reduces values on all processes to a single value using an operation (e.g., sum).
- **Usage**: 
  ```cpp
  MPI_Reduce(&suma, &sumaglob, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  ```
### All Reduce
- **Functionality**: Combines values from all processes and distributes the result to all processes.
- **Usage**: 
  ```cpp
  MPI_Allreduce(&suma, &sumaglob, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  ```
### Gather
- **Functionality**: Gathers data from all processes to one process.
- **Usage**: 
  ```cpp
  if (taskid == 0) {
        sumasParciales = (double *)malloc(numtasks * sizeof(double));
   }
  MPI_Gather(&suma, 1, MPI_DOUBLE, sumasParciales, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  ```
### All Gather
- **Functionality**: Gathers data from all processes and distributes the data to all processes.
- **Usage**: 
  ```cpp
  MPI_Allgather(&suma, 1, MPI_DOUBLE, sumaux, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  ```
## How to Run
To compile and run the program:
```sh
mpicxx communication_method_MPI.cpp -o communication_method_MPI
mpiexec -n 4 ./communication_method_MPI
```

# Exercise 4: MPI Calculate PI

## Overview
This program calculates an approximation of the mathematical constant π (pi) using MPI for parallel computing.

## Explanation
The program divides the calculation of π into multiple intervals and distributes the workload among different processes to compute the sum in parallel. It then reduces the partial sums computed by each process to obtain the final result.

## Code Functionality
1. **Initialization**: 
   - MPI environment is initialized, and process information is retrieved.

2. **Input Interval**: 
   - Process 0 reads the number of intervals (n) from the user.

3. **Broadcast Interval**: 
   - Process 0 broadcasts the value of n to all other processes.

4. **Interval Calculation**: 
   - Each process calculates its portion of the interval and computes the sum of the function `f(x)`.

5. **Printing Results**: 
   - Process 0 prints the calculated value of π and the error compared to the known value of π (π = 3.141592653589793238462643). Each process prints its CPU time for computation.

6. **Finalization**: 
   - MPI environment is finalized.

## How to Run
To compile and run the program:
```sh
mpicxx Exercise4_MPI.cpp -o Exercise4_MPI
mpiexec -n 4 ./Exercise4_MPI
```