# 2d-heat-conduct
Generate temperatures for a 2-dimensional heat conductivity problem.

## Building
To compile the code (requires GNU C compiler) execute:
```
$ gcc -o 2d-heat-conduct -fopenmp 2d-heat-conduct.c
```

To read the manpage execute:
```
$ man ./2d-heat-conduct.1
```

## Running
An example execution of the program is:
```
$ ./2d-heat-conduct -a jacobi -d 50 100 150 200 -s 5 5 -S static -t 0.001 -w
```
NOTE: Please refer to manpage for definition of options required to execute the program.

## Abstract
The goal is to solve a rectangular 2-dimensional heat conductivity problem with two different iterative algorithms, namely Jacobi and Gauss-Seidel, which should, in theory, be faster with the Gauss-Seidel algorithm, from serial to parallel versions with various numbers of threads. It is evident from observations of practical experimentation that the optimised parallel version for the Gauss-Seidel method is consistently faster in execution time as a result of fewer iterations and converging faster with less error.

## Further information
In relation to the original code, it was observed to have poor performance, which was identified as iterating through too many conditional statements relating to the algorithm. Even though this makes for reduced code, it comes with a performance penalty, which was refactored by putting the conditional statements in relation to the algorithm outside of the iterations.

This latest implementation had a dramatic increase in the speed of execution, to an order of magnitude. The main difference between the Jacobi and Gauss-Seidel implementations is that the Jacobi algorithm executes more time-consuming iterations than the Gauss-Seidel method. Further optimisation was pursued by nesting OpenMP parallel directives.

In conclusion, the main difference in the Jacobi and Gauss-Seidel algorithms is that for the Jacobi algorithm it depends on values acquired from a previous point in time from a former iteration, whereas with Gauss-Seidel it uses the most recent values during iteration. However, only the serial version is a true Gauss-Seidel implementation and parallelising the serial version creates a quasi-Gauss-Seidel (a hybrid with Jacobi), which enables the correct iteration to occur, otherwise the implementation of a true Gauss-Seidel with parallel code will not iterate correctly, as was observed during the implementation with the optimisation of the parallel code. Worth noting is the fact that the Gauss-Seidel method should in theory converge faster with less error than that of the Jacobi version, which is evident in the results obtained from benchmarking metrics.

As referenced from benchmark performance results, quasi-Gauss-Seidel truly makes an impact with a low tolerance and a high number of threads, whereby in this scenario it is dramatically faster than the Jacobi method, the reason being that the quasi-Gauss-Seidel method does fewer iterations. However, under circumstances of low tolerance and increased number of threads, the compiler optimisation has a negligible improvement on both the Jacobi and quasi-Gauss-Seidel methods, even with no compiler optimisation “-O0”, which usually is much slower in execution.

Consistently for both algorithms, the optimised parallel code is always progressively faster concerning the increased number of threads. However, further optimisation was researched by using different scheduling, although the dynamic scheduler was consistently ever so slightly slower than that of the static scheduler. To conclude, the Gauss-Seidel algorithm is consistently faster than the Jacobi algorithm in the vast majority of scenarios, therefore proving the theory to support this observation from practical experimentation.
