# Fibonacci

**Name**: Fibonacci  
**Contact Person**: Antoni Navarro, (antoni.navarro@bsc.es)  
**Platform**: OmpSs-2  

## Description
The Fibonacci sequence is a sequence of numbers characterized by the fact that every number after the first two is the sum of the two preceding ones. It usually starts by either the pair (1 1) or (0 1).  

For additional information refer to:  
* https://en.wikipedia.org/wiki/Fibonacci_number

### General Algorithm
This benchmark contains:
* A recursive implementation using taskwaits.  
* A recursive implementation without the usage of taskwaits.  

## Build instructions
Building of the application(s) is as easy as:  
```
make fibonacci|fibonacci_notaskwait|all

```

## Run instructions
```
./fibonacci|fibonacci_notaskwait <R> <N> <F>
```

where,
* R = The number of times the algorithm has to execute (repetitions)
* N = The size of the problem / how many numbers must be generated (fibonacci sequence size)
* F = The maximum depth (numbers generated) in which tasks become final^*

## Others
* ^* Information about the final clause can be found at https://pm.bsc.es/ompss-2-docs/spec/directives/task.html?highlight=final.  
