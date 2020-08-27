# sequence-alignment-openMP
COMP90025 - Parallel and Multicore Computing - 2020S2 - Assignment2A

## how to run
### windows
#### sequential
- ```g++ -std=c++14 seqalign.cpp -o sequential```
- ```sequential.exe < seq.dat```

#### openMP
- ```g++ -fopenmp -std=c++14 seqalign_parallel.cpp -o parallel```
- ```parallel.exe < seq.dat```

### spartan
- ```sbatch sequential.slurm``` 
- ```sbatch parallel.slurm```

#### see jobs under execution
- ```squeue -u xuliny```

#### Spartan Weather Report
- ```spartan-weather```