# sequence-alignment-openMP
COMP90025 - Parallel and Multicore Computing - 2020S2 - Assignment2A

## specification
- [specification](./docs/Project2A.pdf)

## report
- [report](./docs/COMP90025_Assignment_02a.pdf)

## code submission
- [code](./submit/xuliny-seqalignomp.cpp)

## result
- Finally, 8.2 out of 10
- <img src="./docs/code result.png"/>
- <img src="./docs/report result.png"/>
- <img src="./docs/ranking.png"/>
### class ranking distribution
- <img src="./docs/1.png"/>
- <img src="./docs/2.png"/>
- <img src="./docs/3.png"/>
- <img src="./docs/4.png"/>
- <img src="./docs/results cutoff.png"/>

## submit
- ```cd submit```
- ```/data/gpfs/projects/punim0520/pmc/pub/bin/handin seqalignomp xuliny-seqalignomp.cpp```

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

## self test
- ```parallel7.exe < easy.dat > parallel_easy.txt```
- ```sequential.exe < easy.dat > sequential_easy.txt```
- ```parallel7.exe < medium.dat > parallel_medium.txt```
- ```sequential.exe < medium.dat > sequential_medium.txt```
- ```parallel7.exe < hard.dat > parallel_hard.txt```
- ```sequential.exe < hard.dat > sequential_hard.txt```