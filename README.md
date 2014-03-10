# Molecular Dynamics


## General

This code supports the training on Modern Fortran, OpenMP and OpenACC 
that I teach for LSU and LONI HPC users. This code is not for 
learning Molecular Dynamics (though I might use it in the future). 
Can also be used for a tutorial on Makefiles. 

## Objective

Take the original code and rewrite using Modern Fortran concepts 
learned viz,

   * spit code into smaller subunits: modules, functions and subroutines
   * generalize code so that the following parameters are read as input
      + number of atoms or number of unit cells (you can't do both)
      + simulation/equilibration temperature
      + number of times steps
      + you will need to make use of allocatable arrays
   * use object oriented concepts such as
      + derived type data types
      + generic procedures and operator overloading
   * Extra exercise: add another potential e.g. Morse potential for 
     the simulation

Parallelize the code using

   1. OpenMP
   2. OpenACC
   3. CUDA Fortran (when I learn how to do it)

## Code Description
   * Original Code: md-orig.f90
      + This is the original code that students need to to work on. 
        The Objective is to modify this code using as many features
        as you have grasped.
   * Modified Code: md-v1.f90
      + This version splits the original code into many subprograms.
   * Modified Code: md-v2.f90
      + This version generalizes the size of the program, reads input 
        parameters and uses allocatable arrays. Arrays are passed as
        explicit shape arrays.
   * Modified Code: md-v3.f90
      + In this version, arrays are passed as assumed shape arrays 
        with intent attribute creating a need for interface.
      + A module for calculating the potential and force is added.
      + Morse potential is also added
      + If this is final version, the only thing left for you to 
        learn is Object Oriented concepts. Congratulations!
   * Modified Code: md-v4.f90
      + Derived Types for MD variables are introduced.
   * Oject Oriented Code : md-v5.f90
      + Adds generic procedures and operator overloading to the 
        md-v4.f90 code. 
      + If this is your final version, you have learned a lot of
        Modern Fortran Concepts.  Congratulations!
   * OpenMP Code  : md-omp.f90
      + Parallelize the md.f90 code using OpenMP directives.
   * OpenACC Code : md-acc.f90
      + Parallelize the md.f90 code using OpenACC directives.
        For this code to work, you need to have access to a machine 
        with an accelerator such as NVIDIA GPUs and an OpenACC 
        supported compiler such as Portland Group Compilers. 
        It should work with Cray or CAPS Compilers but I don't 
        have access to one.

## Requires

Compilers: 

   * gfortran, 
   * intel fortran or 
   * portland group compilers
 
For OpenACC/CUDA, portland group compilers is required and a machine with
   NVIDIA GPUs (code is tested on NVIDIA Tesla M2090, don't have access to
   AMD GPUs)

## Compiling

```
$ make [option] [COMP=compiler]
```
options: 

   * Serial  : v{0-5}
   * OpenMP  : omp
   * OpenACC : acc

compiler:

   * gfortran : gcc
   * intel    : intel
   * portland : pgi
   * IBM XL   : ibm (Makefile doesn't work on AIX)

## Input Files

 * md.inp : Input file for version 2-5, openmmp and openacc (default potential is lj)
 * md.in : Input file for version 3-5, openmp and openacc
    
## Running

Serial Code
```
$ ./md{0-5} [<input, needed for version 2-5>]
```

OpenMP Code
```
$ OMP_NUM_THREADS=x ./mdo <input>
```
where x is between 1 and number of processor cores

OpenACC Code
```
$ ./mda <input>
```

## Output Files

The code will produce output to screen, you can redirect output to a file and analyze the simulation temperature and energy.
```
$ ./md0 > md.out
```

The code also produces a atom.xyz file with atomic coordinates at each time step.

## Visualization of MD Simulation using VMD

The VMD script, atom.vmd will read in the atom.xyz file and enable one to view the simulation
```
$ vmd -e atom.vmd
```

## Author

Alex Pacheco, 

HPC User Services Consultant,

Louisiana State University

## Credits

None yet 

