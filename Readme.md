press `Ctrl + Shift + V` to view this file (if elsewhere from git)

# Introduction

This repository groups two of my past course reports (2022-2023) that involved some python coding : simulation and modelisation (SIMU) and numerical methods (MENU).

## SIMU
In this folder, you will find the python library I coded (in collaboration with my binome) to simulate molecular dynamics (more specifically, the brownian motion). There are two run scripts : 
 - run.py : runs the molecular dynamics simulation up to 1000 frames with 36 particles inside a box with rigid borders and with perfect collisions. This replicates a large set of gas molecules.
 - run_b.py : runs the simuluation of the Brownian motion of a large particle that collides with a large set of smaller particles, which move with different velocities in different random directions.

A specific ReadMe.md file is inside the file that indicates how to run the code. In short, type the following commands in powershell terminal (in the root) : 

```sh 
pip install -r requirements.txt
```

``` sh 
cd MENU
```

``` sh 
python run.py
```

You can also read the file `rapport.md` that explains the results and the project (written in French). 

## MENU

In this folder, there are many notebooks (written in French as well) that explore many notions in numerical calculus. Among them :
- Rapport_Projet_Poisson.ipynb : explains the approach and the results of Poisson's equation 2D numerical solving using many finite differences methods (Jacobi, Gauss-Seidel, adaptative Gauss-Seidel).
- precision.ipynb : explores the numerical epsilon.
- linear_systems.ipynb : numerical methods to solve linear equations (linear algebra, matrix inversion...)  
- Fourrier_transformation.ipynb : explores numerical fourrier transformation methods.
- integrals.ipynb : displays the different numerical integral computation methods.

