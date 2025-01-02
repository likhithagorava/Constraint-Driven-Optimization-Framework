# Constraint-Driven-Optimization-Framework
## Overview
This project implements a Penalty-Based Constraint Optimization Tool designed to solve constrained optimization problems involving n-variable functions. This framework reduces the constrained optimization problem into an unconstrained function with the help of Quadratic penalty function, which is then solved using advanced optimization techniques like Marquardt's Method, Golden Section Search, and Bounding Phase Method.


## Features
### Penalty Method:
Converts constrained optimization problems into unconstrained optimization by incorporating penalties for constraint violations.
### Marquardt's Method: 
A robust technique for solving single-variable optimization problems.
### Unisearch Implementing Methods:
Bounding Phase Method: Determines the initial inerval for the optimum value.
Golden Section Search: Finds the optimal solution within the identified interval.

## Workflow
### Problem Definition:
Input an n-variable function along with constraints (equality and/or inequality).

### Penalty Method:

Formulate the penalized objective function by adding penalty terms for constraint violations.
### Quadratic Penalty Function
This is one of the most common forms:
The Quadratic Penalty Function is one of the most common forms used in optimization for handling constraints. It can be defined as:

####    p(x) = ∑_i max(0, g_i(x))^2 + ∑_j h_j(x)^2
Where:
g_i(x) represents the inequality constraints g_i(x) ≤ 0.
h_j(x) represents the equality constraints h_j(x) = 0.

This approach ensures penalties are only applied when constraints are violated.

### Apply Marquardt's Method for solving the reduced function.
### Utilize Bounding Phase Method to find a bracket for the minimum.
### Use Golden Section Search to pinpoint the optimal solution within the bracket.

### Iterative Refinement:

Adjust the penalty coefficient iteratively to refine the solution and ensure convergence to the feasible region.
