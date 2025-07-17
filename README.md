# Healthcare Analytics: Test Set Cover Problem

## Problem Description

This project solves the **Test Set Cover Problem** for healthcare analytics. Given a set of medical tests and diseases, the goal is to find the minimum cost collection of tests that can distinguish between every pair of diseases.

### Input Format

- Number of tests (n)
- Number of diseases (m)
- Cost of each test
- Binary matrix A where A[i][j] = 1 if test i is positive for disease j

### Objective

Minimize the total cost of selected tests while ensuring that for every pair of diseases, at least one selected test gives different results for the two diseases.

## Solution Approach

The solution implements a **custom branch-and-bound algorithm** with the following optimizations:

### Core Algorithm

- **Best-First Search**: Uses a priority queue to explore most promising nodes first
- **Linear Programming Relaxation**: Gets lower bounds by solving LP relaxation with CPLEX
- **Intelligent Branching**: Selects variables based on split quality and cost efficiency
- **Symmetry Breaking**: Adds constraints to eliminate equivalent solutions for identical tests

### Key Features

- **Pruning**: Eliminates branches when lower bound exceeds current best solution
- **Greedy Upper Bound**: Initial feasible solution to improve pruning
- **Active Disease Tracking**: Maintains which diseases still need to be distinguished in each branch

## Usage

### Setup

```bash
# Compile (installs dependencies and CPLEX)
./compile.sh

# Run on a single instance
./run.sh <input_file>

# Run all test instances
./runAll.sh
```

### Input Files

Test instances are located in the `input/` directory with naming convention:
`{num_tests}_{num_diseases}_{density}_{instance_id}.ip`

Example: `25_25_0.25_1.ip` = 25 tests, 25 diseases, 25% density, instance #1

## Results

The solver finds optimal solutions efficiently across various problem sizes:

| Instance Size | Avg Time (s) | Example Result |
| ------------- | ------------ | -------------- |
| 25×25         | 0.13-0.14    | Cost: 30-57    |
| 50×50         | 0.36         | Cost: 53       |
| 50×100        | 1.48-1.64    | Cost: 92-103   |
| 100×100       | 2.34-3.06    | Cost: 71-153   |
| 100×200       | 11.0-12.29   | Cost: 69-86    |

All solutions are proven optimal ("OPT") using the branch-and-bound algorithm.

## Technical Implementation

### Dependencies

- **Python 3.9+**
- **NumPy**: Matrix operations and data handling
- **CPLEX/DOcplex**: Linear programming solver for bounds and verification

### Key Files

- `src/main.py`: Entry point and result formatting
- `src/ipinstance.py`: Core IP formulation and branch-and-bound solver
- `src/model_timer.py`: Performance timing utility

### Output Format

Results are logged in JSON format:

```json
{
	"Instance": "25_25_0.25_1.ip",
	"Time": "0.13",
	"Result": "57",
	"Solution": "OPT"
}
```

## Mathematical Formulation

**Variables**: x[i] ∈ {0,1} for each test i

**Objective**: Minimize Σ(cost[i] × x[i])

**Constraints**: For each disease pair (j,k), ensure at least one test distinguishes them:
Σ(x[i] : A[i][j] ≠ A[i][k]) ≥ 1

## Team

- carmstr8

## Course

CSCI 2951-O - Optimization Methods in Computer Science
