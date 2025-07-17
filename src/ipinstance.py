import heapq
import random
from calendar import c
from dataclasses import dataclass
from itertools import combinations, count
from typing import Any, Literal

import numpy as np
from docplex.mp.model import Model


@dataclass(frozen=True)
class IPConfig:
    numTests: int  # number of tests
    numDiseases: int  # number of diseases
    costOfTest: np.ndarray  # [numTests] the cost of each test
    # [numTests][numDiseases] 0/1 matrix if test is positive for disease
    A: np.ndarray


#  * File Format
#  * #Tests (i.e., n)
#  * #Diseases (i.e., m)
#  * Cost_1 Cost_2 . . . Cost_n
#  * A(1,1) A(1,2) . . . A(1, m)
#  * A(2,1) A(2,2) . . . A(2, m)
#  * . . . . . . . . . . . . . .

def data_parse(filename: str):
    try:
        with open(filename, "r") as fl:
            numTests = int(fl.readline().strip())  # n
            numDiseases = int(fl.readline().strip())  # m

            costOfTest = np.array([float(i)
                                  for i in fl.readline().strip().split()])

            A = np.zeros((numTests, numDiseases))
            for i in range(0, numTests):
                A[i, :] = np.array([int(i)
                                   for i in fl.readline().strip().split()])
            return numTests, numDiseases, costOfTest, A
    except Exception as e:
        print(f"Error reading instance file. File format may be incorrect.{e}")
        exit(1)


class IPInstance:

    def __init__(self, filename: str) -> None:
        numT, numD, cst, A = data_parse(filename)
        self.numTests = numT
        self.numDiseases = numD
        self.costOfTest = cst
        self.A = A
        self.model = Model()  # CPLEX solver

    def toString(self):
        out = ""
        out = f"Number of test: {self.numTests}\n"
        out += f"Number of diseases: {self.numDiseases}\n"
        cst_str = " ".join([str(i) for i in self.costOfTest])
        out += f"Cost of tests: {cst_str}\n"
        A_str = "\n".join([" ".join([str(j) for j in self.A[i]])
                          for i in range(0, self.A.shape[0])])
        out += f"A:\n{A_str}"
        return out

    # returns minimum cost of tests nessecary to fully distinguishing diseases
    def solve(self) -> int:
        # Build model
        model = Model()
        x = {i: model.binary_var(name=f"x_{i}") for i in range(self.numTests)}

        # Objective
        model.minimize(
            model.sum(self.costOfTest[i] * x[i] for i in range(self.numTests)))

        # Cover constraints
        for j, k in combinations(range(self.numDiseases), 2):
            model.add_constraint(
                model.sum(x[i]
                          for i in range(self.numTests) if self.A[i][j] != self.A[i][k]) >= 1
            )

        # Symmetry-breaking: pairwise check identical tests
        for i, j in combinations(range(self.numTests), 2):
            if i < j and all(self.A[i][k] == self.A[j][k] for k in range(self.numDiseases)):
                model.add_constraint(x[i] > x[j])

        # Solve with branch-and-bound
        ub = float('inf')
        result = self.best_first_branch_and_bound(model, ub)
        return int(result) if result is not None else -1

    def best_first_branch_and_bound(self, model: Model, ub: float) -> float:
        if not model.solve():
            return ub  # infeasible

        tol = 1e-5
        pq = [(-model.objective_value, 0, model, set(range(self.numDiseases)))]
        cnt = count(1)

        # Get list of variable names that actually exist in the model
        var_names = {var.name: i for i,
                     var in enumerate(model.iter_variables())}

        while pq:
            neg_lb, _, node, active = heapq.heappop(pq)
            lb = -neg_lb
            if lb >= ub - tol:
                continue  # prunable

            # Only consider variables that exist in the model
            scores = []
            for var_name, i in var_names.items():
                if not var_name.startswith('x_'):
                    continue
                var = node.get_var_by_name(var_name)
                if var is None:
                    continue
                if abs(var.solution_value - round(var.solution_value)) > tol:
                    p = sum(self.A[i, j] for j in active)
                    n = len(active)
                    if 0 < p < n:
                        # Calculate how close the split is to 50/50
                        split_quality = 1 - abs(p/n - 0.5)  # Will be 1 for perfect 50/50, 0 for worst split
                        scores.append((split_quality / self.costOfTest[i], i, var_name))

            if not scores:
                ub = min(ub, lb)
                continue

            _, _, name = max(scores)

            for val in (0, 1):
                child = node.clone()
                child.add_constraint(child.get_var_by_name(name) == val)
                # Extract index from variable name
                idx = int(name.split('_')[1])
                new_active = {j for j in active if self.A[idx, j] == val}
                if child.solve() and (new_lb := child.objective_value) < ub - tol:
                    heapq.heappush(pq, (-new_lb, next(cnt), child, new_active))

        return ub

    def greedy_initial_ub(self) -> int:
        pairs = {(i, j) for i in range(self.numDiseases)
                 for j in range(i+1, self.numDiseases)}
        chosen = set()
        while pairs:
            best, best_ratio = None, float("inf")
            for test in range(self.numTests):
                if test in chosen:
                    continue
                count = sum(
                    1 for i, j in pairs if self.A[test][i] != self.A[test][j])
                if count and (ratio := self.costOfTest[test] / count) < best_ratio:
                    best_ratio, best = ratio, test
            if best is None:
                break
            chosen.add(best)
            pairs = {(i, j)
                     for i, j in pairs if self.A[best][i] == self.A[best][j]}
        return int(sum(self.costOfTest[test] for test in chosen))

    def branch_and_bound(self, model: Model, ub: float) -> float:
        tol = 1e-5
        if not model.solve():
            return ub  # infeasible

        lb = model.objective_value
        if lb >= ub - tol:
            return ub  # prunable

        x_vars = [(i, model.get_var_by_name(f"x_{i}"))
                  for i in range(self.numTests)]

        # If all variables are (nearly) integer, update ub
        if all(abs(var.solution_value - round(var.solution_value)) < tol for _, var in x_vars):
            return min(ub, lb)

        # Pick the first fractional variable in the specified order
        for i, var in x_vars:
            if abs(var.solution_value - round(var.solution_value)) > tol:
                frac_index = i
                break

        # Branch: enforce x[frac_index] == 0
        model0 = model.clone()
        model0.add_constraint(model0.get_var_by_name(f"x_{frac_index}") == 0)
        ub = min(ub, self.branch_and_bound(model0, ub))

        # Branch: enforce x[frac_index] == 1
        model1 = model.clone()
        model1.add_constraint(model1.get_var_by_name(f"x_{frac_index}") == 1)
        ub = min(ub, self.branch_and_bound(model1, ub))

        return ub

    def ip_solve(self) -> int:
        # uses CPLEX to verify optimal values
        x = self.model.binary_var_list(self.numTests, name="x")
        self.model.minimize(sum(self.costOfTest[i] * x[i]
                            for i in range(self.numTests)))

        self.model.add_constraints(
            (self.model.sum(x[i] for i in range(self.numTests)
                            if self.A[i][j] != self.A[i][k]) >= 1)
            for j, k in combinations(range(self.numDiseases), 2))

        if solution := self.model.solve():
            return int(self.model.objective_value)
        else:
            return -1
