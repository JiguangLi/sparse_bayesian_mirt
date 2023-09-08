import cvxpy as cp
import numpy as np
import typing


class ProbitSubProblem:
    def __init__(self,
                 y_j: np.ndarray,
                 thetas: np.ndarray,
                 l1_penalty_j: np.ndarray,
                 k: int,
                 num_samples: int,
                 problem_idx: int,
                 loading_constraints: typing.Optional[typing.Dict[typing.Tuple[int, int], float]] = None,
                 ):
        self.thetas = thetas
        self.y_j = y_j
        self.alpha_j = cp.Variable((k, 1))
        self.num_samples = num_samples
        self.d_j = cp.Variable((1, 1))
        self.log_likelihood = cp.sum(cp.multiply(self.y_j, cp.log_normcdf(self.thetas @ self.alpha_j + self.d_j)) +
                                     cp.multiply(1 - self.y_j,
                                                 cp.log_normcdf(-(self.thetas @ self.alpha_j + self.d_j))))
        self.l1_penalty_j = l1_penalty_j
        self.penal_term = cp.norm(cp.multiply(self.alpha_j, self.l1_penalty_j.reshape(-1, 1)), 1)
        self.loading_constraints = loading_constraints
        constraints = []
        for pos, val in self.loading_constraints.items():
            if pos[0] == problem_idx:
                constraints += [self.alpha_j[pos[1]] == val]
        self.problem = cp.Problem(
            cp.Maximize(
                self.log_likelihood / self.num_samples - self.penal_term - 0.5 * cp.square(cp.norm(self.d_j, 2))
            ),
            constraints
        )

    def solve(self):
        self.problem.solve()
        return self.alpha_j.value.flatten(), self.d_j.value[0][0]


class ProbitOptimizationProblem:

    def __init__(self,
                 y: np.ndarray,
                 thetas: np.ndarray,
                 l1_penalty: np.ndarray,
                 num_samples: int,
                 k: int,
                 loading_constraints: typing.Optional[typing.Dict[typing.Tuple[int, int], float]] = None,
    ):
        self.data = y
        self.thetas = thetas
        self.l1_penalty = l1_penalty
        self.n = y.shape[0]
        self.m = y.shape[1]
        self.num_samples = num_samples
        self.k = k
        self.loading_constraints = loading_constraints
        self.problems_dict = {}
        self.construct_sub_problems()

    def constrain_loadings(self, item_idx):
        """Check if there is loading constraints for a given item index"""
        if not self.loading_constraints:
            return False
        else:
            items = [x[0] for x in self.loading_constraints.keys()]
            result = item_idx in items
            return result

    def construct_sub_problems(self):
        for j in range(self.m):
            self.problems_dict[j] = ProbitSubProblem(y_j=np.repeat(self.data[:, j], self.num_samples).reshape(-1, 1),
                                                     thetas=self.thetas,
                                                     l1_penalty_j=self.l1_penalty[j, :],
                                                     k=self.k,
                                                     num_samples=self.num_samples,
                                                     problem_idx=j,
                                                     loading_constraints=self.loading_constraints
                                                     )

    def solve_sub_problem(self, problem_idx):
        return self.problems_dict[problem_idx].solve()