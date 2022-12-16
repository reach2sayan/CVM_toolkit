"""
Optimiser module for SRO Correction
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Generator, Type

import numpy as np
from scipy.optimize import minimize
from scipy.optimize import OptimizeWarning, OptimizeResult

from toolkit.constraints.CorrelationConstraints import CorrelationConstraints
from toolkit.bounds.CorrelationBounds import CorrelationBounds
from toolkit.optimizers.ClusterOptimizer import ClusterOptimizer

from toolkit.io.atatio import get_random_structure

@dataclass(kw_only=True, order=False, eq=False,)
class CVMOptimizer(ClusterOptimizer):
    """
    Class to fit a T-dependent function to Short-Range-order corrections
    to CALPHAD free energy
    """

    options: dict
    optimized_result: Type[OptimizeResult] = None
    _T: float = 100
    _structure_generator: Generator[np.ndarray, None, None] = field(init=False)

    def __post_init__(self) -> None:


        super().__post_init__()
        assert self.cluster.vmatrix_array.shape == (len(self._multconfig_kb), len(self._mults_eci))
        self._constraints = CorrelationConstraints(self.cluster.vmatrix_array,
                                                    self.cluster.ordered_correlations,
                                                    self.cluster.disordered_correlations,
                                                    self.norm_constrained
                                                   ).constraints

        self._bounds = CorrelationBounds(self.cluster.num_clusters,
                                         len(self.cluster.single_point_clusters),
                                         self.cluster.disordered_correlations[self.cluster.single_point_clusters]
                                        ).sro_bounds

        self._structure_generator = get_random_structure(self.cluster.structure, self.cluster)

    @property
    def temperature(self):
        return self._T

    @temperature.setter
    def temperature(self, T):
        self._T = T

    def get_energy(self: CVMOptimizer,
                   correlations: np.ndarray
                  ) -> float:
        return self._F(correlations,
                       self._mults_eci,
                       self._multconfig_kb,
                       self.cluster.vmatrix_array,
                       self._vrhologrho,
                       self._T
                      )

    def fit(self: CVMOptimizer) -> (float, np.ndarray, np.ndarray, float):

        result = None
        result_correlations = self.cluster.disordered_correlations.copy()
        result_value = self.get_energy(result_correlations)
        result_constr_viol = np.float64(0.0)
        result_grad = np.zeros(result_correlations.shape[0])

        earlystop = 0
        trial = 0
        for trial in range(self.num_trials):

            accepted = False
            corrs_attempt = next(self._structure_generator)

            if trial == 0:
                corrs_attempt = self.cluster.disordered_correlations.copy()

            f_attempt = self.get_energy(corrs_attempt)
            args = (self._mults_eci, self._multconfig_kb, self.cluster.vmatrix_array, self._vrhologrho, self.temperature,)
            try:
                temp_results = minimize(self._F,
                                        corrs_attempt,
                                        method='trust-constr',
                                        args=args,
                                        options=self.options,
                                        jac=self._dF,
                                        hess=self._d2F,
                                        constraints=self._constraints,
                                        bounds=self._bounds,
                                       )
            except OptimizeWarning as opt_warn:
                print(opt_warn)
                print(f'WARNING. Optimisation Failure: step {trial}, T = {self.temperature}K')
                print(f'Trial Correlations:\n {corrs_attempt}')
                print(f'Trial Configuration Probabilities:\n self.cluster.print_config_probabilities(corrs_attempt)')

            if temp_results.constr_violation < self.constr_tol and temp_results.fun < result_value:

                earlystop = 0
                result = temp_results.copy()

                result_value = temp_results.fun
                result_grad = temp_results.grad.copy()
                result_correlations = temp_results.x.copy()
                result_constr_viol = temp_results.constr_violation

                accepted = True

            if self.print_output:
                print(f'Trial No.: {trial}')
                print(f'Current attempt correlations: {corrs_attempt}')
                print(f'Trial Validity: {self.cluster.check_correlation_validity(corrs_attempt)}')
                print(f'Attempt Free Energy @ T = {self.temperature}K : {f_attempt/self.cluster.num_lat_atoms}')
                print(f'Current Free Energy @ T = {self.temperature}K : {temp_results.fun/self.cluster.num_lat_atoms}')
                print(f'Current minimum correlations: {temp_results.x}')
                print(f"Gradient: {np.array2string(temp_results.grad)}")
                print(f"Constraint Violation: {temp_results.constr_violation}")
                print(f"Stop Status: {temp_results.status} | {temp_results.message}")
                print(f"Acccepted? : {accepted}")
                print(f"Current Min Free Energy @ T = {self.temperature}K : {result_value/self.cluster.num_lat_atoms}")
                print(f"Current Best Contr. Viol : {result_constr_viol}")

            earlystop += 1
            if earlystop > self.early_stopping_count and trial > self.num_trials/2:
                print(f'No improvement for consecutive {self.early_stopping_count} steps. After half of total steps ({int(self.num_trials/2)}) were done')
                break

        self.optimized_result = result
        return (result_value, result_correlations, result_grad, result_constr_viol)

    def __repr__(self):
        if self.optimized_result is None:
            print('No Optimisation has been performed yet')
            print(f'Fitting parameters:\n{self.options}')
            print(f'Constraints:\n{self._constraints}')
            print(f'Approximating Derivatives:\n{self.approx_deriv}')
            print(self._dF)
            print(self._d2F)
        else:
            print(f'Current minimum correlations: {self.optimized_result.x}')
            print(f"Gradient: {np.array2string(self.optimized_result.grad)}")
            print(f"Constraint Violation: {self.optimized_result.constr_violation}")
            print(
                f"Stop Status: {self.optimized_result.status} | {self.optimized_result.message}")
        return ''
