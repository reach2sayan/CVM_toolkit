"""
Optimiser module for SRO Correction
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Generator
import subprocess

import numpy as np
from scipy.optimize import minimize
from scipy.optimize import OptimizeWarning

from ..constraints.CorrelationConstraints import CorrelationConstraints
from ..bounds.CorrelationBounds import CorrelationBounds
from .ClusterOptimizer import ClusterOptimizer

from ..io.atatio import get_random_structure

@dataclass(kw_only=True, order=False, eq=False, slots=True)
class ShortRangeOrderOptimization(ClusterOptimizer):

    options: dict
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

        self._structure_generator = get_random_structure(self.structure, self.cluster)

    def fit(self: ShortRangeOrderOptimization,
            params: dict
           ) -> (float, np.ndarray, np.ndarray, float):

        temperature = params['temperature']
        result_correlations = self.cluster.disordered_correlations.copy()
        result_value = self._F(result_correlations,
                               self._mults_eci,
                               self._multconfig_kb,
                               self.cluster.vmatrix_array,
                               self._vrhologrho,
                               temperature
                              )
        result_constr_viol = np.float64(0.0)
        result_grad = np.zeros(result_correlations.shape[0])

        earlystop = 0
        trial = 0
        for trial in range(self.num_trials):

            accepted = False
            corrs_attempt = next(self._structure_generator)

            if trial == 0:
                corrs_attempt = self.cluster.disordered_correlations.copy()

            f_attempt = self._F(corrs_attempt,
                                self._mults_eci,
                                self._multconfig_kb,
                                self.cluster.vmatrix_array,
                                self._vrhologrho,
                                temperature
                               )

            args = (self._mults_eci, self._multconfig_kb, self.cluster.vmatrix_array, self._vrhologrho, temperature,)
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
                print(f'WARNING. Optimisation Failure: step {trial}, T = {temperature}K')
                print(f'Trial Correlations:\n {corrs_attempt}')
                print(f'Trial Configuration Probabilities:\n self.cluster.print_config_probabilities(corrs_attempt)')

            if temp_results.constr_violation < self.constr_tol and temp_results.fun < result_value:

                earlystop = 0
                result_value = temp_results.fun
                result_grad = temp_results.grad.copy()
                result_correlations = temp_results.x.copy()
                result_constr_viol = temp_results.constr_violation

                accepted = True

            if self.print_output:
                print(f'Trial No.: {trial}')
                print(f'Current attempt correlations: {corrs_attempt}')
                print(f'Trial Validity: {self.cluster.check_result_validity(corrs_attempt)}')
                print(f'Attempt Free Energy @ T = {temperature}K : {f_attempt/self.num_lat_atoms}')
                print(f'Current Free Energy @ T = {temperature}K : {temp_results.fun/self.num_lat_atoms}')
                print(f'Current minimum correlations: {temp_results.x}')
                print(f"Gradient: {np.array2string(temp_results.grad)}")
                print(f"Constraint Violation: {temp_results.constr_violation}")
                print(
                    f"Stop Status: {temp_results.status} | {temp_results.message}")

                print(f"Acccepted? : {accepted}")
                print(f"Current Min Free Energy @ T = {temperature}K : {result_value/self.num_lat_atoms}")
                print(f"Current Best Contr. Viol : {result_constr_viol}")
                print('\n====================================\n')

            earlystop += 1
            if earlystop > self.early_stopping_count and trial > self.num_trials/2:
                print(
                    f'No improvement for consecutive {self.early_stopping_count} steps. After half of total steps ({int(self.num_trials/2)}) were done')
                break

        return (result_value, result_correlations, result_grad, result_constr_viol)
