"""
Optimiser module for SRO Correction
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Type
import numpy as np
from scipy.optimize import linprog
from scipy.optimize import OptimizeWarning, OptimizeResult

from toolkit.optimizers.ClusterOptimizer import ClusterOptimizer


@dataclass(kw_only=True, order=False, eq=False, slots=True,)
class OrderedStateOptimizer(ClusterOptimizer):
    """
    Finds the ordered structure given cluster information using Linear Programming
    Input:
        cluster_data - ClusterInfo object contatning vmat, eci and cluster information
        corr - sample corr, only used as a guess for the refined simplex method. This also
        specified the point correlations that are kept fixed.
        method - method of the linear programming; simplex or interior point
        options - extra options for the linear programming problem
        Output
    """

    options: dict
    method: str = 'revised-simplex'
    optimized_result: Type[OptimizeResult] = None

    def __post_init__(self: OrderedStateOptimizer) -> None:
        self._bounds = [(self.cluster.disordered_correlations[idx], self.cluster.disordered_correlations[idx]) if cluster['type'] == 1 else (
            1, 1) if cluster['type'] == 0 else (-1, 1) for idx, cluster in self.cluster.clusters.items()]

    def fit(self: OrderedStateOptimizer) -> None:

        obj = self.cluster.eci_array * self.cluster.clusmult_array
        try:
            result = linprog(obj,
                             A_ub=-1 * self.cluster.vmatrix_array,
                             b_ub=np.zeros(self.cluster.vmatrix_array.shape[0]),
                             bounds=self._bounds,
                             options=self.options,
                             method=self.method
                            )
        except OptimizeWarning as opt_warn:
            print(opt_warn)
        if result.success:
            self.cluster.ordered_correlations = result.x
            self.optimized_result = result
            print('Ordered State calculations completed...')
            if self.print_output:
                np.savetxt('ordered_correlations.out', result.x)
                with open('ordered_rho.out', 'w', encoding='utf-8') as frho:
                    for vmat in self.cluster.vmat.values():
                        frho.write(f'{" ".join(map(str,vmat@result.x))}\n')
        else:
            print(
                f'WARNING: linear programming for ordered correlation search failed: {result.status} - {result.message}\nExiting...')

    def __repr__(self) -> str:
        print('\n--------------------------')
        print('Ordered State Description:')
        if self.optimized_result is None:
            print('No Optimisation has been performed yet')
            print(f'Fitting parameters:\n{self.options}')
        else:
            print(f'Ordered Correlations: \n{self.optimized_result.x}')
            print(f'Ordered State Free Energy: {self.optimized_result.fun}')
        print('-----------------------------')
        return ''
