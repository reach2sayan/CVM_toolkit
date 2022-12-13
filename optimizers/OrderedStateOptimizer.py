"""
Optimiser module for SRO Correction
"""

from __future__ import annotations
from dataclasses import dataclass, field
import numpy as np
from scipy.optimize import linprog, OptimizeWarning

from .ClusterOptimizer import ClusterOptimizer


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

    def __post_init__(self) -> None:
        self._bounds = [(self._disordered_correlations[idx], self._disordered_correlations[idx]) if cluster['type'] == 1 else (
            1, 1) if cluster['type'] == 0 else (-1, 1) for idx, cluster in self.cluster._clusters.items()]


    def fit(self):

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
            print('Ordered State calculations completed...')
            if self.print_output:
                np.savetxt('ordered_correlations.out', result.x)
                with open('ordered_rho.out', 'w', encoding='utf-8') as frho:
                    for vmat in self.cluster._vmat.values():
                        frho.write(f'{" ".join(map(str,vmat@result.x))}\n')
            return result
        print(
            f'WARNING: linear programming for ordered correlation search failed: {result.status} - {result.message}\nExiting...')
        return result
