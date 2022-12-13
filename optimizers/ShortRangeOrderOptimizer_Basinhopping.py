"""
Optimiser module for SRO Correction
"""

from __future__ import annotations
from dataclasses import dataclass, field
from typing import Callable, Union

import numpy as np
from scipy.optimize import basinhopping
from scipy.optimize import OptimizeWarning

from ..constraints.CorrelationConstraints import CorrelationConstraints
from ..bounds.CorrelationBounds import CorrelationBounds

from ..basinhopping_features.basinhopping_features import BasinHoppingBounds, BasinHoppingStep
from .ClusterOptimizer import ClusterOptimizer

@dataclass(kw_only=True, order=False, eq=False, slots=True)
class ShortRangeOrderOptimization_Basinhopping(ClusterOptimizer):

    options: dict
    stepsize: float
    bh_temperature: float = 1.0
    bh_interval: int = 50

    _bh_callback: Callable[[float,float,bool],Union[None,bool]] = field(init=False)
    _bh_accept: Callable[[dict],bool] = field(init=False)
    _bh_step: Callable[[np.ndarray],np.ndarray] = field(init=False)

    def __post_init__(self) -> None:

        super().__post_init__()
        self._constraints = CorrelationConstraints(self.cluster.vmatrix_array,
                                                  self.cluster.ordered_correlations,
                                                  self.cluster.disordered_correlations,
                                                  self.norm_constrained
                                                 ).constraints

        self._bounds = CorrelationBounds(self.cluster.num_clusters,
                                         len(self.cluster.single_point_clusters),
                                         self.cluster.disordered_correlations[self.cluster.single_point_clusters]
                                        ).sro_bounds

        assert self.cluster.vmatrix_array.shape == (len(self._multconfig_kb), len(self._mults_eci))

        self._bh_accept = BasinHoppingBounds(self.cluster.vmatrix_array)
        self._bh_step = BasinHoppingStep(self.cluster, self.stepsize)
        self._bh_callback = None #TODO

    def fit(self: ShortRangeOrderOptimization_Basinhopping,
            params: dict,
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

        minimizer_kwargs = {'args' : (self._mults_eci, self._multconfig_kb, self.cluster.vmatrix_array, self._vrhologrho, temperature,),
                            'method': 'trust-constr',
                            'options': self.options,
                            'jac': self._dF, 'hess': self._d2F,
                            'constraints' : self._constraints,
                            'bounds': self._bounds,
                           }

        try:
            result = basinhopping(self._F,
                                  x0 = result_correlations,
                                  minimizer_kwargs = minimizer_kwargs,
                                  niter = self.num_trials,
                                  T = self.bh_temperature,
                                  stepsize = self.stepsize,
                                  take_step = self._bh_step,
                                  accept_test = self._bh_accept,
                                  niter_success = self.early_stopping_count,
                                  interval = self.bh_interval,
                                  seed = self._seed,
                                  callback = self._bh_callback,
                                  disp=self.print_output
                                 )
        except OptimizeWarning as opt_warn:
            print(opt_warn)

        try:
            assert result.lowest_optimization_result.constr_violation < self.constr_tol
        except AssertionError:
            print('Major Constrain violation')

        result_value = result.fun
        result_grad = result.lowest_optimization_result.grad.copy()
        result_correlations = result.x.copy()
        result_constr_viol = result.lowest_optimization_result.constr_violation

        return (result_value, result_correlations, result_grad, result_constr_viol)
