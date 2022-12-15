"""
Optimiser module for SRO Correction
"""

from __future__ import annotations
from typing import Callable, Type, Any
from dataclasses import dataclass, field
from abc import ABC, abstractmethod
import numpy as np
from scipy.optimize import BFGS

from toolkit.cluster.Cluster import Cluster
from toolkit.functions.energyfunctions import F, F_hessian, F_jacobian

@dataclass(kw_only=True, order=False, eq=False, repr=False)
class ClusterOptimizer(ABC):
    """
    Base Optimiser Class
    """

    cluster: Type[Cluster]

    print_output: bool = True
    approx_deriv: bool = True
    num_trials: int = 100
    constr_tol: float = 1e-10
    early_stopping_count: int = 100
    norm_constrained: bool = False

    _bounds: Any = field(init=False)
    _constraints: list = field(init=False,default_factory=list)
    _mults_eci: np.ndarray = field(init=False)
    _multconfig_kb: np.ndarray = field(init=False)
    _vrhologrho: Callable[[np.ndarray],np.ndarray] = np.vectorize(lambda rho: rho * np.log(np.abs(rho)))
    _seed: int = 42
    _F: Callable[[np.ndarray,
                 np.ndarray,
                 np.ndarray,
                 np.ndarray,
                 Callable[[np.ndarray],np.ndarray],
                 float,
                ], float] = F
    _dF: Callable[...,np.ndarray] = field(init=False)
    _d2F: Callable[...,np.ndarray] = field(init=False)

    def __post_init__(self) -> None:

        self._mults_eci = self.cluster.clusmult_array * self.cluster.eci_array
        self._multconfig_kb = self.cluster.configmult_array * self.cluster.kb_array
        
        if self.approx_deriv:
            print('Approximating the derivatives - Jacobian : a 3-point finite diffrence scheme, Hessian : BFGS')
            self._dF = '3-point'
            self._d2F = BFGS()
        else:
            self._dF = F_jacobian
            self._d2F = F_hessian

    @abstractmethod
    def fit(self, **kwargs: Any):
        pass
