import numpy as np

class BasinHoppingBounds:
    """
    Class to constrain the trial correlations of Basin Hopping
    """
    def __init__(self, all_vmat: np.ndarray) -> None:
        self.all_vmat = all_vmat

    def __call__(self, **kwargs: dict) -> bool:
        rho = self.all_vmat @ kwargs['x_new']
        return np.all((rho >= 0.0) & (rho <= 1.0))

class BasinHoppingStep:
    """
    Class to define the step in a Basin Hopping algorithm
    """
    def __init__(self,
                 cluster: dict,
                 stepsize: float = 0.01,
                 seed: int = 42,
                ) -> None:

        self.cluster = cluster
        self.stepsize = stepsize
        self._rng = np.random.default_rng(seed)

    def __call__(self, x: np.ndarray) -> None:

        jitter = np.array([0,
                           *[0]*len(self.cluster.single_point_clusters),
                           *self._rng.normal(0,
                                             self.stepsize,
                                             self.cluster.num_clusters -
                                             len(self.cluster.single_point_clusters) - 1
                                            )
                          ]
                         )
        return x + jitter
