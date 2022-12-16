from __future__ import annotations
from dataclasses import dataclass, field, InitVar
from functools import cached_property
import itertools
import subprocess
import os
from pathlib import Path
import numpy as np
from toolkit.io.atatio import read_clusters, read_kbcoeffs, read_configmult, read_clustermult, read_eci, read_configs, read_vmatrix

EPSILON = 1e-2

@dataclass(kw_only=True, order=False, eq=False,)
class Cluster:
    """
    Class to hold Cluster Description
    """

    _clusters_fname: InitVar[str] = field(default='clusters.out', repr=False)
    _eci_fname: InitVar[str] = field(default='eci.out', repr=False)
    _clustermult_fname: InitVar[str] = field(default='clusmult.out', repr=False)
    _config_fname: InitVar[str] = field(default='config.out', repr=False)
    _configmult_fname: InitVar[str] = field(default='configmult.out', repr=False)
    _kb_fname: InitVar[str] = field(default='configkb.out', repr=False)
    _vmat_fname: InitVar[str] = field(default='vmat.out', repr=False)

    _lattice_fname: str = field(default='lat.in')
    _structure_fname: str = field(default='str.in')
    _sqs_structure_fname: np.ndarray = field(default='str_relax.out')
    _input_structure_fname: np.ndarray = field(default='str.in')

    _ordered_correlations: np.ndarray = None

    clusters: dict = field(init=False)
    kb: dict = field(init=False)
    configmult: dict = field(init=False)
    clustermult: dict = field(init=False)
    configs: dict = field(init=False)
    vmat: dict = field(init=False)
    eci: dict = field(init=False)
    phase: str = field(init=False)
    structure: str = field(init=False)

    def __post_init__(self: Cluster,
                      _clusters_fname: str,
                      _eci_fname: str,
                      _clustermult_fname: str,
                      _config_fname: str,
                      _configmult_fname: str,
                      _kb_fname: str,
                      _vmat_fname: str,
                     ) -> None:

        self.structure = os.getcwd()
        path_ = Path(self.structure)
        self.phase = str(path_.parent.absolute())

        self.clusters = read_clusters(f'{self.structure}/{_clusters_fname}')
        self.kb = read_kbcoeffs(f'{self.structure}/{_kb_fname}')
        self.clustermult = read_clustermult(f'{self.structure}/{_clustermult_fname}')
        self.configmult = read_configmult(f'{self.structure}/{_configmult_fname}')
        self.configs = read_configs(f'{self.structure}/{_config_fname}')
        self.vmat = read_vmatrix(f'{self.structure}/{_vmat_fname}')
        self.eci = read_eci(f'{self.structure}/{_eci_fname}')

        #TODO this part stinks but is useful now. To be removed later.
        #fix the lattice scale difference between input geometry and relax geometry
        try:
            with open(f"{self.structure}/{self._input_structure_fname}", 'r') as str_in:
                strin_lines = str_in.readlines()
        except FileNotFoundError:
            print(f'{self.structure}/{self._input_structure_fname} does not exists')
        try:
            with open(f"{self.structure}/{self._sqs_structure_fname}", 'r') as str_out:
                strout_lines = str_out.readlines()
        except FileNotFoundError:
            print(f'{self.structure}/{self._sqs_structure_fname} does not exists')

        if len(strin_lines[0].split(' ')) > 3:
            strout_lines[0:4] = strin_lines[0:4]
            strout_lines.pop(4)
            strout_lines.pop(5)
        else:
            strout_lines[0:6] = strin_lines[0:6]

        with open(f'{self.structure}/{self._sqs_structure_fname}_temp','w',encoding='utf-8') as scaled_relax:
            for line in strout_lines:
                scaled_relax.write(line)

    @property
    def sqs_correlations(self):
        corr_sqs_ = subprocess.run(['corrdump', '-c', f'-cf={self.structure}/{self._clusters_fname}', f'-s={self.structure}/{self._sqs_structure_fname}_temp', f'-l={self.structure}/{self._lattice_fname}'],
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  check=True
                                 )
        # convert from bytes to string list
        corr_sqs_ = corr_sqs_.stdout.decode('utf-8').split('\t')[:-1]
        corr_sqs_ = np.array(corr_sqs_, dtype=np.float32)  # convert to arrays

        return corr_sqs_

    @property
    def ordered_correlations(self: Cluster) -> np.ndarray:
        return self._ordered_correlations

    @ordered_correlations.setter
    def ordered_correlations(self: Cluster, correlations: np.ndarray) -> np.ndarray:
        self._ordered_correlations = correlations

    @property
    def order2disoder_distance(self: Cluster) -> float:
        return np.linalg.norm(self._ordered_correlations - self.disordered_correlations)

    @cached_property
    def num_lat_atoms(self: Cluster) -> int:
        try:
            with open(f"{self.structure}/{self._lattice_fname}", 'r') as lat_out:
                lat_lines = lat_out.readlines()
        except FileNotFoundError:
            print(f'{self.structure}/{self._lattice_fname} does not exists')

        if len(lat_lines[0].split(' ')) > 3:
            return len(lat_lines) - 4
        else:
            return len(lat_lines) - 6

    @cached_property
    def num_str_atoms(self):
        try:
            with open(f"{self.structure}/{self._input_structure_fname}", 'r') as str_out:
                str_lines = str_out.readlines()
        except FileNotFoundError:
            print(f'{self.structure}/{self._input_structure_fname} does not exists')

        if len(str_lines[0].split(' ')) > 3:
            return len(str_lines) - 4
        else:
            return len(str_lines) - 6

    @cached_property
    def num_configs(self: Cluster) -> int:
        return len(self.configmult)

    @cached_property
    def single_point_clusters(self: Cluster) -> list:
        return [cluster_idx for cluster_idx, cluster in self.clusters.items() if cluster['type'] == 1]

    @cached_property
    def num_clusters(self: Cluster) -> int:
        return len(self.clusters)

    @cached_property
    def clusmult_array(self: Cluster) -> np.ndarray:
        return np.array(list(self.clustermult.values()))

    @cached_property
    def eci_array(self: Cluster) -> np.ndarray:
        return np.array(list(self.eci.values()))

    @cached_property
    def configmult_array(self: Cluster) -> np.ndarray:
        return np.array(list(itertools.chain.from_iterable(list(self.configmult.values()))))

    @cached_property
    def kb_array(self: Cluster) -> np.ndarray:
        return np.array(list(itertools.chain.from_iterable([[kb for _ in range(
            len(self.configmult[idx]))] for idx, kb in self.kb.items()])))

    @cached_property
    def vmatrix_array(self: Cluster) -> np.ndarray:
        return np.vstack(list(self.vmat.values()))

    @cached_property
    def disordered_correlations(self: Cluster) -> np.ndarray:

        corr_rnd = subprocess.run(['corrdump', '-c', f'-cf={self.structure}/{self._clusters_fname}', f'-s={self.structure}/{self._input_structure_fname}', f'-l={self.structure}/{self._lattice_fname}', '-rnd'],
                                  stdout=subprocess.PIPE,
                                  stderr=subprocess.PIPE,
                                  check=True
                                 )
        # convert from bytes to string list
        corr_rnd = corr_rnd.stdout.decode('utf-8').split('\t')[:-1]
        corr_rnd = np.array(corr_rnd, dtype=np.float32)  # convert to arrays

        return corr_rnd

    def check_correlation_validity(self: Cluster,
                                   correlations: np.ndarray,
                                  ) -> bool:
        try:
            for vmat in self.vmat.values():
                rho = vmat @ correlations
                assert np.all([(x >= 0.0 - EPSILON) and (x <= 1.0 + EPSILON) for x in rho])
        except AssertionError:
            return False
        return True

    def print_correlations_to_file(self: Cluster,
                                   correlations: np.ndarray,
                                   filename: str = 'correlations.out'
                                  ) -> None:
        np.savetxt(filename, correlations)

    def print_config_probabilities(self: Cluster,
                                   correlations: np.ndarray,
                                  ) -> str:
        for vmat in self.vmat.values():
            print(np.array2string(vmat@correlations))
        return ''

    def print_config_probabilities_to_file(self: Cluster,
                                           correlations: np.ndarray,
                                           filename: str = 'rho.out'
                                          ) -> None:
        with open(filename, 'w', encoding='utf-8') as frho:
            for vmat in self.vmat.values():
                frho.write(f'{" ".join(map(str,vmat@correlations))}\n')

    def distance_from_ordered(self: Cluster, correlations: np.ndarray) -> float:
        return np.linalg.norm(self._ordered_correlations - correlations)

    def distance_from_disordered(self: Cluster, correlations: np.ndarray) -> float:
        return np.linalg.norm(self.disordered_correlations - correlations)

    def __repr__(self: Cluster) -> str:

        print('\n===========================')
        print('Cluster Description:')
        print('\nClusters:')
        print('{0:<8s}|{1:<8s}|{2:<8s}|{3:<8s}'.format('Index','Type', 'Mult','Radius'))
        for idx, cluster in self.clusters.items():
            assert self.clustermult[idx] == cluster['mult']
            print('{0:<8d}|{1:<8d}|{2:<12d}|{3:<8f}'.format(idx, cluster['type'], cluster['mult'], cluster['length']))

        print('\nConfigs:')
        print('{0:<8s}|{1:<8s}'.format('Index','No. of subconfigs'))
        for idx, config in self.configs.items():
            assert len(self.configmult[idx]) == config['num_of_subclus']
            print('{0:<8d}|{1:<8d}'.format(idx, config['num_of_subclus'],))

        print('\nECIs:')
        print('{0:<8s}|{1:<8s}'.format('Index','ECI'))
        for idx, eci in self.eci.items():
            print('{0:<8d}|{1:<18f}'.format(idx, eci))

        print('--------------------------------------------------')
        print('Disordered Correlations:')
        print(f'{self.disordered_correlations}')
        print('Disordered Configuration Probabilities:')
        print(f'{self.print_config_probabilities(self.disordered_correlations)}')
        print('--------------------------------------------------')

        if self.ordered_correlations is None:
            print('Ordered Correlations not calculated.')
        else:
            print('Ordered Correlations:')
            print(f'{self.ordered_correlations}')
            print('Ordered Configuration Probabilities:')
            print(f'{self.print_config_probabilities(self.ordered_correlations)}')
        print('--------------------------------------------------')

        print('SQS Correlations:')
        print(f'{self.sqs_correlations}')
        print('Ordered Configuration Probabilities:')
        print(f'{self.print_config_probabilities(self.sqs_correlations)}')
        print('--------------------------------------------------')

        print(f'Total Configurations: {self.vmatrix_array.shape[0]}')
        print(f'Total Clusters: {self.vmatrix_array.shape[1]}')
        print(f"Phase: {self.phase.rsplit('/',maxsplit=1)[-1]}")
        print(f'No. of lattice atoms: {self.num_lat_atoms}')
        print(f"Structure: {self.structure.split('/')[-1]}")
        print(f'No. of structure atoms: {self.num_str_atoms}')
        print(f'=========================================')

        return ''
