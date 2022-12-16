from __future__ import annotations
from dataclasses import dataclass, field
from typing import Callable
import sys
import os
import json
import inspect

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from sympy.parsing.sympy_parser import (
    parse_expr,
    standard_transformations,
    implicit_multiplication,
)
from toolkit.io.SROResults import SROResults

eV2J = 96491.5666370759
#return popt, pcov, return_line, return_line_replaced, func.__code__.co_varnames[1:]
@dataclass(kw_only=True, eq=False, order=False)
class SROCorrectionModel:
    """
    Function to fit SRO correction as a function of T
    Inputs:
        func - the fitting function, the first parameter is the independent parameter
        results - json file containing the results of a CVM optimisation
        coeff_in - Optional parameter for initial guesses
        method - method of the fit
        verbose - verbosity
    Output:
        popt - Best fit params
        pcov - covariance of the parameters
    """
    func: Callable[[...],float]
    num_str_atoms: int

    in_Joules: bool = True
    method: str = 'dogbox'
    print_output: bool = True

    _data_output_fname: str = 'sro_correction_data.out'
    _pred_output_fname: str = 'sro_correction_predicted.out'
    popt: np.ndarray = None
    pcov: np.ndarray = None
    _p0: np.ndarray  = None
    _data: pd.DataFrame = None

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, data):
        if isinstance(data, str):
            self._data = pd.read_csv(data)
        elif isinstance(data, pd.DataFrame):
            self._data = data

    @property
    def p0(self: SROCorrectionModel) -> None:
        return self._p0

    @p0.setter
    def p0(self: SROCorrectionModel,
           coeff_in_fname: str
          ) -> None:
        try:
            self._p0 = np.loadtxt(coeff_in_fname)
        except OSError:
            print('File containing initial coefficients not found..taking defaults...')
            self._p0 = None

    @property
    def return_line(self):
        for line in inspect.getsource(self.func).split('\n'):
            if 'return' in line:
                return_line = line.replace('return','').strip()
        return str(return_line)

    @property
    def sro_function(self):
        return_line_replaced = self.return_line
        for index, varname in enumerate(self.func.__code__.co_varnames[1:]):
            return_line_replaced = return_line_replaced.replace(varname,f'{self.popt[index]}')

        func_replace = {'np.exp':'exp',
                        'np.abs':'ABS',
                        'np.sin':'SIN'}

        for old_item, new_item in func_replace.items():
            return_line_replaced = return_line_replaced.replace(old_item,new_item)

        return_line_replaced = str(parse_expr(return_line_replaced,
                                              transformations=standard_transformations+(implicit_multiplication,)).expand(basic=True))
        return_line_replaced = return_line_replaced.replace('/T','*(T**(-1))')
        return_line_replaced = return_line_replaced.replace(' ','')
        return_line_replaced = return_line_replaced.replace('++','+')
        return_line_replaced = return_line_replaced.replace('+-','-')
        return_line_replaced = return_line_replaced.replace('-+','-')
        return_line_replaced = return_line_replaced.replace('--','+')
        return_line_replaced = return_line_replaced.replace('exp','EXP')

        return return_line_replaced

    def fit(self: SROCorrectionModel) -> None:

        xdata = self._data[self._data.temperature != 0].temperature.values
        ydata = self._data[self._data.temperature != 0].apply(lambda x: x.F_opt - x.F_rnd, axis=1).values

        ydata *= self.num_str_atoms
        if self.in_Joules:
            ydata = ydata*eV2J

        self.popt, self.pcov = curve_fit(f=self.func,
                                         xdata=xdata, ydata=ydata,
                                         p0=self.p0,
                                         method=self.method,
                                         maxfev=10_000_000
                                        )

        if self.print_output:
            xcont = np.linspace(min(xdata), max(xdata), 1000)[:, np.newaxis]
            ycont = self.func(xcont,*self.popt)
            np.savetxt(f'{self._pred_output_fname}',np.hstack((xcont, ycont)))
            np.savetxt(f'{self._data_output_fname}',np.hstack((xdata[:,np.newaxis], ydata[:,np.newaxis])))

    def plot_fit(self: SROCorrectionModel,
                 image_name: str = 'cvmfit.svg',
                 **matplotlib_kwargs,
                ) -> None:

        if matplotlib_kwargs is not None:
            plt.rc(matplotlib_kwargs)
        else:
            plt.style.use('classic')
            plt.rc('text', usetex=True)
            plt.rc('font', family='serif',weight='bold',)
            plt.rc('font', family='serif',weight='bold',)
            plt.rc('xtick', labelsize='large')
            plt.rc('ytick', labelsize='large')
        structure = os.getcwd()

        ydata = self._data[self._data.temperature != 0].apply(lambda x: x.F_opt - x.F_rnd, axis=1).values
        ydata *= self.num_str_atoms
        if self.in_Joules:
            ydata = ydata*eV2J

        xdata = self._data[self._data.temperature != 0].temperature.values
        xcont = np.linspace(min(xdata), max(xdata), 1000)[:, np.newaxis]
        ycont = self.func(xcont,*self.popt)

        plt.plot(xcont,ycont,'b-',label='fitted')
        plt.plot(xdata,ydata,'r.',label='data')
        plt.xlabel('temperature (in K)')
        plt.suptitle(f'SRO Correction in {structure.split("/")[-1]}')
        plt.grid()
        if self.in_Joules:
            plt.ylabel(r'CVM Correction (F$_{cvm}$ - F$_{rnd}$) (in J/mol)')
        else:
            plt.ylabel('SRO Correction (in eV/atom)')
        plt.legend()

        plt.savefig(f'{structure}/{image_name}')
