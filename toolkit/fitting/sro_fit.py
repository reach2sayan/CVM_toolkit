from __future__ import annotations
from typing import Callable
import sys
import json
import inspect

import numpy as np
from scipy.optimize import curve_fit
from sympy.parsing.sympy_parser import (
    parse_expr,
    standard_transformations,
    implicit_multiplication,
)

eV2J = 96491.5666370759
#return popt, pcov, return_line, return_line_replaced, func.__code__.co_varnames[1:]
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
    def __init__(self: SROCorrectionModel,
                 func: Callable[[...],float],
                 num_str_atoms: int,
                 in_Joules: bool = True,
                 method: str = 'dogbox',
                 print_output: bool = True,
                ) -> None:

        self._p0 = None
        self.func = func
        self.num_str_atoms = num_str_atoms
        self.in_Joules = in_Joules
        self.method = method
        self.print_output = print_output

        self.popt = None
        self.pcov = None

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
    def return_line_replaced(self):
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

        if self.print_output:
            print(f'Output to func:\n{return_line_replaced}')
            return return_line_replaced

    def sro_fit(self, results):
        try:
            with open(results, 'r', encoding='utf-8') as fhandle:
                data = json.load(fhandle)
        except FileNotFoundError:
            sys.exit(f"Data file {results.split('/')[-1]} not found...")

        xdata = np.array([item.get('temperature') for item in data if item.get('temperature') != 0])
        ydata = np.array([item.get('F_cvm') - item.get('F_rnd') for item in data if item.get('temperature') != 0])*self.num_str_atoms
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
            np.savetxt('sro_fitting_calculated.out',np.hstack((xcont, ycont)))
            np.savetxt('sro_fitting_data.out',np.hstack((xdata[:,np.newaxis], ydata[:,np.newaxis])))

