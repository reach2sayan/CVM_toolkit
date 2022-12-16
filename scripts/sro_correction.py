import os
import sys
import subprocess
from pathlib import Path
from datetime import datetime
from collections.abc import Iterable

from toolkit.io.argparser import SRO_argument_parser
from toolkit.io.SROResults import SROResults
from toolkit.cluster.Cluster import Cluster
from toolkit.optimizers.OrderedStateOptimizer import OrderedStateOptimizer
from toolkit.optimizers.CVMOptimizer import CVMOptimizer
from toolkit.fitting.SROCorrectionModel import SROCorrectionModel
from toolkit.functions.sro_model import sro_model
from toolkit.logger.Logger import Logger

import numpy as np

def custom_linspace(start: float, stop: float, step:float =1) -> Iterable[float]:
    """
    Like np.linspace but uses step instead of num
    This is inclusive to stop, so if start=1, stop=3, step=0.5
    Output is: array([1., 1.5, 2., 2.5, 3.])
    """

    return np.linspace(start, stop, int((stop - start) / step + 1))

if __name__ == '__main__':

    np.set_printoptions(suppress=True, precision=4)

    structure = os.getcwd()
    path = Path(structure)
    phase = str(path.parent.absolute())

    #parse arguments
    args = SRO_argument_parser()
    tdate = datetime.now().strftime('%d%b-%H%m')
    if args.norm_constraint:
        log_fname = f'{args.log}-cons-{tdate}'
        out_fname = f'{args.out}-cons-{tdate}'
    else:
        log_fname = f'{args.log}-nocons-{tdate}'
        out_fname = f'{args.out}-nocons-{tdate}'
    sys.stdout = Logger(sys.stdout, log_fname, args.toscreen)

    results = SROResults(phase = phase,
                         structure = structure,
                         norm_constrained = args.norm_constraint,
                        )

    if args.fit_correction_only:
        cluster = Cluster(_clusters_fname = args.clusters,
                          _eci_fname = args.eci,
                          _clustermult_fname = args.clustermult,
                          _config_fname = args.config,
                          _configmult_fname = args.configmult,
                          _kb_fname = args.kikuchi_barker,
                          _vmat_fname = args.vmat,
                          _lattice_fname = args.lat
                         )
        sro_correction_model = SROCorrectionModel(func=sro_model,
                                                  num_str_atoms = cluster.num_str_atoms,
                                                  in_Joules = args.inJoules,
                                                  print_output = args.disp
                                                 )
        sro_correction_model.data = args.out
        _ = sro_correction_model.fit()
        sro_correction_model.plot_fit()
        with open(f'{structure}/func','w',encoding='utf-8') as correction_func:
            correction_func.write(sro_correction_model.sro_function)
            print('Flag to fit SRO Correction function only. Exiting.')
            sys.exit(0)

    try:
        _ = subprocess.run(['cvmclus', f'-m={phase}/{args.maximal_clusters}', f'-l={structure}/{args.lat}'],
                           stdout=subprocess.PIPE,
                           stderr=subprocess.PIPE,
                           check=True
                           )
    except subprocess.SubprocessError as suberr:
        print('Error in generating cluster configuration files for CVM.')
        print('Continuing to run. Might work if the required files already exists, maybe created by some previous run of cvmclus')

    #Create Cluster Object
    cluster = Cluster(_clusters_fname = args.clusters,
                      _eci_fname = args.eci,
                      _clustermult_fname = args.clustermult,
                      _config_fname = args.config,
                      _configmult_fname = args.configmult,
                      _kb_fname = args.kikuchi_barker,
                      _vmat_fname = args.vmat,
                      _lattice_fname = args.lat
                     )

    #Calculate Ordered State
    options_ordered = {'disp': bool(args.verbose),
                       'maxiter': args.maxiter_linprog,
                      }

    opt_ordered = OrderedStateOptimizer(cluster = cluster,
                                        print_output = args.disp,
                                        num_trials = args.maxiter_linprog,
                                        options = options_ordered,
                                        method = args.method_linprog,
                                       )
    _ = opt_ordered.fit()

    if args.fit_ordered_only:
        print('Flag to fit ordered state only found. Exiting.')
        sys.exit(0)

    if args.basinhopping.title() == 'True':
        raise NotImplementedError('In the works...Consider using the random search')
    else:
        options = {'verbose': args.verbose,
                   'maxiter': args.maxiter,
                   'xtol': args.xtol,
                   'gtol': args.gtol,
                   'barrier_tol': args.barrier_tol,
                   'initial_tr_radius': args.initial_tr_radius,
                   'initial_constr_penalty': args.initial_constr_penalty,
                  }
        opt_sro = CVMOptimizer(cluster = cluster,
                               print_output = args.disp,
                               approx_deriv = args.approx_deriv,
                               num_trials = args.global_iterations,
                               constr_tol = args.constr_tol,
                               early_stopping_count = args.earlystop,
                               norm_constrained = args.norm_constraint,
                               options = options,
                              )

    #MAIN LOOP
    results_ = []
    for T in custom_linspace(start=args.Tmin, stop=args.Tmax, step=args.Tstep):

        opt_sro.temperature = T
        print('=' * 25)
        print(f'Optimising at temperature {opt_sro.temperature}K')

        F_ordered = opt_sro.get_energy(opt_sro.cluster.ordered_correlations)
        print(f'Ordered Correlations @ T = {T}K:')
        print(f'{opt_sro.cluster.ordered_correlations}')
        print(f'Ordered CVM Free Energy (eV/atom) @ T = {opt_sro.temperature}K: {F_ordered/opt_sro.cluster.num_lat_atoms}')

        F_disordered = opt_sro.get_energy(opt_sro.cluster.disordered_correlations)
        print(f'Disordered Correlations @ T = {T}K:')
        print(f'{opt_sro.cluster.ordered_correlations}')
        print(f'Disordered CVM Free Energy (eV/atom) @ T = {opt_sro.temperature}K: {F_disordered/opt_sro.cluster.num_lat_atoms}')

        F_sqs = opt_sro.get_energy(opt_sro.cluster.sqs_correlations)
        print(f'SQS Correlations @ T = {T}K:')
        print(f'{opt_sro.cluster.sqs_correlations}')
        print(f'SQS CVM Free Energy (eV/atom) @ T = {opt_sro.temperature}K: {F_sqs/opt_sro.cluster.num_lat_atoms}')

        opt_F, opt_correlations, opt_grad, opt_constr_viol = opt_sro.fit()

        print(f'Optimised CVM Free energy (eV/atom) @ T = {opt_sro.temperature}K\t: {opt_F/opt_sro.cluster.num_lat_atoms}')
        print(f'Optimised Correlations @ T = {opt_sro.temperature}K\t: {opt_correlations}')
        print(f'Constraint Violation: {opt_constr_viol}')
        print(f'Optimised Free energy Gradient: {np.array2string(opt_grad)}')
        print('Optimised Cluster Configuration Probabilities:')
        opt_sro.cluster.print_config_probabilities(opt_correlations)
        print('=' * 25)

        results_.append({'phase'        : opt_sro.cluster.phase.rsplit('/',maxsplit=1)[-1],
                         'structure'    : opt_sro.cluster.structure.split('/')[-1],
                         'temperature'  : opt_sro.temperature,
                         'F_sqs'        : F_sqs,
                         'F_ord'        : F_ordered,
                         'F_rnd'        : F_disordered,
                         'F_opt'        : opt_F,
                         'constr_tol'   : opt_constr_viol,
                         'correlations' : opt_correlations,
                        }
                       )
        results.result = results_
        results.save_to_file(f'{opt_sro.cluster.structure}/{out_fname}')

    sro_correction_model = SROCorrectionModel(func=sro_model,
                                              num_str_atoms = opt_sro.cluster.num_str_atoms,
                                              in_Joules = args.inJoules,
                                              print_output = args.disp
                                             )
    sro_correction_model.data = results.result
    _ = sro_correction_model.fit()
    sro_correction_model.plot_fit()
    if args.disp:
        if args.inJoules:
            print('SRO Correction function (in J/mol):')
        else:
            print('SRO Correction function (in eV/atom):')
        print(sro_correction_model.sro_function)
    with open(f'{opt_sro.cluster.structure}/func','w',encoding='utf-8') as correction_func:
        correction_func.write(sro_correction_model.sro_function)
