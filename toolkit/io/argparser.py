import argparse

def SRO_argument_parser():
    """
    Argument Parse for SRO correction code
    """
    parser = argparse.ArgumentParser(prog='sro_correction',
                                     description='CVM SRO Error Correction Code by Sayan Samanta and Axel van de Walle',
                                     epilog='The code uses scipy heavily for all numerical optimization. Thanks guys.',
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                    )

    opt_params = parser.add_argument_group("Parameters related to CVM Optimisation")
    sro_fit_params = parser.add_argument_group("Parameters related to SRO fitting")
    clus_params = parser.add_argument_group("Parameters related to CVM description")

    parser.add_argument('--seed', '-s',
                        default=42,
                        type=float,
                        help='set numpy random seed (for reproducible experiments)'
                       )
    parser.add_argument('--disp', '-d',
                        action='store_true',
                        default='False',
                        help="Flag to turn intermediate outputs on",
                       )
    parser.add_argument('--log', '-l',
                        default='log.out',
                        help="Filename for the log file "
                       )
    parser.add_argument('--out', '-o',
                        default='result.csv',
                        help='Indicates the name of the output file'
                       )
    parser.add_argument('--toscreen','-ts',
                        action='store_true',
                        help='Flag to output to screen in addition to log file'
                       )

    clus_params.add_argument('--eci','-e',
                             default='eci.out',
                             help="file containing ECI's ",
                            )
    clus_params.add_argument('--vmat','-vm',
                             default='vmat.out',
                             help="file containing the vmatrix ",
                            )
    clus_params.add_argument('--clusters','-cl',
                             default='clusters.out',
                             help="file containing cluster description",
                            )
    clus_params.add_argument('--maximal_clusters','-mcl',
                             default='maxclus.in',
                             help="file contain the maximal cluster coordinates",
                            )
    clus_params.add_argument('--clustermult','-clm',
                             default='clusmult.out',
                             help="file containing cluster multiplicities",
                            )
    clus_params.add_argument('--kikuchi_barker','-kb',
                             default='configkb.out',
                             help="file containing the kikuchi-barker coefficients ",
                            )
    clus_params.add_argument('--configmult','-com',
                             default='configmult.out',
                             help="file containing cluster configuration multiplicities",
                            )
    clus_params.add_argument('--config','-co',
                             default='config.out',
                             help="file containing cluster configuration descriptions ",
                            )
    clus_params.add_argument('--lat', '-la',
                             default='lat.in',
                             help="contains the lattice description of the phase "
                            )

    sro_fit_params.add_argument('--Tmin','-Tl',
                                default=100,
                                type=int,
                                help="Minimum temperature for SRO correction"
                               )
    sro_fit_params.add_argument('--Tmax','-Tm',
                                default=2000,
                                type=int,
                                help="Maximum temperature for SRO correction"
                               )
    sro_fit_params.add_argument('--Tstep','-Ts',
                                default=100,
                                type=int,
                                help="Temperature increment for SRO correction"
                               )
    sro_fit_params.add_argument('--inJoules','-ij',
                                action='store_true',
                                default=False,
                                help="Sets energy units to Joules/mol"
                               )
    sro_fit_params.add_argument('--sro_method','-sm',
                                default='lm',
                                help="method of fitting SRO correction model",
                               )
    sro_fit_params.add_argument('--coeff_out','-cout',
                                default='sro_coeffs.out',
                                help="filename to store coefficients of SRO correction model ",
                                )
    sro_fit_params.add_argument('--coeff_in','-cin',
                                default='sro_coeffs.in',
                                help="filename containing initial coefficient for fitting SRO correction model",
                                )
    sro_fit_params.add_argument('--fit_correction_only','-fsr',
                                action='store_true',
                                default=False,
                                help="Flag to skip SRO correction optimisation but fit pre-existing data to SRO correction model",
                                )

    opt_params.add_argument('--norm_constraint','-nc',
                            action='store_true',
                            default=False,
                            help="Flag to enable norm constraint"
                            )
    opt_params.add_argument('--constr_tol','-ctol',
                            type=float,
                            default=1e-8,
                            help="Value for maximum allowed constraints violation"
                            )
    opt_params.add_argument('--maxiter','-mit',
                            type=int,
                            default=5000,
                            help="Maximum no. of iterations for local minimizer"
                           )
    opt_params.add_argument('--maxiter_linprog','-mitlin',
                            type=int,
                            default=500000,
                            help="Maximum no. of iterations for linear programming search"
                           )
    opt_params.add_argument('--global_iterations','-git',
                            type=int,
                            default=50,
                            help="No. of global optimisations search steps"
                           )
    opt_params.add_argument('--xtol','-x',
                            type=float,
                            default=1e-8,
                            help="Tolerance for termination of local minimizer by the change of correlations.\
                            The algorithm will terminate when ``tr_radius < xtol``, where\
                            ``tr_radius`` is the radius of the trust region used in the algorithm"
                           )
    opt_params.add_argument('--gtol','-g',
                            type=float,
                            default=1e-8,
                            help="The algorithm will terminate when both the infinity norm (i.e., max abs value)\
                            of the Lagrangian gradient and the constraint violation\
                            are smaller than ``gtol``"
                           )
    opt_params.add_argument('--barrier_tol','-btol',
                            type=float,
                            default=1e-8,
                            help="Threshold on the barrier parameter for the algorithm termination."
                           )
    opt_params.add_argument('--initial_tr_radius','-itr',
                            type=float,
                            default=1,
                            help="Initial trust radius. It reflects the trust the algorithm puts in the\
                            local approximation of the optimization problem. For an accurate local approximation\
                            the trust-region should be large and for an approximation valid\
                            only close to the current point it should be a small one."
                           )
    opt_params.add_argument('--initial_constr_penalty','-icp',
                            type=float,
                            default=1,
                            help="Initial Constraint Penalty. The penalty parameter is used for\
                            balancing the requirements of decreasing the objective function\
                            and satisfying the constraints."
                           )
    opt_params.add_argument('--fit_ordered_only','-foo',
                            action='store_true',
                            default=False,
                            help="Flag to find find_ordered state only and exit."
                            )
    opt_params.add_argument('--method_linprog','-linm',
                            default='highs',
                            help="Method of linear programming for finding ordered correlations",
                            )
    opt_params.add_argument('--basinhopping','-bh',
                            action='store_true',
                            default='False',
                            help="Flag to switch on basinhopping",
                           )
    opt_params.add_argument('--verbose', '-v', action='count', default=0,
                            help="Indicate the verbosity of the fit ",
                            )
    opt_params.add_argument('--approx_deriv','-ad',
                            action='store_true',
                            default=False,
                            help="Flag to enable estimation of derivatives",
                            )
    opt_params.add_argument('--earlystop','-es',
                            default=20,
                            type=int,
                            help="Number of steps to break out of trials if no new minima has been found",
                            )
    opt_params.add_argument('--initial_stepsize','-is',
                            default=0.1,
                            type=float,
                            help="Initial stepsize of the basinhopping algorithm from the disordered phase",
                            )

    args = parser.parse_args()
    return args
