import argparse

def SRO_argument_parser():
    """
    Argument Parse for SRO correction code
    """

    description = r"""
$$$$$$\  $$\    $$\ $$\      $$\        $$$$$$\             $$\     $$\               $$\                     $$\     $$\
$$  __$$\ $$ |   $$ |$$$\    $$$ |      $$  __$$\            $$ |    \__|              \__|                    $$ |    \__|
$$ /  \__|$$ |   $$ |$$$$\  $$$$ |      $$ /  $$ | $$$$$$\ $$$$$$\   $$\ $$$$$$\$$$$\  $$\ $$$$$$$$\ $$$$$$\ $$$$$$\   $$\  $$$$$$\  $$$$$$$\
$$ |      \$$\  $$  |$$\$$\$$ $$ |      $$ |  $$ |$$  __$$\\_$$  _|  $$ |$$  _$$  _$$\ $$ |\____$$  |\____$$\\_$$  _|  $$ |$$  __$$\ $$  __$$\
$$ |       \$$\$$  / $$ \$$$  $$ |      $$ |  $$ |$$ /  $$ | $$ |    $$ |$$ / $$ / $$ |$$ |  $$$$ _/ $$$$$$$ | $$ |    $$ |$$ /  $$ |$$ |  $$ |
$$ |  $$\   \$$$  /  $$ |\$  /$$ |      $$ |  $$ |$$ |  $$ | $$ |$$\ $$ |$$ | $$ | $$ |$$ | $$  _/  $$  __$$ | $$ |$$\ $$ |$$ |  $$ |$$ |  $$ |
\$$$$$$  |   \$  /   $$ | \_/ $$ |       $$$$$$  |$$$$$$$  | \$$$$  |$$ |$$ | $$ | $$ |$$ |$$$$$$$$\\$$$$$$$ | \$$$$  |$$ |\$$$$$$  |$$ |  $$ |
 \______/     \_/    \__|     \__|       \______/ $$  ____/   \____/ \__|\__| \__| \__|\__|\________|\_______|  \____/ \__| \______/ \__|  \__|
                                                  $$ |
                                                  $$ |
                                                  \__|
$$$$$$$$\                  $$\ $$\       $$\   $$\
\__$$  __|                 $$ |$$ |      \__|  $$ |
   $$ | $$$$$$\   $$$$$$\  $$ |$$ |  $$\ $$\ $$$$$$\
   $$ |$$  __$$\ $$  __$$\ $$ |$$ | $$  |$$ |\_$$  _|
   $$ |$$ /  $$ |$$ /  $$ |$$ |$$$$$$  / $$ |  $$ |
   $$ |$$ |  $$ |$$ |  $$ |$$ |$$  _$$<  $$ |  $$ |$$\
   $$ |\$$$$$$  |\$$$$$$  |$$ |$$ | \$$\ $$ |  \$$$$  |
   \__| \______/  \______/ \__|\__|  \__|\__|   \____/


written by Sayan Samanta and Axel van de Walle @ Brown University, RI, USA
"""
    parser = argparse.ArgumentParser(prog='sro_correction',
                                     #description='CVM SRO Error Correction Code by Sayan Samanta and Axel van de Walle',
                                     description=description,
                                     epilog='The code uses scipy heavily for all numerical optimization. Thanks guys.',
                                     #formatter_class=argparse.ArgumentDefaultsHelpFormatter,
                                     formatter_class=argparse.RawDescriptionHelpFormatter
                                    )

    opt_params = parser.add_argument_group("Parameters related to CVM Optimisation")
    sro_fit_params = parser.add_argument_group("Parameters related to SRO fitting")
    clus_params = parser.add_argument_group("Parameters related to CVM description")

    parser.add_argument('--seed', '-s',
                        default=42,
                        type=float,
                        help='set numpy random seed (for reproducible experiments) [default: %(default)s]'
                       )
    parser.add_argument('--disp', '-d',
                        action='store_true',
                        default='False',
                        help="Flag to turn intermediate outputs on [default: %(default)s]",
                       )
    parser.add_argument('--log', '-l',
                        default='log.out',
                        help="Filename for the log file [default: %(default)s]"
                       )
    parser.add_argument('--out', '-o',
                        default='result.csv',
                        help='Indicates the name of the output file [default: %(default)s]'
                       )
    parser.add_argument('--toscreen','-ts',
                        action='store_true',
                        help='Flag to output to screen in addition to log file [default: %(default)s]'
                       )

    clus_params.add_argument('--eci','-e',
                             default='eci.out',
                             help="file containing ECI's [default: %(default)s]",
                            )
    clus_params.add_argument('--vmat','-vm',
                             default='vmat.out',
                             help="file containing the vmatrix [default: %(default)s]",
                            )
    clus_params.add_argument('--clusters','-cl',
                             default='clusters.out',
                             help="file containing cluster description [default: %(default)s]" ,
                            )
    clus_params.add_argument('--maximal_clusters','-mcl',
                             default='maxclus.in',
                             help="file contain the maximal cluster coordinates [default: %(default)s]",
                            )
    clus_params.add_argument('--clustermult','-clm',
                             default='clusmult.out',
                             help="file containing cluster multiplicities [default: %(default)s]",
                            )
    clus_params.add_argument('--kikuchi_barker','-kb',
                             default='configkb.out',
                             help="file containing the kikuchi-barker coefficients [default: %(default)s]",
                            )
    clus_params.add_argument('--configmult','-com',
                             default='configmult.out',
                             help="file containing cluster configuration multiplicities [default: %(default)s]",
                            )
    clus_params.add_argument('--config','-co',
                             default='config.out',
                             help="file containing cluster configuration descriptions [default: %(default)s]",
                            )
    clus_params.add_argument('--lat', '-la',
                             default='lat.in',
                             help="contains the lattice description of the phase [default: %(default)s]"
                            )

    sro_fit_params.add_argument('--Tmin','-Tl',
                                default=100,
                                type=int,
                                help="Minimum temperature for SRO correction [default: %(default)s]"
                               )
    sro_fit_params.add_argument('--Tmax','-Tm',
                                default=2000,
                                type=int,
                                help="Maximum temperature for SRO correction [default: %(default)s]"
                               )
    sro_fit_params.add_argument('--Tstep','-Ts',
                                default=100,
                                type=int,
                                help="Temperature increment for SRO correction [default: %(default)s]"
                               )
    sro_fit_params.add_argument('--inJoules','-ij',
                                action='store_true',
                                default=False,
                                help="Sets energy units to Joules/mol [default: %(default)s]"
                               )
    sro_fit_params.add_argument('--sro_method','-sm',
                                default='lm',
                                help="method of fitting SRO correction model [default: %(default)s]",
                               )
    sro_fit_params.add_argument('--coeff_out','-cout',
                                default='sro_coeffs.out',
                                help="filename to store coefficients of SRO correction model [default: %(default)s]",
                                )
    sro_fit_params.add_argument('--coeff_in','-cin',
                                default='sro_coeffs.in',
                                help="filename containing initial coefficient for fitting SRO correction model [default: %(default)s]",
                                )
    sro_fit_params.add_argument('--fit_correction_only','-fsr',
                                action='store_true',
                                default=False,
                                help="Flag to skip SRO correction optimisation but fit pre-existing data to SRO correction model [default: %(default)s]",
                                )

    opt_params.add_argument('--norm_constraint','-nc',
                            action='store_true',
                            default=False,
                            help="Flag to enable norm constraint [default: %(default)s]"
                            )
    opt_params.add_argument('--constr_tol','-ctol',
                            type=float,
                            default=1e-8,
                            help="Value for maximum allowed constraints violation [default: %(default)s]"
                            )
    opt_params.add_argument('--maxiter','-mit',
                            type=int,
                            default=5000,
                            help="Maximum no. of iterations for local minimizer [default: %(default)s]"
                           )
    opt_params.add_argument('--maxiter_linprog','-mitlin',
                            type=int,
                            default=500000,
                            help="Maximum no. of iterations for linear programming search [default: %(default)s]"
                           )
    opt_params.add_argument('--global_iterations','-git',
                            type=int,
                            default=50,
                            help="No. of global optimisations search steps [default: %(default)s]"
                           )
    opt_params.add_argument('--xtol','-x',
                            type=float,
                            default=1e-8,
                            help="Tolerance for termination of local minimizer by the change of correlations.\
                            The algorithm will terminate when ``tr_radius < xtol``, where\
                            ``tr_radius`` is the radius of the trust region used in the algorithm [default: %(default)s]"
                           )
    opt_params.add_argument('--gtol','-g',
                            type=float,
                            default=1e-8,
                            help="The algorithm will terminate when both the infinity norm (i.e., max abs value)\
                            of the Lagrangian gradient and the constraint violation\
                            are smaller than ``gtol`` [default: %(default)s]"
                           )
    opt_params.add_argument('--barrier_tol','-btol',
                            type=float,
                            default=1e-8,
                            help="Threshold on the barrier parameter for the algorithm termination [default: %(default)s]"
                           )
    opt_params.add_argument('--initial_tr_radius','-itr',
                            type=float,
                            default=1,
                            help="Initial trust radius. It reflects the trust the algorithm puts in the\
                            local approximation of the optimization problem. For an accurate local approximation\
                            the trust-region should be large and for an approximation valid\
                            only close to the current point it should be a small one [default: %(default)s]"
                           )
    opt_params.add_argument('--initial_constr_penalty','-icp',
                            type=float,
                            default=1,
                            help="Initial Constraint Penalty. The penalty parameter is used for\
                            balancing the requirements of decreasing the objective function\
                            and satisfying the constraints [default: %(default)s]"
                           )
    opt_params.add_argument('--fit_ordered_only','-foo',
                            action='store_true',
                            default=False,
                            help="Flag to find find_ordered state only and exit [default: %(default)s]"
                            )
    opt_params.add_argument('--method_linprog','-linm',
                            default='highs',
                            help="Method of linear programming for finding ordered correlations [default: %(default)s]",
                            )
    opt_params.add_argument('--basinhopping','-bh',
                            action='store_true',
                            default='False',
                            help="Flag to switch on basinhopping [default: %(default)s]",
                           )
    opt_params.add_argument('--verbose', '-v', action='count', default=0,
                            help="Indicate the verbosity of the fit [default: %(default)s]",
                            )
    opt_params.add_argument('--approx_deriv','-ad',
                            action='store_true',
                            default=False,
                            help="Flag to enable estimation of derivatives [default: %(default)s]",
                            )
    opt_params.add_argument('--earlystop','-es',
                            default=20,
                            type=int,
                            help="Number of steps to break out of trials if no new minima has been found [default: %(default)s]",
                            )
    opt_params.add_argument('--initial_stepsize','-is',
                            default=0.1,
                            type=float,
                            help="Initial stepsize of the basinhopping algorithm from the disordered phase [default: %(default)s]",
                            )

    args = parser.parse_args()
    return args
