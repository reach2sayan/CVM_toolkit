CVM SRO Error Correction Code by Sayan Samanta and Axel van de Walle

usage: sro_correction [-h] [--seed SEED] [--disp]
                      [--log_level {INFO,DEBUG,WARNING,ERROR}] [--log LOG]
                      [--eci ECI] [--vmat VMAT] [--clusters CLUSTERS]
                      [--maximal_clusters MAXIMAL_CLUSTERS]
                      [--clustermult CLUSTERMULT] [--kb KB]
                      [--configmult CONFIGMULT] [--config CONFIG] [--lat LAT]
                      [--Tmin TMIN] [--Tmax TMAX] [--Tstep TSTEP] [--inJoules]
                      [--sro_method SRO_METHOD] [--coeff_out COEFF_OUT]
                      [--coeff_in COEFF_IN] [--skip_sro_fit] [--fit_sro_only]
                      [--norm_constraint] [--constr_tol CONSTR_TOL]
                      [--maxiter MAXITER] [--maxiter_linprog MAXITER_LINPROG]
                      [--global_iterations GLOBAL_ITERATIONS] [--xtol XTOL]
                      [--gtol GTOL] [--barrier_tol BARRIER_TOL]
                      [--initial_tr_radius INITIAL_TR_RADIUS]
                      [--initial_constr_penalty INITIAL_CONSTR_PENALTY]
                      [--fit_ordered_only] [--method_linprog METHOD_LINPROG]
                      [--basinhopping] [--verbose] [--approx_deriv]
                      [--earlystop EARLYSTOP]
                      [--initial_stepsize INITIAL_STEPSIZE] [--out OUT]


options:
  -h, --help            show this help message and exit
  --seed SEED           set numpy random seed (for reproducible experiments)
                        (default: 42)
  --disp                Flag to turn intermediate outputs on (default: False)
  --log_level {INFO,DEBUG,WARNING,ERROR}
                        Flag to set output level (default: INFO)
  --log LOG             Filename for the log file (default: log.out)

Parameters related to CVM Optimisation:
  --norm_constraint     Flag to enable norm constraint (default: False)
  --constr_tol CONSTR_TOL
                        Value for maximum allowed constraints violation
                        (default: 1e-08)
  --maxiter MAXITER     Maximum no. of iterations for local minimizer
                        (default: 5000)
  --maxiter_linprog MAXITER_LINPROG
                        Maximum no. of iterations for linear programming
                        search (default: 500000)
  --global_iterations GLOBAL_ITERATIONS
                        No. of global optimisations search steps (default: 50)
  --xtol XTOL           Tolerance for termination of local minimizer by the
                        change of correlations. The algorithm will terminate
                        when ``tr_radius < xtol``, where ``tr_radius`` is the
                        radius of the trust region used in the algorithm
                        (default: 1e-08)
  --gtol GTOL           The algorithm will terminate when both the infinity
                        norm (i.e., max abs value) of the Lagrangian gradient
                        and the constraint violation are smaller than ``gtol``
                        (default: 1e-08)
  --barrier_tol BARRIER_TOL
                        Threshold on the barrier parameter for the algorithm
                        termination. (default: 1e-08)
  --initial_tr_radius INITIAL_TR_RADIUS
                        Initial trust radius. It reflects the trust the
                        algorithm puts in the local approximation of the
                        optimization problem. For an accurate local
                        approximation the trust-region should be large and for
                        an approximation valid only close to the current point
                        it should be a small one. (default: 1)
  --initial_constr_penalty INITIAL_CONSTR_PENALTY
                        Initial Constraint Penalty. The penalty parameter is
                        used for balancing the requirements of decreasing the
                        objective function and satisfying the constraints.
                        (default: 1)
  --fit_ordered_only    Flag to find find_ordered state only and exit.
                        (default: False)
  --method_linprog METHOD_LINPROG
                        Method of linear programming for finding ordered
                        correlations (default: highs)
  --basinhopping        Flag to switch on basinhopping (default: False)
  --verbose, -v         Indicate the verbosity of the fit (default: 0)
  --approx_deriv        Flag to enable estimation of derivatives (default:
                        False)
  --earlystop EARLYSTOP
                        Number of steps to break out of trials if no new
                        minima has been found (default: 20)
  --initial_stepsize INITIAL_STEPSIZE
                        Initial stepsize of the basinhopping algorithm from
                        the disordered phase (default: 0.1)
  --out OUT             Indicates the name of the output file (default:
                        result.json)

Parameters related to SRO fitting:
  --Tmin TMIN           Minimum temperature for SRO correction (default: 100)
  --Tmax TMAX           Maximum temperature for SRO correction (default: 2000)
  --Tstep TSTEP         Temperature increment for SRO correction (default:
                        100)
  --inJoules            Sets energy units to Joules/mol (default: False)
  --sro_method SRO_METHOD
                        method of fitting SRO correction model (default: lm)
  --coeff_out COEFF_OUT
                        filename to store coefficients of SRO correction model
                        (default: sro_coeffs.out)
  --coeff_in COEFF_IN   filename containing initial coefficient for fitting
                        SRO correction model (default: sro_coeffs.in)
  --skip_sro_fit        Flag to skip fitting SRO correction model (default:
                        False)
  --fit_sro_only        Flag to skip SRO correction optimisation but fit pre-
                        existing data to SRO correction model (default: False)

Parameters related to CVM description:
  --eci ECI             file containing ECI's (default: eci.out)
  --vmat VMAT           file containing the vmatrix (default: vmat.out)
  --clusters CLUSTERS   file containing cluster description (default:
                        clusters.out)
  --maximal_clusters MAXIMAL_CLUSTERS
                        file contain the maximal cluster coordinates (default:
                        maxclus.in)
  --clustermult CLUSTERMULT
                        file containing cluster multiplicities (default:
                        clusmult.out)
  --kb KB               file containing the kikuchi-barker coefficients
                        (default: configkb.out)
  --configmult CONFIGMULT
                        file containing cluster configuration multiplicities
                        (default: configmult.out)
  --config CONFIG       file containing cluster configuration descriptions
                        (default: config.out)
  --lat LAT, -l LAT     contains the lattice description of the phase
                        (default: lat.in)

The code uses scipy heavily for all numerical optimization. Thanks guys.
