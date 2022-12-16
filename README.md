## CVM Optimization Toolkit
By Sayan Samanta and Axel van de Walle @ Brown University RI USA
---
### Usage
```
usage: sro_correction [-h] [--seed SEED] [--disp] [--log LOG] [--out OUT] [--toscreen] [--eci ECI] [--vmat VMAT] [--clusters CLUSTERS]
                      [--maximal_clusters MAXIMAL_CLUSTERS] [--clustermult CLUSTERMULT] [--kikuchi_barker KIKUCHI_BARKER] [--configmult CONFIGMULT]
                      [--config CONFIG] [--lat LAT] [--Tmin TMIN] [--Tmax TMAX] [--Tstep TSTEP] [--inJoules] [--sro_method SRO_METHOD]
                      [--coeff_out COEFF_OUT] [--coeff_in COEFF_IN] [--fit_correction_only] [--norm_constraint] [--constr_tol CONSTR_TOL]
                      [--maxiter MAXITER] [--maxiter_linprog MAXITER_LINPROG] [--global_iterations GLOBAL_ITERATIONS] [--xtol XTOL] [--gtol GTOL]
                      [--barrier_tol BARRIER_TOL] [--initial_tr_radius INITIAL_TR_RADIUS] [--initial_constr_penalty INITIAL_CONSTR_PENALTY]
                      [--fit_ordered_only] [--method_linprog METHOD_LINPROG] [--basinhopping] [--verbose] [--approx_deriv] [--earlystop EARLYSTOP]
                      [--initial_stepsize INITIAL_STEPSIZE]
```
### Parameters
#### Quick reference table
|Short    |Long                      |Default         |Description                                                                                                                                                                                                                                                                                                                                              |
|---------|--------------------------|----------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|`-h`     |`--help`                  |                |show this help message and exit                                                                                                                                                                                                                                                                                                                          |
|`-s`     |`--seed`                  |`42`            |set numpy random seed (for reproducible experiments)                                                                                                                                                                                                                                                                                                     |
|`-d`     |`--disp`                  |`False`         |Flag to turn intermediate outputs on                                                                                                                                                                                                                                                                                                                     |
|`-l`     |`--log`                   |`log.out`       |Filename for the log file                                                                                                                                                                                                                                                                                                                                |
|`-o`     |`--out`                   |`result.csv`    |Indicates the name of the output file                                                                                                                                                                                                                                                                                                                    |
|`-ts`    |`--toscreen`              |                |Flag to output to screen in addition to log file                                                                                                                                                                                                                                                                                                         |
|`-e`     |`--eci`                   |`eci.out`       |file containing ECI's                                                                                                                                                                                                                                                                                                                                    |
|`-vm`    |`--vmat`                  |`vmat.out`      |file containing the vmatrix                                                                                                                                                                                                                                                                                                                              |
|`-cl`    |`--clusters`              |`clusters.out`  |file containing cluster description                                                                                                                                                                                                                                                                                                                      |
|`-mcl`   |`--maximal_clusters`      |`maxclus.in`    |file contain the maximal cluster coordinates                                                                                                                                                                                                                                                                                                             |
|`-clm`   |`--clustermult`           |`clusmult.out`  |file containing cluster multiplicities                                                                                                                                                                                                                                                                                                                   |
|`-kb`    |`--kikuchi_barker`        |`configkb.out`  |file containing the kikuchi-barker coefficients                                                                                                                                                                                                                                                                                                          |
|`-com`   |`--configmult`            |`configmult.out`|file containing cluster configuration multiplicities                                                                                                                                                                                                                                                                                                     |
|`-co`    |`--config`                |`config.out`    |file containing cluster configuration descriptions                                                                                                                                                                                                                                                                                                       |
|`-la`    |`--lat`                   |`lat.in`        |contains the lattice description of the phase                                                                                                                                                                                                                                                                                                            |
|`-Tl`    |`--Tmin`                  |`100`           |Minimum temperature for SRO correction                                                                                                                                                                                                                                                                                                                   |
|`-Tm`    |`--Tmax`                  |`2000`          |Maximum temperature for SRO correction                                                                                                                                                                                                                                                                                                                   |
|`-Ts`    |`--Tstep`                 |`100`           |Temperature increment for SRO correction                                                                                                                                                                                                                                                                                                                 |
|`-ij`    |`--inJoules`              |                |Sets energy units to Joules/mol                                                                                                                                                                                                                                                                                                                          |
|`-sm`    |`--sro_method`            |`lm`            |method of fitting SRO correction model                                                                                                                                                                                                                                                                                                                   |
|`-cout`  |`--coeff_out`             |`sro_coeffs.out`|filename to store coefficients of SRO correction model                                                                                                                                                                                                                                                                                                   |
|`-cin`   |`--coeff_in`              |`sro_coeffs.in` |filename containing initial coefficient for fitting SRO correction model                                                                                                                                                                                                                                                                                 |
|`-fsr`   |`--fit_correction_only`   |                |Flag to skip SRO correction optimisation but fit pre-existing data to SRO correction model                                                                                                                                                                                                                                                               |
|`-nc`    |`--norm_constraint`       |                |Flag to enable norm constraint                                                                                                                                                                                                                                                                                                                           |
|`-ctol`  |`--constr_tol`            |`1e-08`         |Value for maximum allowed constraints violation                                                                                                                                                                                                                                                                                                          |
|`-mit`   |`--maxiter`               |`5000`          |Maximum no. of iterations for local minimizer                                                                                                                                                                                                                                                                                                            |
|`-mitlin`|`--maxiter_linprog`       |`500000`        |Maximum no. of iterations for linear programming search                                                                                                                                                                                                                                                                                                  |
|`-git`   |`--global_iterations`     |`50`            |No. of global optimisations search steps                                                                                                                                                                                                                                                                                                                 |
|`-x`     |`--xtol`                  |`1e-08`         |Tolerance for termination of local minimizer by the change of correlations.                        The algorithm will terminate when ``tr_radius < xtol``, where                        ``tr_radius`` is the radius of the trust region used in the algorithm                                                                                            |
|`-g`     |`--gtol`                  |`1e-08`         |The algorithm will terminate when both the infinity norm (i.e., max abs value)                        of the Lagrangian gradient and the constraint violation                        are smaller than ``gtol``                                                                                                                                           |
|`-btol`  |`--barrier_tol`           |`1e-08`         |Threshold on the barrier parameter for the algorithm termination.                                                                                                                                                                                                                                                                                        |
|`-itr`   |`--initial_tr_radius`     |`1`             |Initial trust radius. It reflects the trust the algorithm puts in the                        local approximation of the optimization problem. For an accurate local approximation                        the trust-region should be large and for an approximation valid                        only close to the current point it should be a small one.|
|`-icp`   |`--initial_constr_penalty`|`1`             |Initial Constraint Penalty. The penalty parameter is used for                        balancing the requirements of decreasing the objective function                        and satisfying the constraints.                                                                                                                                              |
|`-foo`   |`--fit_ordered_only`      |                |Flag to find find_ordered state only and exit.                                                                                                                                                                                                                                                                                                           |
|`-linm`  |`--method_linprog`        |`highs`         |Method of linear programming for finding ordered correlations                                                                                                                                                                                                                                                                                            |
|`-bh`    |`--basinhopping`          |`False`         |Flag to switch on basinhopping                                                                                                                                                                                                                                                                                                                           |
|`-v`     |`--verbose`               |`0`             |Indicate the verbosity of the fit                                                                                                                                                                                                                                                                                                                        |
|`-ad`    |`--approx_deriv`          |                |Flag to enable estimation of derivatives                                                                                                                                                                                                                                                                                                                 |
|`-es`    |`--earlystop`             |`20`            |Number of steps to break out of trials if no new minima has been found                                                                                                                                                                                                                                                                                   |
|`-is`    |`--initial_stepsize`      |`0.1`           |Initial stepsize of the basinhopping algorithm from the disordered phase                                                                                                                                                                                                                                                                                 |

#### `-h`, `--help`
show this help message and exit

#### `--seed`, `-s` (Default: 42)
set numpy random seed (for reproducible experiments)

#### `--disp`, `-d` (Default: False)
Flag to turn intermediate outputs on

#### `--log`, `-l` (Default: log.out)
Filename for the log file

#### `--out`, `-o` (Default: result.csv)
Indicates the name of the output file

#### `--toscreen`, `-ts`
Flag to output to screen in addition to log file

#### `--eci`, `-e` (Default: eci.out)
file containing ECI's

#### `--vmat`, `-vm` (Default: vmat.out)
file containing the vmatrix

#### `--clusters`, `-cl` (Default: clusters.out)
file containing cluster description

#### `--maximal_clusters`, `-mcl` (Default: maxclus.in)
file contain the maximal cluster coordinates

#### `--clustermult`, `-clm` (Default: clusmult.out)
file containing cluster multiplicities

#### `--kikuchi_barker`, `-kb` (Default: configkb.out)
file containing the kikuchi-barker coefficients

#### `--configmult`, `-com` (Default: configmult.out)
file containing cluster configuration multiplicities

#### `--config`, `-co` (Default: config.out)
file containing cluster configuration descriptions

#### `--lat`, `-la` (Default: lat.in)
contains the lattice description of the phase

#### `--Tmin`, `-Tl` (Default: 100)
Minimum temperature for SRO correction

#### `--Tmax`, `-Tm` (Default: 2000)
Maximum temperature for SRO correction

#### `--Tstep`, `-Ts` (Default: 100)
Temperature increment for SRO correction

#### `--inJoules`, `-ij`
Sets energy units to Joules/mol

#### `--sro_method`, `-sm` (Default: lm)
method of fitting SRO correction model

#### `--coeff_out`, `-cout` (Default: sro_coeffs.out)
filename to store coefficients of SRO correction model

#### `--coeff_in`, `-cin` (Default: sro_coeffs.in)
filename containing initial coefficient for fitting SRO correction model

#### `--fit_correction_only`, `-fsr`
Flag to skip SRO correction optimisation but fit pre-existing data to SRO
correction model

#### `--norm_constraint`, `-nc`
Flag to enable norm constraint

#### `--constr_tol`, `-ctol` (Default: 1e-08)
Value for maximum allowed constraints violation

#### `--maxiter`, `-mit` (Default: 5000)
Maximum no. of iterations for local minimizer

#### `--maxiter_linprog`, `-mitlin` (Default: 500000)
Maximum no. of iterations for linear programming search

#### `--global_iterations`, `-git` (Default: 50)
No. of global optimisations search steps

#### `--xtol`, `-x` (Default: 1e-08)
Tolerance for termination of local minimizer by the change of correlations.
The algorithm will terminate when ``tr_radius < xtol``, where
``tr_radius`` is the radius of the trust region used in the algorithm

#### `--gtol`, `-g` (Default: 1e-08)
The algorithm will terminate when both the infinity norm (i.e., max abs value)
of the Lagrangian gradient and the constraint violation
are smaller than ``gtol``

#### `--barrier_tol`, `-btol` (Default: 1e-08)
Threshold on the barrier parameter for the algorithm termination.

#### `--initial_tr_radius`, `-itr` (Default: 1)
Initial trust radius. It reflects the trust the algorithm puts in the
local approximation of the optimization problem. For an accurate local
approximation                        the trust-region should be large and for
an approximation valid                        only close to the current point
it should be a small one.

#### `--initial_constr_penalty`, `-icp` (Default: 1)
Initial Constraint Penalty. The penalty parameter is used for
balancing the requirements of decreasing the objective function
and satisfying the constraints.

#### `--fit_ordered_only`, `-foo`
Flag to find find ordered state only and exit.

#### `--method_linprog`, `-linm` (Default: highs)
Method of linear programming for finding ordered correlations

#### `--basinhopping`, `-bh` (Default: False)
Flag to switch on basinhopping

#### `--verbose`, `-v` (Default: 0)
Indicate the verbosity of the fit

#### `--approx_deriv`, `-ad`
Flag to enable estimation of derivatives

#### `--earlystop`, `-es` (Default: 20)
Number of steps to break out of trials if no new minima has been found

#### `--initial_stepsize`, `-is` (Default: 0.1)
Initial stepsize of the basinhopping algorithm from the disordered phase
