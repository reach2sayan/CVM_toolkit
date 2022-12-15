## CVM Toolkit
by Sayan Samanta and Axel van de Walle @ Brown University, RI USA
### Usage
```
usage: sro_correction [-h] [--seed SEED] [--disp] [--log_level {INFO,DEBUG,WARNING,ERROR}] [--log LOG] [--eci ECI] [--vmat VMAT]
                      [--clusters CLUSTERS] [--maximal_clusters MAXIMAL_CLUSTERS] [--clustermult CLUSTERMULT] [--kb KB] [--configmult CONFIGMULT]
                      [--config CONFIG] [--lat LAT] [--Tmin TMIN] [--Tmax TMAX] [--Tstep TSTEP] [--inJoules] [--sro_method SRO_METHOD]
                      [--coeff_out COEFF_OUT] [--coeff_in COEFF_IN] [--skip_sro_fit] [--fit_sro_only] [--norm_constraint] [--constr_tol CONSTR_TOL]
                      [--maxiter MAXITER] [--maxiter_linprog MAXITER_LINPROG] [--global_iterations GLOBAL_ITERATIONS] [--xtol XTOL] [--gtol GTOL]
                      [--barrier_tol BARRIER_TOL] [--initial_tr_radius INITIAL_TR_RADIUS] [--initial_constr_penalty INITIAL_CONSTR_PENALTY]
                      [--fit_ordered_only] [--method_linprog METHOD_LINPROG] [--basinhopping] [--verbose] [--approx_deriv] [--earlystop EARLYSTOP]
                      [--initial_stepsize INITIAL_STEPSIZE] [--out OUT]
```
### Parameters
#### Quick reference table
|Short|Long                      |Default                |Description                                                                                                                                                                                                                                                                                                                                              |
|-----|--------------------------|-----------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
|`-h` |`--help`                  |                       |show this help message and exit                                                                                                                                                                                                                                                                                                                          |
|     |`--seed`                  |`42`                   |set numpy random seed (for reproducible experiments)                                                                                                                                                                                                                                                                                                     |
|     |`--disp`                  |`False`                |Flag to turn intermediate outputs on                                                                                                                                                                                                                                                                                                                     |
|     |`--log_level`             |`INFO`                 |Flag to set output level                                                                                                                                                                                                                                                                                                                                 |
|     |`--log`                   |`log.out`              |Filename for the log file                                                                                                                                                                                                                                                                                                                                |
|     |`--eci`                   |`eci.out`              |file containing ECI's                                                                                                                                                                                                                                                                                                                                    |
|     |`--vmat`                  |`vmat.out`             |file containing the vmatrix                                                                                                                                                                                                                                                                                                                              |
|     |`--clusters`              |`clusters.out`         |file containing cluster description                                                                                                                                                                                                                                                                                                                      |
|     |`--maximal_clusters`      |`maxclus.in`           |file contain the maximal cluster coordinates                                                                                                                                                                                                                                                                                                             |
|     |`--clustermult`           |`clusmult.out`         |file containing cluster multiplicities                                                                                                                                                                                                                                                                                                                   |
|     |`--kb`                    |`configkb.out`         |file containing the kikuchi-barker coefficients                                                                                                                                                                                                                                                                                                          |
|     |`--configmult`            |`configmult.out`       |file containing cluster configuration multiplicities                                                                                                                                                                                                                                                                                                     |
|     |`--config`                |`config.out`           |file containing cluster configuration descriptions                                                                                                                                                                                                                                                                                                       |
|`-l` |`--lat`                   |`lat.in`               |contains the lattice description of the phase                                                                                                                                                                                                                                                                                                            |
|     |`--Tmin`                  |`100`                  |Minimum temperature for SRO correction                                                                                                                                                                                                                                                                                                                   |
|     |`--Tmax`                  |`2000`                 |Maximum temperature for SRO correction                                                                                                                                                                                                                                                                                                                   |
|     |`--Tstep`                 |`100`                  |Temperature increment for SRO correction                                                                                                                                                                                                                                                                                                                 |
|     |`--inJoules`              |                       |Sets energy units to Joules/mol                                                                                                                                                                                                                                                                                                                          |
|     |`--sro_method`            |`lm`                   |method of fitting SRO correction model                                                                                                                                                                                                                                                                                                                   |
|     |`--coeff_out`             |`sro_coeffs.out`       |filename to store coefficients of SRO correction model                                                                                                                                                                                                                                                                                                   |
|     |`--coeff_in`              |`sro_coeffs.in`        |filename containing initial coefficient for fitting SRO correction model                                                                                                                                                                                                                                                                                 |
|     |`--skip_sro_fit`          |                       |Flag to skip fitting SRO correction model                                                                                                                                                                                                                                                                                                                |
|     |`--fit_sro_only`          |                       |Flag to skip SRO correction optimisation but fit pre-existing data to SRO correction model                                                                                                                                                                                                                                                               |
|     |`--norm_constraint`       |                       |Flag to enable norm constraint                                                                                                                                                                                                                                                                                                                           |
|     |`--constr_tol`            |`1e-08`                |Value for maximum allowed constraints violation                                                                                                                                                                                                                                                                                                          |
|     |`--maxiter`               |`5000`                 |Maximum no. of iterations for local minimizer                                                                                                                                                                                                                                                                                                            |
|     |`--maxiter_linprog`       |`500000`               |Maximum no. of iterations for linear programming search                                                                                                                                                                                                                                                                                                  |
|     |`--global_iterations`     |`50`                   |No. of global optimisations search steps                                                                                                                                                                                                                                                                                                                 |
|     |`--xtol`                  |`1e-08`                |Tolerance for termination of local minimizer by the change of correlations.                        The algorithm will terminate when ``tr_radius < xtol`` , where    ``tr_radius`` is the radius of the trust region used in the algorithm                                                                                            |
|     |`--gtol`                  |`1e-08`                |The algorithm will terminate when both the infinity norm (i.e., max abs value) of the Lagrangian gradient and the constraint violation                        are smaller than ``gtol``                                                                                                |
|     |`--barrier_tol`           |`1e-08`                |Threshold on the barrier parameter for the algorithm termination.                                                                                                                                                                                                                                                                                        |
|     |`--initial_tr_radius`     |`1`                    |Initial trust radius. It reflects the trust the algorithm puts in the                        local approximation of the optimization problem. For an accurate local approximation the trust-region should be large and for an approximation valid only close to the current point it should be a small one.|
|     |`--initial_constr_penalty`|`1`                    |Initial Constraint Penalty. The penalty parameter is used for balancing the requirements of decreasing the objective function and satisfying the constraints.                                                                                                                                              |
|     |`--fit_ordered_only`      |                       |Flag to find find_ordered state only and exit.                                                                                                                                                                                                                                                                                                           |
|     |`--method_linprog`        |`highs`                |Method of linear programming for finding ordered correlations                                                                                                                                                                                                                                                                                            |
|     |`--basinhopping`          |`False`                |Flag to switch on basinhopping                                                                                                                                                                                                                                                                                                                           |
|`-v` |`--verbose`               |`0`                    |Indicate the verbosity of the fit                                                                                                                                                                                                                                                                                                                        |
|     |`--approx_deriv`          |                       |Flag to enable estimation of derivatives                                                                                                                                                                                                                                                                                                                 |
|     |`--earlystop`             |`20`                   |Number of steps to break out of trials if no new minima has been found                                                                                                                                                                                                                                                                                   |
|     |`--initial_stepsize`      |`0.1`                  |Initial stepsize of the basinhopping algorithm from the disordered phase                                                                                                                                                                                                                                                                                 |
|     |`--out`                   |`result-15Dec-1612.csv`|Indicates the name of the output file. The default name appends the current data and time                                                                                                                                                                                                                                         |

#### `-h`, `--help`
show this help message and exit

#### `--seed` (Default: 42)
set numpy random seed (for reproducible experiments)

#### `--disp` (Default: False)
Flag to turn intermediate outputs on

#### `--log_level` (Default: INFO)
Flag to set output level

#### `--log` (Default: log.out)
Filename for the log file

#### `--eci` (Default: eci.out)
file containing ECI's

#### `--vmat` (Default: vmat.out)
file containing the vmatrix

#### `--clusters` (Default: clusters.out)
file containing cluster description

#### `--maximal_clusters` (Default: maxclus.in)
file contain the maximal cluster coordinates

#### `--clustermult` (Default: clusmult.out)
file containing cluster multiplicities

#### `--kb` (Default: configkb.out)
file containing the kikuchi-barker coefficients

#### `--configmult` (Default: configmult.out)
file containing cluster configuration multiplicities

#### `--config` (Default: config.out)
file containing cluster configuration descriptions

#### `--lat`, `-l` (Default: lat.in)
contains the lattice description of the phase

#### `--Tmin` (Default: 100)
Minimum temperature for SRO correction

#### `--Tmax` (Default: 2000)
Maximum temperature for SRO correction

#### `--Tstep` (Default: 100)
Temperature increment for SRO correction

#### `--inJoules`
Sets energy units to Joules/mol

#### `--sro_method` (Default: lm)
method of fitting SRO correction model

#### `--coeff_out` (Default: sro_coeffs.out)
filename to store coefficients of SRO correction model

#### `--coeff_in` (Default: sro_coeffs.in)
filename containing initial coefficient for fitting SRO correction model

#### `--skip_sro_fit`
Flag to skip fitting SRO correction model

#### `--fit_sro_only`
Flag to skip SRO correction optimisation but fit pre-existing data to SRO
correction model

#### `--norm_constraint`
Flag to enable norm constraint

#### `--constr_tol` (Default: 1e-08)
Value for maximum allowed constraints violation

#### `--maxiter` (Default: 5000)
Maximum no. of iterations for local minimizer

#### `--maxiter_linprog` (Default: 500000)
Maximum no. of iterations for linear programming search

#### `--global_iterations` (Default: 50)
No. of global optimisations search steps

#### `--xtol` (Default: 1e-08)
Tolerance for termination of local minimizer by the change of correlations.
The algorithm will terminate when ``tr_radius < xtol``, where ``tr_radius`` is the radius of the trust region used in the algorithm

#### `--gtol` (Default: 1e-08)
The algorithm will terminate when both the infinity norm (i.e., max abs value)
of the Lagrangian gradient and the constraint violation are smaller than ``gtol``

#### `--barrier_tol` (Default: 1e-08)
Threshold on the barrier parameter for the algorithm termination.

#### `--initial_tr_radius` (Default: 1)
Initial trust radius. It reflects the trust the algorithm puts in the local  approximation of the optimization problem. For an accurate local approximation  the trust-region should be large and for an approximation valid                        only close to the current point it should be a small one.

#### `--initial_constr_penalty` (Default: 1)
Initial Constraint Penalty. The penalty parameter is used for balancing the requirements of decreasing the objective function and satisfying the constraints.

#### `--fit_ordered_only`
Flag to find find_ordered state only and exit.

#### `--method_linprog` (Default: highs)
Method of linear programming for finding ordered correlations

#### `--basinhopping` (Default: False)
Flag to switch on basinhopping

#### `--verbose`, `-v` (Default: 0)
Indicate the verbosity of the fit

#### `--approx_deriv`
Flag to enable estimation of derivatives

#### `--earlystop` (Default: 20)
Number of steps to break out of trials if no new minima has been found

#### `--initial_stepsize` (Default: 0.1)
Initial stepsize of the basinhopping algorithm from the disordered phase

#### `--out` (Default: result-%d%b%H%m.csv)
Indicates the name of the output file. The default name appends the current data and time
