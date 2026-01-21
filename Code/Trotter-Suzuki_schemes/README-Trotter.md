# Trotter-Suzuki decompostion schemes implementation

The code implements the recursive formulae for high order Trotterizations, which allows the evaluation of their efficiency, if the scheme parameters are known.
If the parameters aren't known, the implementation includes a minimizer, that optimizes the theoretical error and provides efficient Trotter-Suzuki schemes.
In its current state the program supports schemes up to order 6.
The code is written in C++, with a dependency on the Eigen linear algebra library.
The building is done by CMake (see `CMakeLists.txt`) and the instructions for this can be found inside `build-run_instructions.txt`.
The class structure can be seen in the schematic below.

## Routines

The program is structured in a way that allows for a symbolic evaluation of schemes as well as a numeric one.
The trade-off between the two is computational cost and accuracy.
However, the loss in accuracy for the numeric part is not very large, and thus we recommend the use of this implementation.
The symbolic part is useful for final accurate results and obtaining the exact polynomial structure of scheme coefficients.
Unfortunately, the memory resources in this implementation scale exponentially, and its usefulness at higher orders is unlikely.

The program also includes an Input/Output component, which makes its usage simple.
One simply needs to give the program an input file with the correct parameters (`/path_to_build/main input_file.in`), and plenty of examples can be found in [Input_files](/Input_files).
More details are described below, but the general parameters needed are:
- routine: 0 - Scheme minimization, 1 - Scheme evaluation 
- mode: 0 - Numeric, 1 - Symbolic
- scalar_type: double, long_double, complex_double and complex_long_double. From this scalar_precision = double or long_double

### Numerical Scheme (routine: 1, mode: 0)

Numerically evaluates a scheme for some evaluation parameters a_eval and b_eval, which are multiplied by 1/2 in relation to those found in the paper (Important distinction!).
Returns the theoretical error and the efficiency per order.

Parameters:
- save_dir (string): Path to save directory
- order (2, 4, 6): Order of the scheme
- no_cycles (1, 2, 3, ...): Number of cycles
- a_eval (vector<scalar_type>): A vector of scheme evaluation parameters a
- b_eval (vector<scalar_type>): A vector of scheme evaluation parameters b
- verbose (0, 1, 2): Option to print values of scheme coefficients between iterations

### Numerical Minimization (routine: 0, mode: 0)

Numerically minimizes a scheme according to a subroutine.
The shared parameters between schemes are:
- save_dir (string): Path to save directory
- order (2, 4, 6): Order of the scheme
- no_cycles (1, 2, 3, ...): Number of cycles
- method (string): Subroutine name

#### Subroutine: minimize

Minimizes once with given initial parameters and returns the evaluated theoretical error and efficiency per order.

Hyperparameters:
- a_init (vector<scalar_type>): A vector of initial scheme parameters a
- b_init (vector<scalar_type>): A vector of initial scheme parameters b
- eps1 (array<scalar_precision, 4>): an array of convergence criteria and acceptance condition
- wi (vector<scalar_precision>): a vector of weights per order
- step (scalar_precision): derivative step size
- n_iter (int): Number of iterations before stopping the minimization
- Ls (array<scalar_precision, 2>): lambda iteration update constants
- lambda (scalar_precision): initial lambda value

#### Subroutine: minimize_origin

Minimizes (with an additional distance from the origin term) once with given initial parameters and returns the evaluated theoretical error and efficiency per order.

Additional hyperparameter:
- ratio (scalar_precision): Ratio between the leading order error and the distance from the origin term

#### Subroutine: min_twostep

Minimizes with two steps with given initial parameters and returns the evaluated theoretical error and efficiency per order.
First step minimizes the full chi2 function, while the second one imposes the constraints.

Hyperparameters:
- a_init (vector<scalar_type>): A vector of initial scheme parameters a
- b_init (vector<scalar_type>): A vector of initial scheme parameters b
- eps1 (array<scalar_type, 4>): an array of convergence criteria and acceptance condition for the first step
- eps2 (array<scalar_type, 4>): an array of convergence criteria and acceptance condition for the second step
- wi (vector<scalar_type>): a vector of weights per order
- step (scalar_type): derivative step size
- n_iter (int): Number of iterations before stopping the minimization
- Ls (array<scalar_type, 2>): lambda iteration update constants
- lambda (scalar_type): initial lambda value
- freeze (0, 1): Boolean to either freeze the Hessian or not(1)

#### Subroutine: min_twostep_origin

Minimizes with two steps with given initial parameters and returns the evaluated theoretical error and efficiency per order.
First step minimizes the full chi2 function (with an additional distance from the origin term), while the second one imposes the constraints.

Additional hyperparameter:
- ratio (scalar_precision): Ratio between the leading order error and the distance from the origin term

#### Subroutine: find

Runs the minimizer many times with random initial conditions to find as many distinct minima.
Returns a table of found minima along with their theoretical efficiencies.

Additional hyperparameters:
- N (int): Number of initial conditions to minimize
- steps (1, 2): How many steps to take for the minimization - runs either minimize(1) or min_twostep(2)
- threshold (scalar_precision): Threshold to cut unconverged schemes
- mu (scalar_precision): Mean of the normal distribution, which samples the initial conditions
- sigma (scalar_precision): The standard deviation, which samples the initial conditions
- tol (scalar_precision): Comparison tolerance between similar converged schemes
- freeze (0, 1): Boolean to either freeze (1) the Hessian or not(0)

#### Subroutine: find_origin

Runs the minimizer (with the additional distance from the origin term) many times with random initial conditions to find as many distinct minima.
Returns a table of found minima along with their theoretical efficiencies.

Additional hyperparameters:
- ratio (scalar_precision): Ratio between the leading order error and the distance from the origin term


## Dependency graph

The diagram below provides an overview of the class structure.
It illustrates the dependencies between the symbolic and numeric components, as well as the shared helper modules that connect them.

![Class structure](Class_structure.png)

