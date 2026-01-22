# Time evolution of quantum symplectic systems

The code for the quantum dynamics of symplectic systems is implemented in `Symplectic_evolution.py`, which includes a class for Trotterized time evolution of a symplectic model, which can be included in the file as a new class.
If this class includes exact time evolution of it's model it is possible to compute the Trotter error.
For this reason, we implement the quantum Harmonic oscillator, which we use to test out the program.

Additionaly, the file script include an IO class, which handles an input file, which needs to be defined by the user, and outputs the simulated data.
An example for such an input file (`template.in`) can be found inside [Harmonic](Input_files/Errors/Harmonic/), with more scheme tests found in [tests](Input_files/Errors/Harmonic/tests).

