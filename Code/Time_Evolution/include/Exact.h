#ifndef _EXACT_H_
#define _EXACT_H_

#include "Model.h"
#include <Eigen/Eigenvalues>

// Class to store the Hamiltonian of the model (1 dimensional chain)
template <typename RealT>
class Exact
{
public:
  Model<RealT> model;                 // The model
  int dim;                            // The dimension of the Hamiltonian
  Operators::MatrixCX<RealT> H;       // The Hamiltonian
  Operators::VectorCX<RealT> eigvals; // The eigenvalues
  Operators::MatrixCX<RealT> eigvecs; // The eigenvectors

  // Exact constructor
  Exact(const Model<RealT> &model_);

  // Exact deconstructor
  ~Exact();

  // Method to build the Hamiltonian
  void buildH();

  // Method to diagonalize the Hamiltonian
  void diagonalize();

  // Method to construct the real time evolution operator
  Operators::MatrixCX<RealT> real_evolve_op(const RealT &t);

  // Method to construct the imaginary time evolution operator
  Operators::MatrixCX<RealT> imag_evolve_op(const RealT &t);

  // Method to real time evolve a state
  Operators::VectorCX<RealT> real_evolve(const RealT &t, const Operators::VectorCX<RealT> &state);

   // Method to imaginary time evolve a state
   Operators::VectorCX<RealT> imag_evolve(const RealT &t, const Operators::VectorCX<RealT> &state);
};

#endif // _EXACT_H_
