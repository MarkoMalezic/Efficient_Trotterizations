#include "Exact.h"

// Exact constructor
template <typename RealT>
Exact<RealT>::Exact(const Model<RealT> &model_)
    : model(model_), dim(pow(2, model_.L))
{
  H = Operators::MatrixCX<RealT>::Zero(dim, dim);
  buildH();
  diagonalize();
}

// Exact deconstructor
template <typename RealT>
Exact<RealT>::~Exact()
{
}

// Method to build the Hamiltonian
template <typename RealT>
void Exact<RealT>::buildH()
{
  // Add the X terms
  for (int i{0}; i < model.L; ++i)
  {
    for (int j{0}; j < model.termsX[i].ops.size(); ++j)
    {
      H += Operators::embed_op2(model.termsX[i].ops[j], model.termsX[i].sites[0], model.L) * model.termsX[i].coefs[j];
    }
  }
  // Add the Y terms
  for (int i{0}; i < model.L; ++i)
  {
    for (int j{0}; j < model.termsY[i].ops.size(); ++j)
    {
      H += Operators::embed_op2(model.termsY[i].ops[j], model.termsY[i].sites[0], model.L) * model.termsY[i].coefs[j];
    }
  }
  // Add the Z terms
  for (int i{0}; i < model.L; ++i)
  {
    for (int j{0}; j < model.termsZ[i].ops.size(); ++j)
    {
      H += Operators::embed_op2(model.termsZ[i].ops[j], model.termsZ[i].sites[0], model.L) * model.termsZ[i].coefs[j];
    }
  }
}

// Method to diagonalize the Hamiltonian
template <typename RealT>
void Exact<RealT>::diagonalize()
{
  Eigen::SelfAdjointEigenSolver<Operators::MatrixCX<RealT>> solver(H);
  eigvals = solver.eigenvalues();
  eigvecs = solver.eigenvectors();
}

// Method to construct the real time evolution operator
template <typename RealT>
Operators::MatrixCX<RealT> Exact<RealT>::real_evolve_op(const RealT &t)
{
  Operators::VectorCX<RealT> exp_eigvals = eigvals.unaryExpr([t](complex<RealT> val)
                                                             { return exp(complex<RealT>(0.0, -1.0) * val * t); });
  return eigvecs * exp_eigvals.asDiagonal() * eigvecs.adjoint();
}

// Method to construct the imaginary time evolution operator
template <typename RealT>
Operators::MatrixCX<RealT> Exact<RealT>::imag_evolve_op(const RealT &t)
{
  Operators::VectorCX<RealT> exp_eigvals = eigvals.unaryExpr([t](complex<RealT> val)
                                                             { return exp(- val * t); });
  return eigvecs * exp_eigvals.asDiagonal() * eigvecs.adjoint();
}

// Method to real time evolve a state
template <typename RealT>
Operators::VectorCX<RealT> Exact<RealT>::real_evolve(const RealT &t, const Operators::VectorCX<RealT> &state)
{
  Operators::VectorCX<RealT> result = eigvecs.adjoint() * state;
  result = eigvals.unaryExpr([t](complex<RealT> val)
                             { return exp(complex<RealT>(0.0, -1.0) * val * t); }).array() * result.array();
  return eigvecs * result;
}

// Method to real time evolve a state
template <typename RealT>
Operators::VectorCX<RealT> Exact<RealT>::imag_evolve(const RealT &t, const Operators::VectorCX<RealT> &state)
{
  Operators::VectorCX<RealT> result = eigvecs.adjoint() * state;
  result = eigvals.unaryExpr([t](complex<RealT> val)
                             { return exp(- val * t); }).array() * result.array();
  result = eigvecs * result;
  return result / result.norm();
}

template class Exact<double>;
template class Exact<long double>;
