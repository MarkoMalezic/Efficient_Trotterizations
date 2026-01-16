#ifndef _MODEL_H_
#define _MODEL_H_

#include <vector>
#include <iostream>
#include <iomanip>

#include "Operators.h"

using namespace std;

using Eigen::VectorXd;

template <typename RealT>
struct Hterm
{
  vector<array<typename Operators::MatrixC2<RealT>, 2>> ops; // The vector of operators
  vector<RealT> coefs;                                       // The vector of coefficients
  array<int, 2> sites;                                       // The operator sites

  // Hterm constructor
  Hterm(const vector<array<typename Operators::MatrixC2<RealT>, 2>> &ops_, vector<RealT> coefs_, const array<int, 2> &sites_)
      : ops(ops_), coefs(coefs_), sites(sites_) {}

  // Hterm default constructor
  Hterm() : ops(), coefs(), sites({0, 0}) {}

  // Hterm deconstructor
  ~Hterm() = default;

  // Method to compute a product of all exponentiated (commuting) operators using a real prefactor
  typename Operators::MatrixC4<RealT> prod_exp(const RealT &prefac, const int &L, const bool &imag_time)
  {
    typename Operators::MatrixC4<RealT> prod = Operators::MatrixC4<RealT>::Identity();
    for (int i{0}; i < ops.size(); ++i)
    {
      prod *= Operators::exp_op2(prefac * coefs[i], ops[i], sites[0], L, imag_time);
    }
    return prod;
  }

  // Method to compute a product of all exponentiated (commuting) operators using a complex prefactor
  typename Operators::MatrixC4<RealT> prod_exp(const complex<RealT> &prefac, const int &L, const bool &imag_time)
  {
    typename Operators::MatrixC4<RealT> prod = Operators::MatrixC4<RealT>::Identity();
    for (int i{0}; i < ops.size(); ++i)
    {
      prod *= Operators::exp_op2(prefac * coefs[i], ops[i], sites[0], L, imag_time);
    }
    return prod;
  }
};

// Class to store the Hamiltonian of the model (1 dimensional chain)
template <typename RealT>
class Model
{
public:
  int L;        // Length of the chain
  bool bound;   // Whether the chain is periodic (0) or open (1)
  string label; // The label of the model
  // Non-commutating terms of the Hamiltonian
  vector<Hterm<RealT>> termsX; // X terms
  vector<Hterm<RealT>> termsY; // Y terms
  vector<Hterm<RealT>> termsZ; // Z terms

  // Model constructor
  Model(const int &L_, const string &label_ = "", const bool &bound_ = false);

  // Model deconstructor
  virtual ~Model();

  // Method to add a single term to the Hamiltonian
  // op: The operator
  // dir: The direction (X = 0, Y = 1, Z = 2)
  // sites: The sites the operator acts on
  void add_term(const int &ind, const int &dir, const vector<array<typename Operators::MatrixC2<RealT>, 2>> &ops,
                vector<RealT> coefs, const array<int, 2> &sites);

  // Virtual method to build the Hamiltonian terms
  virtual void build_model();
};

template <typename RealT>
class ModelXZ : public Model<RealT>
{
public:
  RealT J;                                      // The coupling constant
  Eigen::Matrix<RealT, Eigen::Dynamic, 1> hmag; // The magnetic field

  ModelXZ(int L_, bool bound_ = false) : Model<RealT>(L_, "ModelXZ", bound_), J(1.0), hmag(L_) {}
  ModelXZ(int L_, RealT J_, Eigen::Matrix<RealT, Eigen::Dynamic, 1> hmag_, string label_ = "ModelXZ", bool bound_ = false) : Model<RealT>(L_, label_, bound_), J(J_), hmag(hmag_) {}

  void build_model() override;
};

template <typename RealT>
class ModelXXZ : public Model<RealT>
{
public:
  RealT J;                                      // The coupling constant
  Eigen::Matrix<RealT, Eigen::Dynamic, 1> hmag; // The magnetic field

  ModelXXZ(int L_, bool bound_ = false) : Model<RealT>(L_, "ModelXXZ", bound_), J(1.0), hmag(L_) {}
  ModelXXZ(int L_, RealT J_, Eigen::Matrix<RealT, Eigen::Dynamic, 1> hmag_, string label_ = "ModelXXZ", bool bound_ = false) : Model<RealT>(L_, label_, bound_), J(J_), hmag(hmag_) {}

  void build_model() override;
};

template <typename RealT>
class ModelHeisenberg : public Model<RealT>
{
public:
  array<RealT, 3> J;                            // The coupling constant
  Eigen::Matrix<RealT, Eigen::Dynamic, 1> hmag; // The magnetic field

  ModelHeisenberg(int L_, bool bound_ = false) : Model<RealT>(L_, "ModelHeisenberg", bound_), hmag(L_) {}
  ModelHeisenberg(int L_, array<RealT, 3> J_, Eigen::Matrix<RealT, Eigen::Dynamic, 1> hmag_, string label_ = "ModelHeisenberg", bool bound_ = false) : Model<RealT>(L_, label_, bound_), J(J_), hmag(hmag_) {}

  void build_model() override;
};

#endif // _MODEL_H_
