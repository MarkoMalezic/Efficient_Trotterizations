#ifndef _TROTTER_H_
#define _TROTTER_H_

#include "Model.h"

using Eigen::VectorXd;
using Eigen::VectorXcd;
using VectorXld = Eigen::Matrix<long double, Eigen::Dynamic, 1>;
using VectorXcld = Eigen::Matrix<complex<long double>, Eigen::Dynamic, 1>;

// Class to store the Hamiltonian of the model (1 dimensional chain)
template <typename Vec>
class Trotter
{
public:
  using RealT = typename Eigen::NumTraits<typename Vec::Scalar>::Real;

  Model<RealT> model; // The model
  int dim;            // The dimension of the Hamiltonian
  int q;              // The number of cycles
  bool group;         // The grouping of stages (local = 0, global = 1)
  // Coefficients in the symmetric basis
  Vec a_vec; // The a vector
  Vec b_vec; // The b vector
  // Coefficients in the ramp basis
  Vec c_vec; // The c vector
  Vec d_vec; // The d vector

  // Trotter constructor
  Trotter(const Model<RealT> &model_, const int &q_, const bool &group_, const Vec &a_vec_, const Vec &b_vec_);

  // Trotter deconstructor
  ~Trotter();

  // Method to transfrom from the symmetric basis to the standard basis
  void to_standard();

  // Method to transfrom from the symmetric basis to the ramp basis
  void to_ramp();

  // Method to compute an uphill ramp in the building of the step operator (real prefactor)
  void step_ramp_up(Operators::MatrixCX<RealT> &step_op, const RealT &prefac, const bool &imag_time);

  // Method to compute an uphill ramp in the building of the step operator (complex prefactor)
  void step_ramp_up(Operators::MatrixCX<RealT> &step_op, const complex<RealT> &prefac, const bool &imag_time);

  // Method to compute a downhill ramp in the building of the step operator (real prefactor)
  void step_ramp_down(Operators::MatrixCX<RealT> &step_op, const RealT &prefac, const bool &imag_time);

  // Method to compute a downhill ramp in the building of the step operator (complex prefactor)
  void step_ramp_down(Operators::MatrixCX<RealT> &step_op, const complex<RealT> &prefac, const bool &imag_time);

  // Method to build the step operator S(h)
  Operators::MatrixCX<RealT> build_step_op(const RealT &h, const bool &imag_time = false);

  // Method to build the time evolution operator U(t)
  Operators::MatrixCX<RealT> build_evolve_op(const RealT &t, const RealT &h, const bool &imag_time = false);

  // Method to compute an uphill ramp on the state (real prefactor)
  void take_ramp_up(Operators::VectorCX<RealT> &state, const RealT &prefac, const bool &imag_time);

  // Method to compute an uphill ramp on the state (complex prefactor)
  void take_ramp_up(Operators::VectorCX<RealT> &state, const complex<RealT> &prefac, const bool &imag_time);

  // Method to compute a downhill ramp on the state (real prefactor)
  void take_ramp_down(Operators::VectorCX<RealT> &state, const RealT &prefac, const bool &imag_time);

  // Method to compute a downhill ramp on the state (complex prefactor)
  void take_ramp_down(Operators::VectorCX<RealT> &state, const complex<RealT> &prefac, const bool &imag_time);

  // Method to take a single step in the time evolution of a state
  Operators::VectorCX<RealT> take_step(const RealT &h, Operators::VectorCX<RealT> state, const bool &imag_time = false);

  // Method to evolve a state in time
  Operators::VectorCX<RealT> evolve(const RealT &t, const int &Nt, Operators::VectorCX<RealT> state, const bool &imag_time = false);
};

#endif // _TROTTER_H_
