#include "Trotter.h"

// Trotter constructor
template <typename Vec>
Trotter<Vec>::Trotter(const Model<RealT> &model_, const int &q_, const bool &group_, const Vec &a_vec_, const Vec &b_vec_)
    : model(model_), dim(pow(2, model_.L)), q(q_), group(group_), a_vec(a_vec_), b_vec(b_vec_)
{
  c_vec.resize(q);
  d_vec.resize(q);
  to_ramp();
  cout << "c_vec: ";
  cout << fixed << setprecision(20) << c_vec << endl;
  cout << "d_vec: ";
  cout << fixed << setprecision(20) << d_vec << endl;
}

// Trotter deconstructor
template <typename Vec>
Trotter<Vec>::~Trotter()
{
}

// Method to transfrom from the symmetric basis to the standard basis
template <typename Vec>
void Trotter<Vec>::to_standard()
{
  Vec a_new(a_vec.size());
  Vec b_new(b_vec.size());
  if (q % 2 == 0)
  {
    for (size_t i{0}; i < a_vec.size(); ++i)
    {
      if (i == q / 2)
      {
        a_new[i] = a_vec[0];
      }
      else
      {
        a_new[i] = a_vec[a_vec.size() - 1 - i] / RealT(2.0);
      }
    }
    for (size_t i{0}; i < b_vec.size(); ++i)
    {
      b_new[i] = b_vec[b_vec.size() - 1 - i] / RealT(2.0);
    }
  }
  else
  {
    for (size_t i{0}; i < a_vec.size(); ++i)
    {
      a_new[i] = a_vec[a_vec.size() - 1 - i] / RealT(2.0);
    }
    for (size_t i{0}; i < b_vec.size(); ++i)
    {
      if (i == (q - 1) / 2)
      {
        b_new[i] = b_vec[0];
      }
      else
      {
        b_new[i] = b_vec[b_vec.size() - 1 - i] / RealT(2.0);
      }
    }
  }
  a_vec = a_new;
  b_vec = b_new;
}

// Method to transfrom from the symmetric basis to the ramp basis
template <typename Vec>
void Trotter<Vec>::to_ramp()
{
  to_standard();
  if (q % 2 == 0)
  {
    int back = 0;
    for (int i{0}; i < q; ++i)
    {
      if (i == 0)
      {
        c_vec[0] = a_vec[0];
        d_vec[0] = b_vec[0] - c_vec[0];
      }
      else if (i < q / 2)
      {
        c_vec[i] = a_vec[i] - d_vec[i - 1];
        d_vec[i] = b_vec[i] - c_vec[i];
      }
      else
      {
        c_vec[i] = a_vec[q / 2 - back] - d_vec[i - 1];
        d_vec[i] = b_vec[(q - 2) / 2 - back] - c_vec[i];
        back += 1;
      }
    }
  }
  else
  {
    int back = 0;
    for (int i{0}; i < q; ++i)
    {
      if (i == 0)
      {
        c_vec[0] = a_vec[0];
        d_vec[0] = b_vec[0] - c_vec[0];
      }
      else if (i < (q + 1) / 2)
      {
        c_vec[i] = a_vec[i] - d_vec[i - 1];
        d_vec[i] = b_vec[i] - c_vec[i];
      }
      else
      {
        c_vec[i] = a_vec[(q - 1) / 2 - back] - d_vec[i - 1];
        d_vec[i] = b_vec[(q - 3) / 2 - back] - c_vec[i];
        back += 1;
      }
    }
  }
}

// Method to compute an uphill ramp in the building of the step operator (real prefactor)
template <typename Vec>
void Trotter<Vec>::step_ramp_up(Operators::MatrixCX<RealT> &step_op, const RealT &prefac, const bool &imag_time)
{
  // Local grouping
  if (group)
  {
    // First: Loop over the terms in the Hamiltonian (chain)
    // Second: Loop over the directions
    for (int j{0}; j < model.L; ++j)
    {
      // Loop over the operators in the X direction
      if (model.termsX.size() > 0)
      {
        for (int k{0}; k < model.termsX[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opX{Operators::exp_embed_op2(prefac * model.termsX[j].coefs[k], model.termsX[j].ops[k], model.termsX[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opX;
        }
      }
      // Loop over the operators in the Y direction
      if (model.termsY.size() > 0)
      {
        for (int k{0}; k < model.termsY[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opY{Operators::exp_embed_op2(prefac * model.termsY[j].coefs[k], model.termsY[j].ops[k], model.termsY[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opY;
        }
      }
      // Loop over the operators in the Z direction
      if (model.termsZ.size() > 0)
      {
        for (int k{0}; k < model.termsZ[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opZ{Operators::exp_embed_op2(prefac * model.termsZ[j].coefs[k], model.termsZ[j].ops[k], model.termsZ[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opZ;
        }
      }
    }
  }
  // Global grouping
  else
  {
    // First: Loop over the directions
    // Second: Loop over the terms in the Hamiltonian (chain)
    // Loop over the operators in the X direction
    if (model.termsX.size() > 0)
    {
      for (int j{0}; j < model.L; ++j)
      {
        for (int k{0}; k < model.termsX[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opX{Operators::exp_embed_op2(prefac * model.termsX[j].coefs[k], model.termsX[j].ops[k], model.termsX[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opX;
        }
      }
    }
    // Loop over the operators in the Y direction
    if (model.termsY.size() > 0)
    {
      for (int j{0}; j < model.L; ++j)
      {
        for (int k{0}; k < model.termsY[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opY{Operators::exp_embed_op2(prefac * model.termsY[j].coefs[k], model.termsY[j].ops[k], model.termsY[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opY;
        }
      }
    }
    // Loop over the operators in the Z direction
    if (model.termsZ.size() > 0)
    {
      for (int j{0}; j < model.L; ++j)
      {
        for (int k{0}; k < model.termsZ[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opZ{Operators::exp_embed_op2(prefac * model.termsZ[j].coefs[k], model.termsZ[j].ops[k], model.termsZ[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opZ;
        }
      }
    }
  }
}

// Method to compute an uphill ramp in the building of the step operator (complex prefactor)
template <typename Vec>
void Trotter<Vec>::step_ramp_up(Operators::MatrixCX<RealT> &step_op, const complex<RealT> &prefac, const bool &imag_time)
{
  // Local grouping
  if (group)
  {
    // First: Loop over the terms in the Hamiltonian (chain)
    // Second: Loop over the directions
    for (int j{0}; j < model.L; ++j)
    {
      // Loop over the operators in the X direction
      if (model.termsX.size() > 0)
      {
        for (int k{0}; k < model.termsX[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opX{Operators::exp_embed_op2(prefac * model.termsX[j].coefs[k], model.termsX[j].ops[k], model.termsX[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opX;
        }
      }
      // Loop over the operators in the Y direction
      if (model.termsY.size() > 0)
      {
        for (int k{0}; k < model.termsY[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opY{Operators::exp_embed_op2(prefac * model.termsY[j].coefs[k], model.termsY[j].ops[k], model.termsY[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opY;
        }
      }
      // Loop over the operators in the Z direction
      if (model.termsZ.size() > 0)
      {
        for (int k{0}; k < model.termsZ[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opZ{Operators::exp_embed_op2(prefac * model.termsZ[j].coefs[k], model.termsZ[j].ops[k], model.termsZ[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opZ;
        }
      }
    }
  }
  // Global grouping
  else
  {
    // First: Loop over the directions
    // Second: Loop over the terms in the Hamiltonian (chain)
    // Loop over the operators in the X direction
    if (model.termsX.size() > 0)
    {
      for (int j{0}; j < model.L; ++j)
      {
        for (int k{0}; k < model.termsX[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opX{Operators::exp_embed_op2(prefac * model.termsX[j].coefs[k], model.termsX[j].ops[k], model.termsX[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opX;
        }
      }
    }
    // Loop over the operators in the Y direction
    if (model.termsY.size() > 0)
    {
      for (int j{0}; j < model.L; ++j)
      {
        for (int k{0}; k < model.termsY[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opY{Operators::exp_embed_op2(prefac * model.termsY[j].coefs[k], model.termsY[j].ops[k], model.termsY[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opY;
        }
      }
    }
    // Loop over the operators in the Z direction
    if (model.termsZ.size() > 0)
    {
      for (int j{0}; j < model.L; ++j)
      {
        for (int k{0}; k < model.termsZ[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opZ{Operators::exp_embed_op2(prefac * model.termsZ[j].coefs[k], model.termsZ[j].ops[k], model.termsZ[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opZ;
        }
      }
    }
  }
}

// Method to compute a downhill ramp in the building of the step operator (real prefactor)
template <typename Vec>
void Trotter<Vec>::step_ramp_down(Operators::MatrixCX<RealT> &step_op, const RealT &prefac, const bool &imag_time)
{
  // Local grouping
  if (group)
  {
    // First: Loop over the terms in the Hamiltonian (chain)
    // Second: Loop over the directions
    for (int j{model.L - 1}; j >= 0; --j)
    {
      // Loop over the operators in the Z direction
      if (model.termsZ.size() > 0)
      {
        for (int k{0}; k < model.termsZ[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opZ{Operators::exp_embed_op2(prefac * model.termsZ[j].coefs[k], model.termsZ[j].ops[k], model.termsZ[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opZ;
        }
      }
      // Loop over the operators in the Y direction
      if (model.termsY.size() > 0)
      {
        for (int k{0}; k < model.termsY[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opY{Operators::exp_embed_op2(prefac * model.termsY[j].coefs[k], model.termsY[j].ops[k], model.termsY[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opY;
        }
      }
      // Loop over the operators in the X direction
      if (model.termsX.size() > 0)
      {
        for (int k{0}; k < model.termsX[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opX{Operators::exp_embed_op2(prefac * model.termsX[j].coefs[k], model.termsX[j].ops[k], model.termsX[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opX;
        }
      }
    }
  }
  // Global grouping
  else
  {
    // First: Loop over the directions
    // Second: Loop over the terms in the Hamiltonian (chain)
    // Loop over the operators in the Z direction
    if (model.termsZ.size() > 0)
    {
      for (int j{model.L - 1}; j >= 0; --j)
      {
        for (int k{0}; k < model.termsZ[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opZ{Operators::exp_embed_op2(prefac * model.termsZ[j].coefs[k], model.termsZ[j].ops[k], model.termsZ[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opZ;
        }
      }
    }
    // Loop over the operators in the Y direction
    if (model.termsY.size() > 0)
    {
      for (int j{model.L - 1}; j >= 0; --j)
      {
        for (int k{0}; k < model.termsY[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opY{Operators::exp_embed_op2(prefac * model.termsY[j].coefs[k], model.termsY[j].ops[k], model.termsY[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opY;
        }
      }
    }
    // Loop over the operators in the X direction
    if (model.termsX.size() > 0)
    {
      for (int j{model.L - 1}; j >= 0; --j)
      {
        for (int k{0}; k < model.termsX[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opX{Operators::exp_embed_op2(prefac * model.termsX[j].coefs[k], model.termsX[j].ops[k], model.termsX[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opX;
        }
      }
    }
  }
}

// Method to compute a downhill ramp in the building of the step operator (complex prefactor)
template <typename Vec>
void Trotter<Vec>::step_ramp_down(Operators::MatrixCX<RealT> &step_op, const complex<RealT> &prefac, const bool &imag_time)
{
  // Local grouping
  if (group)
  {
    // First: Loop over the terms in the Hamiltonian (chain)
    // Second: Loop over the directions
    for (int j{model.L - 1}; j >= 0; --j)
    {
      // Loop over the operators in the Z direction
      if (model.termsZ.size() > 0)
      {
        for (int k{0}; k < model.termsZ[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opZ{Operators::exp_embed_op2(prefac * model.termsZ[j].coefs[k], model.termsZ[j].ops[k], model.termsZ[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opZ;
        }
      }
      // Loop over the operators in the Y direction
      if (model.termsY.size() > 0)
      {
        for (int k{0}; k < model.termsY[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opY{Operators::exp_embed_op2(prefac * model.termsY[j].coefs[k], model.termsY[j].ops[k], model.termsY[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opY;
        }
      }
      // Loop over the operators in the X direction
      if (model.termsX.size() > 0)
      {
        for (int k{0}; k < model.termsX[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opX{Operators::exp_embed_op2(prefac * model.termsX[j].coefs[k], model.termsX[j].ops[k], model.termsX[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opX;
        }
      }
    }
  }
  // Global grouping
  else
  {
    // First: Loop over the directions
    // Second: Loop over the terms in the Hamiltonian (chain)
    // Loop over the operators in the Z direction
    if (model.termsZ.size() > 0)
    {
      for (int j{model.L - 1}; j >= 0; --j)
      {
        for (int k{0}; k < model.termsZ[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opZ{Operators::exp_embed_op2(prefac * model.termsZ[j].coefs[k], model.termsZ[j].ops[k], model.termsZ[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opZ;
        }
      }
    }
    // Loop over the operators in the Y direction
    if (model.termsY.size() > 0)
    {
      for (int j{model.L - 1}; j >= 0; --j)
      {
        for (int k{0}; k < model.termsY[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opY{Operators::exp_embed_op2(prefac * model.termsY[j].coefs[k], model.termsY[j].ops[k], model.termsY[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opY;
        }
      }
    }
    // Loop over the operators in the X direction
    if (model.termsX.size() > 0)
    {
      for (int j{model.L - 1}; j >= 0; --j)
      {
        for (int k{0}; k < model.termsX[j].ops.size(); ++k)
        {
          // Exactly exponentiate and embed the operators in the full Hilbert space
          Operators::MatrixCX<RealT> opX{Operators::exp_embed_op2(prefac * model.termsX[j].coefs[k], model.termsX[j].ops[k], model.termsX[j].sites[0], model.L, imag_time)};
          // Multiply the step operator S with the exponential exp (-i h c/d op)
          step_op *= opX;
        }
      }
    }
  }
}

// Method to build the step operator S(h)
template <typename Vec>
Operators::MatrixCX<typename Eigen::NumTraits<typename Vec::Scalar>::Real> Trotter<Vec>::build_step_op(const RealT &h, const bool &imag_time)
{
  Operators::MatrixCX<RealT> step_op = Operators::MatrixCX<RealT>::Identity(dim, dim);
  // Loop over ramps (q)
  for (int i{0}; i < q; ++i)
  {
    // Multiply the step operator with the uphill ramp from the right
    step_ramp_up(step_op, c_vec[i] * h, imag_time);
    // Multiply the step operator with the downhill ramp from the right
    step_ramp_down(step_op, d_vec[i] * h, imag_time);
  }
  return step_op;
}

// Method to build the time evolution operator U(t)
template <typename Vec>
Operators::MatrixCX<typename Eigen::NumTraits<typename Vec::Scalar>::Real> Trotter<Vec>::build_evolve_op(const RealT &t, const RealT &h, const bool &imag_time)
{
  Operators::MatrixCX<RealT> evolve_op = Operators::MatrixCX<RealT>::Identity(dim, dim);
  Operators::MatrixCX<RealT> step_op = build_step_op(h, imag_time);
  int steps = static_cast<int>(t / h);
  for (int i{0}; i < steps; ++i)
  {
    evolve_op *= step_op;
  }
  return evolve_op;
}

// Method to compute an uphill ramp on the state (real prefactor)
template <typename Vec>
void Trotter<Vec>::take_ramp_up(Operators::VectorCX<RealT> &state, const RealT &prefac, const bool &imag_time)
{
  // Local grouping
  if (group)
  {
    // First: Loop over the terms in the Hamiltonian (chain)
    // Second: Loop over the directions
    for (int j{model.L - 1}; j >= 0; --j)
    {
      // Loop over the operators in the Z direction
      if (model.termsZ.size() > 0)
      {
        Operators::MatrixC4<RealT> opZ = model.termsZ[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opZ, model.termsZ[j].sites, model.L);
      }
      // Loop over the operators in the Y direction
      if (model.termsY.size() > 0)
      {
        Operators::MatrixC4<RealT> opY = model.termsY[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opY, model.termsY[j].sites, model.L);
      }
      // Loop over the operators in the X direction
      if (model.termsX.size() > 0)
      {
        Operators::MatrixC4<RealT> opX = model.termsX[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opX, model.termsX[j].sites, model.L);
      }
    }
  }
  // Global grouping
  else
  {
    // First: Loop over the directions
    // Second: Loop over the terms in the Hamiltonian (chain)
    // Loop over the operators in the Z direction
    if (model.termsZ.size() > 0)
    {
      for (int j{model.L - 1}; j >= 0; --j)
      {
        Operators::MatrixC4<RealT> opZ = model.termsZ[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opZ, model.termsZ[j].sites, model.L);
      }
    }
    // Loop over the operators in the Y direction
    if (model.termsY.size() > 0)
    {
      for (int j{model.L - 1}; j >= 0; --j)
      {
        Operators::MatrixC4<RealT> opY = model.termsY[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opY, model.termsY[j].sites, model.L);
      }
    }
    // Loop over the operators in the X direction
    if (model.termsX.size() > 0)
    {
      for (int j{model.L - 1}; j >= 0; --j)
      {
        Operators::MatrixC4<RealT> opX = model.termsX[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opX, model.termsX[j].sites, model.L);
      }
    }
  }
}

// Method to compute an uphill ramp on the state (complex prefactor)
template <typename Vec>
void Trotter<Vec>::take_ramp_up(Operators::VectorCX<RealT> &state, const complex<RealT> &prefac, const bool &imag_time)
{
  // Local grouping
  if (group)
  {
    // First: Loop over the terms in the Hamiltonian (chain)
    // Second: Loop over the directions
    for (int j{model.L - 1}; j >= 0; --j)
    {
      // Loop over the operators in the Z direction
      if (model.termsZ.size() > 0)
      {
        Operators::MatrixC4<RealT> opZ = model.termsZ[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opZ, model.termsZ[j].sites, model.L);
      }
      // Loop over the operators in the Y direction
      if (model.termsY.size() > 0)
      {
        Operators::MatrixC4<RealT> opY = model.termsY[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opY, model.termsY[j].sites, model.L);
      }
      // Loop over the operators in the X direction
      if (model.termsX.size() > 0)
      {
        Operators::MatrixC4<RealT> opX = model.termsX[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opX, model.termsX[j].sites, model.L);
      }
    }
  }
  // Global grouping
  else
  {
    // First: Loop over the directions
    // Second: Loop over the terms in the Hamiltonian (chain)
    // Loop over the operators in the Z direction
    if (model.termsZ.size() > 0)
    {
      for (int j{model.L - 1}; j >= 0; --j)
      {
        Operators::MatrixC4<RealT> opZ = model.termsZ[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opZ, model.termsZ[j].sites, model.L);
      }
    }
    // Loop over the operators in the Y direction
    if (model.termsY.size() > 0)
    {
      for (int j{model.L - 1}; j >= 0; --j)
      {
        Operators::MatrixC4<RealT> opY = model.termsY[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opY, model.termsY[j].sites, model.L);
      }
    }
    // Loop over the operators in the X direction
    if (model.termsX.size() > 0)
    {
      for (int j{model.L - 1}; j >= 0; --j)
      {
        Operators::MatrixC4<RealT> opX = model.termsX[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opX, model.termsX[j].sites, model.L);
      }
    }
  }
}

// Method to compute a downhill ramp on the state (real prefactor)
template <typename Vec>
void Trotter<Vec>::take_ramp_down(Operators::VectorCX<RealT> &state, const RealT &prefac, const bool &imag_time)
{
  // Local grouping
  if (group)
  {
    // First: Loop over the terms in the Hamiltonian (chain)
    // Second: Loop over the directions
    for (int j{0}; j < model.L; ++j)
    {
      // Loop over the operators in the X direction
      if (model.termsX.size() > 0)
      {
        Operators::MatrixC4<RealT> opX = model.termsX[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opX, model.termsX[j].sites, model.L);
      }
      // Loop over the operators in the Y direction
      if (model.termsY.size() > 0)
      {
        Operators::MatrixC4<RealT> opY = model.termsY[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opY, model.termsY[j].sites, model.L);
      }
      // Loop over the operators in the Z direction
      if (model.termsZ.size() > 0)
      {
        Operators::MatrixC4<RealT> opZ = model.termsZ[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opZ, model.termsZ[j].sites, model.L);
      }
    }
  }
  // Global grouping
  else
  {
    // First: Loop over the directions
    // Second: Loop over the terms in the Hamiltonian (chain)
    // Loop over the operators in the X direction
    if (model.termsX.size() > 0)
    {
      for (int j{0}; j < model.L; ++j)
      {
        Operators::MatrixC4<RealT> opX = model.termsX[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opX, model.termsX[j].sites, model.L);
      }
    }
    // Loop over the operators in the Y direction
    if (model.termsY.size() > 0)
    {
      for (int j{0}; j < model.L; ++j)
      {
        Operators::MatrixC4<RealT> opY = model.termsY[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opY, model.termsY[j].sites, model.L);
      }
    }
    // Loop over the operators in the Z direction
    if (model.termsZ.size() > 0)
    {
      for (int j{0}; j < model.L; ++j)
      {
        Operators::MatrixC4<RealT> opZ = model.termsZ[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opZ, model.termsZ[j].sites, model.L);
      }
    }
  }
}

// Method to compute a downhill ramp on the state (complex prefactor)
template <typename Vec>
void Trotter<Vec>::take_ramp_down(Operators::VectorCX<RealT> &state, const complex<RealT> &prefac, const bool &imag_time)
{
  // Local grouping
  if (group)
  {
    // First: Loop over the terms in the Hamiltonian (chain)
    // Second: Loop over the directions
    for (int j{0}; j < model.L; ++j)
    {
      // Loop over the operators in the X direction
      if (model.termsX.size() > 0)
      {
        Operators::MatrixC4<RealT> opX = model.termsX[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opX, model.termsX[j].sites, model.L);
      }
      // Loop over the operators in the Y direction
      if (model.termsY.size() > 0)
      {
        Operators::MatrixC4<RealT> opY = model.termsY[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opY, model.termsY[j].sites, model.L);
      }
      // Loop over the operators in the Z direction
      if (model.termsZ.size() > 0)
      {
        Operators::MatrixC4<RealT> opZ = model.termsZ[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opZ, model.termsZ[j].sites, model.L);
      }
    }
  }
  // Global grouping
  else
  {
    // First: Loop over the directions
    // Second: Loop over the terms in the Hamiltonian (chain)
    // Loop over the operators in the X direction
    if (model.termsX.size() > 0)
    {
      for (int j{0}; j < model.L; ++j)
      {
        Operators::MatrixC4<RealT> opX = model.termsX[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opX, model.termsX[j].sites, model.L);
      }
    }
    // Loop over the operators in the Y direction
    if (model.termsY.size() > 0)
    {
      for (int j{0}; j < model.L; ++j)
      {
        Operators::MatrixC4<RealT> opY = model.termsY[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opY, model.termsY[j].sites, model.L);
      }
    }
    // Loop over the operators in the Z direction
    if (model.termsZ.size() > 0)
    {
      for (int j{0}; j < model.L; ++j)
      {
        Operators::MatrixC4<RealT> opZ = model.termsZ[j].prod_exp(prefac, model.L, imag_time);
        state = Operators::apply_op2(state, opZ, model.termsZ[j].sites, model.L);
      }
    }
  }
}

// Method to take a single step in the time evolution of a state
template <typename Vec>
Operators::VectorCX<typename Eigen::NumTraits<typename Vec::Scalar>::Real> Trotter<Vec>::take_step(const RealT &h, Operators::VectorCX<RealT> state, const bool &imag_time)
{
  // Loop over the ramps (q)
  for (int i{0}; i < q; ++i)
  {
    // Apply the downhill ramp to the state
    take_ramp_down(state, d_vec[q-1-i]*h, imag_time);
    // Apply the uphill ramp to the state
    take_ramp_up(state, c_vec[q-1-i]*h, imag_time);
  }
  return state / state.norm();
  //return state;
}

// Method to evolve a state in time
template <typename Vec>
Operators::VectorCX<typename Eigen::NumTraits<typename Vec::Scalar>::Real> Trotter<Vec>::evolve(const RealT &t, const int &Nt, Operators::VectorCX<RealT> state, const bool &imag_time)
{
  RealT h = static_cast<RealT>(t / Nt);
  for (int i{0}; i < Nt; ++i)
  {
    //cout << "Taking step " << i+1 << " / " << Nt << " with step size h = " << h << endl;
    state = take_step(h, state, imag_time);
  }
  return state;
}

template class Trotter<VectorXd>;
template class Trotter<VectorXcd>;
template class Trotter<VectorXld>;
template class Trotter<VectorXcld>;
