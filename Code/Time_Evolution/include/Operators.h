#include <complex>
#include <array>
#include <Eigen/Dense>
#include <unsupported/Eigen/KroneckerProduct>
#include <unsupported/Eigen/MatrixFunctions>

using namespace std;

namespace Operators
{
  // Define a template for precision type
  template <typename RealT>
  using MatrixC2 = Eigen::Matrix<complex<RealT>, 2, 2>;

  template <typename RealT>
  using MatrixC4 = Eigen::Matrix<complex<RealT>, 4, 4>;

  template <typename RealT>
  using MatrixCX = Eigen::Matrix<complex<RealT>, Eigen::Dynamic, Eigen::Dynamic>;

  template <typename RealT>
  using VectorCX = Eigen::Matrix<complex<RealT>, Eigen::Dynamic, 1>;

  // Identity matrix
  template <typename RealT>
  inline MatrixC2<RealT> Id()
  {
    return MatrixC2<RealT>::Identity();
  }

  // Pauli-X matrix
  template <typename RealT>
  inline MatrixC2<RealT> PauliX()
  {
    MatrixC2<RealT> X;
    X << 0, 1,
         1, 0;
    return X;
  }

  // Pauli-Y matrix
  template <typename RealT>
  inline MatrixC2<RealT> PauliY()
  {
    complex<RealT> I(0, 1);
    MatrixC2<RealT> Y;
    Y << 0, -I,
         I, 0;
    return Y;
  }

  // Pauli-Z matrix
  template <typename RealT>
  inline MatrixC2<RealT> PauliZ()
  {
    MatrixC2<RealT> Z;
    Z << 1, 0,
         0, -1;
    return Z;
  }

  // Kronecker product helper
  template <typename RealT>
  inline MatrixC4<RealT> kron(const MatrixC2<RealT> &A, const MatrixC2<RealT> &B)
  {
    return Eigen::kroneckerProduct(A, B).eval(); // .eval() ensures it's a MatrixC, not expression
  }

  // Embed 2 neighboring operators in the full Hilbert space
  template <typename RealT>
  inline MatrixCX<RealT> embed_op2(const array<MatrixC2<RealT>, 2> &ops, const int &site, const int &L)
  {
    if (L == 2)
    {
      if (site == 0)
      {
        return Eigen::kroneckerProduct(ops[0], ops[1]).eval();
      }
      else
      {
        return Eigen::kroneckerProduct(ops[1], ops[0]).eval();
      }
    }
    else
    {
      if (site == 0)
      {
        int dimr = pow(2, L - 2);
        return Eigen::kroneckerProduct(Eigen::kroneckerProduct(ops[0], ops[1]).eval(), MatrixCX<RealT>::Identity(dimr, dimr)).eval();
      }
      else if (site == L - 2)
      {
        int diml = pow(2, L - 2);
        return Eigen::kroneckerProduct(MatrixCX<RealT>::Identity(diml, diml), Eigen::kroneckerProduct(ops[0], ops[1]).eval()).eval();
      }
      else if (site == L - 1)
      {
        int dim = pow(2, L - 2);
        return Eigen::kroneckerProduct(Eigen::kroneckerProduct(ops[1], MatrixCX<RealT>::Identity(dim, dim)).eval(), ops[0]).eval();
      }
      else
      {
        int diml = pow(2, site);
        int dimr = pow(2, L - site - 2);
        return Eigen::kroneckerProduct(Eigen::kroneckerProduct(MatrixCX<RealT>::Identity(diml, diml), Eigen::kroneckerProduct(ops[0], ops[1]).eval()).eval(), MatrixCX<RealT>::Identity(dimr, dimr)).eval();
      }
    }
  }

  // Exactly exponentiate (Pauli matrices) with real coefficient c using real/imaginary (false/true) time
  // and embed 2 neighbouring operators in the full Hilbert space: exp (- i * c * op) or exp (- c * op)
  template <typename RealT>
  inline MatrixCX<RealT> exp_embed_op2(const RealT &c, const array<MatrixC2<RealT>, 2> &ops, const int &site, const int &L, const bool &imag_time = false)
  {
    MatrixCX<RealT> op_embed{embed_op2(ops, site, L)};
    if (imag_time)
    {
      return cosh(c) * MatrixCX<RealT>::Identity(pow(2, L), pow(2, L)) - sinh(c) * op_embed;

    }
    else
    {
      return cos(c) * MatrixCX<RealT>::Identity(pow(2, L), pow(2, L)) - complex<RealT>(0.0, 1.0) * sin(c) * op_embed;
    }
  }

  // Exactly exponentiate (Pauli matrices) with complex coefficient c using real/imaginary (false/true) time
  // and embed 2 neighbouring operators in the full Hilbert space: exp (- i * (a + ib) * op) or exp (-(a + ib) * op)
  template <typename RealT>
  inline MatrixCX<RealT> exp_embed_op2(const complex<RealT> &c, const array<MatrixC2<RealT>, 2> &ops, const int &site, const int &L, const bool &imag_time = false)
  {
    MatrixCX<RealT> op_embed{embed_op2(ops, site, L)};
    if (imag_time)
    {
      RealT cosb = cos(c.imag());
      RealT sinb = sin(c.imag());
      RealT cosha = cosh(c.real());
      RealT sinha = sinh(c.real());
      complex<RealT> i(0.0, 1.0);
      return (cosha*cosb + i*sinha*sinb) * MatrixCX<RealT>::Identity(pow(2, L), pow(2, L)) - (sinha*cosb + i*cosha*sinb) * op_embed;
    }
    else
    {
      RealT cosa = cos(c.real());
      RealT sina = sin(c.real());
      RealT coshb = cosh(c.imag());
      RealT sinhb = sinh(c.imag());
      complex<RealT> i(0.0, 1.0);
      return (coshb*cosa - i*sinhb*sina) * MatrixCX<RealT>::Identity(pow(2, L), pow(2, L)) + (sinhb*cosa - i*coshb*sina) * op_embed;
    }
  }

  // Exactly exponentiate (Pauli matrices) with real coefficient c
  // exp (- i * c * op) or exp (- c * op)
  template <typename RealT>
  inline MatrixC4<RealT> exp_op2(const RealT &c, const array<MatrixC2<RealT>, 2> &ops, const int &site, const int &L, const bool &imag_time = false)
  {
    MatrixC4<RealT> op;
    if (site == L - 1)
    {
      op = kron(ops[1], ops[0]);
    }
    else
    {
      op = kron(ops[0], ops[1]);
    }
    if (imag_time)
    {
      return cosh(c) * MatrixC4<RealT>::Identity() - sinh(c) * op;
    }
    else
    {
      return cos(c) * MatrixC4<RealT>::Identity() - complex<RealT>(0.0, 1.0) * sin(c) * op;
    }
  }

  // Exactly exponentiate (Pauli matrices) with complex coefficient c
  // exp exp (- i * (a + ib) * op) or exp (-(a + ib) * op)
  template <typename RealT>
  inline MatrixC4<RealT> exp_op2(const complex<RealT> &c, const array<MatrixC2<RealT>, 2> &ops, const int &site, const int &L, const bool &imag_time = false)
  {
    MatrixC4<RealT> op;
    if (site == L - 1)
    {
      op = kron(ops[1], ops[0]);
    }
    else
    {
      op = kron(ops[0], ops[1]);
    }
    if (imag_time)
    {
      RealT cosb = cos(c.imag());
      RealT sinb = sin(c.imag());
      RealT cosha = cosh(c.real());
      RealT sinha = sinh(c.real());
      complex<RealT> i(0.0, 1.0);
      return (cosha*cosb + i*sinha*sinb) * MatrixCX<RealT>::Identity(pow(2, L), pow(2, L)) - (sinha*cosb + i*cosha*sinb) * op;
    }
    else
    {
      RealT cosa = cos(c.real());
      RealT sina = sin(c.real());
      RealT coshb = cosh(c.imag());
      RealT sinhb = sinh(c.imag());
      complex<RealT> i(0.0, 1.0);
      return (coshb*cosa - i*sinhb*sina) * MatrixC4<RealT>::Identity() + (sinhb*cosa - i*coshb*sina) * op;
    }
  }

  // Apply 2 neighbouring operators to a state
  template <typename RealT>
  inline VectorCX<RealT> apply_op2(const VectorCX<RealT> &state, const MatrixC4<RealT> &op, array<int, 2> &sites, const int &L)
  {
    // Bit masks
    int mask1 = 1 << (L - 1 - sites[0]);
    int mask2 = 1 << (L - 1 - sites[1]);

    VectorCX<RealT> result = VectorCX<RealT>::Zero(state.size());

    for (int i{0}; i < state.size(); ++i)
    {
      // Extract bits at site1 and site2
      int bit1 = (i & mask1) >> (L - 1 - sites[0]);
      int bit2 = (i & mask2) >> (L - 1 - sites[1]);

      // Combine bits
      int b;
      if (sites[1] != 0)
      {
        b = (bit1 << 1) | bit2;
      }
      else
      {
        b = (bit2 << 1) | bit1;
      }

      // Loop over possible bp
      for (int bp = 0; bp < 4; ++bp)
      {
        int bit1p;
        int bit2p;
        if (sites[1] != 0)
        {
          bit1p = (bp >> 1) & 1;
          bit2p = bp & 1;
        }
        else
        {
          bit1p = bp & 1;
          bit2p = (bp >> 1) & 1;
        }

        int j = i;
        // Clear bits
        j &= ~mask1;
        j &= ~mask2;
        // Set bits to new values
        j |= (bit1p << (L - 1 - sites[0]));
        j |= (bit2p << (L - 1 - sites[1]));

        // Apply operator
        result[j] += op(bp, b) * state[i];
      }
    }
    return result;
  }
}
