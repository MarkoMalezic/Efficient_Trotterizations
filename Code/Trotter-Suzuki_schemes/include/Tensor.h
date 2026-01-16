#ifndef _TENSOR_H_
#define _TENSOR_H_

#include <string>
#include <iomanip>
#include <fstream>
#include <Eigen/Dense>
#include "Metadata.h"

using Eigen::VectorXcd;
using Eigen::VectorXd;
using VectorXld = Eigen::Matrix<long double, Eigen::Dynamic, 1>;
using VectorXcld = Eigen::Matrix<complex<long double>, Eigen::Dynamic, 1>;

// Class which stores a matrix of dimensions dim_B x dim_A
// using template for the Tensor to take in any precision (float?, double, long double, ...)
template <typename RealT>
class Tensor
{
public:
  // Metadata of the tensor stores its dimension, power and number of parameters
  Metadata metadata;
  // Tensor to write data to
  vector<vector<RealT>> T;

  // Null Tensor constructor
  Tensor(const int &q, const array<int, 2> &coef);

  // Overloaded null Tensor constructor
  Tensor(const int &q, const array<int, 2> &coef, array<int, 2> pows);

  // Tensor constructor with 1 value
  Tensor(const RealT &Q, Powers &powers, const int &q, const array<int, 2> &coef);

  // Overloaded Tensor constructor with 1 value
  Tensor(const RealT &Q, Powers &powers, const int &q, const array<int, 2> &coef, const array<int, 2> &pows);

  // Tensor default constructor
  Tensor();

  // Tensor copy constructor
  Tensor(const Tensor &other);

  // Copy assignment operator
  Tensor &operator=(const Tensor &other);

  // Move constructor for Tensor
  Tensor(Tensor &&other) noexcept;

  // Move assignment operator
  Tensor &operator=(Tensor &&other) noexcept;

  // Tensor deconstructor
  ~Tensor();

  // Method to transfor a Tensor to a higher no. cycles or change coefficient
  Tensor transform(const Metadata &meta_other, Powers &powers, const RealT &Q);

  // Overloaded method to transfor a Tensor to a higher no. cycles
  Tensor transform(const int &nq);

  // Method to multiply this tensor and another one to a higher order
  Tensor multiply(const Tensor &other) &;

  // Method to evaluate the polynomial of the Tensor
  template <typename Vec>
  typename Vec::value_type evaluate(const Vec &a_vec, const Vec &b_vec) const;

  // Method to take a single derivative of the tensor
  void derivative(Tensor &dtensor, const bool &param, const int &k);

  // Method to take multiple derivative of the tensor
  Tensor derivatives(vector<int> a_dev, vector<int> b_dev);

  // Method to display the Tensor
  void display_tensor();

  // Method to display the Tensor to RealT precision
  void display_tensor_precise();

  // Method to display the polynomial of the Tensor
  void display_poly();

  // Method to save the Tensor
  void save(const string &dir, const string &label);

  // Method to load the Tensor
  void load(const string &dir, const string &label);

  // Method to check if the Tensor is null
  bool iszero() const;

  // Overload the += operator for Tensor addition
  Tensor &operator+=(const Tensor &other);

  // Overload the -= operator for Tensor addition
  Tensor &operator-=(const Tensor &other);

  // Overload the * operator for multiplication with RealT
  Tensor operator*(const RealT &value);

  // Overload the *= operator for multiplication of Tensors and integers
  Tensor &operator*=(const int &value);
};

#endif // _TENSOR_H_
