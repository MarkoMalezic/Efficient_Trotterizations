#ifndef _SCHEME_H_
#define _SCHEME_H_

#include <optional>
#include "Coefficients.h"

// Class which constructs a Trotter-Suzuki scheme by recursively iterating from the inside out
template <typename RealT>
class Scheme
{
public:
  int n;                     // Order O(h^(n+1)) of the scheme
  int q;                     // Number of cycles
  Coefficients<RealT> coefs; // Coefficients of all orders

  // indices tracking the steps taken
  int ind_A{0};
  int ind_B{0};

  // Boolean tracking which step to take. false = stepA, true = stepB
  bool swtch{false};

  // Check what to print: - 0 = Nothing
  //                      - 1 = Cycles
  //                      - 2 = Coefficient values
  char verbose;

  // Scheme constructor
  Scheme(const int &order, const int &no_cycles, const char &verbose = 1);

  // Scheme load constructor
  Scheme(const string &directory, const int &order, const int &no_cycles, const char &verbose = 1);

  // Scheme default constructor
  Scheme();

  // Copy constructor
  Scheme(const Scheme& other);

  // Copy assignment operator
  Scheme &operator=(const Scheme &other);

  // Move assignment operator
  Scheme &operator=(Scheme &&other) noexcept;

  // Scheme destructor
  ~Scheme();

  // Methods to update scheme
  void step_A();
  void step_B();
  // Method to iterate through the cycles
  void iterate(optional<int> total_steps_opt = nullopt);

  // Method to transform to a higher number of cycles
  void transform(const int &nq);

  // Method to save the computed Scheme
  void save(const string &directory);

  // Method, which calculates the error for the desired order
  template <typename Vec>
  RealT err2(const int &order, const Vec &a_eval, const Vec &b_eval);

  // Method, which calculates the efficiency for the desired order
  template <typename Vec>
  RealT eff(const int &order, const Vec &a_eval, const Vec &b_eval);

  // Method, which displays squared errors of different order (if known)
  template <typename Vec>
  void display_err2(const Vec &a_eval, const Vec &b_eval);
  
  // Method, which displays efficiencies of different order (if known)
  template <typename Vec>
  void display_eff(const Vec &a_eval, const Vec &b_eval);
};

#endif // _SCHEME_H_
