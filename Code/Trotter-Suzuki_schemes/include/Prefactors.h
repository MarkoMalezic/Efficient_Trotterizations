#ifndef _PREFACTORS_H_
#define _PREFACTORS_H_

#include <vector>
#include <type_traits>

using namespace std;

// Class which stores the prefactors needed to compute Coefficients
// using template for the Tensor to take in any class (double, ...)
template <typename RealT>
class Prefactors
{
public:
  // Order n=3
  static const RealT alpha;
  static const RealT beta;
  // Order n=5
  static const vector<RealT> gammas;
  // Order n=7
  static const vector<RealT> deltas;

  // Prefactors default constructor
  Prefactors();

  // Prefactors deconstructor
  ~Prefactors();
};

#endif // _PREFACTORS_H_
