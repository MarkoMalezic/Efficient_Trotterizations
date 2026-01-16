#ifndef _POWERS_H_
#define _POWERS_H_

#include <vector>
#include <array>
#include <numeric>
#include <iostream>

using namespace std;

// Helper function to get the binomial coefficient
int binom(int n, int k);

// Helper function to transform the vector of powers to a lexicographical index
int pow_to_lex(const vector<int> &pow, int D, int P);

// Helper function to transform a lexicographical index to a vector of powers
vector<int> inv_lex_to_pow(int inv_lex_ind, int P, int N);

// Class which stores powers of a and b
class Powers
{
public:
  // Powers of ai and bi parameters
  vector<int> a_pow;
  vector<int> b_pow;

  // Powers constructor
  Powers(const vector<int> &a_vec, const vector<int> &b_vec);

  // Powers copy constructor
  Powers(const Powers &other);

  // Powers default constructor
  Powers();

  // Powers deconstructor
  ~Powers();

  // Method to display both a_pow and b_pow
  void display() const;

  // Method to check whether the class is empty
  bool empty() const;

  // Method to check whether the class is zero
  bool iszero() const;

  // Method to transform the powers of ai and bi to a lexicographical order
  array<const int, 2> get_lex_inds(array<int, 2> dims, array<int, 2> pows) const;
};

#endif // _POWERS_H_
