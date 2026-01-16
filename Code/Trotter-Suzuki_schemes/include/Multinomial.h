#ifndef _MULTINOMIAL_H_
#define _MULTINOMIAL_H_

#include <vector>
#include <algorithm>
#include <iterator>

using namespace std;

// Function to generate a factorial with possibility to end it earlier
unsigned long factorial(const int &n, const int &end);

// Helper function to recursively generate combinations of powers
// excluding the combinations with parameter > k
vector<vector<int>> generate_combs(const int &m, const int &n,
                                   int current_sum, int index,
                                   vector<int> &current_comb, int k);

// Function to calculate the prefactor of a combination
unsigned long generate_prefac(const int &n, vector<int> comb);

class Multinomial
{
public:
  // Parameters n and m of a multinomial: (x_1 + x_2 + ... + x_m)^n
  int n;
  int m;
  // All combinations of powers and their prefactors
  vector<vector<int>> combs;
  vector<int> prefacs{};

  // Multinomial constructor
  Multinomial(const int &n, const int &m, const int &k);

  // Multinomial default constructor
  Multinomial();

  // Multinomial deconstructor
  ~Multinomial();
};

#endif // _MULTINOMIAL_H_
