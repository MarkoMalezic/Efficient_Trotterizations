#include "Multinomial.h"

// Helper function to generate a factorial with possibility to end it earlier
// e.g. factorial(5, 2) = 5*4*3, while factorial(5, 0) = 5! = 5*4*3*2*1
unsigned long factorial(const int &n, const int &end)
{
  if (n < end)
  {
    throw invalid_argument("Error: n must be greater than end.");
  }
  unsigned long result = 1;
  for (int i = end + 1; i <= n; ++i)
  {
    result *= i;
  }
  return result;
}

// Helper function to recursively generate combinations of powers for only the first 'k' parameters
vector<vector<int>> generate_combs(const int &k, const int &n,
                                   int current_sum, int index,
                                   vector<int> &current_comb)
{
  vector<vector<int>> all_combs;

  if (index == k - 1)
  {
    // Last variable, set its value to make the sum equal to n
    current_comb[index] = n - current_sum;
    all_combs.push_back(current_comb);
    return all_combs;
  }

  // Iterate over all possible values for the current variable
  for (int i = 0; i <= n - current_sum; ++i)
  {
    current_comb[index] = i;
    // Combine results from the recursive call
    auto combs = generate_combs(k, n, current_sum + i, index + 1, current_comb);
    all_combs.insert(all_combs.end(), combs.begin(), combs.end());
  }

  return all_combs;
}

// Helper function to calculate the prefactor of a multinomial combination
unsigned long generate_prefac(const int &n, vector<int> comb)
{
  if (comb.empty())
  {
    throw invalid_argument("Error: comb vector is empty.");
  }
  unsigned long comb_prefac = factorial(n, 0);
  for (int val : comb)
  {
    comb_prefac /= factorial(val, 0);
  }
  return comb_prefac;
}

// Constructor for Multinomial
Multinomial::Multinomial(const int &n, const int &m, const int &k)
    : n(n), m(m)
{
  vector<int> current_comb(m, 0);
  // Generate all combinations for the multinomial, allowing powers only for the first 'k' parameters
  combs = generate_combs(k, n, 0, 0, current_comb);
  for (auto &comb : combs)
  {
    comb.resize(m, 0); // Pad with zeros to make the size equal to 'm'
  }
  sort(combs.rbegin(), combs.rend()); // Sort in reverse lexicographical order
  for (const auto &comb : combs)
  {
    prefacs.push_back(generate_prefac(n, comb)); // Calculate the prefactors
  }
}

// Default constructor for Multinomial
Multinomial::Multinomial()
{
}

// Destructor for Multinomial
Multinomial::~Multinomial()
{
}
