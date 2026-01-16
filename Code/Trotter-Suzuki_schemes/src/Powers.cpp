#include "Powers.h"

// Helper function to calculate the binomial coefficients
int binom(int n, int k)
{
  if (k > n || k < 0)
    return 0;
  if (k == 0 || k == n)
    return 1;
  k = std::min(k, n - k);
  int c = 1;
  for (int i = 1; i <= k; ++i)
  {
    c = c * (n - k + i) / i;
  }
  return c;
}

// Helper function to get a single lexicographical index
// Input: - pow = vector with powers on the parameters
//        - D = dimension/maximal lexicographical index
//        - P = Total power
int pow_to_lex(const vector<int> &pow, int D, int P)
{
  int N{static_cast<int>(pow.size())};
  int i{0};
  // If the total power is 0 the index is simply 0
  if (P == 0)
  {
    return 0;
    // If the total power is 1 the index a sum
  }
  else if (P == 1)
  {
    while (pow[i] != 1)
    {
      i += 1;
    }
    return i;
  }
  else
  {
    while (P > 0 || N < 0)
    {
      // Otherwise calculate the no. combinations, which come before and subtract them
      for (int j{0}; j < pow[i]; ++j)
      {
        D -= binom(P - j + N - 2, N - 2);
      }
      P -= pow[i];
      i += 1;
      N -= 1;
    }
    return D - 1;
  }
}

// Helper function to transform an inverse lexicographical index to a vector powers
// Input: - inv_lex_ind = The inverse of a lexigraphical index: inv_lex_ind = dim - lex_ind
//        - P = Total power
//        - N = Number of parameters
vector<int> inv_lex_to_pow(int inv_lex_ind, int P, int N)
{
  vector<int> pow(N, 0);
  for (int i = 0; i < N; ++i)
  {
    if (i == N - 1)
    {             // Last position
      pow[i] = P; // Assign all remaining power
    }
    else
    {
      // Calculate the no. combinations that come before and add to the power vector
      for (int power = 0; power <= P; ++power)
      {
        int combs = binom(P - power + N - 2 - i, N - 2 - i);
        if (inv_lex_ind < combs)
        {
          pow[i] = power;
          P -= power;
          break;
        }
        else
        {
          inv_lex_ind -= combs;
        }
      }
    }
  }
  return pow;
}

// Constructor for Powers
Powers::Powers(const vector<int> &a_vec, const vector<int> &b_vec)
    : a_pow(a_vec), b_pow(b_vec)
{
}

// Copy constructor for Powers
Powers::Powers(const Powers &other)
    : a_pow(other.a_pow), b_pow(other.b_pow)
{
}

// Default constructor for Powers
Powers::Powers()
{
}

// Destructor for Powers
Powers::~Powers()
{
}

// Method to display both a_pow and b_pow
void Powers::display() const
{
  cout << "a_pow = (";
  for (size_t powi{0}; powi < a_pow.size(); ++powi)
  {
    cout << a_pow[powi] << (powi != a_pow.size() - 1 ? ", " : "");
  }
  cout << ")" << endl;
  cout << "b_pow = (";
  for (size_t powi{0}; powi < b_pow.size(); ++powi)
  {
    cout << b_pow[powi] << (powi != b_pow.size() - 1 ? ", " : "");
  }
  cout << ")" << endl;
}

// Method to check whether the Powers class is empty
// If either powers of A or B are empty then return true
bool Powers::empty() const
{
  if (a_pow.empty() || b_pow.empty())
  {
    return true;
  }
  else
  {
    return false;
  }
}

// Method to check whether the Powers class is zero
// If both powers of A and B are zero then return true
bool Powers::iszero() const
{
  if (accumulate(a_pow.begin(), a_pow.end(), 0) == 0 && accumulate(b_pow.begin(), b_pow.end(), 0) == 0)
  {
    return true;
  }
  else
  {
    return false;
  }
}

// Method to transform the powers of ai and bi to a lexicographical order
array<const int, 2> Powers::get_lex_inds(array<int, 2> dims, array<int, 2> pows) const
{
  // Get the index i - the b part
  int indi = pow_to_lex(b_pow, dims[0], pows[0]);
  // Get the index j - the a part
  int indj = pow_to_lex(a_pow, dims[1], pows[1]);
  return {indi, indj};
}
