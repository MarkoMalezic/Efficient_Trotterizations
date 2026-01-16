#include "Metadata.h"

// Helper functions to get the Tensor powers based on the desired coefficient
array<int, 2> get_pows(const array<int, 2> &coef)
{
  // Map coef to the its powers of ai and bi
  map<array<int, 2>, array<int, 2>> pow_map = {
      // Order n=3:
      {{3, 1}, {1, 2}}, // alpha
      {{3, 2}, {2, 1}}, // beta
      // Order n=5
      {{5, 1}, {1, 4}}, // gamma1
      {{5, 2}, {2, 3}}, // gamma2
      {{5, 3}, {2, 3}}, // gamma3
      {{5, 4}, {3, 2}}, // gamma4
      {{5, 5}, {3, 2}}, // gamma5
      {{5, 6}, {4, 1}}, // gamma6
      // Order n=7
      {{7, 1}, {1, 6}},  // delta1
      {{7, 2}, {2, 5}},  // delta2
      {{7, 3}, {2, 5}},  // delta3
      {{7, 4}, {2, 5}},  // delta4
      {{7, 5}, {3, 4}},  // delta5
      {{7, 6}, {3, 4}},  // delta6
      {{7, 7}, {3, 4}},  // delta7
      {{7, 8}, {3, 4}},  // delta8
      {{7, 9}, {3, 4}},  // delta9
      {{7, 10}, {4, 3}}, // delta10
      {{7, 11}, {4, 3}}, // delta11
      {{7, 12}, {4, 3}}, // delta12
      {{7, 13}, {4, 3}}, // delta13
      {{7, 14}, {4, 3}}, // delta14
      {{7, 15}, {5, 2}}, // delta15
      {{7, 16}, {5, 2}}, // delta16
      {{7, 17}, {5, 2}}, // delta17
      {{7, 18}, {6, 1}},  // delta18
      // Order n=9
      {{9, 1}, {1, 8}},  // epsilon1
      {{9, 2}, {2, 7}},  // epsilon2
      {{9, 3}, {2, 7}},  // epsilon3
      {{9, 4}, {2, 7}},  // epsilon4
      {{9, 5}, {2, 7}},  // epsilon5
      {{9, 6}, {3, 6}},  // epsilon6
      {{9, 7}, {3, 6}},  // epsilon7
      {{9, 8}, {3, 6}},  // epsilon8
      {{9, 9}, {3, 6}},  // epsilon9
      {{9, 10}, {3, 6}}, // epsilon10
      {{9, 11}, {3, 6}}, // epsilon11
      {{9, 12}, {3, 6}}, // epsilon12
      {{9, 13}, {3, 6}}, // epsilon13
      {{9, 14}, {3, 6}}, // epsilon14
      {{9, 15}, {4, 5}}, // epsilon15
      {{9, 16}, {4, 5}}, // epsilon16
      {{9, 17}, {4, 5}}, // epsilon17
      {{9, 18}, {4, 5}},  // epsilon18
      {{9, 19}, {4, 5}},  // epsilon19
      {{9, 20}, {4, 5}},  // epsilon20
      {{9, 21}, {4, 5}},  // epsilon21
      {{9, 22}, {4, 5}},  // epsilon22
      {{9, 23}, {4, 5}},  // epsilon23
      {{9, 24}, {4, 5}},  // epsilon24
      {{9, 25}, {4, 5}},  // epsilon25
      {{9, 26}, {4, 5}},  // epsilon26
      {{9, 27}, {4, 5}},  // epsilon27
      {{9, 28}, {4, 5}},  // epsilon28
      {{9, 29}, {5, 4}},  // epsilon29
      {{9, 30}, {5, 4}},  // epsilon30
      {{9, 31}, {5, 4}},  // epsilon31
      {{9, 32}, {5, 4}},  // epsilon32
      {{9, 33}, {5, 4}},  // epsilon33
      {{9, 34}, {5, 4}},  // epsilon34
      {{9, 35}, {5, 4}},  // epsilon35
      {{9, 36}, {5, 4}},  // epsilon36
      {{9, 37}, {5, 4}},  // epsilon37
      {{9, 38}, {5, 4}},  // epsilon38
      {{9, 39}, {5, 4}},  // epsilon39
      {{9, 40}, {5, 4}},  // epsilon40
      {{9, 41}, {5, 4}},  // epsilon41
      {{9, 42}, {5, 4}},  // epsilon42
      {{9, 43}, {6, 3}},  // epsilon43
      {{9, 44}, {6, 3}},  // epsilon44
      {{9, 45}, {6, 3}},  // epsilon45
      {{9, 46}, {6, 3}},  // epsilon46
      {{9, 47}, {6, 3}},  // epsilon47
      {{9, 48}, {6, 3}},  // epsilon48
      {{9, 49}, {6, 3}},  // epsilon49
      {{9, 50}, {6, 3}},  // epsilon50
      {{9, 51}, {6, 3}},  // epsilon51
      {{9, 52}, {7, 2}},  // epsilon52
      {{9, 53}, {7, 2}},  // epsilon53
      {{9, 54}, {7, 2}},  // epsilon54
      {{9, 55}, {7, 2}},  // epsilon55
      {{9, 56}, {8, 1}}   // epsilon56
  };
  // Find the powers based on the coef pair
  auto it = pow_map.find(coef);
  if (it != pow_map.end())
  {
    return it->second;
  }
  throw invalid_argument("Error: Coefficient not found in the map.");
}

// Constructor for Metadata
Metadata::Metadata(const int &q, const array<int, 2> &coef)
    : q(q), coef(coef)
{
  // The number of parameters is dependent on even/odd number of cycles
  if (q % 2 == 0)
  {
    nps = {q / 2, q / 2 + 1};
  }
  else
  {
    nps = {(q + 1) / 2, (q + 1) / 2};
  }
  pows = get_pows(coef); // Get powers from the corresponding map
  // Calculate the dimensions according to the powers and number of parameters
  dims = {binom(pows[0] + nps[0] - 1, nps[0] - 1), binom(pows[1] + nps[1] - 1, nps[1] - 1)};
}

// Overloaded constructor for Metadata
// Used to initialize the powers directly without using the map
Metadata::Metadata(const int &q, const array<int, 2> &coef, array<int, 2> pows)
    : q(q), coef(coef), pows(pows)
{
  // The number of parameters is dependent on even/odd number of cycles
  if (q % 2 == 0)
  {
    nps = {q / 2, q / 2 + 1};
  }
  else
  {
    nps = {(q + 1) / 2, (q + 1) / 2};
  }
  // Calculate the dimensions according to the powers and number of parameters
  dims = {binom(pows[0] + nps[0] - 1, nps[0] - 1), binom(pows[1] + nps[1] - 1, nps[1] - 1)};
}

// Default constructor for Metadata
Metadata::Metadata()
{
}

// Copy Constructor
Metadata::Metadata(const Metadata &other)
    : q(other.q), coef(other.coef), nps(other.nps), pows(other.pows), dims(other.dims)
{
}

// Copy assignement operator
Metadata &Metadata::operator=(const Metadata &other)
{
  if (this != &other)
  { // Check for self-assignment
    q = other.q;
    coef = other.coef;
    nps = other.nps;
    pows = other.pows;
    dims = other.dims;
  }
  return *this;
}

// Move Constructor
Metadata::Metadata(Metadata &&other) noexcept
    : q(other.q), coef(std::move(other.coef)), nps(std::move(other.nps)),
      pows(std::move(other.pows)), dims(std::move(other.dims))
{
}

// Move assignement operator
Metadata &Metadata::operator=(Metadata &&other) noexcept
{
  if (this != &other)
  { // Check for self-assignment
    q = other.q;
    coef = std::move(other.coef);
    nps = std::move(other.nps);
    pows = std::move(other.pows);
    dims = std::move(other.dims);
  }
  return *this;
}

// Destructor for Metadata
Metadata::~Metadata()
{
}

// Method to display the metadata
void Metadata::display()
{
  cout << "no. cycles q = " << q << ", coefficient label: {" << coef[0] << ", " << coef[1] << "}" << endl;
  cout << "no. parameters {nB, nA} = {" << nps[0] << ", " << nps[1] << "}" << endl;
  cout << "maximal powers {powB, powA} = {" << pows[0] << ", " << pows[1] << "}" << endl;
  cout << "dimensions dimB x dimA = " << dims[0] << " x " << dims[1] << endl;
}
