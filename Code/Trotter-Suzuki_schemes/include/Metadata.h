#ifndef _METADATA_H_
#define _METADATA_H_

#include <map>
#include <stdexcept>
#include "Powers.h"

// Map to find the powers
map<array<int, 2>, array<int, 2>> get_map();

// Map from coef to dimensions and powers of the Tensor class
array<int, 2> get_pows(const array<int, 2> &coef);

// Class which holds important information of the Tensor class
class Metadata
{
public:
  int q;                // Number of cycles
  array<int, 2> coef;   // The coefficient id: e.g. alpha = {3, 1}
  array<int, 2> nps{};  // Number of parameters ai and bi: {nB, nA}
  array<int, 2> pows{}; // powers of the tensor: {pow_B, pow_A}
  array<int, 2> dims{}; // dimensions of the tensor: {dim_B, dim_A}

  // Metadata constructor
  Metadata(const int &q, const array<int, 2> &coef);

  // Overloaded Metadata constructor
  Metadata(const int &q, const array<int, 2> &coef, array<int, 2> pows);

  // Metadata default constructor
  Metadata();

  // Metadata copy constructor
  Metadata(const Metadata &other);

  // Copy assignment operator
  Metadata &operator=(const Metadata &other);

  // Move constructor for Metadata
  Metadata(Metadata &&other) noexcept;

  // Move assignment operator
  Metadata &operator=(Metadata &&other) noexcept;

  // Metadata deconstructor
  ~Metadata();

  // Method to display the metadata
  void display();
};

#endif // _METADATA_H_
