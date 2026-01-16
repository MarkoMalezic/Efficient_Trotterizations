#ifndef _INOUT_H_
#define _INOUT_H_

#include <chrono>
#include <iostream>
#include <fstream>
#include <memory>
#include <random>
#include "H5Cpp.h"
#include "Exact.h"
#include "Trotter.h"

// Helper function to parse a boolean value
bool parse_bool(const string &str);

// Helper function to parse a real or complex value
template <typename Scalar>
Scalar parse_value(const string &str);

// Helper function to parse a list to extract vectors or arrays
template <typename List>
List parse_list(const string &str);

// Helper function to generate a random vector
template <typename RealT>
Eigen::Matrix<RealT, Eigen::Dynamic, 1> generate_random_vector(int size, const RealT &min = -0.1, const RealT &max = 0.1, const int &seed = 42);

// Class responsible for reading the input file
template <typename Scalar>
class InOut
{
public:
  using RealT = typename Eigen::NumTraits<Scalar>::Real;
  // Read file stream
  ifstream rfile;

  int routine;          // Which routine to use (Time evolution: 0, Error estimation: 1)
  string save_file{""}; // File to save to (if not provided no writing occurs)

  // The model object, which is needed for all routines
  unique_ptr<Model<RealT>> model;

  // InOut constructor
  InOut(const string read_file);

  // InOut deconstructor
  ~InOut();

  // Method to build the model object
  void build_model(const string &model_str);

  // *** Routines ***

  // Time evolution
  // Methods: - exact: time evolve a state using exact diagonalization
  //          - trotter: time evolve a state using Trotter decomposition
  void time_evolve();

  // Error estimation
  // Methods: - operator: compare the time evolution operators from exact diagonalization and Trotter decomposition
  //          - state: compare the time evolved states from exact diagonalization and Trotter decomposition
  void error_estimate();
};

#endif // _INOUT_H_
