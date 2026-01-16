#include "Prefactors.h"

// Specialization for double type
template <>
const double Prefactors<double>::alpha = -1.0 / 24.0;

template <>
const double Prefactors<double>::beta = 1.0 / 12.0;

template <>
const vector<double> Prefactors<double>::gammas = {
    7.0 / 5760.0,
    1.0 / 480.0,
    1.0 / 360.0,
    1.0 / 360.0,
    -1.0 / 120.0,
    -1.0 / 720.0};

template <>
const vector<double> Prefactors<double>::deltas = {
    -31.0 / 967680.0,
    -19.0 / 80640.0,
    23.0 / 161280.0,
    -1.0 / 10080.0,
    -1.0 / 10080.0,
    -1.0 / 20160.0,
    13.0 / 60480.0,
    -1.0 / 5040.0,
    -1.0 / 3360.0,
    -19.0 / 40320.0,
    1.0 / 6720.0,
    -1.0 / 7560.0,
    1.0 / 1008.0,
    -1.0 / 10080.0,
    -1.0 / 10080.0,
    1.0 / 10080.0,
    1.0 / 5040.0,
    1.0 / 30240.0};
    
// Specialization for long double type
template <>
const long double Prefactors<long double>::alpha = -1.0 / 24.0L;

template <>
const long double Prefactors<long double>::beta = 1.0 / 12.0L;

template <>
const vector<long double> Prefactors<long double>::gammas = {
    7.0 / 5760.0L,
    1.0 / 480.0L,
    1.0 / 360.0L,
    1.0 / 360.0L,
    -1.0 / 120.0L,
    -1.0 / 720.0L};

template <>
const vector<long double> Prefactors<long double>::deltas = {
    -31.0 / 967680.0L,
    -19.0 / 80640.0L,
    23.0 / 161280.0L,
    -1.0 / 10080.0L,
    -1.0 / 10080.0L,
    -1.0 / 20160.0L,
    13.0 / 60480.0L,
    -1.0 / 5040.0L,
    -1.0 / 3360.0L,
    -19.0 / 40320.0L,
    1.0 / 6720.0L,
    -1.0 / 7560.0L,
    1.0 / 1008.0L,
    -1.0 / 10080.0L,
    -1.0 / 10080.0L,
    1.0 / 10080.0L,
    1.0 / 5040.0L,
    1.0 / 30240.0L};

// Default constructor for Prefactors
template <typename RealT>
Prefactors<RealT>::Prefactors()
{
}

// Destructor for Prefactors
template <typename RealT>
Prefactors<RealT>::~Prefactors()
{
}

// Explicit instantiation for the template class
template class Prefactors<double>;
template class Prefactors<long double>;
