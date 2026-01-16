#include "NCoefficients.h"

// Define the static member for double specialization
template <>
const Prefactors<double> NCoefficients<VectorXd>::prefacs = Prefactors<double>();
template <>
const Prefactors<double> NCoefficients<VectorXcd>::prefacs = Prefactors<double>();

// Define the static member for double specialization
template <>
const Prefactors<long double> NCoefficients<VectorXld>::prefacs = Prefactors<long double>();
template <>
const Prefactors<long double> NCoefficients<VectorXcld>::prefacs = Prefactors<long double>();


// Constructor for NCoefficients
template <typename Vec>
NCoefficients<Vec>::NCoefficients(const Vec &a_eval, const Vec &b_eval)
    : a_eval(a_eval), b_eval(b_eval), gammas(6), deltas(18)
{
}

// Default constructor for NCoefficients
template <typename Vec>
NCoefficients<Vec>::NCoefficients()
{
}

// Destructor for NCoefficients
template <typename Vec>
NCoefficients<Vec>::~NCoefficients()
{
}

// Method to update the alpha coefficient with the A operator: exp(A/2) exp(Phi) exp(A/2)
template <typename Vec>
void NCoefficients<Vec>::alpha_stepA(int &ind_A)
{
  // Prefactor initialization
  RealT alpha_Q = prefacs.alpha;
  RealT beta_Q = prefacs.beta;

  // First term
  alpha += alpha_Q * pow(a_eval[ind_A], 2) * sigma;
  // Second term
  alpha -= beta_Q * a_eval[ind_A] * nu * sigma;
}

// Method to update the alpha coefficient with the B operator: exp(B/2) exp(Phi) exp(B/2)
template <typename Vec>
void NCoefficients<Vec>::alpha_stepB(int &ind_B)
{
  // Prefactor initialization
  RealT beta_Q = prefacs.beta;

  // First term
  alpha += beta_Q * b_eval[ind_B] * pow(nu, 2);
}

// Method to update the beta coefficient with the A operator: exp(A/2) exp(Phi) exp(A/2)
template <typename Vec>
void NCoefficients<Vec>::beta_stepA(int &ind_A)
{
  // Prefactor initialization
  RealT beta_Q = prefacs.beta;

  // First term
  beta += beta_Q * a_eval[ind_A] * pow(sigma, 2);
}

// Method to update the beta coefficient with the B operator: exp(B/2) exp(Phi) exp(B/2)
template <typename Vec>
void NCoefficients<Vec>::beta_stepB(int &ind_B)
{
  // Prefactor initialization
  RealT alpha_Q = prefacs.alpha;
  RealT beta_Q = prefacs.beta;

  // First term
  beta += alpha_Q * pow(b_eval[ind_B], 2) * nu;
  // Second term
  beta -= beta_Q * b_eval[ind_B] * nu * sigma;
}

// Method to update the gamma coefficients by switching between operators A, B
template <typename Vec>
void *NCoefficients<Vec>::gammas_step(const int &gamma_ind, int &ind)
{
  // Prefactor initialization
  RealT alpha_Q = prefacs.alpha;
  RealT beta_Q = prefacs.beta;
  vector<RealT> gammas_Q = prefacs.gammas;

  switch (gamma_ind)
  {
  case 1:
  {
    // First term
    gammas[0] += gammas_Q[0] * pow(a_eval[ind], 4) * sigma;

    // Second term
    gammas[0] += (gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * sigma;

    // Third term
    gammas[0] -= (gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * sigma;

    // Fourth term
    gammas[0] -= gammas_Q[5] * a_eval[ind] * pow(nu, 3) * sigma;

    // Fifth term
    gammas[0] += alpha_Q * pow(a_eval[ind], 2) * alpha;

    // Sixth term
    gammas[0] -= beta_Q * a_eval[ind] * nu * alpha;
    break;
  }
  case 2:
  {
    // First term
    gammas[1] += gammas_Q[1] * pow(a_eval[ind], 3) * pow(sigma, 2);

    // Second term
    gammas[1] -= (2*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * pow(sigma, 2);

    // Third term
    gammas[1] -= 2*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * pow(sigma, 2);

    // Fourth term
    gammas[1] -= alpha_Q * pow(a_eval[ind], 2) * beta;

    // Fifth term
    gammas[1] += beta_Q * a_eval[ind] * nu * beta;

    // Sixth term
    gammas[1] += beta_Q * a_eval[ind] * sigma * alpha;
    break;
  }
  case 3:
  {
    // First term
    gammas[2] += gammas_Q[2] * pow(a_eval[ind], 3) * pow(sigma, 2);

    // Second term
    gammas[2] -= gammas_Q[4] * pow(a_eval[ind], 2) * nu * pow(sigma, 2);

    // Third term
    gammas[2] -= gammas_Q[5] * a_eval[ind] * pow(nu, 2) * pow(sigma, 2);

    // Fourth term
    gammas[2] -= 2*beta_Q * a_eval[ind] * sigma * alpha;
    break;
  }
  case 4:
  {
    // First term
    gammas[3] += gammas_Q[3] * pow(a_eval[ind], 2) * pow(sigma, 3);

    // Second term
    gammas[3] += gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 3);

    // Third term
    gammas[3] += beta_Q * a_eval[ind] * sigma * beta;
    break;
  }
  case 5:
  {
    // First term
    gammas[4] += gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 3);

    // Second term
    gammas[4] += 2*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 3);

    // Third term
    gammas[4] -= 2*beta_Q * a_eval[ind] * sigma * beta;
    break;
  }
  case 6:
  {
    // First term
    gammas[5] += gammas_Q[5] * a_eval[ind] * pow(sigma, 4);
    break;
  }
  default:
    throw invalid_argument("Invalid gamma index");
  }
  return nullptr;
}

// Method to update the gamma coefficients by switching between operators A, B
template <typename Vec>
void *NCoefficients<Vec>::deltas_step(const int &delta_ind, int &ind)
{
  // Prefactor initialization
  RealT alpha_Q = prefacs.alpha;
  RealT beta_Q = prefacs.beta;
  vector<RealT> gammas_Q = prefacs.gammas;
  vector<RealT> deltas_Q = prefacs.deltas;

  switch (delta_ind)
  {
  case 1:
  {
    // First term
    deltas[0] += deltas_Q[0] * pow(a_eval[ind], 6) * sigma;

    // Second term
    deltas[0] += (deltas_Q[1] + deltas_Q[2] + deltas_Q[3]) * pow(a_eval[ind], 5) * nu * sigma;

    // Third term
    deltas[0] += (deltas_Q[4] + deltas_Q[5] + deltas_Q[6] + deltas_Q[7] + deltas_Q[8]) * pow(a_eval[ind], 4) * pow(nu, 2) * sigma;

    // Fourth term
    deltas[0] -= (deltas_Q[9] + deltas_Q[10] + deltas_Q[11] + deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * pow(nu, 3) * sigma;

    // Fifth term
    deltas[0] -= (deltas_Q[14] + deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 4) * sigma;

    // Sixth term
    deltas[0] -= deltas_Q[17] * a_eval[ind] * pow(nu, 5) * sigma;

    // Seventh term
    deltas[0] += gammas_Q[0] * pow(a_eval[ind], 4) * alpha;

    // Eighth term
    deltas[0] += (gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * alpha;

    // Ninth term
    deltas[0] -= (gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * alpha;

    // Tenth term
    deltas[0] -= gammas_Q[5] * a_eval[ind] * pow(nu, 3) * alpha;

    // Eleventh term
    deltas[0] += alpha_Q * pow(a_eval[ind], 2) * gammas[0];

    // Twelfth term
    deltas[0] -= beta_Q * a_eval[ind] * nu * gammas[0];
    break;
  }
  case 2:
  {
    // First term
    deltas[1] += deltas_Q[1] * pow(a_eval[ind], 5) * pow(sigma, 2);

    // Second term
    deltas[1] += (3*deltas_Q[4] + 2*deltas_Q[5] + deltas_Q[6] + deltas_Q[7] + 2*deltas_Q[8]) * pow(a_eval[ind], 4) * nu * pow(sigma, 2);

    // Third term
    deltas[1] -= (3*deltas_Q[9] + 4*deltas_Q[10] + 4*deltas_Q[11] + 3*deltas_Q[12] + 2*deltas_Q[13]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 2);

    // Fourth term
    deltas[1] -= (5*deltas_Q[14] + 5*deltas_Q[15] + 4*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 2);

    // Fifth term
    deltas[1] -= 5*deltas_Q[17] * a_eval[ind] * pow(nu, 4) * pow(sigma, 2);

    // Sixth term
    deltas[1] -= 2*gammas_Q[0] * pow(a_eval[ind], 4) * beta;

    // Seventh term
    deltas[1] += gammas_Q[2] * pow(a_eval[ind], 3) * sigma * alpha;

    // Eighth term
    deltas[1] -= 2*(gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * beta;

    // Ninth term
    deltas[1] -= gammas_Q[4] * pow(a_eval[ind], 2) * nu * sigma * alpha;

    // Tenth term
    deltas[1] += 2*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * beta;

    // Eleventh term
    deltas[1] -= gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * alpha;

    // Twelfth term
    deltas[1] += 2*gammas_Q[5] * a_eval[ind] * pow(nu, 3) * beta;

    // Thirteenth term
    deltas[1] += 2*alpha_Q * pow(a_eval[ind], 2) * gammas[1];

    // Fourtheenth term
    deltas[1] += alpha_Q * pow(a_eval[ind], 2) * gammas[2];

    // Fiftheenth term
    deltas[1] -= 2*beta_Q * a_eval[ind] * nu * gammas[1];

    // Sixteenth term
    deltas[1] -= beta_Q * a_eval[ind] * nu * gammas[2];

    // Seventeenth term
    deltas[1] -= beta_Q * a_eval[ind] * pow(alpha, 2);
    break;
  }
  case 3:
  {
    // First term
    deltas[2] += deltas_Q[2] * pow(a_eval[ind], 5) * pow(sigma, 2);

    // Second term
    deltas[2] -= (deltas_Q[4] - deltas_Q[6] + deltas_Q[8]) * pow(a_eval[ind], 4) * nu * pow(sigma, 2);

    // Third term
    deltas[2] += (deltas_Q[10] + 2*deltas_Q[11] + deltas_Q[12]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 2);

    // Fourth term
    deltas[2] += (deltas_Q[14] + 2*deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 2);

    // Fifth term
    deltas[2] += deltas_Q[17] * a_eval[ind] * pow(nu, 4) * pow(sigma, 2);

    // Sixth term
    deltas[2] += gammas_Q[0] * pow(a_eval[ind], 4) * beta;

    // Seventh term
    deltas[2] += (gammas_Q[1] - 2*gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * alpha;

    // Eighth term
    deltas[2] += (gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * nu * beta;

    // Ninth term
    deltas[2] -= (2*gammas_Q[3] - gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * alpha;

    // Tenth term
    deltas[2] -= (gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(nu, 2) * beta;

    // Eleventh term
    deltas[2] -= gammas_Q[5] * a_eval[ind] * pow(nu, 3) * beta;

    // Twelfth term
    deltas[2] -= alpha_Q * pow(a_eval[ind], 2) * gammas[1];

    // Thirteenth term
    deltas[2] += beta_Q * a_eval[ind] * sigma * gammas[0];

    // Fourtheenth term
    deltas[2] += beta_Q * a_eval[ind] * nu * gammas[1];

    // Fiftheenth term
    deltas[2] += 2*beta_Q * a_eval[ind] * pow(alpha, 2);
    break;
  }
  case 4:
  {
    // First term
    deltas[3] += deltas_Q[3] * pow(a_eval[ind], 5) * pow(sigma, 2);

    // Second term
    deltas[3] += (deltas_Q[7] + deltas_Q[8]) * pow(a_eval[ind], 4) * nu * pow(sigma, 2);

    // Third term
    deltas[3] -= (deltas_Q[11] + deltas_Q[12] + deltas_Q[13]) * pow(a_eval[ind], 3) * pow(nu, 2) * pow(sigma, 2);

    // Fourth term
    deltas[3] -= (deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 3) * pow(sigma, 2);

    // Fifth term
    deltas[3] -= deltas_Q[17] * a_eval[ind] * pow(nu, 4) * pow(sigma, 2);

    // Sixth term
    deltas[3] += 2*gammas_Q[2] * pow(a_eval[ind], 3) * sigma * alpha;

    // Seventh term
    deltas[3] -= 2*gammas_Q[4] * pow(a_eval[ind], 2) * nu * sigma * alpha;

    // Eighth term
    deltas[3] -= 2*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * alpha;

    // Ninth term
    deltas[3] -= 2*beta_Q * a_eval[ind] * sigma * gammas[0];

    // Tenth term
    deltas[3] -= beta_Q * a_eval[ind] * pow(alpha, 2);
    break;
  }
  case 5:
  {
    // First term
    deltas[4] += deltas_Q[4] * pow(a_eval[ind], 4) * pow(sigma, 3);

    // Second term
    deltas[4] -= (deltas_Q[9] + 3*deltas_Q[10] + 3*deltas_Q[11] + deltas_Q[12]) * pow(a_eval[ind], 3) * nu * pow(sigma, 3);

    // Third term
    deltas[4] -= (5*deltas_Q[14] + 5*deltas_Q[15] + 3*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 3);

    // Fourth term
    deltas[4] -= 5*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 3);

    // Fifth term
    deltas[4] += gammas_Q[1] * pow(a_eval[ind], 3) * sigma * beta;

    // Sixth term
    deltas[4] -= (gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * alpha;

    // Seventh term
    deltas[4] -= (2*gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * beta;

    // Eighth term
    deltas[4] -= 3*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * alpha;

    // Ninth term
    deltas[4] -= 2*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * beta;

    // Tenth term
    deltas[4] -= 3*alpha_Q * pow(a_eval[ind], 2) * gammas[3];

    // Eleventh term
    deltas[4] -= alpha_Q * pow(a_eval[ind], 2) * gammas[4];

    // Twelfth term
    deltas[4] += 3*beta_Q * a_eval[ind] * nu * gammas[3];

    // Thirteenth term
    deltas[4] += beta_Q * a_eval[ind] * nu * gammas[4];

    // Fourtheenth term
    deltas[4] += beta_Q * a_eval[ind] * alpha * beta;
    break;
  }
  case 6:
  {
    // First term
    deltas[5] += deltas_Q[5] * pow(a_eval[ind], 4) * pow(sigma, 3);

    // Second term
    deltas[5] -= (deltas_Q[9] - deltas_Q[10] - 3*deltas_Q[11] + deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 3);

    // Third term
    deltas[5] += (deltas_Q[14] + 3*deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 3);

    // Fourth term
    deltas[5] += deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 3);

    // Fifth term
    deltas[5] -= 3*gammas_Q[1] * pow(a_eval[ind], 3) * sigma * beta;

    // Sixth term
    deltas[5] += (3*gammas_Q[3] + 2*gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * alpha;

    // Seventh term
    deltas[5] += (6*gammas_Q[3] + 3*gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * beta;

    // Eighth term
    deltas[5] += 7*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * alpha;

    // Ninth term
    deltas[5] += 6*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * beta;

    // Tenth term
    deltas[5] += 3*alpha_Q * pow(a_eval[ind], 2) * gammas[3];

    // Eleventh term
    deltas[5] += beta_Q * a_eval[ind] * sigma * gammas[1];

    // Twelfth term
    deltas[5] -= 3*beta_Q * a_eval[ind] * nu * gammas[3];

    // Thirteenth term
    deltas[5] -= 2*beta_Q * a_eval[ind] * alpha * beta;
    break;
  }
  case 7:
  {
    // First term
    deltas[6] += deltas_Q[6] * pow(a_eval[ind], 4) * pow(sigma, 3);

    // Second term
    deltas[6] -= (deltas_Q[9] + deltas_Q[10] + deltas_Q[11]) * pow(a_eval[ind], 3) * nu * pow(sigma, 3);

    // Third term
    deltas[6] -= (2*deltas_Q[14] + deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 3);

    // Fourth term
    deltas[6] -= 2*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 3);

    // Fifth term
    deltas[6] += (gammas_Q[1] + gammas_Q[2]) * pow(a_eval[ind], 3) * sigma * beta;

    // Sixth term
    deltas[6] -= 3*gammas_Q[3] * pow(a_eval[ind], 2) * pow(sigma, 2) * alpha;

    // Seventh term
    deltas[6] -= 2*(gammas_Q[3] + gammas_Q[4]) * pow(a_eval[ind], 2) * nu * sigma * beta;

    // Eighth term
    deltas[6] -= 3*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * alpha;

    // Ninth term
    deltas[6] -= 3*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * beta;

    // Tenth term
    deltas[6] -= alpha_Q * pow(a_eval[ind], 2) * gammas[3];

    // Eleventh term
    deltas[6] += beta_Q * a_eval[ind] * sigma * gammas[2];

    // Twelfth term
    deltas[6] += beta_Q * a_eval[ind] * nu * gammas[3];

    // Thirteenth term
    deltas[6] -= beta_Q * a_eval[ind] * alpha * beta;
    break;
  }
  case 8:
  {
    // First term
    deltas[7] += deltas_Q[7] * pow(a_eval[ind], 4) * pow(sigma, 3);

    // Second term
    deltas[7] -= (deltas_Q[12] + 2*deltas_Q[13]) * pow(a_eval[ind], 3) * nu * pow(sigma, 3);

    // Third term
    deltas[7] -= (deltas_Q[15] + 2*deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 3);

    // Fourth term
    deltas[7] -= 3*deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 3);

    // Fifth term
    deltas[7] -= 2*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * alpha;

    // Sixth term
    deltas[7] -= 4*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * alpha;

    // Seventh term
    deltas[7] -= 2*beta_Q * a_eval[ind] * sigma * gammas[2];
    break;
  }
  case 9:
  {
    // First term
    deltas[8] += deltas_Q[8] * pow(a_eval[ind], 4) * pow(sigma, 3);

    // Second term
    deltas[8] -= (2*deltas_Q[11] + deltas_Q[12]) * pow(a_eval[ind], 3) * nu * pow(sigma, 3);

    // Third term
    deltas[8] -= (2*deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * pow(nu, 2) * pow(sigma, 3);

    // Fourth term
    deltas[8] -= deltas_Q[17] * a_eval[ind] * pow(nu, 3) * pow(sigma, 3);

    // Fifth term
    deltas[8] -= 2*gammas_Q[2] * pow(a_eval[ind], 3) * sigma * beta;

    // Sixth term
    deltas[8] += 2*gammas_Q[4] * pow(a_eval[ind], 2) * nu * sigma * beta;

    // Seventh term
    deltas[8] += 2*gammas_Q[5] * a_eval[ind] * pow(nu, 2) * sigma * beta;

    // Eighth term
    deltas[8] -= 2*beta_Q * a_eval[ind] * sigma * gammas[1];

    // Ninth term
    deltas[8] += 2*beta_Q * a_eval[ind] * alpha * beta;
    break;
  }
  case 10:
  {
    // First term
    deltas[9] += deltas_Q[9] * pow(a_eval[ind], 3) * pow(sigma, 4);

    // Second term
    deltas[9] += (deltas_Q[14] - deltas_Q[15]) * pow(a_eval[ind], 2) * nu * pow(sigma, 4);

    // Third term
    deltas[9] += deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    deltas[9] -= (2*gammas_Q[3] - gammas_Q[4]) * pow(a_eval[ind], 2) * pow(sigma, 2) * beta;

    // Fifth term
    deltas[9] -= gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * alpha;

    // Sixth term
    deltas[9] -= alpha_Q * pow(a_eval[ind], 2) * gammas[5];

    // Seventh term
    deltas[9] += beta_Q * a_eval[ind] * sigma * gammas[4];

    // Eighth term
    deltas[9] += beta_Q * a_eval[ind] * nu * gammas[5];

    // Ninth term
    deltas[9] -= beta_Q * a_eval[ind] * pow(beta, 2);
    break;
  }
  case 11:
  {
    // First term
    deltas[10] += deltas_Q[10] * pow(a_eval[ind], 3) * pow(sigma, 4);

    // Second term
    deltas[10] += (3*deltas_Q[14] + 2*deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 4);

    // Third term
    deltas[10] += 3*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    deltas[10] += gammas_Q[3] * pow(a_eval[ind], 2) * pow(sigma, 2) * beta;

    // Fifth term
    deltas[10] += gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * beta;

    // Sixth term
    deltas[10] += 2*alpha_Q * pow(a_eval[ind], 2) * gammas[5];

    // Seventh term
    deltas[10] += beta_Q * a_eval[ind] * sigma * gammas[3];

    // Eighth term
    deltas[10] -= 2*beta_Q * a_eval[ind] * nu * gammas[5];
    break;
  }
  case 12:
  {
    // First term
    deltas[11] += deltas_Q[11] * pow(a_eval[ind], 3) * pow(sigma, 4);

    // Second term
    deltas[11] += (deltas_Q[15] + deltas_Q[16]) * pow(a_eval[ind], 2) * nu * pow(sigma, 4);

    // Third term
    deltas[11] += 2*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    deltas[11] += gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * beta;

    // Fifth term
    deltas[11] += 4*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * alpha;

    // Sixth term
    deltas[11] += 2*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * beta;

    // Seventh term
    deltas[11] -= 2*beta_Q * a_eval[ind] * sigma * gammas[3];
    break;
  }
  case 13:
  {
    // First term
    deltas[12] += deltas_Q[12] * pow(a_eval[ind], 3) * pow(sigma, 4);

    // Second term
    deltas[12] += 2*deltas_Q[15] * pow(a_eval[ind], 2) * nu * pow(sigma, 4);

    // Third term
    deltas[12] -= deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    deltas[12] -= 4*gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * beta;

    // Fifth term
    deltas[12] -= 8*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * alpha;

    // Sixth term
    deltas[12] -= 8*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * beta;

    // Seventh term
    deltas[12] -= 2*beta_Q * a_eval[ind] * sigma * gammas[4];

    // Ninth term
    deltas[12] += 2*beta_Q * a_eval[ind] * pow(beta, 2);
    break;
  }
  case 14:
  {
    // First term
    deltas[13] += deltas_Q[13] * pow(a_eval[ind], 3) * pow(sigma, 4);

    // Second term
    deltas[13] += 2*deltas_Q[16] * pow(a_eval[ind], 2) * nu * pow(sigma, 4);

    // Third term
    deltas[13] += 5*deltas_Q[17] * a_eval[ind] * pow(nu, 2) * pow(sigma, 4);

    // Fourth term
    deltas[13] += gammas_Q[4] * pow(a_eval[ind], 2) * pow(sigma, 2) * beta;

    // Fifth term
    deltas[13] += 6*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * alpha;

    // Sixth term
    deltas[13] += 2*gammas_Q[5] * a_eval[ind] * nu * pow(sigma, 2) * beta;

    // Seventh term
    deltas[13] -= beta_Q * a_eval[ind] * pow(beta, 2);
    break;
  }
  case 15:
  {
    // First term
    deltas[14] += deltas_Q[14] * pow(a_eval[ind], 2) * pow(sigma, 5);

    // Second term
    deltas[14] += deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 5);

    // Third term
    deltas[14] += gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * beta;

    // Fourth term
    deltas[14] += beta_Q * a_eval[ind] * sigma * gammas[5];
    break;
  }
  case 16:
  {
    // First term
    deltas[15] += deltas_Q[15] * pow(a_eval[ind], 2) * pow(sigma, 5);

    // Second term
    deltas[15] -= deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 5);

    // Third term
    deltas[15] -= 2*beta_Q * a_eval[ind] * sigma * gammas[5];
    break;
  }
  case 17:
  {
    // First term
    deltas[16] += deltas_Q[16] * pow(a_eval[ind], 2) * pow(sigma, 5);

    // Second term
    deltas[16] += 5*deltas_Q[17] * a_eval[ind] * nu * pow(sigma, 5);

    // Third term
    deltas[16] -= 2*gammas_Q[5] * a_eval[ind] * pow(sigma, 3) * beta;
    break;
  }
  case 18:
  {
    // First term
    deltas[17] += deltas_Q[17] * a_eval[ind] * pow(sigma, 6);
    break;
  }
  default:
    throw invalid_argument("Invalid gamma index");
  }
  return nullptr;
}


// Explicit instantiation for the template class
template class NCoefficients<VectorXd>;
template class NCoefficients<VectorXcd>;
template class NCoefficients<VectorXld>;
template class NCoefficients<VectorXcld>;
