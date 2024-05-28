/*----------------------------------------------------------------------*/
/*! \file
\brief Implementation of an isotropic coupled material depending on the first Cauchy-Green invariant
and the Jacobi determinant

\level 1
*/
/*----------------------------------------------------------------------*/

#include "4C_matelast_coup13apow.hpp"

#include "4C_material_input_base.hpp"

FOUR_C_NAMESPACE_OPEN

MAT::ELASTIC::PAR::Coup13aPow::Coup13aPow(const Teuchos::RCP<CORE::MAT::PAR::Material>& matdata)
    : Parameter(matdata),
      c_(matdata->Get<double>("C")),
      d_(matdata->Get<int>("D")),
      a_(matdata->Get<double>("A"))
{
}

MAT::ELASTIC::Coup13aPow::Coup13aPow(MAT::ELASTIC::PAR::Coup13aPow* params) : params_(params) {}

void MAT::ELASTIC::Coup13aPow::AddStrainEnergy(double& psi, const CORE::LINALG::Matrix<3, 1>& prinv,
    const CORE::LINALG::Matrix<3, 1>& modinv, const CORE::LINALG::Matrix<6, 1>& glstrain,
    const int gp, const int eleGID)
{
  // material Constants
  const double c = params_->c_;
  const int d = params_->d_;
  const double a = params_->a_;

  // strain energy: Psi = c (I_{\boldsymbol{C}}*(III_{\boldsymbol{C}}^(-a))-3)^d
  // add to overall strain energy
  psi += c * pow((prinv(0) * (pow(prinv(2), -a)) - 3.), d);
}

void MAT::ELASTIC::Coup13aPow::add_derivatives_principal(CORE::LINALG::Matrix<3, 1>& dPI,
    CORE::LINALG::Matrix<6, 1>& ddPII, const CORE::LINALG::Matrix<3, 1>& prinv, const int gp,
    const int eleGID)
{
  const double c = params_->c_;
  const int d = params_->d_;
  const double a = params_->a_;

  // If d<2 the material model is not stress free in the reference configuration
  if (d < 2)
    FOUR_C_THROW(
        "The Elast_Coup13aPow - material only works for positive integer exponents, which are "
        "larger than two.");
  // a==0 not necessary to implement extra as equal to coup1pow
  if (a == 0) FOUR_C_THROW("Use Elast_Coup1Pow.");
  // Material not implemented for negative exponents a
  if (a < 0) FOUR_C_THROW("Use positive values for the exponent a.");
  // Restriction for a (due to polyconvexity)
  if (a > (1. / 3. - 1. / (2. * d)))
    std::cout << "\nWARNING: A should be smaller then 1./3.- 1/(2*D). And D should be odd."
              << std::endl;

  const double I1I3a3 = prinv(0) * pow(prinv(2), -a) - 3.;

  dPI(0) += c * d * pow(prinv(2), -a) * pow(I1I3a3, d - 1.);
  dPI(2) += -a * c * d * prinv(0) * pow(prinv(2), -a - 1.) * pow(I1I3a3, d - 1.);

  if (d == 2)
  {
    ddPII(0) += c * d * (d - 1.) * pow(prinv(2), -2. * a);
    ddPII(2) += a * (a + 1.) * c * d * prinv(0) * pow(prinv(2), -a - 2.) * pow(I1I3a3, d - 1.) +
                a * a * c * d * (d - 1.) * prinv(0) * prinv(0) * pow(prinv(2), (-2. * a - 2.));
    ddPII(4) += -a * c * d * pow(prinv(2), -a - 1.) * pow(I1I3a3, d - 1.) -
                a * c * d * (d - 1.) * prinv(0) * pow(prinv(2), (-2. * a - 1.));
  }
  else
  {
    ddPII(0) += c * d * (d - 1.) * pow(prinv(2), -2. * a) * pow(I1I3a3, d - 2.);
    ddPII(2) += a * (a + 1.) * c * d * prinv(0) * pow(prinv(2), -a - 2.) * pow(I1I3a3, d - 1.) +
                a * a * c * d * (d - 1.) * prinv(0) * prinv(0) * pow(prinv(2), (-2. * a - 2.)) *
                    pow(I1I3a3, d - 2.);
    ddPII(4) +=
        -a * c * d * pow(prinv(2), -a - 1.) * pow(I1I3a3, d - 1.) -
        a * c * d * (d - 1.) * prinv(0) * pow(prinv(2), (-2. * a - 1.)) * pow(I1I3a3, d - 2.);
  }
}
FOUR_C_NAMESPACE_CLOSE