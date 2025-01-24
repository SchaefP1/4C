// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_my_expansion_test.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_global_data.hpp"
#include "4C_linalg_fixedsizematrix_tensor_products.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_mat_elast_summand.hpp"
#include "4C_mat_elasthyper_service.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_service.hpp"

#include <Epetra_SerialDenseSolver.h>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::MyExpansionTest_ElastHyper::MyExpansionTest_ElastHyper(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      nummat_elastiniso_(matdata.parameters.get<int>("NUMMATEL3D")),
      matids_elastiniso_(matdata.parameters.get<std::vector<int>>("MATIDSEL3D")),
      density_(matdata.parameters.get<double>("DENS")),
      exp_rate_(matdata.parameters.get<double>("EXPRATE")),
      rho_0_(matdata.parameters.get<double>("rho_0")),
      rho_ss_(matdata.parameters.get<double>("rho_ss"))
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::Mat::Material> Mat::PAR::MyExpansionTest_ElastHyper::create_material()
{
  return std::make_shared<Mat::MyExpansionTest_ElastHyper>(this);
}

Mat::MyExpansionTest_ElastHyperType Mat::MyExpansionTest_ElastHyperType::instance_;


Core::Communication::ParObject* Mat::MyExpansionTest_ElastHyperType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Mat::MyExpansionTest_ElastHyper* gr_elhy = new Mat::MyExpansionTest_ElastHyper();
  gr_elhy->unpack(buffer);

  return gr_elhy;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::MyExpansionTest_ElastHyper::MyExpansionTest_ElastHyper()
    : params_(nullptr), potsumeliso_(0), rho_average_(0)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::MyExpansionTest_ElastHyper::MyExpansionTest_ElastHyper(
    Mat::PAR::MyExpansionTest_ElastHyper* params)
    : params_(params), potsumeliso_(0), rho_average_(0)
{
  // make sure the referenced materials in material list have quick access parameters
  std::vector<int>::const_iterator m;

  // 3d Elastin matrix
  for (m = params_->matids_elastiniso_.begin(); m != params_->matids_elastiniso_.end(); ++m)
  {
    const int matid = *m;
    std::shared_ptr<Mat::Elastic::Summand> sum = Mat::Elastic::Summand::factory(matid);
    if (sum == nullptr) FOUR_C_THROW("Failed to allocate");
    potsumeliso_.push_back(sum);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::MyExpansionTest_ElastHyper::pack(Core::Communication::PackBuffer& data) const
{
  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);
  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  add_to_pack(data, rho_average_);

  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    // loop map of associated potential summands
    for (const auto& p : potsumeliso_) p->pack_summand(data);
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::MyExpansionTest_ElastHyper::unpack(Core::Communication::UnpackBuffer& buffer)
{
  // make sure we have a pristine material
  potsumeliso_.clear();

  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // matid
  int matid;
  extract_from_pack(buffer, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != nullptr)
  {
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::MyExpansionTest_ElastHyper*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }
  }

  extract_from_pack(buffer, rho_average_);

  if (params_ != nullptr)  // summands are not accessible in postprocessing mode
  {
    // make sure the referenced materials in material list have quick access parameters
    std::vector<int>::const_iterator m;

    for (m = params_->matids_elastiniso_.begin(); m != params_->matids_elastiniso_.end(); ++m)
    {
      const int matid = *m;
      std::shared_ptr<Mat::Elastic::Summand> sum = Mat::Elastic::Summand::factory(matid);
      if (sum == nullptr) FOUR_C_THROW("Failed to allocate");
      potsumeliso_.push_back(sum);
    }

    // loop map of associated potential summands
    for (auto& p : potsumeliso_)
    {
      p->unpack_summand(buffer);
    }

    // in the postprocessing mode, we do not unpack everything we have packed
    // -> position check cannot be done in this case
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::MyExpansionTest_ElastHyper::setup(
    int numgp, const Core::IO::InputParameterContainer& container)
{
  // 3D elastin matrix
  for (auto& p : potsumeliso_) p->setup(numgp, container);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::MyExpansionTest_ElastHyper::update()
{
  // loop map of associated potential summands

  // 3D elastin matrix
  for (auto& p : potsumeliso_) p->update();
}

void Mat::MyExpansionTest_ElastHyper::evaluate_cauchy_n_dir_and_derivatives(
    const Core::LinAlg::Matrix<3, 3>& defgrd, const Core::LinAlg::Matrix<3, 1>& n,
    const Core::LinAlg::Matrix<3, 1>& dir, double& cauchy_n_dir,
    Core::LinAlg::Matrix<3, 1>* d_cauchyndir_dn, Core::LinAlg::Matrix<3, 1>* d_cauchyndir_ddir,
    Core::LinAlg::Matrix<9, 1>* d_cauchyndir_dF, Core::LinAlg::Matrix<9, 9>* d2_cauchyndir_dF2,
    Core::LinAlg::Matrix<9, 3>* d2_cauchyndir_dF_dn,
    Core::LinAlg::Matrix<9, 3>* d2_cauchyndir_dF_ddir, int gp, int eleGID,
    const double* concentration, const double* temp, double* d_cauchyndir_dT,
    Core::LinAlg::Matrix<9, 1>* d2_cauchyndir_dF_dT)
{
  // reset sigma contracted with n and dir
  cauchy_n_dir = 0.0;
  double temperature = 293.15;
  double temperature_0 = 293.15;
  // double rho_s = 0.0;
  double rho_ss = params_->rho_ss_;
  double rho_s0 = params_->rho_0_;

  double alpha_temp = 13.0e-6;  // this is the value for steel just to test
  double alpha_reaction = params_->exp_rate_;
  double scaling_factor = (1.0 + alpha_temp * (temperature - temperature_0)) *
                          (1.0 + alpha_reaction * (rho_average_ - rho_s0) / (rho_ss - rho_s0));

  // std::cout << "evaluate_cauchy_n_dir_and_derivatives " << rho_average_ << " " << scaling_factor
  // << std::endl;
  static Core::LinAlg::Matrix<6, 1> idV(true);
  for (int i = 0; i < 3; ++i) idV(i) = 1.0;
  static Core::LinAlg::Matrix<3, 3> idM(true);
  for (int i = 0; i < 3; ++i) idM(i, i) = 1.0;
  static Core::LinAlg::Matrix<3, 3> iFinM(true);
  for (int i = 0; i < 3; ++i) iFinM(i, i) = 1.0;
  iFinM.scale(1 / scaling_factor);
  // TODO I need my inelasatic stuff here
  // inelastic_->evaluate_inverse_inelastic_def_grad(&defgrd, iFinM);
  static Core::LinAlg::Matrix<3, 3> FeM(true);
  FeM.multiply_nn(1.0, defgrd, iFinM, 0.0);

  // get elastic left cauchy-green tensor and corresponding principal invariants
  static Core::LinAlg::Matrix<3, 3> beM(true);
  beM.multiply_nt(1.0, FeM, FeM, 0.0);
  static Core::LinAlg::Matrix<6, 1> beV_strain(true);
  Core::LinAlg::Voigt::Strains::matrix_to_vector(beM, beV_strain);
  static Core::LinAlg::Matrix<3, 1> prinv(true);
  Core::LinAlg::Voigt::Strains::invariants_principal(prinv, beV_strain);
  static Core::LinAlg::Matrix<6, 1> beV_stress(true);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(beM, beV_stress);

  static Core::LinAlg::Matrix<3, 1> beMdn(true);
  beMdn.multiply(1.0, beM, n, 0.0);
  const double beMdnddir = beMdn.dot(dir);
  static Core::LinAlg::Matrix<3, 1> beMddir(true);
  beMddir.multiply(1.0, beM, dir, 0.0);

  static Core::LinAlg::Matrix<3, 3> ibeM(true);
  ibeM.invert(beM);
  static Core::LinAlg::Matrix<6, 1> ibeV_stress(true);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(ibeM, ibeV_stress);
  static Core::LinAlg::Matrix<3, 1> ibeMdn(true);
  ibeMdn.multiply(1.0, ibeM, n, 0.0);
  const double ibeMdnddir = ibeMdn.dot(dir);
  static Core::LinAlg::Matrix<3, 1> ibeMddir(true);
  ibeMddir.multiply(1.0, ibeM, dir, 0.0);

  // derivatives of principle invariants of elastic left cauchy-green tensor
  static Core::LinAlg::Matrix<3, 1> dPI(true);
  static Core::LinAlg::Matrix<6, 1> ddPII(true);
  evaluate_invariant_derivatives(prinv, gp, eleGID, dPI, ddPII);

  const double detFe = FeM.determinant();
  const double nddir = n.dot(dir);
  const double prefac = 2.0 / detFe;

  // calculate \mat{\sigma} \cdot \vec{n} \cdot \vec{v}
  cauchy_n_dir = prefac * (prinv(1) * dPI(1) * nddir + prinv(2) * dPI(2) * nddir +
                              dPI(0) * beMdnddir - prinv(2) * dPI(1) * ibeMdnddir);

  if (d_cauchyndir_dn)
  {
    d_cauchyndir_dn->update(prinv(1) * dPI(1) + prinv(2) * dPI(2), dir, 0.0);
    d_cauchyndir_dn->update(dPI(0), beMddir, 1.0);
    d_cauchyndir_dn->update(-prinv(2) * dPI(1), ibeMddir, 1.0);
    d_cauchyndir_dn->scale(prefac);
  }

  if (d_cauchyndir_ddir)
  {
    d_cauchyndir_ddir->update(prinv(1) * dPI(1) + prinv(2) * dPI(2), n, 0.0);
    d_cauchyndir_ddir->update(dPI(0), beMdn, 1.0);
    d_cauchyndir_ddir->update(-prinv(2) * dPI(1), ibeMdn, 1.0);
    d_cauchyndir_ddir->scale(prefac);
  }

  if (d_cauchyndir_dF)
  {
    static Core::LinAlg::Matrix<6, 1> d_I1_be(true);
    d_I1_be = idV;
    static Core::LinAlg::Matrix<6, 1> d_I2_be(true);
    d_I2_be.update(prinv(0), idV, -1.0, beV_stress);
    static Core::LinAlg::Matrix<6, 1> d_I3_be(true);
    d_I3_be.update(prinv(2), ibeV_stress, 0.0);

    // calculation of \partial b_{el} / \partial F (elastic left cauchy-green w.r.t. deformation
    // gradient)
    static Core::LinAlg::Matrix<6, 9> d_be_dFe(true);
    d_be_dFe.clear();
    Core::LinAlg::Tensor::add_right_non_symmetric_holzapfel_product_strain_like(
        d_be_dFe, idM, FeM, 1.0);
    static Core::LinAlg::Matrix<9, 9> d_Fe_dF(true);
    d_Fe_dF.clear();
    Core::LinAlg::Tensor::add_non_symmetric_product(1.0, idM, iFinM, d_Fe_dF);
    static Core::LinAlg::Matrix<6, 9> d_be_dF(true);
    d_be_dF.multiply(1.0, d_be_dFe, d_Fe_dF, 0.0);

    // calculation of \partial I_i / \partial F (Invariants of b_{el} w.r.t. deformation gradient)
    static Core::LinAlg::Matrix<9, 1> d_I1_dF(true);
    static Core::LinAlg::Matrix<9, 1> d_I2_dF(true);
    static Core::LinAlg::Matrix<9, 1> d_I3_dF(true);
    d_I1_dF.multiply_tn(1.0, d_be_dF, d_I1_be, 0.0);
    d_I2_dF.multiply_tn(1.0, d_be_dF, d_I2_be, 0.0);
    d_I3_dF.multiply_tn(1.0, d_be_dF, d_I3_be, 0.0);

    // add d_cauchyndir_dI1 \odot d_I1_dF and clear static matrix
    d_cauchyndir_dF->update(prefac * (prinv(1) * ddPII(5) * nddir + prinv(2) * ddPII(4) * nddir +
                                         ddPII(0) * beMdnddir - prinv(2) * ddPII(5) * ibeMdnddir),
        d_I1_dF, 0.0);
    // add d_cauchyndir_dI2 \odot d_I2_dF
    d_cauchyndir_dF->update(
        prefac * (dPI(1) * nddir + prinv(1) * ddPII(1) * nddir + prinv(2) * ddPII(3) * nddir +
                     ddPII(5) * beMdnddir - prinv(2) * ddPII(1) * ibeMdnddir),
        d_I2_dF, 1.0);
    // add d_cauchyndir_dI3 \odot d_I3_dF
    d_cauchyndir_dF->update(
        prefac * (prinv(1) * ddPII(3) * nddir + dPI(2) * nddir + prinv(2) * ddPII(2) * nddir +
                     ddPII(4) * beMdnddir - dPI(1) * ibeMdnddir - prinv(2) * ddPII(3) * ibeMdnddir),
        d_I3_dF, 1.0);

    // next three updates add partial derivative of snt w.r.t. the deformation gradient F for
    // constant invariants first part is term arising from \partial Je^{-1} / \partial F
    static Core::LinAlg::Matrix<3, 3> iFeM(true);
    static Core::LinAlg::Matrix<3, 3> iFeTM(true);
    iFeM.invert(FeM);
    iFeTM.update_t(1.0, iFeM, 0.0);
    static Core::LinAlg::Matrix<9, 1> iFeTV(true);
    Core::LinAlg::Voigt::matrix_3x3_to_9x1(iFeTM, iFeTV);
    static Core::LinAlg::Matrix<1, 9> d_iJe_dFV(true);
    d_iJe_dFV.multiply_tn(1.0, iFeTV, d_Fe_dF, 0.0);
    d_cauchyndir_dF->update_t(-cauchy_n_dir, d_iJe_dFV, 1.0);

    // second part is term arising from \partial b_el * n * v / \partial F
    static Core::LinAlg::Matrix<3, 3> FeMiFinTM(true);
    FeMiFinTM.multiply_nt(1.0, FeM, iFinM, 0.0);
    static Core::LinAlg::Matrix<3, 1> tempvec(true);
    tempvec.multiply_tn(1.0, FeMiFinTM, n, 0.0);
    static Core::LinAlg::Matrix<3, 3> d_bednddir_dF(true);
    d_bednddir_dF.multiply_nt(1.0, dir, tempvec, 0.0);
    // now reuse tempvec
    tempvec.multiply_tn(1.0, FeMiFinTM, dir, 0.0);
    d_bednddir_dF.multiply_nt(1.0, n, tempvec, 1.0);
    static Core::LinAlg::Matrix<9, 1> d_bednddir_dFV(true);
    Core::LinAlg::Voigt::matrix_3x3_to_9x1(d_bednddir_dF, d_bednddir_dFV);
    d_cauchyndir_dF->update(prefac * dPI(0), d_bednddir_dFV, 1.0);

    // third part is term arising from \partial b_el^{-1} * n * v / \partial F
    static Core::LinAlg::Matrix<3, 3> iFM(true);
    iFM.invert(defgrd);
    static Core::LinAlg::Matrix<3, 1> tempvec2(true);
    tempvec.multiply(1.0, ibeM, dir, 0.0);
    tempvec2.multiply(1.0, iFM, n, 0.0);
    static Core::LinAlg::Matrix<3, 3> d_ibednddir_dFM(true);
    d_ibednddir_dFM.multiply_nt(1.0, tempvec, tempvec2, 0.0);
    // now reuse tempvecs
    tempvec.multiply(1.0, ibeM, n, 0.0);
    tempvec2.multiply(1.0, iFM, dir, 0.0);
    d_ibednddir_dFM.multiply_nt(1.0, tempvec, tempvec2, 1.0);
    d_ibednddir_dFM.scale(-1.0);
    static Core::LinAlg::Matrix<9, 1> d_ibednddir_dFV(true);
    Core::LinAlg::Voigt::matrix_3x3_to_9x1(d_ibednddir_dFM, d_ibednddir_dFV);
    d_cauchyndir_dF->update(-prefac * prinv(2) * dPI(1), d_ibednddir_dFV, 1.0);
  }
}

/*--------------------------------------------------------------------*
 *--------------------------------------------------------------------*/
void Mat::MyExpansionTest_ElastHyper::evaluate_invariant_derivatives(
    const Core::LinAlg::Matrix<3, 1>& prinv, const int gp, const int eleGID,
    Core::LinAlg::Matrix<3, 1>& dPI, Core::LinAlg::Matrix<6, 1>& ddPII) const
{
  // clear variables
  dPI.clear();
  ddPII.clear();

  // loop over map of associated potential summands
  // derivatives of strain energy function w.r.t. principal invariants
  for (const auto& p : potsumeliso_)
  {
    p->add_derivatives_principal(dPI, ddPII, prinv, gp, eleGID);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::MyExpansionTest_ElastHyper::evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<6, 1>* glstrain, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<6, 1>* stress, Core::LinAlg::Matrix<6, 6>* cmat, const int gp,
    const int eleGID)
{
  // blank resulting quantities
  // ... even if it is an implicit law that cmat is zero upon input
  stress->clear();
  cmat->clear();
  // std::cout << "evaluate " << gp << std::endl;
  double temperature = 1.0;
  double temperature_0 = 293.15;
  double rho_s = 0.0;
  double rho_ss = params_->rho_ss_;
  double rho_s0 = params_->rho_0_;

  double alpha_temp = 13.0e-6;  // this is the value for steel just to test
  double alpha_reaction = params_->exp_rate_;

  // TODO: do not read from parameter list!
  if (params.isParameter("scalar"))
  {
    std::shared_ptr<std::vector<double>> scalars =
        params.get<std::shared_ptr<std::vector<double>>>("scalar");
    try
    {
      temperature = scalars->at(2);
      rho_s = scalars->at(1);
    }
    catch (const std::out_of_range& e)
    {
      FOUR_C_THROW("Not enough scalars to use my expansion test material");
    }
  }

  rho_average_ = rho_s;  // only works for homogeneous exp. only used for contact TODO change
  double scaling_factor = (1.0 + alpha_temp * (temperature - temperature_0)) *
                          (1.0 + alpha_reaction * (rho_s - rho_s0) / (rho_ss - rho_s0));

  //  set-up inverse expansion tensor
  Core::LinAlg::Matrix<3, 3> iFexpM(true);
  for (int i = 0; i < 3; ++i) iFexpM(i, i) = 1.0;
  iFexpM.scale(1 / scaling_factor);

  // some variables
  static Core::LinAlg::Matrix<6, 1> iCinv(true);
  static Core::LinAlg::Matrix<6, 1> iCinCiCinv(true);  // C_in^-1 C C_in^-1 (voight)
  static Core::LinAlg::Matrix<6, 1> iCv(true);         // C^-1 (voight)
  static Core::LinAlg::Matrix<3, 1> prinv(true);
  static Core::LinAlg::Matrix<3, 3> iCinCM(true);  // C_in^-1 (matrix)
  static Core::LinAlg::Matrix<3, 3> iFinCeM(true);
  static Core::LinAlg::Matrix<9, 1> CiFin9x1(true);
  Core::LinAlg::Matrix<9, 1> CiFinCe9x1(true);
  Core::LinAlg::Matrix<9, 1> CiFiniCe9x1(true);

  EvaluateKinQuantElast(defgrd, iFexpM, iCinv, iCinCiCinv, iCv, iCinCM, iFinCeM, CiFin9x1,
      CiFinCe9x1, CiFiniCe9x1, prinv, gp);

  Core::LinAlg::Matrix<3, 1> dPIe(true);
  Core::LinAlg::Matrix<6, 1> ddPIIe(true);

  static Core::LinAlg::Matrix<3, 1> dPI_i(true);
  static Core::LinAlg::Matrix<6, 1> ddPII_i(true);

  // loop map of associated potential summands
  // derivatives of strain energy function w.r.t. principal invariants
  for (const auto& p : potsumeliso_)
  {
    dPI_i.clear();
    ddPII_i.clear();
    p->add_derivatives_principal(dPI_i, ddPII_i, prinv, gp, eleGID);
    dPIe.update(1.0, dPI_i, 1.0);
    ddPIIe.update(1.0, ddPII_i, 1.0);
  }

  // derivatives of decoupled (volumetric or isochoric) materials weighted with their mass fraction
  // in the constraint mixture
  static Core::LinAlg::Matrix<3, 1> modinv(true);
  Mat::invariants_modified(modinv, prinv);
  Core::LinAlg::Matrix<3, 1> dPmodI(true);
  Core::LinAlg::Matrix<6, 1> ddPmodII(true);
  for (const auto& p : potsumeliso_)
  {
    dPI_i.clear();
    ddPII_i.clear();
    p->add_derivatives_modified(dPI_i, ddPII_i, modinv, gp, eleGID);
    dPmodI.update(1.0, dPI_i, 1.0);
    ddPmodII.update(1.0, ddPII_i, 1.0);
  }

  // convert decoupled derivatives to principal derivatives
  Mat::convert_mod_to_princ(prinv, dPmodI, ddPmodII, dPIe, ddPIIe);

  // 2nd Piola Kirchhoff stress factors (according to Holzapfel-Nonlinear Solid Mechanics p. 216)
  static Core::LinAlg::Matrix<3, 1> gamma(true);
  // constitutive tensor factors (according to Holzapfel-Nonlinear Solid Mechanics p. 261)
  static Core::LinAlg::Matrix<8, 1> delta(true);

  // compose coefficients
  calculate_gamma_delta(gamma, delta, prinv, dPIe, ddPIIe);

  EvaluateIsotropicPrincElast(*stress, *cmat, iCinv, iCinCiCinv, iCv, gamma, delta);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::MyExpansionTest_ElastHyper::EvaluateKinQuantElast(
    Core::LinAlg::Matrix<3, 3> const* const defgrd, Core::LinAlg::Matrix<3, 3> const& iFinM,
    Core::LinAlg::Matrix<6, 1>& iCinv, Core::LinAlg::Matrix<6, 1>& iCinCiCinv,
    Core::LinAlg::Matrix<6, 1>& iCv, Core::LinAlg::Matrix<3, 3>& iCinCM,
    Core::LinAlg::Matrix<3, 3>& iFinCeM, Core::LinAlg::Matrix<9, 1>& CiFin9x1,
    Core::LinAlg::Matrix<9, 1>& CiFinCe9x1, Core::LinAlg::Matrix<9, 1>& CiFiniCe9x1,
    Core::LinAlg::Matrix<3, 1>& prinv, const int gp)
{
  // inverse inelastic right Cauchy-Green
  static Core::LinAlg::Matrix<3, 3> iCinM(true);
  iCinM.multiply_nt(1.0, iFinM, iFinM, 0.0);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(iCinM, iCinv);

  // inverse right Cauchy-Green
  static Core::LinAlg::Matrix<3, 3> iCM(true);
  static Core::LinAlg::Matrix<3, 3> CM(true);
  CM.multiply_tn(1.0, *defgrd, *defgrd, 0.0);
  iCM.invert(CM);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(iCM, iCv);

  // C_{in}^{-1} * C * C_{in}^{-1}
  static Core::LinAlg::Matrix<3, 3> tmp(true);
  static Core::LinAlg::Matrix<3, 3> iCinCiCinM;
  tmp.multiply_nn(1.0, iCinM, CM, 0.0);
  iCinCiCinM.multiply_nn(1.0, tmp, iCinM, 0.0);
  Core::LinAlg::Voigt::Stresses::matrix_to_vector(iCinCiCinM, iCinCiCinv);

  // elastic right Cauchy-Green in strain-like Voigt notation.
  tmp.multiply_nn(1.0, *defgrd, iFinM, 0.0);
  static Core::LinAlg::Matrix<3, 3> CeM(true);
  CeM.multiply_tn(1.0, tmp, tmp, 0.0);
  static Core::LinAlg::Matrix<6, 1> Ce_strain(true);
  Core::LinAlg::Voigt::Strains::matrix_to_vector(CeM, Ce_strain);

  // principal invariants of elastic right Cauchy-Green strain
  Core::LinAlg::Voigt::Strains::invariants_principal(prinv, Ce_strain);

  // C_{in}^{-1} * C
  iCinCM.multiply_nn(1.0, iCinM, CM, 0.0);

  // F_{in}^{-1} * C_e
  iFinCeM.multiply_nn(1.0, iFinM, CeM, 0.0);

  // C * F_{in}^{-1}
  static Core::LinAlg::Matrix<3, 3> CiFinM(true);
  CiFinM.multiply_nn(1.0, CM, iFinM, 0.0);
  Core::LinAlg::Voigt::matrix_3x3_to_9x1(CiFinM, CiFin9x1);

  // C * F_{in}^{-1} * C_e
  static Core::LinAlg::Matrix<3, 3> CiFinCeM(true);
  tmp.multiply_nn(1.0, CM, iFinM, 0.0);
  CiFinCeM.multiply_nn(1.0, tmp, CeM, 0.0);
  Core::LinAlg::Voigt::matrix_3x3_to_9x1(CiFinCeM, CiFinCe9x1);

  // C * F_{in}^{-1} * C_e^{-1}
  static Core::LinAlg::Matrix<3, 3> CiFiniCeM(true);
  static Core::LinAlg::Matrix<3, 3> iCeM(true);
  iCeM.invert(CeM);
  tmp.multiply_nn(1.0, CM, iFinM, 0.0);
  CiFiniCeM.multiply_nn(1.0, tmp, iCeM, 0.0);
  Core::LinAlg::Voigt::matrix_3x3_to_9x1(CiFiniCeM, CiFiniCe9x1);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::MyExpansionTest_ElastHyper::EvaluateIsotropicPrincElast(
    Core::LinAlg::Matrix<6, 1>& stressisoprinc, Core::LinAlg::Matrix<6, 6>& cmatisoprinc,
    Core::LinAlg::Matrix<6, 1> const& iCinv, Core::LinAlg::Matrix<6, 1> const& iCinCiCinv,
    Core::LinAlg::Matrix<6, 1> const& iCv, Core::LinAlg::Matrix<3, 1> const& gamma,
    Core::LinAlg::Matrix<8, 1> const& delta) const
{
  // 2nd Piola Kirchhoff stresses
  stressisoprinc.update(gamma(0), iCinv, 1.0);
  stressisoprinc.update(gamma(1), iCinCiCinv, 1.0);
  stressisoprinc.update(gamma(2), iCv, 1.0);

  // constitutive tensor
  cmatisoprinc.multiply_nt(delta(0), iCinv, iCinv, 1.);
  cmatisoprinc.multiply_nt(delta(1), iCinCiCinv, iCinv, 1.);
  cmatisoprinc.multiply_nt(delta(1), iCinv, iCinCiCinv, 1.);
  cmatisoprinc.multiply_nt(delta(2), iCinv, iCv, 1.);
  cmatisoprinc.multiply_nt(delta(2), iCv, iCinv, 1.);
  cmatisoprinc.multiply_nt(delta(3), iCinCiCinv, iCinCiCinv, 1.);
  cmatisoprinc.multiply_nt(delta(4), iCinCiCinv, iCv, 1.);
  cmatisoprinc.multiply_nt(delta(4), iCv, iCinCiCinv, 1.);
  cmatisoprinc.multiply_nt(delta(5), iCv, iCv, 1.);
  Core::LinAlg::Tensor::add_holzapfel_product(cmatisoprinc, iCv, delta(6));
  Core::LinAlg::Tensor::add_holzapfel_product(cmatisoprinc, iCinv, delta(7));
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::MyExpansionTest_ElastHyper::vis_names(std::map<std::string, int>& names) const
{
  // loop map of associated potential summands
  // 3D elastin matrix
  for (auto& p : potsumeliso_) p->vis_names(names);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool Mat::MyExpansionTest_ElastHyper::vis_data(
    const std::string& name, std::vector<double>& data, int numgp, int eleID) const
{
  int return_val = 0;

  // loop map of associated potential summands
  // 3D elastin matrix
  for (auto& p : potsumeliso_) return_val += p->vis_data(name, data, numgp, eleID);

  return (bool)return_val;
}

FOUR_C_NAMESPACE_CLOSE
