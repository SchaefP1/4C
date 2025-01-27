// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_plasticdruckerprager_exp_test.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_comm_parobject.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_structure.hpp"
#include "4C_linalg_FADmatrix_utils.hpp"
#include "4C_linalg_fixedsizematrix.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_mat_material_factory.hpp"
#include "4C_mat_par_bundle.hpp"
#include "4C_mat_stvenantkirchhoff.hpp"
#include "4C_material_base.hpp"
#include "4C_material_parameter_base.hpp"
#include "4C_utils_local_newton.hpp"

FOUR_C_NAMESPACE_OPEN


Mat::PAR::PlasticDruckerPragerExpTest::PlasticDruckerPragerExpTest(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : Parameter(matdata),
      youngs_(matdata.parameters.get<double>("YOUNG")),
      poissonratio_(matdata.parameters.get<double>("NUE")),
      density_(matdata.parameters.get<double>("DENS")),
      isohard_(matdata.parameters.get<double>("ISOHARD")),
      abstol_(matdata.parameters.get<double>("TOL")),
      cohesion_(matdata.parameters.get<double>("C")),
      eta_(matdata.parameters.get<double>("ETA")),
      xi_(matdata.parameters.get<double>("XI")),
      etabar_(matdata.parameters.get<double>("ETABAR")),
      itermax_(matdata.parameters.get<int>("MAXITER"))
{
}

std::shared_ptr<Core::Mat::Material> Mat::PAR::PlasticDruckerPragerExpTest::create_material()
{
  return std::make_shared<Mat::PlasticDruckerPragerExpTest>(this);
}

Mat::PlasticDruckerPragerExpTestType Mat::PlasticDruckerPragerExpTestType::instance_;

Core::Communication::ParObject* Mat::PlasticDruckerPragerExpTestType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Mat::PlasticDruckerPragerExpTest* plastic = new Mat::PlasticDruckerPragerExpTest();
  plastic->unpack(buffer);
  return plastic;
}

Mat::PlasticDruckerPragerExpTest::PlasticDruckerPragerExpTest() : params_(nullptr), rho_average_(0)
{
}

Mat::PlasticDruckerPragerExpTest::PlasticDruckerPragerExpTest(
    Mat::PAR::PlasticDruckerPragerExpTest* params)
    : params_(params), rho_average_(0)
{
}

void Mat::PlasticDruckerPragerExpTest::pack(Core::Communication::PackBuffer& data) const
{
  int type = unique_par_object_id();
  add_to_pack(data, type);
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();
  add_to_pack(data, matid);
  int histsize = Initialized() ? strainpllast_.size() : 0;
  add_to_pack(data, histsize);
  for (int var = 0; var < histsize; ++var)
  {
    add_to_pack(data, strainpllast_.at(var));
    add_to_pack(data, strainbarpllast_.at(var));
  }
  add_to_pack(data, rho_average_);
}

void Mat::PlasticDruckerPragerExpTest::unpack(Core::Communication::UnpackBuffer& buffer)
{
  isinit_ = true;


  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

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
        params_ = static_cast<Mat::PAR::PlasticDruckerPragerExpTest*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

    int histsize;
    extract_from_pack(buffer, histsize);

    if (histsize == 0) isinit_ = false;

    strainpllast_ = std::vector<Core::LinAlg::Matrix<NUM_STRESS_3D, 1>>();
    strainplcurr_ = std::vector<Core::LinAlg::Matrix<NUM_STRESS_3D, 1>>();
    strainbarpllast_ = std::vector<double>();
    strainbarplcurr_ = std::vector<double>();
    for (int var = 0; var < histsize; ++var)
    {
      Core::LinAlg::Matrix<NUM_STRESS_3D, 1> tmp_vect(true);
      double tmp_scalar = 0.0;

      extract_from_pack(buffer, tmp_vect);
      strainpllast_.push_back(tmp_vect);

      extract_from_pack(buffer, tmp_scalar);
      strainbarpllast_.push_back(tmp_scalar);

      strainplcurr_.push_back(tmp_vect);
      strainbarplcurr_.push_back(tmp_scalar);
    }
  }
}

void Mat::PlasticDruckerPragerExpTest::evaluate_cauchy_n_dir_and_derivatives(
    const Core::LinAlg::Matrix<3, 3>& defgrd, const Core::LinAlg::Matrix<3, 1>& n,
    const Core::LinAlg::Matrix<3, 1>& dir, double& cauchy_n_dir,
    Core::LinAlg::Matrix<3, 1>* d_cauchyndir_dn, Core::LinAlg::Matrix<3, 1>* d_cauchyndir_ddir,
    Core::LinAlg::Matrix<9, 1>* d_cauchyndir_dF, Core::LinAlg::Matrix<9, 9>* d2_cauchyndir_dF2,
    Core::LinAlg::Matrix<9, 3>* d2_cauchyndir_dF_dn,
    Core::LinAlg::Matrix<9, 3>* d2_cauchyndir_dF_ddir, int gp, int eleGID,
    const double* concentration, const double* temp, double* d_cauchyndir_dT,
    Core::LinAlg::Matrix<9, 1>* d2_cauchyndir_dF_dT)
{
  std::cout << "defgrd:" << std::endl;
  defgrd.print(std::cout);
  //  // reset sigma contracted with n and dir
  //  cauchy_n_dir = 0.0;
  //  double temperature = 293.15;
  //  double temperature_0 = 293.15;
  //  double rho_s = 0.0;
  //  double rho_ss = 8520.0;
  //  double rho_s0 = 8400.0;
  //  double scaling_factor = 1.0;
  //
  //  double alpha_temp = 13.0e-6;  // this is the value for steel just to test
  //  double alpha_reaction = 0.065;
  //  scaling_factor = (1.0 + alpha_temp * (temperature - temperature_0)) *
  //                   (1.0 + alpha_reaction * (rho_average_ - rho_s0) / (rho_ss - rho_s0));
  //
  //  //std::cout << "evaluate_cauchy_n_dir_and_derivatives " << rho_average_ << " " <<
  //  scaling_factor << std::endl; static Core::LinAlg::Matrix<6, 1> idV(true); for (int i = 0; i <
  //  3; ++i) idV(i) = 1.0; static Core::LinAlg::Matrix<3, 3> idM(true); for (int i = 0; i < 3; ++i)
  //  idM(i, i) = 1.0; static Core::LinAlg::Matrix<3, 3> iFinM(true); for (int i = 0; i < 3; ++i)
  //  iFinM(i, i) = 1.0; iFinM.scale(1 / scaling_factor);
  //  // TODO I need my inelasatic stuff here
  //  // inelastic_->evaluate_inverse_inelastic_def_grad(&defgrd, iFinM);
  //  static Core::LinAlg::Matrix<3, 3> FeM(true);
  //  FeM.multiply_nn(1.0, defgrd, iFinM, 0.0);
  //
  //  // get elastic left cauchy-green tensor and corresponding principal invariants
  //  static Core::LinAlg::Matrix<3, 3> beM(true);
  //  beM.multiply_nt(1.0, FeM, FeM, 0.0);
  //  static Core::LinAlg::Matrix<6, 1> beV_strain(true);
  //  Core::LinAlg::Voigt::Strains::matrix_to_vector(beM, beV_strain);
  //  static Core::LinAlg::Matrix<3, 1> prinv(true);
  //  Core::LinAlg::Voigt::Strains::invariants_principal(prinv, beV_strain);
  //  static Core::LinAlg::Matrix<6, 1> beV_stress(true);
  //  Core::LinAlg::Voigt::Stresses::matrix_to_vector(beM, beV_stress);
  //
  //  static Core::LinAlg::Matrix<3, 1> beMdn(true);
  //  beMdn.multiply(1.0, beM, n, 0.0);
  //  const double beMdnddir = beMdn.dot(dir);
  //  static Core::LinAlg::Matrix<3, 1> beMddir(true);
  //  beMddir.multiply(1.0, beM, dir, 0.0);
  //
  //  static Core::LinAlg::Matrix<3, 3> ibeM(true);
  //  ibeM.invert(beM);
  //  static Core::LinAlg::Matrix<6, 1> ibeV_stress(true);
  //  Core::LinAlg::Voigt::Stresses::matrix_to_vector(ibeM, ibeV_stress);
  //  static Core::LinAlg::Matrix<3, 1> ibeMdn(true);
  //  ibeMdn.multiply(1.0, ibeM, n, 0.0);
  //  const double ibeMdnddir = ibeMdn.dot(dir);
  //  static Core::LinAlg::Matrix<3, 1> ibeMddir(true);
  //  ibeMddir.multiply(1.0, ibeM, dir, 0.0);
  //
  //  // derivatives of principle invariants of elastic left cauchy-green tensor
  //  static Core::LinAlg::Matrix<3, 1> dPI(true);
  //  static Core::LinAlg::Matrix<6, 1> ddPII(true);
  //  evaluate_invariant_derivatives(prinv, gp, eleGID, dPI, ddPII);
  //
  //  const double detFe = FeM.determinant();
  //  const double nddir = n.dot(dir);
  //  const double prefac = 2.0 / detFe;
  //
  //  // calculate \mat{\sigma} \cdot \vec{n} \cdot \vec{v}
  //  cauchy_n_dir = prefac * (prinv(1) * dPI(1) * nddir + prinv(2) * dPI(2) * nddir +
  //                              dPI(0) * beMdnddir - prinv(2) * dPI(1) * ibeMdnddir);
  //
  //  if (d_cauchyndir_dn)
  //  {
  //    d_cauchyndir_dn->update(prinv(1) * dPI(1) + prinv(2) * dPI(2), dir, 0.0);
  //    d_cauchyndir_dn->update(dPI(0), beMddir, 1.0);
  //    d_cauchyndir_dn->update(-prinv(2) * dPI(1), ibeMddir, 1.0);
  //    d_cauchyndir_dn->scale(prefac);
  //  }
  //
  //  if (d_cauchyndir_ddir)
  //  {
  //    d_cauchyndir_ddir->update(prinv(1) * dPI(1) + prinv(2) * dPI(2), n, 0.0);
  //    d_cauchyndir_ddir->update(dPI(0), beMdn, 1.0);
  //    d_cauchyndir_ddir->update(-prinv(2) * dPI(1), ibeMdn, 1.0);
  //    d_cauchyndir_ddir->scale(prefac);
  //  }
  //
  //  if (d_cauchyndir_dF)
  //  {
  //    static Core::LinAlg::Matrix<6, 1> d_I1_be(true);
  //    d_I1_be = idV;
  //    static Core::LinAlg::Matrix<6, 1> d_I2_be(true);
  //    d_I2_be.update(prinv(0), idV, -1.0, beV_stress);
  //    static Core::LinAlg::Matrix<6, 1> d_I3_be(true);
  //    d_I3_be.update(prinv(2), ibeV_stress, 0.0);
  //
  //    // calculation of \partial b_{el} / \partial F (elastic left cauchy-green w.r.t. deformation
  //    // gradient)
  //    static Core::LinAlg::Matrix<6, 9> d_be_dFe(true);
  //    d_be_dFe.clear();
  //    add_right_non_symmetric_holzapfel_product_strain_like(d_be_dFe, idM, FeM, 1.0);
  //    static Core::LinAlg::Matrix<9, 9> d_Fe_dF(true);
  //    d_Fe_dF.clear();
  //    add_non_symmetric_product(1.0, idM, iFinM, d_Fe_dF);
  //    static Core::LinAlg::Matrix<6, 9> d_be_dF(true);
  //    d_be_dF.multiply(1.0, d_be_dFe, d_Fe_dF, 0.0);
  //
  //    // calculation of \partial I_i / \partial F (Invariants of b_{el} w.r.t. deformation
  //    gradient) static Core::LinAlg::Matrix<9, 1> d_I1_dF(true); static Core::LinAlg::Matrix<9, 1>
  //    d_I2_dF(true); static Core::LinAlg::Matrix<9, 1> d_I3_dF(true); d_I1_dF.multiply_tn(1.0,
  //    d_be_dF, d_I1_be, 0.0); d_I2_dF.multiply_tn(1.0, d_be_dF, d_I2_be, 0.0);
  //    d_I3_dF.multiply_tn(1.0, d_be_dF, d_I3_be, 0.0);
  //
  //    // add d_cauchyndir_dI1 \odot d_I1_dF and clear static matrix
  //    d_cauchyndir_dF->update(prefac * (prinv(1) * ddPII(5) * nddir + prinv(2) * ddPII(4) * nddir
  //    +
  //                                         ddPII(0) * beMdnddir - prinv(2) * ddPII(5) *
  //                                         ibeMdnddir),
  //        d_I1_dF, 0.0);
  //    // add d_cauchyndir_dI2 \odot d_I2_dF
  //    d_cauchyndir_dF->update(
  //        prefac * (dPI(1) * nddir + prinv(1) * ddPII(1) * nddir + prinv(2) * ddPII(3) * nddir +
  //                     ddPII(5) * beMdnddir - prinv(2) * ddPII(1) * ibeMdnddir),
  //        d_I2_dF, 1.0);
  //    // add d_cauchyndir_dI3 \odot d_I3_dF
  //    d_cauchyndir_dF->update(
  //        prefac * (prinv(1) * ddPII(3) * nddir + dPI(2) * nddir + prinv(2) * ddPII(2) * nddir +
  //                     ddPII(4) * beMdnddir - dPI(1) * ibeMdnddir - prinv(2) * ddPII(3) *
  //                     ibeMdnddir),
  //        d_I3_dF, 1.0);
  //
  //    // next three updates add partial derivative of snt w.r.t. the deformation gradient F for
  //    // constant invariants first part is term arising from \partial Je^{-1} / \partial F
  //    static Core::LinAlg::Matrix<3, 3> iFeM(true);
  //    static Core::LinAlg::Matrix<3, 3> iFeTM(true);
  //    iFeM.invert(FeM);
  //    iFeTM.update_t(1.0, iFeM, 0.0);
  //    static Core::LinAlg::Matrix<9, 1> iFeTV(true);
  //    Core::LinAlg::Voigt::matrix_3x3_to_9x1(iFeTM, iFeTV);
  //    static Core::LinAlg::Matrix<1, 9> d_iJe_dFV(true);
  //    d_iJe_dFV.multiply_tn(1.0, iFeTV, d_Fe_dF, 0.0);
  //    d_cauchyndir_dF->update_t(-cauchy_n_dir, d_iJe_dFV, 1.0);
  //
  //    // second part is term arising from \partial b_el * n * v / \partial F
  //    static Core::LinAlg::Matrix<3, 3> FeMiFinTM(true);
  //    FeMiFinTM.multiply_nt(1.0, FeM, iFinM, 0.0);
  //    static Core::LinAlg::Matrix<3, 1> tempvec(true);
  //    tempvec.multiply_tn(1.0, FeMiFinTM, n, 0.0);
  //    static Core::LinAlg::Matrix<3, 3> d_bednddir_dF(true);
  //    d_bednddir_dF.multiply_nt(1.0, dir, tempvec, 0.0);
  //    // now reuse tempvec
  //    tempvec.multiply_tn(1.0, FeMiFinTM, dir, 0.0);
  //    d_bednddir_dF.multiply_nt(1.0, n, tempvec, 1.0);
  //    static Core::LinAlg::Matrix<9, 1> d_bednddir_dFV(true);
  //    Core::LinAlg::Voigt::matrix_3x3_to_9x1(d_bednddir_dF, d_bednddir_dFV);
  //    d_cauchyndir_dF->update(prefac * dPI(0), d_bednddir_dFV, 1.0);
  //
  //    // third part is term arising from \partial b_el^{-1} * n * v / \partial F
  //    static Core::LinAlg::Matrix<3, 3> iFM(true);
  //    iFM.invert(defgrd);
  //    static Core::LinAlg::Matrix<3, 1> tempvec2(true);
  //    tempvec.multiply(1.0, ibeM, dir, 0.0);
  //    tempvec2.multiply(1.0, iFM, n, 0.0);
  //    static Core::LinAlg::Matrix<3, 3> d_ibednddir_dFM(true);
  //    d_ibednddir_dFM.multiply_nt(1.0, tempvec, tempvec2, 0.0);
  //    // now reuse tempvecs
  //    tempvec.multiply(1.0, ibeM, n, 0.0);
  //    tempvec2.multiply(1.0, iFM, dir, 0.0);
  //    d_ibednddir_dFM.multiply_nt(1.0, tempvec, tempvec2, 1.0);
  //    d_ibednddir_dFM.scale(-1.0);
  //    static Core::LinAlg::Matrix<9, 1> d_ibednddir_dFV(true);
  //    Core::LinAlg::Voigt::matrix_3x3_to_9x1(d_ibednddir_dFM, d_ibednddir_dFV);
  //    d_cauchyndir_dF->update(-prefac * prinv(2) * dPI(1), d_ibednddir_dFV, 1.0);
  //  }
}

void Mat::PlasticDruckerPragerExpTest::setup(
    int numgp, const Core::IO::InputParameterContainer& container)
{
  strainpllast_.resize(numgp);
  strainplcurr_.resize(numgp);

  strainbarpllast_.resize(numgp);
  strainbarplcurr_.resize(numgp);

  isinit_ = true;
}

void Mat::PlasticDruckerPragerExpTest::update()
{
  strainpllast_ = strainplcurr_;
  strainbarpllast_ = strainbarplcurr_;

  std::for_each(strainplcurr_.begin(), strainplcurr_.end(), [](auto& item) { item.clear(); });
  std::fill(strainbarplcurr_.begin(), strainbarplcurr_.end(), 0.0);
}

void Mat::PlasticDruckerPragerExpTest::setup_cmat(
    Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat)
{
  double young = params_->youngs_;

  double nu = params_->poissonratio_;

  Mat::StVenantKirchhoff::fill_cmat(cmat, young, nu);
}

void Mat::PlasticDruckerPragerExpTest::setup_cmat_elasto_plastic_cone(
    Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat, double Dgamma, double G, double Kappa,
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& devstrain, double xi, double Hiso, double eta,
    double etabar)
{
  cmat.clear();

  Core::LinAlg::Matrix<NUM_STRESS_3D, 1> id2(true);
  Core::LinAlg::Voigt::identity_matrix(id2);
  Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D> id4sharp(true);
  for (int i = 0; i < 3; i++) id4sharp(i, i) = 1.0;
  for (int i = 3; i < 6; i++) id4sharp(i, i) = 0.5;

  const double normdevstrain =
      sqrt(devstrain(0) * devstrain(0) + devstrain(1) * devstrain(1) + devstrain(2) * devstrain(2) +
           2 * (devstrain(3) * devstrain(3) + devstrain(4) * devstrain(4) +
                   devstrain(5) * devstrain(5)));
  const double epfac = 2 * G * (1 - (Dgamma / sqrt(2) / normdevstrain));

  cmat.update(epfac, id4sharp, 1.0);

  double epfac1 = 0.0;
  double epfac2 = 0.0;
  double epfac3 = 0.0;
  double epfac4 = 0.0;
  epfac1 = epfac / (-3.0);
  cmat.multiply_nt(epfac1, id2, id2, 1.0);

  double A = 0.0;
  A = 1 / (G + Kappa * etabar * eta + xi * xi * Hiso);
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1> D(true);
  D.update(1 / normdevstrain, devstrain);
  epfac2 = 2 * G * (Dgamma / (sqrt(2) * normdevstrain) - G * A);

  for (int k = 0; k < 6; ++k)
  {
    for (int i = 0; i < 6; ++i)
    {
      cmat(i, k) += epfac2 * D(i) * D(k);
    }
  }

  epfac3 = -sqrt(2) * G * A * Kappa;

  for (int k = 0; k < 6; ++k)
  {
    for (int i = 0; i < 6; ++i)
    {
      cmat(k, i) += epfac3 * (eta * D(k) * id2(i) + etabar * id2(k) * D(i));
    }
  }

  epfac4 = Kappa * (1 - Kappa * eta * etabar * A);

  for (int k = 0; k < 6; ++k)
  {
    for (int i = 0; i < 6; ++i)
    {
      cmat(i, k) += epfac4 * id2(k) * id2(i);
    }
  }
}

void Mat::PlasticDruckerPragerExpTest::setup_cmat_elasto_plastic_apex(
    Core::LinAlg::Matrix<NUM_STRESS_3D, NUM_STRESS_3D>& cmat, double Kappa,
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1>& devstrain, double xi, double Hiso, double eta,
    double etabar)
{
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1> id2(true);
  for (int i = 0; i < 3; i++) id2(i) = 1.0;
  double epfac = 0.0;
  epfac = Kappa * (1 - Kappa / (Kappa + xi / eta * xi / etabar * Hiso));
  cmat.clear();
  cmat.multiply_nt(epfac, id2, id2, 0.0);
}

template <typename ScalarT>
void Mat::PlasticDruckerPragerExpTest::EvaluateFAD(const Core::LinAlg::Matrix<3, 3>* defgrd,
    const Core::LinAlg::Matrix<6, 1, ScalarT>* linstrain, Teuchos::ParameterList& params,
    Core::LinAlg::Matrix<6, 1, ScalarT>* stress, Core::LinAlg::Matrix<6, 6>* cmat, const int gp,
    const int eleGID)
{
  // params.print();

  double temperature = 293.15;
  double temperature_0 = 293.15;
  double rho_s = 0.0;
  double rho_ss = 8520.0;
  double rho_s0 = 8400.0;
  double scaling_factor = 1.0;

  double alpha_temp = 13.0e-6;  // this is the value for steel just to test
  double alpha_reaction = 0.065;

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
  // TODO: do not read from parameter list!
  if (params.isParameter("gp_conc"))
  {
    std::shared_ptr<std::vector<std::vector<double>>> scalars =
        params.get<std::shared_ptr<std::vector<std::vector<double>>>>("gp_conc");
    try
    {
      // std::cout << "Test scalars" << std::endl;
      std::vector<double> scalars_at_gp = scalars->at(gp);
      temperature = scalars_at_gp.at(2);
      rho_s = scalars_at_gp.at(1);
    }
    catch (const std::out_of_range& e)
    {
      FOUR_C_THROW("Not enough scalars to use my expansion test material");
    }
  }
  rho_average_ = rho_s;  // only works for homogeneous exp. TODO change
  // std::cout << "temperature:" << temperature  << "rho_s:" << rho_s << std::endl;
  scaling_factor = alpha_temp * (temperature - temperature_0) +
                   alpha_reaction * (rho_s - rho_s0) / (rho_ss - rho_s0);

  Core::LinAlg::Matrix<Mat::NUM_STRESS_3D, 1> plstrain(true);

  ScalarT young = params_->youngs_;
  ScalarT nu = params_->poissonratio_;
  ScalarT Hiso = params_->isohard_;
  ScalarT cohesion = params_->cohesion_;
  ScalarT eta = params_->eta_;
  ScalarT xi = params_->xi_;
  ScalarT etabar = params_->etabar_;
  const int itermax = params_->itermax_;
  ScalarT G = 0.0;
  G = young / (2.0 * (1.0 + nu));
  ScalarT kappa = 0.0;
  kappa = young / (3.0 * (1.0 - 2.0 * nu));
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1> id2(true);
  for (int i = 0; i < 3; i++) id2(i) = 1.0;
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT> strain(*linstrain);

  // subtract inelastic strain so that strain is only elastic and plastic strain
  for (int i = 0; i < 3; ++i) strain(i) -= scaling_factor;

  Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT> strain_p(false);
  for (int i = 0; i < 6; i++) strain_p(i, 0) = strainpllast_.at(gp)(i, 0);
  ScalarT strainbar_p = 0.0;
  strainbar_p = (strainbarpllast_.at(gp));
  if (strainbarpllast_.at(gp) < 0.0)
    FOUR_C_THROW("accumulated plastic strain has to be equal to or greater than zero!");

  for (int i = 3; i < 6; ++i) strain(i) /= 2.0;
  for (int i = 3; i < 6; ++i) strain_p(i) /= 2.0;

  Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT> strain_e(true);
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT> trialstrain_e(false);
  trialstrain_e.update(1.0, strain, (-1.0), strain_p);
  ScalarT tracestrain = trialstrain_e(0) + trialstrain_e(1) + trialstrain_e(2);
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT> volumetricstrain(false);
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT> id2Scalar(true);
  for (int i = 0; i < NUM_STRESS_3D; ++i) id2Scalar(i) = static_cast<ScalarT>(id2(i));
  volumetricstrain.update((tracestrain / 3.0), id2Scalar);
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT> devstrain(false);
  devstrain.update(1.0, trialstrain_e, (-1.0), volumetricstrain);

  //  if (gp == 1 && scaling_factor > 0.0001)
  //  {
  //    std::cout << "scaling_factor:" << scaling_factor << std::endl;
  //    std::cout << "strain:" << strain << std::endl;
  //    std::cout << "trialstrain_e:" << trialstrain_e << std::endl;
  //    std::cout << "strain_p:" << strain_p << std::endl;
  //  }

  ScalarT p = kappa * tracestrain;
  ScalarT p_trial = p;
  Core::LinAlg::Matrix<NUM_STRESS_3D, 1, ScalarT> devstress(false);
  devstress.update((2.0 * G), devstrain);

  ScalarT J2 = 1.0 / 2.0 *
                   (devstress(0) * devstress(0) + devstress(1) * devstress(1) +
                       devstress(2) * devstress(2)) +
               devstress(3) * devstress(3) + devstress(4) * devstress(4) +
               devstress(5) * devstress(5);
  ScalarT Phi_trial = 0.0;
  Phi_trial = std::sqrt(J2) + eta * p - xi * cohesion - xi * Hiso * strainbar_p;
  ScalarT Dgamma = 0.0;
  ScalarT dstrainv = 0.0;
  if (Phi_trial / abs(cohesion) > params_->abstol_)
  {
    auto returnToConeFunctAndDeriv = [this, &G, &kappa, &Phi_trial](ScalarT Dgamma_init)
    { return this->return_to_cone_funct_and_deriv(Dgamma_init, G, kappa, Phi_trial); };

    const double tol = params_->abstol_;
    Dgamma =
        Core::Utils::solve_local_newton(returnToConeFunctAndDeriv, Dgamma, tol * cohesion, itermax);
    strainbar_p = (strainbarpllast_.at(gp)) + xi * Dgamma;
    devstress.scale(1.0 - (G * Dgamma / std::sqrt(J2)));
    p = p_trial - kappa * etabar * Dgamma;
    if ((std::sqrt(J2) - G * Dgamma) / abs(cohesion) < params_->abstol_)
    {
      strainbar_p = (strainbarpllast_.at(gp));
      auto returnToApexFunctAndDeriv = [this, &p_trial, &kappa, &strainbar_p](ScalarT dstrainv_init)
      { return this->return_to_apex_funct_and_deriv(dstrainv_init, p_trial, kappa, strainbar_p); };

      const double tol = params_->abstol_;
      dstrainv = Core::Utils::solve_local_newton(
          returnToApexFunctAndDeriv, dstrainv, tol * cohesion, itermax);
      strainbar_p = (strainbarpllast_.at(gp)) + xi / eta * dstrainv;
      p = p_trial - kappa * dstrainv;
      for (int i = 0; i < 6; i++) devstress(i) = 0.0;
    }
    Stress(p, devstress, *stress);
    strain_e.update(1 / G / 2, devstress, p / kappa / 3, id2Scalar);
    for (int i = 3; i < 6; ++i) strain_e(i) *= 2.0;
    for (int i = 3; i < 6; ++i) strain(i) *= 2.0;
    strain_p.update(1.0, strain, -1.0, strain_e);

    strainplcurr_.at(gp) = Core::FADUtils::cast_to_double(strain_p);
    strainbarplcurr_.at(gp) = Core::FADUtils::cast_to_double(strainbar_p);
  }
  else
  {
    Stress(p, devstress, *stress);
    strain_e.update(trialstrain_e);
    for (int i = 3; i < 6; ++i) strain_e(i) *= 2.0;
    for (int i = 3; i < 6; ++i) strain(i) *= 2.0;

    strainplcurr_.at(gp) = strainpllast_.at(gp);
    strainbarplcurr_.at(gp) = strainbarpllast_.at(gp);
  }
  if (Phi_trial > 0)
  {
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1> devstraindouble =
        Core::FADUtils::cast_to_double(devstrain);
    if (dstrainv != 0.0)
    {
      setup_cmat_elasto_plastic_apex(*cmat, Core::FADUtils::cast_to_double(kappa), devstraindouble,
          Core::FADUtils::cast_to_double(xi), Core::FADUtils::cast_to_double(Hiso),
          Core::FADUtils::cast_to_double(eta), Core::FADUtils::cast_to_double(etabar));
    }
    else
    {
      setup_cmat_elasto_plastic_cone(*cmat, Core::FADUtils::cast_to_double(Dgamma),
          Core::FADUtils::cast_to_double(G), Core::FADUtils::cast_to_double(kappa), devstraindouble,
          Core::FADUtils::cast_to_double(xi), Core::FADUtils::cast_to_double(Hiso),
          Core::FADUtils::cast_to_double(eta), Core::FADUtils::cast_to_double(etabar));
    }
  }
  else
  {
    setup_cmat(*cmat);
  }
}

void Mat::PlasticDruckerPragerExpTest::register_output_data_names(
    std::unordered_map<std::string, int>& names_and_size) const
{
  names_and_size["accumulated_plastic_strain"] = 1;
  names_and_size["plastic_strain"] = 6;
}

bool Mat::PlasticDruckerPragerExpTest::evaluate_output_data(
    const std::string& name, Core::LinAlg::SerialDenseMatrix& data) const
{
  if (name == "accumulated_plastic_strain")
  {
    for (std::size_t gp = 0; gp < strainbarplcurr_.size(); ++gp)
    {
      data(gp, 0) = strainbarplcurr_.at(int(gp));
    }
    return true;
  }
  if (name == "plastic_strain")
  {
    for (std::size_t gp = 0; gp < strainplcurr_.size(); ++gp)
    {
      const double* values = strainplcurr_.at(gp).data();
      for (std::size_t i = 0; i < 6; ++i)
      {
        data(gp, i) = values[i];
      }
    }
    return true;
  }
  return false;
}

template <typename T>
void Mat::PlasticDruckerPragerExpTest::Stress(const T p,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1, T>& devstress,
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1, T>& stress)
{
  stress.update(devstress);
  for (int i = 0; i < 3; ++i) stress(i) += p;
}

template <typename T>
std::pair<T, T> Mat::PlasticDruckerPragerExpTest::return_to_cone_funct_and_deriv(
    T Dgamma, T G, T kappa, T Phi_trial)
{
  double Hiso = params_->isohard_;
  double eta = params_->eta_;
  double xi = params_->xi_;
  double etabar = params_->etabar_;
  T Res = Phi_trial - Dgamma * (G + eta * kappa * etabar) - (xi * xi * Dgamma * Hiso);
  T d = -G - (kappa * etabar * eta) - (xi * xi * Hiso);
  return {Res, d};
}

template <typename T>
std::pair<T, T> Mat::PlasticDruckerPragerExpTest::return_to_apex_funct_and_deriv(
    T dstrainv, T p, T kappa, T strainbar_p)
{
  double Hiso = params_->isohard_;
  double eta = params_->eta_;
  double xi = params_->xi_;
  double cohesion = params_->cohesion_;
  double etabar = params_->etabar_;
  double alpha = xi / eta;
  double beta = xi / etabar;
  T Res =
      beta * cohesion + beta * strainbar_p * Hiso - p + dstrainv * (alpha * beta * Hiso + kappa);
  T d = xi * xi / eta / etabar * Hiso + kappa;
  return {Res, d};
}

template void Mat::PlasticDruckerPragerExpTest::EvaluateFAD<double>(
    const Core::LinAlg::Matrix<3, 3>*, const Core::LinAlg::Matrix<6, 1, double>*,
    Teuchos::ParameterList&, Core::LinAlg::Matrix<6, 1, double>*, Core::LinAlg::Matrix<6, 6>*,
    int gp, int eleGID);
template void Mat::PlasticDruckerPragerExpTest::EvaluateFAD<FAD>(const Core::LinAlg::Matrix<3, 3>*,
    const Core::LinAlg::Matrix<6, 1, FAD>*, Teuchos::ParameterList&,
    Core::LinAlg::Matrix<6, 1, FAD>*, Core::LinAlg::Matrix<6, 6>*, int gp, int eleGID);
template void Mat::PlasticDruckerPragerExpTest::Stress<double>(const double p,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1, double>& devstress,
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1, double>& stress);
template void Mat::PlasticDruckerPragerExpTest::Stress<FAD>(const FAD p,
    const Core::LinAlg::Matrix<NUM_STRESS_3D, 1, FAD>& devstress,
    Core::LinAlg::Matrix<NUM_STRESS_3D, 1, FAD>& stress);
template std::pair<double, double>
Mat::PlasticDruckerPragerExpTest::return_to_cone_funct_and_deriv<double>(
    double Dgamma, double G, double kappa, double Phi_trial);
template std::pair<FAD, FAD> Mat::PlasticDruckerPragerExpTest::return_to_cone_funct_and_deriv<FAD>(
    FAD Dgamma, FAD G, FAD kappa, FAD Phi_trial);
template std::pair<double, double>
Mat::PlasticDruckerPragerExpTest::return_to_apex_funct_and_deriv<double>(
    double dstrainv, double p, double kappa, double strainbar_p);
template std::pair<FAD, FAD> Mat::PlasticDruckerPragerExpTest::return_to_apex_funct_and_deriv<FAD>(
    FAD dstrainv, FAD p, FAD kappa, FAD strainbar_p);

FOUR_C_NAMESPACE_CLOSE
