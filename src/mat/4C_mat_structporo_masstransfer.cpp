// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_mat_structporo_masstransfer.hpp"

#include "4C_comm_pack_helpers.hpp"
#include "4C_global_data.hpp"
#include "4C_mat_my_expansion_test.hpp"
#include "4C_mat_par_bundle.hpp"

#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::PAR::StructPoroMasstransfer::StructPoroMasstransfer(
    const Core::Mat::PAR::Parameter::Data& matdata)
    : StructPoro(matdata),
      functionID_(matdata.parameters.get<int>("FUNCTIONID")),
      rateConstant_(matdata.parameters.get<double>("RATECONSTANT"))
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<Core::Mat::Material> Mat::PAR::StructPoroMasstransfer::create_material()
{
  return std::make_shared<Mat::StructPoroMasstransfer>(this);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::StructPoroMasstransferType Mat::StructPoroMasstransferType::instance_;

Core::Communication::ParObject* Mat::StructPoroMasstransferType::create(
    Core::Communication::UnpackBuffer& buffer)
{
  Mat::StructPoroMasstransfer* struct_poro = new Mat::StructPoroMasstransfer();
  struct_poro->unpack(buffer);
  return struct_poro;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::StructPoroMasstransfer::StructPoroMasstransfer() : params_(NULL) {}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Mat::StructPoroMasstransfer::StructPoroMasstransfer(Mat::PAR::StructPoroMasstransfer* params)
    : StructPoro(params), params_(params)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mat::StructPoroMasstransfer::poro_setup(
    int numgp, const Core::IO::InputParameterContainer& container)
{
  StructPoro::poro_setup(numgp, container);
  masstransferRate_.resize(numgp);
  masstransfer_dp_.resize(numgp);
  density_n_.resize(numgp);
  density_current_.resize(numgp);
  temperature_current_.resize(numgp);
  rho_s_current_.resize(numgp);
  dtemperature_dt_.resize(numgp);
  density_set_flag.resize(numgp);
  std::fill(masstransferRate_.begin(), masstransferRate_.end(), 0.0);
  std::fill(masstransfer_dp_.begin(), masstransfer_dp_.end(), 0.0);
  std::fill(density_n_.begin(), density_n_.end(), 1.0);
  std::fill(density_current_.begin(), density_current_.end(), 1.0);
  std::fill(temperature_current_.begin(), temperature_current_.end(), 1.0);
  std::fill(rho_s_current_.begin(), rho_s_current_.end(), 0.0);
  std::fill(dtemperature_dt_.begin(), dtemperature_dt_.end(), 0.0);
  std::fill(density_set_flag.begin(), density_set_flag.end(), false);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mat::StructPoroMasstransfer::update()
{
  mat_->update();
  density_n_ = density_current_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::StructPoroMasstransfer::pack(Core::Communication::PackBuffer& data) const
{
  // pack type of this instance of ParObject
  int type = unique_par_object_id();
  add_to_pack(data, type);

  // matid
  int matid = -1;
  if (params_ != nullptr) matid = params_->id();  // in case we are in post-process mode
  add_to_pack(data, matid);

  // add base class material
  StructPoro::pack(data);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Mat::StructPoroMasstransfer::unpack(Core::Communication::UnpackBuffer& buffer)
{
  Core::Communication::extract_and_assert_id(buffer, unique_par_object_id());

  // matid
  int matid;
  extract_from_pack(buffer, matid);
  params_ = nullptr;
  if (Global::Problem::instance()->materials() != nullptr)
    if (Global::Problem::instance()->materials()->num() != 0)
    {
      const int probinst = Global::Problem::instance()->materials()->get_read_from_problem();
      Core::Mat::PAR::Parameter* mat =
          Global::Problem::instance(probinst)->materials()->parameter_by_id(matid);
      if (mat->type() == material_type())
        params_ = static_cast<Mat::PAR::StructPoroMasstransfer*>(mat);
      else
        FOUR_C_THROW("Type of parameter material %d does not fit to calling type %d", mat->type(),
            material_type());
    }

  // extract base class material
  StructPoro::unpack(buffer);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mat::StructPoroMasstransfer::ComputeMasstransfer(Teuchos::ParameterList& params, double press,
    int gp, double& masstransferRate, double& masstransfer_dp, double& masstransfer_dphi,
    double& masstransfer_dJ, double& density_current, double& density_n, double& ddensity_dp,
    double& ddensity_dT, double& dddensity_dTdp, double& dtemperature_dt, bool save)
{  // std::cout << "begin " << std::endl;
  double temperature = 1.0;
  double rho_s;
  double rho_ss = 8520.0;
  double t_tot = params.get<double>("total time");
  // TODO: do not read from parameter list!

  //  constant function

  if (params.isParameter("scalardt"))
  {
    std::shared_ptr<std::vector<double>> scalardt =
        params.get<std::shared_ptr<std::vector<double>>>("scalardt");
    dtemperature_dt_[gp] = scalardt->at(2);
  }

  if (params.isParameter("scalar"))
  {
    std::shared_ptr<std::vector<double>> scalars =
        params.get<std::shared_ptr<std::vector<double>>>("scalar");
    if (scalars->size() >
        1)  // in fluid_ele_poro_calc I only have access to one scalar (I know where this happens
            // but don't know why it is done like this) So I just do the calc in solid when I have
            // all scalars and then save the result for the gp to reuse in fluid. Horrible hack I
            // should fix and remains to be seen if it works
    {
      temperature = scalars->at(2);
      rho_s = scalars->at(1);

      // save value to use it later or when called without scalar available
      temperature_current_[gp] = temperature;
      rho_s_current_[gp] = rho_s;
    }
  }

  if (params_->functionID_ == 0)
  {
    masstransferRate = params_->rateConstant_;
    masstransfer_dp = 0.0;
  }
  else if (params_->functionID_ == 1)
  {
    masstransferRate = params_->rateConstant_ * press;
    masstransfer_dp = params_->rateConstant_;
  }
  else if (params_->functionID_ == 2)
  {
    if (t_tot >= 0.1)
    {
      // std::cout << "temperature " << temperature << " rho_s " << rho_s << " press " << press
      // << std::endl;
      masstransferRate_[gp] =
          -59.9 * exp(-21.17e03 / 8.314 / temperature_current_[gp]) *
          (log(press) - 12.919 - log(1.0e05) + 3704.4 / temperature_current_[gp]) *
          (rho_ss - rho_s_current_[gp]);
      masstransfer_dp_[gp] = -59.9 * exp(-21.17e03 / 8.314 / temperature_current_[gp]) *
                             (1 / press) * (rho_ss - rho_s_current_[gp]);
      // std::cout << "passed " << std::endl;
    }

    // This would be for weakly compressible
    if (true)
    {
      double ref_press = 800000.0;
      double MH_over_R = 2.01588e-03 / 8.314;
      density_current = MH_over_R * ref_press / temperature_current_[gp];
      ddensity_dp = 0.0;
      ddensity_dT = -MH_over_R * ref_press / temperature_current_[gp] / temperature_current_[gp];
    }
    // this would be ideal gas
    else
    {
      double MH_over_R = 2.01588e-03 / 8.314;
      density_current = MH_over_R * press / temperature_current_[gp];
      ddensity_dp = MH_over_R / temperature_current_[gp];
      ddensity_dT = -MH_over_R * press / temperature_current_[gp] / temperature_current_[gp];
    }
    density_current_[gp] = density_current;
    if (!density_set_flag[gp])
    {
      density_n_[gp] = density_current;
      density_set_flag[gp] = true;
    }
    density_n = density_n_[gp];

    dtemperature_dt = dtemperature_dt_[gp];
    masstransferRate = masstransferRate_[gp];
    masstransfer_dp = masstransfer_dp_[gp];
  }
  else if (params_->functionID_ == 3 || params_->functionID_ == 4)
  {  // like functionID 0 (constant reaction) just comp
    masstransferRate = params_->rateConstant_;
    masstransfer_dp = 0.0;
    dtemperature_dt = dtemperature_dt_[gp];
    if (params_->functionID_ == 3)
    {
      double ref_press = 1.0;
      double MH_over_R = 1.0;
      density_current = MH_over_R * ref_press / temperature_current_[gp];
      ddensity_dp = 0.0;
      ddensity_dT = -MH_over_R * ref_press / temperature_current_[gp] / temperature_current_[gp];
      dddensity_dTdp = 0.0;
    }
    // this would be ideal gas
    if (params_->functionID_ == 4)
    {
      double MH_over_R = 1.0;
      density_current = MH_over_R * press / temperature_current_[gp];
      ddensity_dp = MH_over_R / temperature_current_[gp];
      ddensity_dT = -MH_over_R * press / temperature_current_[gp] / temperature_current_[gp];
      dddensity_dTdp = -MH_over_R / temperature_current_[gp] / temperature_current_[gp];
    }
    density_current_[gp] = density_current;
    // if the was no previous step so the density there hasn't been set so use current one
    if (!density_set_flag[gp])
    {
      density_n_[gp] = density_current;
      density_set_flag[gp] = true;
    }
    density_n = density_n_[gp];
  }
  else
    FOUR_C_THROW(
        "Type of function specified by functionID %d not implemented", params_->functionID_);
  // std::cout << "end " << std::endl;
  // std::cout << "temperature_current_ " << std::setprecision(10) << temperature_current_[gp] <<
  // std::endl;
  if (save)
  {
    (void)gp;
    // maybe save masstransferRate and masstransfer_dp here, not used yet
  }
}
/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mat::StructPoroMasstransfer::vis_names(std::map<std::string, int>& names)
{
  // call base class
  StructPoro::vis_names(names);
  //  std::string name = "reference_porosity";
  //  names[name] = 1;  // scalar
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool Mat::StructPoroMasstransfer::vis_data(
    const std::string& name, std::vector<double>& data, int numgp, int eleID)
{
  // call base class
  if (StructPoro::vis_data(name, data, numgp, eleID)) return true;
  //  if (name == "reference_porosity")
  //  {
  //    if ((int)data.size() != 1) FOUR_C_THROW("size mismatch");
  //    data[0] = RefPorosityAv();
  //    return true;
  //  }
  return false;
}

FOUR_C_NAMESPACE_CLOSE
