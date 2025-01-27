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
  exp_rate_.resize(numgp);
  std::fill(exp_rate_.begin(), exp_rate_.end(), 0.0);  // TODO see if there are better ways
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void Mat::StructPoroMasstransfer::update() { mat_->update(); }

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
    double& masstransfer_dJ, bool save)
{
  double temperature = 0.0;
  // TODO: do not read from parameter list!
  if (params.isParameter("scalar"))
  {
    Teuchos::RCP<std::vector<double>> scalars =
        params.get<Teuchos::RCP<std::vector<double>>>("scalar");
    temperature = scalars->at(0);
  }
  std::cout << "Temp is: " << temperature << std::endl;
  // constant function
  if (params_->functionID_ == 0)
  {
    masstransferRate = params_->rateConstant_ * temperature;
    masstransfer_dp = 0.0;
    // TODO need something to test a_dus that needs to be set here
  }
  else if (params_->functionID_ == 1)
  {
    masstransferRate = params_->rateConstant_ * press * temperature;
    masstransfer_dp = params_->rateConstant_ * temperature;
  }
  else
    FOUR_C_THROW(
        "Type of function specified by functionID %d not implemented", params_->functionID_);

  exp_rate_[gp] = masstransferRate;  // TODO see if this works well, possibly move into save?
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
