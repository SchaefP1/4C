// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_STRUCTPORO_MASSTRANSFER_HPP
#define FOUR_C_MAT_STRUCTPORO_MASSTRANSFER_HPP

#include "4C_config.hpp"

#include "4C_mat_structporo.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Mat
{
  // forward declaration
  class StructPoroMasstransfer;

  namespace PAR
  {
    class StructPoroMasstransfer : public PAR::StructPoro
    {
      friend class Mat::StructPoroMasstransfer;

     public:
      /// standard constructor
      StructPoroMasstransfer(const Core::Mat::PAR::Parameter::Data& matdata);

      /// destructor
      virtual ~StructPoroMasstransfer() { ; }

      /// create material instance of matching type with my parameters
      std::shared_ptr<Core::Mat::Material> create_material() override;

      /// @name material parameters
      //@{

      /// determine if to use a constant rate function (0) or linear in pressure(1)
      int functionID_;

      /// the constant or linear factor determining the mass transfer rate depending on functionID_
      double rateConstant_;
      //@}

    };  // class StructPoroMasstransfer

  }  // namespace PAR

  class StructPoroMasstransferType : public Core::Communication::ParObjectType
  {
   public:
    std::string name() const { return "StructPoroMasstransferType"; }

    static StructPoroMasstransferType& instance() { return instance_; };

    Core::Communication::ParObject* create(Core::Communication::UnpackBuffer& buffer) override;

   private:
    static StructPoroMasstransferType instance_;
  };

  /*----------------------------------------------------------------------*/
  //! Wrapper for StructPoro material including mass transfer between sceleton and fluid
  //!
  //! This object exists (several times) at every element

  /*!
    This class is a child of the regular StructPoro material, but has additional functions to
    account for mass transfer between the sceleton and fluid phase. This can account for example for
    an exchange of mass between an occulded porosity which moves with the sceleton and the flowing
    fluid or for fluid getting absorbed by the sceleton like in the case of metal hydride storage
    tanks.

    The main methods of this class are:
    ...
   */

  class StructPoroMasstransfer : public StructPoro
  {
   public:
    /// construct empty material object
    StructPoroMasstransfer();

    /// construct the material object given material parameters
    explicit StructPoroMasstransfer(Mat::PAR::StructPoroMasstransfer* params);

    //! @name Packing and Unpacking

    /*!
     \brief Return unique ParObject id

     every class implementing ParObject needs a unique id defined at the
     top of parobject.H (this file) and should return it in this method.
     */
    int unique_par_object_id() const override
    {
      return StructPoroMasstransferType::instance().unique_par_object_id();
    }

    /*!
     \brief Pack this class so it can be communicated

     Resizes the vector data and stores all information of a class in it.
     The first information to be stored in data has to be the
     unique parobject id delivered by unique_par_object_id() which will then
     identify the exact class on the receiving processor.

     \param data (in/out): char vector to store class information
     */
    void pack(Core::Communication::PackBuffer& data) const override;

    /*!
     \brief Unpack data from a char vector into this class

     The vector data contains all information to rebuild the
     exact copy of an instance of a class on a different processor.
     The first entry in data has to be an integer which is the unique
     parobject id defined at the top of this file and delivered by
     unique_par_object_id().

     \param data (in) : vector storing all data to be unpacked into this
     instance.
     */
    void unpack(Core::Communication::UnpackBuffer& buffer) override;

    //@}

    //! @name Access methods

    /// material type
    Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::m_structporomasstransfer;
    }

    /// return copy of this material object
    std::shared_ptr<Core::Mat::Material> clone() const override
    {
      return std::make_shared<StructPoroMasstransfer>(*this);
    }

    /// Initialize internal variables
    void poro_setup(int numgp,  ///< number of Gauss points
        const Core::IO::InputParameterContainer& container);

    void update() override;

    /// Return quick accessible material parameter data
    Core::Mat::PAR::Parameter* parameter() const override { return params_; }

    //! compute current masstransfer and its derivative and save it
    void ComputeMasstransfer(Teuchos::ParameterList& params,  //!< (i) element parameter list
        double press,                                         //!< (i) pressure at gauss point
        int gp,                                               //!< (i) number of current gauss point
        double& masstransferRate,   //!< (o) rate of masstransfer between sceleton and fluid at gp
        double& masstransfer_dp,    //!< (o) derivative of masstransfer wrt pressure
        double& masstransfer_dphi,  //!< (o) derivative of masstransfer wrt porosity
        double& masstransfer_dJ,    //!< (o) derivative of masstransfer wrt jacobian
        bool save = true);

    //@}

    //! @name Visualization methods

    /// Return names of visualization data
    virtual void vis_names(std::map<std::string, int>& names);

    /// Return visualization data
    virtual bool vis_data(const std::string& name, std::vector<double>& data, int numgp, int eleID);

    //@}

   protected:
    /// my material parameters
    Mat::PAR::StructPoroMasstransfer* params_;

    /// current expansion rate at each gp
    std::vector<double> exp_rate_;
    /// current expansion rate at each gp
    std::vector<double> masstransferRate_;
    /// current expansion rate at each gp
    std::vector<double> masstransfer_dp_;
  };

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif /* FOUR_C_MAT_STRUCTPORO_MASSTRANSFER_HPP */
