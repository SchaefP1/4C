// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_MAT_MY_EXPANSION_TEST_HPP
#define FOUR_C_MAT_MY_EXPANSION_TEST_HPP


#include "4C_config.hpp"

#include "4C_comm_parobjectfactory.hpp"
#include "4C_mat_so3_material.hpp"
#include "4C_material_parameter_base.hpp"

FOUR_C_NAMESPACE_OPEN

// forward declaration due to avoid header definition
namespace Mat
{
  namespace Elastic
  {
    class Summand;
  }  // namespace Elastic

  // forward declaration
  class MyExpansionTest_ElastHyper;

  namespace PAR
  {
    class MyExpansionTest_ElastHyper : public Core::Mat::PAR::Parameter
    {
      friend class Mat::MyExpansionTest_ElastHyper;

     public:
      /// standard constructor
      ///
      /// This constructor recursively calls the constructors of the
      /// parameter sets of the hyperelastic summands.
      MyExpansionTest_ElastHyper(const Core::Mat::PAR::Parameter::Data& matdata);

      /// @name material parameters
      //@{

      /// length of 3d elastin matrix material list
      const int nummat_elastiniso_;

      /// the list of 3d elastin matrix material IDs
      const std::vector<int> matids_elastiniso_;

      /// material mass density
      const double density_;

      /// expansion rate
      const double exp_rate_;

      /// rho at zero expansion
      const double rho_0_;

      /// rho at full expansion
      const double rho_ss_;
      //@}

      /// create material instance of matching type with my parameters
      std::shared_ptr<Core::Mat::Material> create_material() override;

    };  // class MyExpansionTest_ElastHyper

  }  // namespace PAR

  class MyExpansionTest_ElastHyperType : public Core::Communication::ParObjectType
  {
   public:
    std::string name() const override { return "MyExpansionTest_ElastHyperType"; }

    static MyExpansionTest_ElastHyperType& instance() { return instance_; };

    Core::Communication::ParObject* create(Core::Communication::UnpackBuffer& buffer) override;

   private:
    static MyExpansionTest_ElastHyperType instance_;
  };


  /*----------------------------------------------------------------------*/
  class Material;

  class MyExpansionTest_ElastHyper : public So3Material
  {
   public:
    /// construct empty material object
    MyExpansionTest_ElastHyper();

    /// construct the material object given material parameters
    explicit MyExpansionTest_ElastHyper(Mat::PAR::MyExpansionTest_ElastHyper* params);

    ///@name Packing and Unpacking
    //@{

    /// \brief Return unique ParObject id
    ///
    /// every class implementing ParObject needs a unique id defined at the
    /// top of drt_parobject.hpp (this file) and should return it in this method.
    int unique_par_object_id() const override
    {
      return MyExpansionTest_ElastHyperType::instance().unique_par_object_id();
    }

    /// \brief Pack this class so it can be communicated
    ///
    /// Resizes the vector data and stores all information of a class in it.
    /// The first information to be stored in data has to be the
    /// unique parobject id delivered by unique_par_object_id() which will then
    /// identify the exact class on the receiving processor.
    ///
    /// \param data (in/out): char vector to store class information
    void pack(Core::Communication::PackBuffer& data) const override;

    /// \brief Unpack data from a char vector into this class
    ///
    /// The vector data contains all information to rebuild the
    /// exact copy of an instance of a class on a different processor.
    /// The first entry in data has to be an integer which is the unique
    /// parobject id defined at the top of this file and delivered by
    /// unique_par_object_id().
    ///
    /// \param data (in) : vector storing all data to be unpacked into this
    ///                    instance.
    void unpack(Core::Communication::UnpackBuffer& buffer) override;

    //@}

    /// check if element kinematics and material kinematics are compatible
    void valid_kinematics(Inpar::Solid::KinemType kinem) override
    {
      if (!(kinem == Inpar::Solid::KinemType::nonlinearTotLag))
        FOUR_C_THROW("element and material kinematics are not compatible");
    }

    /// material type
    Core::Materials::MaterialType material_type() const override
    {
      return Core::Materials::m_my_expansion_test;
    }

    /// return copy of this material object
    std::shared_ptr<Core::Mat::Material> clone() const override
    {
      return std::make_shared<MyExpansionTest_ElastHyper>(*this);
    }

    /// material mass density
    double density() const override { return params_->density_; }

    void evaluate_cauchy_n_dir_and_derivatives(const Core::LinAlg::Matrix<3, 3>& defgrd,
        const Core::LinAlg::Matrix<3, 1>& n, const Core::LinAlg::Matrix<3, 1>& dir,
        double& cauchy_n_dir, Core::LinAlg::Matrix<3, 1>* d_cauchyndir_dn,
        Core::LinAlg::Matrix<3, 1>* d_cauchyndir_ddir, Core::LinAlg::Matrix<9, 1>* d_cauchyndir_dF,
        Core::LinAlg::Matrix<9, 9>* d2_cauchyndir_dF2,
        Core::LinAlg::Matrix<9, 3>* d2_cauchyndir_dF_dn,
        Core::LinAlg::Matrix<9, 3>* d2_cauchyndir_dF_ddir, int gp, int eleGID,
        const double* concentration, const double* temp, double* d_cauchyndir_dT,
        Core::LinAlg::Matrix<9, 1>* d2_cauchyndir_dF_dT) override;

    /*!
     * @brief calculates the derivatives of the hyper-elastic laws with respect to the invariants
     *
     * @param[in] prinv   Principal invariants of the elastic right Cauchy-Green tensor
     * @param[in] gp      current gauss point
     * @param[in] eleGID  Element ID
     * @param[out] dPI    First derivative w.r.t. principle invariants
     * @param[out] ddPII  Second derivative w.r.t. principle invariants
     */
    void evaluate_invariant_derivatives(const Core::LinAlg::Matrix<3, 1>& prinv, int gp, int eleGID,
        Core::LinAlg::Matrix<3, 1>& dPI, Core::LinAlg::Matrix<6, 1>& ddPII) const;

    /// hyperelastic stress response plus elasticity tensor
    void evaluate(const Core::LinAlg::Matrix<3, 3>* defgrd,  ///< Deformation gradient
        const Core::LinAlg::Matrix<6, 1>* glstrain,          ///< Green-Lagrange strain
        Teuchos::ParameterList& params,      ///< Container for additional information
        Core::LinAlg::Matrix<6, 1>* stress,  ///< 2nd Piola-Kirchhoff stresses
        Core::LinAlg::Matrix<6, 6>* cmat,    ///< Constitutive matrix
        int gp,                              ///< Gauss point
        int eleGID) override;                ///< Element ID


    /// Evaluates some kinematic quantities which are used in stress and elasticity tensor
    /// calculation
    static void EvaluateKinQuantElast(
        Core::LinAlg::Matrix<3, 3> const* const defgrd,  ///< Deformation gradient
        Core::LinAlg::Matrix<3, 3> const& iFinM,         ///< Inverse inelastic deformation gradient
        Core::LinAlg::Matrix<6, 1>& iCinv,        ///< Inverse inelastic right Cauchy-Green tensor
        Core::LinAlg::Matrix<6, 1>& iCinCiCinv,   ///< C_{in}^{-1} * C * C_{in}^{-1}
        Core::LinAlg::Matrix<6, 1>& iCv,          ///< Inverse right Cauchy-Green tensor
        Core::LinAlg::Matrix<3, 3>& iCinCM,       ///< C_{in}^{-1} * C
        Core::LinAlg::Matrix<3, 3>& iFinCeM,      ///< F_{in}^{-1} * C_e
        Core::LinAlg::Matrix<9, 1>& CiFin9x1,     ///< C * F_{in}^{-1}
        Core::LinAlg::Matrix<9, 1>& CiFinCe9x1,   ///< C * F_{in}^{-1} * C_e
        Core::LinAlg::Matrix<9, 1>& CiFiniCe9x1,  ///< C * F_{in}^{-1} * C_e^{-1}
        Core::LinAlg::Matrix<3, 1>&
            prinv,      ///< Principal invariants of elastic right Cauchy-Green tensor
        int const gp);  ///< Current Gauss-Point

    /// calculates the isotropic stress and elasticity tensor for coupled configuration
    void EvaluateIsotropicPrincElast(
        Core::LinAlg::Matrix<6, 1>& stressisoprinc,  ///< 2nd Piola Kirchhoff stress
        Core::LinAlg::Matrix<6, 6>& cmatisoprinc,    ///< Elasticity tensor
        Core::LinAlg::Matrix<6, 1> const& iCinv,  ///< Inverse inelastic right Cauchy-Green tensor
        Core::LinAlg::Matrix<6, 1> const& iCinCiCinv,  ///< C_{in}^{-1} * C * C_{in}^{-1}
        Core::LinAlg::Matrix<6, 1> const& iCv,         ///< Inverse right Cauchy-Green tensor
        Core::LinAlg::Matrix<3, 1> const& gamma,       ///< Factors for stress calculation
        Core::LinAlg::Matrix<8, 1> const& delta)
        const;  ///< Factors for elasticity tensor calculation


    /// setup
    void setup(int numgp, const Core::IO::InputParameterContainer& container) override;

    /// update
    void update() override;

    /// Return quick accessible material parameter data
    Core::Mat::PAR::Parameter* parameter() const override { return params_; }

    /// Return names of visualization data
    void vis_names(std::map<std::string, int>& names) const override;

    /// Return visualization data
    bool vis_data(
        const std::string& name, std::vector<double>& data, int numgp, int eleID) const override;

   private:
    /// My material parameters
    Mat::PAR::MyExpansionTest_ElastHyper* params_;

    /// Map to elastin 3d matrix material summands
    std::vector<std::shared_ptr<Mat::Elastic::Summand>> potsumeliso_;

    /// average rho para for exp. Hack and only works for homogeneous exp but ok for test
    /// TODO: find better option
    double rho_average_;
  };

}  // namespace Mat

FOUR_C_NAMESPACE_CLOSE

#endif /* FOUR_C_MAT_MY_EXPANSION_TEST_HPP */
