// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_particle_interaction_sph_artificialviscosity.hpp"

#include "4C_particle_interaction_utils.hpp"
#include "4C_unittest_utils_assertions_test.hpp"

namespace
{
  using namespace FourC;

  class SPHArtificialViscosityTest : public ::testing::Test
  {
   protected:
    std::unique_ptr<ParticleInteraction::SPHArtificialViscosity> artificialviscosity_;

    SPHArtificialViscosityTest()
    {
      // create artificial viscosity handler
      artificialviscosity_ = std::make_unique<ParticleInteraction::SPHArtificialViscosity>();

      // init artificial viscosity handler
      artificialviscosity_->init();

      // setup artificial viscosity handler
      artificialviscosity_->setup();
    }
    // note: the public functions init() and setup() of class SPHEquationOfStateGenTait are called
    // in SetUp() and thus implicitly tested by all following unittests
  };

  TEST_F(SPHArtificialViscosityTest, ArtificialViscosity)
  {
    const double dens_i = 1.01;
    const double dens_j = 0.97;

    double vel_i[3];
    vel_i[0] = 0.2;
    vel_i[1] = 0.3;
    vel_i[2] = -0.12;
    double vel_j[3];
    vel_j[0] = -0.12;
    vel_j[1] = -0.97;
    vel_j[2] = 0.98;

    const double mass_i = 5.85;
    const double mass_j = 10.32;
    const double artvisc_i = 0.1;
    const double artvisc_j = 0.2;
    const double dWdrij = 0.76;
    const double dWdrji = 0.89;
    const double h_i = 0.2;
    const double h_j = 0.25;
    const double c_i = 10.0;
    const double c_j = 12.5;
    const double abs_rij = 0.3;

    double e_ij[3];
    e_ij[0] = 1.0 / std::sqrt(21);
    e_ij[1] = 2.0 / std::sqrt(21);
    e_ij[2] = 4.0 / std::sqrt(21);

    double acc_i[3] = {0.0};
    double acc_j[3] = {0.0};

    // compute reference solution
    const double h_ij = 0.5 * (h_i + h_j);
    const double c_ij = 0.5 * (c_i + c_j);
    const double dens_ij = 0.5 * (dens_i + dens_j);
    const double e_ij_vrel_ij = ((vel_i[0] - vel_j[0]) * e_ij[0] + (vel_i[1] - vel_j[1]) * e_ij[1] +
                                 (vel_i[2] - vel_j[2]) * e_ij[2]);
    const double epsilon = 0.01;
    const double fac = h_ij * c_ij * e_ij_vrel_ij * abs_rij /
                       (dens_ij * (std::pow(abs_rij, 2) + epsilon * std::pow(h_ij, 2)));

    double acc_i_ref[3];
    ParticleInteraction::Utils::vec_set_scale(acc_i_ref, (artvisc_i * mass_j * dWdrij * fac), e_ij);

    double acc_j_ref[3];
    ParticleInteraction::Utils::vec_set_scale(
        acc_j_ref, (-artvisc_j * mass_i * dWdrji * fac), e_ij);

    artificialviscosity_->artificial_viscosity(vel_i, vel_j, &mass_i, &mass_j, artvisc_i, artvisc_j,
        dWdrij, dWdrji, dens_ij, h_ij, c_ij, abs_rij, e_ij, acc_i, acc_j);

    // compare results
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_i, acc_i_ref, 3, 1.0e-14);
    FOUR_C_EXPECT_ITERABLE_NEAR(acc_j, acc_j_ref, 3, 1.0e-14);
  }
}  // namespace
