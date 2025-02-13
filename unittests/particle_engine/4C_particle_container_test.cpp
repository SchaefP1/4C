// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include <gtest/gtest.h>

#include "4C_particle_engine_container.hpp"
#include "4C_unittest_utils_assertions_test.hpp"


namespace
{
  using namespace FourC;

  class ParticleContainerTest : public ::testing::Test
  {
   protected:
    std::unique_ptr<PARTICLEENGINE::ParticleContainer> container_;

    int statesvectorsize_;

    ParticleContainerTest()
    {
      const int size = 7;
      std::set<PARTICLEENGINE::StateEnum> stateEnumSet = {
          PARTICLEENGINE::Position, PARTICLEENGINE::Velocity, PARTICLEENGINE::Mass};

      // create, init and setup container
      container_ = std::make_unique<PARTICLEENGINE::ParticleContainer>();
      container_->init();
      container_->setup(size, stateEnumSet);

      const int maximum_stored_state_enum_set_value{*(--stateEnumSet.end())};
      statesvectorsize_ = maximum_stored_state_enum_set_value + 1;

      // init some particles
      int index(0);
      int globalid(0);

      PARTICLEENGINE::ParticleStates particle;
      particle.assign(statesvectorsize_, std::vector<double>{});

      // first particle
      {
        globalid = 1;
        particle = create_test_particle({1.20, 0.70, 2.10}, {0.23, 1.76, 3.89}, {0.12});
        container_->add_particle(index, globalid, particle);
      }

      // second particle
      {
        globalid = 2;
        particle = create_test_particle({-1.05, 12.6, -8.54}, {0.25, -21.5, 1.0}, {12.34});
        container_->add_particle(index, globalid, particle);
      }

      // third particle
      {
        globalid = 3;
        particle = create_test_particle({61.0, -2.63, 0.11}, {-7.35, -5.98, 1.11}, {0.5});
        container_->add_particle(index, globalid, particle);
      }
    }

    PARTICLEENGINE::ParticleStates create_test_particle(
        std::vector<double> pos, std::vector<double> vel, std::vector<double> mass)
    {
      PARTICLEENGINE::ParticleStates particle;
      particle.assign(statesvectorsize_, std::vector<double>{});

      particle[PARTICLEENGINE::Position] = pos;
      particle[PARTICLEENGINE::Velocity] = vel;
      particle[PARTICLEENGINE::Mass] = mass;

      return particle;
    }

    // note: the public functions init(), setup() and AddParticle() of class ParticleContainer are
    // called in the constructor and thus implicitly tested by all following unittests
  };

  void compare_particle_states(
      PARTICLEENGINE::ParticleStates& particle_reference, PARTICLEENGINE::ParticleStates& particle)
  {
    ASSERT_EQ(particle_reference.size(), particle.size());

    for (std::size_t i = 0; i < particle.size(); ++i)
    {
      std::vector<double>& state_reference = particle_reference[i];
      std::vector<double>& state = particle[i];

      ASSERT_EQ(state_reference.size(), state.size());

      for (std::size_t j = 0; j < state_reference.size(); ++j)
        EXPECT_NEAR(state_reference[j], state[j], 1e-14)
            << "state '"
            << PARTICLEENGINE::enum_to_state_name(static_cast<PARTICLEENGINE::ParticleState>(i))
            << "' j = " << j;
    }
  }

  TEST_F(ParticleContainerTest, increase_container_size)
  {
    container_->increase_container_size();
    EXPECT_EQ(container_->particles_stored(), 3);
    EXPECT_EQ(container_->container_size(), 14);
  }

  TEST_F(ParticleContainerTest, decrease_container_size)
  {
    container_->decrease_container_size();
    EXPECT_EQ(container_->particles_stored(), 3);
    EXPECT_EQ(container_->container_size(), 3);
  }

  TEST_F(ParticleContainerTest, check_and_decrease_container_size)
  {
    container_->check_and_decrease_container_size();
    EXPECT_EQ(container_->particles_stored(), 3);
    EXPECT_EQ(container_->container_size(), 3);
  }

  TEST_F(ParticleContainerTest, ClearContainer)
  {
    container_->clear_container();
    EXPECT_EQ(container_->particles_stored(), 0);
  }

  TEST_F(ParticleContainerTest, AddParticle)
  {
    // init a particle
    int globalid(0);

    PARTICLEENGINE::ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>{});

    globalid = 4;
    particle = create_test_particle({-1.23, 1.70, 9.10}, {6.23, 2.3, 6.9}, {5.12});

    int index(0);
    container_->add_particle(index, globalid, particle);
    EXPECT_EQ(index, 3);
    EXPECT_EQ(container_->particles_stored(), 4);
  }

  TEST_F(ParticleContainerTest, ReplaceParticle)
  {
    // init a particle
    int globalid(0);

    PARTICLEENGINE::ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>{});
    PARTICLEENGINE::ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>{});

    particle_reference = create_test_particle({-1.23, 1.70, 9.10}, {6.23, 2.3, 6.9}, {5.12});

    int index = 0;

    // replace only states and leave global id untouched
    container_->replace_particle(index, -1, particle_reference);
    EXPECT_EQ(container_->particles_stored(), 3);
    container_->get_particle(index, globalid, particle);
    EXPECT_EQ(globalid, 1);
    compare_particle_states(particle_reference, particle);

    // also replace global id
    container_->replace_particle(index, 4, particle_reference);
    EXPECT_EQ(container_->particles_stored(), 3);
    container_->get_particle(index, globalid, particle);
    EXPECT_EQ(globalid, 4);
    compare_particle_states(particle_reference, particle);
  }

  TEST_F(ParticleContainerTest, GetParticle)
  {
    int globalid(0);
    int globalid_reference(0);

    PARTICLEENGINE::ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>{});
    PARTICLEENGINE::ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>{});

    for (int index = 0; index < 3; ++index)
    {
      SCOPED_TRACE("Particle " + std::to_string(index));
      if (index == 0)
      {
        globalid_reference = 1;
        particle_reference = create_test_particle({1.20, 0.70, 2.10}, {0.23, 1.76, 3.89}, {0.12});
      }
      else if (index == 1)
      {
        globalid_reference = 2;
        particle_reference =
            create_test_particle({-1.05, 12.6, -8.54}, {0.25, -21.5, 1.0}, {12.34});
      }
      else if (index == 2)
      {
        globalid_reference = 3;
        particle_reference = create_test_particle({61.0, -2.63, 0.11}, {-7.35, -5.98, 1.11}, {0.5});
      }

      container_->get_particle(index, globalid, particle);
      EXPECT_EQ(globalid_reference, globalid);
      compare_particle_states(particle_reference, particle);
    }
  }

  TEST_F(ParticleContainerTest, RemoveParticle)
  {
    int globalid(0);
    int globalid_reference(0);

    PARTICLEENGINE::ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>{});
    PARTICLEENGINE::ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>{});

    container_->remove_particle(0);
    EXPECT_EQ(container_->particles_stored(), 2);

    for (int index = 0; index < 2; ++index)
    {
      SCOPED_TRACE("Particle " + std::to_string(index));
      if (index == 0)
      {
        globalid_reference = 3;
        particle_reference = create_test_particle({61.0, -2.63, 0.11}, {-7.35, -5.98, 1.11}, {0.5});
      }
      else if (index == 1)
      {
        globalid_reference = 2;
        particle_reference =
            create_test_particle({-1.05, 12.6, -8.54}, {0.25, -21.5, 1.0}, {12.34});
      }

      container_->get_particle(index, globalid, particle);
      EXPECT_EQ(globalid_reference, globalid);
      compare_particle_states(particle_reference, particle);
    }
  }

  TEST_F(ParticleContainerTest, GetStateDim)
  {
    EXPECT_EQ(container_->get_state_dim(PARTICLEENGINE::Position), 3);
    EXPECT_EQ(container_->get_state_dim(PARTICLEENGINE::Velocity), 3);
    EXPECT_EQ(container_->get_state_dim(PARTICLEENGINE::Mass), 1);
  }

  TEST_F(ParticleContainerTest, GetPtrToState)
  {
    std::array<double, 3> pos = {0.0};
    std::array<double, 3> vel = {0.0};
    std::array<double, 1> mass = {0.0};

    for (int index = 0; index < 3; ++index)
    {
      if (index == 0)
      {
        pos[0] = 1.20;
        pos[1] = 0.70;
        pos[2] = 2.10;
        vel[0] = 0.23;
        vel[1] = 1.76;
        vel[2] = 3.89;
        mass[0] = 0.12;
      }
      else if (index == 1)
      {
        pos[0] = -1.05;
        pos[1] = 12.6;
        pos[2] = -8.54;
        vel[0] = 0.25;
        vel[1] = -21.5;
        vel[2] = 1.0;
        mass[0] = 12.34;
      }
      else if (index == 2)
      {
        pos[0] = 61.0;
        pos[1] = -2.63;
        pos[2] = 0.11;
        vel[0] = -7.35;
        vel[1] = -5.98;
        vel[2] = 1.11;
        mass[0] = 0.5;
      }

      const double* currpos = container_->get_ptr_to_state(PARTICLEENGINE::Position, index);
      FOUR_C_EXPECT_ITERABLE_NEAR(currpos, pos.begin(), 3, 1e-14);

      const double* currvel = container_->get_ptr_to_state(PARTICLEENGINE::Velocity, index);
      FOUR_C_EXPECT_ITERABLE_NEAR(currvel, vel.begin(), 3, 1e-14);

      const double* currmass = container_->get_ptr_to_state(PARTICLEENGINE::Mass, index);
      EXPECT_NEAR(currmass[0], mass[0], 1e-14);
    }
  }

  TEST_F(ParticleContainerTest, CondGetPtrToState)
  {
    std::array<double, 3> pos = {0.0};
    std::array<double, 3> vel = {0.0};
    std::array<double, 1> mass = {0.0};

    for (int index = 0; index < 3; ++index)
    {
      if (index == 0)
      {
        pos[0] = 1.20;
        pos[1] = 0.70;
        pos[2] = 2.10;
        vel[0] = 0.23;
        vel[1] = 1.76;
        vel[2] = 3.89;
        mass[0] = 0.12;
      }
      else if (index == 1)
      {
        pos[0] = -1.05;
        pos[1] = 12.6;
        pos[2] = -8.54;
        vel[0] = 0.25;
        vel[1] = -21.5;
        vel[2] = 1.0;
        mass[0] = 12.34;
      }
      else if (index == 2)
      {
        pos[0] = 61.0;
        pos[1] = -2.63;
        pos[2] = 0.11;
        vel[0] = -7.35;
        vel[1] = -5.98;
        vel[2] = 1.11;
        mass[0] = 0.5;
      }

      double* currpos = container_->cond_get_ptr_to_state(PARTICLEENGINE::Position, index);
      FOUR_C_EXPECT_ITERABLE_NEAR(currpos, pos.begin(), 3, 1.0e-14);

      double* currvel = container_->cond_get_ptr_to_state(PARTICLEENGINE::Velocity, index);
      FOUR_C_EXPECT_ITERABLE_NEAR(currvel, vel.begin(), 3, 1.0e-14);

      double* currmass = container_->cond_get_ptr_to_state(PARTICLEENGINE::Mass, index);
      EXPECT_NEAR(currmass[0], mass[0], 1e-14);
    }
  }

  TEST_F(ParticleContainerTest, CondGetPtrToStateNotStored)
  {
    double* currpos = container_->cond_get_ptr_to_state(PARTICLEENGINE::Acceleration, 0);
    EXPECT_EQ(currpos, nullptr);
  }

  TEST_F(ParticleContainerTest, GetPtrToGlobalID)
  {
    int globalid_reference(0);

    for (int index = 0; index < 3; ++index)
    {
      if (index == 0)
        globalid_reference = 1;
      else if (index == 1)
        globalid_reference = 2;
      else if (index == 2)
        globalid_reference = 3;

      int* globalid = container_->get_ptr_to_global_id(index);
      EXPECT_EQ(globalid[0], globalid_reference);
    }
  }

  TEST_F(ParticleContainerTest, ScaleState)
  {
    int globalid(0);

    PARTICLEENGINE::ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>{});
    PARTICLEENGINE::ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>{});

    container_->scale_state(1.5, PARTICLEENGINE::Position);
    container_->scale_state(3.25, PARTICLEENGINE::Velocity);
    container_->scale_state(0.95, PARTICLEENGINE::Mass);

    for (int index = 0; index < 3; ++index)
    {
      SCOPED_TRACE("Particle " + std::to_string(index));
      if (index == 0)
      {
        particle_reference =
            create_test_particle({1.8, 1.05, 3.15}, {0.7475, 5.72, 12.6425}, {0.114});
      }
      else if (index == 1)
      {
        particle_reference =
            create_test_particle({-1.575, 18.9, -12.81}, {0.8125, -69.875, 3.25}, {11.723});
      }
      else if (index == 2)
      {
        particle_reference =
            create_test_particle({91.5, -3.945, 0.165}, {-23.8875, -19.435, 3.6075}, {0.475});
      }

      container_->get_particle(index, globalid, particle);

      compare_particle_states(particle_reference, particle);
    }
  }

  TEST_F(ParticleContainerTest, UpdateState)
  {
    int globalid(0);

    PARTICLEENGINE::ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>{});
    PARTICLEENGINE::ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>{});

    container_->update_state(1.0, PARTICLEENGINE::Position, 0.5, PARTICLEENGINE::Velocity);

    for (int index = 0; index < 3; ++index)
    {
      SCOPED_TRACE("Particle " + std::to_string(index));
      if (index == 0)
      {
        particle_reference = create_test_particle({1.315, 1.58, 4.045}, {0.23, 1.76, 3.89}, {0.12});
      }
      else if (index == 1)
      {
        particle_reference =
            create_test_particle({-0.925, 1.85, -8.04}, {0.25, -21.5, 1.0}, {12.34});
      }
      else if (index == 2)
      {
        particle_reference =
            create_test_particle({57.325, -5.62, 0.665}, {-7.35, -5.98, 1.11}, {0.5});
      }

      container_->get_particle(index, globalid, particle);

      compare_particle_states(particle_reference, particle);
    }
  }

  TEST_F(ParticleContainerTest, set_state)
  {
    int globalid(0);

    PARTICLEENGINE::ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>{});
    PARTICLEENGINE::ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>{});

    std::vector<double> pos(3);
    std::vector<double> vel(3);
    std::vector<double> mass(1);

    pos[0] = 3.15;
    pos[1] = -1.45;
    pos[2] = 9.5;
    vel[0] = -21.30;
    vel[1] = -4.33;
    vel[2] = 0.933;
    mass[0] = 1.234;

    particle_reference = create_test_particle(pos, vel, mass);

    container_->set_state(pos, PARTICLEENGINE::Position);
    container_->set_state(vel, PARTICLEENGINE::Velocity);
    container_->set_state(mass, PARTICLEENGINE::Mass);

    for (int index = 0; index < 3; ++index)
    {
      container_->get_particle(index, globalid, particle);
      compare_particle_states(particle_reference, particle);
    }
  }

  TEST_F(ParticleContainerTest, ClearState)
  {
    int globalid(0);

    PARTICLEENGINE::ParticleStates particle;
    particle.assign(statesvectorsize_, std::vector<double>{});
    PARTICLEENGINE::ParticleStates particle_reference;
    particle_reference.assign(statesvectorsize_, std::vector<double>{});

    particle_reference = create_test_particle({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0});

    container_->clear_state(PARTICLEENGINE::Position);
    container_->clear_state(PARTICLEENGINE::Velocity);
    container_->clear_state(PARTICLEENGINE::Mass);

    for (int index = 0; index < 3; ++index)
    {
      container_->get_particle(index, globalid, particle);
      compare_particle_states(particle_reference, particle);
    }
  }

  TEST_F(ParticleContainerTest, GetMinValueOfState)
  {
    EXPECT_NEAR(container_->get_min_value_of_state(PARTICLEENGINE::Mass), 0.12, 1e-14);
    EXPECT_NEAR(container_->get_min_value_of_state(PARTICLEENGINE::Position), -8.54, 1e-14);
  }

  TEST_F(ParticleContainerTest, GetMaxValueOfState)
  {
    EXPECT_NEAR(container_->get_max_value_of_state(PARTICLEENGINE::Mass), 12.34, 1e-14);
    EXPECT_NEAR(container_->get_max_value_of_state(PARTICLEENGINE::Position), 61.0, 1e-14);
  }

  TEST_F(ParticleContainerTest, GetStoredStates)
  {
    const std::set<PARTICLEENGINE::StateEnum>& particleStates = container_->get_stored_states();
    EXPECT_EQ(particleStates.size(), 3);
    EXPECT_TRUE(particleStates.find(PARTICLEENGINE::Position) != particleStates.end());
    EXPECT_TRUE(particleStates.find(PARTICLEENGINE::Velocity) != particleStates.end());
    EXPECT_TRUE(particleStates.find(PARTICLEENGINE::Mass) != particleStates.end());
  }

  TEST_F(ParticleContainerTest, HaveStoredState)
  {
    EXPECT_TRUE(container_->have_stored_state(PARTICLEENGINE::Position));
    EXPECT_TRUE(container_->have_stored_state(PARTICLEENGINE::Velocity));
    EXPECT_TRUE(container_->have_stored_state(PARTICLEENGINE::Mass));

    EXPECT_FALSE(container_->have_stored_state(PARTICLEENGINE::Acceleration));
    EXPECT_FALSE(container_->have_stored_state(PARTICLEENGINE::Density));
  }

  TEST_F(ParticleContainerTest, ContainerSize) { EXPECT_EQ(container_->container_size(), 7); }

  TEST_F(ParticleContainerTest, ParticlesStored) { EXPECT_EQ(container_->particles_stored(), 3); }
}  // namespace
