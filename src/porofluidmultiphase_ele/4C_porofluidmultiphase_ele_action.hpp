/*----------------------------------------------------------------------*/
/*! \file
 \brief element actions for evaluation of porofluidmultiphase element

   \level 3

 *----------------------------------------------------------------------*/

#ifndef FOUR_C_POROFLUIDMULTIPHASE_ELE_ACTION_HPP
#define FOUR_C_POROFLUIDMULTIPHASE_ELE_ACTION_HPP

#include "4C_config.hpp"

FOUR_C_NAMESPACE_OPEN

namespace POROFLUIDMULTIPHASE
{
  /*--------------------------------------------------------------------------*/
  /*!
   * \brief enum that provides all possible porofluidmultiphase actions
   *///                                                            vuong 08/16
  /*--------------------------------------------------------------------------*/
  enum Action
  {
    // domain action
    set_general_parameter,        // set general parameters for element evaluation
    set_timestep_parameter,       // set time-integration parameters for element evaluation
    calc_mat_and_rhs,             // calc system matrix and residual,
    calc_fluid_struct_coupl_mat,  // calculate off-diagonal fluid-structure coupling matrix
    calc_fluid_scatra_coupl_mat,  // calculate off-diagonal fluid-structure coupling matrix
    calc_initial_time_deriv,      // calculate initial time derivative
    calc_pres_and_sat,            // calculate pressure and saturation
    calc_solidpressure,           // calculate solid pressure
    calc_porosity,                // calculate porosity
    calc_domain_integrals,        // calculate domain integrals
    calc_error,                   // calc_error
    recon_flux_at_nodes,          // reconstruct flux at nodes
    calc_phase_velocities,        // calculate phase velocities
    get_access_from_scatra,       // get access from scatra-framework to evaluate variables
    get_access_from_artcoupling,  // get access from artery-coupling to evaluate variables
    calc_valid_dofs               // calculate the valid volume fraction pressure and species dofs
  };                              // enum Action

  /*--------------------------------------------------------------------------*/
  /*!
   * \brief enum that provides all possible porofluidmultiphase actions on a boundary
   *///                                                            vuong 08/16
  /*--------------------------------------------------------------------------*/
  enum BoundaryAction
  {
    bd_calc_Neumann,  // evaluate neumann loads
  };                  // enum POROFLUIDMULTIPHASE::BoundaryAction
}  // namespace POROFLUIDMULTIPHASE

FOUR_C_NAMESPACE_CLOSE

#endif