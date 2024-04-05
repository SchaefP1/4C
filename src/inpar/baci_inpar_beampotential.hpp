/*----------------------------------------------------------------------*/
/*! \file

\brief input parameter definitions for beam potential-based interactions

\level 3

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_INPAR_BEAMPOTENTIAL_HPP
#define BACI_INPAR_BEAMPOTENTIAL_HPP

#include "baci_config.hpp"

#include "baci_utils_parameter_list.hpp"

BACI_NAMESPACE_OPEN

// forward declaration
namespace INPUT
{
  class ConditionDefinition;
}

/*----------------------------------------------------------------------*/
namespace INPAR
{
  namespace BEAMPOTENTIAL
  {
    /// type of potential interaction
    /// (this enum represents the input file parameter BEAMPOTENTIAL_TYPE)
    enum BeamPotentialType
    {
      beampot_surf,  ///< surface potential
      beampot_vol,   ///< volume potential
      beampot_vague
    };

    /// available strategies/methods to evaluate potential interaction
    /// (this enum represents the input file parameter STRATEGY)
    enum BeamPotentialStrategy
    {
      strategy_doublelengthspec_largesepapprox,         ///< double length specific potential, large
                                                        ///< separations
      strategy_doublelengthspec_smallsepapprox,         ///< double length specific potential, small
                                                        ///< separations
      strategy_singlelengthspec_smallsepapprox,         ///< single length specific potential, small
                                                        ///< separations
      strategy_singlelengthspec_smallsepapprox_simple,  ///< reduced variant of the previous one
      strategy_vague
    };

    /// available types to regularize the force law for separations smaller than
    /// the specified regularization parameter
    enum BeamPotentialRegularizationType
    {
      regularization_linear,    ///< linear extrapolation
      regularization_constant,  ///< constant extrapolation, i.e. f(r)=f(r_reg) for all r<r_reg
      regularization_none       ///< no regularization
    };

    /**
     * \brief rule for how to assign the role of slave and master to beam elements
     */
    enum class MasterSlaveChoice
    {
      smaller_eleGID_is_slave,
      higher_eleGID_is_slave,
      choice_master_slave_vague
    };

    /// set the beam potential parameters
    void SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list);

    /// set beam potential specific conditions
    void SetValidConditions(std::vector<Teuchos::RCP<INPUT::ConditionDefinition>>& condlist);

  }  // namespace BEAMPOTENTIAL

}  // namespace INPAR

/*----------------------------------------------------------------------*/
BACI_NAMESPACE_CLOSE

#endif  // INPAR_BEAMPOTENTIAL_H