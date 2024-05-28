/*-----------------------------------------------------------*/
/*! \file

\brief Base class for modelevaluators in partitioned algorithms.


\level 3

*/
/*-----------------------------------------------------------*/


#include "4C_structure_new_model_evaluator_multiphysics.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::MODELEVALUATOR::Multiphysics::Multiphysics() : active_mt_(mt_none)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Multiphysics::Init(
    const Teuchos::RCP<STR::MODELEVALUATOR::Data>& eval_data_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataGlobalState>& gstate_ptr,
    const Teuchos::RCP<STR::TIMINT::BaseDataIO>& gio_ptr,
    const Teuchos::RCP<STR::Integrator>& int_ptr,
    const Teuchos::RCP<const STR::TIMINT::Base>& timint_ptr, const int& dof_offset)
{
  STR::MODELEVALUATOR::Generic::Init(
      eval_data_ptr, gstate_ptr, gio_ptr, int_ptr, timint_ptr, dof_offset);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Multiphysics::Setup()
{
  check_init();
  issetup_ = true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Multiphysics::Reset(const Epetra_Vector& x)
{
  check_active_model_type();

  get_model_evaluator_from_map(active_mt_)->Reset(x);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Multiphysics::assemble_force(
    Epetra_Vector& f, const double& timefac_np) const
{
  check_active_model_type();

  get_model_evaluator_from_map(active_mt_)->assemble_force(f, timefac_np);
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Multiphysics::assemble_jacobian(
    CORE::LINALG::SparseOperator& jac, const double& timefac_np) const
{
  check_active_model_type();

  get_model_evaluator_from_map(active_mt_)->assemble_jacobian(jac, timefac_np);
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Multiphysics::evaluate_force()
{
  check_active_model_type();

  get_model_evaluator_from_map(active_mt_)->evaluate_force();
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Multiphysics::evaluate_stiff()
{
  check_active_model_type();

  get_model_evaluator_from_map(active_mt_)->evaluate_stiff();
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
bool STR::MODELEVALUATOR::Multiphysics::evaluate_force_stiff()
{
  check_active_model_type();

  get_model_evaluator_from_map(active_mt_)->evaluate_force_stiff();
  return true;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void STR::MODELEVALUATOR::Multiphysics::UpdateStepState(const double& timefac_n)
{
  check_active_model_type();

  get_model_evaluator_from_map(active_mt_)->UpdateStepState(timefac_n);
}

FOUR_C_NAMESPACE_CLOSE