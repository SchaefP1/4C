/*----------------------------------------------------------------------*/
/*! \file

\brief Evaluate boundary conditions for scalar transport problems

\level 2

 *----------------------------------------------------------------------*/
#include "4C_lib_discret.hpp"
#include "4C_mat_elchmat.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_ele_action.hpp"
#include "4C_scatra_ele_boundary_calc.hpp"
#include "4C_scatra_ele_boundary_factory.hpp"
#include "4C_scatra_ele_parameter_elch.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TransportBoundary::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  FOUR_C_THROW("not implemented. Use the Evaluate() method with Location Array instead!");
  return -1;
}

/*----------------------------------------------------------------------*
 |  evaluate the element (public)                             gjb 01/09 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TransportBoundary::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, LocationArray& la,
    CORE::LINALG::SerialDenseMatrix& elemat1, CORE::LINALG::SerialDenseMatrix& elemat2,
    CORE::LINALG::SerialDenseVector& elevec1, CORE::LINALG::SerialDenseVector& elevec2,
    CORE::LINALG::SerialDenseVector& elevec3)
{
  // we assume here, that numdofpernode is equal for every node within
  // the discretization and does not change during the computations
  const int numdofpernode = NumDofPerNode(*(Nodes()[0]));
  int numscal = numdofpernode;

  // perform additional operations specific to implementation type
  switch (parent_element()->ImplType())
  {
    case INPAR::SCATRA::impltype_elch_diffcond:
    case INPAR::SCATRA::impltype_elch_diffcond_thermo:
    case INPAR::SCATRA::impltype_elch_electrode:
    case INPAR::SCATRA::impltype_elch_electrode_growth:
    case INPAR::SCATRA::impltype_elch_electrode_thermo:
    case INPAR::SCATRA::impltype_elch_NP:
    {
      // adapt number of transported scalars for electrochemistry problems
      numscal -= 1;

      // get the material of the first element
      // we assume here, that the material is equal for all elements in this discretization
      // get the parent element including its material
      Teuchos::RCP<CORE::MAT::Material> material = parent_element()->Material();
      if (material->MaterialType() == CORE::Materials::m_elchmat)
        numscal = static_cast<const MAT::ElchMat*>(material.get())->NumScal();

      break;
    }

    case INPAR::SCATRA::impltype_std:
    case INPAR::SCATRA::impltype_advreac:
    case INPAR::SCATRA::impltype_refconcreac:
    case INPAR::SCATRA::impltype_chemo:
    case INPAR::SCATRA::impltype_chemoreac:
    case INPAR::SCATRA::impltype_aniso:
    case INPAR::SCATRA::impltype_cardiac_monodomain:
    case INPAR::SCATRA::impltype_levelset:
    case INPAR::SCATRA::impltype_loma:
    case INPAR::SCATRA::impltype_poro:
    case INPAR::SCATRA::impltype_pororeac:
    case INPAR::SCATRA::impltype_thermo_elch_diffcond:
    case INPAR::SCATRA::impltype_thermo_elch_electrode:
    case INPAR::SCATRA::impltype_multipororeac:
      // do nothing in these cases
      break;

    default:
    {
      // other implementation types are invalid
      FOUR_C_THROW("Invalid implementation type!");
      break;
    }
  }

  // all physics-related stuff is included in the implementation class that can
  // be used in principle inside any element (at the moment: only Transport
  // boundary element)
  // If this element has special features/ methods that do not fit in the
  // generalized implementation class, you have to do a switch here in order to
  // call element-specific routines
  return DRT::ELEMENTS::ScaTraBoundaryFactory::ProvideImpl(
      this, parent_element()->ImplType(), numdofpernode, numscal, discretization.Name())
      ->Evaluate(this, params, discretization, la, elemat1, elemat2, elevec1, elevec2, elevec3);
}


/*----------------------------------------------------------------------*
 | evaluate Neumann boundary condition on boundary element   fang 01/15 |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::TransportBoundary::evaluate_neumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, CORE::Conditions::Condition& condition,
    std::vector<int>& lm, CORE::LINALG::SerialDenseVector& elevec1,
    CORE::LINALG::SerialDenseMatrix* elemat1)
{
  // add Neumann boundary condition to parameter list
  params.set<CORE::Conditions::Condition*>("condition", &condition);

  LocationArray la(discretization.NumDofSets());
  DRT::Element::LocationVector(discretization, la, false);

  // evaluate boundary element
  return Evaluate(params, discretization, la, *elemat1, *elemat1, elevec1, elevec1, elevec1);
}


/*----------------------------------------------------------------------*
 |  Get degrees of freedom used by this element                (public) |
 *----------------------------------------------------------------------*/
void DRT::ELEMENTS::TransportBoundary::LocationVector(const Discretization& dis, LocationArray& la,
    bool doDirichlet, const std::string& condstring, Teuchos::ParameterList& params) const
{
  // check for the action parameter
  const auto action = Teuchos::getIntegralValue<SCATRA::BoundaryAction>(params, "action");
  switch (action)
  {
    case SCATRA::BoundaryAction::calc_weak_Dirichlet:
      // special cases: the boundary element assembles also into
      // the inner dofs of its parent element
      // note: using these actions, the element will get the parent location vector
      //       as input in the respective evaluate routines
      parent_element()->LocationVector(dis, la, doDirichlet);
      break;
    default:
      DRT::Element::LocationVector(dis, la, doDirichlet);
      break;
  }
}

FOUR_C_NAMESPACE_CLOSE