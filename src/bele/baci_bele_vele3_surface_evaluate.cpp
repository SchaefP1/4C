/*----------------------------------------------------------------------*/
/*! \file

\brief volume element


\level 2
*/
/*----------------------------------------------------------------------*/

#include "baci_bele_vele3.H"
#include "baci_lib_discret.H"
#include "baci_utils_exceptions.H"



/*----------------------------------------------------------------------*
 |  evaluate the element (public)                            gammi 04/07|
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Vele3Surface::Evaluate(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, std::vector<int>& lm, Epetra_SerialDenseMatrix& elemat1,
    Epetra_SerialDenseMatrix& elemat2, Epetra_SerialDenseVector& elevec1,
    Epetra_SerialDenseVector& elevec2, Epetra_SerialDenseVector& elevec3)
{
  return 0;
}


/*----------------------------------------------------------------------*
 |  do nothing (public)                                      u.may 05/09|
 |                                                                      |
 |  The function is just a dummy.                                       |
 *----------------------------------------------------------------------*/
int DRT::ELEMENTS::Vele3Surface::EvaluateNeumann(Teuchos::ParameterList& params,
    DRT::Discretization& discretization, DRT::Condition& condition, std::vector<int>& lm,
    Epetra_SerialDenseVector& elevec1, Epetra_SerialDenseMatrix* elemat1)
{
  return 0;
}