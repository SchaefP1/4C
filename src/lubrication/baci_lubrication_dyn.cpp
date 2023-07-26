/*--------------------------------------------------------------------------*/
/*! \file

\brief entry point for the solution of Lubrication problems

\level 3


*/
/*--------------------------------------------------------------------------*/

#include "baci_lubrication_dyn.H"

#include "baci_adapter_lubrication.H"

#include "baci_lib_utils_createdis.H"

#include "baci_lubrication_timint_implicit.H"
#include "baci_lib_dofset_predefineddofnumber.H"


/*----------------------------------------------------------------------*
 | Main control routine for Lubrication problems            wirtz 11/15 |
 *----------------------------------------------------------------------*/
void lubrication_dyn(int restart)
{
  // access the communicator
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("lubrication")->Comm();

  // print problem type
  if (comm.MyPID() == 0)
  {
    std::cout << "###################################################" << std::endl;
    std::cout << "# YOUR PROBLEM TYPE: " << DRT::Problem::Instance()->ProblemName() << std::endl;
    std::cout << "###################################################" << std::endl;
  }

  // access the problem-specific parameter list
  const Teuchos::ParameterList& lubricationdyn =
      DRT::Problem::Instance()->LubricationDynamicParams();

  // access the lubrication discretization
  Teuchos::RCP<DRT::Discretization> lubricationdis =
      DRT::Problem::Instance()->GetDis("lubrication");

  lubricationdis->FillComplete();

  // we directly use the elements from the Lubrication elements section
  if (lubricationdis->NumGlobalNodes() == 0)
    dserror("No elements in the ---LUBRICATION ELEMENTS section");

  // add proxy of velocity related degrees of freedom to lubrication discretization
  Teuchos::RCP<DRT::DofSetInterface> dofsetaux = Teuchos::rcp(
      new DRT::DofSetPredefinedDoFNumber(DRT::Problem::Instance()->NDim(), 0, 0, true));
  if (lubricationdis->AddDofSet(dofsetaux) != 1)
    dserror("lub discretization has illegal number of dofsets!");

  // finalize discretization
  lubricationdis->FillComplete(true, false, false);

  // get linear solver id from LUBRICATION DYNAMIC
  const int linsolvernumber = lubricationdyn.get<int>("LINEAR_SOLVER");
  if (linsolvernumber == (-1))
    dserror(
        "no linear solver defined for LUBRICATION problem. Please set LINEAR_SOLVER in LUBRICATION "
        "DYNAMIC to a valid number!");

  // create instance of Lubrication basis algorithm
  Teuchos::RCP<ADAPTER::LubricationBaseAlgorithm> lubricationonly =
      Teuchos::rcp(new ADAPTER::LubricationBaseAlgorithm());

  // setup Lubrication basis algorithm
  lubricationonly->Setup(
      lubricationdyn, lubricationdyn, DRT::Problem::Instance()->SolverParams(linsolvernumber));

  // read the restart information, set vectors and variables
  if (restart) lubricationonly->LubricationField()->ReadRestart(restart);

  // enter time loop to solve problem
  (lubricationonly->LubricationField())->TimeLoop();

  // perform the result test if required
  DRT::Problem::Instance()->AddFieldTest(lubricationonly->CreateLubricationFieldTest());
  DRT::Problem::Instance()->TestAll(comm);

  return;

}  // end of lubrication_dyn()