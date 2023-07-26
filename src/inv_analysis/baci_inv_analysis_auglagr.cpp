/*----------------------------------------------------------------------*/
/*! \file
\brief Augmetend lagrangian functional control algorithm

\level 3

*/
/*----------------------------------------------------------------------*/
#include <string>
#include <iostream>
#include <fstream>
#include <cstdlib>

#include "baci_inv_analysis_auglagr.H"
#include "baci_inv_analysis_utils.H"

#include "baci_inv_analysis_matpar_manager.H"
#include "baci_inv_analysis_objective_funct.H"
#include "baci_inv_analysis_regularization_base.H"
#include "baci_inv_analysis_optimizer_base.H"

#include "baci_adapter_str_invana.H"
#include "baci_adapter_str_factory.H"
#include "baci_adapter_str_structure_new.H"
#include "baci_adapter_str_structure.H"


#include "baci_inv_analysis_timint_adjoint.H"
#include "baci_inv_analysis_timint_adjoint_prestress.H"

#include "baci_io_control.H"
#include "baci_io_pstream.H"
#include "baci_io_hdf.H"
#include "baci_io.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_prestress_service.H"
#include "baci_timestepping_mstep.H"
#include "baci_lib_parobject.H"
#include "baci_linalg_utils_sparse_algebra_create.H"
#include "baci_inpar_statinvanalysis.H"
#include "baci_inpar_validparameters.H"

#include <Teuchos_ParameterList.hpp>


/*----------------------------------------------------------------------*/
/* standard constructor                                      keh 09/14  */
/*----------------------------------------------------------------------*/
INVANA::InvanaAugLagr::InvanaAugLagr()
    : InvanaBase(),
      inputfile_(Teuchos::null),
      dis_(Teuchos::null),
      elementdata_(Teuchos::null),
      disdual_(Teuchos::null),
      disdualp_(Teuchos::null),
      timestep_(0.0),
      msteps_(0),
      pstype_(INPAR::STR::PreStress::none),
      pstime_(0.0),
      itertopc_(10)
{
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& invp = DRT::Problem::Instance()->StatInverseAnalysisParams();

  // this is supposed to be the number of simulation steps in the primal AND the dual problem
  msteps_ = sdyn.get<int>("NUMSTEP");
  timestep_ = sdyn.get<double>("TIMESTEP");

  // prestress stuff
  pstype_ = Teuchos::getIntegralValue<INPAR::STR::PreStress>(invp, "PRESTRESS");
  if (::UTILS::PRESTRESS::IsMulf(pstype_)) pstime_ = ::UTILS::PRESTRESS::GetPrestressTime();

  // initialize the vector of time steps according to the structural dynamic params
  time_.resize(msteps_, 0.0);
  for (int i = 0; i < msteps_; i++) time_[i] = (i + 1) * timestep_;

  itertopc_ = invp.get<int>("ITERTOPC");

  fpcounter_ = 0;
}

/*----------------------------------------------------------------------*/
/* Setup                                                     keh 09/14  */
/*----------------------------------------------------------------------*/
void INVANA::InvanaAugLagr::Setup()
{
  if (not Discret()->Filled() || not Discret()->HaveDofs())
    dserror("Discretisation is not complete or has no dofs!");

  // initialize "state" vectors
  dis_ = Teuchos::rcp(new Epetra_MultiVector(*(Discret()->DofRowMap()), msteps_, true));
  disdual_ = Teuchos::rcp(new Epetra_MultiVector(*(Discret()->DofRowMap()), msteps_, true));
  disdualp_ = Teuchos::rcp(new Epetra_MultiVector(*(Discret()->DofRowMap()), msteps_, true));

  return;
}

/*----------------------------------------------------------------------*/
/* MStep EpetraVector to EpetraMultiVector                   keh 10/13  */
/*----------------------------------------------------------------------*/
void INVANA::InvanaAugLagr::MStepEpetraToEpetraMulti(
    Teuchos::RCP<std::map<int, Epetra_Vector>> mstepvec, Teuchos::RCP<Epetra_MultiVector> multivec)
{
  multivec->Scale(0.0);
  std::map<int, Epetra_Vector>::iterator curr;
  for (curr = mstepvec->begin(); curr != mstepvec->end(); curr++)
    (*multivec)(curr->first - 1)->Update(1.0, curr->second, 0.0);
}

/*----------------------------------------------------------------------*/
/* Mstep double to std::vector<double>                       keh 10/13  */
/*----------------------------------------------------------------------*/
void INVANA::InvanaAugLagr::MStepDToStdVecD(
    Teuchos::RCP<std::map<int, double>> mstepvec, std::vector<double>* stdvec)
{
  stdvec->resize(msteps_, 0.0);
  std::map<int, double>::iterator curr;
  for (curr = mstepvec->begin(); curr != mstepvec->end(); curr++)
    (*stdvec)[curr->first - 1] = curr->second;
}

/*----------------------------------------------------------------------*/
/* solve primal problem                                      keh 10/13  */
/*----------------------------------------------------------------------*/
int INVANA::InvanaAugLagr::SolveForwardProblem()
{
  int err = 0;

  // use the same control file for every run since usually the last one is of interest
  Discret()->Writer()->OverwriteResultFile();

  // get input lists
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();

  // major switch to different time integrators
  switch (DRT::INPUT::IntegralValue<INPAR::STR::DynamicType>(sdyn, "DYNAMICTYP"))
  {
    case INPAR::STR::dyna_statics:
    {
      // create an adapterbase and adapter
      Teuchos::RCP<ADAPTER::StructureInvana> structadapter = Teuchos::null;
      // FixMe The following switch is just a temporal hack, such we can jump between the new and
      // the old structure implementation. Has to be deleted after the clean-up has been finished!
      const enum INPAR::STR::IntegrationStrategy intstrat =
          DRT::INPUT::IntegralValue<INPAR::STR::IntegrationStrategy>(sdyn, "INT_STRATEGY");
      switch (intstrat)
      {
        // -------------------------------------------------------------------
        // old implementation
        // -------------------------------------------------------------------
        case INPAR::STR::int_old:
        {
          ADAPTER::StructureBaseAlgorithm adapterbase_old(sdyn, sdyn, Discret());
          structadapter =
              Teuchos::rcp_dynamic_cast<ADAPTER::StructureInvana>(adapterbase_old.StructureField());
          structadapter->Setup();

          break;
        }
        // -------------------------------------------------------------------
        // new implementation
        // -------------------------------------------------------------------
        default:
        {
          Teuchos::RCP<ADAPTER::StructureBaseAlgorithmNew> adapterbase_ptr =
              ADAPTER::STR::BuildStructureAlgorithm(sdyn);
          adapterbase_ptr->Init(sdyn, const_cast<Teuchos::ParameterList&>(sdyn), Discret());
          adapterbase_ptr->Setup();
          structadapter = Teuchos::rcp_dynamic_cast<ADAPTER::StructureInvana>(
              adapterbase_ptr->StructureField(), true);

          break;
        }
      }

      // do restart but the one which is explicitly given in the INVERSE ANALYSIS section
      // and only if we are not in parameter continuation mode
      if (FPRestart() and fpcounter_ <= itertopc_)
      {
        DRT::Problem::Instance()->SetInputControlFile(InputControl());
        structadapter->ReadRestart(FPRestart());
      }

      if (fpcounter_ > itertopc_)
      {
        for (int i = 0; i < (int)time_.size(); i++)
        {
          // find whether measurements exist for this simulation timestep
          int mstep = ObjectiveFunct()->FindStep(time_[i]);

          // do this step only in case of measurements or when it is the last prestress step
          if (mstep != -1 or (time_[i] > pstime_ - structadapter->Dt() and time_[i] < pstime_))
          {
            structadapter->SetTimeStepStateOld(time_[i], i + 1, Teuchos::rcp((*dis_)(i), false),
                Teuchos::rcp((*dis_)(i), false));  // veln is not used so far so just put disn

            err = structadapter->Integrate();
          }
        }
      }
      else
      {
        err = structadapter->Integrate();
      }

      // get displacement and time
      MStepEpetraToEpetraMulti(structadapter->DispSteps(), dis_);
      MStepDToStdVecD(structadapter->TimeSteps(), &time_);

      break;
    }
    case INPAR::STR::dyna_genalpha:
    case INPAR::STR::dyna_onesteptheta:
    case INPAR::STR::dyna_gemm:
    case INPAR::STR::dyna_expleuler:
    case INPAR::STR::dyna_centrdiff:
    case INPAR::STR::dyna_ab2:
    case INPAR::STR::dyna_euma:
    case INPAR::STR::dyna_euimsto:
      dserror("inverse analysis for statics only so far");
      break;
    default:
      dserror("unknown time integration scheme '%s'", sdyn.get<std::string>("DYNAMICTYP").c_str());
      break;
  }

  fpcounter_ += 1;
  return err;
}


/*----------------------------------------------------------------------*/
/* solve dual problem                                       keh 10/13   */
/*----------------------------------------------------------------------*/
void INVANA::InvanaAugLagr::SolveAdjointProblem()
{
  // Setup RHS for the adjoints
  Teuchos::RCP<Epetra_Vector> objgrad =
      Teuchos::rcp(new Epetra_Vector(*(Discret()->DofRowMap()), true));

  std::vector<double> mtime = ObjectiveFunct()->MeasuredTime();
  Teuchos::RCP<Epetra_MultiVector> adjrhs =
      Teuchos::rcp(new Epetra_MultiVector(*(Discret()->DofRowMap()), (int)mtime.size(), true));

  for (int i = 0; i < (int)time_.size(); i++)
  {
    // find whether measurements exist for this simulation timestep
    int mstep = ObjectiveFunct()->FindStep(time_[i]);
    if (mstep != -1)
    {
      ObjectiveFunct()->EvaluateGradient(Teuchos::rcp((*dis_)(i), false), time_[i], objgrad);
      (*adjrhs)(mstep)->Scale(1.0, *objgrad);
    }
  }

  // initialize adjoint time integration with RHS as input
  Teuchos::RCP<STR::TimIntAdjoint> timintadj;
  if (::UTILS::PRESTRESS::IsNone(pstype_))
    timintadj = Teuchos::rcp(new STR::TimIntAdjoint(Discret()));
  else if (::UTILS::PRESTRESS::IsMulf(pstype_))
    timintadj = Teuchos::rcp(new STR::TimIntAdjointPrestress(Discret()));

  timintadj->SetupAdjoint(adjrhs, mtime, dis_, time_);

  // adjoint time integration
  timintadj->Integrate();

  // get the solution
  disdual_->Update(1.0, *timintadj->ExtractSolution(), 0.0);

  if (::UTILS::PRESTRESS::IsMulf(pstype_))
  {
    disdualp_->Update(1.0, *timintadj->ExtractPrestressSolution(), 0.0);
  }
}

/*----------------------------------------------------------------------*/
/* evaluate the value and/or gradient of the problem        keh 10/13   */
/*----------------------------------------------------------------------*/
int INVANA::InvanaAugLagr::Evaluate(
    const Epetra_MultiVector& sol, double* val, Teuchos::RCP<Epetra_MultiVector> gradient)
{
  int err = 0;

  // if (Optimizer()->Runc()<=itertopc_)
  Matman()->ReplaceParams(sol);

  if (fpcounter_ <= itertopc_) ResetDiscretization();

  if (gradient != Teuchos::null or val != NULL)
  {
    err = SolveForwardProblem();
    EvaluateError(sol, val);

    if (gradient != Teuchos::null)
    {
      SolveAdjointProblem();
      EvaluateGradient(sol, gradient);
    }
  }

  // do scaling
  double fac = ObjectiveFunct()->GetScaleFac();
  if (gradient != Teuchos::null) gradient->Scale(fac);
  if (val != NULL) *val *= fac;

  return err;
}

/*----------------------------------------------------------------------*/
/* evaluate gradient of the objective function              keh 10/13   */
/*----------------------------------------------------------------------*/
void INVANA::InvanaAugLagr::EvaluateGradient(
    const Epetra_MultiVector& sol, Teuchos::RCP<Epetra_MultiVector> gradient)
{
  // zero out gradient vector initially
  gradient->Scale(0.0);

  Teuchos::RCP<Epetra_Vector> zeros;
  zeros = CORE::LINALG::CreateVector(*(Discret()->DofRowMap()), true);

  // loop the time steps
  for (int j = 0; j < msteps_; j++)
  {
    int mstep = ObjectiveFunct()->FindStep(time_[j]);
    if (mstep == -1) continue;

    Discret()->SetState(0, "displacement", Teuchos::rcp((*dis_)(j), false));
    Discret()->SetState(0, "residual displacement", zeros);
    Discret()->SetState(0, "dual displacement", Teuchos::rcp((*disdual_)(j), false));
    SetTimeStepHistory(j + 1);

    Matman()->AddEvaluate(time_[j], gradient);

    Discret()->ClearState();

    if (::UTILS::PRESTRESS::IsMulf(pstype_))
    {
      int stepps = (int)(pstime_ / (timestep_));
      Discret()->SetState(0, "displacement", Teuchos::rcp((*dis_)(stepps - 1), false));
      Discret()->SetState(0, "residual displacement", zeros);
      Discret()->SetState(0, "dual displacement", Teuchos::rcp((*disdualp_)(j), false));
      SetTimeStepHistory(stepps);

      Matman()->AddEvaluate(time_[stepps - 1], gradient);
    }
  }

  if (Regman() != Teuchos::null) Regman()->EvaluateGradient(sol, gradient);
}

/*----------------------------------------------------------------------*/
/* evaluate gradient of the objective function using FD     keh 12/14   */
/*----------------------------------------------------------------------*/
void INVANA::InvanaAugLagr::EvaluateGradientFD(
    const Epetra_MultiVector& sol, Teuchos::RCP<Epetra_MultiVector> gradient)
{
  if (Discret()->Comm().NumProc() != 1)
    dserror("FD gradient evaluation is only implemented for single processor simulations");

  // zero out gradient vector initially
  gradient->Scale(0.0);

  // get a non const copy of the parameters
  Teuchos::RCP<Epetra_MultiVector> params_0 = Teuchos::rcp(new Epetra_MultiVector(sol));

  // Evaluate the reference solution:
  double val_0;
  Matman()->ReplaceParams(sol);
  ResetDiscretization();
  SolveForwardProblem();
  EvaluateError(*params_0, &val_0);

  // do the perturbation loop
  double alpha = 1.0e-7;
  double beta = 1.0e-7;
  double dp = 0.0;  // perturbation parameter
  double val_p;
  for (int j = 0; j < params_0->MyLength(); j++)
  {
    dp = (*(*params_0)(0))[j] * alpha + beta;
    params_0->SumIntoGlobalValue(j, 0, -dp);
    Matman()->ReplaceParams(*params_0);
    ResetDiscretization();
    SolveForwardProblem();
    EvaluateError(*params_0, &val_p);
    gradient->ReplaceGlobalValue(j, 0, (val_0 - val_p) / dp);
    params_0->Update(1.0, sol, 0.0);
  }
}

/*----------------------------------------------------------------------*/
/* Reset the discretization                                 keh 10/13   */
/*----------------------------------------------------------------------*/
void INVANA::InvanaAugLagr::ResetDiscretization()
{
  Teuchos::ParameterList p;
  p.set("action", "calc_struct_reset_all");
  Discret()->Evaluate(p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
}


/*----------------------------------------------------------------------*/
void INVANA::InvanaAugLagr::SetTimeStepHistory(int timestep)
{
  Teuchos::ParameterList p;
  p.set("timestep", timestep);
  p.set("action", "calc_struct_recover_istep");
  Discret()->Evaluate(p, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null, Teuchos::null);
}

/*----------------------------------------------------------------------*/
/* Evaluate the objective function                          keh 10/13   */
/*----------------------------------------------------------------------*/
void INVANA::InvanaAugLagr::EvaluateError(const Epetra_MultiVector& sol, double* val)
{
  *val = 0.0;

  double toadd;
  for (int i = 0; i < (int)time_.size(); i++)
  {
    // find whether measurements exist for this simulation timestep
    int mstep = ObjectiveFunct()->FindStep(time_[i]);
    if (mstep != -1)
    {
      ObjectiveFunct()->Evaluate(Teuchos::rcp((*dis_)(i), false), time_[i], toadd);
      *val += toadd;
    }
  }

  if (Regman() != Teuchos::null) Regman()->Evaluate(sol, val);
}