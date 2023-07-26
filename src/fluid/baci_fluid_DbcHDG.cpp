/*-----------------------------------------------------------*/
/*! \file

\brief contains utils functions for Dirichlet Boundary Conditions of HDG discretizations

\level 2

*/

#include "baci_elemag_ele_action.H"

#include "baci_fluid_DbcHDG.h"
#include "baci_fluid_ele_action.H"
#include "baci_fluid_ele_calc.H"
#include "baci_fluid_ele_calc_hdg.H"
#include "baci_fluid_ele_hdg.H"

#include "baci_lib_discret_hdg.H"
#include "baci_lib_globalproblem.H"
#include "baci_lib_utils_parameter_list.H"

#include "baci_scatra_ele_action.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::UTILS::DbcHDG_Fluid::ReadDirichletCondition(const DRT::DiscretizationInterface& discret,
    const DRT::Condition& cond, double time, DRT::UTILS::Dbc::DbcInfo& info,
    const Teuchos::RCP<std::set<int>>* dbcgids, int hierarchical_order) const
{
  // no need to check the cast, because it has been done during
  // the build process (see BuildDbc())
  const DRT::DiscretizationFaces& face_discret =
      static_cast<const DRT::DiscretizationFaces&>(discret);

  ReadDirichletCondition(face_discret, cond, time, info, dbcgids, hierarchical_order);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::UTILS::DbcHDG_Fluid::ReadDirichletCondition(const DRT::DiscretizationFaces& discret,
    const DRT::Condition& cond, double time, DRT::UTILS::Dbc::DbcInfo& info,
    const Teuchos::RCP<std::set<int>>* dbcgids, int hierarchical_order) const

{
  // call to corresponding method in base class; safety checks inside
  DRT::UTILS::Dbc::ReadDirichletCondition(discret, cond, time, info, dbcgids, hierarchical_order);

  // say good bye if there are no face elements
  if (discret.FaceRowMap() == NULL) return;

  // get onoff toggles
  const std::vector<int>* onoff = cond.Get<std::vector<int>>("onoff");

  if (discret.NumMyRowFaces() > 0)
  {
    // initialize with true on each proc except proc 0
    bool pressureDone = discret.Comm().MyPID() != 0;

    // loop over all faces
    for (int i = 0; i < discret.NumMyRowFaces(); ++i)
    {
      const DRT::FaceElement* faceele = dynamic_cast<const DRT::FaceElement*>(discret.lRowFace(i));
      const unsigned int dofperface =
          faceele->ParentMasterElement()->NumDofPerFace(faceele->FaceMasterNumber());
      const unsigned int dofpercomponent =
          faceele->ParentMasterElement()->NumDofPerComponent(faceele->FaceMasterNumber());
      const unsigned int component = dofperface / dofpercomponent;

      if (onoff->size() <= component || (*onoff)[component] == 0 ||
          DRT::Problem::Instance(0)->GetProblemType() != ProblemType::fluid)
        pressureDone = true;
      if (!pressureDone)
      {
        if (discret.NumMyRowElements() > 0 && discret.Comm().MyPID() == 0)
        {
          std::vector<int> predof = discret.Dof(0, discret.lRowElement(0));
          const int gid = predof[0];
          const int lid = discret.DofRowMap(0)->LID(gid);

          // set toggle vector
          info.toggle[lid] = 1;
          // amend vector of DOF-IDs which are Dirichlet BCs
          if (dbcgids[set_row] != Teuchos::null) (*dbcgids[set_row]).insert(gid);
          pressureDone = true;
        }
      }

      // do only faces where all nodes are present in the node list
      bool faceRelevant = true;
      int nummynodes = discret.lRowFace(i)->NumNode();
      const int* mynodes = discret.lRowFace(i)->NodeIds();
      for (int j = 0; j < nummynodes; ++j)
        if (!cond.ContainsNode(mynodes[j]))
        {
          faceRelevant = false;
          break;
        }
      if (!faceRelevant) continue;

      // get dofs of current face element
      std::vector<int> dofs = discret.Dof(0, discret.lRowFace(i));

      // loop over dofs
      for (unsigned int j = 0; j < dofperface; ++j)
      {
        // get global id
        const int gid = dofs[j];
        // get corresponding local id
        const int lid = info.toggle.Map().LID(gid);
        if (lid < 0)
          dserror(
              "Global id %d not on this proc %d in system vector", dofs[j], discret.Comm().MyPID());
        // get position of label for this dof in condition line
        int onesetj = j / dofpercomponent;

        if ((*onoff)[onesetj] == 0)
        {
          // no DBC on this dof, set toggle zero
          info.toggle[lid] = 0;
          // get rid of entry in DBC map - if it exists
          if (dbcgids[set_row] != Teuchos::null) (*dbcgids[set_row]).erase(gid);
          continue;
        }
        else  // if ((*onoff)[onesetj]==1)
        {
          // dof has DBC, set toggle vector one
          info.toggle[lid] = 1;
          // amend vector of DOF-IDs which are dirichlet BCs
          if (dbcgids[set_row] != Teuchos::null) (*dbcgids[set_row]).insert(gid);
        }

      }  // loop over DOFs of face
    }    // loop over all faces
  }      // if there are faces

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::UTILS::DbcHDG_Fluid::DoDirichletCondition(const DRT::DiscretizationInterface& discret,
    const DRT::Condition& cond, double time, const Teuchos::RCP<Epetra_Vector>* systemvectors,
    const Epetra_IntVector& toggle, const Teuchos::RCP<std::set<int>>* dbcgids) const
{
  // no need to check the cast, because it has been done during
  // the build process (see BuildDbc())
  const DRT::DiscretizationFaces& face_discret =
      static_cast<const DRT::DiscretizationFaces&>(discret);

  DoDirichletCondition(face_discret, cond, time, systemvectors, toggle);
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void FLD::UTILS::DbcHDG_Fluid::DoDirichletCondition(const DRT::DiscretizationFaces& discret,
    const DRT::Condition& cond, double time, const Teuchos::RCP<Epetra_Vector>* systemvectors,
    const Epetra_IntVector& toggle) const
{
  // call corresponding method from base class; safety checks inside
  DRT::UTILS::Dbc::DoDirichletCondition(discret, cond, time, systemvectors, toggle, NULL);

  // say good bye if there are no face elements
  if (discret.FaceRowMap() == NULL) return;

  // get ids of conditioned nodes
  const std::vector<int>* nodeids = cond.Nodes();
  if (!nodeids) dserror("Dirichlet condition does not have nodal cloud");

  // get curves, functs, vals, and onoff toggles from the condition
  const std::vector<int>* funct = cond.Get<std::vector<int>>("funct");
  const std::vector<double>* val = cond.Get<std::vector<double>>("val");
  const std::vector<int>* onoff = cond.Get<std::vector<int>>("onoff");

  // determine highest degree of time derivative
  // and first existent system vector to apply DBC to
  unsigned deg = 0;  // highest degree of requested time derivative
  Teuchos::RCP<Epetra_Vector> systemvectoraux = Teuchos::null;  // auxiliar system vector
  if (systemvectors[0] != Teuchos::null)
  {
    deg = 0;
    systemvectoraux = systemvectors[0];
  }
  if (systemvectors[1] != Teuchos::null)
  {
    deg = 1;
    if (systemvectoraux == Teuchos::null) systemvectoraux = systemvectors[1];
  }
  if (systemvectors[2] != Teuchos::null)
  {
    deg = 2;
    if (systemvectoraux == Teuchos::null) systemvectoraux = systemvectors[2];
  }

  // do we have faces?
  if (discret.NumMyRowFaces() > 0)
  {
    Epetra_SerialDenseVector elevec1, elevec2, elevec3;
    Epetra_SerialDenseMatrix elemat1, elemat2;
    DRT::Element::LocationArray dummy(1);
    Teuchos::ParameterList initParams;
    if (DRT::Problem::Instance(0)->GetProblemType() == ProblemType::elemag)
      initParams.set<int>("action", ELEMAG::project_dirich_field);
    else if (DRT::Problem::Instance(0)->GetProblemType() == ProblemType::scatra)
      DRT::UTILS::AddEnumClassToParameterList<SCATRA::Action>(
          "action", SCATRA::Action::project_dirich_field, initParams);
    else
      initParams.set<int>(
          "action", FLD::project_fluid_field);  // TODO: Introduce a general action type that is
                                                // valid for all problems
    if (funct != NULL)
    {
      Teuchos::Array<int> functarray(*funct);
      initParams.set("funct", functarray);
    }
    Teuchos::Array<int> onoffarray(*onoff);
    initParams.set("onoff", onoffarray);
    initParams.set("time", time);

    // initialize with true if proc is not proc 0
    bool pressureDone = discret.Comm().MyPID() != 0;

    // loop over all faces
    for (int i = 0; i < discret.NumMyRowFaces(); ++i)
    {
      const DRT::FaceElement* faceele = dynamic_cast<const DRT::FaceElement*>(discret.lRowFace(i));
      const unsigned int dofperface =
          faceele->ParentMasterElement()->NumDofPerFace(faceele->FaceMasterNumber());
      const unsigned int dofpercomponent =
          faceele->ParentMasterElement()->NumDofPerComponent(faceele->FaceMasterNumber());
      const unsigned int component = dofperface / dofpercomponent;

      if (onoff->size() <= component || (*onoff)[component] == 0 ||
          DRT::Problem::Instance(0)->GetProblemType() != ProblemType::fluid)
        pressureDone = true;
      if (!pressureDone)
      {
        if (discret.NumMyRowElements() > 0 && discret.Comm().MyPID() == 0)
        {
          std::vector<int> predof = discret.Dof(0, discret.lRowElement(0));
          const int gid = predof[0];
          const int lid = discret.DofRowMap(0)->LID(gid);

          // amend vector of DOF-IDs which are Dirichlet BCs
          if (systemvectors[0] != Teuchos::null) (*systemvectors[0])[lid] = 0.0;
          if (systemvectors[1] != Teuchos::null) (*systemvectors[1])[lid] = 0.0;
          if (systemvectors[2] != Teuchos::null) (*systemvectors[2])[lid] = 0.0;

          // --------------------------------------------------------------------------------------
          // get parameters
          Teuchos::ParameterList params = DRT::Problem::Instance()->FluidDynamicParams();

          // check whether the imposition of the average pressure is requested
          const int dopressavgbc =
              DRT::INPUT::IntegralValue<INPAR::FLUID::PressAvgBc>(params, "PRESSAVGBC");

          if (dopressavgbc == INPAR::FLUID::yes_pressure_average_bc)
          {
            double pressureavgBC = 0.0;

            // get 1st element
            DRT::Element* ele = discret.lRowElement(0);
            DRT::ELEMENTS::Fluid* fluidele = dynamic_cast<DRT::ELEMENTS::Fluid*>(ele);

            // get material
            Teuchos::RCP<MAT::Material> mat = ele->Material();

            // get discretization type
            const DRT::Element::DiscretizationType distype = ele->Shape();

            // evaluate pressure average     //TODO als make it valid for every discretization type
            Epetra_SerialDenseVector elevec = Epetra_SerialDenseVector(1);
            if (distype == DRT::Element::DiscretizationType::quad4)
              DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::DiscretizationType::quad4>::Instance()
                  ->EvaluatePressureAverage(fluidele, params, mat, elevec);
            else if (distype == DRT::Element::DiscretizationType::quad8)
              DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::DiscretizationType::quad8>::Instance()
                  ->EvaluatePressureAverage(fluidele, params, mat, elevec);
            else if (distype == DRT::Element::DiscretizationType::quad9)
              DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::DiscretizationType::quad9>::Instance()
                  ->EvaluatePressureAverage(fluidele, params, mat, elevec);
            else if (distype == DRT::Element::DiscretizationType::tri3)
              DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::DiscretizationType::tri3>::Instance()
                  ->EvaluatePressureAverage(fluidele, params, mat, elevec);
            else if (distype == DRT::Element::DiscretizationType::tri6)
              DRT::ELEMENTS::FluidEleCalcHDG<DRT::Element::DiscretizationType::tri6>::Instance()
                  ->EvaluatePressureAverage(fluidele, params, mat, elevec);
            else
              dserror("Given distype currently not implemented.");
            pressureavgBC = elevec[0];

            (*systemvectors[0])[lid] = pressureavgBC;

            std::cout << "\n-----------------------------------------------------------------------"
                         "-------------------"
                      << std::endl;
            std::cout << "| Warning: Imposing the analytical average pressure in the first element "
                         "as Dirichlet BC |"
                      << std::endl;
            std::cout << "-------------------------------------------------------------------------"
                         "-----------------\n"
                      << std::endl;
          }
          // --------------------------------------------------------------------------------------
          pressureDone = true;
        }
      }
      int nummynodes = discret.lRowFace(i)->NumNode();
      const int* mynodes = discret.lRowFace(i)->NodeIds();

      // do only faces where all nodes are present in the node list
      bool faceRelevant = true;
      for (int j = 0; j < nummynodes; ++j)
        if (!cond.ContainsNode(mynodes[j]))
        {
          faceRelevant = false;
          break;
        }
      if (!faceRelevant) continue;

      initParams.set<unsigned int>(
          "faceconsider", static_cast<unsigned int>(faceele->FaceMasterNumber()));
      if (static_cast<unsigned int>(elevec1.M()) != dofperface) elevec1.Shape(dofperface, 1);
      std::vector<int> dofs = discret.Dof(0, discret.lRowFace(i));

      bool do_evaluate = false;
      if (funct != NULL)
        for (unsigned int i = 0; i < component; ++i)
          if ((*funct)[i] > 0) do_evaluate = true;

      if (do_evaluate)
      {
        // cast the const qualifier away, thus the Evaluate routine can be called.
        DRT::DiscretizationFaces& non_const_dis = const_cast<DRT::DiscretizationFaces&>(discret);
        faceele->ParentMasterElement()->Evaluate(
            initParams, non_const_dis, dummy, elemat1, elemat2, elevec1, elevec2, elevec3);
      }
      else
        for (unsigned int i = 0; i < dofperface; ++i) elevec1(i) = 1.;

      // loop over face dofs
      for (unsigned int j = 0; j < dofperface; ++j)
      {
        // get global id
        const int gid = dofs[j];
        // get corresponding local id
        const int lid = toggle.Map().LID(gid);
        if (lid < 0)
          dserror(
              "Global id %d not on this proc %d in system vector", dofs[j], discret.Comm().MyPID());
        // get position of label for this dof in condition line
        int onesetj = j / dofpercomponent;

        // check whether dof gid is a dbc gid
        if (toggle[lid] == 0) continue;

        std::vector<double> value(deg + 1, (*val)[onesetj]);

        // assign value
        if (systemvectors[0] != Teuchos::null) (*systemvectors[0])[lid] = value[0] * elevec1(j);
        if (systemvectors[1] != Teuchos::null) (*systemvectors[1])[lid] = value[1] * elevec1(j);
        if (systemvectors[2] != Teuchos::null) (*systemvectors[2])[lid] = value[2] * elevec1(j);

      }  // loop over all DOFs
    }    // loop over all faces

  }  // if there are faces

  return;
}