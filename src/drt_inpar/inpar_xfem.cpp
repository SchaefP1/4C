/*----------------------------------------------------------------------*/
/*!
\file inpar_xfem.cpp

\brief Input parameters for XFEM

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
</pre>
*/

/*----------------------------------------------------------------------*/



#include "drt_validparameters.H"
#include "inpar_xfem.H"
#include "inpar_cut.H"
#include "../drt_lib/drt_conditiondefinition.H"



void INPAR::XFEM::SetValidParameters(Teuchos::RCP<Teuchos::ParameterList> list)
{
  using namespace DRT::INPUT;
  using Teuchos::tuple;
  using Teuchos::setStringToIntegralParameter;

  Teuchos::Array<std::string> yesnotuple = tuple<std::string>("Yes","No","yes","no","YES","NO");
  Teuchos::Array<int> yesnovalue = tuple<int>(true,false,true,false,true,false);

  Teuchos::ParameterList& xfem_general = list->sublist("XFEM GENERAL",false,"");

  // OUTPUT options
  setStringToIntegralParameter<int>("GMSH_DEBUG_OUT","No","Do you want to write extended Gmsh output for each timestep?",
                               yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("GMSH_DEBUG_OUT_SCREEN","No","Do you want to be informed, if Gmsh output is written?",
                                 yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("GMSH_SOL_OUT","No","Do you want to write extended Gmsh output for each timestep?",
                               yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("GMSH_EOS_OUT","No","Do you want to write extended Gmsh output for each timestep?",
                               yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("GMSH_DISCRET_OUT","No","Do you want to write extended Gmsh output for each timestep?",
                               yesnotuple,yesnovalue,&xfem_general);
  setStringToIntegralParameter<int>("GMSH_CUT_OUT","No","Do you want to write extended Gmsh output for each timestep?",
                               yesnotuple,yesnovalue,&xfem_general);


  IntParameter("MAX_NUM_DOFSETS",3,"Maximum number of volumecells in the XFEM element",&xfem_general);


  setStringToIntegralParameter<int>("NODAL_DOFSET_STRATEGY","full","which strategy used for the nodal dofset management per node?",
                                 tuple<std::string>("OneDofset_PerNodeAndPosition","full"),
                                 tuple<int>(
                                     INPAR::CUT::NDS_Strategy_OneDofset_PerNodeAndPosition,
                                     INPAR::CUT::NDS_Strategy_full
                                     ),
                                 &xfem_general);


  // Integration options
  setStringToIntegralParameter<int>("VOLUME_GAUSS_POINTS_BY","Tessellation","how to find Gauss Points for the cut volumes",
                               tuple<std::string>("Tessellation","MomentFitting","DirectDivergence"),
                               tuple<int>(
                                   INPAR::CUT::VCellGaussPts_Tessellation,
                                   INPAR::CUT::VCellGaussPts_MomentFitting,
                                   INPAR::CUT::VCellGaussPts_DirectDivergence
                                   ),
                               &xfem_general);

  setStringToIntegralParameter<int>("BOUNDARY_GAUSS_POINTS_BY","Tessellation","how to find Gauss Points for the boundary cells",
                                 tuple<std::string>("Tessellation","MomentFitting","DirectDivergence"),
                                 tuple<int>(
                                     INPAR::CUT::BCellGaussPts_Tessellation,
                                     INPAR::CUT::BCellGaussPts_MomentFitting,
                                     INPAR::CUT::BCellGaussPts_DirectDivergence
                                     ),
                                 &xfem_general);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfluid_dyn = list->sublist("XFLUID DYNAMIC",false,"");

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfluid_general = xfluid_dyn.sublist("GENERAL",false,"");

  // Do we use more than one fluid discretization?
  BoolParameter("XFLUIDFLUID","no","Use an embedded fluid patch.", &xfluid_general);

  // How many monolithic steps we keep the fluidfluid-boundary fixed
  IntParameter("RELAXING_ALE_EVERY",1,"Relaxing Ale after how many monolithic steps",&xfluid_general);

  BoolParameter("RELAXING_ALE","yes","switch on/off for relaxing Ale in monolithic fluid-fluid-fsi",&xfluid_general);

  DoubleParameter("XFLUIDFLUID_SEARCHRADIUS",1.0,"Radius of the search tree",&xfluid_general);

  // xfluidfluid-fsi-monolithic approach
  setStringToIntegralParameter<int>("MONOLITHIC_XFFSI_APPROACH","xffsi_fixedALE_partitioned","The monolithic apporach for xfluidfluid-fsi",
                                    tuple<std::string>("xffsi_full_newton", "xffsi_fixedALE_interpolation", "xffsi_fixedALE_partitioned"),
                                    tuple<int>(
                                      INPAR::XFEM::XFFSI_Full_Newton,    //xffsi with no fixed xfem-coupling
                                      INPAR::XFEM::XFFSI_FixedALE_Interpolation,  // xffsi with fixed xfem-coupling in every newtonstep
                                                                                  // and interpolations for embedded-dis afterwards
                                      INPAR::XFEM::XFFSI_FixedALE_Partitioned      // xffsi with fixed xfem-coupling in every newtonstep
                                                                                  // and solving fluid-field again
                                      ),
                                    &xfluid_general);

  // xfluidfluid time integration approach
  setStringToIntegralParameter<int>("XFLUIDFLUID_TIMEINT","Xff_TimeInt_FullProj","The xfluidfluid-timeintegration approach",
                                    tuple<std::string>("Xff_TimeInt_FullProj", "Xff_TimeInt_ProjIfMoved","Xff_TimeInt_KeepGhostValues","Xff_TimeInt_IncompProj"),
                                    tuple<int>(
                                      INPAR::XFEM::Xff_TimeInt_FullProj   ,      //always project nodes from embedded to background nodes
                                      INPAR::XFEM::Xff_TimeInt_ProjIfMoved,      //project nodes just if the status of background nodes changed
                                      INPAR::XFEM::Xff_TimeInt_KeepGhostValues,  //always keep the ghost values of the background discretization
                                      INPAR::XFEM::Xff_TimeInt_IncompProj        //after projecting nodes do a incompressibility projection
                                      ),
                                    &xfluid_general);

  setStringToIntegralParameter<int>("XFLUID_TIMEINT","STD=COPY/SL_and_GHOST=COPY/GP","The xfluid time integration approach",
                               tuple<std::string>("STD=COPY_and_GHOST=COPY/GP", "STD=COPY/SL_and_GHOST=COPY/GP", "STD=SL(boundary-zone)_and_GHOST=GP", "STD=COPY/PROJ_and_GHOST=COPY/PROJ/GP"),
                               tuple<int>(
                                   INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_AND_GHOST_by_Copy_or_GP,       // STD= only copy, GHOST= copy or ghost penalty reconstruction
                                   INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_SL_AND_GHOST_by_Copy_or_GP, // STD= copy or SL, GHOST= copy or ghost penalty reconstruction
                                   INPAR::XFEM::Xf_TimeIntScheme_STD_by_SL_cut_zone_AND_GHOST_by_GP,        // STD= only SL on whole boundary zone, GHOST= ghost penalty reconstruction
                                   INPAR::XFEM::Xf_TimeIntScheme_STD_by_Copy_or_Proj_AND_GHOST_by_Proj_or_Copy_or_GP
                                   ),
                               &xfluid_general);

  BoolParameter("ALE_XFluid","no","XFluid is Ale Fluid?",&xfluid_general);

  // for new OST-implementation: which interface terms to be evaluated for previous time step
  setStringToIntegralParameter<int>("INTERFACE_TERMS_PREVIOUS_STATE","PreviousState_only_consistency","how to treat interface terms from previous time step (new OST)",
                               tuple<std::string>("PreviousState_only_consistency", "PreviousState_full"),
                               tuple<int>(
                                   INPAR::XFEM::PreviousState_only_consistency, /// evaluate only consistency terms for previous time step
                                   INPAR::XFEM::PreviousState_full             /// evaluate consistency, adjoint consistency and penalty terms or previous time step
                                   ),
                               &xfluid_general);

  /*----------------------------------------------------------------------*/
  Teuchos::ParameterList& xfluid_stab = xfluid_dyn.sublist("STABILIZATION",false,"");

  // Boundary-Coupling options
  setStringToIntegralParameter<int>("COUPLING_METHOD","Nitsche","method how to enforce embedded boundary/coupling conditions at the interface",
                               tuple<std::string>("Hybrid_LM_Cauchy_stress", "Hybrid_LM_viscous_stress", "Nitsche"),
                               tuple<int>(
                                   INPAR::XFEM::Hybrid_LM_Cauchy_stress,  // Cauchy stress-based mixed/hybrid formulation
                                   INPAR::XFEM::Hybrid_LM_viscous_stress, // viscous stress-based mixed/hybrid formulation
                                   INPAR::XFEM::Nitsche                   // Nitsche's formulation
                                   ),
                               &xfluid_stab);

  setStringToIntegralParameter<int>("HYBRID_LM_L2_PROJ","part_ele_proj","perform the L2 projection between stress fields on whole element or on fluid part?",
                               tuple<std::string>("full_ele_proj", "part_ele_proj"),
                               tuple<int>(
                                   INPAR::XFEM::Hybrid_LM_L2_Proj_full,   // L2 stress projection on whole fluid element
                                   INPAR::XFEM::Hybrid_LM_L2_Proj_part    // L2 stress projection on partial fluid element volume
                                   ),
                               &xfluid_stab);

  BoolParameter("VISC_ADJOINT_SYMMETRY","yes","viscous and adjoint viscous interface terms with matching sign?",&xfluid_stab);

  // viscous and convective Nitsche/MSH stabilization parameter
  DoubleParameter("NIT_STAB_FAC", 35.0, " ( stabilization parameter for Nitsche's penalty term",&xfluid_stab);

  setStringToIntegralParameter<int>("VISC_STAB_TRACE_ESTIMATE","CT_div_by_hk","how to estimate the scaling from the trace inequality in Nitsche's method",
                               tuple<std::string>("CT_div_by_hk", "eigenvalue"),
                               tuple<int>(
                                   INPAR::XFEM::ViscStab_TraceEstimate_CT_div_by_hk,   // estimate the trace inequality by a trace-constant CT and hk: CT/hk
                                   INPAR::XFEM::ViscStab_TraceEstimate_eigenvalue      // estimate the trace inequality by solving a local element-wise eigenvalue problem
                                   ),
                               &xfluid_stab);

  setStringToIntegralParameter<int>("VISC_STAB_HK","ele_vol_div_by_max_ele_surf","how to define the characteristic element length in cut elements",
                                 tuple<std::string>(
                                     "vol_equivalent",
                                     "cut_vol_div_by_cut_surf",
                                     "ele_vol_div_by_cut_surf",
                                     "ele_vol_div_by_ele_surf",
                                     "ele_vol_div_by_max_ele_surf"
                                     ),
                                 tuple<int>(
                                     INPAR::XFEM::ViscStab_hk_vol_equivalent,             /// volume equivalent element diameter
                                     INPAR::XFEM::ViscStab_hk_cut_vol_div_by_cut_surf,    /// physical partial/cut volume divided by physical partial/cut surface measure ( used to estimate the cut-dependent inverse estimate on cut elements, not useful for sliver and/or dotted cut situations)
                                     INPAR::XFEM::ViscStab_hk_ele_vol_div_by_cut_surf,    /// full element volume divided by physical partial/cut surface measure ( used to estimate the cut-dependent inverse estimate on cut elements, however, avoids problems with sliver cuts, not useful for dotted cuts)
                                     INPAR::XFEM::ViscStab_hk_ele_vol_div_by_ele_surf,    /// full element volume divided by surface measure ( used for uncut situations, standard weak Dirichlet boundary/coupling conditions)
                                     INPAR::XFEM::ViscStab_hk_ele_vol_div_by_max_ele_surf /// default: full element volume divided by maximal element surface measure ( used to estimate the trace inequality for stretched elements in combination with ghost-penalties)
                                 ),
                                 &xfluid_stab);


  setStringToIntegralParameter<int>("CONV_STAB_SCALING","none","scaling factor for viscous interface stabilization (Nitsche, MSH)",
                                    tuple<std::string>("inflow", "abs_inflow", "none"),
                                    tuple<int>(
                                      INPAR::XFEM::ConvStabScaling_inflow,                // scaling with max(0,-u*n)
                                      INPAR::XFEM::ConvStabScaling_abs_inflow,            // scaling with |u*n|
                                      INPAR::XFEM::ConvStabScaling_none                   // no convective stabilization
                                      ),
                               &xfluid_stab);

  setStringToIntegralParameter<int>("XFF_CONV_STAB_SCALING","none","scaling factor for convective interface stabilization of fluid-fluid Coupling",
                                    tuple<std::string>("inflow", "averaged", "none"),
                                    tuple<int>(
                                      INPAR::XFEM::XFF_ConvStabScaling_upwinding,          // one-sided inflow stabilization
                                      INPAR::XFEM::XFF_ConvStabScaling_only_averaged,      // averaged inflow stabilization
                                      INPAR::XFEM::XFF_ConvStabScaling_none                // no convective stabilization
                                      ),
                               &xfluid_stab);

  setStringToIntegralParameter<int>("MASS_CONSERVATION_COMBO","max","choose the maximum from viscous and convective contributions or just sum both up",
                                    tuple<std::string>("max", "sum"),
                                    tuple<int>(
                                      INPAR::XFEM::MassConservationCombination_max,        /// use the maximum contribution
                                      INPAR::XFEM::MassConservationCombination_sum         /// sum viscous and convective contributions
                                      ),
                               &xfluid_stab);

  setStringToIntegralParameter<int>("MASS_CONSERVATION_SCALING","only_visc","apply additional scaling of penalty term to enforce mass conservation for convection-dominated flow",
                                    tuple<std::string>("full", "only_visc"),
                                    tuple<int>(
                                      INPAR::XFEM::MassConservationScaling_full,           /// apply mass-conserving convective scaling additionally
                                      INPAR::XFEM::MassConservationScaling_only_visc       /// use only the viscous scaling
                                      ),
                               &xfluid_stab);

  BoolParameter("GHOST_PENALTY_STAB","no","switch on/off ghost penalty interface stabilization",&xfluid_stab);

  BoolParameter("GHOST_PENALTY_TRANSIENT_STAB","no","switch on/off ghost penalty transient interface stabilization",&xfluid_stab);

  BoolParameter("GHOST_PENALTY_2nd_STAB","no","switch on/off ghost penalty interface stabilization for 2nd order derivatives",&xfluid_stab);

  DoubleParameter("GHOST_PENALTY_FAC",       0.1, "define stabilization parameter ghost penalty interface stabilization",&xfluid_stab);

  DoubleParameter("GHOST_PENALTY_TRANSIENT_FAC",       0.001, "define stabilization parameter ghost penalty transient interface stabilization",&xfluid_stab);

  BoolParameter("XFF_EOS_PRES_EMB_LAYER","no","switch on/off edge-based pressure stabilization on interface-contributing elements of the embedded fluid",&xfluid_stab);

  BoolParameter("IS_PSEUDO_2D","no","modify viscous interface stabilization due to the vanishing polynomial in third dimension when using strong Dirichlet conditions to block polynomials in one spatial dimension",&xfluid_stab);
}



void INPAR::XFEM::SetValidConditions(const std::vector<Teuchos::RCP<DRT::INPUT::ConditionComponent> > &dirichletbundcomponents,
                                     const std::vector<Teuchos::RCP<DRT::INPUT::ConditionComponent> > &neumanncomponents,
                                     std::vector<Teuchos::RCP<DRT::INPUT::ConditionDefinition> >& condlist)
{
  using namespace DRT::INPUT;

  std::vector<Teuchos::RCP<ConditionComponent> > xfemcomponents;

  xfemcomponents.push_back(Teuchos::rcp(new IntConditionComponent("label")));

  Teuchos::RCP<ConditionDefinition> movingfluid =
        Teuchos::rcp(new ConditionDefinition("DESIGN FLUID MESH VOL CONDITIONS",
                                             "FluidMesh",
                                             "Fluid Mesh",
                                             DRT::Condition::FluidMesh,
                                             true,
                                             DRT::Condition::Volume));
  Teuchos::RCP<ConditionDefinition> fluidfluidcoupling =
      Teuchos::rcp(new ConditionDefinition("DESIGN FLUID FLUID COUPLING SURF CONDITIONS",
                                           "FluidFluidCoupling",
                                           "FLUID FLUID Coupling",
                                           DRT::Condition::FluidFluidCoupling,
                                           true,
                                           DRT::Condition::Surface));
  Teuchos::RCP<ConditionDefinition> ALEfluidcoupling =
      Teuchos::rcp(new ConditionDefinition("DESIGN ALE FLUID COUPLING SURF CONDITIONS",
                                           "ALEFluidCoupling",
                                           "ALE FLUID Coupling",
                                           DRT::Condition::ALEFluidCoupling,
                                           true,
                                           DRT::Condition::Surface));

  for (unsigned i=0; i<xfemcomponents.size(); ++i)
  {
    movingfluid->AddComponent(xfemcomponents[i]);
    fluidfluidcoupling->AddComponent(xfemcomponents[i]);
    ALEfluidcoupling->AddComponent(xfemcomponents[i]);
  }

  condlist.push_back(fluidfluidcoupling);
  condlist.push_back(movingfluid);
  condlist.push_back(ALEfluidcoupling);


  /*--------------------------------------------------------------------*/
  // XFEM coupling conditions

  //*----------------*/
  // Displacement surface condition for XFEM WDBC and Neumann boundary conditions

  Teuchos::RCP<ConditionDefinition> xfem_surf_displacement =
      Teuchos::rcp(new ConditionDefinition(
          "DESIGN XFEM DISPLACEMENT SURF CONDITIONS",
          "XFEMSurfDisplacement",
          "XFEM Surf Displacement",
          DRT::Condition::XFEM_Surf_Displacement,
          true,
          DRT::Condition::Surface));

  xfem_surf_displacement->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("COUPLINGID")));
  xfem_surf_displacement->AddComponent(Teuchos::rcp(new IntVectorConditionComponent("couplingID", 1, false, false, false)));
  xfem_surf_displacement->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("EVALTYPE")));
  xfem_surf_displacement->AddComponent(
    Teuchos::rcp(
      new StringConditionComponent(
        "evaltype","funct",
        Teuchos::tuple<std::string>("zero","funct","implementation"),
        Teuchos::tuple<std::string>("zero","funct","implementation"),
        true)));


  for (unsigned i=0; i<dirichletbundcomponents.size(); ++i)
  {
    xfem_surf_displacement->AddComponent(dirichletbundcomponents[i]);
  }

  condlist.push_back(xfem_surf_displacement);



  //*----------------*/
  // Levelset field condition components

  std::vector<Teuchos::RCP<ConditionComponent> > levelsetfield_components;

  levelsetfield_components.push_back(Teuchos::rcp(new SeparatorConditionComponent("LEVELSETFIELDNO")));
  levelsetfield_components.push_back(Teuchos::rcp(new IntVectorConditionComponent("levelsetfieldno", 1, false, false, false)));

  levelsetfield_components.push_back(Teuchos::rcp(new SeparatorConditionComponent("LEVELSETCURVE")));
  levelsetfield_components.push_back(Teuchos::rcp(new IntVectorConditionComponent("levelsetcurve", 1, true, true, false)));


  //*----------------*/
  // Levelset based Weak Dirichlet conditions

  Teuchos::RCP<ConditionDefinition> xfem_levelset_wdbc =
      Teuchos::rcp(new ConditionDefinition(
          "DESIGN XFEM LEVELSET WEAK DIRICHLET VOL CONDITIONS",
          "XFEMLevelsetWeakDirichlet",
          "XFEM Levelset Weak Dirichlet",
          DRT::Condition::XFEM_Levelset_Weak_Dirichlet,
          true,
          DRT::Condition::Volume));

  for (unsigned i=0; i<levelsetfield_components.size(); ++i)
  {
    xfem_levelset_wdbc->AddComponent(levelsetfield_components[i]);
  }
  for (unsigned i=0; i<dirichletbundcomponents.size(); ++i)
  {
    xfem_levelset_wdbc->AddComponent(dirichletbundcomponents[i]);
  }

  condlist.push_back(xfem_levelset_wdbc);

  //*----------------*/
  // Levelset based Neumann conditions

  Teuchos::RCP<ConditionDefinition> xfem_levelset_neumann =
      Teuchos::rcp(new ConditionDefinition(
          "DESIGN XFEM LEVELSET NEUMANN VOL CONDITIONS",
          "XFEMLevelsetNeumann",
          "XFEM Levelset Neumann",
          DRT::Condition::XFEM_Levelset_Neumann,
          true,
          DRT::Condition::Volume));

  for (unsigned i=0; i<levelsetfield_components.size(); ++i)
  {
    xfem_levelset_neumann->AddComponent(levelsetfield_components[i]);
  }
  for (unsigned i=0; i<neumanncomponents.size(); ++i)
  {
    xfem_levelset_neumann->AddComponent(neumanncomponents[i]);
  }

  condlist.push_back(xfem_levelset_neumann);

  //*----------------*/
  // Levelset based Twophase conditions

  Teuchos::RCP<ConditionDefinition> xfem_levelset_twophase =
      Teuchos::rcp(new ConditionDefinition(
          "DESIGN XFEM LEVELSET TWOPHASE VOL CONDITIONS",
          "XFEMLevelsetTwophase",
          "XFEM Levelset Twophase",
          DRT::Condition::XFEM_Levelset_Twophase,
          true,
          DRT::Condition::Volume));

  for (unsigned i=0; i<levelsetfield_components.size(); ++i)
  {
    xfem_levelset_twophase->AddComponent(levelsetfield_components[i]);
  }

  condlist.push_back(xfem_levelset_twophase);

  //*----------------*/
  // Levelset based Combustion conditions

  Teuchos::RCP<ConditionDefinition> xfem_levelset_combustion =
      Teuchos::rcp(new ConditionDefinition(
          "DESIGN XFEM LEVELSET COMBUSTION VOL CONDITIONS",
          "XFEMLevelsetCombustion",
          "XFEM Levelset Combustion",
          DRT::Condition::XFEM_Levelset_Combustion,
          true,
          DRT::Condition::Volume));

  for (unsigned i=0; i<levelsetfield_components.size(); ++i)
  {
    xfem_levelset_combustion->AddComponent(levelsetfield_components[i]);
  }

  condlist.push_back(xfem_levelset_combustion);

  //*----------------*/
  // Surface Fluid-Fluid coupling conditions


  std::vector<Teuchos::RCP<ConditionComponent> > xfluidfluidsurfcomponents;

  xfluidfluidsurfcomponents.push_back(Teuchos::rcp(new IntConditionComponent("label")));
  xfluidfluidsurfcomponents.push_back(
      Teuchos::rcp(
           new StringConditionComponent(
             "COUPSTRATEGY",
             "xfluid",
             Teuchos::tuple<std::string>("xfluid","embedded","mean"),
             Teuchos::tuple<int>(
                 INPAR::XFEM::Xfluid_Sided,
                 INPAR::XFEM::Embedded_Sided,
                 INPAR::XFEM::Mean)
             )));


  Teuchos::RCP<ConditionDefinition> xfem_surf_fluidfluid =
      Teuchos::rcp(new ConditionDefinition(
          "DESIGN XFEM FLUIDFLUID SURF CONDITIONS",
          "XFEMSurfFluidFluid",
          "XFEM Surf FluidFluid",
          DRT::Condition::XFEM_Surf_FluidFluid,
          true,
          DRT::Condition::Surface));

  for (unsigned i=0; i<xfluidfluidsurfcomponents.size(); ++i)
    xfem_surf_fluidfluid->AddComponent(xfluidfluidsurfcomponents[i]);

  condlist.push_back(xfem_surf_fluidfluid);

  //*----------------*/
  // Surface partitioned XFSI boundary conditions

  Teuchos::RCP<ConditionDefinition> xfem_surf_fsi_part =
      Teuchos::rcp(new ConditionDefinition(
          "DESIGN XFEM FSI PARTITIONED SURF CONDITIONS",
          "XFEMSurfFSIPart",
          "XFEM Surf FSI Part",
          DRT::Condition::XFEM_Surf_FSIPart,
          true,
          DRT::Condition::Surface));

  condlist.push_back(xfem_surf_fsi_part);

  //*----------------*/
  // Surface monolithic XFSI coupling conditions

  Teuchos::RCP<ConditionDefinition> xfem_surf_fsi_mono =
      Teuchos::rcp(new ConditionDefinition(
          "DESIGN XFEM FSI MONOLITHIC SURF CONDITIONS",
          "XFEMSurfFSIMono",
          "XFEM Surf FSI Mono",
          DRT::Condition::XFEM_Surf_FSIMono,
          true,
          DRT::Condition::Surface));

  condlist.push_back(xfem_surf_fsi_mono);


  //*----------------*/
  // Surface partitioned CRACK XFSI boundary conditions

  Teuchos::RCP<ConditionDefinition> xfem_surf_crfsi_part =
      Teuchos::rcp(new ConditionDefinition(
          "DESIGN XFEM CRACK FSI PARTITIONED SURF CONDITIONS",
          "XFEMSurfCrackFSIPart",
          "XFEM Surf Crack FSI Part",
          DRT::Condition::XFEM_Surf_CrackFSIPart,
          true,
          DRT::Condition::Surface));

  condlist.push_back(xfem_surf_crfsi_part);


  //*----------------*/
  // Surface Weak Dirichlet conditions

  Teuchos::RCP<ConditionDefinition> xfem_surf_wdbc =
      Teuchos::rcp(new ConditionDefinition(
          "DESIGN XFEM WEAK DIRICHLET SURF CONDITIONS",
          "XFEMSurfWeakDirichlet",
          "XFEM Surf Weak Dirichlet",
          DRT::Condition::XFEM_Surf_Weak_Dirichlet,
          true,
          DRT::Condition::Surface));

  xfem_surf_wdbc->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("COUPLINGID")));
  xfem_surf_wdbc->AddComponent(Teuchos::rcp(new IntVectorConditionComponent("couplingID", 1, false, false, false)));
  xfem_surf_wdbc->AddComponent(Teuchos::rcp(new SeparatorConditionComponent("EVALTYPE")));
  xfem_surf_wdbc->AddComponent(
    Teuchos::rcp(
      new StringConditionComponent(
        "evaltype","funct_interpolated",
        Teuchos::tuple<std::string>(
            "zero",
            "funct_interpolated",
            "funct_gausspoint",
            "displacement_1storder_wo_initfunct",
            "displacement_2ndorder_wo_initfunct",
            "displacement_1storder_with_initfunct",
            "displacement_2ndorder_with_initfunct"),
        Teuchos::tuple<std::string>(
            "zero",
            "funct_interpolated",
            "funct_gausspoint",
            "displacement_1storder_wo_initfunct",
            "displacement_2ndorder_wo_initfunct",
            "displacement_1storder_with_initfunct",
            "displacement_2ndorder_with_initfunct"),
        true)));

  for (unsigned i=0; i<dirichletbundcomponents.size(); ++i)
  {
    xfem_surf_wdbc->AddComponent(dirichletbundcomponents[i]);
  }

  condlist.push_back(xfem_surf_wdbc);


  //*----------------*/
  // Surface Neumann conditions

  Teuchos::RCP<ConditionDefinition> xfem_surf_neumann =
      Teuchos::rcp(new ConditionDefinition(
          "DESIGN XFEM NEUMANN SURF CONDITIONS",
          "XFEMSurfNeumann",
          "XFEM Surf Neumann",
          DRT::Condition::XFEM_Surf_Neumann,
          true,
          DRT::Condition::Surface));

  for (unsigned i=0; i<neumanncomponents.size(); ++i)
  {
    xfem_surf_neumann->AddComponent(neumanncomponents[i]);
  }

  condlist.push_back(xfem_surf_neumann);


}
