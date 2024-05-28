/*----------------------------------------------------------------------*/
/*! \file

\brief Declaration

\level 1

*/
/*----------------------------------------------------------------------*/


#include "4C_linear_solver_amgnxn_hierarchies.hpp"

#include "4C_linalg_multiply.hpp"
#include "4C_linear_solver_amgnxn_vcycle.hpp"
#include "4C_utils_exceptions.hpp"

#include <EpetraExt_RowMatrixOut.h>
#include <MueLu_EpetraOperator.hpp>
#include <MueLu_MLParameterListInterpreter_decl.hpp>
#include <MueLu_ParameterListInterpreter.hpp>
#include <Teuchos_PtrDecl.hpp>
#include <Teuchos_Time.hpp>
#include <Teuchos_XMLParameterListHelpers.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include <iostream>

FOUR_C_NAMESPACE_OPEN


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::Hierarchies(Teuchos::RCP<AMGNXN::BlockedMatrix> A,
    std::vector<Teuchos::ParameterList> muelu_params, std::vector<int> num_pdes,
    std::vector<int> null_spaces_dim,
    std::vector<Teuchos::RCP<std::vector<double>>> null_spaces_data, int NumLevelAMG,
    std::string verbosity)
    : a_(A),
      muelu_params_(muelu_params),
      num_pdes_(num_pdes),
      null_spaces_dim_(null_spaces_dim),
      null_spaces_data_(null_spaces_data),
      num_blocks_(A->GetNumRows()),
      num_level_amg_(NumLevelAMG),
      verbosity_(verbosity)
{
  // Plausibility checks
  if (a_->GetNumRows() != a_->GetNumCols() or (int)(muelu_params_.size()) != num_blocks_ or
      (int)(num_pdes_.size()) != num_blocks_ or (int)(null_spaces_dim_.size()) != num_blocks_ or
      (int)(null_spaces_data_.size()) != num_blocks_)
    FOUR_C_THROW("Something wrong");

  // Setput
  setup();
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::GetNumLevels(int block)
{
  if (h_block_[block] == Teuchos::null)
    return num_level_max_;
  else
    return h_block_[block]->GetNumLevels();
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>>
CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::GetA(int block)
{
  return a_block_level_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>>
CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::GetP(int block)
{
  return p_block_level_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>>
CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::GetR(int block)
{
  return r_block_level_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::vector<Teuchos::RCP<CORE::LINEAR_SOLVER::AMGNXN::MueluSmootherWrapper>>
CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::GetSPre(int block)
{
  return s_pre_block_level_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

std::vector<Teuchos::RCP<CORE::LINEAR_SOLVER::AMGNXN::MueluSmootherWrapper>>
CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::GetSPos(int block)
{
  return s_pos_block_level_[block];
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::GetNumPDEs(int block) { return num_pdes_[block]; }

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::GetNullSpaceDim(int block)
{
  return null_spaces_dim_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<std::vector<double>> CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::GetNullSpaceData(
    int block)
{
  return null_spaces_data_[block];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::setup()
{
  TEUCHOS_FUNC_TIME_MONITOR("CORE::LINALG::SOLVER::AMGNXN::Hierarchies::Setup");

  using MueLuUtils = MueLu::Utilities<double, int, int, Node>;

  // ===========================================================
  // Build up MueLu Hierarchies of each one of the blocks
  // ===========================================================

  h_block_.assign(num_blocks_, Teuchos::null);
  std::vector<int> offsets(num_level_amg_ - 1, 0);
  for (int block = 0; block < num_blocks_; block++)
  {
    int offsetFineLevel = a_->GetMatrix(block, block)->RowMap().MinAllGID();
    Teuchos::RCP<Epetra_Operator> A_eop = a_->GetMatrix(block, block)->EpetraOperator();
    h_block_[block] =
        build_mue_lu_hierarchy(muelu_params_[block], num_pdes_[block], null_spaces_dim_[block],
            null_spaces_data_[block], A_eop, block, num_blocks_, offsets, offsetFineLevel);
  }

  // ===========================================================
  // Determine number of levels
  // ===========================================================

  num_level_max_ = -10000000;
  num_level_min_ = 10000000;
  for (int block = 0; block < num_blocks_; block++)
  {
    if (h_block_[block] == Teuchos::null) continue;
    int NumLevel_this_block = h_block_[block]->GetNumLevels();
    if (NumLevel_this_block > num_level_max_) num_level_max_ = NumLevel_this_block;
    if (NumLevel_this_block < num_level_min_) num_level_min_ = NumLevel_this_block;
  }


  // ===========================================================
  // Extract matrices, transfer operators and smoothers from the hierarchies
  // ===========================================================

  for (int block = 0; block < num_blocks_; block++)
  {
    // create a dummy hierarchy by repeating the same matrix
    if (h_block_[block] == Teuchos::null)
    {
      std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> A_level(num_level_max_, Teuchos::null);
      std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> P_level(
          num_level_max_ - 1, Teuchos::null);
      std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> R_level(
          num_level_max_ - 1, Teuchos::null);
      std::vector<Teuchos::RCP<AMGNXN::MueluSmootherWrapper>> SPre_level(
          num_level_max_, Teuchos::null);
      std::vector<Teuchos::RCP<AMGNXN::MueluSmootherWrapper>> SPos_level(
          num_level_max_ - 1, Teuchos::null);

      Teuchos::RCP<CORE::LINALG::SparseMatrix> Abb = a_->GetMatrix(block, block);
      Teuchos::RCP<CORE::LINALG::SparseMatrix> Peye = CORE::LINALG::Eye(Abb->DomainMap());
      Teuchos::RCP<CORE::LINALG::SparseMatrix> Reye = CORE::LINALG::Eye(Abb->RangeMap());

      for (int level = 0; level < num_level_max_; level++) A_level[level] = Abb;

      for (int level = 0; level < num_level_max_ - 1; level++)
      {
        P_level[level] = Peye;
        R_level[level] = Reye;
      }

      a_block_level_.push_back(A_level);
      p_block_level_.push_back(P_level);
      r_block_level_.push_back(R_level);
      s_pre_block_level_.push_back(SPre_level);
      s_pos_block_level_.push_back(SPos_level);
    }
    else  // Recover objects created by muelu
    {
      int NumLevel_this_block = h_block_[block]->GetNumLevels();
      std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> A_level(
          NumLevel_this_block, Teuchos::null);
      std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> P_level(
          NumLevel_this_block - 1, Teuchos::null);
      std::vector<Teuchos::RCP<CORE::LINALG::SparseMatrix>> R_level(
          NumLevel_this_block - 1, Teuchos::null);
      std::vector<Teuchos::RCP<AMGNXN::MueluSmootherWrapper>> SPre_level(
          NumLevel_this_block, Teuchos::null);
      std::vector<Teuchos::RCP<AMGNXN::MueluSmootherWrapper>> SPos_level(
          NumLevel_this_block, Teuchos::null);

      // some local typedefs
      typedef Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> Matrix;
      typedef MueLu::SmootherBase<Scalar, LocalOrdinal, GlobalOrdinal, Node> SmootherBase;

      Teuchos::RCP<Matrix> myA = Teuchos::null;
      Teuchos::RCP<Epetra_CrsMatrix> myAcrs = Teuchos::null;
      Teuchos::RCP<CORE::LINALG::SparseMatrix> myAspa = Teuchos::null;
      Teuchos::RCP<SmootherBase> myS = Teuchos::null;
      Teuchos::RCP<CORE::LINEAR_SOLVER::AMGNXN::MueluSmootherWrapper> mySWrap = Teuchos::null;

      bool explicitdirichlet = a_->GetMatrix(0, 0)->ExplicitDirichlet();
      bool savegraph = a_->GetMatrix(0, 0)->SaveGraph();

      for (int level = 0; level < NumLevel_this_block; level++)
      {
        Teuchos::RCP<MueLu::Level> this_level = h_block_[block]->GetLevel(level);
        if (this_level->IsAvailable("A"))
        {
          myA = this_level->Get<Teuchos::RCP<Matrix>>("A");
          myAcrs = MueLuUtils::Op2NonConstEpetraCrs(myA);
          myAspa = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
              myAcrs, CORE::LINALG::Copy, explicitdirichlet, savegraph));
          A_level[level] = myAspa;
        }
        else
          FOUR_C_THROW("Error in extracting A");

        if (this_level->IsAvailable("PreSmoother"))
        {
          myS = this_level->Get<Teuchos::RCP<SmootherBase>>("PreSmoother");
          mySWrap = Teuchos::rcp(new CORE::LINEAR_SOLVER::AMGNXN::MueluSmootherWrapper(myS));
          SPre_level[level] = mySWrap;
        }
        else
          FOUR_C_THROW("Error in extracting PreSmoother");

        if (level < NumLevel_this_block - 1)
        {
          if (this_level->IsAvailable("PostSmoother"))
          {
            myS = this_level->Get<Teuchos::RCP<SmootherBase>>("PostSmoother");
            mySWrap = Teuchos::rcp(new CORE::LINEAR_SOLVER::AMGNXN::MueluSmootherWrapper(myS));
            SPos_level[level] = mySWrap;
          }
          else
            FOUR_C_THROW("Error in extracting PostSmoother");
        }

        if (level != 0)
        {
          if (this_level->IsAvailable("P"))
          {
            myA = this_level->Get<Teuchos::RCP<Matrix>>("P");
            myAcrs = MueLuUtils::Op2NonConstEpetraCrs(myA);
            myAspa = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
                myAcrs, CORE::LINALG::Copy, explicitdirichlet, savegraph));
            P_level[level - 1] = myAspa;
          }
          else
            FOUR_C_THROW("Error in extracting P");

          if (this_level->IsAvailable("R"))
          {
            myA = this_level->Get<Teuchos::RCP<Matrix>>("R");
            myAcrs = MueLuUtils::Op2NonConstEpetraCrs(myA);
            myAspa = Teuchos::rcp(new CORE::LINALG::SparseMatrix(
                myAcrs, CORE::LINALG::Copy, explicitdirichlet, savegraph));
            R_level[level - 1] = myAspa;
          }
          else
            FOUR_C_THROW("Error in extracting R");
        }

      }  // loop in levels

      a_block_level_.push_back(A_level);
      p_block_level_.push_back(P_level);
      r_block_level_.push_back(R_level);
      s_pre_block_level_.push_back(SPre_level);
      s_pos_block_level_.push_back(SPos_level);

    }  // else


  }  // loop in blocks

  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::build_mue_lu_hierarchy(
    Teuchos::ParameterList paramListFromXml, int numdf, int dimns,
    Teuchos::RCP<std::vector<double>> nsdata, Teuchos::RCP<Epetra_Operator> A_eop, int block,
    int NumBlocks, std::vector<int>& offsets, int offsetFineLevel)
{
  TEUCHOS_FUNC_TIME_MONITOR("CORE::LINALG::SOLVER::AMGNXN::Hierarchies::build_mue_lu_hierarchy");

  using MueLuUtils = MueLu::Utilities<double, int, int, Node>;

  Teuchos::Time timer("", true);
  timer.reset();

  Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> H = Teuchos::null;
  bool create_uncoarsened_hierarchy =
      paramListFromXml.get<bool>("create un-coarsened hierarchy", false);
  bool fix_coarse_maps = paramListFromXml.get<bool>(
      "fix coarse maps", false);  // this is required in all the fields if we want to merge and
                                  // solve a coarse level block matrix

  if (not create_uncoarsened_hierarchy)
  {
    // Some cheks
    if (numdf < 1 or dimns < 1) FOUR_C_THROW("Error: PDE equations or null space dimension wrong.");
    if (nsdata == Teuchos::null) FOUR_C_THROW("Error: null space data is empty");

    // Hack making TSI work with the trilinos Q1_2015. The Q1_2014 version worked without this
    if (numdf == 1) offsetFineLevel = 0;

    // Prepare operator for MueLu
    Teuchos::RCP<Epetra_CrsMatrix> A_crs = Teuchos::rcp_dynamic_cast<Epetra_CrsMatrix>(A_eop);
    if (A_crs == Teuchos::null)
      FOUR_C_THROW("Make sure that the input matrix is a Epetra_CrsMatrix (or derived)");
    Teuchos::RCP<Xpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mueluA =
        Teuchos::rcp(new Xpetra::EpetraCrsMatrixT<int, Xpetra::EpetraNode>(A_crs));

    Teuchos::RCP<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mueluA_wrap =
        Teuchos::rcp(new Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>(mueluA));
    Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> mueluOp =
        Teuchos::rcp_dynamic_cast<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(
            mueluA_wrap);
    mueluOp->SetFixedBlockSize(numdf, offsetFineLevel);

    // Prepare null space vector for MueLu
    Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> rowMap = mueluA->getRowMap();
    Teuchos::RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> nspVector =
        Xpetra::MultiVectorFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(
            rowMap, dimns, true);
    for (size_t i = 0; i < Teuchos::as<size_t>(dimns); i++)
    {
      Teuchos::ArrayRCP<Scalar> nspVectori = nspVector->getDataNonConst(i);
      const size_t myLength = nspVector->getLocalLength();
      for (size_t j = 0; j < myLength; j++)
      {
        nspVectori[j] = (*nsdata)[i * myLength + j];
      }
    }



    if (fix_coarse_maps)
    {
      // Add information about maps offsets
      std::string offsets_str("{");
      for (int i = 0; i < (int)offsets.size(); i++)
      {
        offsets_str = offsets_str + convert_int(offsets[i]);
        if (i < (int)offsets.size() - 1) offsets_str = offsets_str + ", ";
      }
      offsets_str = offsets_str + "}";
      if (paramListFromXml.sublist("Factories", true).isSublist("myCoarseMapFactory123"))
        FOUR_C_THROW(
            "We are going to overwrite the factory 'myCoarseMapFactory123'. Please use an other "
            "name");
      Teuchos::ParameterList& myCoarseMapFactoryList =
          paramListFromXml.sublist("Factories", true).sublist("myCoarseMapFactory123");
      myCoarseMapFactoryList.set("factory", "CoarseMapFactory");
      myCoarseMapFactoryList.set("Domain GID offsets", offsets_str);
      if (paramListFromXml.sublist("Hierarchy", true).sublist("All").isParameter("CoarseMap"))
        FOUR_C_THROW("We are going to overwrite 'CoarseMap'. Don't use 'CoarseMap' here.");
      Teuchos::ParameterList& AllList = paramListFromXml.sublist("Hierarchy").sublist("All");
      AllList.set("CoarseMap", "myCoarseMapFactory123");

      if (A_eop->Comm().MyPID() == 0) std::cout << "offsets_str " << offsets_str << std::endl;
    }

    // Add offset for the finest level
    Teuchos::ParameterList& MatrixList = paramListFromXml.sublist("Matrix");
    MatrixList.set<int>("DOF offset", offsetFineLevel);
    MatrixList.set<int>("number of equations", numdf);

    if (verbosity_ == "on" and A_eop->Comm().MyPID() == 0)
    {
      std::cout << "offsetFineLevel " << offsetFineLevel << std::endl;
    }

    // Build up hierarchy
    MueLu::ParameterListInterpreter<Scalar, LocalOrdinal, GlobalOrdinal, Node> mueLuFactory(
        paramListFromXml);
    H = mueLuFactory.CreateHierarchy();
    H->SetDefaultVerbLevel(MueLu::Extreme);  // TODO sure?
    H->GetLevel(0)->Set("A", mueluOp);
    H->GetLevel(0)->Set("Nullspace", nspVector);
    H->GetLevel(0)->setlib(Xpetra::UseEpetra);
    H->setlib(Xpetra::UseEpetra);
    mueLuFactory.SetupHierarchy(*H);

    // Recover information about the maps
    if (fix_coarse_maps)
    {
      int NumLevel_block = H->GetNumLevels();
      Teuchos::RCP<MueLu::Level> this_level = Teuchos::null;
      Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> myA = Teuchos::null;
      Teuchos::RCP<Epetra_CrsMatrix> myAcrs = Teuchos::null;
      for (int level = 1; (level < NumLevel_block) and (level < (int)offsets.size() + 1); level++)
      {
        this_level = H->GetLevel(level);
        if (this_level->IsAvailable("A"))
        {
          myA = this_level
                    ->Get<Teuchos::RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>>>(
                        "A");  // Matrix
          myAcrs = MueLuUtils::Op2NonConstEpetraCrs(myA);
        }
        else
          FOUR_C_THROW("Error in extracting A");

        offsets[level - 1] = offsets[level - 1] + myAcrs->RangeMap().MaxAllGID() +
                             1;  // TODO I think we don't have to overwrite previous result
      }
    }
  }
  else  // when we use a dummy hierarchy we sill have to compute the offsets
  {
    if (fix_coarse_maps)
    {
      for (int level = 0; level < (int)offsets.size(); level++)
        offsets[level] = offsets[level] + A_eop->OperatorRangeMap().MaxAllGID() +
                         1;  // TODO I think we don't have to overwrite previous result
    }
  }

  double elaptime = timer.totalElapsedTime(true);
  if (verbosity_ == "on" and A_eop->Comm().MyPID() == 0)
    std::cout
        << "       Calling CORE::LINALG::SOLVER::AMGNXN::Hierarchies::build_mue_lu_hierarchy takes "
        << std::setw(16) << std::setprecision(6) << elaptime << " s" << std::endl;
  return H;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::GetNumLevelMin() { return num_level_min_; }

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::GetNumBlocks() { return num_blocks_; }

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<CORE::LINEAR_SOLVER::AMGNXN::BlockedMatrix>
CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::GetBlockMatrix()
{
  if (a_ == Teuchos::null) FOUR_C_THROW("No data");
  return a_;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>>
CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::GetH(int block)
{
  Teuchos::RCP<MueLu::Hierarchy<Scalar, LocalOrdinal, GlobalOrdinal, Node>> H = h_block_[block];
  if (H == Teuchos::null) FOUR_C_THROW("No data");
  return H;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<CORE::LINALG::SparseMatrix> CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::GetA(
    int block, int level)
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> A = a_block_level_[block][level];
  if (A == Teuchos::null) FOUR_C_THROW("No data");
  return A;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<CORE::LINALG::SparseMatrix> CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::GetP(
    int block, int level)
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> P = p_block_level_[block][level];
  if (P == Teuchos::null) FOUR_C_THROW("No data");
  return P;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<CORE::LINALG::SparseMatrix> CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::GetR(
    int block, int level)
{
  Teuchos::RCP<CORE::LINALG::SparseMatrix> R = r_block_level_[block][level];
  if (R == Teuchos::null) FOUR_C_THROW("No data");
  return R;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<CORE::LINEAR_SOLVER::AMGNXN::MueluSmootherWrapper>
CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::GetSPre(int block, int level)
{
  Teuchos::RCP<AMGNXN::MueluSmootherWrapper> SPre = s_pre_block_level_[block][level];
  if (SPre == Teuchos::null) FOUR_C_THROW("No data");
  return SPre;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<CORE::LINEAR_SOLVER::AMGNXN::MueluSmootherWrapper>
CORE::LINEAR_SOLVER::AMGNXN::Hierarchies::GetSPos(int block, int level)
{
  Teuchos::RCP<AMGNXN::MueluSmootherWrapper> SPos = s_pos_block_level_[block][level];
  if (SPos == Teuchos::null) FOUR_C_THROW("No data");
  return SPos;
}



/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

CORE::LINEAR_SOLVER::AMGNXN::MonolithicHierarchy::MonolithicHierarchy(
    Teuchos::RCP<AMGNXN::Hierarchies> H, const Teuchos::ParameterList& params,
    const Teuchos::ParameterList& params_smoothers)
    : h_(H), params_(params), params_smoothers_(params_smoothers)
{
  // Expected parameters in params (example)
  //<ParameterList name="params">
  //
  //  <Parameter name="number of levels"                 type="int"  value="..."/>
  //
  //  <Parameter name="smoother: all but coarsest level" type="string"  value="myFinestSmoother"/>
  //
  //  <Parameter name="smoother: coarsest level"         type="string"  value="myCoarsestSmoother"/>
  //
  //  <Parameter name="verbosity"                        type="string"  value="on"/>
  //
  //</ParameterList>

  // Expected parameters in params_smoothers (example)
  //<ParameterList name="params_smoothers">
  //
  //  <ParameterList name="myFinestSmoother">
  //
  //   ...    ...    ...    ...    ...
  //
  //  </ParameterList>
  //
  //  <ParameterList name="myCoarsestSmoother">
  //
  //   ...    ...    ...    ...    ...
  //
  //  </ParameterList>
  //
  //</ParameterList>

  setup();
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

int CORE::LINEAR_SOLVER::AMGNXN::MonolithicHierarchy::GetNumLevels() { return num_levels_; }


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

void CORE::LINEAR_SOLVER::AMGNXN::MonolithicHierarchy::setup()
{
  TEUCHOS_FUNC_TIME_MONITOR("CORE::LINALG::SOLVER::AMGNXN::MonolithicHierarchy::Setup()");

  // ====================================================
  // Create block transfer operators
  // ====================================================

  num_levels_ = params_.get<int>("number of levels", -1);
  if (num_levels_ == -1) FOUR_C_THROW("Missing \"number of levels\" in your xml file");
  num_levels_ = std::min(num_levels_, h_->GetNumLevelMin());
  num_blocks_ = h_->GetNumBlocks();


  // ====================================================
  // Some output
  // ====================================================
  std::string verbosity = params_.get<std::string>("verbosity", "off");

  if (h_->GetBlockMatrix()->GetMatrix(0, 0)->Comm().MyPID() != 0) verbosity = "off";

  if (verbosity == "on")
  {
    // std::cout << "===============================================" << std::endl;
    // std::cout << "AMGnxn preconditioner: debug info  (begin)" << std::endl;
    // std::cout << std::endl;
    // std::cout << "===============================================" << std::endl;
    std::cout << "number of blocks = " << num_blocks_ << std::endl;
    std::cout << "number of levels = " << num_levels_ << std::endl;
    for (int i = 0; i < num_blocks_; i++)
      std::cout << "block " << i << ": number of levels = " << GetHierarchies()->GetNumLevels(i)
                << std::endl;
  }



  // ====================================================
  // Create block transfer operators
  // ====================================================

  p_.assign(num_levels_ - 1, Teuchos::null);
  r_.assign(num_levels_ - 1, Teuchos::null);
  for (int level = 0; level < num_levels_ - 1; level++)
  {
    p_[level] = Teuchos::rcp(new AMGNXN::DiagonalBlockedMatrix(num_blocks_));
    r_[level] = Teuchos::rcp(new AMGNXN::DiagonalBlockedMatrix(num_blocks_));
    for (int block = 0; block < num_blocks_; block++)
    {
      p_[level]->SetMatrix(h_->GetP(block, level), block, block);
      r_[level]->SetMatrix(h_->GetR(block, level), block, block);
    }
  }

  // ====================================================
  // Create coarser matrices
  // ====================================================
  //

  a_.assign(num_levels_, Teuchos::null);
  for (int level = 0; level < num_levels_; level++)
  {
    if (level == 0)
      a_[level] = h_->GetBlockMatrix();
    else
    {
      a_[level] = Teuchos::rcp(new AMGNXN::BlockedMatrix(num_blocks_, num_blocks_));

      for (int block = 0; block < num_blocks_; block++)
        a_[level]->SetMatrix(h_->GetA(block, level), block, block);

      for (int row = 0; row < num_blocks_; row++)
        for (int col = 0; col < num_blocks_; col++)
        {
          if (row != col)  // The diagonal blocks have been already assigned
          {
            Teuchos::RCP<CORE::LINALG::SparseMatrix> A_spa = a_[level - 1]->GetMatrix(row, col);
            Teuchos::RCP<CORE::LINALG::SparseMatrix> P_spa = p_[level - 1]->GetMatrix(col, col);
            Teuchos::RCP<CORE::LINALG::SparseMatrix> R_spa = r_[level - 1]->GetMatrix(row, row);

            Teuchos::RCP<CORE::LINALG::SparseMatrix> AP_spa = Teuchos::null;
            AP_spa = CORE::LINALG::MLMultiply(*A_spa, *P_spa, true);
            if (AP_spa == Teuchos::null) FOUR_C_THROW("Error in AP");

            Teuchos::RCP<CORE::LINALG::SparseMatrix> RAP_spa = Teuchos::null;
            RAP_spa = CORE::LINALG::MLMultiply(*R_spa, *AP_spa, true);
            if (RAP_spa == Teuchos::null) FOUR_C_THROW("Error in RAP");

            a_[level]->SetMatrix(RAP_spa, row, col);
          }
        }
    }
  }

  // ====================================================
  // Create smoothers
  // ====================================================

  spre_.assign(num_levels_, Teuchos::null);
  spos_.assign(num_levels_ - 1, Teuchos::null);
  for (int level = 0; level < num_levels_; level++)
  {
    if (level < num_levels_ - 1)
    {
      spre_[level] = build_smoother(level);
      spos_[level] = spre_[level];  // build_smoother(level);
    }
    else
    {
      spre_[level] = build_smoother(level);
    }
  }

  // if (verbosity=="on")
  //{
  //  //std::cout << "===============================================" << std::endl;
  //  std::cout << std::endl;
  //  std::cout << "AMGnxn preconditioner: debug info  (end)" << std::endl;
  //  std::cout << "===============================================" << std::endl;
  //}

  return;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<CORE::LINEAR_SOLVER::AMGNXN::BlockedMatrix>
CORE::LINEAR_SOLVER::AMGNXN::MonolithicHierarchy::GetA(int level)
{
  return a_[level];
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<CORE::LINEAR_SOLVER::AMGNXN::Hierarchies>
CORE::LINEAR_SOLVER::AMGNXN::MonolithicHierarchy::GetHierarchies()
{
  return h_;
}

/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<CORE::LINEAR_SOLVER::AMGNXN::GenericSmoother>
CORE::LINEAR_SOLVER::AMGNXN::MonolithicHierarchy::build_smoother(int level)
{
  TEUCHOS_FUNC_TIME_MONITOR("CORE::LINALG::SOLVER::AMGNXN::MonolithicHierarchy::build_smoother");

  std::string smother_name;
  if (level < num_levels_ - 1)
  {
    smother_name = params_.get<std::string>("smoother: all but coarsest level", "none");
    if (smother_name == "none")
      FOUR_C_THROW("You have to set the fine level smoother. Fix your xml file.");
  }
  else
  {
    smother_name = params_.get<std::string>("smoother: coarsest level", "none");
    if (smother_name == "none")
      FOUR_C_THROW("You have to set the coarse level smoother. Fix your xml file.");
  }

  std::string verbosity = params_.get<std::string>("verbosity", "off");

  if (GetA(0)->GetMatrix(0, 0)->Comm().MyPID() != 0) verbosity = "off";

  std::vector<int> blocks(h_->GetNumBlocks(), 0);
  for (int i = 0; i < h_->GetNumBlocks(); i++) blocks[i] = i;

  AMGNXN::SmootherFactory mySmootherCreator;
  mySmootherCreator.SetOperator(GetA(level));
  mySmootherCreator.SetParamsSmoother(params_smoothers_);
  mySmootherCreator.SetHierarchies(GetHierarchies());
  mySmootherCreator.SetLevel(level);
  mySmootherCreator.SetBlocks(blocks);
  mySmootherCreator.SetSmootherName(smother_name);
  mySmootherCreator.SetVerbosity(verbosity);

  // Recover null spaces from the hierarchies and give it to the smoother creator
  std::vector<AMGNXN::NullSpaceInfo> null_space_blocks;
  for (int i = 0; i < num_blocks_; i++)
  {
    AMGNXN::NullSpaceInfo myNS(h_->GetNumPDEs(i), h_->GetNullSpaceDim(i), h_->GetNullSpaceData(i));
    null_space_blocks.push_back(myNS);
  }
  mySmootherCreator.set_null_space_all_blocks(null_space_blocks);

  Teuchos::RCP<AMGNXN::GenericSmoother> Sbase = mySmootherCreator.Create();
  Teuchos::RCP<AMGNXN::BlockedSmoother> Sblock =
      Teuchos::rcp_dynamic_cast<AMGNXN::BlockedSmoother>(Sbase);
  if (Sblock == Teuchos::null)
    FOUR_C_THROW("We expect a block smoother. Fix the xml file defining the smoother");
  return Sbase;
}


/*------------------------------------------------------------------------------*/
/*------------------------------------------------------------------------------*/

Teuchos::RCP<CORE::LINEAR_SOLVER::AMGNXN::Vcycle>
CORE::LINEAR_SOLVER::AMGNXN::MonolithicHierarchy::BuildVCycle()
{
  int NumSweeps = 1;  // Hard coded
  int FirstLevel = 0;
  Teuchos::RCP<Vcycle> V = Teuchos::rcp(new Vcycle(num_levels_, NumSweeps, FirstLevel));

  V->SetOperators(a_);
  V->SetProjectors(p_);
  V->SetRestrictors(r_);
  V->SetPreSmoothers(spre_);
  V->SetPosSmoothers(spos_);

  return V;
}

FOUR_C_NAMESPACE_CLOSE