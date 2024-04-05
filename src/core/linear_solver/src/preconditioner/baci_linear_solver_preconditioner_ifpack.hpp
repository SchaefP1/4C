/*----------------------------------------------------------------------------*/
/*! \file

\brief CORE::LINALG::SOLVER wrapper around Trilinos' IFPACK preconditioner

\level 0

*/
/*----------------------------------------------------------------------------*/
#ifndef BACI_LINEAR_SOLVER_PRECONDITIONER_IFPACK_HPP
#define BACI_LINEAR_SOLVER_PRECONDITIONER_IFPACK_HPP

#include "baci_config.hpp"

#include "baci_linear_solver_preconditioner_type.hpp"

#include <Ifpack.h>

BACI_NAMESPACE_OPEN

namespace CORE::LINEAR_SOLVER
{
  /*! \brief  IFPACK preconditioners
   *
   *  Set of standard single-matrix preconditioners.
   */
  class IFPACKPreconditioner : public LINEAR_SOLVER::PreconditionerType
  {
   public:
    //! Constructor (empty)
    IFPACKPreconditioner(Teuchos::ParameterList& ifpacklist, Teuchos::ParameterList& solverlist);

    //! Setup
    void Setup(bool create, Epetra_Operator* matrix, Epetra_MultiVector* x,
        Epetra_MultiVector* b) override;

    /// linear operator used for preconditioning
    Teuchos::RCP<Epetra_Operator> PrecOperator() const override { return prec_; }

    //! return name of sublist in parameter list which contains parameters for preconditioner
    std::string getParameterListName() const override { return "IFPACK Parameters"; }

   private:
    //! IFPACK parameter list
    Teuchos::ParameterList& ifpacklist_;

    //! solver parameter list
    Teuchos::ParameterList& solverlist_;

    //! system of equations used for preconditioning used by P_ only
    Teuchos::RCP<Epetra_RowMatrix> Pmatrix_;

    //! preconditioner
    Teuchos::RCP<Ifpack_Preconditioner> prec_;

  };  // class IFPACKPreconditioner
}  // namespace CORE::LINEAR_SOLVER

BACI_NAMESPACE_CLOSE

#endif