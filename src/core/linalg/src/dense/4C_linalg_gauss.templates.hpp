// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_LINALG_GAUSS_TEMPLATES_HPP
#define FOUR_C_LINALG_GAUSS_TEMPLATES_HPP

#include "4C_config.hpp"

#include "4C_linalg_gauss.hpp"
#include "4C_utils_mathoperations.hpp"

FOUR_C_NAMESPACE_OPEN

namespace Core::LinAlg
{
  template <bool do_piv, unsigned dim, typename Valtype>
  Valtype gauss_elimination(Core::LinAlg::Matrix<dim, dim, Valtype>& A,
      Core::LinAlg::Matrix<dim, 1, Valtype>& b, Core::LinAlg::Matrix<dim, 1, Valtype>& x)
  {
    if (dim > 1)
    {
      bool changesign = false;
      if (not do_piv)
      {
        for (unsigned k = 0; k < dim; ++k)
        {
          A(k, k) = 1. / A(k, k);

          for (unsigned i = k + 1; i < dim; ++i)
          {
            A(i, k) *= A(k, k);
            x(i) = A(i, k);

            for (unsigned j = k + 1; j < dim; ++j)
            {
              A(i, j) -= A(i, k) * A(k, j);
            }
          }

          for (unsigned i = k + 1; i < dim; ++i)
          {
            b(i) -= x(i) * b(k);
          }
        }
      }
      else
      {
        for (unsigned k = 0; k < dim; ++k)
        {
          unsigned pivot = k;

          // search for pivot element
          for (unsigned i = k + 1; i < dim; ++i)
          {
            pivot = (Core::MathOperations<Valtype>::abs(A(pivot, k)) <
                        Core::MathOperations<Valtype>::abs(A(i, k)))
                        ? i
                        : pivot;
          }

          // exchange pivot row and current row
          if (pivot != k)
          {
            for (unsigned j = 0; j < dim; ++j)
            {
              std::swap(A(k, j), A(pivot, j));
            }
            std::swap(b(k, 0), b(pivot, 0));
            changesign = not changesign;
          }

          if (A(k, k) == 0.0)
          {
            return 0.0;
          }

          A(k, k) = 1. / A(k, k);

          for (unsigned i = k + 1; i < dim; ++i)
          {
            A(i, k) *= A(k, k);
            x(i, 0) = A(i, k);

            for (unsigned j = k + 1; j < dim; ++j)
            {
              A(i, j) -= A(i, k) * A(k, j);
            }
          }

          for (unsigned i = k + 1; i < dim; ++i)
          {
            b(i, 0) -= x(i, 0) * b(k, 0);
          }
        }
      }

      // back substitution
      x(dim - 1, 0) = b(dim - 1, 0) * A(dim - 1, dim - 1);

      for (int i = dim - 2; i >= 0; --i)
      {
        for (int j = dim - 1; j > i; --j)
        {
          b(i, 0) -= A(i, j) * x(j, 0);
        }
        x(i, 0) = b(i, 0) * A(i, i);
      }
      Valtype det = 1.0;
      for (unsigned i = 0; i < dim; ++i) det *= 1.0 / A(i, i);

      if (changesign) det *= -1.0;

      return det;
    }
    else
    {
      x(0, 0) = b(0, 0) / A(0, 0);
      return x(0, 0);
    }
  }

  /*!
    \brief computes a Gaussian elimination for a linear system of equations after infnorm scaling

    \tparam dim      (in)    : dimension of the matrix
    \return determinant of system matrix
  */
  template <unsigned dim>
  double scaled_gauss_elimination(Core::LinAlg::Matrix<dim, dim>& A,  ///< (in)    : system matrix
      Core::LinAlg::Matrix<dim, 1>& b,                                ///< (in)    : right-hand-side
      Core::LinAlg::Matrix<dim, 1>& x                                 ///< (out)   : solution vector
  )
  {
    // infnorm scaling
    for (unsigned i = 0; i < dim; ++i)
    {
      // find norm of max entry in row
      double max = std::abs(A(i, 0));
      for (unsigned j = 1; j < dim; ++j)
      {
        const double norm = std::abs(A(i, j));
        if (norm > max) max = norm;
      }

      // close to zero row detected -> matrix does probably not have full rank
      if (max < 1.0e-14)
      {
        return Core::LinAlg::gauss_elimination<true, dim>(A, b, x);
      }

      // scale row with inv of max entry
      const double scale = 1.0 / max;
      for (unsigned j = 0; j < dim; ++j)
      {
        A(i, j) *= scale;
      }
      b(i) *= scale;
    }

    // solve scaled system using pivoting
    return Core::LinAlg::gauss_elimination<true, dim>(A, b, x);
  }

}  // namespace Core::LinAlg

FOUR_C_NAMESPACE_CLOSE

#endif
