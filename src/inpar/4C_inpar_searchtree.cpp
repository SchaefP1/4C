// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_inpar_searchtree.hpp"

#include "4C_utils_parameter_list.hpp"

FOUR_C_NAMESPACE_OPEN



void Inpar::Geo::set_valid_parameters(std::map<std::string, Core::IO::InputSpec>& list)
{
  using Teuchos::tuple;
  using namespace Core::IO::InputSpecBuilders;

  Core::Utils::SectionSpecs search_tree{"SEARCH TREE"};

  search_tree.specs.emplace_back(deprecated_selection<Inpar::Geo::TreeType>("TREE_TYPE",
      {
          {"notree", Inpar::Geo::Notree},
          {"octree3d", Inpar::Geo::Octree3D},
          {"quadtree3d", Inpar::Geo::Quadtree3D},
          {"quadtree2d", Inpar::Geo::Quadtree2D},
      },
      {.description = "set tree type", .default_value = Inpar::Geo::Notree}));

  search_tree.move_into_collection(list);
}

FOUR_C_NAMESPACE_CLOSE
