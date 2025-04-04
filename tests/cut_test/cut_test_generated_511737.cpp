// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_cut_combintersection.hpp"
#include "4C_cut_levelsetintersection.hpp"
#include "4C_cut_meshintersection.hpp"
#include "4C_cut_options.hpp"
#include "4C_cut_side.hpp"
#include "4C_cut_tetmeshintersection.hpp"
#include "4C_cut_volumecell.hpp"
#include "4C_fem_general_utils_local_connectivity_matrices.hpp"

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include "cut_test_utils.hpp"

void test_generated_511737()
{
  Cut::MeshIntersection intersection;
  intersection.GetOptions().Init_for_Cuttests();  // use full cln
  std::vector<int> nids;

  int sidecount = 0;
  std::vector<double> lsvs(8);
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.28396741972129468934;
    tri3_xyze(1, 0) = -0.0033120421506342137205;
    tri3_xyze(2, 0) = -0.60999998945792144323;
    nids.push_back(18884);
    tri3_xyze(0, 1) = -0.29999999991185211101;
    tri3_xyze(1, 1) = 4.4804108985455043434e-09;
    tri3_xyze(2, 1) = -0.60999998982646530532;
    nids.push_back(18886);
    tri3_xyze(0, 2) = -0.29149274647931572302;
    tri3_xyze(1, 2) = -0.012406562203330070981;
    tri3_xyze(2, 2) = -0.60413566292079168285;
    nids.push_back(-2058);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.29999999991185211101;
    tri3_xyze(1, 0) = 4.4804108985455043434e-09;
    tri3_xyze(2, 0) = -0.60999998982646530532;
    nids.push_back(18886);
    tri3_xyze(0, 1) = -0.29911429578782666727;
    tri3_xyze(1, 1) = -0.022083329691285034924;
    tri3_xyze(2, 1) = -0.59827110838966446327;
    nids.push_back(18525);
    tri3_xyze(0, 2) = -0.29149274647931572302;
    tri3_xyze(1, 2) = -0.012406562203330070981;
    tri3_xyze(2, 2) = -0.60413566292079168285;
    nids.push_back(-2058);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.29999999991185211101;
    tri3_xyze(1, 0) = 4.4804108985455043434e-09;
    tri3_xyze(2, 0) = -0.60999998982646530532;
    nids.push_back(18886);
    tri3_xyze(0, 1) = -0.31557688599351790826;
    tri3_xyze(1, 1) = -0.0034890170675985938964;
    tri3_xyze(2, 1) = -0.60999999006249572275;
    nids.push_back(18888);
    tri3_xyze(0, 2) = -0.30727350317533619339;
    tri3_xyze(1, 2) = -0.013063178333627072092;
    tri3_xyze(2, 2) = -0.60413564025141719416;
    nids.push_back(-2059);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix tri3_xyze(3, 3);

    nids.clear();
    tri3_xyze(0, 0) = -0.29911429578782666727;
    tri3_xyze(1, 0) = -0.022083329691285034924;
    tri3_xyze(2, 0) = -0.59827110838966446327;
    nids.push_back(18525);
    tri3_xyze(0, 1) = -0.29999999991185211101;
    tri3_xyze(1, 1) = 4.4804108985455043434e-09;
    tri3_xyze(2, 1) = -0.60999998982646530532;
    nids.push_back(18886);
    tri3_xyze(0, 2) = -0.30727350317533619339;
    tri3_xyze(1, 2) = -0.013063178333627072092;
    tri3_xyze(2, 2) = -0.60413564025141719416;
    nids.push_back(-2059);
    intersection.AddCutSide(++sidecount, nids, tri3_xyze, Core::FE::CellType::tri3);
  }
  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = -0.30047066068049621546;
    hex8_xyze(1, 0) = 0.010854792244467368484;
    hex8_xyze(2, 0) = -0.59026315789473871032;
    nids.push_back(522674);
    hex8_xyze(0, 1) = -0.30066666666666663765;
    hex8_xyze(1, 1) = -1.3586452032633526686e-16;
    hex8_xyze(2, 1) = -0.59026315789473871032;
    nids.push_back(522675);
    hex8_xyze(0, 2) = -0.2933333333333332793;
    hex8_xyze(1, 2) = -1.3264801666129576442e-16;
    hex8_xyze(2, 2) = -0.59026315789473871032;
    nids.push_back(522704);
    hex8_xyze(0, 3) = -0.29314210798097189992;
    hex8_xyze(1, 3) = 0.01059004121411450533;
    hex8_xyze(2, 3) = -0.59026315789473871032;
    nids.push_back(522703);
    hex8_xyze(0, 4) = -0.30047066068049621546;
    hex8_xyze(1, 4) = 0.01085479224446742573;
    hex8_xyze(2, 4) = -0.57052631578947510249;
    nids.push_back(510886);
    hex8_xyze(0, 5) = -0.30066666666666663765;
    hex8_xyze(1, 5) = -7.895761457417655958e-17;
    hex8_xyze(2, 5) = -0.57052631578947510249;
    nids.push_back(510887);
    hex8_xyze(0, 6) = -0.2933333333333332793;
    hex8_xyze(1, 6) = -7.5741110909137057138e-17;
    hex8_xyze(2, 6) = -0.57052631578947510249;
    nids.push_back(510916);
    hex8_xyze(0, 7) = -0.29314210798097189992;
    hex8_xyze(1, 7) = 0.010590041214114562576;
    hex8_xyze(2, 7) = -0.57052631578947510249;
    nids.push_back(510915);

    intersection.add_element(500037, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = -0.30066666666666669316;
    hex8_xyze(1, 0) = -1.8084769459994069553e-16;
    hex8_xyze(2, 0) = -0.60999999999999987566;
    nids.push_back(534463);
    hex8_xyze(0, 1) = -0.30047066068049627097;
    hex8_xyze(1, 1) = -0.010854792244467692877;
    hex8_xyze(2, 1) = -0.60999999999999987566;
    nids.push_back(536399);
    hex8_xyze(0, 2) = -0.29314210798097195543;
    hex8_xyze(1, 2) = -0.010590041214114822785;
    hex8_xyze(2, 2) = -0.60999999999999987566;
    nids.push_back(536428);
    hex8_xyze(0, 3) = -0.29333333333333333481;
    hex8_xyze(1, 3) = -1.7763119093490121774e-16;
    hex8_xyze(2, 3) = -0.60999999999999987566;
    nids.push_back(534492);
    hex8_xyze(0, 4) = -0.30066666666666663765;
    hex8_xyze(1, 4) = -1.3586452032633526686e-16;
    hex8_xyze(2, 4) = -0.59026315789473871032;
    nids.push_back(522675);
    hex8_xyze(0, 5) = -0.30047066068049621546;
    hex8_xyze(1, 5) = -0.010854792244467652979;
    hex8_xyze(2, 5) = -0.59026315789473871032;
    nids.push_back(524611);
    hex8_xyze(0, 6) = -0.29314210798097189992;
    hex8_xyze(1, 6) = -0.010590041214114782886;
    hex8_xyze(2, 6) = -0.59026315789473871032;
    nids.push_back(524640);
    hex8_xyze(0, 7) = -0.2933333333333332793;
    hex8_xyze(1, 7) = -1.3264801666129576442e-16;
    hex8_xyze(2, 7) = -0.59026315789473871032;
    nids.push_back(522704);

    intersection.add_element(513659, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = -0.30779921338002058651;
    hex8_xyze(1, 0) = 0.011119543274820191739;
    hex8_xyze(2, 0) = -0.60999999999999987566;
    nids.push_back(534433);
    hex8_xyze(0, 1) = -0.30800000000000005151;
    hex8_xyze(1, 1) = -1.8406419826498017332e-16;
    hex8_xyze(2, 1) = -0.60999999999999987566;
    nids.push_back(534434);
    hex8_xyze(0, 2) = -0.30066666666666669316;
    hex8_xyze(1, 2) = -1.8084769459994069553e-16;
    hex8_xyze(2, 2) = -0.60999999999999987566;
    nids.push_back(534463);
    hex8_xyze(0, 3) = -0.30047066068049627097;
    hex8_xyze(1, 3) = 0.010854792244467328585;
    hex8_xyze(2, 3) = -0.60999999999999987566;
    nids.push_back(534462);
    hex8_xyze(0, 4) = -0.30779921338002053099;
    hex8_xyze(1, 4) = 0.011119543274820231638;
    hex8_xyze(2, 4) = -0.59026315789473871032;
    nids.push_back(522645);
    hex8_xyze(0, 5) = -0.307999999999999996;
    hex8_xyze(1, 5) = -1.3908102399137472e-16;
    hex8_xyze(2, 5) = -0.59026315789473871032;
    nids.push_back(522646);
    hex8_xyze(0, 6) = -0.30066666666666663765;
    hex8_xyze(1, 6) = -1.3586452032633526686e-16;
    hex8_xyze(2, 6) = -0.59026315789473871032;
    nids.push_back(522675);
    hex8_xyze(0, 7) = -0.30047066068049621546;
    hex8_xyze(1, 7) = 0.010854792244467368484;
    hex8_xyze(2, 7) = -0.59026315789473871032;
    nids.push_back(522674);

    intersection.add_element(511708, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = -0.29988289827639691421;
    hex8_xyze(1, 0) = 0.021695431910718442764;
    hex8_xyze(2, 0) = -0.60999999999999987566;
    nids.push_back(534461);
    hex8_xyze(0, 1) = -0.30047066068049627097;
    hex8_xyze(1, 1) = 0.010854792244467328585;
    hex8_xyze(2, 1) = -0.60999999999999987566;
    nids.push_back(534462);
    hex8_xyze(0, 2) = -0.29314210798097195543;
    hex8_xyze(1, 2) = 0.010590041214114465432;
    hex8_xyze(2, 2) = -0.60999999999999987566;
    nids.push_back(534491);
    hex8_xyze(0, 3) = -0.29256868124526524966;
    hex8_xyze(1, 3) = 0.021166275034847258779;
    hex8_xyze(2, 3) = -0.60999999999999987566;
    nids.push_back(534490);
    hex8_xyze(0, 4) = -0.29988289827639685869;
    hex8_xyze(1, 4) = 0.021695431910718494806;
    hex8_xyze(2, 4) = -0.59026315789473871032;
    nids.push_back(522673);
    hex8_xyze(0, 5) = -0.30047066068049621546;
    hex8_xyze(1, 5) = 0.010854792244467368484;
    hex8_xyze(2, 5) = -0.59026315789473871032;
    nids.push_back(522674);
    hex8_xyze(0, 6) = -0.29314210798097189992;
    hex8_xyze(1, 6) = 0.01059004121411450533;
    hex8_xyze(2, 6) = -0.59026315789473871032;
    nids.push_back(522703);
    hex8_xyze(0, 7) = -0.29256868124526519415;
    hex8_xyze(1, 7) = 0.021166275034847310821;
    hex8_xyze(2, 7) = -0.59026315789473871032;
    nids.push_back(522702);

    intersection.add_element(511736, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = -0.30047066068049627097;
    hex8_xyze(1, 0) = 0.010854792244467328585;
    hex8_xyze(2, 0) = -0.60999999999999987566;
    nids.push_back(534462);
    hex8_xyze(0, 1) = -0.30066666666666669316;
    hex8_xyze(1, 1) = -1.8084769459994069553e-16;
    hex8_xyze(2, 1) = -0.60999999999999987566;
    nids.push_back(534463);
    hex8_xyze(0, 2) = -0.29333333333333333481;
    hex8_xyze(1, 2) = -1.7763119093490121774e-16;
    hex8_xyze(2, 2) = -0.60999999999999987566;
    nids.push_back(534492);
    hex8_xyze(0, 3) = -0.29314210798097195543;
    hex8_xyze(1, 3) = 0.010590041214114465432;
    hex8_xyze(2, 3) = -0.60999999999999987566;
    nids.push_back(534491);
    hex8_xyze(0, 4) = -0.30047066068049621546;
    hex8_xyze(1, 4) = 0.010854792244467368484;
    hex8_xyze(2, 4) = -0.59026315789473871032;
    nids.push_back(522674);
    hex8_xyze(0, 5) = -0.30066666666666663765;
    hex8_xyze(1, 5) = -1.3586452032633526686e-16;
    hex8_xyze(2, 5) = -0.59026315789473871032;
    nids.push_back(522675);
    hex8_xyze(0, 6) = -0.2933333333333332793;
    hex8_xyze(1, 6) = -1.3264801666129576442e-16;
    hex8_xyze(2, 6) = -0.59026315789473871032;
    nids.push_back(522704);
    hex8_xyze(0, 7) = -0.29314210798097189992;
    hex8_xyze(1, 7) = 0.01059004121411450533;
    hex8_xyze(2, 7) = -0.59026315789473871032;
    nids.push_back(522703);

    intersection.add_element(511737, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  {
    Core::LinAlg::SerialDenseMatrix hex8_xyze(3, 8);

    nids.clear();
    hex8_xyze(0, 0) = -0.29314210798097195543;
    hex8_xyze(1, 0) = 0.010590041214114465432;
    hex8_xyze(2, 0) = -0.60999999999999987566;
    nids.push_back(534491);
    hex8_xyze(0, 1) = -0.29333333333333333481;
    hex8_xyze(1, 1) = -1.7763119093490121774e-16;
    hex8_xyze(2, 1) = -0.60999999999999987566;
    nids.push_back(534492);
    hex8_xyze(0, 2) = -0.28600000000000003197;
    hex8_xyze(1, 2) = -1.744146872698617646e-16;
    hex8_xyze(2, 2) = -0.60999999999999987566;
    nids.push_back(534521);
    hex8_xyze(0, 3) = -0.28581355528144763989;
    hex8_xyze(1, 3) = 0.010325290183761602278;
    hex8_xyze(2, 3) = -0.60999999999999987566;
    nids.push_back(534520);
    hex8_xyze(0, 4) = -0.29314210798097189992;
    hex8_xyze(1, 4) = 0.01059004121411450533;
    hex8_xyze(2, 4) = -0.59026315789473871032;
    nids.push_back(522703);
    hex8_xyze(0, 5) = -0.2933333333333332793;
    hex8_xyze(1, 5) = -1.3264801666129576442e-16;
    hex8_xyze(2, 5) = -0.59026315789473871032;
    nids.push_back(522704);
    hex8_xyze(0, 6) = -0.28599999999999997646;
    hex8_xyze(1, 6) = -1.2943151299625631128e-16;
    hex8_xyze(2, 6) = -0.59026315789473871032;
    nids.push_back(522733);
    hex8_xyze(0, 7) = -0.28581355528144758438;
    hex8_xyze(1, 7) = 0.010325290183761642177;
    hex8_xyze(2, 7) = -0.59026315789473871032;
    nids.push_back(522732);

    intersection.add_element(511766, nids, hex8_xyze, Core::FE::CellType::hex8);
  }

  intersection.CutTest_Cut(
      true, Cut::VCellGaussPts_DirectDivergence, Cut::BCellGaussPts_Tessellation);
  intersection.Cut_Finalize(
      true, Cut::VCellGaussPts_DirectDivergence, Cut::BCellGaussPts_Tessellation, false, true);
}
