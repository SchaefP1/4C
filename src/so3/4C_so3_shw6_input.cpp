/*----------------------------------------------------------------------*/
/*! \file
\brief
\level 1


*----------------------------------------------------------------------*/


#include "4C_so3_shw6.hpp"  //**
#include "4C_mat_so3_material.hpp"
#include "4C_io_linedefinition.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ELEMENTS::SoShw6::ReadElement(
    const std::string& eletype, const std::string& distype, INPUT::LineDefinition* linedef)
{
  // read number of material model
  int material = 0;
  linedef->ExtractInt("MAT", material);
  SetMaterial(0, MAT::Factory(material));

  SolidMaterial()->Setup(NUMGPT_WEG6, linedef);

  std::string buffer;
  linedef->ExtractString("KINEM", buffer);


  // geometrically non-linear with Total Lagrangean approach
  if (buffer == "nonlinear")
  {
    kintype_ = INPAR::STR::KinemType::nonlinearTotLag;
  }
  // geometrically linear
  else if (buffer == "linear")
  {
    kintype_ = INPAR::STR::KinemType::linear;
    FOUR_C_THROW("Reading of SOLIDSHW6 element failed onlz nonlinear kinetmatics implemented");
  }

  // geometrically non-linear with Updated Lagrangean approach
  else
    FOUR_C_THROW("Reading of SOLIDSHW6 element failed KINEM unknown");

  // check if material kinematics is compatible to element kinematics
  SolidMaterial()->ValidKinematics(kintype_);

  // Validate that materials doesn't use extended update call.
  if (SolidMaterial()->UsesExtendedUpdate())
    FOUR_C_THROW("This element currently does not support the extended update call.");

  linedef->ExtractString("EAS", buffer);

  // full sohw6 EAS technology
  if (buffer == "soshw6")
  {
    eastype_ = soshw6_easpoisthick;  // EAS to allow linear thickness strain
    neas_ = soshw6_easpoisthick;     // number of eas parameters
    soshw6_easinit();
  }
  // no EAS technology
  else if (buffer == "none")
  {
    eastype_ = soshw6_easnone;
    neas_ = 0;  // number of eas parameters
    std::cout << "Warning: Solid-Shell Wegde6 without EAS" << std::endl;
  }
  else
    FOUR_C_THROW("Reading of SOLIDSHW6 EAS technology failed");

  // check for automatically align material space optimally with parameter space
  optimal_parameterspace_map_ = false;
  nodes_rearranged_ = false;
  if (linedef->HaveNamed("OPTORDER")) optimal_parameterspace_map_ = true;

  return true;
}

FOUR_C_NAMESPACE_CLOSE