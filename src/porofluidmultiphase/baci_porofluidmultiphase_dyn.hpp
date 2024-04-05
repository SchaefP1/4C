/*----------------------------------------------------------------------*/
/*! \file
 \brief entry point (global control routine) for porous multiphase flow

   \level 3

 *----------------------------------------------------------------------*/

#ifndef BACI_POROFLUIDMULTIPHASE_DYN_HPP
#define BACI_POROFLUIDMULTIPHASE_DYN_HPP

#include "baci_config.hpp"

BACI_NAMESPACE_OPEN


/*! entry point for the solution of Lubrication problems */
void porofluidmultiphase_dyn(int restart /* do we have to perform a restart?  */
);


BACI_NAMESPACE_CLOSE

#endif  // POROFLUIDMULTIPHASE_DYN_H