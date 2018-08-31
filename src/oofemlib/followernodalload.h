/*
 *
 *                 #####    #####   ######  ######  ###   ###
 *               ##   ##  ##   ##  ##      ##      ## ### ##
 *              ##   ##  ##   ##  ####    ####    ##  #  ##
 *             ##   ##  ##   ##  ##      ##      ##     ##
 *            ##   ##  ##   ##  ##      ##      ##     ##
 *            #####    #####   ##      ######  ##     ##
 *
 *
 *             OOFEM : Object Oriented Finite Element Code
 *
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This library is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU Lesser General Public
 *  License as published by the Free Software Foundation; either
 *  version 2.1 of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *  Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public
 *  License along with this library; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 */

#ifndef followernodalload_h
#define followernodalload_h

#include "activebc.h"
#include "floatarray.h"
#include "intarray.h"
#include "chartype.h"
#include "valuemodetype.h"
#include "dofmanager.h"
#include "error.h"

#define _IFT_FollowerNodalLoad_Name   "followernodalload"

///@name Input fields for active boundary condition
//@{
#define _IFT_FollowerNodalLoad_node "node"
#define _IFT_FollowerNodalLoad_components "components"
#define _IFT_FollowerNodalLoad_ltf "ltf"
//@}

namespace oofem {
/**
 * Class implementing follower type of nodal load, where the load follows the displacements and rotations 
 * of the associated node. The load components assembled for original undeformed configuration are continuously
 * updated based on displacements and rotations of the associated node.
 * This bc does not reuire any internal DOF
 */
class OOFEM_EXPORT FollowerNodalLoad : public ActiveBoundaryCondition
{
protected:
  int node; // associated node
  FloatArray components; // components defined on deformed configuration
  int ltf; // load time function

public:
    FollowerNodalLoad(int n, Domain * d);
    /// Destructor.
    virtual ~FollowerNodalLoad();

    IRResultType initializeFrom(InputRecord *ir);
    virtual const char *giveInputRecordName() const { return _IFT_FollowerNodalLoad_Name; }
    virtual void assembleVector(FloatArray &answer, TimeStep *tStep,
                                CharType type, ValueModeType mode,
                                const UnknownNumberingScheme &s, FloatArray *eNorms = NULL);

    virtual void giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols,
                                    CharType type, const UnknownNumberingScheme &r_s,
                                    const UnknownNumberingScheme &c_s);

    /// Gives the number of internal dof managers.
    virtual int giveNumberOfInternalDofManagers() { return 0; }

    virtual contextIOResultType saveContext(DataStream &stream, ContextMode mode, void *obj = NULL);
    virtual contextIOResultType restoreContext(DataStream &stream, ContextMode mode, void *obj = NULL);

    virtual const char *giveClassName() const { return "FollowerNodalLoad"; }

protected:
};
} //end of oofem namespace
#endif // followernodalload_h
