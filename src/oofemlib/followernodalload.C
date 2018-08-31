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

#include "followernodalload.h"
#include "classfactory.h"
#include "masterdof.h"
#include "floatmatrix.h"
#include "sparsemtrx.h"
#include "unknownnumberingscheme.h"
#include "function.h"
#include "timestep.h"
#include "datastream.h"
#include "contextioerr.h"
#include "node.h"
#include "domain.h"

namespace oofem {
REGISTER_BoundaryCondition(FollowerNodalLoad);

FollowerNodalLoad :: FollowerNodalLoad(int n, Domain *d) : ActiveBoundaryCondition(n, d)
{
    this->node = 0;
}


FollowerNodalLoad :: ~FollowerNodalLoad()
{
}


IRResultType FollowerNodalLoad :: initializeFrom(InputRecord *ir)
{
    ActiveBoundaryCondition :: initializeFrom(ir);
    IRResultType result;
    rhsTf = 0;

    IR_GIVE_FIELD(ir, node, _IFT_FollowerNodalLoad_node);
    IR_GIVE_FIELD(ir, components, _IFT_FollowerNodalLoad_components);
    IR_GIVE_FIELD(ir, ltf, _IFT_FollowerNodalLoad_ltf);

    return IRRT_OK;
}


void FollowerNodalLoad :: giveLocArray(const UnknownNumberingScheme &r_s,  IntArray &locr)
{
  return this->domain->giveDofManager(this->node)->giveCompleteLocationArray(locr, r_s);
}


void LinearConstraintBC :: assembleVector(FloatArray &answer, TimeStep *tStep,
                                          CharType type, ValueModeType mode,
                                          const UnknownNumberingScheme &s, FloatArray *eNorms)
{
    if ( type == ExternalForcesVector ) {

      DofManager *dman = this->domain->giveDofManager(this->node);

      for (Dof *dof: *dman) {
        if (dof->giveDofID() == D_u) {
          
      }



        // compute true residual
        int size = this->weights.giveSize();
        Dof *mdof = *md->begin();
        Dof *idof;

        // assemble location array
        for ( int _i = 1; _i <= size; _i++ ) {
            factor = 1.;
            if ( weightsTf.giveSize() ) {
                factor = domain->giveFunction( weightsTf.at(_i) )->evaluateAtTime( tStep->giveIntrinsicTime() );
            }
            idof = this->domain->giveDofManager( this->dofmans.at(_i) )->giveDofWithID( this->dofs.at(_i) );
            if ( s.giveDofEquationNumber(idof) ) {
                answer.at( s.giveDofEquationNumber(idof) ) += mdof->giveUnknown(mode, tStep) * this->weights.at(_i) * factor;
            }
            if ( s.giveDofEquationNumber( mdof ) ) {
                answer.at( s.giveDofEquationNumber( mdof ) ) += idof->giveUnknown(mode, tStep) * this->weights.at(_i) * factor;
            }
        }
    } else {
        // use rhs value

        if ( rhsTf ) {
            factor = domain->giveFunction(rhsTf)->evaluateAtTime( tStep->giveIntrinsicTime() );
        }
        this->giveLocArray( s, loc, lambdaeq.at(1) );
        vec.at(1) = rhs * factor;
        answer.assemble(vec, lambdaeq);
    }
}

void LinearConstraintBC :: giveLocationArrays(std :: vector< IntArray > &rows, std :: vector< IntArray > &cols, CharType type, const UnknownNumberingScheme &r_s, const UnknownNumberingScheme &c_s)
{
    rows.resize(3);
    cols.resize(3);

    IntArray loc, lambdaeq(1);
    this->giveLocArray( r_s, loc, lambdaeq.at(1) );
    // column block
    rows [ 0 ] = loc;
    cols [ 0 ] = lambdaeq;
    // row block
    cols [ 1 ] = loc;
    rows [ 1 ] = lambdaeq;
    // diagonal enry (some sparse mtrx implementation requaire this)
    rows [ 2 ] = lambdaeq;
    cols [ 2 ] = lambdaeq;
}


contextIOResultType
LinearConstraintBC :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( mode & CM_Definition ) {
        if ( ( iores = weights.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = weightsTf.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = dofmans.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = dofs.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( !stream.write(rhs) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.write(rhsTf) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( ( iores = lhsType.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = rhsType.storeYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    if ( ( iores = md->saveContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}


contextIOResultType
LinearConstraintBC :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    contextIOResultType iores;
    if ( mode & CM_Definition ) {
        if ( ( iores = weights.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = weightsTf.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = dofmans.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = dofs.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( !stream.read(rhs) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( !stream.read(rhsTf) ) {
            THROW_CIOERR(CIO_IOERR);
        }
        if ( ( iores = lhsType.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
        if ( ( iores = rhsType.restoreYourself(stream) ) != CIO_OK ) {
            THROW_CIOERR(iores);
        }
    }

    if ( ( iores = md->restoreContext(stream, mode) ) != CIO_OK ) {
        THROW_CIOERR(iores);
    }

    return CIO_OK;
}
} //end of oofem namespace
