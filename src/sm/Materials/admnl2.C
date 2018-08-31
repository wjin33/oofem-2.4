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

#include "admnl2.h"
#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "../sm/Elements/structuralelement.h"
#include "sparsemtrx.h"
#include "error.h"
#include "nonlocalmaterialext.h"
#include "contextioerr.h"
#include "stressvector.h"
#include "strainvector.h"
#include "classfactory.h"
#include "dynamicinputrecord.h"
#include "datastream.h"
#include "unknownnumberingscheme.h"

#ifdef __OOFEG
 #include "oofeggraphiccontext.h"
 #include "connectivitytable.h"
#endif

namespace oofem {
REGISTER_Material(ADNLMaterial2)

ADNLMaterial2 :: ADNLMaterial2(int n, Domain *d) : AnisotropicDamageMaterial2(n, d), StructuralNonlocalMaterialExtensionInterface(d)
    //
    // constructor
    //
{ }


ADNLMaterial2 :: ~ADNLMaterial2()
//
// destructor
//
{ }

void
ADNLMaterial2 :: updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep)
{
    /* Implements the service updating local variables in given integration points,
     * which take part in nonlocal average process. Actually, no update is necessary,
     * because the value used for nonlocal averaging is strain vector used for nonlocal secant stiffness
     * computation. It is therefore necessary only to store local strain in corresponding status.
     * This service is declared at StructuralNonlocalMaterial level.
     */
    FloatArray SDstrainVector;

    ADNLMaterial2Status *nlstatus = static_cast< ADNLMaterial2Status * >( this->giveStatus(gp) );

    //this->initTempStatus(gp);

    // subtract stress-independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigenstrain value
    this->giveStressDependentPartOfStrainVector(SDstrainVector, gp, strainVector, tStep, VM_Total);

    // compute and store the local variable to be averaged
    // (typically the local equivalent strain)
    nlstatus->letTempStrainVectorBe(SDstrainVector);

    FloatMatrix Trans(4,4);
    FloatArray GStrain(4),LocStrain(4);
    GStrain.at(1)=SDstrainVector.at(1);
    GStrain.at(2)=SDstrainVector.at(2);
    GStrain.at(3)=SDstrainVector.at(3);
    GStrain.at(4)=SDstrainVector.at(4)/2;

    Trans.zero();
    Trans.at(1,1)=cos(this->theta)*cos(this->theta);
    Trans.at(1,2)=sin(this->theta)*sin(this->theta);
    Trans.at(1,3)=0;
    Trans.at(1,4)=2*sin(this->theta)*cos(this->theta);
    Trans.at(2,1)=sin(this->theta)*sin(this->theta);
    Trans.at(2,2)=cos(this->theta)*cos(this->theta);
    Trans.at(2,3)=0;
    Trans.at(2,4)=-2*sin(this->theta)*cos(this->theta);
    Trans.at(3,1)=0;
    Trans.at(3,2)=0;
    Trans.at(3,3)=1;
    Trans.at(3,4)=0;
    Trans.at(4,1)=-sin(this->theta)*cos(this->theta);
    Trans.at(4,2)=sin(this->theta)*cos(this->theta);
    Trans.at(4,3)=0;
    Trans.at(4,4)=cos(this->theta)*cos(this->theta)-sin(this->theta)*sin(this->theta);


    LocStrain.zero();
    LocStrain.beProductOf(Trans, GStrain);

    FloatArray TempStress;
    TempStress=nlstatus->giveTempStressVector();
    double TrSig;
    TrSig = TempStress.at(1) + TempStress.at(2) + TempStress.at(3);

    FloatArray equivStrain(4);

    this->computeLocalEquivalentStrain(equivStrain, LocStrain,TrSig, gp, tStep);

    // standard formulation based on averaging of equivalent strain
    nlstatus->setLocalEquivalentStrainForAverage(equivStrain);
    return;
}


void
ADNLMaterial2 :: computeEquivalentStrain(FloatArray &answer, const FloatArray &strain,const double &TrSig, GaussPoint *gp, TimeStep *tStep)
{
    FloatArray nonlocalContribution(4), nonlocalEquivalentStrain(4);
    nonlocalContribution.zero();
    nonlocalEquivalentStrain.zero();

    ADNLMaterial2Status *nonlocStatus, *status = static_cast< ADNLMaterial2Status * >( this->giveStatus(gp) );

    this->buildNonlocalPointTable(gp);
    this->updateDomainBeforeNonlocAverage(tStep);

    // compute nonlocal equivalent strain

    std :: list< localIntegrationRecord > *list = this->giveIPIntegrationList(gp); // !

    //double sigmaRatio = 0.; //ratio sigma2/sigma1 used for stress-based averaging
    //double nx, ny; //components of the first principal stress direction (for stress-based averaging)
    //double updatedIntegrationVolume = 0.; //new integration volume. Sum of all new weights used for stress-based averaging


    //Loop over all Gauss points which are in gp's integration domain
    for ( auto &lir : *list ) {
        GaussPoint *neargp = lir.nearGp;
        nonlocStatus = static_cast< ADNLMaterial2Status * >( neargp->giveMaterialStatus() );
        nonlocalContribution = nonlocStatus->giveLocalEquivalentStrainForAverage();
        nonlocalContribution *= lir.weight;
        nonlocalEquivalentStrain += nonlocalContribution;
    }

    if ( scaling == ST_Standard ) { // standard rescaling
        nonlocalEquivalentStrain *= 1. / status->giveIntegrationScale();
    } else if ( scaling == ST_Borino ) { // Borino modification
        double scale = status->giveIntegrationScale();
        if ( scale > 1. ) {
            nonlocalEquivalentStrain *= 1. / scale;
        } else {
            nonlocalEquivalentStrain += ( 1. - scale ) * status->giveLocalEquivalentStrainForAverage();
        }
    }

    // undernonlocal or overnonlocal formulation
    if ( mm != 1. ) {
        FloatArray localEquivalentStrain = status->giveLocalEquivalentStrainForAverage();
        if ( mm >= 0. ) { // harmonic averaging
            if ( localEquivalentStrain.at(1) > 0. && nonlocalEquivalentStrain.at(1) > 0. ) {
                nonlocalEquivalentStrain.at(1) = 1./ ( mm / nonlocalEquivalentStrain.at(1) + ( 1. - mm ) / localEquivalentStrain.at(1) );
            } else {
                nonlocalEquivalentStrain.at(1) = 0;
            }
            if ( localEquivalentStrain.at(2) > 0. && nonlocalEquivalentStrain.at(2) > 0. ) {
                nonlocalEquivalentStrain.at(2) = 1./ ( mm / nonlocalEquivalentStrain.at(2) + ( 1. - mm ) / localEquivalentStrain.at(2) );
            } else {
                nonlocalEquivalentStrain.at(2) = 0;
            }
            if ( localEquivalentStrain.at(3) > 0. && nonlocalEquivalentStrain.at(3) > 0. ) {
                nonlocalEquivalentStrain.at(3) = 1./ ( mm / nonlocalEquivalentStrain.at(3) + ( 1. - mm ) / localEquivalentStrain.at(3) );
            } else {
                nonlocalEquivalentStrain.at(3) = 0;
            }
            if ( localEquivalentStrain.at(4) > 0. && nonlocalEquivalentStrain.at(4) > 0. ) {
                nonlocalEquivalentStrain.at(4) = 1./ ( mm / nonlocalEquivalentStrain.at(4) + ( 1. - mm ) / localEquivalentStrain.at(4) );
            } else {
                nonlocalEquivalentStrain.at(4) = 0;
            }
        } else {   // arithmetic averaging, -mm is used instead of mm
            nonlocalEquivalentStrain.at(1) = -mm * nonlocalEquivalentStrain.at(1) + ( 1. + mm ) * localEquivalentStrain.at(1);
            nonlocalEquivalentStrain.at(2) = -mm * nonlocalEquivalentStrain.at(2) + ( 1. + mm ) * localEquivalentStrain.at(2);
            nonlocalEquivalentStrain.at(3) = -mm * nonlocalEquivalentStrain.at(3) + ( 1. + mm ) * localEquivalentStrain.at(3);
            nonlocalEquivalentStrain.at(4) = -mm * nonlocalEquivalentStrain.at(4) + ( 1. + mm ) * localEquivalentStrain.at(4);
        }
    }

    this->endIPNonlocalAverage(gp);  // I don't know whether this line is needed or not?

    answer = nonlocalEquivalentStrain;
}

Interface *
ADNLMaterial2 :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialExtensionInterfaceType ) {
        return static_cast< StructuralNonlocalMaterialExtensionInterface * >(this);
    } else {
        return NULL;
    }
}


IRResultType
ADNLMaterial2 :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = AnisotropicDamageMaterial2 :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }
    result = StructuralNonlocalMaterialExtensionInterface :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }
/*
*    averType = 0;
*    IR_GIVE_OPTIONAL_FIELD(ir, averType, _IFT_ADNLMaterial2_averagingtype);
*    if ( averType == 2 ) {
*        exponent = 0.5; // default value for averaging type 2
*    }
*
*    if ( averType == 3 ) {
*        exponent = 1.; // default value for averaging type 3
*    }
*
*    if ( averType == 2 || averType == 3 ) {
*        IR_GIVE_OPTIONAL_FIELD(ir, exponent, _IFT_ADNLMaterial2_exp);
*    }
*
*    if ( averType >= 2 && averType <= 5 ) {
*        IR_GIVE_OPTIONAL_FIELD(ir, Rf, _IFT_ADNLMaterial2_rf);
*    }
*/
    return IRRT_OK;
}


int
ADNLMaterial2 :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
  ADNLMaterial2Status *status = static_cast< ADNLMaterial2Status * >( this->giveStatus(gp) );
  if ( type == IST_LocalEquivalentStrain ) {
    answer.resize(4);
    answer.zero();
    answer  = status->giveLocalEquivalentStrainForAverage();
  } else {
    return AnisotropicDamageMaterial2 :: giveIPValue(answer, gp, type, tStep);
  }

  return 1; // to make the compiler happy

}

double
ADNLMaterial2 :: predictRelativeComputationalCost(GaussPoint *gp)
{
    //
    // The values returned come from measurement
    // do not change them unless you know what are you doing
    //
    double cost = 1.2;


    if ( gp->giveMaterialMode() == _3dMat ) {
        cost = 1.5;
    }

    ADNLMaterial2Status *status = static_cast< ADNLMaterial2Status * >( this->giveStatus(gp) );
    int size = status->giveIntegrationDomainList()->size();
    // just a guess (size/10) found optimal
    // cost *= (1.0 + (size/10)*0.5);
    cost *= ( 1.0 + size / 15.0 );

    return cost;
}


/* *********************************************************************************************   */




ADNLMaterial2Status :: ADNLMaterial2Status(int n, Domain *d, GaussPoint *g) :
AnisotropicDamageMaterial2Status(n, d, g), StructuralNonlocalMaterialStatusExtensionInterface()
{
    localEquivalentStrainForAverage.zero();
}


ADNLMaterial2Status :: ~ADNLMaterial2Status()
{ }


void
ADNLMaterial2Status :: printOutputAt(FILE *file, TimeStep *tStep)
{
    MaterialMode mode = gp->giveMaterialMode();
    if ( mode == _PlaneStrain ) { // special treatment of the out-of-plane stress
        FloatArray helpVec;
        MaterialStatus :: printOutputAt(file, tStep);
        fprintf(file, "  strains ");
        StructuralMaterial :: giveFullSymVectorForm(helpVec, strainVector, mode);
        for ( auto &v : helpVec ) {
            fprintf( file, " %.4e", v );
        }
        fprintf(file, "\n              stresses");
        StructuralMaterial :: giveFullSymVectorForm(helpVec, stressVector, mode);
        for ( auto &v : helpVec ) {
            fprintf( file, " %.4e", v );
        }
        fprintf(file, "\n");
    } else {
        StructuralMaterialStatus :: printOutputAt(file, tStep); // standard treatment of strains and stresses
    }
    fprintf(file, "status { ");
    fprintf(file, "\n      Tensile Damage");
    fprintf( file, " %.4e %.4e %.4e",damage.at(1,1),damage.at(1,2),damage.at(1,3));
    fprintf(file, "\n   Compresive Damage");
    fprintf( file, " %.4e %.4e %.4e",damage.at(2,1),damage.at(2,2),damage.at(2,3));
#ifdef keep_track_of_dissipated_energy
    fprintf(file, "\n              dissW %f, freeE %f, stressW %f ", this->dissWork, ( this->stressWork ) - ( this->dissWork ), this->stressWork);
#endif
    fprintf(file, "\n}\n");
}

void
ADNLMaterial2Status :: initTempStatus()
//
// initializes temp variables according to variables form previous equilibrium state.
// builds new crackMap
//
{
    AnisotropicDamageMaterial2Status :: initTempStatus();

}


void
ADNLMaterial2Status :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    AnisotropicDamageMaterial2Status :: updateYourself(tStep);
}


Interface *
ADNLMaterial2Status :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialStatusExtensionInterfaceType ) {
        return static_cast< StructuralNonlocalMaterialStatusExtensionInterface * >(this);
    } else {
        return AnisotropicDamageMaterial2Status :: giveInterface(type);
    }
}

} // end namespace oofem
