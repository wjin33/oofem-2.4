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

#include "admnl1.h"
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
REGISTER_Material(ADNLMaterial);

ADNLMaterial :: ADNLMaterial(int n, Domain *d) : AnisotropicDamageMaterial1(n, d), StructuralNonlocalMaterialExtensionInterface(d)
    //
    // constructor
    //
{ }


ADNLMaterial :: ~ADNLMaterial()
//
// destructor
//
{ }

void
ADNLMaterial :: updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep)
{
    /* Implements the service updating local variables in given integration points,
     * which take part in nonlocal average process. Actually, no update is necessary,
     * because the value used for nonlocal averaging is strain vector used for nonlocal secant stiffness
     * computation. It is therefore necessary only to store local strain in corresponding status.
     * This service is declared at StructuralNonlocalMaterial level.
     */
    FloatArray SDstrainVector;
    double equivStrain;
    ADNLMaterialStatus *nlstatus = static_cast< ADNLMaterialStatus * >( this->giveStatus(gp) );

    this->initTempStatus(gp);

    // subtract stress-independent part
    // note: eigenStrains (temperature) is not contained in mechanical strain stored in gp
    // therefore it is necessary to subtract always the total eigenstrain value
    this->giveStressDependentPartOfStrainVector(SDstrainVector, gp, strainVector, tStep, VM_Total);

    // compute and store the local variable to be averaged
    // (typically the local equivalent strain)
    nlstatus->letTempStrainVectorBe(SDstrainVector);
    this->computeLocalEquivalentStrain(equivStrain, SDstrainVector, gp, tStep);

    // standard formulation based on averaging of equivalent strain
    nlstatus->setLocalEquivalentStrainForAverage(equivStrain);
    return;
}


void
ADNLMaterial :: computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    double nonlocalContribution, nonlocalEquivalentStrain = 0.0;
    ADNLMaterialStatus *nonlocStatus, *status = static_cast< ADNLMaterialStatus * >( this->giveStatus(gp) );

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
        nonlocStatus = static_cast< ADNLMaterialStatus * >( neargp->giveMaterialStatus() );
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
        double localEquivalentStrain = status->giveLocalEquivalentStrainForAverage();
        if ( mm >= 0. ) { // harmonic averaging
            if ( localEquivalentStrain > 0. && nonlocalEquivalentStrain > 0. ) {
                nonlocalEquivalentStrain = 1. / ( mm / nonlocalEquivalentStrain + ( 1. - mm ) / localEquivalentStrain );
            } else {
                nonlocalEquivalentStrain = 0.;
            }
        } else {   // arithmetic averaging, -mm is used instead of mm
            nonlocalEquivalentStrain = -mm * nonlocalEquivalentStrain + ( 1. + mm ) * localEquivalentStrain;
        }
    }

    this->endIPNonlocalAverage(gp);  // I don't know this line is needed?

    kappa = nonlocalEquivalentStrain;

    status->setNonLocalEquivalentStrain(nonlocalEquivalentStrain);
}

Interface *
ADNLMaterial :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialExtensionInterfaceType ) {
        return static_cast< StructuralNonlocalMaterialExtensionInterface * >(this);
    } else {
        return NULL;
    }
}

void
ADNLMaterial :: giveInputRecord(DynamicInputRecord &input)
{
    AnisotropicDamageMaterial1 :: giveInputRecord(input);
    NonlocalMaterialExtensionInterface :: giveInputRecord(input);
}

IRResultType
ADNLMaterial :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    result = AnisotropicDamageMaterial1 :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }
    result = StructuralNonlocalMaterialExtensionInterface :: initializeFrom(ir);
    if ( result != IRRT_OK ) {
        return result;
    }
    return IRRT_OK;
}


int
ADNLMaterial :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
  ADNLMaterialStatus *status = static_cast< ADNLMaterialStatus * >( this->giveStatus(gp) );
  if ( type == IST_LocalEquivalentStrain ) {
      answer.resize(1);
      answer.at(1) = status->giveLocalEquivalentStrainForAverage();
  }else if (type == IST_NonlocalEquivalentStrain){
      answer.resize(1);
      answer.at(1) = status->giveNonLocalEquivalentStrain();
  }else {
    return AnisotropicDamageMaterial1 :: giveIPValue(answer, gp, type, tStep);
  }

  return 1; // to make the compiler happy

}

double
ADNLMaterial :: predictRelativeComputationalCost(GaussPoint *gp)
{
    //
    // The values returned come from measurement
    // do not change them unless you know what are you doing
    //
    double cost = 1.2;


    if ( gp->giveMaterialMode() == _3dMat ) {
        cost = 1.5;
    }

    ADNLMaterialStatus *status = static_cast< ADNLMaterialStatus * >( this->giveStatus(gp) );
    int size = status->giveIntegrationDomainList()->size();
    // just a guess (size/10) found optimal
    // cost *= (1.0 + (size/10)*0.5);
    cost *= ( 1.0 + size / 15.0 );

    return cost;
}


/* *********************************************************************************************   */




ADNLMaterialStatus :: ADNLMaterialStatus(int n, Domain *d, GaussPoint *g) :
AnisotropicDamageMaterial1Status(n, d, g), StructuralNonlocalMaterialStatusExtensionInterface()
{
    localEquivalentStrainForAverage = 0.0;
}


ADNLMaterialStatus :: ~ADNLMaterialStatus()
{ }


void
ADNLMaterialStatus :: printOutputAt(FILE *file, TimeStep *tStep)
{
    MaterialMode mode = gp->giveMaterialMode();
    double Volume = gp->givegausspointvolume();
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
    fprintf(file, "              Damage");
    fprintf( file, " %.4e %.4e %.4e",damage.at(1,1),damage.at(2,2),damage.at(1,2));
    //fprintf(file, "\n              Nl-EqEps-t %f, Nl-EqEps-c %f ", this->equivalenstrain_t, this->equivalenstrain_c);
#ifdef keep_track_of_dissipated_energy
    fprintf(file, "\n              dissW %f freeE %f stressW %f \n", this->dissWork*Volume, ( this->stressWork - this->dissWork )*Volume, this->stressWork*Volume);
#endif
}

void
ADNLMaterialStatus :: initTempStatus()
//
// initializes temp variables according to variables form previous equilibrium state.
// builds new crackMap
//
{
    AnisotropicDamageMaterial1Status :: initTempStatus();

}


void
ADNLMaterialStatus :: updateYourself(TimeStep *tStep)
//
// updates variables (nonTemp variables describing situation at previous equilibrium state)
// after a new equilibrium state has been reached
// temporary variables are having values corresponding to newly reched equilibrium.
//
{
    AnisotropicDamageMaterial1Status :: updateYourself(tStep);
}


Interface *
ADNLMaterialStatus :: giveInterface(InterfaceType type)
{
    if ( type == NonlocalMaterialStatusExtensionInterfaceType ) {
        return static_cast< StructuralNonlocalMaterialStatusExtensionInterface * >(this);
    } else {
        return AnisotropicDamageMaterial1Status :: giveInterface(type);
    }
}

/// Functions for MaterialStatusMapperInterface
void ADNLMaterialStatus ::  copyStateVariables(const MaterialStatus &iStatus)
{
    MaterialStatus &tmpStat = const_cast< MaterialStatus & >(iStatus);
    const ADNLMaterialStatus &nlStatus = dynamic_cast< ADNLMaterialStatus & >(tmpStat);

    this->localEquivalentStrainForAverage = nlStatus.giveLocalEquivalentStrainForAverage();
    this->nonlocalEquivalentStrain =  nlStatus.giveNonLocalEquivalentStrain();

    AnisotropicDamageMaterial1Status :: copyStateVariables(iStatus);


}

} // end namespace oofem
