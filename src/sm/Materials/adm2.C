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
 *---
 *               Copyright (C) 1993 - 2013   Borek Patzak
 *
 *
 *
 *       Czech Technical University, Faculty of Civil Engineering,
 *   Department of Structural Mechanics, 166 29 Prague, Czech Republic
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

// This code is based on the anisotropic damage model proposed by Desmorat, Gatuingt and Ragueneau in
// their paper "Nonlocal anisotropic damage model and related computational aspects for quasi-brittle material"
// published in Engineering Fracture Mechanics 74 (2007) 1539-1560.

#include "adm2.h"
#include "../sm/Materials/structuralmaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "dynamicinputrecord.h"
#include "gausspoint.h"
#include "classfactory.h"
#include <math.h>

namespace oofem {
REGISTER_Material(AnisotropicDamageMaterial2)

AnisotropicDamageMaterial2 :: AnisotropicDamageMaterial2(int n, Domain *d) : StructuralMaterial(n, d),
    RandomMaterialExtensionInterface()
    //
    // constructor
    //
{
    E11 = 0.;
    E22 = 0.;
    nu12 = 0.;
    nu23 = 0.;
    G12 = 0.;
    eqeps_1t = 0;
    eqeps_1s = 0;
    alpha_1t = 0.;
    eqeps_1c = 0.;
    beta_1c = 0.;
    alpha_1c = 0.;
    eqeps_2t = 0;
    eqeps_2s = 0;
    alpha_2t = 0.;
    eqeps_2c = 0.;
    beta_2c = 0.;
    alpha_2c = 0.;
    eta = 0.0;
    theta = 0.;
}

AnisotropicDamageMaterial2 :: ~AnisotropicDamageMaterial2()
//
// destructor
//
{}

/*
Interface *
AnisotropicDamageMaterial2 :: giveInterface(InterfaceType type)
{
    if ( type == MaterialModelMapperInterfaceType ) {
        return static_cast< MaterialModelMapperInterface * >(this);
    } else {
        return NULL;
    }
}
*/

MaterialStatus *
AnisotropicDamageMaterial2 :: CreateStatus(GaussPoint *gp) const
{
    return new AnisotropicDamageMaterial2Status(1, domain, gp);
}


MaterialStatus *
AnisotropicDamageMaterial2 :: giveStatus(GaussPoint *gp) const
{
    MaterialStatus *status = static_cast< MaterialStatus * >( gp->giveMaterialStatus() );
    if ( status == NULL ) {
        // create a new one
        status = this->CreateStatus(gp);

        if ( status != NULL ) {
            gp->setMaterialStatus( status, this->giveNumber() );
            this->_generateStatusVariables(gp);
        }
    }

    return status;
}

int
AnisotropicDamageMaterial2 :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports the given mode
//
{
    return mode == _PlaneStrain;
    //return mode == _3dMat || mode == _PlaneStress || mode == _PlaneStrain || mode == _1dMat;
}

void
AnisotropicDamageMaterial2 :: giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *atTime)
{
    ///////////////////////////////////////////////////////
    //Does not account for 3D unloading case /////////////
    //The 3D case is not complete  ////////////////////////
    #define TOLERANCE 1.e-6 // convergence tolerance for the internal iteration
    //this->initTempStatus(gp); //updata the status, so the current iteration starts from the previous converged results;
    // subtract the stress-independent part of strains (e.g. due to temperature)
    FloatArray SDstrainVector;

    this->giveStressDependentPartOfStrainVector(SDstrainVector, gp, totalStrain, atTime, VM_Total);

            // make sure strain components are correctly stored as the following--WJ

            FloatMatrix Trans(6,6),TransPrim(6,6);
            FloatArray  GStrain(6),LocStrain(6);
            GStrain.at(1)=SDstrainVector.at(1);
            GStrain.at(2)=SDstrainVector.at(2);
            GStrain.at(3)=SDstrainVector.at(3);
            GStrain.at(4)=SDstrainVector.at(4)/2;
            GStrain.at(5)=SDstrainVector.at(5)/2;
            GStrain.at(6)=SDstrainVector.at(6)/2;
            Trans.zero();
            Trans.at(1,1)=cos(this->theta)*cos(this->theta);
            Trans.at(1,2)=sin(this->theta)*sin(this->theta);
            Trans.at(1,6)=2*sin(this->theta)*cos(this->theta);
            Trans.at(2,1)=sin(this->theta)*sin(this->theta);
            Trans.at(2,2)=cos(this->theta)*cos(this->theta);
            Trans.at(2,6)=-2*sin(this->theta)*cos(this->theta);
            Trans.at(3,3)=1;
            Trans.at(4,4)=cos(this->theta);
            Trans.at(4,5)=-sin(this->theta);
            Trans.at(5,4)=sin(this->theta);
            Trans.at(5,5)=cos(this->theta);
            Trans.at(6,1)=-sin(this->theta)*cos(this->theta);
            Trans.at(6,2)=sin(this->theta)*cos(this->theta);
            Trans.at(6,6)=cos(this->theta)*cos(this->theta)-sin(this->theta)*sin(this->theta);
            ///TransPrim.beInverseOf( Trans );
            TransPrim.zero();
            TransPrim.at(1,1)=cos(-this->theta)*cos(-this->theta);
            TransPrim.at(1,2)=sin(-this->theta)*sin(-this->theta);
            TransPrim.at(1,6)=2*sin(-this->theta)*cos(-this->theta);
            TransPrim.at(2,1)=sin(-this->theta)*sin(-this->theta);
            TransPrim.at(2,2)=cos(-this->theta)*cos(-this->theta);
            TransPrim.at(2,6)=-2*sin(-this->theta)*cos(-this->theta);
            TransPrim.at(3,3)=1;
            TransPrim.at(4,4)=cos(-this->theta);
            TransPrim.at(4,5)=-sin(-this->theta);
            TransPrim.at(5,4)=sin(-this->theta);
            TransPrim.at(5,5)=cos(-this->theta);
            TransPrim.at(6,1)=-sin(-this->theta)*cos(-this->theta);
            TransPrim.at(6,2)=sin(-this->theta)*cos(-this->theta);
            TransPrim.at(6,6)=cos(-this->theta)*cos(-this->theta)-sin(-this->theta)*sin(-this->theta);

            LocStrain.zero();
            LocStrain.beProductOf(Trans, GStrain);

            AnisotropicDamageMaterial2Status *status = static_cast< AnisotropicDamageMaterial2Status * >( this->giveStatus(gp) );

            FloatArray TempStress;
            TempStress=status->giveTempStressVector();
            FloatArray LocPreStress(6);
            LocPreStress.zero();
            LocPreStress.beProductOf(Trans, TempStress);

            double  TrEps, NeweqStrain1c, NeweqStrain2c, NeweqStrain3c;
            //TrSig = TempStress.at(1) + TempStress.at(2) + TempStress.at(3);
            TrEps = LocStrain.at(1) + LocStrain.at(2) + LocStrain.at(3);

            int flag;
            if ( TrEps >= 0){
                flag = 1;
            }else{
                flag = 0;
            }

            FloatArray EquivStrain(6);
            EquivStrain.zero();
            //Compute equivalent strain
            if (TrEps > 0) {
                EquivStrain.at(1) =  sqrt(pow(LocStrain.at(1), 2) + (pow(LocStrain.at(5),2)+pow(LocStrain.at(6),2))*pow(this->eqeps_1t/this->eqeps_1s,2));
                EquivStrain.at(3) =  sqrt(pow(LocStrain.at(2), 2) + pow(LocStrain.at(4),2)*pow(this->eqeps_2t/this->eqeps_2s,2)+pow(LocStrain.at(6),2)*pow(this->eqeps_2t/this->eqeps_1s,2));
                EquivStrain.at(5) =  sqrt(pow(LocStrain.at(3), 2) + pow(LocStrain.at(4),2)*pow(this->eqeps_2t/this->eqeps_2s,2)+pow(LocStrain.at(5),2)*pow(this->eqeps_2t/this->eqeps_1s,2));
            }else {
                EquivStrain.at(2) =  sqrt(pow(LocStrain.at(1), 2) + (pow(LocStrain.at(5),2)+pow(LocStrain.at(6),2))*pow(this->eqeps_1c/this->eqeps_1s,2));
                EquivStrain.at(4) =  sqrt(pow(LocStrain.at(2), 2) + pow(LocStrain.at(4),2)*pow(this->eqeps_2c/this->eqeps_2s,2)+pow(LocStrain.at(6),2)*pow(this->eqeps_2c/this->eqeps_1s,2));
                EquivStrain.at(6) =  sqrt(pow(LocStrain.at(3), 2) + pow(LocStrain.at(4),2)*pow(this->eqeps_2c/this->eqeps_2s,2)+pow(LocStrain.at(5),2)*pow(this->eqeps_2c/this->eqeps_1s,2));
            };

            FloatArray DamageT(3),DamageC(3);
            DamageT.zero();
            DamageC.zero();
            DamageT.at(1) = status->giveTempDamage().at(1,1);
            DamageT.at(2) = status->giveTempDamage().at(1,2);
            DamageT.at(3) = status->giveTempDamage().at(1,3);
            DamageC.at(1) = status->giveTempDamage().at(2,1);
            DamageC.at(2) = status->giveTempDamage().at(2,2);
            DamageC.at(3) = status->giveTempDamage().at(2,3);

            if (status->giveFlag_1st() == 0){
                status->setTempEquivalenstrain_1t(this->eqeps_1t);
                status->setTempEquivalenstrain_1c(this->eqeps_1c);
                status->setTempEquivalenstrain_2t(this->eqeps_2t);
                status->setTempEquivalenstrain_2c(this->eqeps_2c);
                status->setTempEquivalenstrain_3t(this->eqeps_2t);
                status->setTempEquivalenstrain_3c(this->eqeps_2c);
                status->setFlag_1st(1);
            }

            /*
            if (EquivStrain.at(1) < this->eqeps_1t){
                status->setTempEquivalenstrain_1t(this->eqeps_1t);
            };
            if (EquivStrain.at(2) < this->eqeps_1c){
                status->setTempEquivalenstrain_1c(this->eqeps_1c);
            };
            if (EquivStrain.at(3) < this->eqeps_2t){
                status->setTempEquivalenstrain_2t(this->eqeps_2t);
            };
            if (EquivStrain.at(4) < this->eqeps_2c){
                status->setTempEquivalenstrain_2c(this->eqeps_2c);
            };
            if (EquivStrain.at(5) < this->eqeps_2t){
                status->setTempEquivalenstrain_3t(this->eqeps_2t);
            };
            if (EquivStrain.at(6) < this->eqeps_2c){
                status->setTempEquivalenstrain_3c(this->eqeps_2c);
            };
            */
            ///
            if (EquivStrain.at(1) > status->giveTempEquivalenstrain_1t()){
                DamageT.at(1) = 1-exp(-(EquivStrain.at(1) - this->eqeps_1t)/this->alpha_1t);
                status->setTempEquivalenstrain_1t(EquivStrain.at(1));
                flag=1;
            }else{
                 DamageT.at(1) = DamageT.at(1);
            };
            ///
            if (EquivStrain.at(3) > status->giveTempEquivalenstrain_2t()){
                DamageT.at(2) = 1-exp(-(EquivStrain.at(3) - this->eqeps_2t)/this->alpha_2t);
                status->setTempEquivalenstrain_2t(EquivStrain.at(3));
                flag=1;
            }else{
                 DamageT.at(2) = DamageT.at(2);
            };
            ///
            if (EquivStrain.at(5) > status->giveTempEquivalenstrain_3t()){
                DamageT.at(3) = 1-exp(-(EquivStrain.at(5) - this->eqeps_2t)/this->alpha_2t);
                status->setTempEquivalenstrain_3t(EquivStrain.at(5));
                flag=1;
            }else{
                 DamageT.at(3) = DamageT.at(3);
            };
            ///
            double confinement = (LocPreStress.at(3)+LocPreStress.at(2))*0.5;
            if (confinement < 0){
                NeweqStrain1c = EquivStrain.at(2)+this->eta*confinement;
            } else {
                NeweqStrain1c = EquivStrain.at(2);
            }
            if ( NeweqStrain1c > status->giveTempEquivalenstrain_1c()){
                DamageC.at(1) = exp((NeweqStrain1c -this->beta_1c)/this->alpha_1c)/(1+exp((NeweqStrain1c-this->beta_1c)/this->alpha_1c));
                status->setTempEquivalenstrain_1c(NeweqStrain1c);
            }else{
                DamageC.at(1) = DamageC.at(1);
            };
            ///
            confinement = LocPreStress.at(1);
            if (confinement < 0){
                NeweqStrain2c = EquivStrain.at(4)+this->eta*confinement;
            } else {
                NeweqStrain2c = EquivStrain.at(4);
            }
            if (NeweqStrain2c > status->giveTempEquivalenstrain_2c()){
                DamageC.at(2) = exp((NeweqStrain2c-this->beta_2c)/this->alpha_2c)/(1+exp((NeweqStrain2c-this->beta_2c)/this->alpha_2c));
                status->setTempEquivalenstrain_2c(NeweqStrain2c);
            }else{
                DamageC.at(2) = DamageC.at(2);
            };
            ///
            confinement = LocPreStress.at(1);
            if (confinement < 0){
                NeweqStrain3c = EquivStrain.at(6)+this->eta*confinement;
            } else {
                NeweqStrain3c = EquivStrain.at(6);
            }
            if (NeweqStrain3c > status->giveTempEquivalenstrain_3c()){
                DamageC.at(3) = exp((NeweqStrain3c-this->beta_2c)/this->alpha_2c)/(1+exp((NeweqStrain3c-this->beta_2c)/this->alpha_2c));
                status->setTempEquivalenstrain_2c(NeweqStrain3c);
            }else{
                DamageC.at(3) = DamageC.at(3);
            };

            status->setTempFlag(flag);
            //status->setTempFlag2(flag2);
            //status->setTempFlag3(flag3);

            FloatMatrix MATS(6,6), MATC(6,6);
            MATS.zero();
            MATC.zero();

            double d1, d2, d3, nu21, G23;

             if (flag == 1){
                d1 = DamageT.at(1);
                d2 = DamageT.at(2);
                d3 = DamageT.at(3);
             }else{
                d1 = DamageC.at(1);
                d2 = DamageC.at(2);
                d3 = DamageC.at(3);
             };



            nu21 = this->E22*this->nu12/this->E11;
            G23  = this->E22/(2*(1+this->nu23));

            MATS.at(1,1)= 1/((1-d1)*this->E11);
            MATS.at(1,2)=-this->nu12/this->E11;
            MATS.at(1,3)=-this->nu12/this->E11;
            MATS.at(2,1)=-nu21/this->E22;
            MATS.at(2,2)= 1/((1-d2)*this->E22);
            MATS.at(2,3)=-this->nu23/this->E22;
            MATS.at(3,1)=-nu21/this->E22;
            MATS.at(3,2)=-this->nu23/this->E22;
            MATS.at(3,3)=1/((1-d3)*this->E22);
            MATS.at(4,4)=1/((1-d2)*(1-d3)*G23);
            MATS.at(5,5)=1/((1-d1)*(1-d3)*this->G12);
            MATS.at(6,6)=1/((1-d1)*(1-d2)*this->G12);

            MATC.beInverseOf( MATS);

           FloatArray LocStress(6),Stress(6);

           LocStress.zero();
           LocStrain.at(4) = 2*LocStrain.at(4);
           LocStrain.at(5) = 2*LocStrain.at(5);
           LocStrain.at(6) = 2*LocStrain.at(6);

           LocStress.beProductOf(MATC, LocStrain);

           Stress.zero();
           Stress.beProductOf(TransPrim, LocStress);

           /*
           TempDamageT.zero();
           TempDamageC.zero();
           TempDamageT.beProductOf(TransPrim, DamageT);
           TempDamageC.beProductOf(TransPrim, DamageC);

           FloatMatrix Damage(3,3);
           Damage.zero();
           Damage.at(1,1)=TempDamageT.at(1,1);
           Damage.at(1,2)=TempDamageT.at(2,1);
           Damage.at(1,3)=TempDamageT.at(4,1);
           Damage.at(2,1)=TempDamageC.at(1,1);
           Damage.at(2,2)=TempDamageC.at(2,1);
           Damage.at(2,3)=TempDamageC.at(4,1);
           */
           FloatMatrix Damage(3,3);
           Damage.zero();
           Damage.at(1,1)=DamageT.at(1);
           Damage.at(1,2)=DamageT.at(2);
           Damage.at(1,3)=DamageT.at(3);
           Damage.at(2,1)=DamageC.at(1);
           Damage.at(2,2)=DamageC.at(2);
           Damage.at(2,3)=DamageC.at(3);
           status->setTempDamage( Damage );
           status->letTempStrainVectorBe(totalStrain);
           answer.resize(6);
           answer.at(1)=Stress.at(1);
           answer.at(2)=Stress.at(2);
           answer.at(3)=Stress.at(3);
           answer.at(4)=Stress.at(4);
           answer.at(5)=Stress.at(5);
           answer.at(6)=Stress.at(6);

           status->letTempStressVectorBe(answer);
        #ifdef keep_track_of_dissipated_energy
           status->computeWork(gp);
        #endif
           return;
}

void
AnisotropicDamageMaterial2 :: giveRealStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp,
                                                              const FloatArray &totalStrain, TimeStep *atTime)
//
// uses the following special methods:
//   computePrincValDir2D  ...  evaluation of eigenvalues and eigenvector of a symmetric 2x2 matrix
//
{
#define TOLERANCE 1.e-6 // convergence tolerance for the internal iteration
    //this->initTempStatus(gp); //updata the status, so the current iteration starts from the previous converged results;
    // subtract the stress-independent part of strains (e.g. due to temperature)
    FloatArray SDstrainVector;

    this->giveStressDependentPartOfStrainVector(SDstrainVector, gp, totalStrain, atTime, VM_Total);

    // make sure strain components are correctly stored as the following--WJ

    FloatMatrix Trans(4,4),TransPrim(4,4);
    FloatArray GStrain(4),LocStrain(4);
    GStrain.at(1)=SDstrainVector.at(1);
    GStrain.at(2)=SDstrainVector.at(2);
    GStrain.at(3)=SDstrainVector.at(3);
    GStrain.at(4)=SDstrainVector.at(4)/2;

    TransPrim.zero();
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
    ///TransPrim.beInverseOf( Trans );
    TransPrim.at(1,1)=cos(-this->theta)*cos(-this->theta);
    TransPrim.at(1,2)=sin(-this->theta)*sin(-this->theta);
    TransPrim.at(1,3)=0;
    TransPrim.at(1,4)=2*sin(-this->theta)*cos(-this->theta);
    TransPrim.at(2,1)=sin(-this->theta)*sin(-this->theta);
    TransPrim.at(2,2)=cos(-this->theta)*cos(-this->theta);
    TransPrim.at(2,3)=0;
    TransPrim.at(2,4)=-2*sin(-this->theta)*cos(-this->theta);
    TransPrim.at(3,1)=0;
    TransPrim.at(3,2)=0;
    TransPrim.at(3,3)=1;
    TransPrim.at(3,4)=0;
    TransPrim.at(4,1)=-sin(-this->theta)*cos(-this->theta);
    TransPrim.at(4,2)=sin(-this->theta)*cos(-this->theta);
    TransPrim.at(4,3)=0;
    TransPrim.at(4,4)=cos(-this->theta)*cos(-this->theta)-sin(-this->theta)*sin(-this->theta);

    LocStrain.zero();
    LocStrain.beProductOf(Trans, GStrain);
    AnisotropicDamageMaterial2Status *status = static_cast< AnisotropicDamageMaterial2Status * >( this->giveStatus(gp) );
    FloatArray TempStress;
    TempStress=status->giveTempStressVector();
    double TrEps, NeweqStrain1c, NeweqStrain2c;
    //TrSig = TempStress.at(1) + TempStress.at(2) + TempStress.at(3);
    TrEps = LocStrain.at(1) + LocStrain.at(2) ;
    FloatArray LocPreStress(4);
    LocPreStress.zero();
    LocPreStress.beProductOf(Trans, TempStress);
    int flag;
    if ( TrEps >= 0){
        flag = 1;
    }else{
        flag = 0;
    }

    FloatArray EquivStrain(4);
    EquivStrain.zero();
    //Compute equivalent strain
    this->computeEquivalentStrain(EquivStrain, LocStrain, TrEps, gp, atTime);

    //Should check whether use tempdamage or damage value
    /*
    FloatMatrix TempDamageT(4,1),TempDamageC(4,1);
    TempDamageT.zero();
    TempDamageC.zero();
    TempDamageT.at(1,1) = status->giveTempDamage().at(1,1);
    TempDamageT.at(2,1) = status->giveTempDamage().at(1,2);
    TempDamageT.at(4,1) = status->giveTempDamage().at(1,3);
    TempDamageC.at(1,1) = status->giveTempDamage().at(2,1);
    TempDamageC.at(2,1) = status->giveTempDamage().at(2,2);
    TempDamageC.at(4,1) = status->giveTempDamage().at(2,3);
    FloatMatrix DamageT(4,1),DamageC(4,1);
    DamageT.zero();
    DamageT.beProductOf(Trans, TempDamageT);
    DamageC.zero();
    DamageC.beProductOf(Trans, TempDamageC);
    */

    FloatArray DamageT(2),DamageC(2);
    DamageT.zero();
    DamageC.zero();
    DamageT.at(1) = status->giveTempDamage().at(1,1);
    DamageT.at(2) = status->giveTempDamage().at(1,2);
    DamageC.at(1) = status->giveTempDamage().at(2,1);
    DamageC.at(2) = status->giveTempDamage().at(2,2);

    if (status->giveFlag_1st() == 0){
        status->setTempEquivalenstrain_1t(this->eqeps_1t);
        status->setTempEquivalenstrain_1c(this->eqeps_1c);
        status->setTempEquivalenstrain_2t(this->eqeps_2t);
        status->setTempEquivalenstrain_2c(this->eqeps_2c);
        status->setFlag_1st(1);
    }
    /*
    if (EquivStrain.at(1) < this->eqeps_1t){
        status->setTempEquivalenstrain_1t(this->eqeps_1t);
    };
    if (EquivStrain.at(2) < this->eqeps_1c){
        status->setTempEquivalenstrain_1c(this->eqeps_1c);
    };
    if (EquivStrain.at(3) < this->eqeps_2t){
        status->setTempEquivalenstrain_2t(this->eqeps_2t);
    };
    if (EquivStrain.at(4) < this->eqeps_2c){
        status->setTempEquivalenstrain_2c(this->eqeps_2c);
    };
    */
    if (EquivStrain.at(1) > status->giveTempEquivalenstrain_1t()){
        DamageT.at(1) = 1-exp(-(EquivStrain.at(1) - this->eqeps_1t)/this->alpha_1t);
        /*
        if (TrStress > 0){
            DamageT.at(1) = 1-exp(-(EquivStrain.at(1) - this->eqeps_1t)/this->alpha_1t);
        }else{
            DamageT.at(1) = (EquivStrain.at(1) - this->eqeps_1t)/(this->alpha_1t+(EquivStrain.at(1) - this->eqeps_1t));
        };*/
        status->setTempEquivalenstrain_1t(EquivStrain.at(1));
        flag=1;
    }else{
         DamageT.at(1) = DamageT.at(1);
    };
    double confinement = (LocPreStress.at(3)+LocPreStress.at(2))*0.5;
    if (confinement < 0){
        NeweqStrain1c = EquivStrain.at(2)+this->eta*confinement;
    } else {
        NeweqStrain1c = EquivStrain.at(2);
    }
    if ( NeweqStrain1c >= status->giveTempEquivalenstrain_1c()){
        DamageC.at(1) = exp((NeweqStrain1c -beta_1c)/this->alpha_1c)/(1+exp((NeweqStrain1c-beta_1c)/this->alpha_1c));
        status->setTempEquivalenstrain_1c(NeweqStrain1c);
        flag = 0;
    }else{
        DamageC.at(1) = DamageC.at(1);
   };

    if (EquivStrain.at(3) > status->giveTempEquivalenstrain_2t()){
        DamageT.at(2) = 1-exp(-(EquivStrain.at(3) - this->eqeps_2t)/this->alpha_2t);
        /*
        if (TrStress > 0){
            DamageT.at(2) = 1-exp(-(EquivStrain.at(3) - this->eqeps_2t)/this->alpha_2t);
        }else{
            DamageT.at(2) = (EquivStrain.at(3) - this->eqeps_2t)/(this->alpha_2t+(EquivStrain.at(3) - this->eqeps_2t));
        };*/
        status->setTempEquivalenstrain_2t(EquivStrain.at(3));
        flag=1;
    }else{
         DamageT.at(2) = DamageT.at(2);
    };

    confinement = LocPreStress.at(1);
    if (confinement < 0){
        NeweqStrain2c = EquivStrain.at(4)+this->eta*confinement;
    } else {
        NeweqStrain2c = EquivStrain.at(4);
    }
    if (NeweqStrain2c >= status->giveTempEquivalenstrain_2c()){
        DamageC.at(2) = exp((NeweqStrain2c-beta_2c)/this->alpha_2c)/(1+exp((NeweqStrain2c-this->beta_2c)/this->alpha_2c));
        status->setTempEquivalenstrain_2c(NeweqStrain2c);
        flag = 0;
    }else{
        DamageC.at(2) = DamageC.at(2);
   };

    status->setTempFlag(flag);

    FloatMatrix MATC(4,4);
    MATC.zero();

    double d1, d2, D, nu21;

     if (flag == 1){
        d1=DamageT.at(1);
        d2=DamageT.at(2);
    }else{
        d1=DamageC.at(1);
        d2= DamageC.at(2);
    };


    nu21=this->E22*this->nu12/this->E11;
    D= (1-d2)*pow(this->nu23,2)+2*(1-d1)*(1-d2)*this->nu12*nu21*this->nu23+(1-d1)*(2-d2)*this->nu12*nu21-1;
    MATC.at(1,1)=this->E11*(1-d1)*((1-d2)*pow(this->nu23,2)-1)/D;
    MATC.at(1,2)=-this->E11*nu21*(1-d1)*(1-d2)*(1+this->nu23)/D;
    MATC.at(1,3)=-this->E11*nu21*(1-d1)*(1+(1-d2)*this->nu23)/D;
    MATC.at(1,4)=0;
    MATC.at(2,1)=MATC.at(1,2);
    MATC.at(2,2)=this->E22*(1-d2)*((1-d1)*this->nu12*nu21-1)/D;
    MATC.at(2,3)=-this->E22*(1-d2)*(this->nu23+(1-d1)*this->nu12*nu21)/D;
    MATC.at(2,4)=0;
    MATC.at(3,1)=MATC.at(1,3);
    MATC.at(3,2)=MATC.at(2,3);
    MATC.at(3,3)=this->E22*(1-d2)*(1-d1)*(this->nu12*nu21-1)/D;
    MATC.at(3,4)=0;
    MATC.at(4,1)=0;
    MATC.at(4,2)=0;
    MATC.at(4,3)=0;
    MATC.at(4,4)=this->G12*(1-d1)*(1-d2);

   FloatArray LocStress(4),Stress(4);
   LocStress.zero();
   LocStrain.at(4) = 2*LocStrain.at(4);
   LocStress.beProductOf(MATC, LocStrain);

   Stress.zero();
   Stress.beProductOf(TransPrim, LocStress);

   /*
   TempDamageT.zero();
   TempDamageC.zero();
   TempDamageT.beProductOf(TransPrim, DamageT);
   TempDamageC.beProductOf(TransPrim, DamageC);

   FloatMatrix Damage(3,3);
   Damage.zero();
   Damage.at(1,1)=TempDamageT.at(1,1);
   Damage.at(1,2)=TempDamageT.at(2,1);
   Damage.at(1,3)=TempDamageT.at(4,1);
   Damage.at(2,1)=TempDamageC.at(1,1);
   Damage.at(2,2)=TempDamageC.at(2,1);
   Damage.at(2,3)=TempDamageC.at(4,1);
   */
   FloatMatrix Damage(3,3);
   Damage.zero();
   Damage.at(1,1)=DamageT.at(1);
   Damage.at(1,2)=DamageT.at(2);
   Damage.at(2,1)=DamageC.at(1);
   Damage.at(2,2)=DamageC.at(2);

   status->setdamaget1(DamageT.at(1));
   status->setdamaget2(DamageT.at(2));
   status->setdamagec1(DamageC.at(1));
   status->setdamagec2(DamageC.at(2));

   status->setTempDamage( Damage );
   status->letTempStrainVectorBe(totalStrain);
   answer.resize(4);
   answer.at(1)=Stress.at(1);
   answer.at(2)=Stress.at(2);
   answer.at(3)=Stress.at(3);
   answer.at(4)=Stress.at(4);
   status->letTempStressVectorBe(answer);
#ifdef keep_track_of_dissipated_energy
   status->computeWork(gp);
#endif
   return;
}

void
AnisotropicDamageMaterial2 :: computeEquivalentStrain(FloatArray &EquivStrain, const FloatArray &strain,const double &TrEps, GaussPoint *gp, TimeStep *tStep)
{

    if (TrEps > 0) {
        EquivStrain.at(1) =  pow(strain.at(1), 2) + pow(strain.at(4),2)*pow(this->eqeps_1t/this->eqeps_1s,2);
        EquivStrain.at(1) = sqrt(EquivStrain.at(1));
        EquivStrain.at(3) =  pow(strain.at(2), 2) + pow(strain.at(4),2)*pow(this->eqeps_2t/this->eqeps_2s,2);
        EquivStrain.at(3) = sqrt(EquivStrain.at(3));
    }else {
        EquivStrain.at(2) =  pow(strain.at(1), 2) + pow(strain.at(4),2)*pow(this->eqeps_1c/this->eqeps_1s,2);
        EquivStrain.at(2) = sqrt(EquivStrain.at(2));
        EquivStrain.at(4) =  pow(strain.at(2), 2) + pow(strain.at(4),2)*pow(this->eqeps_2c/this->eqeps_2s,2);
        EquivStrain.at(4) = sqrt(EquivStrain.at(4));
    };
    return;
}

void AnisotropicDamageMaterial2 :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,MatResponseMode mode,
                                           GaussPoint *gp,TimeStep *tStep)
{

    this->initTempStatus(gp);
    if ( mode == ElasticStiffness ) {

        FloatMatrix MATS(6,6);

        MATS.zero();

        double nu21, G23;
        nu21 = this->E22*this->nu12/this->E11;
        G23  = this->E22/(2*(1+this->nu23));

        MATS.at(1,1)= 1/(this->E11);
        MATS.at(1,2)=-this->nu12/this->E11;
        MATS.at(1,3)=-this->nu12/this->E11;
        MATS.at(2,1)=-nu21/this->E22;
        MATS.at(2,2)= 1/(this->E22);
        MATS.at(2,3)=-this->nu23/this->E22;
        MATS.at(3,1)=-nu21/this->E22;
        MATS.at(3,2)=-this->nu23/this->E22;
        MATS.at(3,3)=1/(this->E22);
        MATS.at(4,4)=1/(G23);
        MATS.at(5,5)=1/(this->G12);
        MATS.at(6,6)=1/(this->G12);

        answer.resize(6,6);
        answer.beInverseOf( MATS );
        return;

    } else {

        FloatMatrix Trans(6,6),TransPrim(6,6);
        Trans.zero();
        Trans.at(1,1)=cos(this->theta)*cos(this->theta);
        Trans.at(1,2)=sin(this->theta)*sin(this->theta);
        Trans.at(1,6)=2*sin(this->theta)*cos(this->theta);
        Trans.at(2,1)=sin(this->theta)*sin(this->theta);
        Trans.at(2,2)=cos(this->theta)*cos(this->theta);
        Trans.at(2,6)=-2*sin(this->theta)*cos(this->theta);
        Trans.at(3,3)=1;
        Trans.at(4,4)=cos(this->theta);
        Trans.at(4,5)=-sin(this->theta);
        Trans.at(5,4)=sin(this->theta);
        Trans.at(5,5)=cos(this->theta);
        Trans.at(6,1)=-sin(this->theta)*cos(this->theta);
        Trans.at(6,2)=sin(this->theta)*cos(this->theta);
        Trans.at(6,6)=cos(this->theta)*cos(this->theta)-sin(this->theta)*sin(this->theta);
        ///TransPrim.beInverseOf( Trans );
        TransPrim.zero();
        TransPrim.at(1,1)=cos(-this->theta)*cos(-this->theta);
        TransPrim.at(1,2)=sin(-this->theta)*sin(-this->theta);
        TransPrim.at(1,6)=2*sin(-this->theta)*cos(-this->theta);
        TransPrim.at(2,1)=sin(-this->theta)*sin(-this->theta);
        TransPrim.at(2,2)=cos(-this->theta)*cos(-this->theta);
        TransPrim.at(2,6)=-2*sin(-this->theta)*cos(-this->theta);
        TransPrim.at(3,3)=1;
        TransPrim.at(4,4)=cos(-this->theta);
        TransPrim.at(4,5)=-sin(-this->theta);
        TransPrim.at(5,4)=sin(-this->theta);
        TransPrim.at(5,5)=cos(-this->theta);
        TransPrim.at(6,1)=-sin(-this->theta)*cos(-this->theta);
        TransPrim.at(6,2)=sin(-this->theta)*cos(-this->theta);
        TransPrim.at(6,6)=cos(-this->theta)*cos(-this->theta)-sin(-this->theta)*sin(-this->theta);
        //
        AnisotropicDamageMaterial2Status *status = static_cast< AnisotropicDamageMaterial2Status * >( this->giveStatus(gp) );
        /*
        FloatMatrix TempDamageT(4,1),TempDamageC(4,1);
        TempDamageT.zero();
        TempDamageC.zero();
        TempDamageT.at(1,1) = status->giveTempDamage().at(1,1);
        TempDamageT.at(2,1) = status->giveTempDamage().at(1,2);
        TempDamageT.at(4,1) = status->giveTempDamage().at(1,3);
        TempDamageC.at(1,1) = status->giveTempDamage().at(2,1);
        TempDamageC.at(2,1) = status->giveTempDamage().at(2,2);
        TempDamageC.at(4,1) = status->giveTempDamage().at(2,3);
        FloatMatrix DamageT(4,1),DamageC(4,1);
        DamageT.zero();
        DamageT.beProductOf(Trans, TempDamageT);
        DamageC.zero();
        DamageC.beProductOf(Trans, TempDamageC);
        */

        FloatArray DamageT(3),DamageC(3);
        DamageT.zero();
        DamageC.zero();
        DamageT.at(1) = status->giveTempDamage().at(1,1);
        DamageT.at(2) = status->giveTempDamage().at(1,2);
        DamageT.at(3) = status->giveTempDamage().at(1,3);
        DamageC.at(1) = status->giveTempDamage().at(2,1);
        DamageC.at(2) = status->giveTempDamage().at(2,2);
        DamageC.at(3) = status->giveTempDamage().at(2,3);


        int flag;
        flag=status->giveTempFlag();


        FloatMatrix MATS(6,6), MATC(6,6);
        MATS.zero();
        MATC.zero();

        double d1, d2, d3, nu21, G23;

         if (flag == 1){
            d1 = DamageT.at(1);
            d2 = DamageT.at(2);
            d3 = DamageT.at(3);
         }else{
            d1 = DamageC.at(1);
            d2 = DamageC.at(2);
            d3 = DamageC.at(3);
         };

        nu21 = this->E22*this->nu12/this->E11;
        G23  = this->E22/(2*(1+this->nu23));

        MATS.at(1,1)= 1/((1-d1)*this->E11);
        MATS.at(1,2)=-this->nu12/this->E11;
        MATS.at(1,3)=-this->nu12/this->E11;
        MATS.at(2,1)=-nu21/this->E22;
        MATS.at(2,2)= 1/((1-d2)*this->E22);
        MATS.at(2,3)=-this->nu23/this->E22;
        MATS.at(3,1)=-nu21/this->E22;
        MATS.at(3,2)=-this->nu23/this->E22;
        MATS.at(3,3)=1/((1-d3)*this->E22);
        MATS.at(4,4)=1/((1-d2)*(1-d3)*G23);
        MATS.at(5,5)=1/((1-d1)*(1-d3)*this->G12);
        MATS.at(6,6)=1/((1-d1)*(1-d2)*this->G12);

        MATC.beInverseOf( MATS);

        FloatMatrix R(6,6),RPrim(6,6),Temp1(6,6),Temp2(6,6),Temp3(6,6);
        R.zero();
        R.at(1,1)=1;R.at(2,2)=1;R.at(3,3)=1;
        R.at(4,4)=2;R.at(5,5)=2;R.at(6,6)=2;
        RPrim.zero();
        RPrim.at(1,1)=1;RPrim.at(2,2)=1;RPrim.at(3,3)=1;
        RPrim.at(4,4)=0.5;RPrim.at(5,5)=0.5;RPrim.at(6,6)=0.5;

        Temp1.zero();
        Temp2.zero();
        Temp1.beProductOf(Trans,RPrim);
        Temp2.beProductOf(R,Temp1);
        Temp3.beProductOf(TransPrim,MATC);
        answer.resize(6,6);
        answer.beProductOf(Temp3,Temp2);
        return;
    };
}
void AnisotropicDamageMaterial2 :: givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseMode mode,
                                                           GaussPoint *gp, TimeStep *atTime)
{
    //this->initTempStatus(gp);
    if ( mode == ElasticStiffness ) {

        double nu21;
        FloatMatrix MATS(4,4);
        nu21=this->E22*this->nu12/this->E11;
        MATS.at(1,1)=1/this->E11;
        MATS.at(1,2)=-nu21/this->E22;
        MATS.at(1,3)=-nu21/this->E22;
        MATS.at(1,4)=0;
        MATS.at(2,1)=-this->nu12/this->E11;
        MATS.at(2,2)=1/this->E22;
        MATS.at(2,3)=-this->nu23/this->E22;
        MATS.at(2,4)=0;
        MATS.at(3,1)=-this->nu12/this->E11;
        MATS.at(3,2)=-this->nu23/this->E22;
        MATS.at(3,3)=1/this->E22;
        MATS.at(3,4)=0;
        MATS.at(4,1)=0;
        MATS.at(4,2)=0;
        MATS.at(4,3)=0;
        MATS.at(4,4)=1/this->G12;
        answer.resize(4,4);
        answer.beInverseOf( MATS );
        return;

    } else if (mode == SecantStiffness) {

        FloatMatrix Trans(4,4),TransPrim(4,4);
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
            //TransPrim.beInverseOf( Trans );
        TransPrim.zero();
        TransPrim.at(1,1)=cos(-this->theta)*cos(-this->theta);
        TransPrim.at(1,2)=sin(-this->theta)*sin(-this->theta);
        TransPrim.at(1,3)=0;
        TransPrim.at(1,4)=2*sin(-this->theta)*cos(-this->theta);
        TransPrim.at(2,1)=sin(-this->theta)*sin(-this->theta);
        TransPrim.at(2,2)=cos(-this->theta)*cos(-this->theta);
        TransPrim.at(2,3)=0;
        TransPrim.at(2,4)=-2*sin(-this->theta)*cos(-this->theta);
        TransPrim.at(3,1)=0;
        TransPrim.at(3,2)=0;
        TransPrim.at(3,3)=1;
        TransPrim.at(3,4)=0;
        TransPrim.at(4,1)=-sin(-this->theta)*cos(-this->theta);
        TransPrim.at(4,2)=sin(-this->theta)*cos(-this->theta);
        TransPrim.at(4,3)=0;
        TransPrim.at(4,4)=cos(-this->theta)*cos(-this->theta)-sin(-this->theta)*sin(-this->theta);

        AnisotropicDamageMaterial2Status *status = static_cast< AnisotropicDamageMaterial2Status * >( this->giveStatus(gp) );
        /*
        FloatMatrix TempDamageT(4,1),TempDamageC(4,1);
        TempDamageT.zero();
        TempDamageC.zero();
        TempDamageT.at(1,1) = status->giveTempDamage().at(1,1);
        TempDamageT.at(2,1) = status->giveTempDamage().at(1,2);
        TempDamageT.at(4,1) = status->giveTempDamage().at(1,3);
        TempDamageC.at(1,1) = status->giveTempDamage().at(2,1);
        TempDamageC.at(2,1) = status->giveTempDamage().at(2,2);
        TempDamageC.at(4,1) = status->giveTempDamage().at(2,3);
        FloatMatrix DamageT(4,1),DamageC(4,1);
        DamageT.zero();
        DamageT.beProductOf(Trans, TempDamageT);
        DamageC.zero();
        DamageC.beProductOf(Trans, TempDamageC);
        */
        FloatArray DamageT(2),DamageC(2);
        DamageT.zero();
        DamageC.zero();
        DamageT.at(1) = status->giveTempDamage().at(1,1);
        DamageT.at(2) = status->giveTempDamage().at(1,2);
        DamageC.at(1) = status->giveTempDamage().at(2,1);
        DamageC.at(2) = status->giveTempDamage().at(2,2);

        int flag;
        flag=status->giveTempFlag();

        double d1, d2, nu21, D;
        if (flag ==1 ){
            d1=DamageT.at(1);
            d2=DamageT.at(2);
        }else{
            d1=DamageC.at(1);
            d2=DamageC.at(2);
        };

        FloatMatrix MATC(4,4);
        nu21=this->E22*this->nu12/this->E11;
        D= (1-d2)*pow(this->nu23,2)+2*(1-d1)*(1-d2)*this->nu12*nu21*this->nu23+(1-d1)*(2-d2)*this->nu12*nu21-1;
        MATC.at(1,1)=this->E11*(1-d1)*((1-d2)*pow(this->nu23,2)-1)/D;
        MATC.at(1,2)=-this->E11*nu21*(1-d1)*(1-d2)*(1+this->nu23)/D;
        MATC.at(1,3)=-this->E11*nu21*(1-d1)*(1+(1-d2)*this->nu23)/D;
        MATC.at(1,4)=0;
        MATC.at(2,1)=MATC.at(1,2);
        MATC.at(2,2)=this->E22*(1-d2)*((1-d1)*this->nu12*nu21-1)/D;
        MATC.at(2,3)=-this->E22*(1-d2)*(this->nu23+(1-d1)*this->nu12*nu21)/D;
        MATC.at(2,4)=0;
        MATC.at(3,1)=MATC.at(1,3);
        MATC.at(3,2)=MATC.at(2,3);
        MATC.at(3,3)=this->E22*(1-d2)*(1-d1)*(this->nu12*nu21-1)/D;
        MATC.at(3,4)=0;
        MATC.at(4,1)=0;
        MATC.at(4,2)=0;
        MATC.at(4,3)=0;
        MATC.at(4,4)=this->G12*(1-d1)*(1-d2);

        FloatMatrix R(4,4),RPrim(4,4),Temp1(4,4),Temp2(4,4),Temp3(4,4);
        R.zero();
        R.at(1,1)=1;
        R.at(2,2)=1;
        R.at(3,3)=1;
        R.at(4,4)=2;
        RPrim.zero();
        RPrim.at(1,1)=1;
        RPrim.at(2,2)=1;
        RPrim.at(3,3)=1;
        RPrim.at(4,4)=0.5;

        Temp1.zero();
        Temp2.zero();
        Temp1.beProductOf(Trans,RPrim);
        Temp2.beProductOf(R,Temp1);
        Temp3.beProductOf(TransPrim,MATC);
        answer.resize(4,4);
        answer.beProductOf(Temp3,Temp2);
        return;

    }else {

        FloatMatrix Trans(4,4),TransPrim(4,4);
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
            //TransPrim.beInverseOf( Trans );
        TransPrim.zero();
        TransPrim.at(1,1)=cos(-this->theta)*cos(-this->theta);
        TransPrim.at(1,2)=sin(-this->theta)*sin(-this->theta);
        TransPrim.at(1,3)=0;
        TransPrim.at(1,4)=2*sin(-this->theta)*cos(-this->theta);
        TransPrim.at(2,1)=sin(-this->theta)*sin(-this->theta);
        TransPrim.at(2,2)=cos(-this->theta)*cos(-this->theta);
        TransPrim.at(2,3)=0;
        TransPrim.at(2,4)=-2*sin(-this->theta)*cos(-this->theta);
        TransPrim.at(3,1)=0;
        TransPrim.at(3,2)=0;
        TransPrim.at(3,3)=1;
        TransPrim.at(3,4)=0;
        TransPrim.at(4,1)=-sin(-this->theta)*cos(-this->theta);
        TransPrim.at(4,2)=sin(-this->theta)*cos(-this->theta);
        TransPrim.at(4,3)=0;
        TransPrim.at(4,4)=cos(-this->theta)*cos(-this->theta)-sin(-this->theta)*sin(-this->theta);

        AnisotropicDamageMaterial2Status *status = static_cast< AnisotropicDamageMaterial2Status * >( this->giveStatus(gp) );
        /*
        FloatMatrix TempDamageT(4,1),TempDamageC(4,1);
        TempDamageT.zero();
        TempDamageC.zero();
        TempDamageT.at(1,1) = status->giveTempDamage().at(1,1);
        TempDamageT.at(2,1) = status->giveTempDamage().at(1,2);
        TempDamageT.at(4,1) = status->giveTempDamage().at(1,3);
        TempDamageC.at(1,1) = status->giveTempDamage().at(2,1);
        TempDamageC.at(2,1) = status->giveTempDamage().at(2,2);
        TempDamageC.at(4,1) = status->giveTempDamage().at(2,3);
        FloatMatrix DamageT(4,1),DamageC(4,1);
        DamageT.zero();
        DamageT.beProductOf(Trans, TempDamageT);
        DamageC.zero();
        DamageC.beProductOf(Trans, TempDamageC);
        */
        FloatArray DamageT(2),DamageC(2);
        DamageT.zero();
        DamageC.zero();
        DamageT.at(1) = status->giveTempDamage().at(1,1);
        DamageT.at(2) = status->giveTempDamage().at(1,2);
        DamageC.at(1) = status->giveTempDamage().at(2,1);
        DamageC.at(2) = status->giveTempDamage().at(2,2);

        FloatArray GStrain(4), LocStrain(4);
        GStrain = status->giveTempStrainVector();
        GStrain.at(4) = GStrain.at(4)/2;
        LocStrain.zero();
        LocStrain.beProductOf(Trans, GStrain);

        int flag;
        flag=status->giveTempFlag();

        double d1, d2, nu21, D, eqeps1, eqeps2, Domega_Deqeps1, Domega_Deqeps2;
        FloatArray Domega1_Dstrain(4), Domega2_Dstrain(4);
        Domega1_Dstrain.zero();
        Domega2_Dstrain.zero();

        if (flag ==1 ){
            d1=DamageT.at(1);
            d2=DamageT.at(2);
            eqeps1 = status->giveTempEquivalenstrain_1t();
            eqeps2 = status->giveTempEquivalenstrain_2t();
            if (d1 > 0){
                Domega_Deqeps1 = exp(-(eqeps1 - this->eqeps_1t)/this->alpha_1t)/this->alpha_1t;
            }else{
                Domega_Deqeps1 = 0;
            }
            if (d2 > 0){
                Domega_Deqeps2 = exp(-(eqeps2 - this->eqeps_2t)/this->alpha_2t)/this->alpha_2t;
            }else{
                Domega_Deqeps2 = 0;
            }
            if (eqeps1 > 0){
                Domega1_Dstrain.at(1)= Domega_Deqeps1*2*LocStrain.at(1)/eqeps1;
                Domega1_Dstrain.at(4)= Domega_Deqeps1*2*LocStrain.at(4)*pow(this->eqeps_1t/this->eqeps_1s,2)/eqeps1;
            }
            if (eqeps2 > 0){
                Domega2_Dstrain.at(1)= Domega_Deqeps2*2*LocStrain.at(2)/eqeps2;
                Domega2_Dstrain.at(4)= Domega_Deqeps2*2*LocStrain.at(4)*pow(this->eqeps_2t/this->eqeps_1s,2)/eqeps2;
            }
        }else{
            d1=DamageC.at(1);
            d2=DamageC.at(2);
            eqeps1 = status->giveTempEquivalenstrain_1c();
            eqeps2 = status->giveTempEquivalenstrain_2c();
            if (d1 > 0){
                Domega_Deqeps1 = exp((eqeps1 -this->beta_1c)/this->alpha_1c)/this->alpha_1c/pow((1+exp((eqeps1-this->beta_1c)/this->alpha_1c)),2);
            }else{
                Domega_Deqeps1 = 0;
            }
            if (d2 > 0){
                Domega_Deqeps2 = exp((eqeps2 -this->beta_2c)/this->alpha_2c)/this->alpha_2c/pow((1+exp((eqeps2-this->beta_2c)/this->alpha_2c)),2);
            }else{
                Domega_Deqeps2 = 0;
            }
            if (eqeps1 > 0){
                Domega1_Dstrain.at(1)= Domega_Deqeps1*2*LocStrain.at(1)/eqeps1;
                Domega1_Dstrain.at(4)= Domega_Deqeps1*2*LocStrain.at(4)*pow(this->eqeps_1c/this->eqeps_1s,2)/eqeps1;
            }
            if (eqeps2 > 0){
                Domega2_Dstrain.at(1)= Domega_Deqeps2*2*LocStrain.at(2)/eqeps2;
                Domega2_Dstrain.at(4)= Domega_Deqeps2*2*LocStrain.at(4)*pow(this->eqeps_2c/this->eqeps_1s,2)/eqeps2;
            }
        };

        FloatMatrix MATC(4,4);
        nu21=this->E22*this->nu12/this->E11;
        D= (1-d2)*pow(this->nu23,2)+2*(1-d1)*(1-d2)*this->nu12*nu21*this->nu23+(1-d1)*(2-d2)*this->nu12*nu21-1;
        MATC.at(1,1)=this->E11*(1-d1)*((1-d2)*pow(this->nu23,2)-1)/D;
        MATC.at(1,2)=-this->E11*nu21*(1-d1)*(1-d2)*(1+this->nu23)/D;
        MATC.at(1,3)=-this->E11*nu21*(1-d1)*(1+(1-d2)*this->nu23)/D;
        MATC.at(1,4)=0;
        MATC.at(2,1)=MATC.at(1,2);
        MATC.at(2,2)=this->E22*(1-d2)*((1-d1)*this->nu12*nu21-1)/D;
        MATC.at(2,3)=-this->E22*(1-d2)*(this->nu23+(1-d1)*this->nu12*nu21)/D;
        MATC.at(2,4)=0;
        MATC.at(3,1)=MATC.at(1,3);
        MATC.at(3,2)=MATC.at(2,3);
        MATC.at(3,3)=this->E22*(1-d2)*(1-d1)*(this->nu12*nu21-1)/D;
        MATC.at(3,4)=0;
        MATC.at(4,1)=0;
        MATC.at(4,2)=0;
        MATC.at(4,3)=0;
        MATC.at(4,4)=this->G12*(1-d1)*(1-d2);


        FloatMatrix DMATC_Domega1(4,4);
        DMATC_Domega1.zero();
        double DDDomega1 = -2*(1-d2)*this->nu12*nu21*this->nu23-(2-d2)*this->nu12*nu21;
        DMATC_Domega1.at(1,1)=(-this->E11*((1-d2)*pow(this->nu23,2)-1)*D-DDDomega1*MATC.at(1,1))/pow(D,2);
        DMATC_Domega1.at(2,2)=(-this->E22*(1-d2)*this->nu12*nu21*D-DDDomega1*MATC.at(2,2))/pow(D,2);
        DMATC_Domega1.at(3,3)=(-this->E22*(1-d2)*(this->nu12*nu21-1)*D-DDDomega1*MATC.at(3,3))/pow(D,2);
        DMATC_Domega1.at(4,4)=-this->G12*(1-d2);
        DMATC_Domega1.at(1,2)=(this->E11*nu21*(1-d2)*(1+this->nu23)*D-DDDomega1*MATC.at(1,2))/pow(D,2);
        DMATC_Domega1.at(1,3)=(this->E11*nu21*(1+(1-d2)*this->nu23)*D-DDDomega1*MATC.at(1,3))/pow(D,2);
        DMATC_Domega1.at(2,3)= (this->E22*(1-d2)*this->nu12*nu21*D-DDDomega1*MATC.at(2,3))/pow(D,2);
        DMATC_Domega1.at(2,1)=DMATC_Domega1.at(1,2);
        DMATC_Domega1.at(3,1)=DMATC_Domega1.at(1,3);
        DMATC_Domega1.at(3,2)=DMATC_Domega1.at(2,3);

        FloatMatrix DMATC_Domega2(4,4);
        DMATC_Domega2.zero();
        double DDDomega2 = -pow(this->nu23,2)-2*(1-d1)*this->nu12*nu21*this->nu23+(1-d1)*this->nu12*nu21;
        DMATC_Domega2.at(1,1)=(-this->E11*(1-d1)*pow(this->nu23,2)*D-DDDomega2*MATC.at(1,1))/pow(D,2);
        DMATC_Domega2.at(2,2)=(-this->E22*((1-d1)*this->nu12*nu21-1)*D-DDDomega2*MATC.at(2,2))/pow(D,2);
        DMATC_Domega2.at(3,3)=(-this->E22*(1-d1)*(this->nu12*nu21-1)*D-DDDomega2*MATC.at(3,3))/pow(D,2);
        DMATC_Domega2.at(4,4)=-this->G12*(1-d1);
        DMATC_Domega2.at(1,2)=(this->E11*nu21*(1-d1)*(1+this->nu23)*D-DDDomega2*MATC.at(1,2))/pow(D,2);
        DMATC_Domega2.at(1,3)=(this->E11*nu21*(1-d1)*this->nu23*D-DDDomega2*MATC.at(1,3))/pow(D,2);
        DMATC_Domega2.at(2,3)= (this->E22*(this->nu23+(1-d1)*this->nu12*nu21)*D-DDDomega2*MATC.at(2,3))/pow(D,2);
        DMATC_Domega2.at(2,1)=DMATC_Domega2.at(1,2);
        DMATC_Domega2.at(3,1)=DMATC_Domega2.at(1,3);
        DMATC_Domega2.at(3,2)=DMATC_Domega2.at(2,3);

        FloatArray DcDeps1(4), DcDeps2(4);
        DcDeps1.beProductOf(DMATC_Domega1, LocStrain);
        DcDeps2.beProductOf(DMATC_Domega2, LocStrain);

        for (int i = 1; i < 4; i++){
            for(int j = 1; i < 4; i++){
                MATC.at(i,j) = MATC.at(i,j) + DcDeps1.at(i)*Domega1_Dstrain.at(j) + DcDeps2.at(i)*Domega2_Dstrain.at(j);
            }
        }

        FloatMatrix R(4,4),RPrim(4,4),Temp1(4,4),Temp2(4,4),Temp3(4,4);
        R.zero();
        R.at(1,1)=1;
        R.at(2,2)=1;
        R.at(3,3)=1;
        R.at(4,4)=2;
        RPrim.zero();
        RPrim.at(1,1)=1;
        RPrim.at(2,2)=1;
        RPrim.at(3,3)=1;
        RPrim.at(4,4)=0.5;

        Temp1.zero();
        Temp2.zero();
        Temp1.beProductOf(Trans,RPrim);
        Temp2.beProductOf(R,Temp1);
        Temp3.beProductOf(TransPrim,MATC);
        answer.resize(4,4);
        answer.beProductOf(Temp3,Temp2);
        return;
    }

}

int AnisotropicDamageMaterial2 ::giveDamageState(GaussPoint *gp)
{
    AnisotropicDamageMaterial2Status *status = static_cast< AnisotropicDamageMaterial2Status * >( this->giveStatus(gp) );

    FloatMatrix Damage(3,3);
    Damage = status->giveTempDamage();

    if (Damage.at(1,1) > 0.0 || Damage.at(1,2) > 0.0 ||Damage.at(2,1) > 0.0 || Damage.at(2,2) > 0.0){
        return 1;
    } else {
        return 0;
    };
}

IRResultType
AnisotropicDamageMaterial2 :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, E11, _IFT_AnisotropicDamageMaterial2_exx);
    IR_GIVE_FIELD(ir, E22, _IFT_AnisotropicDamageMaterial2_eyy);
    IR_GIVE_FIELD(ir, nu12, _IFT_AnisotropicDamageMaterial2_nuxy);
    IR_GIVE_FIELD(ir, nu23, _IFT_AnisotropicDamageMaterial2_nuyz);
    IR_GIVE_FIELD(ir, G12, _IFT_AnisotropicDamageMaterial2_gxy);
    IR_GIVE_FIELD(ir, eqeps_1t, _IFT_AnisotropicDamageMaterial2_eqeps_1t);
    IR_GIVE_FIELD(ir, eqeps_1s, _IFT_AnisotropicDamageMaterial2_eqeps_1s);
    IR_GIVE_FIELD(ir, alpha_1t, _IFT_AnisotropicDamageMaterial2_alpha_1t);
    IR_GIVE_FIELD(ir, eqeps_1c, _IFT_AnisotropicDamageMaterial2_eqeps_1c);
    IR_GIVE_FIELD(ir, beta_1c, _IFT_AnisotropicDamageMaterial2_beta_1c);
    IR_GIVE_FIELD(ir, alpha_1c, _IFT_AnisotropicDamageMaterial2_alpha_1c);
    IR_GIVE_FIELD(ir, eqeps_2t, _IFT_AnisotropicDamageMaterial2_eqeps_2t);
    IR_GIVE_FIELD(ir, eqeps_2s, _IFT_AnisotropicDamageMaterial2_eqeps_2s);
    IR_GIVE_FIELD(ir, alpha_2t, _IFT_AnisotropicDamageMaterial2_alpha_2t);
    IR_GIVE_FIELD(ir, eqeps_2c, _IFT_AnisotropicDamageMaterial2_eqeps_2c);
    IR_GIVE_FIELD(ir, beta_2c, _IFT_AnisotropicDamageMaterial2_beta_2c);
    IR_GIVE_FIELD(ir, alpha_2c, _IFT_AnisotropicDamageMaterial2_alpha_2c);
    IR_GIVE_FIELD(ir, eta, _IFT_AnisotropicDamageMaterial2_eta);
    IR_GIVE_FIELD(ir, theta, _IFT_AnisotropicDamageMaterial2_theta);

    return StructuralMaterial :: initializeFrom(ir);
}

int
AnisotropicDamageMaterial2 :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime)
{
    AnisotropicDamageMaterial2Status *status = static_cast< AnisotropicDamageMaterial2Status * >( this->giveStatus(gp) );
    if  (type == IST_DamageT1 ){
        answer.resize(1);
        answer.at(1) = status->givedamaget1();
        return 1;
    } else if ( type == IST_DamageT2 ) {
        answer.resize(1);
        answer.at(1) = status->givedamaget2();
        return 1;
    } else if ( type == IST_DamageC1 ) {
        answer.resize(1);
        answer.at(1) = status->givedamagec1();
        return 1;
    } else if ( type == IST_DamageC2 ) {
        answer.resize(1);
        answer.at(1) = status->givedamagec2();
        return 1;
    }else if ( type == IST_DamageTensorTemp ) {
        answer.resize(6);
        answer.zero();
        answer.at(1) = status->giveTempDamage().at(1, 1);
        answer.at(2) = status->giveTempDamage().at(1, 2);
        answer.at(3) = status->giveTempDamage().at(1, 3);
        answer.at(4) = status->giveTempDamage().at(2, 1);
        answer.at(5) = status->giveTempDamage().at(2, 2);
        answer.at(6) = status->giveTempDamage().at(2, 3);
        return 1;

    } else if ( type == IST_DamageTensor ) {
        answer.resize(6);
        answer.zero();
        answer.at(1) = status->giveDamage().at(1, 1);
        answer.at(2) = status->giveDamage().at(1, 2);
        answer.at(3) = status->giveDamage().at(1, 3);
        answer.at(4) = status->giveDamage().at(2, 1);
        answer.at(5) = status->giveDamage().at(2, 2);
        answer.at(6) = status->giveDamage().at(2, 3);
        return 1;

#ifdef keep_track_of_dissipated_energy
    } else if ( type == IST_StressWorkDensity ) {
        answer.resize(1);
        answer.at(1) = status->giveStressWork();
        return 1;
    } else if ( type == IST_DissWorkDensity ) {
        answer.resize(1);
        answer.at(1) = status->giveDissWork();
    } else if ( type == IST_FreeEnergyDensity ) {
        answer.resize(1);
        answer.at(1) = status->giveStressWork() - status->giveDissWork();
        return 1;

#endif
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, atTime);
    }

    return 1; // to make the compiler happy
}


/*/\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//\\\\\\\Material status functions\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/


AnisotropicDamageMaterial2Status :: AnisotropicDamageMaterial2Status(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g),
    RandomMaterialStatusExtensionInterface()
{
    equivalenstrain_1t = tempequivalenstrain_1t = 0.0;
    equivalenstrain_1c = tempequivalenstrain_1c = 0.0;
    equivalenstrain_2t = tempequivalenstrain_2t = 0.0;
    equivalenstrain_2c = tempequivalenstrain_2c = 0.0;
    equivalenstrain_3t = tempequivalenstrain_3t = 0.0;
    equivalenstrain_3c = tempequivalenstrain_3c = 0.0;

    damaget1 = 0.0;
    damaget2 = 0.0;
    damagec1 = 0.0;
    damagec2 = 0.0;

    damage.resize(3, 3);
    damage.zero();
    tempDamage.resize(3, 3);
    tempDamage.zero();

    flag = tempFlag = 0;

    flag_1st = 0;

#ifdef keep_track_of_dissipated_energy
    stressWork = tempStressWork = 0.0;
    dissWork = tempDissWork = 0.0;
#endif
}

AnisotropicDamageMaterial2Status :: ~AnisotropicDamageMaterial2Status()
{ }


void
AnisotropicDamageMaterial2Status :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    this->tempequivalenstrain_1t = this->equivalenstrain_1t;
    this->tempequivalenstrain_1c = this->equivalenstrain_1c;
    this->tempequivalenstrain_2t = this->equivalenstrain_2t;
    this->tempequivalenstrain_2c = this->equivalenstrain_2c;
    this->tempequivalenstrain_3t = this->equivalenstrain_3t;
    this->tempequivalenstrain_3c = this->equivalenstrain_3c;
    this->tempDamage = this->damage;
    this->tempFlag = this->flag;

#ifdef keep_track_of_dissipated_energy
    this->tempStressWork = this->stressWork;
    this->tempDissWork = this->dissWork;
#endif
}

void
AnisotropicDamageMaterial2Status :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);
    this->equivalenstrain_1t = this->tempequivalenstrain_1t;
    this->equivalenstrain_1c = this->tempequivalenstrain_1c;
    this->equivalenstrain_2t = this->tempequivalenstrain_2t;
    this->equivalenstrain_2c = this->tempequivalenstrain_2c;
    this->equivalenstrain_3t = this->tempequivalenstrain_3t;
    this->equivalenstrain_3c = this->tempequivalenstrain_3c;
    this->damage = this->tempDamage;
    this->flag = this->tempFlag;
#ifdef keep_track_of_dissipated_energy
    this->stressWork = this->tempStressWork;
    this->dissWork = this->tempDissWork;
#endif
}

Interface *
AnisotropicDamageMaterial2Status :: giveInterface(InterfaceType type)
{
    if ( type == RandomMaterialStatusExtensionInterfaceType ) {
        return static_cast< RandomMaterialStatusExtensionInterface * >(this);
    } else {
        return NULL;
    }
}

void
AnisotropicDamageMaterial2Status :: printOutputAt(FILE *file, TimeStep *tStep)
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

#ifdef keep_track_of_dissipated_energy
void
AnisotropicDamageMaterial2Status :: computeWork(GaussPoint *gp)
{

    // strain increment
     FloatArray deps;
     deps.beDifferenceOf(tempStrainVector, strainVector);
     
     // increment of stress work density
     double dSW = ( tempStressVector.dotProduct(deps) + stressVector.dotProduct(deps) ) / 2.;
     tempStressWork = stressWork + dSW;
     
     // elastically stored energy density
     double We = tempStressVector.dotProduct(tempStrainVector) / 2.;
     
     // dissipative work density
     tempDissWork = tempStressWork - We;
     
}
#endif

} // end namespace oofems
