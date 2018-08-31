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

#include "gausspoint.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "classfactory.h"
#include "intmatpprexcz.h"
#include "dynamicinputrecord.h"
#include <cmath>

namespace oofem {
REGISTER_Material(IntMatPPREXCZ);


IntMatPPREXCZ :: IntMatPPREXCZ(int n, Domain *d) : StructuralInterfaceMaterial(n, d)
  //
  // constructor
  //
{
    mPenaltyStiffness=0;
    mGIc=0;    // fracture energy, mode I
    mGIIc=0;   // fracture energy, mode II
    mSigma=0;  // max tensile stress
    mTau=0;    // max shear
    malpha=0;  //shape factors for tensile crack
    mbeta=0;   //shape factors for shearing crack
    mdeltan=0;    // Final normal crack open widths
    mdeltat=0;    // Final normal crack open widths
    mGamman=0; // Energy constant 1
    mGammat=0; // Energy constant 2
    dGnt=0;
    dGtn=0;
    mdeltan_conj=0;
    mdeltat_conj=0;
}


IntMatPPREXCZ :: ~IntMatPPREXCZ()
//
// destructor
//
{ }

int IntMatPPREXCZ :: checkConsistency()
{
    return 1;
}

void IntMatPPREXCZ :: giveFirstPKTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump,
                                                const FloatMatrix &F, TimeStep *tStep)
{
    this->initTempStatus(gp); //updata the status, so the current iteration starts from the previous converged results;

    // Note that the traction vector is assumed to be decomposed into its tangential and normal part,
    // with tangential directions (1,2) and normal direction (3).

    IntMatPPREXCZStatus *status = static_cast< IntMatPPREXCZStatus * >( this->giveStatus(gp) );

    FloatArray histo_jump=status->giveTempJumpMax();
    double delt_max = histo_jump.at(1);
    double deln_max = histo_jump.at(2);

    double delt= fabs(jump.at(2));
    double deln= jump.at(3);

    double sign_dt;
    if (jump.at(2) >= 0) {
        sign_dt=1;
    }else{
        sign_dt=-1;
    }

    FloatArray tractionTrial(3);
    tractionTrial.zero();

    if ( deln < 0 ){
        tractionTrial.at(3)=mPenaltyStiffness*deln;
    } else if ( deln >= 0 && deln <= mdeltan  && delt <= mdeltat_conj  ) {
        if (deln>=deln_max) {
            tractionTrial.at(3)=-1*malpha*(mGamman/mdeltan)*pow(1-deln/mdeltan, malpha-1)*
                    (mGammat*pow(1-delt/mdeltat, mbeta)+dGtn);
        }else{
            tractionTrial.at(3)=-1*malpha*(mGamman/mdeltan)*pow(1-deln_max/mdeltan, malpha-1)*
                    (mGammat*pow(1-delt/mdeltat, mbeta)+dGtn)*deln/deln_max;
        }
    } else if ( deln > mdeltan || delt > mdeltat_conj ){
        tractionTrial.at(3)=0;
    }

    if (delt == 0 ){
         tractionTrial.at(2)=0;
    }else if ( deln >= 0 && deln <= mdeltan_conj && delt <= mdeltat  ) {
        if (delt>=delt_max) {
            tractionTrial.at(2)=-1*mbeta*(mGammat/mdeltat)*pow(1-delt/mdeltat, mbeta-1)*
                    (mGamman*pow(1-deln/mdeltan, malpha)+dGnt)*sign_dt;
        }else{
            tractionTrial.at(2)=-1*mbeta*(mGammat/mdeltat)*pow(1-delt_max/mdeltat, mbeta-1)*
                    (mGamman*pow(1-deln/mdeltan, malpha)+dGnt)*delt*sign_dt/delt_max;
        }
    } else if ( delt > mdeltat || deln > mdeltan_conj ){
        tractionTrial.at(2)=0;
    }

    if (deln>=deln_max && deln<=mdeltan){
        histo_jump.at(2)=deln;
    }
    if (delt>=delt_max && delt<=mdeltat){
        histo_jump.at(1)=delt;
    }

    answer.resize(3);
    answer.zero();
    answer = tractionTrial;

    status->letTempJumpMaxBe(histo_jump);

    status->letTempJumpBe(jump);
    status->letTempFirstPKTractionBe(answer);
    status->letTempTractionBe(answer);

}

void IntMatPPREXCZ :: give3dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    //OOFEM_WARNING("not implemented. Use numerical Jacobian instead.");

    //answer={ss st sn;
    //        ts tt tn;
    //        ns nt nn;}
    IntMatPPREXCZStatus *status = static_cast< IntMatPPREXCZStatus * >( this->giveStatus(gp) );

    FloatArray histo_jump=status->giveTempJumpMax();
    double delt_max = histo_jump.at(1);
    double deln_max = histo_jump.at(2);

    FloatArray jump=status->giveTempJump();
    double delt= fabs(jump.at(2));
    double deln= jump.at(3);

    double sign_dt;
    if (jump.at(2) >= 0) {
        sign_dt=1;
    }else{
        sign_dt=-1;
    }

    FloatMatrix T_stiff;
    T_stiff.resize(2,2);
    T_stiff.zero();

    if (deln < 0){
        T_stiff.at(1,1)=mPenaltyStiffness;
        T_stiff.at(1,2)=0;
    }else if (deln >= 0 && deln <= mdeltan  && delt <= mdeltat_conj  ) {
        if (deln>=deln_max) {
            T_stiff.at(1,1)=malpha*(malpha-1)*(mGamman/pow(mdeltan,2))*pow(1-deln/mdeltan, malpha-2)*
                    (mGammat*pow(1-delt/mdeltat, mbeta)+dGtn);
            T_stiff.at(1,2)=malpha*mbeta*(mGamman/mdeltan)*(mGammat/mdeltat)*pow(1-deln/mdeltan, malpha-1)*
                    pow(1-delt/mdeltat, mbeta-1)*sign_dt;

        }else{
            T_stiff.at(1,1)=-1*malpha*(mGamman/mdeltan)*pow(1-deln_max/mdeltan, malpha-1)*
                    (mGammat*pow(1-delt/mdeltat, mbeta)+dGtn)/deln_max;
            T_stiff.at(1,2)=malpha*mbeta*(mGamman/mdeltan)*(mGammat/mdeltat)*pow(1-deln/mdeltan, malpha-1)*
                    pow(1-delt/mdeltat, mbeta-1)*sign_dt*deln/deln_max;
        }
    } else {
        T_stiff.at(1,1)=0;
        T_stiff.at(1,2)=0;
    }

    if ( deln >= 0 && deln <= mdeltan_conj && delt <= mdeltat  ) {
        if (delt>=delt_max) {
            T_stiff.at(2,1)=malpha*mbeta*(mGamman/mdeltan)*(mGammat/mdeltat)*pow(1-deln/mdeltan, malpha-1)*
                    pow(1-delt/mdeltat, mbeta-1)*sign_dt;
            T_stiff.at(2,2)=mbeta*(mbeta-1)*(mGammat/pow(mdeltat,2))*pow(1-delt/mdeltat, mbeta-2)*
                    (mGamman*pow(1-deln/mdeltan, malpha)+dGnt);
        }else{
            T_stiff.at(2,1)=malpha*mbeta*(mGamman/mdeltan)*(mGammat/mdeltat)*pow(1-deln/mdeltan, malpha-1)*
                    pow(1-delt/mdeltat, mbeta-1)*sign_dt*delt/delt_max;
            T_stiff.at(2,2)=-1*mbeta*(mGammat/mdeltat)*pow(1-delt_max/mdeltat, mbeta-1)*
                    (mGamman*pow(1-deln/mdeltan, malpha)+dGnt)*sign_dt/delt_max;
        }
    } else if ( delt > mdeltat || deln > mdeltan_conj ){
        T_stiff.at(2,1)=0;
        T_stiff.at(2,2)=0;
    }

    answer.resize(3,3);
    answer.zero();
    answer.at(2,2)=T_stiff.at(2,2);
    answer.at(2,3)=T_stiff.at(2,1);
    answer.at(3,2)=T_stiff.at(1,2);
    answer.at(3,3)=T_stiff.at(1,1);

}

int IntMatPPREXCZ :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    IntMatPPREXCZStatus *status = static_cast< IntMatPPREXCZStatus * >( this->giveStatus(gp) );
    if ( type == IST_DamageScalar ) {
        double damage1=status->giveTempJumpMax().at(1)/mdeltat;
        double damage2=status->giveTempJumpMax().at(2)/mdeltan;
        answer.resize(1);
        if (damage1>=damage2){
            answer.at(1) = damage1;
        }else{
            answer.at(1) = damage2;
        }
        return 1;
    } else {
        return StructuralInterfaceMaterial :: giveIPValue(answer, gp, type, tStep);
    }
}

IRResultType IntMatPPREXCZ :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                    // Required by IR_GIVE_FIELD macro

    //IR_GIVE_FIELD(ir, mPenaltyStiffness, _IFT_IntMatPPREXCZ_PenaltyStiffness);

    IR_GIVE_FIELD(ir, mGIc, _IFT_IntMatPPREXCZ_phi1);

    mGIIc = mGIc;                                   //Defaults to GIc
    IR_GIVE_OPTIONAL_FIELD(ir, mGIIc, _IFT_IntMatPPREXCZ_phi2);

    IR_GIVE_FIELD(ir, mSigma, _IFT_IntMatPPREXCZ_sigm);

    IR_GIVE_FIELD(ir, mTau, _IFT_IntMatPPREXCZ_taum);

    IR_GIVE_FIELD(ir, malpha, _IFT_IntMatPPREXCZ_alpha);

    IR_GIVE_FIELD(ir, mbeta, _IFT_IntMatPPREXCZ_beta);

    mdeltan=malpha*mGIc/mSigma;

    mdeltat=mbeta*mGIIc/mTau;

    mPenaltyStiffness=10*mSigma/mdeltan;

    if ( mGIc < mGIIc ){
        dGnt=0;
        dGtn=mGIIc - mGIc;
    }else if(mGIc > mGIIc){
        dGnt=mGIc - mGIIc;
        dGtn=0;
    }else{
        dGnt=0;
        dGtn=0;
    };

    mdeltan_conj=mdeltan-mdeltan*pow(dGnt/mGIc,1/malpha);
    mdeltat_conj=mdeltat-mdeltat*pow(dGtn/mGIIc,1/mbeta);

    if ( mGIc == mGIIc ){
        mGamman = -mGIc;
        mGammat =1;
    } else {
        mGamman =pow(-mGIc, dGnt/(mGIc - mGIIc));
        mGammat =pow(-mGIIc, dGtn/(mGIIc - mGIc));
    };
    return StructuralInterfaceMaterial :: initializeFrom(ir);
}

void IntMatPPREXCZ :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralInterfaceMaterial :: giveInputRecord(input);

    input.setField(mGIc, _IFT_IntMatPPREXCZ_phi1);
    input.setField(mGIIc, _IFT_IntMatPPREXCZ_phi2);
    input.setField(mSigma, _IFT_IntMatPPREXCZ_sigm);
    input.setField(mTau, _IFT_IntMatPPREXCZ_taum);
    input.setField(malpha, _IFT_IntMatPPREXCZ_alpha);
    input.setField(mbeta, _IFT_IntMatPPREXCZ_beta);

}

void IntMatPPREXCZ :: printYourself()
{
    printf("\nInitializing IntMatPPREXCZ:\n");
    printf("mPenaltyStiffness: %e\n", mPenaltyStiffness);
    printf("mGIc: %e\n", mGIc);
    printf("mGIIc: %e\n", mGIIc);
    printf("mSigma: %e\n", mSigma);
    printf("mTau: %e\n", mTau);;
    printf("malpha: %e\n\n", malpha);
    printf("mbeta: %e\n\n", mbeta);
}

/*/\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//\\\\\\\Material status functions\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/


IntMatPPREXCZStatus :: IntMatPPREXCZStatus(int n, Domain *d, GaussPoint *g) : StructuralInterfaceMaterialStatus(n, d, g)
  //
  // constructor
  //
{
   JumpMax.resize(2);
   JumpMax.zero();
   tempJumpMax.resize(2);
   tempJumpMax.zero();
}

IntMatPPREXCZStatus :: ~IntMatPPREXCZStatus()
//
// destructor
//
{ }

void IntMatPPREXCZStatus :: initTempStatus()
{
    StructuralInterfaceMaterialStatus :: initTempStatus();

    this->tempJumpMax = this->JumpMax;
}

void IntMatPPREXCZStatus :: updateYourself(TimeStep *tStep)
{
    JumpMax = tempJumpMax;;

    StructuralInterfaceMaterialStatus ::updateYourself(tStep);
}


void IntMatPPREXCZStatus :: copyStateVariables(const MaterialStatus &iStatus)
{
    StructuralInterfaceMaterialStatus :: copyStateVariables(iStatus);

    MaterialStatus &tmpStat = const_cast< MaterialStatus & >(iStatus);
    const IntMatPPREXCZStatus &structStatus = dynamic_cast< IntMatPPREXCZStatus & >(tmpStat);

    JumpMax   = structStatus.JumpMax;
    tempJumpMax   = structStatus.tempJumpMax;

}

void IntMatPPREXCZStatus :: addStateVariables(const MaterialStatus &iStatus)
{
    OOFEM_ERROR("not implemented.");
}


} /* namespace oofem */
