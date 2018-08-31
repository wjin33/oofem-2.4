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
#include "intmatpprincz.h"
#include "dynamicinputrecord.h"
#include <cmath>

namespace oofem {
REGISTER_Material(IntMatPPRINCZ);


IntMatPPRINCZ :: IntMatPPRINCZ(int n, Domain *d) : StructuralInterfaceMaterial(n, d)
  //
  // constructor
  //
{
    PenaltyStiffness=0;
    GIc=0;    // fracture energy, mode I
    GIIc=0;   // fracture energy, mode II
    Sigma=0;  // max tensile stress
    Tau=0;    // max shear
    alpha=0;  //shape factors for tensile crack
    beta=0;   //shape factors for shearing crack
    deltan=0;    // Final normal crack open widths
    deltat=0;    // Final normal crack open widths
    Gamman=0; // Energy constant 1
    Gammat=0; // Energy constant 2
    dGnt=0;
    dGtn=0;
    deltan_conj=0;
    deltat_conj=0;
    lambdan=0; //normal penalty stiffness
    lambdat=0; //tangent penalty stiffness
    m=0;
    n=0;
}


IntMatPPRINCZ :: ~IntMatPPRINCZ()
//
// destructor
//
{ }

int IntMatPPRINCZ :: checkConsistency()
{
    return 1;
}

void IntMatPPRINCZ :: giveFirstPKTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump,
                                                const FloatMatrix &F, TimeStep *tStep)
{
    this->initTempStatus(gp); //updata the status, so the current iteration starts from the previous converged results;

    // Note that the traction vector is assumed to be decomposed into its tangential and normal part,
    // with tangential directions (1,2) and normal direction (3).

    IntMatPPRINCZStatus *status = static_cast< IntMatPPRINCZStatus * >( this->giveStatus(gp) );

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
        tractionTrial.at(3)=PenaltyStiffness*deln;
    } else if ( deln >= 0 && deln <= deltan  && delt <= deltat_conj  ) {
        if (deln>=deln_max) {
            tractionTrial.at(3)=(Gamman/deltan)*(m*pow(1-deln/deltan, alpha)*pow(m/alpha+deln/deltan, m-1)-
                                alpha*pow(1-deln/deltan, alpha-1)*pow(m/alpha+deln/deltan, m))*
                                (Gammat*pow(1-delt/deltat, beta)*pow(n/beta+delt/deltat, n)+dGtn);
        }else{
            tractionTrial.at(3)=(Gamman/deltan)*(m*pow(1-deln_max/deltan, alpha)*pow(m/alpha+deln_max/deltan, m-1)-
                                alpha*pow(1-deln_max/deltan, alpha-1)*pow(m/alpha+deln_max/deltan, m))*
                                (Gammat*pow(1-delt/deltat, beta)*pow(n/beta+delt/deltat, n)+dGtn)*deln/deln_max;
        }
    } else if ( deln > deltan || delt > deltat_conj ){
        tractionTrial.at(3)=0;
    }

    if ( deln >= 0 && deln <= deltan_conj && delt <= deltat  ) {
        if (delt>=delt_max) {
            tractionTrial.at(2)=(Gammat/deltat)*(n*pow(1-delt/deltat, beta)*pow(n/beta+delt/deltat, n-1)-
                                beta*pow(1-delt/deltat, beta-1)*pow(n/beta+delt/deltat, n))*
                                (Gamman*pow(1-deln/deltan, alpha)*pow(m/alpha+deln/deltan, m)+dGnt)*sign_dt;
        }else{
            tractionTrial.at(2)=(Gammat/deltat)*(n*pow(1-delt_max/deltat, beta)*pow(n/beta+delt_max/deltat, n-1)-
                                beta*pow( 1-delt_max/deltat, beta-1)*pow(n/beta+delt_max/deltat, n))*
                                (Gamman*pow(1-deln/deltan, alpha)*pow(m/alpha+deln/deltan, m)+dGnt)*delt*sign_dt/delt_max;
        }
    } else if ( delt > deltat || deln > deltan_conj ){
        tractionTrial.at(2)=0;
    }

    if (deln>=deln_max && deln<=deltan){
        histo_jump.at(2)=deln;
    }
    if (delt>=delt_max && delt<=deltat){
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

void IntMatPPRINCZ :: give3dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    //OOFEM_WARNING("not implemented. Use numerical Jacobian instead.");

    //answer={ss st sn;
    //        ts tt tn;
    //        ns nt nn;}
    IntMatPPRINCZStatus *status = static_cast< IntMatPPRINCZStatus * >( this->giveStatus(gp) );

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
        T_stiff.at(1,1)=PenaltyStiffness;
        T_stiff.at(1,2)=0;
    }else if (deln >= 0 && deln <= deltan  && delt <= deltat_conj ) {
        if (deln>=deln_max) {
            T_stiff.at(1,1) = Gamman/pow(deltan,2)*((pow(m,2)-m)*pow(1-deln/deltan,alpha)*pow(m/alpha+deln/deltan,m-2)+
                (pow(alpha,2)-alpha)*pow(1-deln/deltan,alpha-2)*pow(m/alpha+deln/deltan,m)-
                2*alpha*m*pow(1-deln/deltan,alpha-1)*pow(m/alpha+deln/deltan,m-1))*
                (Gammat*pow(1-delt/deltat,beta)*pow(n/beta+delt/deltat,n)+dGtn);
            T_stiff.at(1,2)=Gamman*Gammat/(deltat*deltan)*(m*pow(1-deln/deltan,alpha)*pow(m/alpha+deln/deltan,m-1)-
            	alpha*pow(1-deln/deltan, alpha-1)*pow(m/alpha+deln/deltan,m))*
                (n*pow(1-delt/deltat,beta)*pow(n/beta+delt/deltat,n-1)-
                beta*pow(1-delt/deltat,beta-1)*pow(n/beta+delt/deltat,n))*sign_dt;
        }else{
            T_stiff.at(1,1)=(Gamman/deltan)*(m*pow(1-deln_max/deltan, alpha)*pow(m/alpha+deln_max/deltan, m-1)-
                alpha*pow(1-deln_max/deltan, alpha-1)*pow(m/alpha+deln_max/deltan, m))*
                (Gammat*pow(1-delt/deltat, beta)*pow(n/beta+delt/deltat, n)+dGtn)/deln_max;
            T_stiff.at(1,2)=Gamman*Gammat/(deltat*deltan)*(m*pow(1-deln_max/deltan,alpha)*pow(m/alpha+deln_max/deltan,m-1)-
            	alpha*pow(1-deln_max/deltan, alpha-1)*pow(m/alpha+deln_max/deltan,m))*
                (n*pow(1-delt/deltat,beta)*pow(n/beta+delt/deltat,n-1)-
                beta*pow(1-delt/deltat,beta-1)*pow(n/beta+delt/deltat,n))*sign_dt*deln/deln_max;
        }
    } else {
        T_stiff.at(1,1)=0;
        T_stiff.at(1,2)=0;
    }

    if ( deln >= 0 && deln <= deltan_conj && delt <= deltat  ) {
        if (delt>=delt_max) {
            T_stiff.at(2,1) = T_stiff.at(1,2);
            T_stiff.at(2,2) = Gammat/pow(deltat,2)*((pow(n,2)-n)*pow(1-delt/deltat,beta)*pow(n/beta+delt/deltat,n-2)+
            	(pow(beta,2)-beta)*pow(1-delt/deltat, beta-2)*pow(n/beta+delt/deltat,n)-
            	2*beta*n*pow(1-delt/deltat,beta-1)*pow(n/beta+delt/deltat,n-1))*
            	(Gamman*pow(1-deln/deltan, alpha)*pow(m/alpha+deln/deltan,m)+dGnt);
        }else{
            T_stiff.at(2,1)=Gamman*Gammat/(deltat*deltan)*(m*pow(1-deln/deltan,alpha)*pow(m/alpha+deln/deltan,m-1)-
            	alpha*pow(1-deln/deltan, alpha-1)*pow(m/alpha+deln/deltan,m))*
                (n*pow(1-delt_max/deltat,beta)*pow(n/beta+delt_max/deltat,n-1)-
                beta*pow(1-delt_max/deltat,beta-1)*pow(n/beta+delt_max/deltat,n))*sign_dt*delt/delt_max;
            T_stiff.at(2,2)= (Gammat/deltat)*(n*pow(1-delt_max/deltat, beta)*pow(n/beta+delt_max/deltat, n-1)-
                beta*pow( 1-delt_max/deltat, beta-1)*pow(n/beta+delt_max/deltat, n))*
                (Gamman*pow(1-deln/deltan, alpha)*pow(m/alpha+deln/deltan, m)+dGnt)/delt_max;
        }
    } else if ( delt > deltat || deln > deltan_conj ){
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

int IntMatPPRINCZ :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    IntMatPPRINCZStatus *status = static_cast< IntMatPPRINCZStatus * >( this->giveStatus(gp) );
    if ( type == IST_DamageScalar ) {
        double damage1=status->giveTempJumpMax().at(1)/deltat;
        double damage2=status->giveTempJumpMax().at(2)/deltan;
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

IRResultType IntMatPPRINCZ :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                    // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, lambdan, _IFT_IntMatPPRINCZ_lambdan);
    IR_GIVE_FIELD(ir, lambdat, _IFT_IntMatPPRINCZ_lambdat);

    IR_GIVE_FIELD(ir, GIc, _IFT_IntMatPPRINCZ_phi1);
    GIIc = GIc;                                   //Defaults to GIc
    IR_GIVE_OPTIONAL_FIELD(ir, GIIc, _IFT_IntMatPPRINCZ_phi2);

    IR_GIVE_FIELD(ir, Sigma, _IFT_IntMatPPRINCZ_sigm);
    IR_GIVE_FIELD(ir, Tau, _IFT_IntMatPPRINCZ_taum);

    IR_GIVE_FIELD(ir, alpha, _IFT_IntMatPPRINCZ_alpha);
    IR_GIVE_FIELD(ir, beta, _IFT_IntMatPPRINCZ_beta);

    m = alpha*(alpha-1)*pow(lambdan, 2)/(1-alpha*pow(lambdan,2));
    n = beta*(beta-1)*pow(lambdat, 2)/(1-beta*pow(lambdat,2));

    deltan = GIc/Sigma*alpha*lambdan*pow(1-lambdan, alpha-1)*(alpha/m+1)*pow(alpha*lambdan/m+1,m-1);
    deltat = GIIc/Tau*beta*lambdat*pow(1-lambdat, beta-1)*(beta/n+1)*pow(beta*lambdat/n+1,n-1);

    PenaltyStiffness=10*Sigma/deltan;

    if ( GIc < GIIc ){
        dGnt=0;
        dGtn=GIIc - GIc;
    }else if(GIc > GIIc){
        dGnt=GIc - GIIc;
        dGtn=0;
    }else{
        dGnt=0;
        dGtn=0;
    };

    if ( GIc == GIIc ){
        Gamman = -GIc*pow(alpha/m, m);
        Gammat = pow(beta/n, n);
    } else {
        Gamman =pow(-GIc, dGnt/(GIc - GIIc))*pow(alpha/m, m);
        Gammat =pow(-GIIc, dGtn/(GIIc - GIc))*pow(beta/n, n);
    };

    double x0=0;
    double x1=deltat;
    double tol=0.5*10e-7;
    double xm, f_x0, f_xm;

    if (GIc < GIIc){
        while ((x1-x0)>tol){
            xm=x0+(x1-x0)/2;
            f_x0=Gammat*(pow(1-fabs(x0)/deltat,beta))*(pow((n/beta)+fabs(x0)/deltat,n))+dGtn;
            f_xm=Gammat*(pow(1-fabs(xm)/deltat,beta))*(pow((n/beta)+fabs(xm)/deltat,n))+dGtn;
            if (sign(f_x0)==sign(f_xm)){
                x0=xm;
            }else{
                x1=xm;
            }
        }
        deltat_conj=xm;
    }else {
        deltat_conj=deltat;
    }

    x0=0;
    x1=deltan;
    if (GIc > GIIc){
        while ((x1-x0)>tol){
            xm=x0+(x1-x0)/2;
            f_x0=Gamman*(pow(1-fabs(x0)/deltan,alpha))*(pow((m/alpha)+fabs(x0)/deltan,m))+dGnt;
            f_xm=Gamman*(pow(1-fabs(xm)/deltan,alpha))*(pow((m/alpha)+fabs(xm)/deltan,m))+dGnt;
            if (sign(f_x0)==sign(f_xm)){
                x0=xm;
            }else{
                x1=xm;
            }
        }
        deltan_conj=xm;
    }else{
        deltan_conj=deltan;
    }


    return StructuralInterfaceMaterial :: initializeFrom(ir);
}

double IntMatPPRINCZ :: sign (double x ){
    if (x>0){
        return 1;
    }else if (x<0){
        return -1;
    }else{
        return 0;
    }
}

void IntMatPPRINCZ :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralInterfaceMaterial :: giveInputRecord(input);

    input.setField(GIc, _IFT_IntMatPPRINCZ_phi1);
    input.setField(GIIc, _IFT_IntMatPPRINCZ_phi2);
    input.setField(Sigma, _IFT_IntMatPPRINCZ_sigm);
    input.setField(Tau, _IFT_IntMatPPRINCZ_taum);
    input.setField(alpha, _IFT_IntMatPPRINCZ_alpha);
    input.setField(beta, _IFT_IntMatPPRINCZ_beta);
    input.setField(lambdan, _IFT_IntMatPPRINCZ_lambdan);
    input.setField(lambdat, _IFT_IntMatPPRINCZ_lambdat);

}

void IntMatPPRINCZ :: printYourself()
{
    printf("\nInitializing IntMatPPRINCZ:\n");
    printf("PenaltyStiffness: %e\n", PenaltyStiffness);
    printf("GIc: %e\n", GIc);
    printf("GIIc: %e\n", GIIc);
    printf("Sigma: %e\n", Sigma);
    printf("Tau: %e\n", Tau);;
    printf("alpha: %e\n\n", alpha);
    printf("beta: %e\n\n", beta);
}

/*/\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//\\\\\\\Material status functions\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/


IntMatPPRINCZStatus :: IntMatPPRINCZStatus(int n, Domain *d, GaussPoint *g) : StructuralInterfaceMaterialStatus(n, d, g)
  //
  // constructor
  //
{
   JumpMax.resize(2);
   JumpMax.zero();
   tempJumpMax.resize(2);
   tempJumpMax.zero();
}

IntMatPPRINCZStatus :: ~IntMatPPRINCZStatus()
//
// destructor
//
{ }

void IntMatPPRINCZStatus :: initTempStatus()
{
    StructuralInterfaceMaterialStatus :: initTempStatus();

    this->tempJumpMax = this->JumpMax;
}

void IntMatPPRINCZStatus :: updateYourself(TimeStep *tStep)
{
    JumpMax = tempJumpMax;;

    StructuralInterfaceMaterialStatus ::updateYourself(tStep);
}


void IntMatPPRINCZStatus :: copyStateVariables(const MaterialStatus &iStatus)
{
    StructuralInterfaceMaterialStatus :: copyStateVariables(iStatus);

    MaterialStatus &tmpStat = const_cast< MaterialStatus & >(iStatus);
    const IntMatPPRINCZStatus &structStatus = dynamic_cast< IntMatPPRINCZStatus & >(tmpStat);

    JumpMax   = structStatus.JumpMax;
    tempJumpMax   = structStatus.tempJumpMax;

}

void IntMatPPRINCZStatus :: addStateVariables(const MaterialStatus &iStatus)
{
    OOFEM_ERROR("not implemented.");
}


} /* namespace oofem */
