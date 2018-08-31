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
#include "intmatpprinadm1.h"
#include "dynamicinputrecord.h"
#include <cmath>

namespace oofem {
REGISTER_Material(IntMatPPRINADM1);


IntMatPPRINADM1 :: IntMatPPRINADM1(int n, Domain *d) : StructuralInterfaceMaterial(n, d)
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
    TotalG=0;   //Total energy release rate to be dissipated
}


IntMatPPRINADM1 :: ~IntMatPPRINADM1()
//
// destructor
//
{ }

int IntMatPPRINADM1 :: checkConsistency()
{
    return 1;
}

void IntMatPPRINADM1 :: giveFirstPKTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump,
                                                const FloatMatrix &F, TimeStep *tStep)
{
    this->initTempStatus(gp); //updata the status, so the current iteration starts from the previous converged results;

    // Note that the traction vector is assumed to be decomposed into its tangential and normal part,
    // with tangential directions (1,2) and normal direction (3).

    IntMatPPRINADM1Status *status = static_cast< IntMatPPRINADM1Status * >( this->giveStatus(gp) );

    FloatArray MaterialPara = status->delivertomaterial();
    PenaltyStiffness = MaterialPara.at(1);
    GIc = MaterialPara.at(2);
    GIIc = MaterialPara.at(3);
    Sigma = MaterialPara.at(4);
    Tau = MaterialPara.at(5);
    deltan = MaterialPara.at(6);
    deltat = MaterialPara.at(7);
    Gamman = MaterialPara.at(8);
    Gammat = MaterialPara.at(9);
    dGnt = MaterialPara.at(10);
    dGtn = MaterialPara.at(11);
    deltan_conj = MaterialPara.at(12);
    deltat_conj = MaterialPara.at(13);
    m = MaterialPara.at(14);
    n = MaterialPara.at(15);

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

void IntMatPPRINADM1 :: give3dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep)
{
    //OOFEM_WARNING("not implemented. Use numerical Jacobian instead.");

    //answer={ss st sn;
    //        ts tt tn;
    //        ns nt nn;}
    IntMatPPRINADM1Status *status = static_cast< IntMatPPRINADM1Status * >( this->giveStatus(gp) );

    FloatArray MaterialPara = status->delivertomaterial();
    PenaltyStiffness = MaterialPara.at(1);
    GIc = MaterialPara.at(2);
    GIIc = MaterialPara.at(3);
    Sigma = MaterialPara.at(4);
    Tau = MaterialPara.at(5);
    deltan = MaterialPara.at(6);
    deltat = MaterialPara.at(7);
    Gamman = MaterialPara.at(8);
    Gammat = MaterialPara.at(9);
    dGnt = MaterialPara.at(10);
    dGtn = MaterialPara.at(11);
    deltan_conj = MaterialPara.at(12);
    deltat_conj = MaterialPara.at(13);
    m = MaterialPara.at(14);
    n = MaterialPara.at(15);

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

int IntMatPPRINADM1 :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep)
{
    IntMatPPRINADM1Status *status = static_cast< IntMatPPRINADM1Status * >( this->giveStatus(gp) );
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

IRResultType IntMatPPRINADM1 :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                    // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, lambdan, _IFT_IntMatPPRINADM1_lambdan);
    IR_GIVE_FIELD(ir, lambdat, _IFT_IntMatPPRINADM1_lambdat);

    IR_GIVE_FIELD(ir, alpha, _IFT_IntMatPPRINADM1_alpha);
    IR_GIVE_FIELD(ir, beta, _IFT_IntMatPPRINADM1_beta);

    IR_GIVE_FIELD(ir, TotalG, _IFT_IntMatPPRINADM1_totalg);

    return StructuralInterfaceMaterial :: initializeFrom(ir);
}

void IntMatPPRINADM1 :: giveInputRecord(DynamicInputRecord &input)
{
    StructuralInterfaceMaterial :: giveInputRecord(input);

    input.setField(alpha, _IFT_IntMatPPRINADM1_alpha);
    input.setField(beta, _IFT_IntMatPPRINADM1_beta);
    input.setField(lambdan, _IFT_IntMatPPRINADM1_lambdan);
    input.setField(lambdat, _IFT_IntMatPPRINADM1_lambdat);
    input.setField(TotalG, _IFT_IntMatPPRINADM1_totalg);

}

void IntMatPPRINADM1 :: printYourself()
{
    printf("\nInitializing IntMatPPRINADM1:\n");
    printf("PenaltyStiffness: %e\n", PenaltyStiffness);
    printf("alpha: %e\n\n", alpha);
    printf("beta: %e\n\n", beta);
}


/*/\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//\\\\\\\Material status functions\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/


IntMatPPRINADM1Status :: IntMatPPRINADM1Status(int n, Domain *d, GaussPoint *g, const double &a1, const double &a2, const double &a3, const double &a4, const double &a5) : StructuralInterfaceMaterialStatus(n, d, g)
  //
  // constructor
  //
{
    JumpMax.resize(2);
    JumpMax.zero();
    tempJumpMax.resize(2);
    tempJumpMax.zero();
    PenaltyStiffness=0;
    GIc=0;    // fracture energy, mode I
    GIIc=0;   // fracture energy, mode II
    Sigma=0;  // max tensile stress
    Tau=0;    // max shear
    alpha=a1;  //shape factors for tensile crack
    beta=a2;   //shape factors for shearing crack
    deltan=0;    // Final normal crack open widths
    deltat=0;    // Final normal crack open widths
    Gamman=0; // Energy constant 1
    Gammat=0; // Energy constant 2
    dGnt=0;
    dGtn=0;
    deltan_conj=0;
    deltat_conj=0;
    lambdan=a3; //normal penalty stiffness
    lambdat=a4; //tangent penalty stiffness
    m=0;
    n=0;
    TotalG=a5;   //Total energy release rate to be dissipated
}

IntMatPPRINADM1Status :: ~IntMatPPRINADM1Status()
//
// destructor
//
{ }

FloatArray IntMatPPRINADM1Status ::delivertomaterial()
{
    FloatArray materialparameter(15);
    materialparameter.at(1) = PenaltyStiffness;
    materialparameter.at(2) = GIc;
    materialparameter.at(3) = GIIc;
    materialparameter.at(4) = Sigma;
    materialparameter.at(5) = Tau;
    materialparameter.at(6) = deltan;
    materialparameter.at(7) = deltat;
    materialparameter.at(8) = Gamman;
    materialparameter.at(9) = Gammat;
    materialparameter.at(10) = dGnt;
    materialparameter.at(11) = dGtn;
    materialparameter.at(12) = deltan_conj;
    materialparameter.at(13) = deltat_conj;
    materialparameter.at(14) = m;
    materialparameter.at(15) = n;
    return materialparameter;
}

void IntMatPPRINADM1Status :: initTempStatus()
{
    StructuralInterfaceMaterialStatus :: initTempStatus();

    this->tempJumpMax = this->JumpMax;
}

void IntMatPPRINADM1Status :: updateYourself(TimeStep *tStep)
{
    JumpMax = tempJumpMax;;

    StructuralInterfaceMaterialStatus ::updateYourself(tStep);
}


void IntMatPPRINADM1Status :: copyStateVariables(const MaterialStatus &iStatus)
{
    StructuralInterfaceMaterialStatus :: copyStateVariables(iStatus);

    MaterialStatus &tmpStat = const_cast< MaterialStatus & >(iStatus);
    const IntMatPPRINADM1Status &structStatus = dynamic_cast< IntMatPPRINADM1Status & >(tmpStat);

    JumpMax   = structStatus.JumpMax;
    tempJumpMax   = structStatus.tempJumpMax;
    PenaltyStiffness =structStatus.PenaltyStiffness;
    GIc = structStatus.GIc;
    GIIc = structStatus.GIIc;
    Sigma = structStatus.Sigma;
    Tau = structStatus.Tau;
    alpha = structStatus.alpha;
    beta = structStatus.beta;
    deltan = structStatus.deltan;
    deltat = structStatus.deltat;
    Gamman = structStatus.Gamman;
    Gammat = structStatus.Gammat;
    dGnt = structStatus.dGnt;
    dGtn = structStatus.dGtn;
    deltan_conj = structStatus.deltan_conj;
    deltat_conj = structStatus.deltat_conj;
    lambdan = structStatus.lambdan;
    lambdat = structStatus.lambdat;
    m = structStatus.m;
    n = structStatus.n;
    TotalG = structStatus.TotalG;

}

void IntMatPPRINADM1Status :: addStateVariables(const MaterialStatus &iStatus)
{
    OOFEM_ERROR("not implemented.");
} 

void IntMatPPRINADM1Status :: initializeFrom(const double &NeedDissipated_G)
{
    GIc= NeedDissipated_G;
    GIIc = GIc;                    //Defaults to GIc


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

}

double IntMatPPRINADM1Status :: sign (double x ){
    if (x>0){
        return 1;
    }else if (x<0){
        return -1;
    }else{
        return 0;
    }
}


} /* namespace oofem */
