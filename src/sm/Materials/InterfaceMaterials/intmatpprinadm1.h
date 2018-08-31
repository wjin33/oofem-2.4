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

#ifndef INTMatPPRINADM1_H_
#define INTMatPPRINADM1_H_

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

///@name Input fields for IntMatPPRINADM1
//@{
#define _IFT_IntMatPPRINADM1_Name "intmatpprinadm1"
//#define _IFT_IntMatPPRINADM1_PenaltyStiffness "kn"
#define _IFT_IntMatPPRINADM1_totalg  "totalg"
#define _IFT_IntMatPPRINADM1_alpha "alpha"
#define _IFT_IntMatPPRINADM1_beta  "beta"
#define _IFT_IntMatPPRINADM1_lambdan "lambdan"
#define _IFT_IntMatPPRINADM1_lambdat  "lambdat"
//@}

namespace oofem {
/**
 * This class implements associated Material Status for IntMatPPRINADM1
 */
class IntMatPPRINADM1Status : public StructuralInterfaceMaterialStatus
{
public:
    IntMatPPRINADM1Status(int n, Domain * d, GaussPoint * g, const double &a1, const double &a2, const double &a3, const double &a4, const double &a5);
    virtual ~IntMatPPRINADM1Status();

    /// history variable
    /// Vector measure of the largest displacement jump ever reached in material.
    FloatArray JumpMax;
    /// Non-equilibrated vector measure of the largest displacement jump ever reached in material.
    FloatArray tempJumpMax;

    double PenaltyStiffness;
    double GIc;    // fracture energy, mode I
    double GIIc;   // fracture energy, mode II
    double Sigma;  // max tensile stress
    double Tau;    // max shear
    double alpha;  //shape factors for tensile crack
    double beta;   //shape factors for shearing crack
    double lambdan; //normal penalty stiffness
    double lambdat; //tangent penalty stiffness
    double m;
    double n;

    double deltan;       // Final normal crack open widths
    double deltat;       // Final shear crack displacement
    double deltan_conj;  // Conjugated final normal crack widths which shear traction is zero
    double deltat_conj;  // Conjugated final shear crack widths which tensile traction is zero

    double Gamman; // Energy constant 1
    double Gammat; // Energy constant 2
    double dGtn;    // Energy constant 3
    double dGnt;    // Energy constant 4
    double TotalG;   //Total energy release rate to be dissipated

    /// Returns the last equilibrated vector measure largest displacement jump.
    const FloatArray &giveJumpMax() { return JumpMax; }
    /// Returns the temp. vector measure largest displacement jump.
    const FloatArray &giveTempJumpMax() { return tempJumpMax; }
    /// Sets the vector measure of the displacement jump to given value.
    void letJumpMaxBe(FloatArray newJumpMax) { JumpMax = std :: move(newJumpMax); }
    /// Sets the temp vector measure of the displacement jump to given value.
    void letTempJumpMaxBe(FloatArray newTempJumpMax) { tempJumpMax = std :: move(newTempJumpMax); }
    void letTranctionBe(double sigmaMax){Sigma=sigmaMax; Tau=sigmaMax; }
    FloatArray delivertomaterial();
    void initializeFrom(const double &NeedDissipated_G);
    double sign(double x );

    virtual const char *giveClassName() const { return "IntMatPPRINADM1Status"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    /// Functions for MaterialStatusMapperInterface
    virtual void copyStateVariables(const MaterialStatus &iStatus);
    virtual void addStateVariables(const MaterialStatus &iStatus);
};


/**
 * PPR cohesive zone model incorporating to XFEM.
 * Created on: July 17, 2017
 * @author Wencheng Jin
 */
class IntMatPPRINADM1 : public StructuralInterfaceMaterial
{
public:
    IntMatPPRINADM1(int n, Domain * d);
    virtual ~IntMatPPRINADM1();

protected:
    /// Material parameters
    double PenaltyStiffness;
    double GIc;    // fracture energy, mode I
    double GIIc;   // fracture energy, mode II
    double Sigma;  // max tensile stress
    double Tau;    // max shear
    double alpha;  //shape factors for tensile crack
    double beta;   //shape factors for shearing crack
    double lambdan; //normal penalty stiffness
    double lambdat; //tangent penalty stiffness
    double m;
    double n;

    double deltan;       // Final normal crack open widths
    double deltat;       // Final shear crack displacement
    double deltan_conj;  // Conjugated final normal crack widths which shear traction is zero
    double deltat_conj;  // Conjugated final shear crack widths which tensile traction is zero

    double Gamman; // Energy constant 1
    double Gammat; // Energy constant 2
    double dGtn;    // Energy constant 3
    double dGnt;    // Energy constant 4
    double TotalG;   //Total energy release rate to be dissipated

    virtual int checkConsistency();

public:

    virtual int hasNonLinearBehaviour()   { return 1; }

    virtual const char *giveClassName() const { return "IntMatPPRINADM1"; }
    virtual const char *giveInputRecordName() const { return _IFT_IntMatPPRINADM1_Name; }

    virtual void giveFirstPKTraction_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &jump,
                                        const FloatMatrix &F, TimeStep *tStep);

    // Dummy implementation, we must rely on numerical computation of the tangent.
    virtual void give3dStiffnessMatrix_dTdj(FloatMatrix &answer, MatResponseMode rMode, GaussPoint *gp, TimeStep *tStep);

    virtual bool hasAnalyticalTangentStiffness() const { return true; }

private:
    // Help functions
    //double computeYieldFunction(const double &iTractionNormal, const double &iTractionTang) const;
    //void computeTraction(FloatArray &oT, const FloatArray &iTTrial, const double &iPlastMultInc) const;

public:
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new IntMatPPRINADM1Status(1, domain, gp, alpha, beta, lambdan, lambdat, TotalG); }
    virtual void printYourself();

};
} /* namespace oofem */
#endif /* INTMatPPRINADM1_H_ */
