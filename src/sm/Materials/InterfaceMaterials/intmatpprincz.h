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

#ifndef INTMatPPRINCZ_H_
#define INTMatPPRINCZ_H_

#include "structuralinterfacematerial.h"
#include "structuralinterfacematerialstatus.h"

///@name Input fields for IntMatPPRINCZ
//@{
#define _IFT_IntMatPPRINCZ_Name "intmatpprincz"
//#define _IFT_IntMatPPRINCZ_PenaltyStiffness "kn"
#define _IFT_IntMatPPRINCZ_phi1  "phi1"
#define _IFT_IntMatPPRINCZ_phi2  "phi2"
#define _IFT_IntMatPPRINCZ_sigm  "sigm"
#define _IFT_IntMatPPRINCZ_taum  "taum"
#define _IFT_IntMatPPRINCZ_alpha "alpha"
#define _IFT_IntMatPPRINCZ_beta  "beta"
#define _IFT_IntMatPPRINCZ_lambdan "lambdan"
#define _IFT_IntMatPPRINCZ_lambdat  "lambdat"
//@}

namespace oofem {
/**
 * This class implements associated Material Status for IntMatPPRINCZ
 */
class IntMatPPRINCZStatus : public StructuralInterfaceMaterialStatus
{
public:
    IntMatPPRINCZStatus(int n, Domain * d, GaussPoint * g);
    virtual ~IntMatPPRINCZStatus();

    /// history variable
    /// Vector measure of the largest displacement jump ever reached in material.
    FloatArray JumpMax;
    /// Non-equilibrated vector measure of the largest displacement jump ever reached in material.
    FloatArray tempJumpMax;

    /// Returns the last equilibrated vector measure largest displacement jump.
    const FloatArray &giveJumpMax() { return JumpMax; }
    /// Returns the temp. vector measure largest displacement jump.
    const FloatArray &giveTempJumpMax() { return tempJumpMax; }
    /// Sets the vector measure of the displacement jump to given value.
    void letJumpMaxBe(FloatArray newJumpMax) { JumpMax = std :: move(newJumpMax); }
    /// Sets the temp vector measure of the displacement jump to given value.
    void letTempJumpMaxBe(FloatArray newTempJumpMax) { tempJumpMax = std :: move(newTempJumpMax); }


    virtual const char *giveClassName() const { return "IntMatPPRINCZStatus"; }

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
class IntMatPPRINCZ : public StructuralInterfaceMaterial
{
public:
    IntMatPPRINCZ(int n, Domain * d);
    virtual ~IntMatPPRINCZ();

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

    virtual int checkConsistency();

public:

    virtual int hasNonLinearBehaviour()   { return 1; }

    virtual const char *giveClassName() const { return "IntMatPPRINCZ"; }
    virtual const char *giveInputRecordName() const { return _IFT_IntMatPPRINCZ_Name; }

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

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new IntMatPPRINCZStatus(1, domain, gp); }
    virtual void printYourself();
    double sign(double x );
};
} /* namespace oofem */
#endif /* INTMatPPRINCZ_H_ */
