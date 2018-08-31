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

#ifndef adm2_h
#define adm2_h

// this turns on or off a bunch of internal variables
// that allow tracing the distribution of dissipated energy
// (can be turned off if such information is not needed)
#define keep_track_of_dissipated_energy

#include "material.h"
#include "../sm/Materials/structuralmaterial.h"
#include "isolinearelasticmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "randommaterialext.h"

///@name Input fields for Anisotropic Damage model with intrinsic anisotropic material
//@{
// @todo add input parametres which are read from input file here
#define _IFT_AnisotropicDamageMaterial2_Name "adm2"
#define _IFT_AnisotropicDamageMaterial2_exx "exx"
#define _IFT_AnisotropicDamageMaterial2_eyy "eyy"
#define _IFT_AnisotropicDamageMaterial2_nuxy "nuxy"
#define _IFT_AnisotropicDamageMaterial2_nuyz "nuyz"
#define _IFT_AnisotropicDamageMaterial2_gxy "gxy"
#define _IFT_AnisotropicDamageMaterial2_eqeps_1t "eqeps_1t"
#define _IFT_AnisotropicDamageMaterial2_eqeps_1s "eqeps_1s"
#define _IFT_AnisotropicDamageMaterial2_alpha_1t "alpha_1t"
#define _IFT_AnisotropicDamageMaterial2_eqeps_1c "eqeps_1c"
#define _IFT_AnisotropicDamageMaterial2_beta_1c "beta_1c"
#define _IFT_AnisotropicDamageMaterial2_alpha_1c "alpha_1c"
#define _IFT_AnisotropicDamageMaterial2_eqeps_2t "eqeps_2t"
#define _IFT_AnisotropicDamageMaterial2_eqeps_2s "eqeps_2s"
#define _IFT_AnisotropicDamageMaterial2_alpha_2t "alpha_2t"
#define _IFT_AnisotropicDamageMaterial2_eqeps_2c "eqeps_2c"
#define _IFT_AnisotropicDamageMaterial2_beta_2c "beta_2c"
#define _IFT_AnisotropicDamageMaterial2_alpha_2c "alpha_2c"
#define _IFT_AnisotropicDamageMaterial2_eta "eta"
#define _IFT_AnisotropicDamageMaterial2_theta "theta"
//@}


namespace oofem {
class GaussPoint;

/**
 * This class implements associated Material Status to AnisotropicDamageMaterial.
 * Stores a damage tensor and hardening variable (and possible extra information).
 */
class AnisotropicDamageMaterial2Status : public StructuralMaterialStatus,
        public RandomMaterialStatusExtensionInterface
{
protected:
    /// Scalar measure of the largest tensile strain level ever reached in material in directioin 1
    double equivalenstrain_1t;
    /// Non-equilibrated scalar measure of the largest tensile strain level.
    double tempequivalenstrain_1t;
    /// Scalar measure of the largest compressive strain level ever reached in materialin directioin 1
    double equivalenstrain_1c;
    /// Non-equilibrated scalar measure of the largest compressive strain level.
    double tempequivalenstrain_1c;
    /// Scalar measure of the largest tensile strain level ever reached in material in directioin 1
    double equivalenstrain_2t;
    /// Non-equilibrated scalar measure of the largest tensile strain level.
    double tempequivalenstrain_2t;
    /// Scalar measure of the largest compressive strain level ever reached in materialin directioin 1
    double equivalenstrain_2c;
    /// Non-equilibrated scalar measure of the largest compressive strain level.
    double tempequivalenstrain_2c;
    /// Scalar measure of the largest tensile strain level ever reached in material in directioin 1
    double equivalenstrain_3t;
    /// Non-equilibrated scalar measure of the largest tensile strain level.
    double tempequivalenstrain_3t;
    /// Scalar measure of the largest compressive strain level ever reached in materialin directioin 1
    double equivalenstrain_3c;
    /// Non-equilibrated scalar measure of the largest compressive strain level.
    double tempequivalenstrain_3c;
    /// Second order damage tensor
    /* damage components
     * damage11=omega_1t
     * damage12=omega_2t
     * damage13=omega_12s
     * damage21=omega_1c
     * damage22=omega_2c
     * damage23=omega_12s
     */
    FloatMatrix damage;
    /// Non-equilibrated second order damage tensor
    FloatMatrix tempDamage;
    double damaget1;
    double damaget2;
    double damagec1;
    double damagec2;
    /// This flag turns into 1 and remains 1 when the trace of the damage tensor is >1 in compression (tr(strainTensor)<0)
    int flag;
    int tempFlag;

    int flag_1st;

#ifdef keep_track_of_dissipated_energy
    /// Density of total work done by stresses on strain increments.
    double stressWork;
    /// Non-equilibrated density of total work done by stresses on strain increments.
    double tempStressWork;
    /// Density of dissipated work.
    double dissWork;
    /// Non-equilibrated density of dissipated work.
    double tempDissWork;
#endif

public:
    /// Constructor
    AnisotropicDamageMaterial2Status(int n, Domain *d, GaussPoint *g);
    /// Destructor
    virtual ~AnisotropicDamageMaterial2Status();

    /// Returns the last equilibrated scalar measure of the largest tensile strain level.
    double giveEquivalenstrain_1t() { return equivalenstrain_1t; }
    /// Returns the temp. scalar measure of the largest tensile strain level.
    double giveTempEquivalenstrain_1t() { return tempequivalenstrain_1t; }
    /// Returns the last equilibrated scalar measure of the largest tensile strain level.
    double giveEquivalenstrain_1c() { return equivalenstrain_1c; }
    /// Returns the temp. scalar measure of the largest tensile strain level.
    double giveTempEquivalenstrain_1c() { return tempequivalenstrain_1c; }
    /// Returns the last equilibrated scalar measure of the largest tensile strain level.
    double giveEquivalenstrain_2t() { return equivalenstrain_2t; }
    /// Returns the temp. scalar measure of the largest tensile strain level.
    double giveTempEquivalenstrain_2t() { return tempequivalenstrain_2t; }
    /// Returns the last equilibrated scalar measure of the largest tensile strain level.
    double giveEquivalenstrain_2c() { return equivalenstrain_2c; }
    /// Returns the temp. scalar measure of the largest tensile strain level.
    double giveTempEquivalenstrain_2c() { return tempequivalenstrain_2c; }
    /// Returns the last equilibrated scalar measure of the largest tensile strain level.
    double giveEquivalenstrain_3t() { return equivalenstrain_3t; }
    /// Returns the temp. scalar measure of the largest tensile strain level.
    double giveTempEquivalenstrain_3t() { return tempequivalenstrain_3t; }
    /// Returns the last equilibrated scalar measure of the largest tensile strain level.
    double giveEquivalenstrain_3c() { return equivalenstrain_3c; }
    /// Returns the temp. scalar measure of the largest tensile strain level.
    double giveTempEquivalenstrain_3c() { return tempequivalenstrain_3c; }
    /// Sets the temp scalar measure of the largest tensile strain level to given value.
    void setTempEquivalenstrain_1t(double newEquivalenstrain_1t) { tempequivalenstrain_1t = newEquivalenstrain_1t; }
    /// Sets the temp scalar measure of the largest compressive strain level to given value.
    void setTempEquivalenstrain_1c(double newEquivalenstrain_1c) { tempequivalenstrain_1c = newEquivalenstrain_1c; }
    /// Sets the temp scalar measure of the largest tensile strain level to given value.
    void setTempEquivalenstrain_2t(double newEquivalenstrain_2t) { tempequivalenstrain_2t = newEquivalenstrain_2t; }
    /// Sets the temp scalar measure of the largest compressive strain level to given value.
    void setTempEquivalenstrain_2c(double newEquivalenstrain_2c) { tempequivalenstrain_2c = newEquivalenstrain_2c; }
    /// Sets the temp scalar measure of the largest tensile strain level to given value.
    void setTempEquivalenstrain_3t(double newEquivalenstrain_3t) { tempequivalenstrain_3t = newEquivalenstrain_3t; }
    /// Sets the temp scalar measure of the largest compressive strain level to given value.
    void setTempEquivalenstrain_3c(double newEquivalenstrain_3c) { tempequivalenstrain_3c = newEquivalenstrain_3c; }
    /// Returns the last equilibrated second order damage tensor.
    const FloatMatrix &giveDamage() { return damage; }
    /// Returns the temp. second order damage tensor.
    const FloatMatrix &giveTempDamage() { return tempDamage; }
    /// Assigns temp. damage tensor to given tensor d
    void setTempDamage(const FloatMatrix &d) { tempDamage = d; }

    void setdamaget1(double newdamaget1) { damaget1 = newdamaget1; }
    void setdamaget2(double newdamaget2) { damaget2 = newdamaget2; }
    void setdamagec1(double newdamagec1) { damagec1 = newdamagec1; }
    void setdamagec2(double newdamagec2) { damagec2 = newdamagec2; }

    double givedamaget1(){return damaget1; }
    double givedamaget2(){return damaget2; }
    double givedamagec1(){return damagec1; }
    double givedamagec2(){return damagec2; }

    /// Returns the value of the flag.
    int giveFlag() { return flag; }
    /// Sets the value of the temporary value of flag.
    void setTempFlag(int newflag) { tempFlag = newflag; }
    /// Returns the value of the temporary value of flag.
    int giveTempFlag() { return tempFlag; }

    int giveFlag_1st() { return flag_1st; }
    void setFlag_1st(int newflag) { flag_1st = newflag; }

    /*
    /// Returns the value of the flag.
    int giveFlag2() { return flag2; }
    /// Sets the value of the temporary value of flag.
    void setTempFlag2(int newflag) { tempFlag2 = newflag; }
    /// Returns the value of the temporary value of flag.
    int giveTempFlag2() { return tempFlag2; }
    */

#ifdef keep_track_of_dissipated_energy
    /// Returns the density of total work of stress on strain increments.
    double giveStressWork() { return stressWork; }
    /// Returns the temp density of total work of stress on strain increments.
    double giveTempStressWork() { return tempStressWork; }
    /// Sets the density of total work of stress on strain increments to given value.
    void setTempStressWork(double w) { tempStressWork = w; }
    /// Returns the density of dissipated work.
    double giveDissWork() { return dissWork; }
    /// Returns the density of temp dissipated work.
    double giveTempDissWork() { return tempDissWork; }
    /// Sets the density of dissipated work to given value.
    void setTempDissWork(double w) { tempDissWork = w; }
    /// Computes the increment of total stress work and of dissipated work.
    void computeWork(GaussPoint *gp);
#endif

    // definition
    virtual const char *giveClassName() const { return "AnisotropicDamageMaterial2Status"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    Interface * giveInterface(InterfaceType type);

};


/**
 * Class representing anisotropic damage model for intrinsic transverse isotropic material
 * @author : Wencheng Jin : wencheng.jin@gatech.edu
 */
class AnisotropicDamageMaterial2 : public StructuralMaterial,
        public RandomMaterialExtensionInterface
{
protected:

    /// Young's modulus
    double E11, E22;
    /// Poisson's ratio
    double nu12, nu23;
    /// Shear modulus
    double G12;
    /// Tensile equivalent strain threshold in axial direction
    double eqeps_1t;
    /// Shear equivalent strain threshold in 12 direction
    double eqeps_1s;
    /// Hardening/softening parameter for tension
    double alpha_1t;
    /// Compressive equivalent strain threshold in axial direction
    double eqeps_1c;
    /// Hardening/softening parameters for compression
    double beta_1c, alpha_1c;
    /// Tensile equivalent strain threshold in transverse direction
    double eqeps_2t;
    /// Shear equivalent strain threshold in 12 direction
    double eqeps_2s;
    /// Hardening/softening parameter for tension
    double alpha_2t;
    /// Compressive equivalent strain threshold in transverse direction
    double eqeps_2c;
    /// Hardening/softening parameters for compression
    double beta_2c, alpha_2c;
    /// confinement dependent parameter
    double eta;
    /// local Angle
    double theta;

public:
    /// Constructor
    AnisotropicDamageMaterial2(int n, Domain *d);
    /// Destructor
    virtual ~AnisotropicDamageMaterial2();

    virtual int hasMaterialModeCapability(MaterialMode mode);

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "AnisotropicDamageMaterial2"; }
    virtual const char *giveInputRecordName() const { return _IFT_AnisotropicDamageMaterial2_Name; }

    //Initialization
    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const ;
    virtual MaterialStatus *giveStatus(GaussPoint *gp) const;

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp,
                            InternalStateType type,TimeStep *atTime);

    virtual void computeEquivalentStrain(FloatArray &EquivStrain, const FloatArray &strain,
                                         const double &TrSig, GaussPoint *gp, TimeStep *tStep);

    void giveRealStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp,
                                          const FloatArray &totalStrain, TimeStep *atTime);

    void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep);

    int giveDamageState(GaussPoint *gp);
/////---------------------------------------------
/*
 *   /// Plain-strain version of stress evaluation algorithm-by Wencheng Jin
 *  void Aijkl_Bkl(const double A[3][3][3][3], const FloatMatrix &B, FloatMatrix &C);
    void computePrincValDir2D(double &D1, double &D2, double &c, double &s, double Dx, double Dy, double Dxy);
    void computeSecantOperator(FloatMatrix &answer, double MATD[3][3][3][3],  const double &Dc, const FloatMatrix &damageTensor);
    //void computeSecantOperator(FloatMatrix &answer, FloatMatrix & MATD,  double Dc, FloatMatrix damageTensor, GaussPoint *gp );
    void MAT2_MAT4(const FloatMatrix &DMATRIX, double TENSOR[3][3][3][3], int ICOE);
    void MAT4_MAT2(const double TENSOR[3][3][3][3], FloatMatrix &DMATRIX, int ICOE);
*/
/////---------------------------------------------

protected:

    virtual void givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseMode mmode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,MatResponseMode mmode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);
};
} // end namespace oofem
#endif // adm2_h
