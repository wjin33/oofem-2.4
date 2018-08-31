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

#ifndef adm1_h
#define adm1_h

// this turns on or off a bunch of internal variables
// that allow tracing the distribution of dissipated energy
// (can be turned off if such information is not needed)
#define keep_track_of_dissipated_energy

#include "material.h"
#include "../sm/Materials/structuralmaterial.h"
#include "isolinearelasticmaterial.h"
#include "../sm/Materials/structuralms.h"
#include "randommaterialext.h"

///@name Input fields for AnisotropicDamageMaterial
//@{
// @todo add input parametres which are read from input file here
#define _IFT_AnisotropicDamageMaterial1_Name "adm1"
#define _IFT_AnisotropicDamageMaterial1_kappat "kappa_t"
#define _IFT_AnisotropicDamageMaterial1_alphat "alpha_t"
#define _IFT_AnisotropicDamageMaterial1_kappac "kappa_c"
#define _IFT_AnisotropicDamageMaterial1_alphac "alpha_c"
#define _IFT_AnisotropicDamageMaterial1_eta "eta"
//@}

namespace oofem {
class GaussPoint;

/**
 * This class implements associated Material Status to AnisotropicDamageMaterial.
 * Stores a damage tensor and hardening variable (and possible extra information).
 */
class AnisotropicDamageMaterial1Status : public StructuralMaterialStatus,
        public RandomMaterialStatusExtensionInterface
{
protected:
    /// Scalar measure of the largest Mazars' tensile strain level ever reached in material.
    double equivalenstrain_t;
    /// Non-equilibrated scalar measure of the largest tensile strain level.
    double tempequivalenstrain_t;
    /// Scalar measure of the largest compressive strain level ever reached in material.
    double equivalenstrain_c;
    /// Non-equilibrated scalar measure of the largest compressive strain level.
    double tempequivalenstrain_c;
    /// Second order damage tensor
    FloatMatrix damage;
    /// Non-equilibrated second order damage tensor
    FloatMatrix tempDamage;
    /// This flag turns into 1 and remains 1 when the trace of the damage tensor is >1 in compression (tr(strainTensor)<0)
    int flag;
    int tempFlag;

    //Bonus parameter
    FloatMatrix permeability;


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
    AnisotropicDamageMaterial1Status(int n, Domain *d, GaussPoint *g);
    /// Destructor
    virtual ~AnisotropicDamageMaterial1Status();

    /// Returns the last equilibrated scalar measure of the largest tensile strain level.
    const double &giveEquivalenstrain_t() const { return equivalenstrain_t; }
    /// Returns the temp. scalar measure of the largest tensile strain level.
    const double &giveTempEquivalenstrain_t() const { return tempequivalenstrain_t; }

    void setEquivalenstrain_t(double newEquivalenstrain_t) { equivalenstrain_t = newEquivalenstrain_t; }
    void setEquivalenstrain_c(double newEquivalenstrain_c) { equivalenstrain_c = newEquivalenstrain_c; }

    /// Returns the last equilibrated scalar measure of the largest tensile strain level.
    const double &giveEquivalenstrain_c() const { return equivalenstrain_c; }
    /// Returns the temp. scalar measure of the largest tensile strain level.
    const double &giveTempEquivalenstrain_c() const { return tempequivalenstrain_c; }
    /// Sets the temp scalar measure of the largest tensile strain level to given value.
    void setTempEquivalenstrain_t(double newEquivalenstrain_t) { tempequivalenstrain_t = newEquivalenstrain_t; }
    /// Sets the temp scalar measure of the largest compressive strain level to given value.
    void setTempEquivalenstrain_c(double newEquivalenstrain_c) { tempequivalenstrain_c = newEquivalenstrain_c; }
    /// Returns the last equilibrated second order damage tensor.
    const FloatMatrix &giveDamage() const { return damage; }
    /// Returns the temp. second order damage tensor.
    const FloatMatrix &giveTempDamage() const { return tempDamage; }
    /// Assigns temp. damage tensor to given tensor d
    void setTempDamage(const FloatMatrix &d) { tempDamage = d; }
    void setDamage(const FloatMatrix &d) { damage = d; }

    const FloatMatrix &givePermeability() const { return permeability; }

    /// Returns the value of the flag.
    const int &giveFlag() const { return flag; }
    /// Sets the value of the temporary value of flag.
    void setTempFlag(int newflag) { tempFlag = newflag; }
    /// Returns the value of the temporary value of flag.
    const int &giveTempFlag() const { return tempFlag; }

#ifdef keep_track_of_dissipated_energy
    /// Returns the density of total work of stress on strain increments.
    const double &giveStressWork() const { return stressWork; }
    /// Returns the temp density of total work of stress on strain increments.
    const double &giveTempStressWork() const { return tempStressWork; }
    /// Sets the density of total work of stress on strain increments to given value.
    void setStressWork(double w) { stressWork = w; }
    void setTempStressWork(double w) { tempStressWork = w; }
    /// Returns the density of dissipated work.
    const double &giveDissWork() const { return dissWork; }
    /// Returns the density of temp dissipated work.
    const double &giveTempDissWork() const { return tempDissWork; }
    /// Sets the density of dissipated work to given value.
    void setDissWork(double w) { tempDissWork = w; }
    void setTempDissWork(double w) { tempDissWork = w; }
    /// Computes the increment of total stress work and of dissipated work.
    void computeWork(GaussPoint *gp);
    void giveInternalStatetype(IntArray &G);
#endif

    // definition
    virtual const char *giveClassName() const { return "AnisotropicDamageMaterial1Status"; }

    void computePrincValDir2D(double &D1, double &D2, double &c, double &s, double Dx, double Dy, double Dxy);

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);
    virtual void printOutputAt(FILE *file, TimeStep *tStep);
    Interface * giveInterface(InterfaceType type);

    /// Functions for MaterialStatusMapperInterface
    virtual void copyStateVariables(const MaterialStatus &iStatus);
    //virtual void addStateVariables(const MaterialStatus &iStatus);


};


/**
 * Class representing anisotropic damage model.
 * @author : Wencheng Jin : wencheng.jin@gatech.edu
 */
class AnisotropicDamageMaterial1 : public IsotropicLinearElasticMaterial,
        public RandomMaterialExtensionInterface
{
protected:

    /// Reference to bulk (undamaged) material
    //IsotropicLinearElasticMaterial *linearElasticMaterial;
    /// Young's modulus
    //double E;
    /// Poisson's ratio
    //double nu;
    /// Tensile damage threshold kappa_t
    double kappa_t;
    /// Compressive damage threshold kappa_c
    double kappa_c;
    /// Tensile hardening parameter
    double alpha_t;
    /// Compressive hardening parameter
    double alpha_c;
    /// confining influence parameter
    double eta;

public:
    /// Constructor
    AnisotropicDamageMaterial1(int n, Domain *d);
    /// Destructor
    virtual ~AnisotropicDamageMaterial1();

    virtual int hasNonLinearBehaviour() { return 1; }
    virtual int hasMaterialModeCapability(MaterialMode mode);

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "AnisotropicDamageMaterial1"; }
    virtual const char *giveInputRecordName() const { return _IFT_AnisotropicDamageMaterial1_Name; }

    // Returns reference to undamaged (bulk) material
/*
    /// Plane-stress version of the stress evaluation algorithm
    virtual void giveRealStressVector_PlaneStress(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep);
//    void computePrincValDir2D(double &D1, double &D2, double &c, double &s, double Dx, double Dy, double Dxy);
    bool checkPrincVal2D(double Dx, double Dy, double Dxy);
    void computeDamage(FloatMatrix &tempDamage, const FloatMatrix &damage, double kappa, double eps1, double eps2, double ceps, double seps, double epsZ);
    //double computeTraceD(double equivStrain);
    double computeOutOfPlaneStrain(const FloatArray &inplaneStrain, const FloatMatrix &dam, bool tens_flag);
    double computeDimensionlessOutOfPlaneStress(const FloatArray &inplaneStrain, double epsZ, const FloatMatrix &dam);
    void computeInplaneStress(FloatArray &inplaneStress, const FloatArray &inplaneStrain, double epsZ, const FloatMatrix &dam);
*/
/////---------------------------------------------
    /// Plain-strain version of stress evaluation algorithm-by Wencheng Jin
    void giveRealStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &totalStrain, TimeStep *atTime);
    void Aijkl_Bkl(const double A[3][3][3][3], const FloatMatrix &B, FloatMatrix &C);
	void computePrincValDir2D(double &D1, double &D2, double &c, double &s, double Dx, double Dy, double Dxy);
    void computeSecantOperator(FloatMatrix &answer, double MATD[3][3][3][3],  const double &Dc, const FloatMatrix &damageTensor);
    //void computeSecantOperator(FloatMatrix &answer, FloatMatrix & MATD,  double Dc, FloatMatrix damageTensor, GaussPoint *gp );
    void MAT2_MAT4(const FloatMatrix &DMATRIX, double TENSOR[3][3][3][3], int ICOE);
    void MAT4_MAT2(const double TENSOR[3][3][3][3], FloatMatrix &DMATRIX, int ICOE);
//	void givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseMode mode, GaussPoint *gp, TimeStep *atTime);
    //initialize the form
    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const ;
    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual MaterialStatus * giveStatus(GaussPoint *gp) const;
    // Used to check error, for errorcheck exportmodule use
    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime);
    virtual void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);
    //Interface * giveInterface(InterfaceType type);
    int giveDamageState(GaussPoint *gp);
/////---------------------------------------------
/*
    void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep);
    void computePrincValDir3D(double &D1, double &D2, double &c, double &s, double Dx, double Dy, double Dxy);
    void jacobi_eigenvalue ( int n, double a[], int it_max, double v[], double d[], int &it_num, int &rot_num );
    void r8mat_diag_get_vector ( int n, double a[], double v[] );
    void r8mat_identity ( int n, double a[] );
    */
/*
    /// Obtains the proportion of the damage tensor that is needed to get the first eigenvalue equal to the damage threshold
    double obtainAlpha1(FloatMatrix tempDamageTensor, double deltaLambda, FloatMatrix positiveStrainTensor, double damageThreshold);
    /// Obtains the proportion of the damage tensor that is needed to get the second eigenvalue equal to the damage threshold
    double obtainAlpha2(FloatMatrix tempDamageTensor, double deltaLambda, FloatMatrix positiveStrainTensor, FloatMatrix ProjMatrix, double damageThreshold);
    /// Obtains the proportion of the damage tensor that is needed to get the third eigenvalue equal to the damage threshold
    double obtainAlpha3(FloatMatrix tempDamageTensor, double deltaLambda, FloatMatrix positiveStrainTensor, FloatArray vec3, double damageThreshold);

    double checkSymmetry(FloatMatrix matrix);

    void correctBigValues(FloatMatrix &matrix);

    double computeTraceD(FloatMatrix tempDamageTensor, FloatMatrix strainTensor, GaussPoint *gp);

    double computeCorrectionFactor(FloatMatrix tempDamageTensor, FloatMatrix strainTensor, GaussPoint *gp);

    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,
                                               MatResponseMode mode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);

    virtual void giveRealStressVector(FloatArray &answer,  GaussPoint *gp,
                                      const FloatArray &reducedStrain, TimeStep *tStep);

    virtual void giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_StressControl(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, const IntArray &strainControl, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }
    virtual void giveRealStressVector_1d(FloatArray &answer, GaussPoint *gp, const FloatArray &reducedE, TimeStep *tStep)
    { this->giveRealStressVector(answer, gp, reducedE, tStep); }

    

    //    virtual InternalStateValueType giveIPValueType(InternalStateType type);
    //virtual void giveThermalDilatationVector(FloatArray &answer, GaussPoint *, TimeStep *);


     * Computes the equivalent strain measure from given strain vector (full form).
     * @param[out] kappa Return parameter, containing the corresponding equivalent strain.
     * @param strain Total strain vector in full form.
     * @param gp Integration point.
     * @param tStep Time step.

    virtual void computeDamageTensor(FloatMatrix &answer, GaussPoint *gp,
                                     const FloatArray &totalStrain, double equivStrain,
                                     TimeStep *atTime);
    */
    virtual void giveInputRecord(DynamicInputRecord &input);

protected:


    virtual void givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseMode mmode,
                                          GaussPoint *gp,
                                         TimeStep *tStep);
/*
    virtual void give3dMaterialStiffnessMatrix(FloatMatrix &answer,MatResponseMode mmode,
                                               GaussPoint *gp,
                                               TimeStep *tStep);

    virtual void givePlaneStressStiffMtrx(FloatMatrix &answer, MatResponseMode mmode,
                                          GaussPoint *gp,
                                          TimeStep *tStep);

    virtual void give1dStressStiffMtrx(FloatMatrix &answer, MatResponseMode mmode,
                                       GaussPoint *gp,
                                       TimeStep *tStep);
    virtual void computePlaneStressStrain(FloatMatrix &answer, FloatMatrix damageTensor, FloatArray totalStrain, GaussPoint *gp,
                                          TimeStep *atTime);
/    virtual void computePlaneStressSigmaZ(double &answer, FloatMatrix damageTensor, FloatArray reducedTotalStrainVector,
                                          double epsZ, GaussPoint *gp, TimeStep *atTime);

*/

};
} // end namespace oofem
#endif // adm1_h
