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

#ifndef admnl1_h
#define admnl1_h

#include "Materials/adm1.h"
#include "Materials/structuralnonlocalmaterialext.h"
//#include "nonlocmatstiffinterface.h"

///@name Input fields for ADNLMaterial
//@{
#define _IFT_ADNLMaterial_Name "admnl1"
#define _IFT_ADNLMaterial_r "r"
#define _IFT_ADNLMaterial_averagingtype "averagingtype"
#define _IFT_ADNLMaterial_exp "exp"
#define _IFT_ADNLMaterial_rf "rf"
//@}

namespace oofem {
class GaussPoint;

/**
 * This class implements associated Material Status to ADNLMaterial (Nonlocal anisotropic damage).
 * Stores local equivalent strain for averaging.
 */
class ADNLMaterialStatus : public AnisotropicDamageMaterial1Status, public StructuralNonlocalMaterialStatusExtensionInterface
{
protected:
    /// Equivalent strain for averaging.
    double localEquivalentStrainForAverage;
    double nonlocalEquivalentStrain;

    /* // Variables used to track loading/reloading
     * public:
     * enum LastStateType {LST_elastic, LST_loading, LST_unloading};
     * LastStateType lst;
     */

public:
    /// Constructor.
    ADNLMaterialStatus(int n, Domain *d, GaussPoint *g);
    /// Destructor.
    virtual ~ADNLMaterialStatus();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    /// Returns the local  equivalent strain to be averaged.
    const double &giveLocalEquivalentStrainForAverage()const { return localEquivalentStrainForAverage; }
    /// Sets the localEquivalentStrainForAverage to given value.
    void setLocalEquivalentStrainForAverage(double ls) { localEquivalentStrainForAverage = ls; }
     /// Returns the nonlocal  equivalent strain.
    void setNonLocalEquivalentStrain(double ls) { nonlocalEquivalentStrain = ls; }
     /// Sets the nonlocalEquivalentStrain to given value.
    const double &giveNonLocalEquivalentStrain()const { return nonlocalEquivalentStrain; }

    // definition
    virtual const char *giveClassName() const { return "ADNLMaterialStatus"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

    /// Functions for MaterialStatusMapperInterface
    virtual void copyStateVariables(const MaterialStatus &iStatus);

    /**
     * Interface requesting service.
     * In the case of nonlocal constitutive models,
     * the use of multiple inheritance is assumed. Typically, the class representing nonlocal
     * constitutive model status is derived both from class representing local status and from class
     * NonlocalMaterialStatusExtensionInterface or from one of its derived classes
     * (which declare services and variables corresponding to specific analysis type).
     * @return In both cases, this function returns pointer to this object, obtained by
     * returning address of component or using pointer conversion from receiver to base class
     * NonlocalMaterialStatusExtensionInterface.
     */
    virtual Interface *giveInterface(InterfaceType);
};


/**
 * This class implements a Nonlocal Anisotropic Damage Model for Concrete in Tension and compression
 * Model based on nonlocal averaging of equivalent strain.
 */
class ADNLMaterial : public AnisotropicDamageMaterial1, public StructuralNonlocalMaterialExtensionInterface
{
    /**
protected:
    /// Final value of interaction radius, for a model with evolving characteristic length.
    double Rf;
    /// Parameter used as an exponent by models with evolving characteristic length.
    double exponent;
    /// Parameter specifying how the weight function should be adjusted due to damage.
    int averType;
     */

public:
    /// Constructor
    ADNLMaterial(int n, Domain *d);
    /// Destructor
    virtual ~ADNLMaterial();

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "ADNLMaterial"; }
    virtual const char *giveInputRecordName() const { return _IFT_ADNLMaterial_Name; }

    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual Interface *giveInterface(InterfaceType it);

    virtual void computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep);

    void computeLocalEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
    { AnisotropicDamageMaterial1 :: computeEquivalentStrain(kappa, strain, gp, tStep); }

    virtual void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new ADNLMaterialStatus(1, AnisotropicDamageMaterial1::domain, gp); }

    virtual double predictRelativeComputationalCost(GaussPoint *gp);

};
} // end namespace oofem
#endif // admnl1_h
