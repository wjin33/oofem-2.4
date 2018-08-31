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

#ifndef admnl2_h
#define admnl2_h

#include "Materials/adm2.h"
#include "Materials/structuralnonlocalmaterialext.h"
//#include "nonlocmatstiffinterface.h"

///@name Input fields for ADNLMaterial2
//@{
#define _IFT_ADNLMaterial2_Name "admnl2"
#define _IFT_ADNLMaterial2_r1 "r1"
#define _IFT_ADNLMaterial2_r2 "r2"
#define _IFT_ADNLMaterial2_averagingtype "averagingtype"
#define _IFT_ADNLMaterial2_exp "exp"
#define _IFT_ADNLMaterial2_rf "rf"
//@}

namespace oofem {
class GaussPoint;

/**
 * This class implements associated Material Status to ADNLMaterial2 (Nonlocal anisotropic damage for transverse isotropic material).
 * Stores local equivalent strains for averaging.
 */
class ADNLMaterial2Status : public AnisotropicDamageMaterial2Status, public StructuralNonlocalMaterialStatusExtensionInterface
{
protected:
    /// Equivalent strain for averaging.
    FloatArray localEquivalentStrainForAverage;

    /* // Variables used to track loading/reloading
     * public:
     * enum LastStateType {LST_elastic, LST_loading, LST_unloading};
     * LastStateType lst;
     */

public:
    /// Constructor.
    ADNLMaterial2Status(int n, Domain *d, GaussPoint *g);
    /// Destructor.
    virtual ~ADNLMaterial2Status();

    virtual void printOutputAt(FILE *file, TimeStep *tStep);

    /// Returns the local  equivalent strain to be averaged.
    FloatArray giveLocalEquivalentStrainForAverage() { return localEquivalentStrainForAverage; }
    /// Sets the localEquivalentStrainForAverage to given value.
    void setLocalEquivalentStrainForAverage(FloatArray ls) { localEquivalentStrainForAverage = ls; }

    // definition
    virtual const char *giveClassName() const { return "ADNLMaterial2Status"; }

    virtual void initTempStatus();
    virtual void updateYourself(TimeStep *tStep);

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
 * This class implements a Nonlocal Anisotropic Damage Model for transverse isotropic material (shale) in Tension and Compression
 * Model based on nonlocal averaging of equivalent strain.
 */
class ADNLMaterial2 : public AnisotropicDamageMaterial2, public StructuralNonlocalMaterialExtensionInterface
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
    ADNLMaterial2(int n, Domain *d);
    /// Destructor
    virtual ~ADNLMaterial2();

    // identification and auxiliary functions
    virtual const char *giveClassName() const { return "ADNLMaterial2"; }
    virtual const char *giveInputRecordName() const { return _IFT_ADNLMaterial2_Name; }

    virtual IRResultType initializeFrom(InputRecord *ir);

    virtual Interface *giveInterface(InterfaceType it);

    virtual void computeEquivalentStrain(FloatArray &EquivStrain, const FloatArray &strain,
                                         const double &TrSig, GaussPoint *gp, TimeStep *tStep);

    void computeLocalEquivalentStrain(FloatArray &kappa, const FloatArray &strain,
                                      const double &TrSig, GaussPoint *gp, TimeStep *tStep)
    { AnisotropicDamageMaterial2 :: computeEquivalentStrain(kappa, strain,TrSig, gp, tStep); }

    virtual void updateBeforeNonlocAverage(const FloatArray &strainVector, GaussPoint *gp, TimeStep *tStep);

    virtual int giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *tStep);

    virtual MaterialStatus *CreateStatus(GaussPoint *gp) const { return new ADNLMaterial2Status(1, AnisotropicDamageMaterial2::domain, gp); }

    virtual double predictRelativeComputationalCost(GaussPoint *gp);

};
} // end namespace oofem
#endif // admnl1_h
