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

#ifndef PLNonlocalDamage_H_
#define PLNonlocalDamage_H_

#include "xfem/propagationlaw.h"

#define _IFT_PLNonlocalDamage_Name "propagationlawnonlocaldamage"
#define _IFT_PLNonlocalDamage_AngleInc "angleinc" ///< Angle between sampling points on the circle
#define _IFT_PLNonlocalDamage_IncLength "incrementlength" ///< Increment length per time step
#define _IFT_PLNonlocalDamage_DamageThreshold "damagethreshold" ///< Threshold for macro fracture initiation
#define _IFT_PLNonlocalDamage_RadialBasisFunc "useradialbasisfunc" ///< If radial basis functions should be used for stress interpolation

namespace oofem {
class Domain;
class EnrichmentDomain;
class DynamicInputRecord;


/**
 * Propagation law that propagates the crack in the direction
 * that gives $ \damage_{r\theta} = 0 && \damage_{\theta|\theta} >= Threshold$.
 * Based on
 * T.P. Fries and M. Baydoun:
 * "Crack propagation with the extended finite element method
 * and a hybrid explicit-implicit crack description",
 * Internat. J. Numer. Methods Engrg 89,
 * pp. 1527--1558 (2012)
 *
 * The damage is evaluated in several points on a circle
 * surrounding the crack tip.
 *
 * Compared to the paper above, the implementation has been extended
 * with a criterion for crack propagation instead of always
 * propagating a predefined increment length. Two options are
 * currently available for stress interpolation:
 * 1) Take stress of closest Gauss point
 * 2) Interpolate with radial basis functions
 *
 * @author Wencheng Jin & Erik Svenning
 */
class OOFEM_EXPORT PLNonlocalDamage : public PropagationLaw
{
public:
    PLNonlocalDamage() :  mAngleInc(0.0), mIncrementLength(0.0), mDamageThreshold(0.0), mUseRadialBasisFunc(false) { }
    virtual ~PLNonlocalDamage() { }

    virtual const char *giveClassName() const { return "PLNonlocalDamage"; }
    virtual const char *giveInputRecordName() const { return _IFT_PLNonlocalDamage_Name; }

    virtual IRResultType initializeFrom(InputRecord *ir);
    virtual void giveInputRecord(DynamicInputRecord &input);

    virtual bool hasPropagation() const { return mIncrementLength > 0.; } ///
    virtual bool propagateInterface(Domain &iDomain, EnrichmentFront &iEnrFront, TipPropagation &oTipProp);

protected:
    double  mAngleInc, mIncrementLength, mDamageThreshold;
    bool mUseRadialBasisFunc;
};
} // end namespace oofem


#endif /* PLNonlocalDamage_H_ */
