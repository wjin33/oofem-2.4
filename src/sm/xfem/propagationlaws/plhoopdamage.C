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

#include "xfem/propagationlaws/plhoopdamage.h"

#include "xfem/propagationlaw.h"
#include "xfem/tipinfo.h"
#include "classfactory.h"
#include "mathfem.h"
#include "dynamicinputrecord.h"
#include "spatiallocalizer.h"
#include "floatmatrix.h"
#include "gausspoint.h"
#include "../sm/Materials/structuralms.h"
#include "../sm/Materials/adm1.h"
#include "xfem/enrichmentitem.h"
#include "feinterpol.h"
#include "xfem/xfemmanager.h"

#include "xfem/XFEMDebugTools.h"

namespace oofem {
REGISTER_PropagationLaw(PLHoopDamage)

/////////////////////////////////////////////
IRResultType PLHoopDamage :: initializeFrom(InputRecord *ir)
{
    IRResultType result;

    IR_GIVE_FIELD(ir, mRadius,                          _IFT_PLHoopDamage_Radius);
    IR_GIVE_FIELD(ir, mAngleInc,                        _IFT_PLHoopDamage_AngleInc);
    IR_GIVE_FIELD(ir, mIncrementLength,         _IFT_PLHoopDamage_IncLength);
    IR_GIVE_FIELD(ir, mDamageThreshold, _IFT_PLHoopDamage_DamageThreshold);

    int useRadialBasisFunc = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, useRadialBasisFunc, _IFT_PLHoopDamage_RadialBasisFunc);
    if ( useRadialBasisFunc == 1 ) {
        mUseRadialBasisFunc = true;
    }

    return IRRT_OK;
}

void PLHoopDamage :: giveInputRecord(DynamicInputRecord &input)
{
    int number = 1;
    input.setRecordKeywordField(this->giveInputRecordName(), number);

    input.setField(mRadius,                             _IFT_PLHoopDamage_Radius);
    input.setField(mAngleInc,                           _IFT_PLHoopDamage_AngleInc);
    input.setField(mIncrementLength,            _IFT_PLHoopDamage_IncLength);
    input.setField(mDamageThreshold,        _IFT_PLHoopDamage_DamageThreshold);

    if ( mUseRadialBasisFunc ) {
        input.setField(1,       _IFT_PLHoopDamage_RadialBasisFunc);
    }
}

bool PLHoopDamage :: propagateInterface(Domain &iDomain, EnrichmentFront &iEnrFront, TipPropagation &oTipProp)
{
    if ( !iEnrFront.propagationIsAllowed() ) {
        return false;
    }

    // Fetch crack tip data
    const TipInfo &tipInfo = iEnrFront.giveTipInfo();

    SpatialLocalizer *localizer = iDomain.giveSpatialLocalizer();

    // Construct circle points on an arc from -90 to 90 degrees
    double angle = -90.0 + mAngleInc;
    std :: vector< double >angles;
    while ( angle <= ( 90.0 - mAngleInc ) ) {
        angles.push_back(angle * M_PI / 180.0);
        angle += mAngleInc;
    }

    const FloatArray &xT    = tipInfo.mGlobalCoord;
    const FloatArray &t     = tipInfo.mTangDir;
    const FloatArray &n     = tipInfo.mNormalDir;

    // Transfrom damage to fracture
    //FloatArray Tip  = xT;
    //Tip.at(1)+=1;
    Element *el = localizer->giveElementContainingPoint(tipInfo.mGlobalCoord);
    //Element *el = localizer->giveElementContainingPoint(tipInfo.mGlobalCoord);
    if ( el != NULL ) {
        std :: vector< FloatArray >circPoints;

        for ( size_t i = 0; i < angles.size(); i++ ) {
            FloatArray tangent(2);
            tangent.zero();
            tangent.add(cos(angles [ i ]), t);
            tangent.add(sin(angles [ i ]), n);
            tangent.normalize();

            FloatArray x(xT);
            x.add(mRadius, tangent);
            circPoints.push_back(x);
        }



        std :: vector< double >omegaTTArray, omegaRTArray;

        // Loop over circle points
        for ( size_t pointIndex = 0; pointIndex < circPoints.size(); pointIndex++ ) {
            FloatMatrix damageTen;

            if ( mUseRadialBasisFunc ) {
                // Interpolate stress with radial basis functions

                // Choose a cut-off length l:
                // take the distance between two nodes in the element containing the
                // crack tip multiplied by a constant factor.
                // ( This choice implies that we hope that the element has reasonable
                // aspect ratio.)
                const FloatArray &x1 = * ( el->giveDofManager(1)->giveCoordinates() );
                const FloatArray &x2 = * ( el->giveDofManager(2)->giveCoordinates() );
                const double l = 1.0 * x1.distance(x2);

                // Use the octree to get all elements that have
                // at least one Gauss point in a certain region around the tip.
                const double searchRadius = 3.0 * l;
                std :: set< int >elIndices;
                localizer->giveAllElementsWithIpWithinBox(elIndices, circPoints [ pointIndex ], searchRadius);


                // Loop over the elements and Gauss points obtained.
                // Evaluate the interpolation.
                FloatMatrix sumQiWiVi;
                double sumWiVi = 0.0;
                for ( int elIndex: elIndices ) {
                    Element *gpEl = iDomain.giveElement(elIndex);

                    for ( GaussPoint *gp_i: *gpEl->giveDefaultIntegrationRulePtr() ) {
                        ////////////////////////////////////////
                        // Compute global gp coordinates
                        FloatArray N;
                        FEInterpolation *interp = gpEl->giveInterpolation();
                        interp->evalN( N, gp_i->giveNaturalCoordinates(), FEIElementGeometryWrapper(gpEl) );


                        // Compute global coordinates of Gauss point
                        FloatArray globalCoord(2);
                        globalCoord.zero();

                        for ( int i = 1; i <= gpEl->giveNumberOfDofManagers(); i++ ) {
                            DofManager *dMan = gpEl->giveDofManager(i);
                            globalCoord.at(1) += N.at(i) * dMan->giveCoordinate(1);
                            globalCoord.at(2) += N.at(i) * dMan->giveCoordinate(2);
                        }


                        ////////////////////////////////////////
                        // Compute weight of kernel function

                        FloatArray tipToGP;
                        tipToGP.beDifferenceOf(globalCoord, xT);
                        bool inFrontOfCrack = true;
                        if ( tipToGP.dotProduct(t) < 0.0 ) {
                            inFrontOfCrack = false;
                        }

                        double r = circPoints [ pointIndex ].distance(globalCoord);

                        if ( r < l && inFrontOfCrack ) {
                            double w = ( ( l - r ) / ( pow(2.0 * M_PI, 1.5) * pow(l, 3) ) ) * exp( -0.5 * pow(r, 2) / pow(l, 2) );

                            // Compute gp volume
                            double V = gpEl->computeVolumeAround(gp_i);

                            // Get stress
                            AnisotropicDamageMaterial1Status *ms = dynamic_cast< AnisotropicDamageMaterial1Status * >( gp_i->giveMaterialStatus() );
                            if ( ms == NULL ) {
                                OOFEM_ERROR("failed to fetch MaterialStatus.");
                            }

                            FloatMatrix damageTenGP = ms->giveDamage();


                            if ( sumQiWiVi.giveNumberOfRows() != damageTenGP.giveNumberOfRows() || sumQiWiVi.giveNumberOfColumns() != damageTenGP.giveNumberOfColumns()  ) {
                                sumQiWiVi.resize( damageTenGP.giveNumberOfRows(),damageTenGP.giveNumberOfColumns()  );
                                sumQiWiVi.zero();
                            }

                            // Add to numerator
                            sumQiWiVi.add(w * V, damageTenGP);

                            // Add to denominator
                            sumWiVi += w * V;
                        }
                    }
                }


                if ( fabs(sumWiVi) > 1.0e-12 ) {
                    damageTen.add(1.0 / sumWiVi, sumQiWiVi);
                } else {
                    // Take stress from closest Gauss point
                    int region = 1;
                    bool useCZGP = false;
                    GaussPoint &gp = * ( localizer->giveClosestIP(circPoints [ pointIndex ], region, useCZGP) );


                    // Compute stresses
                    AnisotropicDamageMaterial1Status *ms = dynamic_cast< AnisotropicDamageMaterial1Status * >( gp.giveMaterialStatus() );
                    if ( ms == NULL ) {
                        OOFEM_ERROR("failed to fetch MaterialStatus.");
                    }

                    damageTen = ms->giveDamage();
                }
            } else {
                // Take stress from closest Gauss point
                int region = 1;
                bool useCZGP = false;
                GaussPoint &gp = * ( localizer->giveClosestIP(circPoints [ pointIndex ], region, useCZGP) );


                // Compute stresses
                AnisotropicDamageMaterial1Status *ms = dynamic_cast< AnisotropicDamageMaterial1Status * >( gp.giveMaterialStatus() );
                if ( ms == NULL ) {
                    OOFEM_ERROR("failed to fetch MaterialStatus.");
                }

                damageTen = ms->giveDamage();
            }

            //FloatMatrix stress(2, 2);

            //int shearPos = stressVec.giveSize();

            //stress.at(1, 1) = stressVec.at(1);
            // stress.at(1, 2) = stressVec.at(shearPos);
            //stress.at(2, 1) = stressVec.at(shearPos);
            //stress.at(2, 2) = stressVec.at(2);


            // Rotation matrix
            FloatMatrix rot(2, 2);
            rot.at(1, 1) =  cos(angles [ pointIndex ]);
            rot.at(1, 2) = -sin(angles [ pointIndex ]);
            rot.at(2, 1) =  sin(angles [ pointIndex ]);
            rot.at(2, 2) =  cos(angles [ pointIndex ]);

            FloatArray tRot, nRot;
            tRot.beProductOf(rot, t);
            nRot.beProductOf(rot, n);

            FloatMatrix rotTot(2, 2);
            rotTot.setColumn(tRot, 1);
            rotTot.setColumn(nRot, 2);


            FloatMatrix tmp, damageRot;

            tmp.beTProductOf(rotTot, damageTen);
            damageRot.beProductOf(tmp, rotTot);


            const double omegaThetaTheta      =               damageRot.at(2, 2);
            omegaTTArray.push_back(omegaThetaTheta);

            const double omegaRTheta          =               damageRot.at(1, 2);
            omegaRTArray.push_back(omegaRTheta);
        }

        //////////////////////////////
        // Compute propagation angle

        // Find angles that fulfill sigRT = 0
        //const double damageTol = 1.0e-9;
        double maxOmegaTT = 0.0, maxAngle = 0.0;
        //bool foundZeroLevel = false;
        /*
        for ( size_t segIndex = 0; segIndex < ( circPoints.size() - 1 ); segIndex++ ) {
            // If the shear damage damageRT changes sign over the segment
            if ( omegaRTArray [ segIndex ] * omegaRTArray [ segIndex + 1 ] < damageTol ) {
                // Compute location of zero level
                double xi = EnrichmentItem :: calcXiZeroLevel(omegaRTArray[segIndex], omegaRTArray[segIndex + 1]);

                double theta                    = 0.5 * ( 1.0 - xi ) * angles [ segIndex ]         + 0.5 * ( 1.0 + xi ) * angles [ segIndex + 1 ];
                double omegaThetaTheta    = 0.5 * ( 1.0 - xi ) * omegaTTArray [ segIndex ] + 0.5 * ( 1.0 + xi ) * omegaTTArray [ segIndex + 1 ];

                //printf("Found candidate: theta: %e sigThetaTheta: %e\n", theta, sigThetaTheta);

                if ( omegaThetaTheta > maxOmegaTT ) {
                    foundZeroLevel = true;
                    maxOmegaTT = omegaThetaTheta;
                    maxAngle = theta;
                }
            }
        }
        */

        for ( size_t segIndex = 0; segIndex < ( circPoints.size() ); segIndex++ ) {
            if ( omegaTTArray[ segIndex ] > maxOmegaTT ) {
                    //foundZeroLevel = true;
                    maxOmegaTT = omegaTTArray[ segIndex ];
                    maxAngle = angles [ segIndex ];
                }
        }

        //if ( !foundZeroLevel ) {
        //    return false;
        //}

        if ( iDomain.giveXfemManager()->giveVtkDebug() ) {
            XFEMDebugTools :: WriteArrayToMatlab("sigTTvsAngle.m", angles, omegaTTArray);
            XFEMDebugTools :: WriteArrayToMatlab("sigRTvsAngle.m", angles, omegaRTArray);

            //XFEMDebugTools :: WriteArrayToGnuplot("sigTTvsAngle.dat", angles, omegaTTArray);
            //XFEMDebugTools :: WriteArrayToGnuplot("sigRTvsAngle.dat", angles, omegaRTArray);
        }

        // Compare with threshold
        printf("maxOmegaTT: %e mDamageThreshold: %e\n", maxOmegaTT, mDamageThreshold);
        //if ( maxOmegaTT > mDamageThreshold && foundZeroLevel ) {
        if ( maxOmegaTT > mDamageThreshold ) {
            // Rotation matrix
            FloatMatrix rot(2, 2);
            rot.at(1, 1) =  cos(maxAngle);
            rot.at(1, 2) = -sin(maxAngle);
            rot.at(2, 1) =  sin(maxAngle);
            rot.at(2, 2) =  cos(maxAngle);

            FloatArray dir;
            dir.beProductOf(rot, tipInfo.mTangDir);

            // Fill up struct
            oTipProp.mTipIndex = tipInfo.mTipIndex;
            oTipProp.mPropagationDir = dir;
            oTipProp.mPropagationLength = mIncrementLength;

            return true;
        }
    }


    return false;
}
} // end namespace oofem
