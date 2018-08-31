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

#include "calmls.h"
#include "verbose.h"
#include "timestep.h"
#include "floatmatrix.h"
#include "datastream.h"
#include "mathfem.h"
#include "element.h"
#include "Elements/structuralelement.h"
#include "../sm/Materials/structuralmaterial.h"
#include "../sm/CrossSections/simplecrosssection.h"

#include "classfactory.h"
#include "engngm.h"
// includes for HPC - not very clean (NumMethod knows what is "node" and "dof")
#include "node.h"
#include "dof.h"
#include "contextioerr.h"
#include "exportmodulemanager.h"
#include "parallelcontext.h"
#include "unknownnumberingscheme.h"
#include <cmath>
#include <list>

namespace oofem {
#define CALM_RESET_STEP_REDUCE 0.05
#define CALM_TANGENT_STIFF_TRESHOLD 0.1
#define CALM_DEFAULT_NRM_TICKS 2
#define CALM_MAX_REL_ERROR_BOUND 1.e10

REGISTER_SparseNonLinearSystemNM(CylindricalALM)

CylindricalALM :: CylindricalALM(Domain *d, EngngModel *m) :
    SparseNonLinearSystemNM(d, m), calm_HPCIndirectDofMask(), calm_HPCDmanDofSrcArray(), ccDofGroups()
{
    nsmax  = 60;       // default maximum number of sweeps allowed
    numberOfRequiredIterations = 3;
    //rtol   = 10.E-3  ;   // convergence tolerance
    //Psi    = 0.1;       // displacement control on
    solved = 0;
    calm_NR_Mode = calm_NR_OldMode = calm_modifiedNRM;
    calm_NR_ModeTick = -1; // do not swith to calm_NR_OldMode
    calm_MANRMSteps = 0;

    //Bergan_k0 = 0.;    // value used for computing Bergan's parameter
    // of current stiffness.

    deltaL    = -1.0;
    // TangenStiffnessTreshold = 0.10;

    // Variables for Hyper Plane Control
    calm_Control = calm_hpc_off; // HPControl is not default
    // linesearch default off

    // number of convergence_criteria dof groups, set to 0 (default behavior)
    nccdg = 0;

    minIterations = 0;

    calm_hpc_init = 0;

    // Maximum number of restarts when convergence not reached during maxiter
    maxRestarts = 3;

    ProElesize=0;

    parallel_context = engngModel->giveParallelContext( d->giveNumber() );
}


CylindricalALM :: ~CylindricalALM()
{
}


NM_Status
CylindricalALM :: solve(SparseMtrx &k, FloatArray &R, FloatArray *R0, FloatArray *iR,
                        FloatArray &X, FloatArray &dX, FloatArray &F,
                        const FloatArray &internalForcesEBENorm, double &ReachedLambda, referenceLoadInputModeType rlm,
                        int &nite, TimeStep *tStep)
{
    FloatArray rhs, DeltaX1, deltaX_P, deltaX_F, XInitial,dXm1,deltaX_;
    FloatArray ddX; // total increment of displacements in iteration
    //double Bergan_k0 = 1.0, bk;
    double XX, RR, RR0, p = 0.0;
    double deltaLambda, Lambda, DeltaLambdam1, DeltaLambda = 0.0;
    double drProduct = 0.0;
    int neq = R.giveSize();
    int irest = 0;
    int i, NodePerEle=3;
    double _RR, _XX;
    NM_Status status;
    bool converged, errorOutOfRangeFlag;
    // print iteration header


    OOFEM_LOG_INFO("CALMLS:       Initial step length: %-15e\n", deltaL);
    if ( nccdg == 0 ) {
        OOFEM_LOG_INFO("CALMLS:       Iteration       LoadLevel       ForceError      DisplError    \n");
        OOFEM_LOG_INFO("----------------------------------------------------------------------------\n");
    } else {
        OOFEM_LOG_INFO("Iter  LoadLevel       ");
        for ( i = 1; i <= nccdg; i++ ) {
            OOFEM_LOG_INFO("ForceError(%02d)  DisplError(%02d)  ", i, i);
        }

        OOFEM_LOG_INFO("\n__________________________________________________________\n");
    }

    //
    // Now smarter method is default:
    // after convergence troubles for (calm_NR_ModeTick subsequent
    // steps newly set calm_NR_Mode will be used. After the calm_NR_OldMode will be restored
    //
    if ( calm_NR_ModeTick == 0 ) {
        calm_NR_Mode = calm_NR_OldMode;
    }

    if ( calm_NR_ModeTick > 0 ) {
        calm_NR_ModeTick--;
    }

    XInitial = X;
    ddX.resize(neq);
    ddX.zero();
    deltaX_P.resize(neq);
    deltaX_P.zero();
    deltaX_F.resize(neq);
    deltaX_F.zero();

    status = NM_None;
    this->giveLinearSolver();


    /*
    std::list<int> mylist;
    int nelem = domain->giveNumberOfElements();
    for ( int ielem = 1; ielem <= nelem; ielem++ ) {
        Element *element = domain->giveElement(ielem);
        int answer=0;
        StructuralElement* Ele=static_cast<StructuralElement*>(domain->giveElement(ielem));
        SimpleCrossSection* Cros=static_cast<SimpleCrossSection*>(Ele->giveStructuralCrossSection());
        for ( auto &gp : *element->giveDefaultIntegrationRulePtr()){
            Cros->giveDamageState_PlaneStrain(answer,gp);
            if (answer){
                IntArray elenodes=element->giveDofManArray();
                for(int i=1; i<=elenodes.giveSize(); i++){
                    mylist.push_back(elenodes.at(i));
                }
                NodePerEle = elenodes.giveSize();
                break;
            }
        }
    } 
    if ( mylist.size() != 0 ){
        ProElesize=mylist.size()/NodePerEle;
        calm_HPCDmanDofSrcArray.resize(4*4*ProElesize);
        int j=1;
        for (std::list<int>::const_iterator iterator = mylist.begin(), end = mylist.end(); iterator != end; ++iterator) {
            calm_HPCDmanDofSrcArray.at(4*j-3)=*iterator;
            calm_HPCDmanDofSrcArray.at(4*j-2)=1;
            calm_HPCDmanDofSrcArray.at(4*j-1)=*iterator;
            calm_HPCDmanDofSrcArray.at(4*j-0)=2;
            j++;
        }
        this->convertHPCMap();
    }
    */

    // create HPC Map if needed
    if ( calm_hpc_init ) {
        this->convertHPCMap();
        calm_hpc_init = 0;
    }
    if (calm_HPCDmanDofSrcArray.isEmpty()){
        calm_Control = calm_hpc_off;
    }else{
        calm_Control = calm_hpc_on;
    }

    if ( R0 ) {
        RR0 = parallel_context->localNorm(* R0);
        RR0 *= RR0;
    } else {
        RR0 = 0.0;
    }

    //
    // A  initial step (predictor)
    //
    //

    //
    // A.2.   We assume positive-definite (0)Kt (tangent stiffness mtrx).
    //
    RR = parallel_context->localNorm(R);
    RR *= RR;

restart:
    //
    // A.1. calculation of (0)VarRt
    //
    dX.zero();
    //engngModel->updateComponent(tStep, InternalRhs, domain); // By not updating this, one obtains the old equilibrated tangent.
    engngModel->updateComponent(tStep, NonLinearLhs, domain);
    linSolver->solve(k, R, deltaX_P);

    // If desired by the user, the solution is (slightly) perturbed, so that various symmetries can be broken.
    // This is useful e.g. to trigger localization in a homogeneous material under uniform stress without
    // the need to introduce material imperfections. The problem itself remains symmetric but the iterative
    // solution is brought to a nonsymmetric state and it gets a chance to converge to a nonsymmetric solution.
    // Parameters of the perturbation technique are specified by the user and by default no perturbation is done. 
    // Milan Jirasek
    //SparseNonLinearSystemNM :: applyPerturbation(&deltaX_P);

    if ( calm_Control == calm_hpc_off ) {
        XX = parallel_context->localNorm(deltaX_P);
        XX *= XX;
        p = sqrt(XX);
    } else if ( calm_Control == calm_hpc_on ) {
        _XX = 0;
        for(int i=1; i<=ProElesize; i++){
            IntArray indxy(6);
            indxy.at(1)=calm_HPCIndirectDofMask.at(i*NodePerEle*2-5);
            indxy.at(2)=calm_HPCIndirectDofMask.at(i*NodePerEle*2-4);
            indxy.at(3)=calm_HPCIndirectDofMask.at(i*NodePerEle*2-3);
            indxy.at(4)=calm_HPCIndirectDofMask.at(i*NodePerEle*2-2);
            indxy.at(5)=calm_HPCIndirectDofMask.at(i*NodePerEle*2-1);
            indxy.at(6)=calm_HPCIndirectDofMask.at(i*NodePerEle*2-0);
            FloatArray Temp_X_P(6);
            for (int j=1; j<=6; j++){
                if (indxy.at(j)){
                    Temp_X_P.at(j)=deltaX_P.at(indxy.at(j));
                }else{
                    Temp_X_P.at(j)=0;
                }
            }
            _XX += pow((Temp_X_P.at(1)-Temp_X_P.at(5)),2);
            _XX += pow((Temp_X_P.at(3)-Temp_X_P.at(1)),2);
            _XX += pow((Temp_X_P.at(5)-Temp_X_P.at(3)),2);
            _XX += pow((Temp_X_P.at(2)-Temp_X_P.at(6)),2);
            _XX += pow((Temp_X_P.at(4)-Temp_X_P.at(2)),2);
            _XX += pow((Temp_X_P.at(6)-Temp_X_P.at(4)),2);
        }
        // In case of paralllel analysis:
        FloatArray collected_XXRR;
        parallel_context->accumulate({_XX, _RR}, collected_XXRR);
        _XX = collected_XXRR(0);
        _RR = collected_XXRR(1);

        p = sqrt(_XX);
    }

    //XR = parallel_context->localDotProduct(deltaX_P, R);
    /* XR is unscaled Bergan's param of current stiffness XR = deltaX_P^T k deltaX_P
     * this is used to test whether k has negative or positive slope */

    Lambda = ReachedLambda;
    DeltaLambda = deltaLambda =  deltaL / p;
    Lambda += DeltaLambda;
    //
    // A.3.
    //

    rhs = R;
    rhs.times(DeltaLambda);
    if ( R0 ) {
        rhs.add(* R0);
    }
    linSolver->solve(k, rhs, dX);
    X.add(dX);

    DeltaX1=dX;

    nite = 0;

    // update solution state counter
    tStep->incrementStateCounter();
    engngModel->updateComponent(tStep, InternalRhs, domain);

    do {
        nite++;
        tStep->incrementSubStepNumber();

        dXm1 = dX;
        DeltaLambdam1 = DeltaLambda;

        //
        // B  - iteration MNRM is used
        //
        // B.1. is ommited because MNRM is used instead of NRM.
        //
        if ( ( calm_NR_Mode == calm_fullNRM ) || ( ( calm_NR_Mode == calm_accelNRM ) && ( nite % calm_MANRMSteps == 0 ) ) ) {
            //
            // ALM with full NRM
            //
            // we assemble new tangent stiffness and compute new deltaX_P
            // internal state of elements is updated by previous calling
            // InternalRhs
            //
            engngModel->updateComponent(tStep, NonLinearLhs, domain);
            //
            // compute deltaX_P for i-th iteration
            //
            linSolver->solve(k, R, deltaX_P);
        }

        // B.2.
        //

        rhs =  R;
        rhs.times(Lambda);
        if ( R0 ) {
            rhs.add(* R0);
        }

        rhs.subtract(F);
        deltaX_F.resize(neq);
        linSolver->solve(k, rhs, deltaX_F);



        //
        // B.3. Compute -- deltaLambda ---
        //

        if ( calm_Control == calm_hpc_off ) {
              // this two lines are necessary if NRM is used
              // (for MNRM they can be computed at startup A1).
              deltaX_.clear();
              deltaX_=dX + deltaX_F;
              double X1X_ = parallel_context->localDotProduct(DeltaX1, deltaX_);
              double X1Xp = parallel_context->localDotProduct(DeltaX1, deltaX_P);
              deltaLambda =  (deltaL*deltaL -X1X_) / X1Xp;

          } else if ( calm_Control == calm_hpc_on ) {
            _XX = 0;
            double X1Xp=0, X1X_=0;
            FloatArray deltaX_b(6);
            deltaX_b.zero();
            for(int i=1; i<=ProElesize; i++){
                IntArray indxy(6);
                indxy.at(1)=calm_HPCIndirectDofMask.at(i*NodePerEle*2-5);
                indxy.at(2)=calm_HPCIndirectDofMask.at(i*NodePerEle*2-4);
                indxy.at(3)=calm_HPCIndirectDofMask.at(i*NodePerEle*2-3);
                indxy.at(4)=calm_HPCIndirectDofMask.at(i*NodePerEle*2-2);
                indxy.at(5)=calm_HPCIndirectDofMask.at(i*NodePerEle*2-1);
                indxy.at(6)=calm_HPCIndirectDofMask.at(i*NodePerEle*2-0);

                FloatArray Temp_X_P(6),Temp_X1(6),Temp_dX(6), Temp_X_F(6);
                for (int j=1; j<=6; j++){
                    if (indxy.at(j)){
                        Temp_X_P.at(j)=deltaX_P.at(indxy.at(j));
                        Temp_X1.at(j)=DeltaX1.at(indxy.at(j));
                        Temp_dX.at(j)=dX.at(indxy.at(j));
                        Temp_X_F.at(j)=deltaX_F.at(indxy.at(j));
                    }else{
                        Temp_X_P.at(j)=0;
                        Temp_X1.at(j)=0;
                        Temp_dX.at(j)=0;
                        Temp_X_F.at(j)=0;
                    }
                }

                X1Xp += (Temp_X_P.at(1)-Temp_X_P.at(5))*(Temp_X1.at(1)-Temp_X1.at(5));
                X1Xp += (Temp_X_P.at(3)-Temp_X_P.at(1))*(Temp_X1.at(3)-Temp_X1.at(1));
                X1Xp += (Temp_X_P.at(5)-Temp_X_P.at(3))*(Temp_X1.at(5)-Temp_X1.at(3));
                X1Xp += (Temp_X_P.at(2)-Temp_X_P.at(6))*(Temp_X1.at(2)-Temp_X1.at(6));
                X1Xp += (Temp_X_P.at(4)-Temp_X_P.at(2))*(Temp_X1.at(4)-Temp_X1.at(2));
                X1Xp += (Temp_X_P.at(6)-Temp_X_P.at(4))*(Temp_X1.at(6)-Temp_X1.at(4));

                deltaX_b.at(1)=(Temp_dX.at(1)-Temp_dX.at(5))+(Temp_X_F.at(1)-Temp_X_F.at(5));
                deltaX_b.at(2)=(Temp_dX.at(3)-Temp_dX.at(5))+(Temp_X_F.at(3)-Temp_X_F.at(5));
                deltaX_b.at(3)=(Temp_dX.at(5)-Temp_dX.at(3))+(Temp_X_F.at(5)-Temp_X_F.at(3));
                deltaX_b.at(4)=(Temp_dX.at(2)-Temp_dX.at(6))+(Temp_X_F.at(2)-Temp_X_F.at(6));
                deltaX_b.at(5)=(Temp_dX.at(4)-Temp_dX.at(2))+(Temp_X_F.at(4)-Temp_X_F.at(2));
                deltaX_b.at(6)=(Temp_dX.at(6)-Temp_dX.at(4))+(Temp_X_F.at(6)-Temp_X_F.at(4));

                X1X_ += (Temp_X1.at(1)-Temp_X1.at(5))*deltaX_b.at(1);
                X1X_ += (Temp_X1.at(3)-Temp_X1.at(1))*deltaX_b.at(2);
                X1X_ += (Temp_X1.at(5)-Temp_X1.at(3))*deltaX_b.at(3);
                X1X_ += (Temp_X1.at(2)-Temp_X1.at(6))*deltaX_b.at(4);
                X1X_ += (Temp_X1.at(4)-Temp_X1.at(2))*deltaX_b.at(5);
                X1X_ += (Temp_X1.at(6)-Temp_X1.at(4))*deltaX_b.at(6);
            }
            // In case of paralllel analysis:
            FloatArray collected_XXRR;
            parallel_context->accumulate({X1Xp, X1X_}, collected_XXRR);
            X1Xp = collected_XXRR(0);
            X1X_ = collected_XXRR(1);
            deltaLambda =  (deltaL*deltaL -X1X_) / X1Xp;
        }

        //
        // B.4.+ B.5.
        //
            //
            // update solution vectors
            //

            ddX.clear();
            ddX=deltaX_P;
            ddX.times(deltaLambda);
            ddX.add(deltaX_F);
            dX = dXm1;
            dX.add(ddX);
            X = XInitial;
            X.add(dX);

            drProduct = parallel_context->localNorm(ddX);
            drProduct *= drProduct;

            tStep->incrementStateCounter();     // update solution state counter
            //
            // B.6.
            //
            DeltaLambda = DeltaLambdam1 + deltaLambda;
            Lambda = ReachedLambda + DeltaLambda;

            tStep->incrementStateCounter();     // update solution state counter
            engngModel->updateComponent(tStep, InternalRhs, domain);

        //
        // B.7.
        //
        // convergence check
        //

        converged = this->checkConvergence(R, R0, F, X, ddX, Lambda, RR0, RR, drProduct,
                                           internalForcesEBENorm, nite, errorOutOfRangeFlag);
        if ( ( nite >= nsmax ) || errorOutOfRangeFlag ) {
            irest++;
            if ( irest <= maxRestarts ) {
                // convergence problems
                // there must be step restart followed by decrease of step length
                // status |= NM_ForceRestart;
                // reduce step length
                deltaL =  deltaL * CALM_RESET_STEP_REDUCE;
                if ( deltaL < minStepLength ) {
                    deltaL = minStepLength;
                }

                // restore previous total displacement vector
                X = XInitial;
                // reset all changes from previous equilibrium state
                engngModel->initStepIncrements();
                dX.zero();
                // restore initial stiffness
                engngModel->updateComponent(tStep, NonLinearLhs, domain);

                OOFEM_LOG_INFO("CALMLS:       Iteration Reset ...\n");

                //calm_NR_OldMode  = calm_NR_Mode;
                calm_NR_Mode     = calm_fullNRM;
                //calm_NR_ModeTick = CALM_DEFAULT_NRM_TICKS;
                goto restart;
            } else {
                status = NM_NoSuccess;
                OOFEM_WARNING("Convergence not reached after %d iterations", nsmax);
                // exit(1);
                break;
            }
        }

        // output of per iteration data
        engngModel->giveExportModuleManager()->doOutput(tStep, true);
    } while ( !converged || ( nite < minIterations ) );

    //
    // update dofs,nodes,Elemms and print result
    //
#ifdef VERBOSE
    // printf ("\nCALM - step iteration finished") ;
#endif

    /* we have computed Bergan's  parameter -> you can adjust deltaL according
     * to this value */
    //
    // there has been restart already - set nite to maxiter
    //
    // if (irest > 0) nite = nsmax;

    if ( nite > numberOfRequiredIterations ) {
        deltaL =  deltaL * numberOfRequiredIterations / nite;
    } else {
        deltaL =  deltaL * sqrt( sqrt( ( double ) numberOfRequiredIterations / ( double ) nite ) );
    }

    if ( deltaL > maxStepLength ) {
        deltaL = maxStepLength;
    }

    if ( deltaL < minStepLength ) {
        deltaL = minStepLength;
    }

    OOFEM_LOG_INFO("CALMLS:       Adjusted step length: %-15e\n", deltaL);

    status = NM_Success;
    solved = 1;
    ReachedLambda = Lambda;

    return status;
}

bool
CylindricalALM :: checkConvergence(const FloatArray &R, const FloatArray *R0, const FloatArray &F,
                                   const FloatArray &X, const FloatArray &ddX,
                                   double Lambda, double RR0, double RR, double drProduct,
                                   const FloatArray &internalForcesEBENorm, int nite, bool &errorOutOfRange)
{

    int _ng = nccdg;
    double forceErr, dispErr;
    FloatArray rhs; // residual of momentum balance eq (unbalanced nodal forces)
    FloatArray dg_forceErr(nccdg), dg_dispErr(nccdg), dg_totalLoadLevel(nccdg), dg_totalDisp(nccdg);
    bool answer;
    EModelDefaultEquationNumbering dn;

    answer = true;
    errorOutOfRange = false;

    // compute residual vector
    rhs =  R;
    rhs.times(Lambda);
    if ( R0 ) {
        rhs.add(* R0);
    }

    rhs.subtract(F);

    if ( _ng > 0 ) {
        forceErr = dispErr = 0.0;
        // zero error norms per group
        dg_forceErr.zero();
        dg_dispErr.zero();
        dg_totalLoadLevel.zero();
        dg_totalDisp.zero();
        // loop over dof managers
        for ( auto &dman : domain->giveDofManagers() ) {
            if ( !dman->isLocal() ) {
                continue;
            }

            // loop over individual dofs
            for ( Dof *_idofptr: *dman ) {
                // loop over dof groups
                for ( int _dg = 1; _dg <= _ng; _dg++ ) {
                    // test if dof ID is in active set
                    if ( ccDofGroups.at(_dg - 1).find( _idofptr->giveDofID() ) != ccDofGroups.at(_dg - 1).end() ) {
                        int _eq = _idofptr->giveEquationNumber(dn);

                        if ( _eq ) {
                            continue;
                        }
                        dg_forceErr.at(_dg) += rhs.at(_eq) * rhs.at(_eq);
                        dg_dispErr.at(_dg)  += ddX.at(_eq) * ddX.at(_eq);
                        // missing - compute norms of total displacement and load vectors (but only for selected dofs)!
                        if ( R0 ) {
                            dg_totalLoadLevel.at(_dg) += R0->at(_eq) * R0->at(_eq);
                        }

                        dg_totalLoadLevel.at(_dg) += R.at(_eq) * R.at(_eq) * Lambda * Lambda;
                        dg_totalDisp.at(_dg) += X.at(_eq) * X.at(_eq);
                    }
                } // end loop over dof groups
            } // end loop over DOFs
        } // end loop over dof managers

        // loop over elements and their DOFs
        for ( auto &elem : domain->giveElements() ) {
            if ( elem->giveParallelMode() != Element_local ) {
                continue;
            }

            // loop over element internal Dofs
            for ( int _idofman = 1; _idofman <= elem->giveNumberOfInternalDofManagers(); _idofman++ ) {
                // loop over individual dofs
                for ( Dof *_idofptr: *elem->giveInternalDofManager(_idofman) ) {
                    // loop over dof groups
                    for ( int _dg = 1; _dg <= _ng; _dg++ ) {
                        // test if dof ID is in active set
                        if ( ccDofGroups.at(_dg - 1).find( _idofptr->giveDofID() ) != ccDofGroups.at(_dg - 1).end() ) {
                            int _eq = _idofptr->giveEquationNumber(dn);

                            if ( _eq ) {
                                continue;
                            }
                            dg_forceErr.at(_dg) += rhs.at(_eq) * rhs.at(_eq);
                            dg_dispErr.at(_dg)  += ddX.at(_eq) * ddX.at(_eq);
                            // missing - compute norms of total displacement and load vectors (but only for selected dofs)!
                            if ( R0 ) {
                                dg_totalLoadLevel.at(_dg) += R0->at(_eq) * R0->at(_eq);
                            }

                            dg_totalLoadLevel.at(_dg) += R.at(_eq) * R.at(_eq) * Lambda * Lambda;
                            dg_totalDisp.at(_dg) += X.at(_eq) * X.at(_eq);
                        }
                    } // end loop over dof groups
                } // end loop over DOFs
            } // end loop over internal element dofmans
        } // end loop over elements

        // exchange individual partition contributions (simultaneously for all groups)
        FloatArray collectiveErr;
        parallel_context->accumulate(dg_forceErr, collectiveErr);
        dg_forceErr = collectiveErr;
        parallel_context->accumulate(dg_dispErr, collectiveErr);
        dg_dispErr = collectiveErr;
        parallel_context->accumulate(dg_totalLoadLevel, collectiveErr);
        dg_totalLoadLevel = collectiveErr;
        parallel_context->accumulate(dg_totalDisp, collectiveErr);
        dg_totalDisp = collectiveErr;

        OOFEM_LOG_INFO("CALMLS:       %-15d %-15e ", nite, Lambda);
        // loop over dof groups
        for ( int _dg = 1; _dg <= _ng; _dg++ ) {
            //  compute a relative error norm
            if ( ( dg_totalLoadLevel.at(_dg) ) < calm_SMALL_ERROR_NUM ) {
                dg_forceErr.at(_dg) = sqrt( dg_forceErr.at(_dg) );
            } else {
                dg_forceErr.at(_dg) = sqrt( dg_forceErr.at(_dg) / dg_totalLoadLevel.at(_dg) );
            }

            //
            // compute displacement error
            //
            if ( dg_totalDisp.at(_dg) < calm_SMALL_ERROR_NUM ) {
                dg_dispErr.at(_dg) = sqrt( dg_dispErr.at(_dg) );
            } else {
                dg_dispErr.at(_dg) = sqrt( dg_dispErr.at(_dg) / dg_totalDisp.at(_dg) );
            }

            if ( ( fabs( dg_forceErr.at(_dg) ) > rtolf.at(_dg) * CALM_MAX_REL_ERROR_BOUND ) ||
                ( fabs( dg_dispErr.at(_dg) )  > rtold.at(_dg) * CALM_MAX_REL_ERROR_BOUND ) ) {
                errorOutOfRange = true;
            }

            if ( ( fabs( dg_forceErr.at(_dg) ) > rtolf.at(_dg) ) || ( fabs( dg_dispErr.at(_dg) ) > rtold.at(_dg) ) ) {
                answer = false;
            }


            OOFEM_LOG_INFO( "%-15e %-15e ", dg_forceErr.at(_dg), dg_dispErr.at(_dg) );
        }

        OOFEM_LOG_INFO("\n");
    } else {
        //
        // _ng==0 (errors computed for all dofs - this is the default)
        //

        //
        // compute force error(s)
        //
        double dXX;
        forceErr = parallel_context->localNorm(rhs);
        forceErr *= forceErr;
        dXX = parallel_context->localNorm(X);
        dXX *= dXX;

        double eNorm = internalForcesEBENorm.sum();
        // we compute a relative error norm
        if ( ( RR0 + RR * Lambda * Lambda ) > calm_SMALL_ERROR_NUM ) {
            forceErr = sqrt( forceErr / ( RR0 + RR * Lambda * Lambda ) );
        } else if ( eNorm > calm_SMALL_ERROR_NUM ) {
            forceErr = sqrt(forceErr / eNorm);
        } else {
            forceErr = sqrt(forceErr);
        }

        //
        // compute displacement error
        //
        if ( dXX < calm_SMALL_ERROR_NUM ) {
            dispErr = drProduct;
        } else {
            dispErr = drProduct / dXX;
            dispErr = sqrt(dispErr);
        }

        if ( ( fabs(forceErr) > rtolf.at(1) * CALM_MAX_REL_ERROR_BOUND ) ||
            ( fabs(dispErr)  > rtold.at(1) * CALM_MAX_REL_ERROR_BOUND ) ) {
            errorOutOfRange = true;
        }

        if ( ( fabs(forceErr) > rtolf.at(1) ) || ( fabs(dispErr) > rtold.at(1) ) ) {
            answer = false;
        }

        OOFEM_LOG_INFO("CALMLS:       %-15d %-15e %-15e %-15e\n", nite, Lambda, forceErr, dispErr);
    } // end default case (all dofs conributing)

    return answer;
}



IRResultType
CylindricalALM :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                   // Required by IR_GIVE_FIELD macro

    double initialStepLength, forcedInitialStepLength;
    int hpcMode;

    nsmax = 30;
    IR_GIVE_OPTIONAL_FIELD(ir, nsmax, _IFT_CylindricalALM_maxiter);
    if ( nsmax < 30 ) {
        nsmax = 30;
    }

    IR_GIVE_OPTIONAL_FIELD(ir, maxRestarts, _IFT_CylindricalALM_maxrestarts);

    minStepLength = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, minStepLength, _IFT_CylindricalALM_minsteplength);

    IR_GIVE_FIELD(ir, maxStepLength, _IFT_CylindricalALM_steplength);
    initialStepLength = maxStepLength;
    IR_GIVE_OPTIONAL_FIELD(ir, initialStepLength, _IFT_CylindricalALM_initialsteplength);

    //if (deltaL <= 0.0)  deltaL=maxStepLength;
    // This method (instanciate) is called not only at the beginning but also
    // after restart from engngModel updateAttributes -> and in this case
    // we want to keep restored deltaL)
    if ( ( deltaL <= 0.0 ) || ( deltaL > maxStepLength ) ) {
        deltaL = initialStepLength;
    }

    // using this, one can enforce deltaL from the input file after restart
    forcedInitialStepLength = 0.;
    IR_GIVE_OPTIONAL_FIELD(ir, forcedInitialStepLength, _IFT_CylindricalALM_forcedinitialsteplength);
    if ( forcedInitialStepLength > 0. ) {
        deltaL = forcedInitialStepLength;
    }

    numberOfRequiredIterations = 3;
    IR_GIVE_OPTIONAL_FIELD(ir, numberOfRequiredIterations, _IFT_CylindricalALM_reqiterations);
    if ( numberOfRequiredIterations < 3 ) {
        numberOfRequiredIterations = 3;
    }

    if ( numberOfRequiredIterations > 1000 ) {
        numberOfRequiredIterations = 1000;
    }

    IR_GIVE_OPTIONAL_FIELD(ir, minIterations, _IFT_CylindricalALM_miniterations);
    if ( result == IRRT_OK ) {
        if ( minIterations > 3 && minIterations < 1000 ) {
            numberOfRequiredIterations = minIterations;
        }

        if ( nsmax <= minIterations ) {
            nsmax = minIterations + 1;
        }
    }

    // read if MANRM method is used
    calm_MANRMSteps = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, calm_MANRMSteps, _IFT_CylindricalALM_manrmsteps);
    if ( calm_MANRMSteps > 0 ) {
        calm_NR_Mode = calm_NR_OldMode = calm_accelNRM;
    } else {
        calm_NR_Mode = calm_modifiedNRM;
    }

    // read if HPC is requested
    hpcMode = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, hpcMode, _IFT_CylindricalALM_hpcmode);

    calm_HPCDmanDofSrcArray.clear();
    IR_GIVE_OPTIONAL_FIELD(ir, calm_HPCDmanDofSrcArray, _IFT_CylindricalALM_hpc);
    ProElesize=calm_HPCDmanDofSrcArray.giveSize()/12;


    // in calm_HPCIndirectDofMask are stored pairs with following meaning:
    // inode idof
    // example HPC 4 1 2 6 1
    // will yield to HPC for node 4 dof 1 and node 6 dof 1
    // calm_HPCIndirectDofMask must be converted to indirect map
    // -> because dof eqs. are not known now, we derefer this to
    // solveYourselfAt() subroutine. The need for converting is indicated by
    // calm_HPControl = hpc_init
    if ( calm_HPCDmanDofSrcArray.giveSize() != 0 ) {
        if ( hpcMode == 1 ) {
            calm_Control = calm_hpc_on; // default is to use hpc_on
        }

        if ( ( calm_HPCDmanDofSrcArray.giveSize() % 2 ) != 0 ) {
            OOFEM_ERROR("HPC Map size must be even number, it contains pairs <node, nodeDof>");
        }

        calm_hpc_init = 1;
    } else {
        if ( hpcMode ) {
            OOFEM_ERROR("HPC Map must be specified");
        }
    }

    int _value = 0;
    IR_GIVE_OPTIONAL_FIELD(ir, _value, _IFT_CylindricalALM_lstype);
    solverType = ( LinSystSolverType ) _value;

    /* initialize optional dof groups for convergence criteria evaluation */
    this->nccdg = 0; // default, no dof cc group, all norms evaluated for all dofs
    IR_GIVE_OPTIONAL_FIELD(ir, nccdg, _IFT_CylindricalALM_nccdg);

    if ( nccdg >= 1 ) {
        IntArray _val;
        char name [ 12 ];
        // create an empty set
        __DofIDSet _set;
        // resize dof group vector
        this->ccDofGroups.resize(nccdg, _set);
        for ( int _i = 0; _i < nccdg; _i++ ) {
            sprintf(name, "%s%d", _IFT_CylindricalALM_ccdg, _i + 1);
            // read dof group as int array under ccdg# keyword
            IR_GIVE_FIELD(ir, _val, name);
            // convert aray into set
            for ( int _j = 1; _j <= _val.giveSize(); _j++ ) {
                ccDofGroups.at(_i).insert( ( DofIDItem ) _val.at(_j) );
            }
        }

        // read relative error tolerances of the solver fo each cc
        // if common rtolv provided, set to this tolerace both rtolf and rtold
        IR_GIVE_OPTIONAL_FIELD(ir, rtolf, _IFT_CylindricalALM_rtolv);
        rtold = rtolf;
        // read optional force and displacement tolerances
        IR_GIVE_OPTIONAL_FIELD(ir, rtolf, _IFT_CylindricalALM_rtolf);
        IR_GIVE_OPTIONAL_FIELD(ir, rtold, _IFT_CylindricalALM_rtold);

        if ( ( rtolf.giveSize() != nccdg ) || ( rtold.giveSize() != nccdg ) ) {
            OOFEM_ERROR("Incompatible size of rtolf or rtold params, expected size %d (nccdg)", nccdg);
        }
    } else {
        nccdg = 0;
        double _rtol = 1.e-3; // default tolerance
        rtolf.resize(1);
        rtold.resize(1);
        // read relative error tolerances of the solver
        // if common rtolv provided, set to this tolerace both rtolf and rtold
        IR_GIVE_OPTIONAL_FIELD(ir, _rtol, _IFT_CylindricalALM_rtolv);
        rtolf.at(1) = rtold.at(1) = _rtol;
        IR_GIVE_OPTIONAL_FIELD(ir, _rtol, _IFT_CylindricalALM_rtolf);
        rtolf.at(1) = _rtol;
        IR_GIVE_OPTIONAL_FIELD(ir, _rtol, _IFT_CylindricalALM_rtold);
        rtold.at(1) = _rtol;
    }

    this->giveLinearSolver()->initializeFrom(ir);

    SparseNonLinearSystemNM :: initializeFrom(ir);

    return IRRT_OK;
}


void CylindricalALM :: convertHPCMap()
//
// converts HPC map from user input map to HPC indirect Map
//
// indirect map:
// size of indirect map is number of controlled DOFs
// on i-th position this map contain equation number of i-th controlled dof
//
// user input map;
// user input map has size 2*controlled DOFs
// and contain pairs (node, nodeDof);
//
// This is used in order to hide equation numbering from user
//
{
    IntArray indirectMap;
    int size;
    EModelDefaultEquationNumbering dn;

    int count = 0;
    size = calm_HPCDmanDofSrcArray.giveSize() / 2;
    indirectMap.resize(size);
    for ( int i = 1; i <= size; i++ ) {
        int inode = calm_HPCDmanDofSrcArray.at(2 * i - 1);
        int idofid = calm_HPCDmanDofSrcArray.at(2 * i);
        for ( auto &dman : domain->giveDofManagers() ) {
            int jglobnum = dman->giveLabel();
            if ( inode == jglobnum ) {
                if ( parallel_context->isLocal( dman.get() ) ) {
                    indirectMap.at(++count) = dman->giveDofWithID(idofid)->giveEquationNumber(dn);
                }
                break;
            }
        }
    }

    if ( count != size ) {
        OOFEM_WARNING("some dofmans/Dofs in HPCarray not recognized");
    }

    calm_HPCIndirectDofMask.resize(count);

    for ( int i = 1; i <= count; i++ ) {
        calm_HPCIndirectDofMask.at(i) = indirectMap.at(i);
    }
}


contextIOResultType
CylindricalALM :: saveContext(DataStream &stream, ContextMode mode, void *obj)
{
    // write current deltaL
    if ( !stream.write(deltaL) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}


contextIOResultType
CylindricalALM :: restoreContext(DataStream &stream, ContextMode mode, void *obj)
{
    // read last deltaL
    if ( !stream.read(deltaL) ) {
        THROW_CIOERR(CIO_IOERR);
    }

    return CIO_OK;
}


SparseLinearSystemNM *
CylindricalALM :: giveLinearSolver()
{
    if ( linSolver ) {
        if ( linSolver->giveLinSystSolverType() == solverType ) {
            return linSolver.get();
        }
    }

    linSolver.reset( classFactory.createSparseLinSolver(solverType, domain, engngModel) );
    if ( !linSolver ) {
        OOFEM_ERROR("linear solver creation failed for lstype %d", solverType);
    }

    return linSolver.get();
}


} // end namespace oofem
