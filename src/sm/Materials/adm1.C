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
 *---
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

// This code is based on the anisotropic damage model proposed by Desmorat, Gatuingt and Ragueneau in
// their paper "Nonlocal anisotropic damage model and related computational aspects for quasi-brittle material"
// published in Engineering Fracture Mechanics 74 (2007) 1539-1560.

#include "adm1.h"
#include "../sm/Materials/structuralmaterial.h"
#include "floatmatrix.h"
#include "floatarray.h"
#include "mathfem.h"
#include "datastream.h"
#include "contextioerr.h"
#include "dynamicinputrecord.h"
#include "gausspoint.h"
#include "classfactory.h"
#include <cmath>

namespace oofem {
REGISTER_Material(AnisotropicDamageMaterial1)

AnisotropicDamageMaterial1 :: AnisotropicDamageMaterial1(int n, Domain *d) : IsotropicLinearElasticMaterial(n, d),
    RandomMaterialExtensionInterface()
    //
    // constructor
    //
{
    //linearElasticMaterial = new IsotropicLinearElasticMaterial(n, d);
    //E = 0.;
    //nu = 0.;
    kappa_t = 0;
    alpha_t = 0;
    kappa_c = 0.;
    alpha_c = 0.;
    eta = 0.;
}

AnisotropicDamageMaterial1 :: ~AnisotropicDamageMaterial1()
//
// destructor
//
{
    //delete linearElasticMaterial;
}

/*
Interface *
AnisotropicDamageMaterial1 :: giveInterface(InterfaceType type)
{
    if ( type == MaterialModelMapperInterfaceType ) {
        return static_cast< MaterialModelMapperInterface * >(this);
    } else {
        return NULL;
    }
}
*/

MaterialStatus *
AnisotropicDamageMaterial1 :: CreateStatus(GaussPoint *gp) const
{
    return new AnisotropicDamageMaterial1Status(1, domain, gp);
}


MaterialStatus *
AnisotropicDamageMaterial1 :: giveStatus(GaussPoint *gp) const
{
    MaterialStatus *status = static_cast< MaterialStatus * >( gp->giveMaterialStatus() );
    if ( status == NULL ) {
        // create a new one
        status = this->CreateStatus(gp);
        if ( status != NULL ) {
            gp->setMaterialStatus( status, this->giveNumber() );
            this->_generateStatusVariables(gp);
        }
    }
    return status;
}

int
AnisotropicDamageMaterial1 :: hasMaterialModeCapability(MaterialMode mode)
//
// returns whether receiver supports the given mode
//
{
    return mode == _PlaneStrain ||  mode == _3dMat;
    //return mode == _3dMat || mode == _PlaneStress || mode == _PlaneStrain || mode == _1dMat;
}

//********************************************************
// Implementation by Wencheng Jin
//********************************************************


void
AnisotropicDamageMaterial1 :: giveRealStressVector_PlaneStrain(FloatArray &answer, GaussPoint *gp,
                                                              const FloatArray &totalStrain, TimeStep *atTime)
//
// uses the following special methods:
//   computePrincValDir2D  ...  evaluation of eigenvalues and eigenvector of a symmetric 2x2 matrix
//
{
#define AD_TOLERANCE 1.e-6 // convergence tolerance for the internal iteration used under plane stress
    this->initTempStatus(gp); //updata the status, so the current iteration starts from the previous converged results;
    // subtract the stress-independent part of strains (e.g. due to temperature)
    FloatArray SDstrainVector;
    this->giveStressDependentPartOfStrainVector(SDstrainVector, gp, totalStrain, atTime, VM_Total);

    // make sure strain components are correctly stored as the following--WJ

    // compute in-plane principal strains: epsilon1 and epsilon2
    // and the components of the first principal strain direction: ceps and seps
    double epsilon1, epsilon2, ceps, seps,treps,Dc;
    this->computePrincValDir2D(epsilon1, epsilon2, ceps, seps, SDstrainVector.at(1), SDstrainVector.at(2), SDstrainVector.at(4) / 2.);

    FloatMatrix EigenVic(2,2);
    EigenVic.at(1,1)=ceps;
    EigenVic.at(1,2)=seps;
    EigenVic.at(2,1)=-seps;
    EigenVic.at(2,2)=ceps;

    treps=epsilon1+epsilon2;

    //added part
    //int FALG0=0;
    //if (epsilon1>0 || epsilon2>0){
    //    FALG0=1;
    //}

    AnisotropicDamageMaterial1Status *status = static_cast< AnisotropicDamageMaterial1Status * >( this->giveStatus(gp) );
    //Should check whether use tempdamage or damage value
    FloatMatrix tempDamage = status->giveTempDamage();

    double trtempomega,Yd=0, Fd, equivStrain;
    trtempomega=tempDamage.at(1,1)+tempDamage.at(2,2);
  
    //Compute equivalent strain
    this->computeEquivalentStrain(equivStrain, SDstrainVector, gp, atTime);

    if (treps > 0) {
        status->setTempFlag(0);
        //if (equivStrain > this->kappa_t){
        //    Fd = equivStrain - status->giveTempEquivalenstrain_t();
        //}else{
        //    Fd = -1;
        //}
        if (epsilon1 > 0) {
            Yd = Yd + epsilon1 * epsilon1;
        }
        if (epsilon2 > 0) {
            Yd = Yd + epsilon2 * epsilon2;
        }
        Fd = equivStrain - ( this->kappa_t + this->alpha_t * trtempomega );
        Dc = -this->nu;
    } else {
        status->setTempFlag(1);
        Fd = equivStrain + this->eta*treps - ( this->kappa_c + this->alpha_c * trtempomega );
        Dc = -2;
        Yd=equivStrain*equivStrain;
    }

    FloatMatrix MATD_2, StressTensor, StrainTensor;
    MATD_2.resize(6,6);
    MATD_2.zero();
    double MATD[3][3][3][3]={0};
    StressTensor.resize(3,3);
    StressTensor.zero();
    StrainTensor.resize(3,3);
    StrainTensor.zero();
    StrainTensor.at(1,1)=SDstrainVector.at(1);
    StrainTensor.at(1,2)=SDstrainVector.at(4)/2;
    StrainTensor.at(2,1)=SDstrainVector.at(4)/2;
    StrainTensor.at(2,2)=SDstrainVector.at(2);
    FloatMatrix Omega(2,2);
    Omega.zero();


    if ( Fd < AD_TOLERANCE ){
        if (treps > 0) {
            status->setTempEquivalenstrain_t(equivStrain);
        } else {
            status->setTempEquivalenstrain_c(equivStrain);
        }
        this->computeSecantOperator( MATD_2, MATD, Dc, tempDamage);

        this->Aijkl_Bkl(MATD, StrainTensor, StressTensor);
    } else {
        if (treps > 0) {

            double lambda = ( equivStrain - status->giveTempEquivalenstrain_t() ) / this->alpha_t;

            //double lambda = 1-exp(-(equivStrain - this->kappa_t)/this->alpha_t);

            status->setTempEquivalenstrain_t(equivStrain);
            if (epsilon1 > 0){
                Omega.at(1,1) = lambda*epsilon1*epsilon1/Yd;
            }
            if (epsilon2 > 0){
                Omega.at(2,2) = lambda*epsilon2*epsilon2/Yd;
            }
            // Omega.at(1,1) = lambda*((epsilon1 + abs(epsilon1))/2)*(( epsilon1+abs(epsilon1))/ 2)/Yd ;
            // Omega.at(2,2) = lambda*((epsilon2 + abs(epsilon2))/2)*(( epsilon2+abs(epsilon2))/ 2)/Yd ;
            FloatMatrix EigenVicInv = EigenVic;
            for (int i=1; i<=2; i++){
                EigenVic.at(1,i)=EigenVic.at(1,i)*Omega.at(i,i);
                EigenVic.at(2,i)=EigenVic.at(2,i)*Omega.at(i,i);
            }

            EigenVicInv.beInverseOf( EigenVicInv );
            Omega.beProductOf( EigenVic, EigenVicInv );

            tempDamage.add(Omega);
            status->setTempDamage( tempDamage );

            //status->setTempDamage( Omega );

        } else {
            double lambda = ( ( equivStrain - status->giveTempEquivalenstrain_c()) + this->eta* treps) / this->alpha_c;
            status->setTempEquivalenstrain_c(equivStrain);
            if ((epsilon1-treps/2) > 0){
                Omega.at(1,1) = lambda*(epsilon1-treps/2)*(epsilon1-treps/2)/Yd;
            }
            if ((epsilon2-treps/2) > 0){
                Omega.at(2,2) = lambda*(epsilon2-treps/2)*(epsilon2-treps/2)/Yd;
            }
            //Omega.at(1,1) = lambda*(((epsilon1-treps/2)+abs(epsilon1-treps/2))/2)*(((epsilon1-treps/2)+abs(epsilon1-treps/2))/2)/Yd;
            //Omega.at(2,2) = lambda*(((epsilon2-treps/2)+abs(epsilon2-treps/2))/2)*(((epsilon2-treps/2)+abs(epsilon2-treps/2))/2)/Yd;
            FloatMatrix EigenVicInv = EigenVic;
            for (int i=1; i<=2; i++){
                EigenVic.at(1,i)=EigenVic.at(1,i)*Omega.at(i,i);
                EigenVic.at(2,i)=EigenVic.at(2,i)*Omega.at(i,i);
            }
            EigenVicInv.beInverseOf( EigenVicInv );
            Omega.beProductOf( EigenVic, EigenVicInv );
            tempDamage.add(Omega);
            status->setTempDamage( tempDamage );
        }

        this->computeSecantOperator( MATD_2, MATD, Dc, tempDamage);

        this->Aijkl_Bkl(MATD, StrainTensor, StressTensor);
    }
    status->letTempStrainVectorBe(totalStrain);
    answer.resize(4);
    answer.at(1)=StressTensor.at(1,1);
    answer.at(2)=StressTensor.at(2,2);
    answer.at(3)=StressTensor.at(3,3);
    answer.at(4)=(StressTensor.at(1,2)+StressTensor.at(2,1))/2;
    status->letTempStressVectorBe(answer);
#ifdef keep_track_of_dissipated_energy
    status->computeWork(gp);
#endif
    return;
}

void
AnisotropicDamageMaterial1 :: computeEquivalentStrain(double &kappa, const FloatArray &strain, GaussPoint *gp, TimeStep *tStep)
{
    if ( strain.isEmpty() ) {
        kappa = 0.;
        return;
    }
    // compute in-plane principal strains: epsilon1 and epsilon2
    // and the components of the first principal strain direction: ceps and seps
    double epsilon1, epsilon2, ceps, seps, treps;
    this->computePrincValDir2D(epsilon1, epsilon2, ceps, seps, strain.at(1), strain.at(2), strain.at(4) / 2.);
    treps=epsilon1+epsilon2;

    double Yd=0;

    //if (treps>0) {
        if (epsilon1 > 0) {
            Yd = Yd + epsilon1 * epsilon1;
        }
        if (epsilon2 > 0) {
            Yd = Yd + epsilon2 * epsilon2;
        }
        kappa = sqrt( Yd );
        /*
    } else {
        if ( (epsilon1- treps /2) > 0 ) {
            Yd = Yd + (epsilon1- treps /2) * (epsilon1- treps /2);
        }
        if ( (epsilon2- treps /2) > 0 ){
            Yd = Yd + (epsilon2- treps /2) * (epsilon2- treps /2);
        }
        kappa = sqrt( Yd );
    }*/
    return;
}

int AnisotropicDamageMaterial1 ::giveDamageState(GaussPoint *gp)
{
    AnisotropicDamageMaterial1Status *status = static_cast< AnisotropicDamageMaterial1Status * >( this->giveStatus(gp) );

    FloatMatrix Damage(3,3);
    Damage = status->giveTempDamage();

    if (Damage.at(1,1) > 0.0 || Damage.at(1,2) > 0.0 ||Damage.at(2,1) > 0.0 || Damage.at(2,2) > 0.0){
        return 1;
    } else {
        return 0;
    };
}


void
AnisotropicDamageMaterial1 ::  Aijkl_Bkl(const double A[3][3][3][3], const FloatMatrix &B, FloatMatrix &C)
{

    C.resize(3,3);
    C.zero();
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            for ( int k = 1; k <= 3; k++ ) {
                for ( int l = 1; l <= 3; l++  ) {
                    C.at(i,j)= C.at(i,j) + A[i-1][j-1][k-1][l-1]*B.at(k,l);
                }
            }
        }
    }
    return;
}

void
AnisotropicDamageMaterial1 :: computePrincValDir2D(double &D1, double &D2, double &c, double &s, double Dx, double Dy, double Dxy)
//
// computes the principal values and directions of a symmetric second-order tensor in 2D
// input: Dx, Dy, Dxy ... components of the tensor wrt global coordinates
// output: D1, D2 ... ordered principal values, D1>=D2
// output: c, s ... components of the unit principal vector associated with D1
//                  (cosine and sine of the angle between the major principal direction and the global x-axis)
{
    double aux1 = ( Dx + Dy ) / 2.;
    double aux2 = ( Dx - Dy ) / 2.;
    double aux3 = sqrt(aux2 * aux2 + Dxy * Dxy);
    D1 = aux1 + aux3;
    D2 = aux1 - aux3;
    // formulae (90)-(92) and the two cases preceding them
    c = 1.;
    s = 0.;         // cases 1 and 2a
    if ( Dxy != 0. ) {  // Consider equal side triangle, with middle line trough
        double t = ( D1 - Dx ) / Dxy;
        c = 1. / sqrt(1. + t * t);
        s = c * t;
    } else if ( Dx < Dy ) { // case 2b
        c = 0.;
        s = 1.;
    }
    return;
}

/*
void
AnisotropicDamageMaterial1 :: giveRealStressVector_3d(FloatArray &answer, GaussPoint *gp,
                                                      const FloatArray &totalStrain, TimeStep *tStep)
{
}
void
AnisotropicDamageMaterial1 :: computePrincValDir3D(double &D1, double &D2, double &c, double &s, double Dx, double Dy, double Dxy)
//
// computes the principal values and directions of a symmetric second-order tensor in 2D
// input: Dx, Dy, Dxy ... components of the tensor wrt global coordinates
{
}

void AnisotropicDamageMaterial1 ::jacobi_eigenvalue ( int n, double a[], int it_max, double v[],
  double d[], int &it_num, int &rot_num )
//
//  Purpose:
//
//    JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
//
//  Discussion:
//
//    This function computes the eigenvalues and eigenvectors of a
//    real symmetric matrix, using Rutishauser's modfications of the classical
//    Jacobi rotation method with threshold pivoting.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    17 September 2013
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int N, the order of the matrix.
//
//    Input, double A[N*N], the matrix, which must be square, real,
//    and symmetric.
//
//    Input, int IT_MAX, the maximum number of iterations.
//
//    Output, double V[N*N], the matrix of eigenvectors.
//
//    Output, double D[N], the eigenvalues, in descending order.
//
//    Output, int &IT_NUM, the total number of iterations.
//
//    Output, int &ROT_NUM, the total number of rotations.
//
{
  double *bw;
  double c;
  double g;
  double gapq;
  double h;
  int i;
  int j;
  int k;
  int l;
  int m;
  int p;
  int q;
  double s;
  double t;
  double tau;
  double term;
  double termp;
  double termq;
  double theta;
  double thresh;
  double w;
  double *zw;

  r8mat_identity ( n, v );

  r8mat_diag_get_vector ( n, a, d );

  bw = new double[n];
  zw = new double[n];

  for ( i = 0; i < n; i++ )
  {
    bw[i] = d[i];
    zw[i] = 0.0;
  }
  it_num = 0;
  rot_num = 0;

  while ( it_num < it_max )
  {
    it_num = it_num + 1;
//
//  The convergence threshold is based on the size of the elements in
//  the strict upper triangle of the matrix.
//
    thresh = 0.0;
    for ( j = 0; j < n; j++ )
    {
      for ( i = 0; i < j; i++ )
      {
        thresh = thresh + a[i+j*n] * a[i+j*n];
      }
    }

    thresh = sqrt ( thresh ) / ( double ) ( 4 * n );

    if ( thresh == 0.0 )
    {
      break;
    }

    for ( p = 0; p < n; p++ )
    {
      for ( q = p + 1; q < n; q++ )
      {
        gapq = 10.0 * fabs ( a[p+q*n] );
        termp = gapq + fabs ( d[p] );
        termq = gapq + fabs ( d[q] );
//
//  Annihilate tiny offdiagonal elements.
//
        if ( 4 < it_num &&
             termp == fabs ( d[p] ) &&
             termq == fabs ( d[q] ) )
        {
          a[p+q*n] = 0.0;
        }
//
//  Otherwise, apply a rotation.
//
        else if ( thresh <= fabs ( a[p+q*n] ) )
        {
          h = d[q] - d[p];
          term = fabs ( h ) + gapq;

          if ( term == fabs ( h ) )
          {
            t = a[p+q*n] / h;
          }
          else
          {
            theta = 0.5 * h / a[p+q*n];
            t = 1.0 / ( fabs ( theta ) + sqrt ( 1.0 + theta * theta ) );
            if ( theta < 0.0 )
            {
              t = - t;
            }
          }
          c = 1.0 / sqrt ( 1.0 + t * t );
          s = t * c;
          tau = s / ( 1.0 + c );
          h = t * a[p+q*n];
//
//  Accumulate corrections to diagonal elements.
//
          zw[p] = zw[p] - h;
          zw[q] = zw[q] + h;
          d[p] = d[p] - h;
          d[q] = d[q] + h;

          a[p+q*n] = 0.0;
//
//  Rotate, using information from the upper triangle of A only.
//
          for ( j = 0; j < p; j++ )
          {
            g = a[j+p*n];
            h = a[j+q*n];
            a[j+p*n] = g - s * ( h + g * tau );
            a[j+q*n] = h + s * ( g - h * tau );
          }

          for ( j = p + 1; j < q; j++ )
          {
            g = a[p+j*n];
            h = a[j+q*n];
            a[p+j*n] = g - s * ( h + g * tau );
            a[j+q*n] = h + s * ( g - h * tau );
          }

          for ( j = q + 1; j < n; j++ )
          {
            g = a[p+j*n];
            h = a[q+j*n];
            a[p+j*n] = g - s * ( h + g * tau );
            a[q+j*n] = h + s * ( g - h * tau );
          }
//
//  Accumulate information in the eigenvector matrix.
//
          for ( j = 0; j < n; j++ )
          {
            g = v[j+p*n];
            h = v[j+q*n];
            v[j+p*n] = g - s * ( h + g * tau );
            v[j+q*n] = h + s * ( g - h * tau );
          }
          rot_num = rot_num + 1;
        }
      }
    }

    for ( i = 0; i < n; i++ )
    {
      bw[i] = bw[i] + zw[i];
      d[i] = bw[i];
      zw[i] = 0.0;
    }
  }
//
//  Restore upper triangle of input matrix.
//
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < j; i++ )
    {
      a[i+j*n] = a[j+i*n];
    }
  }
//
//  Ascending sort the eigenvalues and eigenvectors.
//
  for ( k = 0; k < n - 1; k++ )
  {
    m = k;
    for ( l = k + 1; l < n; l++ )
    {
      if ( d[l] < d[m] )
      {
        m = l;
      }
    }

    if ( m != k )
    {
      t    = d[m];
      d[m] = d[k];
      d[k] = t;
      for ( i = 0; i < n; i++ )
      {
        w        = v[i+m*n];
        v[i+m*n] = v[i+k*n];
        v[i+k*n] = w;
      }
    }
  }

  delete [] bw;
  delete [] zw;

  return;
}

void AnisotropicDamageMaterial1 ::r8mat_diag_get_vector ( int n, double a[], double v[] )
//
//  Purpose:
//
//    R8MAT_DIAG_GET_VECTOR gets the value of the diagonal of an R8MAT.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//
//  Parameters:
//
//    Input, int N, the number of rows and columns of the matrix.
//
//    Input, double A[N*N], the N by N matrix.
//
//    Output, double V[N], the diagonal entries
//    of the matrix.
//
{
  int i;

  for ( i = 0; i < n; i++ )
  {
    v[i] = a[i+i*n];
  }

  return;
}

void AnisotropicDamageMaterial1 ::r8mat_identity ( int n, double a[] )
//
//  Purpose:
//
//    R8MAT_IDENTITY sets the square matrix A to the identity.
//
//  Discussion:
//
//    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
//    in column-major order.
//
//  Parameters:
//
//    Input, int N, the order of A.
//
//    Output, double A[N*N], the N by N identity matrix.
//
{
  int i;
  int j;
  int k;

  k = 0;
  for ( j = 0; j < n; j++ )
  {
    for ( i = 0; i < n; i++ )
    {
      if ( i == j )
      {
        a[k] = 1.0;
      }
      else
      {
        a[k] = 0.0;
      }
      k = k + 1;
    }
  }

  return;
}
void AnisotropicDamageMaterial1 :: give3dMaterialStiffnessMatrix(FloatMatrix &answer,MatResponseMode mmode,
                                                                 GaussPoint *gp, TimeStep *tStep)
{
}
*/

void AnisotropicDamageMaterial1 :: computeSecantOperator(FloatMatrix &answer, double MATD[3][3][3][3],  const double &Dc, const FloatMatrix &damageTensor)
// Implementation of the 3D stiffness matrix, according to the equations 56 and 57 of the reference paper.
{
    //    Provide all the varible needed, no Gauss point value are accessed directly. 
    double B1, B2, H, A1, A2, A3, A4 ;
    double traceD ;
    FloatMatrix IE, MATS_2, Omega3D ;

    B1 = ( 1.0 + this->nu )/ this->E / 2.0 ;
    B2 = this->nu / this->E ;
    H  = 16 * ( 1 - pow( this->nu, 2 ) ) / ( 3 * this->E * ( 2 - this->nu ) );
    A1 = -Dc * H / 70.0 ;
    A2 = ( 7 + 2 * Dc )* H / 7.0 ;
    A3 = Dc * H / 7.0 ;
    A4 = -Dc * H / 35.0 ;

    traceD = damageTensor.at(1,1)+damageTensor.at(2,2);
    Omega3D.resize(3,3);
    Omega3D.zero();
    Omega3D.at(1,1)=damageTensor.at(1,1);
    Omega3D.at(1,2)=damageTensor.at(1,2);
    Omega3D.at(2,1)=damageTensor.at(2,1);
    Omega3D.at(2,2)=damageTensor.at(2,2);

    IE.resize(3, 3);
    IE.zero();
    IE.at(1, 1) = IE.at(2, 2) = IE.at(3, 3) = 1.0;

    double MATS[3][3][3][3]={0};
    for ( int i = 1; i <= 3; i++ ) {
        for ( int j = 1; j <= 3; j++ ) {
            for ( int k = 1; k <= 3; k++ ) {
                for ( int l = 1; l <= 3; l++  ) {
                    MATS[i-1][j-1][k-1][l-1] = B1*(IE.at(i,k)*IE.at(j,l)+IE.at(i,l)*IE.at(j,k))-B2*IE.at(i,j)*IE.at(k,l)+2.0*A1*traceD*IE.at(i,j)*IE.at(k,l)\
                                                + 0.5*A2*(IE.at(i,k)*Omega3D.at(j,l)+IE.at(i,l)*Omega3D.at(j,k)+Omega3D.at(i,k)*IE.at(j,l)+Omega3D.at(i,l)*IE.at(j,k))\
                                                + A3*(IE.at(i,j)*Omega3D.at(k,l)+Omega3D.at(i,j)*IE.at(k,l))+A4*traceD*(IE.at(i,k)*IE.at(j,l)+IE.at(i,l)*IE.at(j,k));
                }
            }
        }
    }
    MATS_2.resize(6,6);
    MATS_2.zero();
    MAT4_MAT2( MATS, MATS_2, 2 );
    MATS_2.beInverseOf( MATS_2 );

    // The resulting material stiffness matrix is built
    answer.resize(6, 6);
    answer.at(1, 1) = MATS_2.at(1, 1);
    answer.at(1, 2) = MATS_2.at(1, 2);
    answer.at(1, 3) = MATS_2.at(1, 3);
    answer.at(1, 4) = MATS_2.at(1, 4);
    answer.at(1, 5) = MATS_2.at(1, 5);
    answer.at(1, 6) = MATS_2.at(1, 6);
    answer.at(2, 1) = MATS_2.at(2, 1);
    answer.at(2, 2) = MATS_2.at(2, 2);
    answer.at(2, 3) = MATS_2.at(2, 3);
    answer.at(2, 4) = MATS_2.at(2, 4);
    answer.at(2, 5) = MATS_2.at(2, 5);
    answer.at(2, 6) = MATS_2.at(2, 6);
    answer.at(3, 1) = MATS_2.at(3, 1);
    answer.at(3, 2) = MATS_2.at(3, 2);
    answer.at(3, 3) = MATS_2.at(3, 3);
    answer.at(3, 4) = MATS_2.at(3, 4);
    answer.at(3, 5) = MATS_2.at(3, 5);
    answer.at(3, 6) = MATS_2.at(3, 6);
    answer.at(4, 1) = MATS_2.at(4, 1);
    answer.at(4, 2) = MATS_2.at(4, 2);
    answer.at(4, 3) = MATS_2.at(4, 3);
    answer.at(4, 4) = MATS_2.at(4, 4);
    answer.at(4, 5) = MATS_2.at(4, 5);
    answer.at(4, 6) = MATS_2.at(4, 6);
    answer.at(5, 1) = MATS_2.at(5, 1);
    answer.at(5, 2) = MATS_2.at(5, 2);
    answer.at(5, 3) = MATS_2.at(5, 3);
    answer.at(5, 4) = MATS_2.at(5, 4);
    answer.at(5, 5) = MATS_2.at(5, 5);
    answer.at(5, 6) = MATS_2.at(5, 6);
    answer.at(6, 1) = MATS_2.at(6, 1);
    answer.at(6, 2) = MATS_2.at(6, 2);
    answer.at(6, 3) = MATS_2.at(6, 3);
    answer.at(6, 4) = MATS_2.at(6, 4);
    answer.at(6, 5) = MATS_2.at(6, 5);
    answer.at(6, 6) = MATS_2.at(6, 6);
    MAT2_MAT4( answer, MATD, 1);

    return;
}



void AnisotropicDamageMaterial1 :: MAT2_MAT4(const FloatMatrix &DMATRIX, double TENSOR[3][3][3][3], int ICOE)
{
    int COE1, COE2;
    //double TENSOR[3][3][3][3]={0};
    if (ICOE == 1){
         COE1=1.;
         COE2=1.;
    }else if (ICOE == 2){
         COE1=2.;
         COE2=4.;
    }      
      TENSOR[0][0][0][0] = DMATRIX.at(1,1);
      TENSOR[0][0][1][1] = DMATRIX.at(1,2);
      TENSOR[0][0][2][2] = DMATRIX.at(1,3);
      TENSOR[0][0][0][1] = DMATRIX.at(1,4)/COE1;
      TENSOR[0][0][1][0] = DMATRIX.at(1,4)/COE1;
      TENSOR[0][0][1][2] = DMATRIX.at(1,5)/COE1;
      TENSOR[0][0][2][1] = DMATRIX.at(1,5)/COE1;
      TENSOR[0][0][0][2] = DMATRIX.at(1,6)/COE1;
      TENSOR[0][0][2][0] = DMATRIX.at(1,6)/COE1;

      TENSOR[1][1][0][0] = DMATRIX.at(2,1);
      TENSOR[1][1][1][1] = DMATRIX.at(2,2);
      TENSOR[1][1][2][2] = DMATRIX.at(2,3);
      TENSOR[1][1][0][1] = DMATRIX.at(2,4)/COE1;
      TENSOR[1][1][1][0] = DMATRIX.at(2,4)/COE1;
      TENSOR[1][1][1][2] = DMATRIX.at(2,5)/COE1;
      TENSOR[1][1][2][1] = DMATRIX.at(2,5)/COE1;
      TENSOR[1][1][0][2] = DMATRIX.at(2,6)/COE1;
      TENSOR[1][1][2][0] = DMATRIX.at(2,6)/COE1;

      TENSOR[2][2][0][0] = DMATRIX.at(3,1);
      TENSOR[2][2][1][1] = DMATRIX.at(3,2);
      TENSOR[2][2][2][2] = DMATRIX.at(3,3);
      TENSOR[2][2][0][1] = DMATRIX.at(3,4)/COE1;
      TENSOR[2][2][1][0] = DMATRIX.at(3,4)/COE1;
      TENSOR[2][2][1][2] = DMATRIX.at(3,5)/COE1;
      TENSOR[2][2][2][1] = DMATRIX.at(3,5)/COE1;
      TENSOR[2][2][0][2] = DMATRIX.at(3,6)/COE1;
      TENSOR[2][2][2][0] = DMATRIX.at(3,6)/COE1;

      TENSOR[0][1][0][0] = DMATRIX.at(4,1)/COE1;
      TENSOR[0][1][1][1] = DMATRIX.at(4,2)/COE1;
      TENSOR[0][1][2][2] = DMATRIX.at(4,3)/COE1;
      TENSOR[0][1][0][1] = DMATRIX.at(4,4)/COE2;
      TENSOR[0][1][1][0] = DMATRIX.at(4,4)/COE2;
      TENSOR[0][1][1][2] = DMATRIX.at(4,5)/COE2;
      TENSOR[0][1][2][1] = DMATRIX.at(4,5)/COE2;
      TENSOR[0][1][0][2] = DMATRIX.at(4,6)/COE2;
      TENSOR[0][1][2][0] = DMATRIX.at(4,6)/COE2;

      TENSOR[1][2][0][0] = DMATRIX.at(5,1)/COE1;
      TENSOR[1][2][1][1] = DMATRIX.at(5,2)/COE1;
      TENSOR[1][2][2][2] = DMATRIX.at(5,3)/COE1;
      TENSOR[1][2][0][1] = DMATRIX.at(5,4)/COE2;
      TENSOR[1][2][1][0] = DMATRIX.at(5,4)/COE2;
      TENSOR[1][2][1][2] = DMATRIX.at(5,5)/COE2;
      TENSOR[1][2][2][1] = DMATRIX.at(5,5)/COE2;
      TENSOR[1][2][0][2] = DMATRIX.at(5,6)/COE2;
      TENSOR[1][2][2][0] = DMATRIX.at(5,6)/COE2;

      TENSOR[0][2][0][0] = DMATRIX.at(6,1)/COE1;
      TENSOR[0][2][1][1] = DMATRIX.at(6,2)/COE1;
      TENSOR[0][2][2][2] = DMATRIX.at(6,3)/COE1;
      TENSOR[0][2][0][1] = DMATRIX.at(6,4)/COE2;
      TENSOR[0][2][1][0] = DMATRIX.at(6,4)/COE2;
      TENSOR[0][2][1][2] = DMATRIX.at(6,5)/COE2;
      TENSOR[0][2][2][1] = DMATRIX.at(6,5)/COE2;
      TENSOR[0][2][0][2] = DMATRIX.at(6,6)/COE2;
      TENSOR[0][2][2][0] = DMATRIX.at(6,6)/COE2;
     
      TENSOR[1][0][0][0] = DMATRIX.at(4,1)/COE1;
      TENSOR[1][0][1][1] = DMATRIX.at(4,2)/COE1;
      TENSOR[1][0][2][2] = DMATRIX.at(4,3)/COE1;
      TENSOR[1][0][0][1] = DMATRIX.at(4,4)/COE2;
      TENSOR[1][0][1][0] = DMATRIX.at(4,4)/COE2;
      TENSOR[1][0][1][2] = DMATRIX.at(4,5)/COE2;
      TENSOR[1][0][2][1] = DMATRIX.at(4,5)/COE2;
      TENSOR[1][0][0][2] = DMATRIX.at(4,6)/COE2;
      TENSOR[1][0][2][0] = DMATRIX.at(4,6)/COE2;

      TENSOR[2][1][0][0] = DMATRIX.at(5,1)/COE1;
      TENSOR[2][1][1][1] = DMATRIX.at(5,2)/COE1;
      TENSOR[2][1][2][2] = DMATRIX.at(5,3)/COE1;
      TENSOR[2][1][0][1] = DMATRIX.at(5,4)/COE2;
      TENSOR[2][1][1][0] = DMATRIX.at(5,4)/COE2;
      TENSOR[2][1][1][2] = DMATRIX.at(5,5)/COE2;
      TENSOR[2][1][2][1] = DMATRIX.at(5,5)/COE2;
      TENSOR[2][1][0][2] = DMATRIX.at(5,6)/COE2;
      TENSOR[2][1][2][0] = DMATRIX.at(5,6)/COE2;

      TENSOR[2][0][0][0] = DMATRIX.at(6,1)/COE1;
      TENSOR[2][0][1][1] = DMATRIX.at(6,2)/COE1;
      TENSOR[2][0][2][2] = DMATRIX.at(6,3)/COE1;
      TENSOR[2][0][0][1] = DMATRIX.at(6,4)/COE2;
      TENSOR[2][0][1][0] = DMATRIX.at(6,4)/COE2;
      TENSOR[2][0][1][2] = DMATRIX.at(6,5)/COE2;
      TENSOR[2][0][2][1] = DMATRIX.at(6,5)/COE2;
      TENSOR[2][0][0][2] = DMATRIX.at(6,6)/COE2;
      TENSOR[2][0][2][0] = DMATRIX.at(6,6)/COE2;
      return;
}

void AnisotropicDamageMaterial1 :: MAT4_MAT2(const double TENSOR[3][3][3][3], FloatMatrix &DMATRIX, int ICOE)
{
    int COE1, COE2;
    DMATRIX.resize(6,6);
    DMATRIX.zero(); 
    if (ICOE == 1){
         COE1=1.;
         COE2=1.;
    }else if (ICOE == 2){
         COE1=2.;
         COE2=4.;
    } 
      DMATRIX.at(1,1)=TENSOR[0][0][0][0];
      DMATRIX.at(1,2)=TENSOR[0][0][1][1];
      DMATRIX.at(1,3)=TENSOR[0][0][2][2];
      DMATRIX.at(1,4)=TENSOR[0][0][0][1]*COE1;
      DMATRIX.at(1,5)=TENSOR[0][0][1][2]*COE1;
      DMATRIX.at(1,6)=TENSOR[0][0][0][2]*COE1;

      DMATRIX.at(2,1)=TENSOR[1][1][0][0];
      DMATRIX.at(2,2)=TENSOR[1][1][1][1];
      DMATRIX.at(2,3)=TENSOR[1][1][2][2];
      DMATRIX.at(2,4)=TENSOR[1][1][0][1]*COE1;
      DMATRIX.at(2,5)=TENSOR[1][1][1][2]*COE1;
      DMATRIX.at(2,6)=TENSOR[1][1][0][2]*COE1;

      DMATRIX.at(3,1)=TENSOR[2][2][0][0];
      DMATRIX.at(3,2)=TENSOR[2][2][1][1];
      DMATRIX.at(3,3)=TENSOR[2][2][2][2];
      DMATRIX.at(3,4)=TENSOR[2][2][0][1]*COE1;
      DMATRIX.at(3,5)=TENSOR[2][2][1][2]*COE1;
      DMATRIX.at(3,6)=TENSOR[2][2][0][2]*COE1;

     // here, engineering shear strain is used

      DMATRIX.at(4,1)=TENSOR[0][1][0][0]*COE1;
      DMATRIX.at(4,2)=TENSOR[0][1][1][1]*COE1;
      DMATRIX.at(4,3)=TENSOR[0][1][2][2]*COE1;
      DMATRIX.at(4,4)=TENSOR[0][1][0][1]*COE2;
      DMATRIX.at(4,5)=TENSOR[0][1][1][2]*COE2;
      DMATRIX.at(4,6)=TENSOR[0][1][0][2]*COE2;

      DMATRIX.at(5,1)=TENSOR[1][2][0][0]*COE1;
      DMATRIX.at(5,2)=TENSOR[1][2][1][1]*COE1;
      DMATRIX.at(5,3)=TENSOR[1][2][2][2]*COE1;
      DMATRIX.at(5,4)=TENSOR[1][2][0][1]*COE2;
      DMATRIX.at(5,5)=TENSOR[1][2][1][2]*COE2;
      DMATRIX.at(5,6)=TENSOR[1][2][0][2]*COE2;

      DMATRIX.at(6,1)=TENSOR[0][2][0][0]*COE1;
      DMATRIX.at(6,2)=TENSOR[0][2][1][1]*COE1;
      DMATRIX.at(6,3)=TENSOR[0][2][2][2]*COE1;
      DMATRIX.at(6,4)=TENSOR[0][2][0][1]*COE2;
      DMATRIX.at(6,5)=TENSOR[0][2][1][2]*COE2;
      DMATRIX.at(6,6)=TENSOR[0][2][0][2]*COE2;
      return;
}



void AnisotropicDamageMaterial1 :: givePlaneStrainStiffMtrx(FloatMatrix &answer, MatResponseMode mode,
                                                           GaussPoint *gp, TimeStep *atTime)
{
    AnisotropicDamageMaterial1Status *status = static_cast< AnisotropicDamageMaterial1Status * >( this->giveStatus(gp) );
    if ( mode == ElasticStiffness ) {
        this->giveStiffnessMatrix(answer, mode, gp, atTime);
        return;
    } else {
        FloatMatrix damageTensor;
        // The damage tensor is read
        damageTensor.resize(2, 2);
        damageTensor.zero();
        damageTensor = status->giveTempDamage();

        FloatArray totalStrain = status->giveTempStrainVector();
        FloatArray reducedTotalStrainVector;
        this->giveStressDependentPartOfStrainVector(reducedTotalStrainVector, gp, totalStrain, atTime, VM_Total);
        double treps, Dc;
        treps = reducedTotalStrainVector.at(1)+reducedTotalStrainVector.at(2) + reducedTotalStrainVector.at(3);
        if ( treps > 0 ){
            Dc = - this->nu;
        } else {
            Dc = - 2;
        }

        FloatMatrix MATD_2;
        MATD_2.resize(6,6);
        MATD_2.zero();
        double MATD[3][3][3][3]={0};

        this->computeSecantOperator( MATD_2, MATD, Dc, damageTensor);

        answer.resize(4,4);
        answer.at(1,1)=MATD_2.at(1,1);
        answer.at(1,2)=MATD_2.at(1,2);
        answer.at(1,3)=MATD_2.at(1,3);
        answer.at(1,4)=MATD_2.at(1,4);
        answer.at(2,1)=MATD_2.at(2,1);
        answer.at(2,2)=MATD_2.at(2,2);
        answer.at(2,3)=MATD_2.at(2,3);
        answer.at(2,4)=MATD_2.at(2,4);
        answer.at(3,1)=MATD_2.at(3,1);
        answer.at(3,2)=MATD_2.at(3,2);
        answer.at(3,3)=MATD_2.at(3,3);
        answer.at(3,4)=MATD_2.at(3,4);
        answer.at(4,1)=MATD_2.at(4,1);
        answer.at(4,2)=MATD_2.at(4,2);
        answer.at(4,3)=MATD_2.at(4,3);
        answer.at(4,4)=MATD_2.at(4,4);
        return;
    }
}

IRResultType
AnisotropicDamageMaterial1 :: initializeFrom(InputRecord *ir)
{
    IRResultType result;                // Required by IR_GIVE_FIELD macro

    IR_GIVE_FIELD(ir, kappa_t, _IFT_AnisotropicDamageMaterial1_kappat);
    IR_GIVE_FIELD(ir, kappa_c, _IFT_AnisotropicDamageMaterial1_kappac);
    IR_GIVE_FIELD(ir, alpha_t, _IFT_AnisotropicDamageMaterial1_alphat);
    IR_GIVE_FIELD(ir, alpha_c, _IFT_AnisotropicDamageMaterial1_alphac);
    IR_GIVE_FIELD(ir, eta, _IFT_AnisotropicDamageMaterial1_eta);

    return IsotropicLinearElasticMaterial::initializeFrom(ir);
}

void
AnisotropicDamageMaterial1 :: giveInputRecord(DynamicInputRecord &input)
{
    IsotropicLinearElasticMaterial :: giveInputRecord(input);
    input.setField(this->kappa_t, _IFT_AnisotropicDamageMaterial1_kappat);
    input.setField(this->alpha_t, _IFT_AnisotropicDamageMaterial1_alphat);
    input.setField(this->kappa_c, _IFT_AnisotropicDamageMaterial1_kappac);
    input.setField(this->alpha_c, _IFT_AnisotropicDamageMaterial1_alphac);
    input.setField(this->eta, _IFT_AnisotropicDamageMaterial1_eta);

}

int
AnisotropicDamageMaterial1 :: giveIPValue(FloatArray &answer, GaussPoint *gp, InternalStateType type, TimeStep *atTime)
{
    AnisotropicDamageMaterial1Status *status = static_cast< AnisotropicDamageMaterial1Status * >( this->giveStatus(gp) );
    if ( type == IST_DamageScalar ) { // returning the trace of the damage tensor
        answer.resize(1);
        answer.at(1) = status->giveDamage().at(1, 1) + status->giveDamage().at(2, 2) ;
        return 1;
    } else if ( type == IST_DamageTensor ) {
        answer.resize(6);
        answer.zero();
        answer.at(1) = status->giveDamage().at(1, 1);
        answer.at(2) = status->giveDamage().at(2, 2);
        answer.at(3) = 0.0;
        answer.at(4) = 0.0;
        answer.at(5) = 0.0;
        answer.at(6) = status->giveDamage().at(1, 2);
        return 1;
    } else if ( type == IST_Permeability ) {
        answer.resize(6);
        answer.zero();
        answer.at(1) = status->givePermeability().at(1, 1);
        answer.at(2) = status->givePermeability().at(2, 2);
        answer.at(3) = 0.0;
        answer.at(4) = 0.0;
        answer.at(5) = 0.0;
        answer.at(6) = status->givePermeability().at(1, 2);
        return 1;
    } else if ( type == IST_PrincipalDamageTensor ) {
        //int checker12=this->checkSymmetry(status->giveDamage());
        FloatMatrix dam = status->giveDamage();
        FloatMatrix eVecs;
        dam.jaco_(answer, eVecs, 20);
        return 1;
    } else if ( type == IST_DamageTensorTemp ) {
        answer.resize(6);
        answer.zero();
        answer.at(1) = status->giveTempDamage().at(1, 1);
        answer.at(2) = status->giveTempDamage().at(2, 2);
        answer.at(3) = 0.0;
        answer.at(4) = 0.0;
        answer.at(5) = 0.0;
        answer.at(6) = status->giveTempDamage().at(1, 2);
        return 1;
    } else if ( type == IST_PrincipalDamageTempTensor ) {
        //int checker13=this->checkSymmetry(status->giveTempDamage());
        FloatMatrix dam = status->giveTempDamage();
        FloatMatrix eVecs;
        dam.jaco_(answer, eVecs, 20);
        return 1;
    } else if ( type == IST_EquivalentStrainC ) {
        answer.resize(1);
        answer.at(1) = status->giveEquivalenstrain_c();
        return 1;
    } else if ( type == IST_EquivalentStrainT ) {
        answer.resize(1);
        answer.at(1) = status->giveEquivalenstrain_t();
        return 1;
#ifdef keep_track_of_dissipated_energy
    } else if ( type == IST_StressWorkDensity ) {
        answer.resize(1);
        answer.at(1) = status->giveStressWork();
        return 1;
    } else if ( type == IST_DissWorkDensity ) {
        answer.resize(1);
        answer.at(1) = status->giveDissWork();
        return 1;
    } else if ( type == IST_FreeEnergyDensity ) {
        answer.resize(1);
        answer.at(1) = status->giveStressWork() - status->giveDissWork();
        return 1;

#endif
    } else {
        return StructuralMaterial :: giveIPValue(answer, gp, type, atTime);
    }

    return 1; // to make the compiler happy
}


/*/\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//\\\\\\\Material status functions\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\
//\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\\*/


AnisotropicDamageMaterial1Status :: AnisotropicDamageMaterial1Status(int n, Domain *d, GaussPoint *g) : StructuralMaterialStatus(n, d, g),
    RandomMaterialStatusExtensionInterface()
{
    equivalenstrain_t = tempequivalenstrain_t = 0.0;
    equivalenstrain_c = tempequivalenstrain_c = 0.0;

    damage.resize(2, 2);
    damage.zero();
    tempDamage.resize(2, 2);
    tempDamage.zero();
    permeability.resize(2, 2);
    permeability.at(1,1)=1e-20;
    permeability.at(2,2)=1e-20;
    flag = tempFlag = 0;

#ifdef keep_track_of_dissipated_energy
    stressWork = tempStressWork = 0.0;
    dissWork = tempDissWork = 0.0;
#endif
}

AnisotropicDamageMaterial1Status :: ~AnisotropicDamageMaterial1Status()
{ }


void
AnisotropicDamageMaterial1Status :: initTempStatus()
{
    StructuralMaterialStatus :: initTempStatus();
    this->tempequivalenstrain_t = this->equivalenstrain_t;
    this->tempequivalenstrain_c = this->equivalenstrain_c;
    this->tempDamage = this->damage;
    this->tempFlag = flag;

#ifdef keep_track_of_dissipated_energy
    this->tempStressWork = this->stressWork;
    this->tempDissWork = this->dissWork;
#endif
}

void
AnisotropicDamageMaterial1Status :: updateYourself(TimeStep *atTime)
{
    StructuralMaterialStatus :: updateYourself(atTime);
    this->equivalenstrain_t = this->tempequivalenstrain_t;
    this->equivalenstrain_c = this->tempequivalenstrain_c;
    this->damage = this->tempDamage;

    double damage1, damage2, ceps, seps;

    this->computePrincValDir2D(damage1, damage2, ceps, seps,this->damage.at(1,1), this->damage.at(2,2), this->damage.at(1,2));

    FloatMatrix EigenVic(2,2), deltaperm(2,2);
    deltaperm.zero();
    EigenVic.at(1,1)=ceps;
    EigenVic.at(1,2)=seps;
    EigenVic.at(2,1)=-seps;
    EigenVic.at(2,2)=ceps;

    double V=1;//m
    int N=2e6;
    double chi=0.005;//
    //double mu=0.001002; //N/m^2.s
    //double gamma = 9800; //N/m^3
    double t=1e-5;
    double pi=3.1415926;
    double coefficient=t*V*chi/12/pi/N;
    deltaperm.at(1,1)= coefficient*(pow(damage1,5/3)*(1-EigenVic.at(1,1)*EigenVic.at(1,1))+pow(damage2,5/3)*(1-EigenVic.at(2,1)*EigenVic.at(2,1)));
    deltaperm.at(2,2)= coefficient*(pow(damage1,5/3)*(1-EigenVic.at(1,2)*EigenVic.at(1,2))+pow(damage2,5/3)*(1-EigenVic.at(2,2)*EigenVic.at(2,2)));

    FloatMatrix EigenVicInv = EigenVic;
    for (int i=1; i<=2; i++){
        EigenVic.at(1,i)=EigenVic.at(1,i)*deltaperm.at(i,i);
        EigenVic.at(2,i)=EigenVic.at(2,i)*deltaperm.at(i,i);
    }
    EigenVicInv.beInverseOf( EigenVicInv );
    deltaperm.beProductOf( EigenVic, EigenVicInv );
    permeability.add(deltaperm);

#ifdef keep_track_of_dissipated_energy
    this->stressWork = this->tempStressWork;
    this->dissWork = this->tempDissWork;
#endif
}

Interface *
AnisotropicDamageMaterial1Status :: giveInterface(InterfaceType type)
{
    if ( type == RandomMaterialStatusExtensionInterfaceType ) {
        return static_cast< RandomMaterialStatusExtensionInterface * >(this);
    } else {
        return NULL;
    }
}

void
AnisotropicDamageMaterial1Status :: printOutputAt(FILE *file, TimeStep *tStep)
{
    MaterialMode mode = gp->giveMaterialMode();
    double Volume = gp->givegausspointvolume();
    if ( mode == _PlaneStrain ) { // special treatment of the out-of-plane stress
        FloatArray helpVec;
        MaterialStatus :: printOutputAt(file, tStep);
        fprintf(file, "  strains ");
        StructuralMaterial :: giveFullSymVectorForm(helpVec, strainVector, mode);
        for ( auto &v : helpVec ) {
            fprintf( file, " %.4e", v );
        }
        fprintf(file, "\n              stresses");
        StructuralMaterial :: giveFullSymVectorForm(helpVec, stressVector, mode);
        for ( auto &v : helpVec ) {
            fprintf( file, " %.4e", v );
        }
        fprintf(file, "\n");
    } else {
        StructuralMaterialStatus :: printOutputAt(file, tStep); // standard treatment of strains and stresses
    }

    fprintf(file, "              Damage");
    fprintf( file, " %.4e %.4e %.4e",damage.at(1,1),damage.at(2,2),damage.at(1,2));
#ifdef keep_track_of_dissipated_energy
    fprintf(file, "\n              dissW %f freeE %f stressW %f\n", this->dissWork*Volume, ( this->stressWork - this->dissWork )*Volume, this->stressWork*Volume);
#endif
}

#ifdef keep_track_of_dissipated_energy
void
AnisotropicDamageMaterial1Status :: computeWork(GaussPoint *gp)
{

    // strain increment
     FloatArray deps;
     deps.beDifferenceOf(tempStrainVector, strainVector);
     
     // increment of stress work density
     double dSW = ( tempStressVector.dotProduct(deps) + stressVector.dotProduct(deps) ) / 2.;
     tempStressWork = stressWork + dSW;
     
     // elastically stored energy density
     double We = tempStressVector.dotProduct(tempStrainVector) / 2.;
     
     // dissipative work density
     tempDissWork = tempStressWork - We;
     
}
#endif


void AnisotropicDamageMaterial1Status ::giveInternalStatetype(IntArray &G){

    G.resize(7);
    G.at(1)=IST_StressTensor;
    G.at(2)=IST_StrainTensor;
    G.at(3)=IST_DamageTensor;
    G.at(4)=IST_StressWorkDensity;
    G.at(5)=IST_DissWorkDensity;
    G.at(6)=IST_EquivalentStrainC;
    G.at(7)=IST_EquivalentStrainT;
    //G.at(8)=IST_NonlocalEquivalentStrain;
}



void AnisotropicDamageMaterial1Status :: copyStateVariables(const MaterialStatus &iStatus)
{
    MaterialStatus &tmpStat = const_cast< MaterialStatus & >(iStatus);
    const AnisotropicDamageMaterial1Status &adm1Status = dynamic_cast< AnisotropicDamageMaterial1Status & >(tmpStat);

    equivalenstrain_t = adm1Status.giveEquivalenstrain_t();
    tempequivalenstrain_t = adm1Status.giveTempEquivalenstrain_t();
    equivalenstrain_c = adm1Status.giveEquivalenstrain_c();
    tempequivalenstrain_c = adm1Status.giveTempEquivalenstrain_c();
    damage = adm1Status.giveDamage();
    tempDamage = adm1Status.giveTempDamage();
    flag = adm1Status.giveFlag();
    tempFlag = adm1Status.giveTempFlag();

    permeability = adm1Status.givePermeability();

#ifdef keep_track_of_dissipated_energy
    stressWork = adm1Status.giveStressWork();
    tempStressWork = adm1Status.giveTempStressWork();
    dissWork = adm1Status.giveDissWork();
    tempDissWork = adm1Status.giveTempDissWork();
#endif

    StructuralMaterialStatus :: copyStateVariables (iStatus);
}


void
 AnisotropicDamageMaterial1Status :: computePrincValDir2D(double &D1, double &D2, double &c, double &s, double Dx, double Dy, double Dxy)
//
// computes the principal values and directions of a symmetric second-order tensor in 2D
// input: Dx, Dy, Dxy ... components of the tensor wrt global coordinates
// output: D1, D2 ... ordered principal values, D1>=D2
// output: c, s ... components of the unit principal vector associated with D1
//                  (cosine and sine of the angle between the major principal direction and the global x-axis)
{
    double aux1 = ( Dx + Dy ) / 2.;
    double aux2 = ( Dx - Dy ) / 2.;
    double aux3 = sqrt(aux2 * aux2 + Dxy * Dxy);
    D1 = aux1 + aux3;
    D2 = aux1 - aux3;
    // formulae (90)-(92) and the two cases preceding them
    c = 1.;
    s = 0.;         // cases 1 and 2a
    if ( Dxy != 0. ) {  // Consider equal side triangle, with middle line trough
        double t = ( D1 - Dx ) / Dxy;
        c = 1. / sqrt(1. + t * t);
        s = c * t;
    } else if ( Dx < Dy ) { // case 2b
        c = 0.;
        s = 1.;
    }
    return;
}



} // end namespace oofems
