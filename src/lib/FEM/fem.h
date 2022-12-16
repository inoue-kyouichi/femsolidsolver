#ifndef _FEM_H_
#define _FEM_H_

//##################################################################################
//
// FEM Base
//
// Copyright (c) 2020 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   fem.h
 * @brief  FEMBase Definition Header
 * @author T. Otani
 */

#include <iostream>
#include <cmath>
#include <string>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <strings.h>
#include <sys/stat.h>
#include <omp.h>
#include <chrono>

#include "allocation.h"
#include "domain.h"
#include "gauss.h"
#include "math_tools.h"
#include "fem_define.h"
#include "fem_base_mathTool.h"
#include "ShapeFunction.h"
#include "fileIO.h"
#include "pardiso_solver.h"
#include "TextParser.h"

class Fem : public Domain{

  //fem.cpp
 public:
  Fem(){
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) I2[i][j]=0e0;
        I2[i][i]=1e0;
    }

    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        for(int k=0;k<3;k++){
          for(int l=0;l<3;l++){
            I4[i][j][k][l]=5e-1*(I2[i][k]*I2[j][l]+I2[i][l]*I2[j][k]);
          }
        }
      }
    }
  }

  double I2[3][3],I4[3][3][3][3];

  double totalVolume;

  double rho;
  ARRAY1D<double> volume,volume0,volumeChangeRatio;
  ARRAY2D<double> U, innerForce, externalForce, RHS;

  std::vector<ARRAY2D<double>> Qu;
  std::vector<ARRAY2D<double>> Mass;
  std::vector<ARRAY4D<double>> Ku;

  void rotationalDirichlet(const int loop);
  void set_rhs_statics();
  void calc_MassMatrix();
  void corrector_statics(const double *u,const double relaxation);
  void calcVolume_hexa(const int &ic,ARRAY1D<double> &elementVolume,const int &numOfNodeInElm,const int &numOfGaussPoint,const bool option);

  void stress_tensor_initialize();
  void setInnerForce();

  void calcStiffnessMatrix(std::function<void(ARRAY2D<double>&,ARRAY2D<double>&,ARRAY2D<double>&,ARRAY2D<double>&)> func,ARRAY2D<double> &U_tmp);

 private:
  void exportRestartData(const int loop);

  //line search
  double line_search(const double *u);
  double line_search_innerProduct(ARRAY2D<double> &Q,const double *u);
  void calc_Q(ARRAY2D<double> &innerForce_tmp,ARRAY2D<double> &U_tmp);

  //fem_preprocessing.cpp
 public:
  void initialize(TextParser &tp);
  void inputDomainInfo(TextParser &tp);
  void allocate();
 private:
  void inputDirichletBoundaryInfo(TextParser &tp);
  void inputNeumannBoundaryInfo(TextParser &tp);
  void inputFiberInfo(TextParser &tp);
  void restart_setting(const int dataNumber,const bool Restart,TextParser &tp);

  //fem_postprocessing.cpp
  public:
  ARRAY1D<double> Mises,Mises_strain;
  ARRAY2D<double> AEigen_Ave,sigmaEigen_Ave,G_strainEigen_Ave;
  ARRAY3D<double> AEigenVector_Ave,sigmaEigenVector_Ave,G_strainEigenVector_Ave;

  void calcEigen(const double (&A)[3][3],double (&AEigen)[3],double (&AEigenVector)[3][3]);
  void normalize(ARRAY3D<double> &AEigenVector_Ave,const int ic);

  //fem_boundary.cpp
 public:
  double boundaryPressure;
  ARRAY2D<double> externalSurfaceForce;
  std::vector<ARRAY2D<double>> BFe;

  void inputSurfaceBoundary(TextParser &tp);
  void calc_externalSurfaceForce_prescribedTraction(std::vector<ElementType> &element,ARRAY2D<double> &Traction);
 private:
  void calc_Traction_element_quad(const int ic,const int numOfNodeInBdElm,const int numOfGaussPoint,ElementType &belement,const double (&TractionForce)[3]);
  void calc_Traction_element_triangle(const int ic,const int numOfNodeInBdElm,const int numOfGaussPoint,ElementType &belement,const double (&TractionForce)[3]);
  // void calc_TractionByPressure_element(const int &ic,const Element &boundaryElement);
  // void calcBFe_inGaussIntegral(ARRAY3D<double> &BFe,ARRAY2D<double> &dNdr,ARRAY2D<double> &X,const int numOfNodeInBdElm,const double weight,const int ic);
};

#endif //_FEM_H_
