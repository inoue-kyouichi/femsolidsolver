#ifndef _SignedDistanceFunction_H_
#define _SignedDistanceFunction_H_

//##################################################################################
//
// FEM Base
//
// Copyright (c) 2019 Mechanical and Bioengineering Systems Lab.,
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

// Nagrath et al., Compt. Methods Appl. Mech. Engrg., 194 (2005), 4565-4587

#include "fem.h"

class SignedDistanceFunction{
  public:
  Fem FEM,FEMorg;
  ARRAY1D<double> SDF;
  VECTOR1D<int> elementToSolve;
  VECTOR2D<int> knownNodes;
  VECTOR2D<int> elementsPerNodes;

  void set1stTetra();
  void calcKnownNodes();
  void calcElementsPerNodes();
  void exportVTU(std::string file);

  void calcSDF();

  private:
  ARRAY1D<int> mask;
  double dt;
  double calcSDFofUnknownNode(const int elm,int &unknown);
  double calcEikonalEquation();
  void calcPredictionMatrix(VECTOR1D<ARRAY2D<double>> &LHS,VECTOR1D<ARRAY2D<double>> &RHS,const int ic);
  void Prediction_Galerkin_inGaussIntegral(VECTOR1D<ARRAY2D<double>> &LHS,VECTOR1D<ARRAY2D<double>> &RHS,ARRAY1D<double> &N,ARRAY2D<double> &dNdr,ARRAY2D<double> &x_current,ARRAY2D<double> &dNdx,const int numOfNodeInElm,const double weight,const int ic);
  void Prediction_SUPG_inGaussIntegral(VECTOR1D<ARRAY2D<double>> &LHS,VECTOR1D<ARRAY2D<double>> &RHS,ARRAY1D<double> &N,ARRAY2D<double> &dNdr,ARRAY2D<double> &x_current,ARRAY2D<double> &dNdx,const int numOfNodeInElm,const double weight,const int ic);
  double SUPG_stabilizationParameter(const double (&dxdr)[3][3],const double (&advel)[3],VTKCellType &CellType);

  void calcSourceTerm(ARRAY1D<double> &b,const int ic);
  void calc_dt();
};


#endif