/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */
#include "eye_mechanics.h"
#include <ostream>
#include <fstream>

using namespace std;

// #################################################################
/**
 * @brief tentative
 */
void EYE_Mechanics::EyeMech::calcBoundaryForce()
{
  VECTOR1D<ARRAY2D<double>> Qb(numOfBoundaryElm);
  for(int ic=0;ic<numOfBoundaryElm;ic++) Qb[ic].allocate(boundaryElement[ic].node.size(),3);

  #pragma omp parallel for
  for(int ic=0;ic<numOfBoundaryElm;ic++){
    calcBoundaryPressure_spatialForm(ic,Qb,boundaryPressure);
  }

  #pragma omp parallel for
  for(int ic=0;ic<numOfNode;ic++){
    for(int j=0;j<3;j++)  boundaryForce(ic,j) = 0e0;
  }

  for(int ic=0;ic<numOfBoundaryElm;ic++){
    for(int p=0;p<boundaryElement[ic].node.size();p++){
      for(int i=0;i<3;i++){
        boundaryForce(boundaryElement[ic].node[p],i) += Qb[ic](p,i);
      }
    }
  }

}

// #################################################################
/**
 * @brief tentative
 */
void EYE_Mechanics::EyeMech::calcBoundaryPressure_spatialForm(const int ic,VECTOR1D<ARRAY2D<double>> &Qb,const double boundaryPressure)
{
  int numOfNodeInElm=boundaryElement[ic].node.size();
  ARRAY2D<double> x_current(numOfNodeInElm,3);
  ARRAY2D<double> x_ref(numOfNodeInElm,3);

  ARRAY1D<double> N(numOfNodeInElm);
  ARRAY2D<double> dNdr(numOfNodeInElm,3);

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      x_current(p,i) = x0(boundaryElement[ic].node[p],i)+U(boundaryElement[ic].node[p],i);
      x_ref(p,i)     = x0(boundaryElement[ic].node[p],i);
    }
  }

  Gauss g(1),g2(2);
  GaussTriangle gTri(1),gTri2(2);


  for(int i=0;i<3;i++){
    for(int p=0;p<numOfNodeInElm;p++){
      Qb[ic](p,i) = 0e0;
    }
  }


  switch(boundaryElement[ic].meshType){
    case VTK_TRIANGLE:
      ShapeFunction2D::C2D3_N(N,gTri.point[0][0],gTri.point[0][1],gTri.point[0][2]);
      ShapeFunction2D::C2D3_dNdr(dNdr,gTri.point[0][0],gTri.point[0][1],gTri.point[0][2]);
      // cout << gTri.weight[0]*5e-1 << endl;
      boundaryPressure_inGaussIntegral(Qb,N,dNdr,x_current,x_ref,numOfNodeInElm,boundaryPressure,gTri.weight[0]*5e-1,ic);
      break;
    case VTK_QUAD:
      for(int i1=0;i1<2;i1++){
        for(int i2=0;i2<2;i2++){
            ShapeFunction2D::C2D4_N(N,g.point[i1],g.point[i2]);
            ShapeFunction2D::C2D4_dNdr(dNdr,g.point[i1],g.point[i2]);
            boundaryPressure_inGaussIntegral(Qb,N,dNdr,x_current,x_ref,numOfNodeInElm,boundaryPressure,g.weight[i1]*g.weight[i2],ic);
        }
      }
      break;
    default:
      cout << "undefined mesh type" << endl;
  }

}

// #################################################################
/**
 * @brief tentative
 */
void EYE_Mechanics::EyeMech::boundaryPressure_inGaussIntegral(VECTOR1D<ARRAY2D<double>> &Qb,ARRAY1D<double> &N,ARRAY2D<double> &dNdr,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,const int numOfNodeInElm,const double boundaryPressure,const int weight,const int ic)
{
  double dxdr1[3],dxdr2[3],dxdr3[3],Jacobian;

  for(int i=0;i<3;i++){
    dxdr1[i] = 0e+0;
    dxdr2[i] = 0e+0;
    for(int p=0;p<numOfNodeInElm;p++){
      dxdr1[i] += dNdr(p,0) * x_current(p,i);
      dxdr2[i] += dNdr(p,1) * x_current(p,i);
    }
  }

  mathTool::crossProduct(dxdr1,dxdr2,dxdr3,Jacobian);

  double normal[3];
  for(int j=0;j<3;j++) normal[j] = dxdr3[j] / Jacobian;

  for(int i=0;i<3;i++){
    for(int p=0;p<numOfNodeInElm;p++){
      Qb[ic](p,i) += -1e0 * boundaryPressure *  N(p) * normal[i] * Jacobian * weight;
    }
  }

}