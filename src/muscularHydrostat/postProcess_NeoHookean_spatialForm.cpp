/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "muscularHydrostat.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <algorithm>

#include <Eigen/Core>
#include <Eigen/Eigen>

using namespace std;
using namespace Eigen;

// #################################################################
/**
 * @brief calc stress tensor of SantVenant material
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] option           true or faluse: calculate tangential stiffness matrix or not.
 */
void muscularHydrostat::Muscle::postProcess_element_spatialForm(const int ic,ARRAY2D<double> &U_tmp,const bool option)
{
  //cout <<ic<<endl;
  //cout <<"1"<<endl;
  double nu = 123.0;
  double Poisson = 1999.0/4080.0;
  int numOfNodeInElm=element[ic].node.size();
  ARRAY2D<double> x_current(numOfNodeInElm,3);
  ARRAY2D<double> x_ref(numOfNodeInElm,3);
  ARRAY2D<double> dNdr(numOfNodeInElm,3);

  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      x_current(p,i) = x0(element[ic].node[p],i)+U_tmp(element[ic].node[p],i);
      x_ref(p,i)     = x0(element[ic].node[p],i);
    }
  }

  double stress[3][3];

  Gauss g(1),g2(2);
  GaussTetra gTet(1),gTet2(2);
  GaussTriangle gTri(1),gTri2(2);
  double G_strainEigen[3],G_strainEigenVector[3][3];
//cout <<"2"<<endl;

for(int i=0;i<3;i++){
        G_strainEigen_Ave(ic,i)=0e0;
        for(int j=0;j<3;j++) G_strainEigenVector_Ave(ic,i,j)=0e0;
      }
  switch(element[ic].meshType){
    case VTK_TETRA:
      ShapeFunction3D::C3D4_dNdr(dNdr,gTet.point[0][0],gTet.point[0][1],gTet.point[0][2],gTet.point[0][3]);
      postProcess_inGaussIntegral(ic,numOfNodeInElm,x_current,x_ref,dNdr,gTet.weight[0]*1e0/6e0,stress,option,G_strainEigen,G_strainEigenVector);
      break;
    //case VTK_HEXAHEDRON:
    //  for(int i1=0;i1<2;i1++){
    //    for(int i2=0;i2<2;i2++){
    //      for(int i3=0;i3<2;i3++){
    //        ShapeFunction3D::C3D8_dNdr(dNdr,g.point[i1],g.point[i2],g.point[i3]);
    //        Takenaka2019_inGaussIntegral(ic,numOfNodeInElm,x_current,x_ref,dNdr,g.weight[i1]*g.weight[i2]*g.weight[i3],stress,option);
    //      }
    //    }
    //  }
    //  break;
    //case VTK_QUADRATIC_TETRA:
    //  for(int i1=0;i1<4;i1++){
    //    ShapeFunction3D::C3D10_dNdr(dNdr,gTet2.point[i1][0],gTet2.point[i1][1],gTet2.point[i1][2],gTet2.point[i1][3]);
    //    Takenaka2019_inGaussIntegral(ic,numOfNodeInElm,x_current,x_ref,dNdr,gTet2.weight[i1]*1e0/6e0,stress,option);
    //  }
    //  break;
    //case VTK_QUADRATIC_HEXAHEDRON:
    //  break;
    //case VTK_WEDGE:
    //  for(int i1=0;i1<2;i1++){
    //    ShapeFunction3D::C3D6_dNdr(dNdr,gTri.point[0][0],gTri.point[0][1],gTri.point[0][2],g.point[i1]);
    //    Takenaka2019_inGaussIntegral(ic,numOfNodeInElm,x_current,x_ref,dNdr,5e-1*gTri.weight[0]*g.weight[i1],stress,option);
    //  }
    //  break;
    default:
      cout << "undefined mesh type" << endl;
      exit(1);
  }
  //cout <<"4"<<endl;
  for(int i=0;i<3;i++){
    G_strainEigen_Ave(ic,i)+=G_strainEigen[i];
    for(int j=0;j<3;j++){
      G_strainEigenVector_Ave(ic,i,j)+=G_strainEigenVector[i][j];
    }
  }
  for(int i=0;i<3;i++){
    G_strainEigen_Ave(ic,i)/=4e0;
  }
  Mises_strain(ic)=1/(1+Poisson)*sqrt(5e-1*(pow(G_strainEigen_Ave(ic,0)-G_strainEigen_Ave(ic,1),2e0)+pow(G_strainEigen_Ave(ic,1)-G_strainEigen_Ave(ic,2),2e0)+pow(G_strainEigen_Ave(ic,2)-G_strainEigen_Ave(ic,0),2e0)));
  normalize(G_strainEigenVector_Ave,ic);
  //cout <<"5"<<endl;
}

// #################################################################
/**
 * @brief calc stress tensor of SantVenant material
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] numOfNodeInElm   number of node in each element
 * @param [in] numOfGaussPoint  number of Gauss point set in each element
 * @param [in] option           true or faluse: calculate tangential stiffness matrix or not.
 */
double muscularHydrostat::Muscle::postProcess_inGaussIntegral(const int &ic,
const int &numOfNodeInElm,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,ARRAY2D<double> &dNdr,const double weight,double (&stress)[3][3],const bool option, double (&G_strainEigen)[3],double (&G_strainEigenVector)[3][3])
{
  //cout <<"start ingauss"<<endl;
  ARRAY2D<double> dNdx(numOfNodeInElm,3);
  double detJ,volume=0e0,J;
  double dXdr[3][3],dxdr[3][3],drdX[3][3],drdx[3][3];
  double sigma[3][3],S[3][3],C[3][3],invC[3][3],F[3][3],B[3][3];
  double elasticityTensor_ref[3][3][3][3],elasticityTensor_current[3][3][3][3],tangentCoefficient[3][3][3][3];

  //PDL
  double Ic4bar;
  int fiberNum=0;
  double lambda;
  double a0[3],a[3];

  //sigma an-isotropic term
  //for(int i=0;i<3;i++) a0[i]=fiberDirection_elm(ic,i);

  FEM_MathTool::calc_dxdr(dxdr,dNdr,x_current,numOfNodeInElm);
  FEM_MathTool::calc_dNdx(dNdx,dNdr,dxdr,numOfNodeInElm);
  detJ = mathTool::calcDeterminant_3x3(dxdr);
  volume += detJ * weight;

  FEM_MathTool::calc_dXdr(dXdr,dNdr,x_ref,numOfNodeInElm);
  mathTool::calcInverseMatrix_3x3(drdX,dXdr);

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      F[i][j]=0e0;
      for(int k=0;k<3;k++) F[i][j] += dxdr[i][k]*drdX[k][j];
    }
  }

  //todo initialStretchRatio：内部・外部入力の違いにより結果が異なる？
  //-------------initial stretch-------------------
  for(int ik=0;ik<fibers[ic].fiber.size();ik++){
    double F_initial[3][3]={},Ftmp[3][3]={};

    int fiberNumber = static_cast<int>(fibers[ic].fiber[ik].group);
    double initialStretchRatio = Material[fiberNumber].initialStretch;

    calc_F_initial(F_initial,fibers[ic].fiber[ik].a0,initialStretchRatio);
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        Ftmp[i][j]=0e0;
        for(int k=0;k<3;k++) Ftmp[i][j]+=F[i][k]*F_initial[k][j];
      }
    }
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) F[i][j]=Ftmp[i][j];
    }
  }
  //--------------------------------------------

  J = mathTool::calcDeterminant_3x3(F);

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      C[i][j]=0e0;
      for(int k=0;k<3;k++){
        C[i][j]+=F[k][i]*F[k][j];
      }
    }
  }
  double G_strain[3][3];
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      G_strain[i][j] = 5e-1*(C[i][j]-I2[i][j]);

      // if(ic==0){
      // cout << "check_Greenstrain" << endl;
      // cout << G_strain[i][j] << endl;
      // }

    }
  }
  calcEigen(G_strain,G_strainEigen,G_strainEigenVector);
  //cout <<"end ingauss"<<endl;
 //for(int i=0;i<3;i++){
 //  for(int j=0;j<3;j++){  
 //    // if(ic==0){
 //    // cout << "check_Green_principal_strain" << endl;
 //    // cout << G_strainEigen[i][j] << endl;
 //    // }

 //  }
 //}

  if(option==false) return volume;
  return volume;
}

  // double sigma[3][3];

  // for(int i=0;i<3;i++){
  //   for(int j=0;j<3;j++){
  //      sigma[i][j] = nu/J*(B[i][j]-I2[i][j])+lambda/J*(log(J))*I2[i][j];
  //   }
  // }

  //       //calc_internal force vector
  // for(int p=0;p<numOfNodeInElm;p++){
  //   for(int i=0;i<3;i++){
  //     for(int j=0;j<3;j++){
  //       Qu[ic](p,i) += sigma[i][j] * dNdx(p,j) * detJ * weight;
  //     }
  //   }
  // }

  // if(option==false) return;

  // double c4[3][3][3][3];
  // for(int i=0;i<3;i++){
  //   for(int j=0;j<3;j++){
  //     for(int k=0;k<3;k++){
  //       for(int l=0;l<3;l++) c4[i][j][k][l] = lambda/J*I2[i][j]*I2[k][l]+2e0/J*(nu-lambda*log(J))*I4[i][j][k][l];
  //     }
  //   }
  // }

//   //calc_tangential_stiffness_matrix
//   for(int p=0;p<numOfNodeInElm;p++){
//     for(int q=0;q<numOfNodeInElm;q++){
//       for(int i=0;i<3;i++){
//         for(int j=0;j<3;j++){
//           for(int k=0;k<3;k++){
//             for(int l=0;l<3;l++){
//               Ku[ic](p,q,i,j) += dNdx(p,k)*(c4[i][k][j][l]+sigma[k][l]*I2[i][j])*dNdx(q,l) * detJ * weight;
//             }
//           }
//         }
//       }
//     }
//   }
// }
