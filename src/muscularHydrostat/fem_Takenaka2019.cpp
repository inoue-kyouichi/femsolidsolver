/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "muscularHydrostat.h"

using namespace std;

// #################################################################
/**
 * @brief calc stress tensor of SantVenant material
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] option           true or false: calculate tangential stiffness matrix or not.
 */
void muscularHydrostat::Muscle::calcStressTensor_Takenaka2019_element_spatialForm(const int ic,ARRAY2D<double> &U_tmp,const bool option)
{
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

  switch(element[ic].meshType){
    case VTK_TETRA:
      ShapeFunction3D::C3D4_dNdr(dNdr,gTet.point[0][0],gTet.point[0][1],gTet.point[0][2],gTet.point[0][3]);
      Takenaka2019_inGaussIntegral(ic,numOfNodeInElm,x_current,x_ref,dNdr,gTet.weight[0]*1e0/6e0,stress,option);
      break;
    case VTK_HEXAHEDRON:
      for(int i1=0;i1<2;i1++){
        for(int i2=0;i2<2;i2++){
          for(int i3=0;i3<2;i3++){
            ShapeFunction3D::C3D8_dNdr(dNdr,g.point[i1],g.point[i2],g.point[i3]);
            Takenaka2019_inGaussIntegral(ic,numOfNodeInElm,x_current,x_ref,dNdr,g.weight[i1]*g.weight[i2]*g.weight[i3],stress,option);
          }
        }
      }
      break;
    case VTK_QUADRATIC_TETRA:
      for(int i1=0;i1<4;i1++){
        ShapeFunction3D::C3D10_dNdr(dNdr,gTet2.point[i1][0],gTet2.point[i1][1],gTet2.point[i1][2],gTet2.point[i1][3]);
        Takenaka2019_inGaussIntegral(ic,numOfNodeInElm,x_current,x_ref,dNdr,gTet2.weight[i1]*1e0/6e0,stress,option);
      }
      break;
    case VTK_QUADRATIC_HEXAHEDRON:
      break;
    case VTK_WEDGE:
      for(int i1=0;i1<2;i1++){
        ShapeFunction3D::C3D6_dNdr(dNdr,gTri.point[0][0],gTri.point[0][1],gTri.point[0][2],g.point[i1]);
        Takenaka2019_inGaussIntegral(ic,numOfNodeInElm,x_current,x_ref,dNdr,5e-1*gTri.weight[0]*g.weight[i1],stress,option);
      }
      break;
    default:
      cout << "undefined mesh type" << endl;
      exit(1);
  }
}

// #################################################################
/**
 * @brief calc stress tensor of PDL with using selective reduced integration
 * @param [in] ic               element number
 * @param [in] U_tmp            displacement vector
 * @param [in] numOfNodeInElm   number of node in each element
 * @param [in] numOfGaussPoint  number of Gauss point set in each element
 * @param [in] option           true or faluse: calculate tangential stiffness matrix or not.
 */
void muscularHydrostat::Muscle::Takenaka2019_inGaussIntegral(const int &ic,
const int &numOfNodeInElm,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,ARRAY2D<double> &dNdr,const double weight,double (&stress)[3][3],const bool mainLoop)
{
  ARRAY2D<double> dNdx(numOfNodeInElm,3);
  double detJ,volume=0e0,J;
  double dXdr[3][3],dxdr[3][3],drdX[3][3],drdx[3][3];
  double sigma[3][3],S[3][3],C[3][3],invC[3][3],F[3][3],B[3][3];
  double elasticityTensor_ref[3][3][3][3],elasticityTensor_current[3][3][3][3],tangentCoefficient[3][3][3][3];

  //nearly-incompressible material
  double pressure,dpressure,Siso[3][3],S_aniso[3][3],Sbar[3][3],P4[3][3][3][3],P4bar[3][3][3][3];
  double invC_odot[3][3][3][3],C4iso[3][3][3][3],C4vol[3][3][3][3],C4bar[3][3][3][3],term2,term4,term4_2;

  //PDL
  double Ic4bar;
  int fiberNum=0;
  double lambda;

  const double c10 = 1.037e3*1e-6; //[MPa]
  const double c20 = 4.86e2*1e-6; //[MPa]
  const double af = 1.0e3*1e-6;//1.0e6*1e-3; //[MPa]
  const double K = 1e5*1e-6; //[MPa]
  //const double c4 = 4.86e2*1e-3;
  //const double c3 = 1e5*1e-3;
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

  mathTool::calcInverseMatrix_3x3(invC,C);

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) P4[i][j][k][l]=I4[i][j][k][l]-1e0/3e0*invC[i][j]*C[k][l];
      }
    }
  }

  //------------------------------specific routine-----------------------------
  //S_bar
  double I_C1 = C[0][0]+C[1][1]+C[2][2];
  double I_C1_mod = I_C1 * pow(J,-2e0/3e0);
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) Sbar[i][j]=2e0*c10*I2[i][j] + 4e0*c20*(I_C1_mod-3e0)*I2[i][j];
    //for(int j=0;j<3;j++) Sbar[i][j]=2e0*c10*I2[i][j]; //an-iso
  }

  //sigma an-isotropic term
  /*Ic4bar=0e0;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) Ic4bar += a0[i]*C[i][j]*a0[j]*pow(J,-2e0/3e0);
    // for(int j=0;j<3;j++) Ic4bar += a0[i]*C[i][j]*a0[j];
  }

  lambda=sqrt(Ic4bar);
  for(int i=0;i<3;i++){
    a[i]=0e0;
    for(int j=0;j<3;j++) a[i] += F[i][j] * a0[j];
  }
  for(int i=0;i<3;i++) a[i] = a[i]/lambda;

   if(Ic4bar<1e0){
     term4=0e0;
     term4_2=0e0;
   }else{
     term4=c4*(exp(Ic4bar-1e0)-1e0);
     term4_2=c4*exp(Ic4bar-1e0);
   }
   
   for(int i=0;i<3;i++){
     for(int j=0;j<3;j++) Sbar[i][j]+=2e0*term4*a0[i]*a0[j];
   }*/

  //S_iso
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      Siso[i][j]=0e0;
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) Siso[i][j] += pow(J,-2e0/3e0)*P4[i][j][k][l]*Sbar[k][l];
      }
    }
  }

  //S_vol or hydrostatic pressure
  // pressure=c3*(J-1e0/J);            //
  // dpressure=c3*(1e0+1e0/(J*J));     //an-iso
  pressure=K*(J-1e0);
  dpressure=K;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      sigma[i][j] = pressure*I2[i][j];
    }
  }

  //transformation
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) sigma[i][j] += F[i][k]*Siso[k][l]*F[j][l]/J;
      }
    }
  }

  if(mainLoop==false){
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) stress[i][j]=sigma[i][j];
    }
    return;
  }

  //--------muscle contraction----------------------
  double contraction[3][3]={};
  for(int ik=0;ik<fibers[ic].fiber.size();ik++){

    for(int i=0;i<3;i++){
      a[i]=0e0;
      for(int j=0;j<3;j++) a[i] += F[i][j] * fibers[ic].fiber[ik].a0[j];
    }
    double lambda = sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]);
    for(int i=0;i<3;i++) a[i] = a[i]/lambda;
    int fiberNumber = static_cast<int>(fibers[ic].fiber[ik].group);
    double contractionCoefficient = Material[fiberNumber].contractionCoefficient;
    //contractionCoefficient *= 5e-1*(c10+c20)/(fibers[ic].fiber.size());
    contractionCoefficient *= af/(fibers[ic].fiber.size());

    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++){
        contraction[i][j] += contractionCoefficient*a[i]*a[j];///J;
      }
    }
  }
  //--------end muscle contraction----------------------

  //calc_internal force vector
  for(int p=0;p<numOfNodeInElm;p++){
    for(int i=0;i<3;i++){
      for(int j=0;j<3;j++) Qu[ic](p,i) += (sigma[i][j]+contraction[i][j]) * dNdx(p,j) * detJ * weight;
    }
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) invC_odot[i][j][k][l]=5e-1*(invC[i][k]*invC[j][l]+invC[i][l]*invC[j][k]);
      }
    }
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) P4bar[i][j][k][l]=invC_odot[i][j][k][l]-1e0/3e0*(invC[i][j]*invC[k][l]);
      }
    }
  }

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++) C4vol[i][j][k][l]=J*(pressure+J*dpressure)*invC[i][j]*invC[k][l]-2e0*J*pressure*invC_odot[i][j][k][l];
      }
    }
  }

  //------------------------------specific routine-----------------------------
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++){
        // C4bar[i][j][k][l]=term4_2*a0[i]*a0[j]*a0[k]*a0[l]; //an-iso
          C4bar[i][j][k][l]=2e0*c20*I2[i][j]*I2[k][l];
        }
      }
    }
  }

  term2=0e0;
  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++) term2+=Sbar[i][j]*C[i][j];
  }
  term2=term2*pow(J,-2e0/3e0);

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++){
          C4iso[i][j][k][l]=0e0;
          for(int m=0;m<3;m++){
            for(int n=0;n<3;n++){
              for(int p=0;p<3;p++){
                for(int q=0;q<3;q++) C4iso[i][j][k][l]+=P4[i][j][m][n]*C4bar[m][n][p][q]*P4[k][l][p][q];
              }
            }
          }
          C4iso[i][j][k][l] += 2e0/3e0*term2*P4bar[i][j][k][l]-2e0/3e0*(invC[i][j]*Siso[k][l]+Siso[i][j]*invC[k][l]);
          C4iso[i][j][k][l] += C4vol[i][j][k][l];
                // C4iso[i][j][k][l]+=term4_2*a0[i]*a0[j]*a0[k]*a0[l]; //an-iso
        }
      }
    }
  }

        //------------------------------end specific routine-----------------------------

  FEM_MathTool::tensorPushForward_4order(elasticityTensor_current,C4iso,F,J);

  for(int i=0;i<3;i++){
    for(int j=0;j<3;j++){
      for(int k=0;k<3;k++){
        for(int l=0;l<3;l++){
          tangentCoefficient[i][j][k][l] = elasticityTensor_current[i][k][j][l]+sigma[k][l]*I2[i][j];
        }
      }
    }
  }

        //-------------------------------------------------------

        //calc_tangential_stiffness_matrix
  for(int p=0;p<numOfNodeInElm;p++){
    for(int q=0;q<numOfNodeInElm;q++){
      for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
          for(int k=0;k<3;k++){
            for(int l=0;l<3;l++){
              Ku[ic](p,q,i,j) += dNdx(p,k)*(tangentCoefficient[i][j][k][l]*dNdx(q,l)) * detJ * weight;
            }
          }
        }
      }
    }
  }
}
