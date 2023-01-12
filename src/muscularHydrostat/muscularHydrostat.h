#ifndef _MUSCULAR_HYDROSTAT_H_
#define _MUSCULAR_HYDROSTAT_H_

//##################################################################################
//
// Eye mechanics simulation
//
// Copyright (c) 2020 Biomechanics Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   eye_Mechanics.h
 * @brief  PDL Header
 * @author T. Otani
 */

#include "fem.h"
#include "pardiso_solver.h"

namespace muscularHydrostat {

class MaterialModel{
  public:
  double contractionCoefficient;
  double initialStretch;
};

enum class  MaterialType{
  M0 = 0,
  M1 = 1,
  M2 = 2
};

//enum class FiberGroup{
//  GGa = 1,
//  GGm = 2,
//  GGp = 3,
//  HG = 4,
//  SG = 5,
//  SLa = 6,
//  SLp = 7,
//  IL = 8,
//  T = 9,
//  V = 10
//};
enum class FiberGroup{
  BD = 1,
  GGa = 2,
  GGm = 3,
};

class Fiber{
  public:
    double a0[3];
    FiberGroup group;
};

class FibersInElement{
  public:
  std::vector<Fiber> fiber;
};

class Muscle : public Fem {
 public:
  Muscle(){
    //Material.resize(10);
    Material.resize(2);
  };
  ~Muscle(){};
  TextParser tp;
  std::string outputDir,fileName;

  // int numOfBoundaryElm;
  // VECTOR1D<ElementType> boundaryElement;
  std::vector<MaterialModel> Material;
  std::vector<FibersInElement> fibers;

 private:
  int dataNumber;
  int Restart;
  int OMPnumThreads;
  int maxIteration,NRiteration;
  double NRtolerance;
  double relaxation;
  PARDISO_solver PARDISO;
  ARRAY2D<double> fiberDirection_elm;

 public: 
  void preprocess();
  void femSolidAnalysis();
  void calcStressTensor();
  void postProcess_element_spatialForm(const int ic,ARRAY2D<double> &U_tmp,const bool option);
  double postProcess_inGaussIntegral(const int &ic,
const int &numOfNodeInElm,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,ARRAY2D<double> &dNdr,const double weight,double (&stress)[3][3],const bool option, double (&G_strainEigen)[3],double (&G_strainEigenVector)[3][3]);

 private:
  bool NRscheme();
  void set_rhs_statics();
  void set_forceddisplacement(int ic);
  void calcStressTensor_SantVenant_element_spatialForm(const int ic,ARRAY2D<double> &U_tmp,const bool option);
  double SantVenant_inGaussIntegral(ARRAY2D<double> &dNdr,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,
        ARRAY2D<double> &dNdx,const int numOfNodeInElm,const double weight,const int ic,const bool option);

  void calcStressTensor_muscle_element_spatialForm(const int ic,ARRAY2D<double> &U_tmp,const bool option);
  void muscle_inGaussIntegral(const int &ic,const int &numOfNodeInElm,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,ARRAY2D<double> &dNdr,const double weight,double (&stress)[3][3],const bool mainLoop);
  void calc_F_initial(double (&F_initial)[3][3],const double (&a0)[3],const double lambda);

  void calcStressTensor_Takenaka2019_element_spatialForm(const int ic,ARRAY2D<double> &U_tmp,const bool option);
  void Takenaka2019_inGaussIntegral(const int &ic,const int &numOfNodeInElm,ARRAY2D<double> &x_current,ARRAY2D<double> &x_ref,ARRAY2D<double> &dNdr,const double weight,double (&stress)[3][3],const bool mainLoop);
  void move_hyoid_bone();

  void inputMaterialParameters(TextParser &tp);
  void inputSolverInfo(TextParser &tp);
  void inputOutputInfo(TextParser &tp);
  void inputFiberInfo(TextParser &tp);
  void inputFiberInfo_cal(TextParser &tp);
  void inputMaterialInfo(TextParser &tp);
  void inputDirichletInfo(TextParser &tp);
  void forceddisplacement(TextParser &tp);
  void inputhyoidboneinfo(TextParser &tp);

  void export_vtu(const std::string &file);
  void export_vtu_strain(const std::string &file,ARRAY2D<double> &strainEigen_ave,ARRAY3D<double> &strainEigenVector_ave,ARRAY1D<double> &Mises_strain);
  void export_vtu_fiber(const std::string &file);
  void export_vtu_hyoidbone(ARRAY2D<double> &x,const std::vector<ElementType> &element,const int &numOfNode,const int &numOfElm,const std::string &file);
  void export_vtu_extention_rate(const std::string &file,ARRAY1D<double> &extention);
};

}

#endif