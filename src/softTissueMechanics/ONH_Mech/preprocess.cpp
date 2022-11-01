/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "eye_mechanics.h"
#include <ostream>
#include <fstream>

using namespace std;

void EYE_Mechanics::EyeMech::preprocess()
{
  inputDomainInfo(tp);
  allocate();

  boundaryForce.allocate(numOfNode,3);

  inputMaterialInfo(tp);

  inputSurfaceInfo(tp);
  inputDirichletInfo(tp);

  // initialize(tp);
  inputSolverInfo(tp);
  inputOutputInfo(tp);

  omp_set_num_threads(OMPnumThreads);

  mkdir(outputDir.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);

  // pressureStiffnessMatrix.resize(numOfBoundaryElm);
  // for(int ic=0;ic<numOfBoundaryElm;ic++){
  //   pressureStiffness[ic].allocate(boundaryElement[ic].node.size(),boundaryElement[ic].node.size(),3,3);
  // }

  //CSR setting
  PARDISO.initialize(3*numOfNode);
  PARDISO.CSR_initialize(inb,numOfNode,3);

  string output = outputDir + "/" + fileName + "_boundary" + ".vtu";
  fileIO::export_vtu_boundary(x,element,numOfNode,numOfElm,ibd,bd,output);
}

// #################################################################
/**
 * @brief solver information from TP file
 */
void EYE_Mechanics::EyeMech::inputSolverInfo(TextParser &tp)
{
  string str,base_label,label;
  int tmp;

  base_label = "/Solver";

  label = base_label + "/boundaryPressure";
  if ( !tp.getInspectedValue(label, boundaryPressure)){
    cout << label << " is not set" << endl;
    exit(0);
  }

  label = base_label + "/dataNumber";
  if ( !tp.getInspectedValue(label, dataNumber)){
    cout << "maxiteration is not set" << endl;
    exit(0);
  }

  label = base_label + "/maxIteration";
  if ( !tp.getInspectedValue(label, maxIteration)){
    cout << "maxiteration is not set" << endl;
    exit(0);
  }

  label = base_label + "/NR_iteration";
  if ( !tp.getInspectedValue(label, NRiteration)){
    cout << "NLiteration is not set" << endl;
    exit(0);
  }

  label = base_label + "/NR_tolerance";
  if ( !tp.getInspectedValue(label, tmp)){
    cout << "NRtolerance is not set" << endl;
    exit(0);
  }
  NRtolerance = pow(1e-1,tmp);

  label = base_label + "/Restart";
  if ( !tp.getInspectedValue(label, Restart)){
    cout << label <<" is not set" << endl;
    exit(0);
  }

  label = base_label + "/OMPnumThreads";
  if ( !tp.getInspectedValue(label, OMPnumThreads)){
    cout << "OMPnumThreads is not set" << endl;
    exit(0);
  }

  label = base_label + "/relaxation";
  if ( !tp.getInspectedValue(label,relaxation)){
    cout << "NLiteration is not set" << endl;
    exit(0);
  }
}

// #################################################################
/**
 * @brief tentative
 */
void EYE_Mechanics::EyeMech::inputMaterialInfo(TextParser &tp)
{
  string str,base_label,label,inputDir;
  base_label = "/Domain";

  label = base_label + "/inputDir";
  if ( !tp.getInspectedValue(label,inputDir)){
    cout << "data format is not set" << endl;
    exit(0);
  }
  string D_file,Dvalue_file;
  label = base_label + "/materialTypeFile";
  if ( !tp.getInspectedValue(label, D_file)){
    cout << label << " is not found" << endl;
    exit(1);
  }
  D_file=inputDir+"/"+D_file;

  ifstream file(D_file);
  if(!file){
    cout << "Error:Input "<< D_file << " not found" << endl;
    exit(1);
  }

  string tmp,dtmp;
  for(int i=0;i<numOfElm;i++){
    getline(file,str);
    element[i].materialType = stoi(str);
  }

}

// #################################################################
/**
 * @brief tentative
 */
void EYE_Mechanics::EyeMech::inputSurfaceInfo(TextParser &tp)
{
  string str,base_label,label,inputDir;
  base_label = "/Domain";

  label = base_label + "/inputDir";
  if ( !tp.getInspectedValue(label,inputDir)){
    cout << "data format is not set" << endl;
    exit(0);
  }
  string D_file,Dvalue_file;
  label = base_label + "/innerSurfaceFile";
  if ( !tp.getInspectedValue(label, D_file)){
    cout << label << " is not found" << endl;
    exit(1);
  }
  D_file=inputDir+"/"+D_file;

  numOfBoundaryElm = fileIO::CountNumbersOfTextLines(D_file);
  boundaryElement.resize(numOfBoundaryElm);

  for(int ic=0;ic<numOfBoundaryElm;ic++){ 
    // meshtype_立方体ボリュームメッシュ、サーフェスが四角形
    // boundaryElement[ic].node.resize(4);
    // boundaryElement[ic].meshType = VTK_QUAD;
    // meshtype_四面体ボリュームメッシュ、サーフェスが三角形
    boundaryElement[ic].node.resize(3);
    boundaryElement[ic].meshType = VTK_TRIANGLE;
  }

  ifstream file(D_file);
  if(!file){
    cout << "Error:Input "<< D_file << " not found" << endl;
    exit(1);
  }

  string tmp,dtmp;
  for(int i=0;i<numOfBoundaryElm;i++){
    getline(file,str);
    istringstream stream(str);
    for(int j=0;j<boundaryElement[i].node.size();j++){
      getline(stream,tmp,' ');
      boundaryElement[i].node[j] = stoi(tmp);
    }
  }


}

// #################################################################
/**
 * @brief Dirichlet information from TP file
 */
void EYE_Mechanics::EyeMech::inputDirichletInfo(TextParser &tp)
{
  string str,base_label,label,inputDir;
  base_label = "/Domain";

  label = base_label + "/inputDir";
  if ( !tp.getInspectedValue(label,inputDir)){
    cout << "data format is not set" << endl;
    exit(0);
  }
  string D_file,Dvalue_file;
  label = base_label + "/dirichletFile";
  if ( !tp.getInspectedValue(label, D_file)){
    cout << label << " is not found" << endl;
  }
  D_file=inputDir+"/"+D_file;

  int numOfData = fileIO::CountNumbersOfTextLines(D_file);

  int node;
  string tmp,dtmp;

  ifstream file(D_file);
  if(!file){
    cout << "Error:Input "<< D_file << " not found" << endl;
    exit(1);
  }

  ibd.allocate(numOfNode,3);
  bd.allocate(numOfNode,3);

  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++){
      bd(i,j) = 0e0;
      ibd(i,j) = 1;
    }
  }

  for(int i=0;i<numOfData;i++){
    getline(file,str);
    istringstream stream(str);
    getline(stream,tmp,' ');
    int number = stoi(tmp);
    for(int j=0;j<3;j++){
      getline(stream,tmp,' ');
      ibd(number,j) = stoi(tmp);
    }
  }

}


// #################################################################
/**
 * @brief output information from TP file
 */
void EYE_Mechanics::EyeMech::inputOutputInfo(TextParser &tp)
{
  string str,base_label,label;

  base_label = "/Output";
  label = base_label + "/outputDir";
  if ( !tp.getInspectedValue(label, outputDir)){
    cout << "outputDir is not set" << endl;
    exit(0);
  }

  label = base_label + "/fileName";
  if ( !tp.getInspectedValue(label, fileName)){
    cout << "fileName is not set" << endl;
    exit(0);
  }
}