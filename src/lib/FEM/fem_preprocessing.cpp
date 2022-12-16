/**
 * @file fem_preprocessing.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "fem.h"

using namespace std;

// #################################################################
/**
 * @brief initialize FEM class
 */
void Fem::initialize(TextParser &tp)
{
  inputDomainInfo(tp);
  // string restartDirName="Restart_"+to_string(dataNumber);
  // mkdir(restartDirName.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);
  allocate();
  // restart_setting(dataNumber,Restart);
  inputDirichletBoundaryInfo(tp);
}

// #################################################################
/**
 * @brief output information from TP file
 */
void Fem::restart_setting(const int dataNumber,const bool Restart,TextParser &tp)
{
  FILE *fp;
  string output;

  string str,base_label,label;

  base_label = "/Domain";
  if(Restart==false){
    output = "Restart_"+to_string(dataNumber)+"/x0.dat";
    if ((fp = fopen(output.c_str(), "w")) == NULL) {
      cout << output << " open error" << endl;
      exit(1);
    }
    for(int i=0;i<numOfNode;i++){
      fprintf(fp,"%e %e %e\n",x0(i,0),x0(i,1),x0(i,2));
    }
    fclose(fp);
  }else{
    label = base_label + "/referenceNodeFile";
    if ( !tp.getInspectedValue(label, output)){
      cout << output <<" is not set" << endl;
      exit(0);
    }
    if ((fp = fopen(output.c_str(), "r")) == NULL) {
      cout << output << " open error" << endl;
      exit(1);
    }
    for(int i=0;i<numOfNode;i++){
      fscanf(fp,"%lf %lf %lf\n",&x0(i,0),&x0(i,1),&x0(i,2));
      for(int j=0;j<3;j++) x(i,j)=x0(i,j);
    }
    fclose(fp);

    label = base_label + "/restartNodeFile";
    if ( !tp.getInspectedValue(label, output)){
      cout << "outputDir is not set" << endl;
      exit(0);
    }
    if ((fp = fopen(output.c_str(), "r")) == NULL) {
      cout << output << " open error" << endl;
      exit(1);
    }
    for(int i=0;i<numOfNode;i++){
      fscanf(fp,"%lf %lf %lf\n",&U(i,0),&U(i,1),&U(i,2));
    }
    fclose(fp);
  }
}

// #################################################################
/**
 * @brief domain information from tp file
 */
void Fem::inputDomainInfo(TextParser &tp)
{
  string str,base_label,label,inputDir;

  base_label = "/Domain";

  label = base_label + "/inputDir";
  if ( !tp.getInspectedValue(label,inputDir)){
    cout << "data format is not set" << endl;
    exit(0);
  }
  string file1,file2,file3,file4,file5;
  label = base_label + "/nodeFile";
  if ( !tp.getInspectedValue(label, file1)){
    cout << label << " is not found" << endl;
    exit(0);
  }
  label = base_label + "/elementFile";
  if ( !tp.getInspectedValue(label,file2)){
    cout << label << " is not found" << endl;
    exit(0);
  }
  label = base_label + "/meshTypeFile";
  if ( !tp.getInspectedValue(label,file3)){
    cout << label << " is not found" << endl;
    exit(0);
  }

  file1=inputDir+"/"+file1;
  file2=inputDir+"/"+file2;
  file3=inputDir+"/"+file3;
  set_geometry(file1,file2,file3);

}

// #################################################################
/**
 * @brief Dirichlet information from TP file
 */
void Fem::inputDirichletBoundaryInfo(TextParser &tp)
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
  label = base_label + "/dirichletValueFile";
  if ( !tp.getInspectedValue(label, Dvalue_file)){
    cout << label << " is not found" << endl;
  }
  D_file=inputDir+"/"+D_file;
  Dvalue_file=inputDir+"/"+Dvalue_file;
  set_dirichlet(D_file,Dvalue_file);
}

// #################################################################
/**
 * @brief Neumann information from TP file
 */
void Fem::inputNeumannBoundaryInfo(TextParser &tp)
{
  string str,base_label,label,inputDir;
  string N_file,Nvalue_file,NmeshType_file;

  base_label = "/Domain";
  label = base_label + "/inputDir";
  if ( !tp.getInspectedValue(label,inputDir)){
    cout << "data format is not set" << endl;
    exit(0);
  }

  label = base_label + "/neumannFile";
  if ( !tp.getInspectedValue(label, N_file)){
    cout << label << " is not found" << endl;
  }
  label = base_label + "/neumannValueFile";
  if ( !tp.getInspectedValue(label, Nvalue_file)){
    cout << label << " is not found" << endl;
  }
  label = base_label + "/neumannMeshTypeFile";
  if ( !tp.getInspectedValue(label, NmeshType_file)){
    cout << label << " is not found" << endl;
  }
  set_neumann(N_file,NmeshType_file,Nvalue_file);
}

// #################################################################
/**
 * @brief allocation.
 */
void Fem::allocate()
{
  Mass.resize(numOfElm);
  Qu.resize(numOfElm);
  Ku.resize(numOfElm);

  for(int ic=0;ic<numOfElm;ic++){
    int numOfNodeInElm = element[ic].node.size();
    Mass[ic].allocate(numOfNodeInElm,numOfNodeInElm);
    Qu[ic].allocate(numOfNodeInElm,3);
    Ku[ic].allocate(numOfNodeInElm,numOfNodeInElm,3,3);
  }

  volume.allocate(numOfElm);
  volume0.allocate(numOfElm);
  volumeChangeRatio.allocate(numOfElm);

  Mises.allocate(numOfElm);
  Mises_strain.allocate(numOfElm);
  AEigen_Ave.allocate(numOfElm,3);
  sigmaEigen_Ave.allocate(numOfElm,3);
  G_strainEigen_Ave.allocate(numOfElm,3);
  AEigenVector_Ave.allocate(numOfElm,3,3);
  sigmaEigenVector_Ave.allocate(numOfElm,3,3);
  G_strainEigenVector_Ave.allocate(numOfElm,3,3);

  U.allocate(numOfNode,3);
  RHS.allocate(numOfNode,3);
  innerForce.allocate(numOfNode,3);
  externalForce.allocate(numOfNode,3);

  for(int ic=0;ic<numOfElm;ic++) volumeChangeRatio(ic)=1e0;

  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++){
      RHS(i,j) = 0e0;
      externalForce(i,j) = 0e0;
      U(i,j) = 0e0;
    }
  }

  //dirichlet boundary conditions
  ibd.allocate(numOfNode,3);
  bd.allocate(numOfNode,3);
  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++){
      ibd(i,j)=1;
      bd(i,j)=0e0;
    }
  }

}
