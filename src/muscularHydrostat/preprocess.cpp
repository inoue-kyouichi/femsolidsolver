/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "muscularHydrostat.h"
#include <ostream>
#include <fstream>

using namespace std;

void muscularHydrostat::Muscle::preprocess()
{
  inputDomainInfo(tp);
  allocate();

  fibers.resize(numOfElm);
  //fiberDirection_elm.allocate(numOfElm,3);
  inputFiberInfo_cal(tp);
  //inputFiberInfo(tp);
  //inputMaterialInfo(tp);

  inputMaterialParameters(tp);

  inputDirichletInfo(tp);

  inputSolverInfo(tp);
  inputOutputInfo(tp);

  omp_set_num_threads(OMPnumThreads);

  mkdir(outputDir.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);

  string file = outputDir + "/test_fiber.vtu";
  export_vtu_fiber(file);

  //CSR setting
  PARDISO.initialize(3*numOfNode);
  PARDISO.CSR_initialize(inb,numOfNode,3);
}

// #################################################################
/**
 * @brief solver information from TP file
 */
void muscularHydrostat::Muscle::inputMaterialParameters(TextParser &tp)
{
  string str,base_label,label;
  int tmp;

  base_label = "/Fibers";

  int numOfMaterials=3; //input yourself
  // label = base_label + "/numOfMaterials";
  // if ( !tp.getInspectedValue(label, numOfMaterials)){
  //   cout << label << " is not set" << endl;
  //   exit(0);
  // }

  Material.resize(numOfMaterials+1);

  for(int ic=1;ic<=numOfMaterials;ic++){

    string base_label2 = base_label + "/F"+to_string(ic);
    label = base_label2 + "/contractionCoefficient";
    if ( !tp.getInspectedValue(label, Material[ic].contractionCoefficient)){
      cout << label << " is not set" << endl;
      exit(0);
    }

    label = base_label2 + "/initialStretch";
    if ( !tp.getInspectedValue(label, Material[ic].initialStretch)){
      cout << label << " is not set" << endl;
      exit(0);
    }
  }

}

// #################################################################
/**
 * @brief solver information from TP file
 */
void muscularHydrostat::Muscle::inputSolverInfo(TextParser &tp)
{
  string str,base_label,label;
  int tmp;

  base_label = "/Solver";

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
void muscularHydrostat::Muscle::inputFiberInfo_cal(TextParser &tp)
{
  string str,base_label,label,inputDir;
  base_label = "/Domain";

  label = base_label + "/inputDir";
  if ( !tp.getInspectedValue(label,inputDir)){
    cout << "data format is not set" << endl;
    exit(0);
  }
  string D_file,Dvalue_file;
  label = base_label + "/fiberFile";
  if ( !tp.getInspectedValue(label, D_file)){
    cout << label << " is not found" << endl;
    exit(1);
  }
  D_file=inputDir+"/"+D_file;

  FILE *fp;
  char buf[50];
  if ((fp = fopen(D_file.c_str(), "r")) == NULL) {
    cout << "file open error" << endl;
    exit(1); 
  }
  fgets(buf,30,fp);
  for(int ic=0;ic<numOfElm;ic++){
    int number;
    fscanf(fp,"%d",&number);
    fibers[ic].fiber.resize(number);
  }

  for(int ic=0;ic<numOfElm;ic++){
    for(int i=0;i<fibers[ic].fiber.size();i++){
      int number;
      fscanf(fp,"%d",&number);
      fibers[ic].fiber[i].group = static_cast<FiberGroup>(number);
    }
  }

  for(int ic=0;ic<numOfElm;ic++){
    for(int i=0;i<fibers[ic].fiber.size();i++){
      double tmp;
      fscanf(fp,"%lf",&tmp);
      fibers[ic].fiber[i].a0[0] = tmp;
    }
  }

  for(int ic=0;ic<numOfElm;ic++){
    for(int i=0;i<fibers[ic].fiber.size();i++){
      double tmp;
      fscanf(fp,"%lf",&tmp);
      fibers[ic].fiber[i].a0[1] = tmp;
    }
  }

  for(int ic=0;ic<numOfElm;ic++){
    for(int i=0;i<fibers[ic].fiber.size();i++){
      double tmp;
      fscanf(fp,"%lf",&tmp);
      fibers[ic].fiber[i].a0[2] = tmp;
    }
  }

  fclose(fp);

  // ifstream file(D_file);
  // if(!file){
  //   cout << "Error:Input "<< D_file << " not found" << endl;
  //   exit(1);
  // }

  // string tmp,dtmp;
  // for(int i=0;i<numOfElm;i++){
  //   getline(file,str);
  //   istringstream stream(str);
  //   for(int j=0;j<3;j++){
  //     getline(stream,tmp,' ');
  //     fiberDirection_elm(i,j) = stod(tmp);
  //   }
  // }
}

// #################################################################
/**
 * @brief tentative
 */
void muscularHydrostat::Muscle::inputFiberInfo(TextParser &tp)
{
  string str,base_label,label,inputDir;
  base_label = "/Domain";

  label = base_label + "/inputDir";
  if ( !tp.getInspectedValue(label,inputDir)){
    cout << "data format is not set" << endl;
    exit(0);
  }
  string D_file,Dvalue_file;
  label = base_label + "/fiberFile";
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
    istringstream stream(str);
    for(int j=0;j<3;j++){
      getline(stream,tmp,' ');
      fiberDirection_elm(i,j) = stod(tmp);
    }
  }
}

// #################################################################
/**
 * @brief tentative
 */
void muscularHydrostat::Muscle::inputMaterialInfo(TextParser &tp)
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

  //string tmp,dtmp;
  //for(int i=0;i<numOfElm;i++){
  //  getline(file,str);
  //  element[i].materialType = stoi(str);
  //}
  string tmp,dtmp;
  for(int i=0;i<numOfElm;i++){
    getline(file,str);
    int mat = stoi(str);
    switch(mat){
      case static_cast<int>(MaterialType::M0):
        element[i].materialType = static_cast<int>(MaterialType::M0);
        break;
      case static_cast<int>(MaterialType::M1):
        element[i].materialType = static_cast<int>(MaterialType::M1);
        break;
      case static_cast<int>(MaterialType::M2):
        element[i].materialType = static_cast<int>(MaterialType::M2);
        break;
      default:
        cout << "undefined material. Exit..." << endl;
        exit(1);
    }
  }

}

// #################################################################
/**
 * @brief Dirichlet information from TP file
 */
void muscularHydrostat::Muscle::inputDirichletInfo(TextParser &tp)
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
void muscularHydrostat::Muscle::inputOutputInfo(TextParser &tp)
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

// #################################################################
/**
 * @brief calc boundary conditions
 * @param [in] stress
 */
void muscularHydrostat::Muscle::export_vtu_fiber(const string &file)
{
  FILE *fp;
  if ((fp = fopen(file.c_str(), "w")) == NULL) {
    cout << file << " open error" << endl;
    exit(1); 
  }

  fprintf(fp,"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp,"<UnstructuredGrid>\n");
  fprintf(fp,"<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n",numOfNode,numOfElm);
  fprintf(fp,"<Points>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",x(i,0),x(i,1),x(i,2));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Points>\n");
  fprintf(fp,"<Cells>\n");
  fprintf(fp,"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    for(int j=0;j<element[i].node.size();j++) fprintf(fp,"%d ",element[i].node[j]);
    fprintf(fp,"\n");
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  int num=0;
  for(int i=0;i<numOfElm;i++){
    num += element[i].node.size();
    fprintf(fp,"%d\n",num);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++) fprintf(fp,"%d\n",element[i].meshType);
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Cells>\n");

  fprintf(fp,"<PointData Vectors=\"displacement[m]\">\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"displacement[m]\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",U(i,0),U(i,1),U(i,2));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");

  fprintf(fp,"<CellData>");
  for(int ic=1;ic<=3;ic++){
    fprintf(fp,"<DataArray type=\"Float64\" Name=\"fiber_%d\" NumberOfComponents=\"3\" format=\"ascii\">\n",ic);
    for(int i=0;i<numOfElm;i++){
      bool test = false;
      for(int j=0;j<fibers[i].fiber.size();j++){
        if(static_cast<int>(fibers[i].fiber[j].group)==ic){
          test = true;
          fprintf(fp,"%e %e %e\n",fibers[i].fiber[j].a0[0],fibers[i].fiber[j].a0[1],fibers[i].fiber[j].a0[2]);
        }
      }
      if(test==false){
        fprintf(fp,"0e0 0e0 0e0\n");
      }
    }
    fprintf(fp,"</DataArray>\n");
  }
  fprintf(fp,"</CellData>\n");
  fprintf(fp,"</Piece>");
  fprintf(fp,"</UnstructuredGrid>");
  fprintf(fp,"</VTKFile>");
  fclose(fp);
}

