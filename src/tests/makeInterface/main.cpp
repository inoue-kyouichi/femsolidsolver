//##################################################################################
//
// test: calc surface curvature
//
// Copyright (c) 2020  Mechanical and Bioengineering Systems Lab.,
//                     Department of Mechanical Science and Bioengineering,
//                     Graduate School of Engineering Science,
//                     Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   main_curvature.cpp
 * @author T. Otani
 * @detail Fem solid analysis with displacement control.
 */

#include "fem.h"
#include "rigidBody.h"

#include "glog/logging.h"

using namespace std;

void calcInterface(ARRAY1D<int> &interface,Fem &FEM,TriangleSet &obj);
void exportVTU(const std::string &output_file,Fem &FEM,ARRAY1D<int> &interface);

int main(int argc,char *argv[])
{
  google::InitGoogleLogging(argv[0]);
  google::InstallFailureSignalHandler();

  int numOfNode,numOfElm;
  Fem FEM;
  RigidBody RBD;
  TextParser tp;

  //read triangle surface
  if(argc!=2){
    cout << "please set input file name. Exit..." << endl;
    exit(1);
  }
  std::string input_file = argv[1];
  int ierror;
  if ((ierror = tp.read(input_file)) != TP_NO_ERROR) {
    printf("\tError at reading '%s' file\n", input_file.c_str());
    return 1;
  }

  omp_set_num_threads(2);

  FEM.inputDomainInfo(tp);
  RBD.initialize(tp);

  ARRAY1D<int> interface(FEM.numOfNode);
  for(int ic=0;ic<FEM.numOfNode;ic++) interface(ic) = 1;

  calcInterface(interface,FEM,RBD.obj);

  std::string output_file = "test.vtu";
  exportVTU(output_file,FEM,interface);

  return 0;
}

void calcInterface(ARRAY1D<int> &interface,Fem &FEM,TriangleSet &obj)
{
  for(int ic=0;ic<obj.numOfNode;ic++){
    cout << ic << endl;
    double minLengh = 1e10;
    int number;
    for(int i=0;i<FEM.numOfNode;i++){
      double length=0e0;
      for(int j=0;j<3;j++) length += pow(FEM.x(i,j)-obj.x(ic,j),2e0);
      length = sqrt(length);
      if(length<minLengh){
        minLengh = length;
        number = i;
      }
    }

    if(number==156){
      cout << minLengh << endl;
      exit(1);
    }


    if(minLengh<1e-12){
      interface(number) = 0;
    }
  }

}

void exportVTU(const std::string &output_file,Fem &FEM,ARRAY1D<int> &interface)
{
  FILE *fp;
  if ((fp = fopen(output_file.c_str(), "w")) == NULL) {
    cout << output_file << " open error" << endl;
    exit(1); 
  }

  fprintf(fp,"<VTKFile type=\"UnstructuredGrid\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n");
  fprintf(fp,"<UnstructuredGrid>\n");
  fprintf(fp,"<Piece NumberOfPoints= \"%d\" NumberOfCells= \"%d\" >\n",FEM.numOfNode,FEM.numOfElm);
  fprintf(fp,"<Points>\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"Points\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<FEM.numOfNode;i++){
    fprintf(fp,"%e %e %e\n",FEM.x(i,0),FEM.x(i,1),FEM.x(i,2));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Points>\n");
  fprintf(fp,"<Cells>\n");
  fprintf(fp,"<DataArray type=\"Int64\" Name=\"connectivity\" format=\"ascii\">\n");
  for(int i=0;i<FEM.numOfElm;i++){
    for(int j=0;j<FEM.element[i].node.size();j++) fprintf(fp,"%d ",FEM.element[i].node[j]);
    fprintf(fp,"\n");
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"Int64\" Name=\"offsets\" format=\"ascii\">\n");
  int num=0;
  for(int i=0;i<FEM.numOfElm;i++){
    num += FEM.element[i].node.size();
    fprintf(fp,"%d\n",num);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"<DataArray type=\"UInt8\" Name=\"types\" format=\"ascii\">\n");
  for(int i=0;i<FEM.numOfElm;i++) fprintf(fp,"%d\n",FEM.element[i].meshType);
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</Cells>\n");

  fprintf(fp,"<PointData>\n");

  fprintf(fp,"<DataArray type=\"UInt32\" Name=\"interface\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<FEM.numOfNode;i++){
    fprintf(fp,"%d\n",interface(i));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");

  fprintf(fp,"<CellData>");
  fprintf(fp,"</CellData>\n");
  fprintf(fp,"</Piece>");
  fprintf(fp,"</UnstructuredGrid>");
  fprintf(fp,"</VTKFile>");
  fclose(fp);

}