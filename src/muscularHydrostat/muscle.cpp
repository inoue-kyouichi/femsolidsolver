/**
 * @file fem.cpp
 * @brief Fem class
 * @author T. Otani
 */

#include "muscularHydrostat.h"
#include <ostream>
#include <fstream>

using namespace std;

// #################################################################
/**
 * @brief fem solid analysis routine
 * @param [inout] PARDISO   PARDISO class
 */
void muscularHydrostat::Muscle::femSolidAnalysis()
{
  int output_iter=1;
  string output;

  if(Restart==0){
    output = outputDir + "/test_" + to_string(0) + ".vtu";
    fileIO::export_vtu(x,element,numOfNode,numOfElm,U,volumeChangeRatio,output);
  }

  // for(int ic=0;ic<numOfElm;ic++) calcVolume_hexa(ic,volume0,8,2,0);

  for(int loop=Restart+1;loop<=maxIteration;loop++){

    if(NRscheme()==true){
      printf("NaN detected in loop=%d. Exit...\n",loop);
      exit(1);
    }

    #pragma omp parallel for
    for(int i=0;i<numOfNode;i++){
      for(int j=0;j<3;j++) x(i,j) = x0(i,j) + U(i,j);
    }

    // FILE *fp;
    // output = "Restart/U_" + to_string(loop) + ".dat";
    // if ((fp = fopen(output.c_str(), "w")) == NULL) {
    //   cout << "file open error" << endl;
    //   exit(1);
    // }
    // for(int i=0;i<numOfNode;i++){
    //   fprintf(fp,"%e %e %e\n",U(i,0),U(i,1),U(i,2));
    // }
    // fclose(fp);

    // if(loop%output_iter==0){
    //   output = outputDir + "/test_" + to_string(loop/output_iter) + ".vtu";
    //   fileIO::export_vtu(x,element,numOfNode,numOfElm,U,volumeChangeRatio,lambda_ave,output);
    // }
  }

  output = outputDir + "/test.vtu";
  export_vtu(output);

}

// #################################################################
/**
 * @brief fem solid analysis routine
 * @param [inout] PARDISO   PARDISO class
 */
bool muscularHydrostat::Muscle::NRscheme()
{
    double residual,residual0,norm,norm0,min1=0e0,min2=0e0;
    string output;
    int num1,num2;

    for(int ic=1;ic<=NRiteration;ic++){

      calcStressTensor();  //calc K and Q
      //calcBoundaryForce();
      set_rhs_statics();
      set_forceddisplacement(ic);

      PARDISO.set_CSR_value3D(Ku,element,numOfNode,numOfElm,inb);
      PARDISO.set_CSR_dirichlet_boundary_condition3D(numOfNode,ibd);

      for(int i=0;i<numOfNode;i++){
        for(int j=0;j<3;j++) PARDISO.b[i+j*numOfNode]=RHS(i,j);
      }

      if(ic==1) residual0 = PARDISO.vector_norm(numOfNode*3,PARDISO.b);

      PARDISO.main(3*numOfNode,OMPnumThreads);
      norm = PARDISO.vector_norm(numOfNode*3,PARDISO.x);
      if(ic==1) norm0 = norm;
      corrector_statics(PARDISO.x,relaxation);
      // residual=line_search(PARDISO.x);

      //find minimum
      if(ic == 1) min1 = norm/norm0;
      if(norm/norm0 <= min1){
         min1 = norm/norm0;
         num1 = ic;
      }
      
      
      if(isnan(norm)){
        cout << "residual is nan. Exit..." << endl;
        cout <<  "min NR iter = " << num1 << endl;
        cout << "min norm/norm0 = " << min1 << endl;
        return true;
      }

      // for(int i=0;i<numOfElm;i++){
      //   calcVolume_hexa(i,volume,8,2,1);
      //   volumeChangeRatio[i]=volume[i]/volume0[i];
      // }
      for(int i=0;i<numOfNode;i++){
        for(int j=0;j<3;j++) x(i,j) = x0(i,j) + U(i,j);
      }
      // output = outputDir + "/test_NR_" + to_string(ic) + ".vtu";
      // fileIO::export_vtu(x,element,numOfNode,numOfElm,U,volumeChangeRatio,lambda_ave,output);

      printf("NR iter.=%d norm/norm0=%e\n",ic,norm/norm0);

      //find minimum
      if(ic == 1) min2 = norm/norm0;
      if(norm/norm0 <= min2){
         min2 = norm/norm0;
         num2 = ic;
         if(min2 < 5e-3){
          output = outputDir + "/min_NR_"+to_string(num2)+".vtu";
          export_vtu(output);
         }
      }
      //

      if(ic%10==0){
        output = outputDir + "/NR_"+to_string(ic/10)+".vtu";
        export_vtu(output);
      }
      if(norm/norm0<NRtolerance)  break;
       
      // if(test!=1 && ic>50) break;
    }
    return false;
}

// #################################################################
/**
 * @brief calc stress tensor
 */
void muscularHydrostat::Muscle::calcStressTensor()
{
  double elementVolume,volume;

  stress_tensor_initialize();

  volume = 0e0;

  #pragma omp parallel for
  for(int ic=0;ic<numOfElm;ic++){
    // calcStressTensor_muscle_element_spatialForm(ic,U,true);
    calcStressTensor_Takenaka2019_element_spatialForm(ic,U,true);
  }

  for(int ic=0;ic<numOfElm;ic++){
    for(int p=0;p<element[ic].node.size();p++){
      for(int i=0;i<3;i++){
        innerForce(element[ic].node[p],i) += Qu[ic](p,i);
      }
    }
  }
}

// #################################################################
/**
 * @brief assembly right-hand side
 */
void muscularHydrostat::Muscle::set_rhs_statics()
{
 #pragma omp parallel for
  for(int ic=0;ic<numOfNode;ic++){
    for(int j=0;j<3;j++){
      // RHS(ic,j) = (double)ibd(ic,j)*( - innerForce(ic,j)+boundaryForce(ic,j));
      RHS(ic,j) = (double)ibd(ic,j)*( - innerForce(ic,j));
    }
  }
}

// #################################################################
/**
 * @brief set forced displacement
 */
void muscularHydrostat::Muscle::set_forceddisplacement(int ic)
{
  ARRAY2D<double> Forceddisplacement;
  Forceddisplacement.allocate(numOfNode,3);
  #pragma omp parallel for
  for(int i=0;i<numOfNode;i++){
    for(int j=0;j<3;j++) Forceddisplacement(i,j) = 0e0;
  }

  for(int i=0; i<numOfFD; i++){
    int o = FD(i,0);
    Forceddisplacement(o,0) = 0e0;  //forceddisplacement for x 0.1[mm]
    Forceddisplacement(o,1) = -1e2;  //forceddisplacement for y 0.1[mm]
    Forceddisplacement(o,2) = 0e0;  //forceddisplacement for z 0.1[mm]
  }

  if(ic == 1){
    for(int i=0;i<numOfFD;i++){
      int o = FD(i,0);
      for(int j=0;j<3;j++){
        RHS(o,j) = Forceddisplacement(o,j);
      }
    }
  }

}

// #################################################################
/**
 * @brief calc boundary conditions
 * @param [in] stress
 */
void muscularHydrostat::Muscle::export_vtu(const string &file)
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

  fprintf(fp,"<PointData Vectors=\"displacement[mm]\">\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"displacement[mm]\" NumberOfComponents=\"3\" format=\"ascii\">\n");
  for(int i=0;i<numOfNode;i++){
    fprintf(fp,"%e %e %e\n",U(i,0),U(i,1),U(i,2));
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</PointData>\n");

  fprintf(fp,"<CellData>");
  fprintf(fp,"<DataArray type=\"Int64\" Name=\"Material\" NumberOfComponents=\"1\" format=\"ascii\">\n");
  for(int i=0;i<numOfElm;i++){
    fprintf(fp,"%d\n",element[i].materialType);
  }
  fprintf(fp,"</DataArray>\n");
  fprintf(fp,"</CellData>\n");
  fprintf(fp,"</Piece>");
  fprintf(fp,"</UnstructuredGrid>");
  fprintf(fp,"</VTKFile>");
  fclose(fp);
}
