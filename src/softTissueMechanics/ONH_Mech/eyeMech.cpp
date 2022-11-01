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
 * @brief fem solid analysis routine
 * @param [inout] PARDISO   PARDISO class
 */
void EYE_Mechanics::EyeMech::femSolidAnalysis()
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
bool EYE_Mechanics::EyeMech::NRscheme()
{
    double residual,residual0,norm,norm0;
    string output;

    for(int ic=1;ic<=NRiteration;ic++){

      calcStressTensor();  //calc K and Q
      calcBoundaryForce();
      set_rhs_statics();

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

      if(isnan(norm)){
        cout << "residual is nan. Exit..." << endl;
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

      if(norm/norm0<NRtolerance) break;
      // if(test!=1 && ic>50) break;
    }
    return false;
}

// #################################################################
/**
 * @brief calc stress tensor
 */
void EYE_Mechanics::EyeMech::calcStressTensor()
{
  double elementVolume,volume;

  stress_tensor_initialize();

  volume = 0e0;

  #pragma omp parallel for
  for(int ic=0;ic<numOfElm;ic++){
    switch(element[ic].materialType){
      case MaterialType::LC:
        calcStressTensor_NeoHookean_element_spatialForm(ic,10e3,0.49e0,U,true);
      break;
      case MaterialType::SCLERA:
        calcStressTensor_NeoHookean_element_spatialForm(ic,100e3,0.49e0,U,true);
      break;
      case MaterialType::PRELAMINA:
        calcStressTensor_NeoHookean_element_spatialForm(ic,10000/2.98,0.49e0,U,true);
      break;
      default:
        cout << "undefined material. Exit.." << endl;
        exit(1);
    }
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
void EYE_Mechanics::EyeMech::set_rhs_statics()
{
 #pragma omp parallel for
  for(int ic=0;ic<numOfNode;ic++){
    for(int j=0;j<3;j++){
      RHS(ic,j) = (double)ibd(ic,j)*( - innerForce(ic,j)+boundaryForce(ic,j));
    }
  }
}


// #################################################################
/**
 * @brief calc boundary conditions
 * @param [in] stress
 */
void EYE_Mechanics::EyeMech::export_vtu(const string &file)
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

  fprintf(fp,"<PointData Vectors=\"displacement[m/s]\">\n");
  fprintf(fp,"<DataArray type=\"Float64\" Name=\"displacement[m/s]\" NumberOfComponents=\"3\" format=\"ascii\">\n");
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
