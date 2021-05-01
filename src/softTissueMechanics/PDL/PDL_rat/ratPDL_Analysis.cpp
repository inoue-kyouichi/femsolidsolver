/**
 * @file rigidElasticInteraction.cpp
 * @brief RigidElasticInteraction class
 * @author T. Otani
 */

#include "PDL_analysis.h"

using namespace std;

// #################################################################
/**
 * @brief rigid body-elastic body interaction problem
 */
void ratPDL::RigidElasticInteraction::RigidElasticInteraction::mainLoop()
{
  int output_iter=1;
  string output;

  //------------------------------------------------------------
  FILE *fp;
  string outputFile_tooth = outputDir + "/toothDisplacement.dat";
  if ((fp = fopen(outputFile_tooth.c_str(), "w")) == NULL) {
    exit(1);
  }
  fprintf(fp,"loop F Ux Uy Uz Qx Qy Qz ARMx ARMy ARMz\n");
  fclose(fp);
  //------------------------------------------------------------

  for(int j=0;j<3;j++) initialMomentArm[j] = FUpoint[j]-RBdy.xg[j];
  // for(int ic=0;ic<numOfElm;ic++) calcVolume_hexa(ic,volume0,8,2,false);

  //linear elastic material only
  ElasticBody.stress_tensor_initialize();
  for(int ic=0;ic<ElasticBody.numOfElm;ic++) ElasticBody.calcStressTensor_LinearElastic_element_spatialForm(ic,true);

  for(int loop=1;loop<=maxIteration;loop++){

    for(int j=0;j<3;j++) FU[j] = (double)loop/(double)maxIteration * FU_input[j];

    if(NRscheme()==1){
      printf("NR scheme is wrong. Loop=%d\n",loop);
      exit(1);
    }

    for(int i=0;i<ElasticBody.numOfNode;i++){
      for(int j=0;j<3;j++) ElasticBody.x(i,j) = ElasticBody.x0(i,j) + ElasticBody.U(i,j);
    }
    // for(int i=0;i<numOfElm;i++){
    //   calcVolume_hexa(i,volume,8,2,true);
    //   volumeChangeRatio(i)=volume(i)/volume0(i);
    // }


    for(int ic=0;ic<ElasticBody.numOfElm;ic++) ElasticBody.postProcess_LinearElastic_element_spatialForm(ic,true);
    // exportRestartData(loop);

    RBdy.updateShape();
    output = outputDir + "/tooth_"+to_string(dataNumber)+"_"+to_string(loop)+".ply";
    RBdy.exportPLY(output);

    double angle[3];
    mathTool::calc_thetaFromRotationMatrix(angle,RBdy.R);

    double momentArm[3];
    for(int i=0;i<3;i++){
      momentArm[i]=0e0;
      for(int j=0;j<3;j++) momentArm[i] += RBdy.R[i][j] * initialMomentArm[j];
      momentArm[i]-=initialMomentArm[i];
    }

    if ((fp = fopen(outputFile_tooth.c_str(), "a")) == NULL) {
      exit(1);
    }
    // double ForceMagnitude= (double)loop/(double)maxIteration * sqrt(FU_input[0]*FU_input[0]+FU_input[1]*FU_input[1]+FU_input[2]*FU_input[2]);
    // fprintf(fp,"%d %e %e %e %e %e %e %e %e %e %e\n",loop,ForceMagnitude,RBdy.U[0],RBdy.U[1],RBdy.U[2],angle[0],angle[1],angle[2],momentArm[0],momentArm[1],momentArm[2]);
    // fclose(fp);

    output = outputDir + "/PDL_"+to_string(dataNumber)+"_"+to_string(loop)+".vtu";
    fileIO::export_vtu_Mises(ElasticBody.x,ElasticBody.element,ElasticBody.numOfNode,ElasticBody.numOfElm,ElasticBody.U,ElasticBody.Mises,output);
    // fileIO::export_vtu(x,element,numOfNode,numOfElm,U,volumeChangeRatio,lambda_ave,sigmaEigen_Ave,sigmaEigenVector_Ave,output);
  }
}

// #################################################################
/**
 * @brief fem solid analysis routine
 */
int ratPDL::RigidElasticInteraction::NRscheme()
{
  double residual,residual0,norm,norm0;
  string output;

  for(int ic=1;ic<=NRiteration;ic++){

    // ElasticBody.calcStressTensor();  //calc K and Q

    //linear elastic material only
    #pragma omp parallel for
    for(int i=0;i<ElasticBody.numOfNode;i++){
      for(int j=0;j<3;j++) ElasticBody.innerForce(i,j) = 0e0;
    }

    for (int ia = 0; ia < ElasticBody.numOfElm; ia++){
      for (int p = 0; p < ElasticBody.element[ia].node.size(); p++){
       for (int q = 0; q < ElasticBody.element[ia].node.size(); q++){
        for (int i = 0; i < 3; i++){
          for (int j = 0; j < 3; j++) ElasticBody.innerForce(ElasticBody.element[ia].node[p],i) += ElasticBody.Ku[ia](p,q,i,j) * ElasticBody.U(ElasticBody.element[ia].node[q],j);
        }
      }
    }
  }

    calcTemporalFw(RBdy);

    //elastic body-rigid body interaction
    for(int i=0;i<numOfCP;i++){
      for(int j=0;j<3;j++) ElasticBody.innerForce(CP(i),j)+=LAMBDA(i,j);
    }

    ElasticBody.set_rhs_statics();
    PARDISO.set_CSR_value3D(ElasticBody.Ku,ElasticBody.element,ElasticBody.numOfNode,ElasticBody.numOfElm,ElasticBody.inb);

    //-----rigid body interaction term-------
    calcRigidBodyInteractionTerm(ElasticBody.U,RBdy);
    PARDISO.set_CSR_value_rigidBodyInteraction(ElasticBody.numOfNode,iCP,Rb,Kqq,numOfCP);
    PARDISO.set_CSR_dirichlet_boundary_condition3D(ElasticBody.numOfNode,ElasticBody.ibd);

    //rhs term
    for(int i=0;i<ElasticBody.numOfNode;i++){
      for(int j=0;j<3;j++) PARDISO.b[i+j*ElasticBody.numOfNode]=ElasticBody.RHS(i,j);
    }
    for(int i=0;i<numOfCP;i++){
      for(int j=0;j<3;j++) PARDISO.b[3*ElasticBody.numOfNode+i+j*numOfCP]=-Qlambda(i,j);
    }
    for(int j=0;j<3;j++){
      PARDISO.b[3*ElasticBody.numOfNode+3*numOfCP+j]  = FU[j]-QU[j];
      PARDISO.b[3*ElasticBody.numOfNode+3*numOfCP+3+j]= Fw[j]-Qw[j];
    }

    PARDISO.main(3*ElasticBody.numOfNode+3*numOfCP+6,OMPnumThreads);
    norm = PARDISO.vector_norm(3*ElasticBody.numOfNode+3*numOfCP+6,PARDISO.x);

    if(isnan(norm)){
      cout << "norm is nan " << endl;
      return 1;
    }
    if(ic==1) norm0 = norm;
    corrector_statics(ElasticBody.U,PARDISO.x,RBdy,ElasticBody.numOfNode,relaxation);

    for(int i=0;i<ElasticBody.numOfNode;i++){
      for(int j=0;j<3;j++) ElasticBody.x(i,j) = ElasticBody.x0(i,j) + ElasticBody.U(i,j);
    }

    //for debug
    // output = outputDir + "/test_NR_" + to_string(ic) + ".vtu";
    // fileIO::export_vtu(x,element,numOfNode,numOfElm,U,volumeChangeRatio,lambda_ave,output);

    //for debug
    // RBdy.updateShape();
    // output = outputDir + "/rigidBody_NR_" + to_string(ic) + ".ply";
    // RBdy.exportPLY(output);
    printf("NR iter.=%d norm/norm0=%e\n",ic,norm/norm0);
    //printf("NR iter=%d Tooth displacement=(%e %e %e)\n",ic,RBdy.U[0],RBdy.U[1],RBdy.U[2]);

    if(norm/norm0<NRtolerance) break;
    // if(test!=1 && ic>50) break;
  }
  return 0;
}

// #################################################################
/**
 * @brief calc b0
 * @param [in] RBdy          rigid body class
 */
void ratPDL::RigidElasticInteraction::initialize_rigidBodyInteraction()
{
  ElasticBody.initialize(tp);
  ElasticBody.inputMaterialInfo(tp);
  RBdy.initialize(tp);
  inputSolverInfo(tp);
  inputOutputInfo(tp);
  inputRigidBodyInterface();

  LAMBDA.allocate(numOfCP,3);
  Rb.allocate(numOfCP,3,3);
  b0.allocate(numOfCP,3);
  b.allocate(numOfCP,3);
  Qlambda.allocate(numOfCP,3);

  for(int i=0;i<numOfCP;i++){
    for(int j=0;j<3;j++) LAMBDA(i,j)=0e0;
  }

  for(int ic=0;ic<numOfCP;ic++){
    for(int j=0;j<3;j++) b0(ic,j)=ElasticBody.x(CP(ic),j)-RBdy.xg[j];
  }

  omp_set_num_threads(OMPnumThreads);

  mkdir(outputDir.c_str(),S_IRWXU | S_IRWXG | S_IRWXO);

  for(int i=0;i<ElasticBody.numOfNode;i++){
      if(ElasticBody.ibd(i,0)==0 && iCP(i)!=-1) cout << i << endl;
  }

  //CSR setting
  PARDISO.initialize(3*ElasticBody.numOfNode,3*numOfCP);
  PARDISO.CSR_initialize(ElasticBody.inb,ElasticBody.numOfNode,iCP,CP,numOfCP,3);

  string output = outputDir + "/" + fileName + "_boundary" + ".vtu";
  fileIO::export_vtu_boundary(ElasticBody.x,ElasticBody.element,ElasticBody.numOfNode,ElasticBody.numOfElm,ElasticBody.ibd,ElasticBody.bd,output);
  output = outputDir + "/start.ply";
  RBdy.exportPLY(output);
}


// #################################################################
/**
 * @brief domain information from tp file
 */
void ratPDL::RigidElasticInteraction::inputRigidBodyInterface()
{
  string str,base_label,label,inputDir;

  base_label = "/Domain";

  label = base_label + "/inputDir";
  if ( !tp.getInspectedValue(label,inputDir)){
    cout << "data format is not set" << endl;
    exit(0);
  }
  string file;
  label = base_label + "/interface";
  if ( !tp.getInspectedValue(label, file)){
    cout << label << " is not found" << endl;
    exit(0);
  }

  file=inputDir + "/" + file;

  numOfCP = fileIO::CountNumbersOfTextLines(file);
  CP.allocate(numOfCP);
  iCP.allocate(ElasticBody.numOfNode);

  FILE *fp;
  if ((fp = fopen(file.c_str(), "r")) == NULL) {
    cout << "file open error" << endl;
    exit(1);
  }
  for(int i=0;i<numOfCP;i++) fscanf(fp,"%d\n",&CP(i));
  fclose(fp);

  for(int i=0;i<ElasticBody.numOfNode;i++) iCP(i)=-1;
  for(int i=0;i<numOfCP;i++) iCP(CP(i))=i;

  base_label = "/RigidBody";
  // label = base_label + "/Force";
  // if ( !tp.getInspectedVector(label,FU_input,3)){
  //   cout << label << " is not found" << endl;
  //   exit(0);
  // }

  label = base_label + "/ForcePoint";
  if ( !tp.getInspectedVector(label,FUpoint,3)){
    cout << label << " is not found" << endl;
    exit(0);
  }
  // label = base_label + "/Moment";
  // if ( !tp.getInspectedVector(label,Fw,3)){
  //   cout << label << " is not found" << endl;
  //   exit(0);
  // }
}

// #################################################################
/**
 * @brief solver information from TP file
 */
void ratPDL::RigidElasticInteraction::inputSolverInfo(TextParser &tp)
{
  string str,base_label,label;
  int tmp;

  base_label = "/Solver";

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
 * @brief output information from TP file
 */
void ratPDL::RigidElasticInteraction::inputOutputInfo(TextParser &tp)
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
 * @brief assembly right-hand side
 */
void ratPDL::RigidElasticInteraction::set_rhs_statics()
{
 // #pragma omp parallel for
  for(int ic=0;ic<ElasticBody.numOfNode;ic++){
    for(int j=0;j<3;j++){
      ElasticBody.RHS(ic,j) = (double)ElasticBody.ibd(ic,j)*( - ElasticBody.innerForce(ic,j));
    }
  }
}