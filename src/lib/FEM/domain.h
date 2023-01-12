#ifndef _DOMAIN_H_
#define _DOMAIN_H_

//##################################################################################
//
// FEM Base
//
// Copyright (c) 2016 Mechanical and Bioengineering Systems Lab.,
//                    Graduate School of Engineering Science and Bioengineering,
//                    Osaka University.
// All rights reserved.
//
//##################################################################################

/**
 * @file   domain.h
 * @brief  FEMBase Definition Header
 * @author T. Otani
 */

#include <string>
#include "fem_define.h"
#include "allocation.h"

class Domain{
public:
  int numOfNode, numOfElm, numOfDirichlet,numOfFD, numOfNode_hb, numOfElm_hb;
  int numOfNeumann;
  int meshType;
  ARRAY2D<double> x,x0,x_hb;
  std::vector<ElementType> element,belement,element_hb;
  ARRAY2D<int> ibd;
  ARRAY2D<double> bd;
  ARRAY2D<double> bn;
  VECTOR2D<int> inb;
  VECTOR2D<int> ieb;
  ARRAY2D<int> FD;


  void set_geometry(const std::string &file1,const std::string &file2,const std::string &file3);
  void set_dirichlet(const std::string &D_file);
  void set_dirichlet(const std::string &D_file,const std::string &Dvalue_file);
  void set_neumann(const std::string &N_file,const std::string &meshTypeFile,const std::string &NvalueFile);
  void export_vtk(const std::string &file);

private:
  void calc_adjacent_nodes();
  void calc_adjacent_elements();
private:

};

#endif //_DOMAIN_H_
