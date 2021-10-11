/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "Model/Convolution.hpp"
#include "Model/Tapering.hpp"
#include "Model/EModelProperty.hpp"
#include "Anamorphosis/Anam.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"

class ModTrans : public AStringable
{
public:
  ModTrans();
  ModTrans(const ModTrans &m);
  ModTrans& operator= (const ModTrans &m);
  virtual ~ModTrans();

  void cancelProperty();
  int addConvolution(int conv_type,
                     int conv_idir,
                     int conv_ndisc,
                     double conv_range);
  int addAnamorphosis(int anam_type,
                      int anam_nclass,
                      int anam_iclass,
                      int anam_var,
                      double anam_coefr,
                      double anam_coefs,
                      VectorDouble& anam_strcnt,
                      VectorDouble& anam_stats);
  int addTapering(int tape_type,double tape_range);

  virtual String toString(int level = 0) const override;

  const EModelProperty& getModTransMode() const { return _modTransMode; }
  int getAnamIClass() const                     { return _anamIClass; }
  int getAnamNClass() const                     { return _anamNClass; }
  int getAnamPointBlock() const                 { return _anamPointBlock;  }
  void setAnamIClass(int iclass)                { _anamIClass = iclass; }
  void setAnamVar(int var)                      { _anamPointBlock = var; }
  const VectorDouble& getAnamStrCount() const   { return _anamStrCount; }
  const VectorDouble& getAnamMeans() const      { return _anamMeans; }
  double getAnamMeans(int iclass) const         { return _anamMeans[iclass]; }

  Anam* getAnam()        const { return _anam; }
  Convolution* getConv() const { return _conv; }
  Tapering* getTape()    const { return _tape; }

private:
  EModelProperty _modTransMode;
  Convolution* _conv;
  Tapering*    _tape;
  Anam*        _anam;
  int    _anamIClass;         /* Target factor (-1: discretized grade) */
  int    _anamNClass;         /* Number of indicator classes */
  int    _anamPointBlock;     /* Type of point / block covariance */
  VectorDouble _anamStrCount; /* Array of structure count per model (IR)  */
  VectorDouble _anamMeans;    /* Array of statistics per class */
};
