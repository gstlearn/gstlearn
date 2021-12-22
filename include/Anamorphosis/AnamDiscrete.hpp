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

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Anamorphosis/Anam.hpp"
#include "Matrix/MatrixRectangular.hpp"

class GSTLEARN_EXPORT AnamDiscrete: public Anam
{
public:
  AnamDiscrete(const EAnam& type = EAnam::UNDEFINED);
  AnamDiscrete(const AnamDiscrete &m);
  AnamDiscrete& operator= (const AnamDiscrete &m);
  virtual ~AnamDiscrete();

  virtual String toString(int level) const override;
  virtual VectorDouble z2f(int nfact, const VectorInt& ifacs, double z) const = 0;
  virtual void calculateMeanAndVariance();

  double getMean() const { return _mean; }
  int getNCut() const { return _nCut; }
  int getNClass() const { return _nCut + 1; }
  int getNElem() const { return _nElem; }
  const VectorDouble& getZCut() const { return _zCut; }
  double getZCut(int i) const { return _zCut[i]; }
  double getVariance() const { return _variance; }

  void setMean(double mean) { _mean = mean; }
  void setVariance(double variance) { _variance = variance; }
  virtual void setNCut(int ncut);
  virtual void setZCut(const VectorDouble& zcut);
  void setStats(const VectorDouble& stats);

  // Function for using Stats in DD anamorphosis
  double getDDStatProp  (int iclass) const;
  double getDDStatZmoy  (int iclass) const;
  double getDDStatCnorm (int iclass) const;
  double getDDStatLambda(int iclass) const;
  double getDDStatU     (int iclass) const;
  double getDDStatMul   (int iclass) const;
  void   setDDStatProp  (int iclass, double value);
  void   setDDStatZmoy  (int iclass, double value);
  void   setDDStatCnorm (int iclass, double value);
  void   setDDStatLambda(int iclass, double value);
  void   setDDStatU     (int iclass, double value);
  void   setDDStatMul   (int iclass, double value);

  // Function for using Stats in IR anamorphosis
  double getIRStatT   (int iclass) const;
  double getIRStatQ   (int iclass) const;
  double getIRStatZ   (int iclass) const;
  double getIRStatB   (int iclass) const;
  double getIRStatR   (int iclass) const;
  double getIRStatRV  (int iclass) const;
  void setIRStatT(int iclass, double value);
  void setIRStatQ(int iclass, double value);
  void setIRStatZ(int iclass, double value);
  void setIRStatB(int iclass, double value);
  void setIRStatR(int iclass, double value);
  void setIRStatRV(int iclass, double value);

  const MatrixRectangular& getStats() const { return _stats; }

private:
  bool _isClassValid(int iclass) const;
  void _resize();

private:
  int _nCut;
  int _nElem;
  double _mean;
  double _variance;
  VectorDouble _zCut;
  MatrixRectangular _stats;
};
