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
#include "Covariances/CovLMC.hpp"
#include "Matrix/MatrixRectangular.hpp"

class ASpace;
class SpacePoint;
class CovAniso;
class Model;

typedef struct
{
  std::string convName;
  double convScale;
  double (*convFunc)(double);
} Def_Convolution;

/* Prototyping the internal covariance functions */
double _conv_uniform(double v);
double _conv_exponential(double v);
double _conv_gaussian(double v);
double _conv_sincard(double v);

GSTLEARN_EXPORT Def_Convolution& D_CONV(int rank);

class GSTLEARN_EXPORT CovLMCConvolution : public CovLMC
{
public:
  CovLMCConvolution(int conv_type,
                    int conv_dir,
                    double conv_range,
                    int conv_ndisc = 10,
                    const ASpace* space = nullptr);
  CovLMCConvolution(const CovLMCConvolution &r);
  CovLMCConvolution& operator= (const CovLMCConvolution &r);
  virtual ~CovLMCConvolution();

  virtual IClonable* clone() const override { return new CovLMCConvolution(*this); };
  virtual String toString(int /*level*/) const;

  virtual double eval0(int ivar,
                       int jvar,
                       const CovCalcMode& mode = CovCalcMode()) const override;
  virtual double eval(int ivar,
                      int jvar,
                      const SpacePoint& p1,
                      const SpacePoint& p2,
                      const CovCalcMode& mode = CovCalcMode()) const override;

  int init(int conv_type, int conv_idir, double conv_range, int conv_ndisc);

  int getConvDir() const { return _convDir; }
  int getConvDiscNumber() const { return _convDiscNumber; }
  double getConvRange() const { return _convRange; }
  int getConvType() const { return _convType; }
  const VectorDouble& getConvWeight() const { return _convWeight; }
  const MatrixRectangular& getConvIncr() const { return _convIncr; }
  VectorDouble getConvIncr(int rank) const { return _convIncr.getColumn(rank); }
  int getConvNumber() const { return _convNumber; }

private:
  int _getConvNumber();

private:
  int _convType; /* Convolution type */
  int _convDir; /* Convolution direction: 0:X, 1:Y, 2:Z, 3:XY, 4:XYZ */
  int _convDiscNumber; /* Number of discretization per direction */
  double _convRange; /* Convolution Range */
  int _convNumber;
  MatrixRectangular _convIncr; /* Discretization lags */
  VectorDouble      _convWeight; /* Weights for convolution */
};
