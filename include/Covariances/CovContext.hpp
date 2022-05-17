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
#include "Basic/Vector.hpp"
#include "Space/ASpaceObject.hpp"

class ASpace;
class Vario;
class Db;

class GSTLEARN_EXPORT CovContext : public ASpaceObject
{
public:
  CovContext(int nvar = 1,
             const ASpace* space = nullptr,
             double field = 1);
  CovContext(int nvar,
             int ndim,
             double field = 1,
             const VectorDouble& mean = VectorDouble(),
             const VectorDouble& covar0 = VectorDouble());
  CovContext(const Db *db,
             const ASpace* space = nullptr);
  CovContext(const Vario* vario,
             const ASpace* space = nullptr);
  CovContext(const CovContext &r);
  CovContext& operator= (const CovContext &r);
  virtual ~CovContext();

  /// AStringable interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Indicate if I am consistent with the provided space
  virtual bool isConsistent(const ASpace* space) const override;

  static CovContext* create(int nvar, int ndim, double field = 1.);

  bool isEqual(const CovContext &r) const;

  int                 getNVar()         const { return _nVar; }
  double              getField()        const { return _field; }
  const VectorDouble& getMean()         const { return _mean; }
  const VectorDouble& getCovar0()       const { return _covar0; }

  double getMean(int ivar) const;
  double getCovar0(int ivar, int jvar) const;

  void setNVar(int nvar)                 { _nVar = nvar; _update(); }
  void setField(double field)            { _field = field; }

  void setMean(const VectorDouble& mean);
  void setMean(int ivar, const double mean);

  void setCovar0(const VectorDouble& covar0);
  void setCovar0(int ivar, int jvar, double covar0);

  void copyCovContext(const CovContext& ctxt);

private:
  int           _nVar;         /*! Number of variables */
  double        _field;        /*! Field maximum size */
  VectorDouble  _mean;         /*! Array of Variable Mean */
  VectorDouble  _covar0;       /*! Variance-Covariance matrix (used for covariances) */

private:
  int _getIndex(int ivar, int jvar) const;
  void _update();
};
