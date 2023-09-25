/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES Paris / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Drifts/ADrift.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/VectorHelper.hpp"


class ASpace;
class SpacePoint;
class Db;

class GSTLEARN_EXPORT DriftList : public ADrift, public ICloneable
{
public:
  DriftList(const ASpace* space = nullptr);
  DriftList(const DriftList &r);
  DriftList& operator= (const DriftList &r);
  virtual ~DriftList();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftList)

  /// ASpaceObject Interface
  virtual bool isConsistent(const ASpace* space) const override;

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ADrift interface */
  virtual double eval(const Db* db, int iech) const override;
  int getNVariables() const override;

  int getDriftNumber() const { return static_cast<int>(_drifts.size()); }

  // Add one elementary drift structure
  void addDrift(const ADriftElem* drift);
  // Remove an elementary drift structure
  void delDrift(unsigned int i);
  // Remove all elementary drift structures
  void delAllDrifts();

  const VectorBool& getFiltered() const { return _filtered; }
  void setFiltered(const VectorBool& filtered) { _filtered = filtered; }
  bool isFiltered(int i) const;
  void setFiltered(int i, bool filter);
  int  getDriftEquationNumber() const;
  bool hasExternalDrift() const;
  bool isValid() const;

  /// TODO : to be removed (encapsulation)
  ////////////////////////////////////////////////
  const ADriftElem*  getDrift(int il) const;
  ADriftElem*        getDrift(int il); /// beurk :(
  int                getRankFex(int il) const;
  String             getDriftName(int il) const;
  ////////////////////////////////////////////////

  const VectorDouble& getDriftCoef() const { return _driftCoef; }

  /**
   *
   * @param ivar Rank of the variable (_nVar)
   * @param il Rank of the drift function
   * @param ib Rank of the drift equation (_driftEquationNumber)
   * @return
   */
  double getDriftCoef(int ivar, int il, int ib) const { return _driftCoef[_getAddress(ivar,il,ib)]; }
  void setDriftCoef(int ivar, int il, int ib, double value) { _driftCoef[_getAddress(ivar,il,ib)] = value; }
  void resetDriftCoeff() { VectorHelper::fill(_driftCoef, 0.); }
  VectorDouble getDriftCoefByPart(int ivar, int ib) const;
  void setDriftCoefByPart(int ivar, int ib, const VectorDouble& coef);

  double getDrift(const Db* db, int ib, int iech) const;
  VectorDouble getDriftByColumn(const Db* db, int ib, bool useSel = true) const;
  VectorDouble getDriftBySample(const Db* db, int iech) const;
  VectorVectorDouble getDrifts(const Db* db, bool useSel = true) const;
  bool isFlagLinked() const { return _flagLinked; }
  VectorDouble evalDrifts(const Db* db,
                          const VectorDouble& coeffs,
                          bool useSel = false) const;
  int getDriftMaxIRFOrder(void) const;
  bool isDriftDefined(const VectorInt &powers, int rank_fex = 0) const;
  bool isDriftDifferentDefined(const VectorInt &powers, int rank_fex = -1) const;

  void copyCovContext(const CovContext& ctxt);

  void setFlagLinked(bool flagLinked) { _flagLinked = flagLinked; }

  void updateDriftList();

private:
  bool _isDriftIndexValid(int i) const;
  bool _isDriftEquationValid(int ib) const;
  int  _getAddress(int ivar, int il, int ib) const
  {
    return (ib + getDriftEquationNumber() * (il + getDriftNumber() * ivar));
  }

#ifndef SWIG
protected:
  bool _flagLinked;
  VectorDouble             _driftCoef; /* Array of Drift Coefficients */
  std::vector<ADriftElem*> _drifts;    /* Vector of elementary drift functions */
  VectorBool               _filtered;  /* Vector of filtered flags (Dimension: as _drifts) */
#endif
};
