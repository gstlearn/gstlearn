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

#include "Enum/ECalcMember.hpp"
#include "Drifts/ADrift.hpp"
#include "Basic/ICloneable.hpp"
#include "Covariances/CovContext.hpp"
#include "Matrix/MatrixRectangular.hpp"

class ASpace;
class SpacePoint;
class Db;
class ELoc;

/**
 * \brief
 * This class provides the information on **Drift** part of the Model. The drift plays the role of the average of the
 * target Random Function which may be constant or vary as a function with low frequency variations (by opposition to the
 * complementary part of the Spatial Characteristics which is described by its Covariance)
 *
 * This class essentially contains a list of basic (active)e drift functions: see ADrift.hpp for details.
 *
 * This class also carry other important informations:
 * - a vector giving the status of each basic drift functions: it may be *active* or *filtered*
 * - some additional information defining some relationship between the basic drift function: this is used for the special
 * case where the different Random Functions obey to algebraic relations.
 */
class GSTLEARN_EXPORT DriftList : public AStringable, public ICloneable
{
public:
  DriftList(const CovContext &ctxt = CovContext());
  DriftList(const DriftList &r);
  DriftList& operator= (const DriftList &r);
  virtual ~DriftList();

  /// ICloneable interface
  IMPLEMENT_CLONING(DriftList)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int getNVariables() const { return _ctxt.getNVar(); }
  int getDriftNumber() const { return static_cast<int>(_drifts.size()); }
  bool hasDrift() const { return !_drifts.empty(); }

  // Add one elementary drift structure
  void addDrift(const ADrift* drift);
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
  const ADrift*  getDrift(int il) const;
  int            getRankFex(int il) const;
  String         getDriftName(int il) const;
  ////////////////////////////////////////////////

  const VectorDouble& getBetaHats() const { return _betaHat; }

  void setDriftCLByPart(int ivar, int ib, const VectorDouble& coef);
  void resetDriftList();

  bool isDriftSampleDefined(const Db *db,
                            int ib,
                            int nech,
                            const VectorInt &nbgh,
                            const ELoc &loctype) const;
  double getDrift(const Db* db, int ib, int iech) const;

  VectorVectorDouble getDrifts(const Db* db, bool useSel = true) const;
  bool isFlagLinked() const { return _flagLinked; }
  bool isFlagCombined() const { return _flagCombined; }
  int getDriftMaxIRFOrder(void) const;
  bool isDriftDefined(const VectorInt &powers, int rank_fex = 0) const;
  bool isDriftDifferentDefined(const VectorInt &powers, int rank_fex = -1) const;

  void copyCovContext(const CovContext& ctxt) { _ctxt.copyCovContext(ctxt); }

  void setFlagLinked(bool flagLinked) { _flagLinked = flagLinked; }
  void setFlagCombined(bool flagCombined) { _flagCombined = flagCombined; }
  void setBetaHat(const VectorDouble &betaHat) { _betaHat = betaHat; }

  double evalDrift(const Db* db,
                   int iech,
                   int il,
                   const ECalcMember& member = ECalcMember::fromKey("LHS")) const;
  double evalDriftCoef(const Db *db, int iech, const VectorDouble &coeffs) const;
  VectorDouble evalDriftCoefs(const Db *db,
                              const VectorDouble &coeffs,
                              bool useSel = false) const;
  VectorDouble evalDriftBySample(const Db *db,
                                 int iech,
                                 const ECalcMember &member = ECalcMember::fromKey("LHS")) const;
  void evalDriftBySampleInPlace(const Db *db,
                                int iech,
                                const ECalcMember &member,
                                VectorDouble &drftab) const;
  MatrixRectangular evalDriftMatrix(const Db *db,
                                    int ivar0 = -1,
                                    const VectorInt &nbgh = VectorInt(),
                                    const ECalcMember &member = ECalcMember::fromKey("LHS"));
  double evalDriftValue(const Db *db,
                        int iech,
                        int ivar,
                        int ib,
                        const ECalcMember &member = ECalcMember::fromKey("LHS")) const;

private:
  bool _isDriftIndexValid(int i) const;
  bool _isDriftEquationValid(int ib) const;
  int  _getAddress(int ivar, int il, int ib) const
  {
    return (ib + getDriftEquationNumber() * (il + getDriftNumber() * ivar));
  }
  double _getDriftCL(int ivar, int il, int ib) const { return _driftCL[_getAddress(ivar,il,ib)]; }
  void   _setDriftCL(int ivar, int il, int ib, double value) { _driftCL[_getAddress(ivar,il,ib)] = value; }
  VectorInt _getActiveVariables(int ivar0) const;

#ifndef SWIG
protected:
  bool                 _flagLinked;   /* True when Drift equations are linked */
  bool                 _flagCombined; /* True when Drift is obtained by Linear Combination */
  VectorDouble         _driftCL;      /* Linear combination of Drift Coefficients */
  std::vector<ADrift*> _drifts;       /* Vector of elementary drift functions */
  VectorDouble         _betaHat;      /* Drift coefficients by ML */
  VectorBool           _filtered;     /* Vector of filtered flags (Dimension: as _drifts) */
  CovContext           _ctxt;         /* Context (space, number of variables, ...) */
#endif
};
