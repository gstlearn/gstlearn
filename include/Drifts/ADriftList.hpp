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

#include "Drifts/ADrift.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Drifts/EDrift.hpp"

#include "Basic/Vector.hpp"
#include "Basic/IClonable.hpp"

class ASpace;
class SpacePoint;
class Db;

class ADriftList : public ADrift, public IClonable
{
public:
  ADriftList(bool flagLinked = false, const ASpace* space = nullptr);
  ADriftList(const ADriftList &r);
  ADriftList& operator= (const ADriftList &r);
  virtual ~ADriftList();

  ///////////////////////////////////////////////////
  /// IClonable interface */
  virtual IClonable* clone() const override;
  ///////////////////////////////////////////////////

  ///////////////////////////////////////////////////
  /// ASpaceObject Interface
  virtual bool isConsistent(const ASpace* space) const override;
  ///////////////////////////////////////////////////

  ///////////////////////////////////////////////////
  /// AStringable Interface
  virtual String toString(int level = 0) const override;
  ///////////////////////////////////////////////////

  virtual double eval(const Db* db, int iech) const override;
  int getNVariables() const override;

  int getDriftNumber() const { return static_cast<int>(_drifts.size()); }

  // Add an elementary drift structure
  void addDrift(const ADriftElem* drift);
  // Remove an elementary drift structure
  void delDrift(unsigned int i);
  // Remove all elementary drift structures
  void delAllDrift();

  const std::vector<ADriftElem*>& getDrifts() const { return _drifts; }
  void setDrifts(const std::vector<ADriftElem*>& drifts) { _drifts = drifts; }
  const VectorBool& getFiltered() const { return _filtered; }
  void setFiltered(const VectorBool& filtered) { _filtered = filtered; }
  bool isFiltered(int i) const;
  void setFiltered(int i, bool filter);
  int  getDriftEquationNumber() const { return _driftEquationNumber; }

  /// TODO : to be removed (encapsulation)
  ////////////////////////////////////////////////
  const ADriftElem*  getDrift(int il) const;
  ADriftElem*        getDrift(int il); /// beurk :(
  const EDrift&      getType(int il) const;
  int                getRankFex(int il) const;
  String             getDriftName(int il) const;
  void               setType(int il, const EDrift& type);
  ////////////////////////////////////////////////

  const VectorDouble& getCoefDrift() const { return _coefDrift; }

  /**
   *
   * @param ivar Rank of the variable (_nVar)
   * @param il Rank of the drift function
   * @param ib Rank of the drift equation (_driftEquationNumber)
   * @return
   */
  double getCoefDrift(int ivar, int il, int ib) const { return _coefDrift[_getAddress(ivar,il,ib)]; }
  void setCoefDrift(int ivar, int il, int ib, double value) { _coefDrift[_getAddress(ivar,il,ib)] = value; }
  void setCoefDrift(int rank, double coeff) { _coefDrift[rank] = coeff; }

  VectorDouble getDrift(const Db* db, int ib, bool useSel = true);
  VectorVectorDouble getDrifts(const Db* db, bool useSel = true);

private:
  bool _isDriftIndexValid(int i) const;
  bool _isDriftEquationValid(int ib) const;
  int  _getAddress(int ivar, int il, int ib) const
  {
    return (ib + getDriftEquationNumber() * (il + getDriftNumber() * ivar));
  }
  void _updateCoefDrift();
  void _setDriftEquationNumber();

#ifndef SWIG
protected:
  bool _flagLinked;
  int  _driftEquationNumber;
  VectorDouble             _coefDrift; /* Array of Drift Coefficients */
  std::vector<ADriftElem*> _drifts;    /* Vector of elementary drift functions */
  VectorBool               _filtered;  /* Vector of filtered flags */
#endif
};
