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

class Db;

class GSTLEARN_EXPORT Selectivity
{
public:
  Selectivity(int nclass = 0);
  Selectivity(const Selectivity &m);
  Selectivity& operator= (const Selectivity &m);
  virtual ~Selectivity();

  static Selectivity* createFromDb(const VectorDouble& cuts, const Db* db);
  static Selectivity* createFromTab(const VectorDouble& cuts,
                                    const VectorDouble& tab,
                                    const VectorDouble& weights = VectorDouble());

  void setBest(int iclass, double best);
  void setMest(int iclass, double mest);
  void setQest(int iclass, double qest);
  void setQstd(int iclass, double qstd);
  void setTest(int iclass, double test);
  void setTstd(int iclass, double tstd);
  void setZcut(int iclass, double zcut);

  int    getNClass() const { return _nClass; }
  double getBest(int iclass) const;
  double getMest(int iclass) const;
  double getQest(int iclass) const;
  double getQstd(int iclass) const;
  double getTest(int iclass) const;
  double getTstd(int iclass) const;
  double getZcut(int iclass) const;

  VectorDouble& getTest() { return _Test; }
  void calculateBenefitGrade();
  void dumpGini();
  void correctTonnageOrder();

private:
  bool _isValid(int iclass) const;

private:
  int _nClass;
  VectorDouble _Zcut;
  VectorDouble _Test;
  VectorDouble _Qest;
  VectorDouble _Best;
  VectorDouble _Mest;
  VectorDouble _Tstd;
  VectorDouble _Qstd;
};
