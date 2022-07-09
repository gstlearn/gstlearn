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
#include "Stats/ESelectivity.hpp"

class Db;

class GSTLEARN_EXPORT Selectivity
{
public:
  Selectivity(int ncut);
  Selectivity(const VectorDouble& zcuts);
  Selectivity(const Selectivity &m);
  Selectivity& operator= (const Selectivity &m);
  virtual ~Selectivity();

  static Selectivity* create(int ncut);
  static Selectivity* createByCodes(const std::vector<ESelectivity>& codes,
                                    bool flag_est,
                                    bool flag_std,
                                    bool flag_inter,
                                    int ncut,
                                    double proba = TEST,
                                    bool verbose = false);
  static Selectivity* createByKeys(const VectorString& scodes,
                                   bool flag_est,
                                   bool flag_std,
                                   bool flag_inter,
                                   int ncut,
                                   double proba = TEST,
                                   bool verbose = false);

  void   setZcut(int icut, double zcut);
  int    getNCuts() const { return _nCut; }
  double getZcut(int icut) const;
  int    getNQT() const { return ESelectivity::getSize(); }
  int    getVariableNumber() const;

  void defineRecoveries(const std::vector<ESelectivity>& codes,
                        bool flag_est,
                        bool flag_std,
                        bool flag_inter,
                        double proba = TEST,
                        bool verbose = false);

protected:
  bool _isValid(int icut) const;

private:
  void _printQTvars(const char *title, int type, int number) const;
  void _defineVariableRanks();

private:
  int _nCut;
  VectorDouble _Zcut;
  VectorInt _numberQTEst;
  VectorInt _numberQTStd;
  VectorInt _rankQTEst;
  VectorInt _rankQTStd;

};
