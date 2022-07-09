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
  static Selectivity* create(const VectorDouble& zcut);
  static Selectivity* createByCodes(const std::vector<ESelectivity>& codes,
                                    const VectorDouble& zcuts,
                                    bool flag_est,
                                    bool flag_std,
                                    bool flag_inter,
                                    double proba = TEST,
                                    bool verbose = false);
  static Selectivity* createByKeys(const VectorString& scodes,
                                   const VectorDouble& zcuts,
                                   bool flag_est,
                                   bool flag_std,
                                   bool flag_inter,
                                   double proba = TEST,
                                   bool verbose = false);

  void   setZcut(int icut, double zcut);
  int    getNCuts() const { return _nCut; }
  double getZcut(int icut) const;
  const VectorDouble& getZcut() const { return _Zcut; }
  int    getNQT() const { return ESelectivity::getSize(); }
  int    getVariableNumber() const;

  void defineRecoveries(const std::vector<ESelectivity>& codes,
                        bool flag_est,
                        bool flag_std,
                        bool flag_inter,
                        double proba = TEST,
                        bool verbose = false);

  bool isUsed(const ESelectivity& code) const;
  bool isUsedEst(const ESelectivity& code) const;
  bool isUsedStD(const ESelectivity& code) const;
  bool isNeededT() const;
  bool isNeededQ() const;
  int  getAddressQTEst(const ESelectivity& code, int iptr0, int rank=0) const;
  int  getAddressQTStD(const ESelectivity& code, int iptr0, int rank=0) const;
  int  getNumberQTEst(const ESelectivity& code) const { return _numberQTEst[code.getValue()]; }
  int  getNumberQTStd(const ESelectivity& code) const { return _numberQTStd[code.getValue()]; }
  const VectorInt getNumberQTEst() const { return _numberQTEst; }
  const VectorInt getNumberQTStd() const { return _numberQTStd; }

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
