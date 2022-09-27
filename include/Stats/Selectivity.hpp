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

#include "Enum/ESelectivity.hpp"

#include "Basic/Vector.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/AStringable.hpp"
#include "Matrix/MatrixInt.hpp"
#include "Basic/Table.hpp"

class Db;
class AAnam;

class GSTLEARN_EXPORT Selectivity: public ICloneable, public AStringable
{
public:
  Selectivity(int ncut = 0);
  Selectivity(const VectorDouble &zcuts,
              double zmax = TEST,
              double proba = TEST,
              bool flag_tonnage_correct = false);
  Selectivity(const Selectivity &m);
  Selectivity& operator= (const Selectivity &m);
  virtual ~Selectivity();

  /// ICloneable interface
  IMPLEMENT_CLONING(Selectivity)

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static Selectivity* create(int ncut);
  static Selectivity* create(const VectorDouble& zcut);
  static Selectivity* createByCodes(const std::vector<ESelectivity>& codes,
                                    const VectorDouble& zcuts,
                                    bool flag_est,
                                    bool flag_std,
                                    double proba = TEST,
                                    bool verbose = false);
  static Selectivity* createByKeys(const VectorString& keys,
                                   const VectorDouble& zcuts,
                                   bool flag_est,
                                   bool flag_std,
                                   double zmax = TEST,
                                   bool flag_correct = false,
                                   double proba = TEST,
                                   bool verbose = false);
  static Selectivity* createInterpolation(const VectorDouble& zcuts,
                                          const Selectivity& selecin,
                                          bool verbose);

  int calculateFromDb(const Db* db, bool autoCuts = false);
  int calculateFromArray(const VectorDouble& tab,
                         const VectorDouble& weights = VectorDouble(),
                         bool autoCuts = false);
  int calculateFromAnam(AAnam* anam);

  const Table eval(const Db *db, bool autoCuts = false);
  const Table eval(const VectorDouble &tab,
                   const VectorDouble &weights = VectorDouble(),
                   bool autoCuts = false);
  const Table eval(AAnam *anam);

  void   resetCuts(const VectorDouble& zcuts);
  int    getNCuts() const { return static_cast<int>(_Zcut.size()); }
  int    getNQT() const { return static_cast<int>(ESelectivity::getSize()); }
  int    getVariableNumber() const;
  String getVariableName(const ESelectivity& code, int icut, int mode) const;
  String getVariableName(int rank0) const;
  VectorString getVariableNames() const;

  void   setZcut(int icut, double zcut);
  void   setBest(int icut, double best);
  void   setMest(int icut, double mest);
  void   setQest(int icut, double qest);
  void   setQstd(int icut, double qstd);
  void   setTest(int icut, double test);
  void   setTstd(int icut, double tstd);

  double getZcut(int icut) const;
  double getBest(int icut) const;
  double getMest(int icut) const;
  double getQest(int icut) const;
  double getQstd(int icut) const;
  double getTest(int icut) const;
  double getTstd(int icut) const;
  const VectorDouble& getZcut() const { return _Zcut; }

  void calculateBenefitAndGrade();
  void dumpGini();
  void correctTonnageOrder();
  void defineRecoveries(const std::vector<ESelectivity>& codes,
                        bool flag_est,
                        bool flag_std,
                        double proba = TEST,
                        bool verbose = false);

  bool isUsed(const ESelectivity& code) const;
  bool isUsedEst(const ESelectivity& code) const;
  bool isUsedStD(const ESelectivity& code) const;
  bool isNeededT() const;
  bool isNeededQ() const;
  int  getAddressQTEst(const ESelectivity& code, int iptr0, int rank=0) const;
  int  getAddressQTStd(const ESelectivity& code, int iptr0, int rank=0) const;
  int  getNumberQTEst(const ESelectivity& code) const;
  int  getNumberQTStd(const ESelectivity& code) const;
  const VectorInt getNumberQTEst() const;
  const VectorInt getNumberQTStd() const;
  void storeInDb(Db *db, int iech0, int iptr, double zestim, double zstdev);
  void interpolateSelectivity(const Selectivity* selecin);

  void setFlagTonnageCorrect(bool flagTonnageCorrect) { _flagTonnageCorrect = flagTonnageCorrect; }
  void setZmax(double zmax) { _zmax = zmax; }
  bool isFlagTonnageCorrect() const { return _flagTonnageCorrect; }
  double getZmax() const { return _zmax; }
  bool isOnlyZDefined() const { return _flagOnlyZDefined; }

  const Table getStats() const { return _stats; }

private:
  VectorString _getAllNames() const;
  void _printQTvars(const char *title, int type, int number) const;
  void _defineVariableRanks();
  bool _isRecoveryDefined() const;
  bool _isValidCut(int icut) const;
  void _interpolateInterval(double zval,
                            double zi0,
                            double zi1,
                            double ti0,
                            double ti1,
                            double qi0,
                            double qi1,
                            double *tval,
                            double *qval,
                            double tol = EPSILON3);
  void _concatenate(VectorString& names,
                    const ESelectivity& code,
                    int mode) const;
  bool _isMultiplied(const ESelectivity& code) const;
  void _defineAutomaticCutoffs(const VectorDouble& tab, double eps = EPSILON3);

private:
  VectorDouble _Zcut;
  Table _stats;
  double _zmax;
  double _proba;
  bool   _flagTonnageCorrect;
  MatrixInt _numberQT;
  MatrixInt _rankQT;
  bool _flagOnlyZDefined;
};
