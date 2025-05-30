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

#include "../Matrix/Table.hpp"
#include "Enum/ESelectivity.hpp"

#include "Basic/VectorNumT.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/NamingConvention.hpp"
#include "Matrix/MatrixInt.hpp"

class Db;
class AAnam;

class GSTLEARN_EXPORT Selectivity: public AStringable, public ICloneable
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
  static Selectivity* createByCuts(const VectorDouble& zcut);
  static Selectivity* createByCodes(const std::vector<ESelectivity>& codes,
                                    const VectorDouble& zcuts = VectorDouble(),
                                    bool flag_est = true,
                                    bool flag_std = true,
                                    double proba = TEST,
                                    bool verbose = false);
  static Selectivity* createByKeys(const VectorString& keys,
                                   const VectorDouble& zcuts = VectorDouble(),
                                   bool flag_est = true,
                                   bool flag_std = true,
                                   double zmax = TEST,
                                   bool flag_tonnage_correct = false,
                                   double proba = TEST,
                                   bool verbose = false);
  static Selectivity* createInterpolation(const VectorDouble& zcuts,
                                          const Selectivity& selecin,
                                          bool verbose);

  int calculateFromDb(const Db* db, bool autoCuts = false);
  int calculateFromArray(const VectorDouble& tab,
                         const VectorDouble& weights = VectorDouble(),
                         bool autoCuts = false);
  int calculateFromAnamorphosis(AAnam* anam);

  Table eval(const Db *db, bool autoCuts = false);
  Table evalFromArray(const VectorDouble &tab,
                            const VectorDouble &weights = VectorDouble(),
                            bool autoCuts = false);
  Table evalFromAnamorphosis(AAnam *anam);

  void   resetCuts(const VectorDouble& zcuts);
  int    getNCuts() const { return static_cast<int>(_Zcut.size()); }
  static int getNQT() { return static_cast<int>(ESelectivity::getSize()); }
  int    getNVar() const;
  String getVariableName(const ESelectivity& code, int icut, int mode) const;
  String getVariableName(int rank0) const;
  VectorString getVariableNames() const;

  void   setZcut(int iclass, double zcut);
  void   setBest(int iclass, double best);
  void   setMest(int iclass, double mest);
  void   setQest(int iclass, double qest);
  void   setQstd(int iclass, double qstd);
  void   setTest(int iclass, double test);
  void   setTstd(int iclass, double tstd);

  double getZcut(int iclass) const;
  double getBest(int iclass) const;
  double getMest(int iclass) const;
  double getQest(int iclass) const;
  double getQstd(int iclass) const;
  double getTest(int iclass) const;
  double getTstd(int iclass) const;
  const VectorDouble& getZcut() const { return _Zcut; }

  void calculateBenefitAndGrade();
  void dumpGini() const;
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
  int  getNQTEst(const ESelectivity& code) const;
  int  getNQTStd(const ESelectivity& code) const;
  VectorInt getNQTEst() const;
  VectorInt geNQTStd() const;
  void storeInDb(Db *db, int iech0, int iptr, double zestim, double zstdev) const;
  void interpolateSelectivity(const Selectivity* selecin);

  void setFlagTonnageCorrect(bool flagTonnageCorrect) { _flagTonnageCorrect = flagTonnageCorrect; }
  void setZmax(double zmax) { _zmax = zmax; }
  bool isFlagTonnageCorrect() const { return _flagTonnageCorrect; }
  double getZmax() const { return _zmax; }
  bool isOnlyZDefined() const { return _flagOnlyZDefined; }

  Table getStats() const;
  Table getAllStats() const { return _stats; }

  const MatrixInt& getCountQT() const { return _numberQT; }
  const MatrixInt& getRankQt() const { return _rankQT; }

private:
  static VectorString _getAllNames();
  static void _printQTvars(const char *title, int type, int number);
  void _defineVariableRanks();
  bool _isRecoveryDefined() const;
  bool _isValidCut(int iclass) const;
  static void _interpolateInterval(double zval,
                                   double zi0,
                                   double zi1,
                                   double ti0,
                                   double ti1,
                                   double qi0,
                                   double qi1,
                                   double* tval,
                                   double* qval,
                                   double tol = EPSILON3);
  void _concatenate(VectorString& names,
                    const ESelectivity& code,
                    int mode) const;
  static bool _isMultiplied(const ESelectivity& code);
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

GSTLEARN_EXPORT int dbSelectivity(Db *db,
                                  const String &name,
                                  const VectorDouble& zcuts,
                                  const NamingConvention &namconv = NamingConvention(
                                      "Selectivity"));
