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

// WARNING: Make this include list as small as possible!
#include "Basic/AStringable.hpp"

class GSTLEARN_EXPORT Option_AutoFit : public AStringable
{
private:
  bool _verbose;                   /* Verbose option */
  int _wmode;                      /* Weighting option (used in Goulard) */
  int _maxiter;                    /* Maximum number of iterations */
  int _flag_intrinsic;             /* Ask for an intrinsic model */
  double _tolstop;                 /* Tolerance for the stopping criterion */
  double _tolred;                  /* Scaled tolerance (used in calculations) */
  double _epsdelta;                /* Tolerance for the search */
  double _tolsigma;                /* Percentage of variance below which a structure is discarded */
  double _initdelta;               /* Initial radius of the trusting area */
  double _constantSillValue;       /* Constant Sill as a constraint */
  VectorDouble _constantSills;     /* Array of constant Sills (expanded to the variables) */

 public:
  Option_AutoFit();
  Option_AutoFit(const Option_AutoFit &m);
  Option_AutoFit& operator= (const Option_AutoFit &m);
  virtual ~Option_AutoFit();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  double getConstantSillValue() const { return _constantSillValue; }
  double getEpsdelta() const { return _epsdelta; }
  int getFlagIntrinsic() const { return _flag_intrinsic; }
  double getInitdelta() const { return _initdelta; }
  int getMaxiter() const { return _maxiter; }
  double getTolred() const { return _tolred; }
  double getTolsigma() const { return _tolsigma; }
  double getTolstop() const { return _tolstop; }
  bool getVerbose() const { return _verbose; }
  int getWmode() const { return _wmode; }
  const VectorDouble& getConstantSills() const { return _constantSills; }
  double getConstantSills(int ivar) const { return _constantSills[ivar]; }

  void setConstantSillValue(double value) { _constantSillValue = value; }
  void setEpsdelta(double epsdelta) { _epsdelta = epsdelta; }
  void setFlagIntrinsic(int flagIntrinsic) { _flag_intrinsic = flagIntrinsic; }
  void setInitdelta(double initdelta) { _initdelta = initdelta; }
  void setMaxiter(int maxiter) { _maxiter = maxiter; }
  void setTolred(double tolred) { _tolred = tolred; }
  void setTolsigma(double tolsigma) { _tolsigma = tolsigma; }
  void setTolstop(double tolstop) { _tolstop = tolstop; }
  void setVerbose(bool verbose) { _verbose = verbose; }
  void setWmode(int wmode) { _wmode = wmode; }
  void setConstantSills(int nvar);
};
