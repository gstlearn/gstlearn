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

#include "Basic/AStringable.hpp"

/**
 * This class defines the options and parameters used during the Variogram Fitting.
 * All the parameters described hereafter are either available in the construction,
 * or can be set using a specific get() function.
 * - verbose: Ask for a verbose option during the Automatic Model fitting
 * - wmode: Weighting option (see comments on the setWMode() function)
 * - maxiter: Maximum number of iterations
 * - flag_intrinsic: When True, fit a Model which includes at least one Intrinsic basic structure
 * - tolstop: Define an absolute criterion used for stopping the iterations
 * - tolred: Define the relative criterion used for stopping the iterations
 * - epsdelta: Define the tolerance used for the search
 * - tolsigma: Percentage of the variance below which a basic structure is discarded
 * - initdelta: Initial radius of the trusting area
 */
class GSTLEARN_EXPORT Option_AutoFit : public AStringable
{
 public:
  Option_AutoFit();
  Option_AutoFit(const Option_AutoFit &m);
  Option_AutoFit& operator= (const Option_AutoFit &m);
  virtual ~Option_AutoFit();

  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  double getEpsdelta() const { return _epsdelta; }
  int getFlagIntrinsic() const { return _flag_intrinsic; }
  double getInitdelta() const { return _initdelta; }
  int getMaxiter() const { return _maxiter; }
  double getTolred() const { return _tolred; }
  double getTolsigma() const { return _tolsigma; }
  double getTolstop() const { return _tolstop; }
  bool getVerbose() const { return _verbose; }
  int getWmode() const { return _wmode; }
  bool isUseEigenLibrary() const { return _useEigenLibrary; }

  void setEpsdelta(double epsdelta) { _epsdelta = epsdelta; }
  void setFlagIntrinsic(int flagIntrinsic) { _flag_intrinsic = flagIntrinsic; }
  void setInitdelta(double initdelta) { _initdelta = initdelta; }
  void setMaxiter(int maxiter) { _maxiter = maxiter; }
  void setTolred(double tolred) { _tolred = tolred; }
  void setTolsigma(double tolsigma) { _tolsigma = tolsigma; }
  void setTolstop(double tolstop) { _tolstop = tolstop; }
  void setVerbose(bool verbose) { _verbose = verbose; }
  /**
 * Set the type of the weighting function used in the fitting procedure.
 * This function is defined in the case of several directional experimental variograms,
 * calculated in a multivariate case:
 * 0: The weight is constant  
 * 1: The weight is proportional to the number of pairs  
 * 2: The weight is proportional to the number of pairs and inverse proportional to the distance  
 * 3: The weight is inverse proportional to the number of lags for each direction  
 * @param wmode       type of weighting function (0, 1, 2 or 3, see above) 
 * @note The default value for wmode is 2 
 */
  void setWmode(int wmode) { _wmode = wmode; }
  void setUseEigenLibrary(bool useEigenLibrary) { _useEigenLibrary = useEigenLibrary; }

 private:
   bool   _verbose;                 /* Verbose option */
   int    _wmode;                   /* Weighting option (used in Goulard) */
   int    _maxiter;                 /* Maximum number of iterations */
   int    _flag_intrinsic;          /* Ask for an intrinsic model */
   double _tolstop;                 /* Tolerance for the stopping criterion */
   double _tolred;                  /* Scaled tolerance (used in calculations) */
   double _epsdelta;                /* Tolerance for the search */
   double _tolsigma;                /* Percentage of variance below which a structure is discarded */
   double _initdelta;               /* Initial radius of the trusting area */
   bool   _useEigenLibrary;         /* Force/discard use of Eigen library for computing eigen factors */
};
