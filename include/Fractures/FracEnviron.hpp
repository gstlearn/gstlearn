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
#include "Basic/ASerializable.hpp"
#include "Fractures/FracFault.hpp"
#include "Fractures/FracFamily.hpp"

class GSTLEARN_EXPORT FracEnviron: public AStringable, public ASerializable
{
public:
  FracEnviron(double xmax = 0.,
              double ymax = 0.,
              double deltax = 0.,
              double deltay = 0.,
              double mean = 0.,
              double stdev = 0.);
  FracEnviron(const FracEnviron& r);
  FracEnviron& operator=(const FracEnviron& r);
  virtual ~FracEnviron();

  /// Interface to AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  static FracEnviron* createFromNF(const String& neutralFilename,
                                   bool verbose = true);
  static FracEnviron* create(double xmax = 0.,
                             double ymax = 0.,
                             double deltax = 0.,
                             double deltay = 0,
                             double mean = 0.,
                             double stdev = 0.);

  int getNFamilies() const { return (int) _families.size(); }
  int getNFaults() const { return (int) _faults.size(); }

  double getDeltax() const { return _deltax; }
  double getDeltay() const { return _deltay; }
  double getMean() const { return _mean; }
  double getStdev() const { return _stdev; }
  double getXmax() const { return _xmax; }
  double getYmax() const { return _ymax; }
  double getXextend() const;

  const FracFault& getFault(int i) const { return _faults[i]; }
  const FracFamily& getFamily(int i) const { return _families[i]; }

  void addFamily(const FracFamily& family) { _families.push_back(family); }
  void addFault(const FracFault& fault) { _faults.push_back(fault); }

protected:
  /// Interface for ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os,bool verbose = false) const override;
  String _getNFName() const override { return "Fracture Environ"; }

private:
  double _xmax;                 //!< Maximum horizontal distance
  double _ymax;                 //!< Maximum vertical distance
  double _deltax;               //!< Dilation along the horizontal axis
  double _deltay;               //!< Dilation along the vertical axis
  double _mean;                 //!< Mean of thickness distribution
  double _stdev;                //!< Standard deviation of thickness distribution
  std::vector<FracFamily> _families; //!< Family definition
  std::vector<FracFault>  _faults;   //!< Fault definition
};
