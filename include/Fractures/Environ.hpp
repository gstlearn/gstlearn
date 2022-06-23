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

#include "FracFault.hpp"

#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/Vector.hpp"
#include "Fractures/Family.hpp"

class GSTLEARN_EXPORT Environ: public AStringable, public ASerializable
{
public:
  Environ(double xmax = 0.,
          double ymax = 0.,
          double deltax = 0.,
          double deltay = 0,
          double xextend = 0.,
          double mean = 0.,
          double stdev = 0.);
  Environ(const Environ& r);
  Environ& operator=(const Environ& r);
  virtual ~Environ();

  static Environ* createFromNF(const String& neutralFilename,
                               bool verbose = false);
  static Environ* create(double xmax = 0.,
                         double ymax = 0.,
                         double deltax = 0.,
                         double deltay = 0,
                         double xextend = 0.,
                         double mean = 0.,
                         double stdev = 0.);
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int getNFamilies() const { return (int) _families.size(); }
  int getNFaults() const { return (int) _faults.size(); }

  double getDeltax() const { return _deltax; }
  double getDeltay() const { return _deltay; }
  double getMean() const { return _mean; }
  double getStdev() const { return _stdev; }
  double getXextend() const { return _xextend; }
  double getXmax() const { return _xmax; }
  double getYmax() const { return _ymax; }

  const FracFault& getFault(int i) const { return _faults[i]; }
  const Family& getFamily(int i) const { return _families[i]; }

  void addFamily(const Family& family) { _families.push_back(family); }
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
  double _xextend;              //!< Field extension along horizontal axis
  double _mean;                 //!< Mean of thickness distribution
  double _stdev;                //!< Standard deviation of thickness distribution
  std::vector<Family> _families; //!< Family definition
  std::vector<FracFault>  _faults;   //!< Fault definition
};
