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

#include "Basic/AStringable.hpp"
#include <limits>

class GSTLEARN_EXPORT ParamInfo : public AStringable
{
  public:
  ParamInfo(double valuser, double valmin = -std::numeric_limits<double>::infinity(),
                            double valmax = std::numeric_limits<double>::infinity(), 
                            bool mute = true);
  ParamInfo(const ParamInfo &m);
  ParamInfo& operator= (const ParamInfo &m);
  virtual ~ParamInfo();
  void setValMin(double valmin);
  void setValMax(double valmax);
  void setMutable(bool mute){_mutable = mute;}
  String toString(const AStringFormat* strfmt) const override;

private :
  double _valUser;
  double _valMin;
  double _valMax;
  double _valMinUser;
  double _valMaxUser;
  bool   _mutable;
};
