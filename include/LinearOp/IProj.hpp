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

#include <Basic/VectorNumT.hpp>

namespace gstlrn
{
class GSTLEARN_EXPORT IProj
{
public:
  IProj() { }
  virtual ~IProj() { }
  VectorDouble point2mesh(const VectorDouble& inv) const;
  VectorDouble mesh2point(const VectorDouble& inv) const;
  Id point2mesh(const VectorDouble& inv, VectorDouble& outv) const;
  Id mesh2point(const VectorDouble& inv, VectorDouble& outv) const;
#ifndef SWIG
  Id point2mesh(const constvect inv, vect out) const;
  Id mesh2point(const constvect inv, vect out) const;
#endif

  virtual Id getNApex() const = 0;
  virtual Id getNPoint() const = 0;


#ifndef SWIG
  Id addMesh2point(const constvect inv, vect outv) const;
  Id addPoint2mesh(const constvect inv, vect outv) const;
  void mesh2point2mesh(const constvect inv, vect outv) const;
  void point2mesh2point(const constvect inv, vect outv) const;

protected:
  virtual Id _addPoint2mesh(const constvect inv, vect outv) const
  {
    DECLARE_UNUSED(inv);
    DECLARE_UNUSED(outv);
    return 1;
  }
  virtual Id _addMesh2point(const constvect inv, vect outv) const
  {
    DECLARE_UNUSED(inv);
    DECLARE_UNUSED(outv);
    return 1;
  }
  #endif

private:
  mutable VectorDouble _workpoint;
  mutable VectorDouble _workmesh;

};
}