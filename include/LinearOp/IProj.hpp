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

class GSTLEARN_EXPORT IProj
{
public:
  IProj() { }
  virtual ~IProj() { }
  int point2mesh(const VectorDouble& inv, VectorDouble& outv) const;
  int mesh2point(const VectorDouble& inv, VectorDouble& outv) const;
#ifndef SWIG
  int point2mesh(const constvect inv, vect out) const;
  int mesh2point(const constvect inv, vect out) const;
#endif

  virtual int getApexNumber() const = 0;
  virtual int getPointNumber() const = 0;

#ifndef SWIG
  int addMesh2point(const constvect inv, vect outv) const;
  int addPoint2mesh(const constvect inv, vect outv) const;

protected:
  virtual int _addPoint2mesh(const constvect inv, vect outv) const
  {
    DECLARE_UNUSED(inv);
    DECLARE_UNUSED(outv);
    return 1;
  }
  virtual int _addMesh2point(const constvect inv, vect outv) const
  {
    DECLARE_UNUSED(inv);
    DECLARE_UNUSED(outv);
    return 1;
  }
  #endif
};
