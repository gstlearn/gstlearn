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
#include "Basic/ICloneable.hpp"
#include "Basic/VectorNumT.hpp"

/* Elementary Drift function
 * */

class Db;

/**
 * \brief
 * This class describes one basic Drift Function.
 *
 * It is the uppermost class of the Drift Tree and is conceived as simple as possible on purpose
 * (in order to let the user defined its own version if necessary): it must simply be able to return its value
 * at the location of one sample from a Db.
 *
 * This returned value depends on the implementation of this basic drift function and mainly depends upon:
 * - the space dimension
 * - the coordinates of the target
 *
 * If NDIM represents the space dimension, each basic drift function belongs to one of the following categories:
 *
 * - an internal drift function (**DRIFTM**) characterized by a vector of coefficients P (of dimension NDIM):
 * it returns the numerical expression elaborated starting from the coordinates (X) of one sample such as:
 *
 *      X_1**P_1 * X_2**P_2 * ... * X_NDIM**P_NDIM
 *
 * - an **external drift** function (**DRIFTF**) identified by its rank (corresponding locator ELoc::F)
 * for one sample.
 */
class GSTLEARN_EXPORT ADrift : public AStringable, public ICloneable
{
public:
  ADrift();
  ADrift(const ADrift &r);
  ADrift& operator= (const ADrift &r);
  virtual ~ADrift();

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// Interface for daughter classes
  virtual VectorInt getPowers() const { return VectorInt(); }
  virtual int       getRankFex() const { return 0; }

  // ADriftelem Interface
  virtual String getDriftName() const = 0;
  virtual int    getOrderIRF() const = 0;
  virtual int    getOrderIRFIdim(int idim) const = 0;
  virtual double eval(const Db* db,int iech) const = 0;
  virtual int    getDriftNDimMax() const { return 0; }
  virtual bool   isDriftExternal() const { return false; }

};
