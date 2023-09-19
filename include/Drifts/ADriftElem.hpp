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

#include "Enum/EDrift.hpp"

#include "Drifts/ADrift.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CovContext.hpp"

/* Elementary Drift function
 * */

class Db;

class GSTLEARN_EXPORT ADriftElem : public ADrift, public ICloneable
{
public:
  ADriftElem(const EDrift &type,
             const CovContext &ctxt = CovContext());
  ADriftElem(const ADriftElem &r);
  ADriftElem& operator= (const ADriftElem &r);
  virtual ~ADriftElem();

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ASpaceObject Interface
  virtual bool isConsistent(const ASpace* space) const override;

  /// ADrift Interface
  virtual int getNVariables() const override { return _ctxt.getNVar(); }

  /// Interface for daughter classes
  virtual VectorInt getPowers() const { return VectorInt(); }
  virtual int       getRankFex() const { return 0; }

  // ADriftelem Interface
  virtual String getDriftName() const = 0;
  virtual int    getOrderIRF() const = 0;
  virtual int    getOrderIRFIdim(int idim) const = 0;
  virtual int    getNDim() const { return 0; }
  virtual bool   isDriftExternal() const { return false; }
  virtual double eval(const Db* db,int iech) const override = 0;

  const EDrift&  getDriftType() const { return _type; }

  void copyCovContext(const CovContext& ctxt) { _ctxt.copyCovContext(ctxt); }

private:
  CovContext  _ctxt;  /* Context (space, number of variables, ...) */
  EDrift _type;       /* Drift function type */
};
