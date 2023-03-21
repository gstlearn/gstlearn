/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Authors: <authors>                                                         */
/* Website: <website>                                                         */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "Enum/EDrift.hpp"

#include "Drifts/ADrift.hpp"
#include "Drifts/ADriftElem.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/ASerializable.hpp"
#include "Basic/ICloneable.hpp"
#include "Covariances/CovContext.hpp"

/* Elementary Drift function
 * */

class Db;

class GSTLEARN_EXPORT ADriftElem : public ADrift, public ASerializable, public ICloneable
{
public:
  ADriftElem(const EDrift &type,
             const CovContext &ctxt = CovContext(),
             int rankFex = 0);
  ADriftElem(const ADriftElem &r);
  ADriftElem& operator= (const ADriftElem &r);
  virtual ~ADriftElem();

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ASpaceObject Interface
  virtual bool isConsistent(const ASpace* space) const override;

  /// ADrift Interface
  virtual int getNVariables() const override { return _ctxt.getNVar(); }

  virtual String getDriftSymbol() const = 0;
  virtual String getDriftName() const = 0;
  virtual int    getOrderIRF() const = 0;
  virtual int    getNDim() const { return 0; }
  virtual bool   getDriftExternal() const { return false; }
  virtual double eval(const Db* db,int iech) const override = 0;

  int  getRankFex() const { return _rankFex; }
  void setRankFex(int rankFex) { _rankFex = rankFex; }
  const EDrift& getType() const { return _type; }
  void setType(const EDrift& type) { _type = type; }

  void copyCovContext(const CovContext& ctxt) { _ctxt.copyCovContext(ctxt); }

protected:
  /// Interface to ASerializable
  virtual bool _deserialize(std::istream& is, bool verbose = false) override;
  virtual bool _serialize(std::ostream& os, bool verbose = false) const override;
  String _getNFName() const override { return "Drift"; }

private:
  CovContext  _ctxt;  /* Context (space, number of variables, ...) */
  EDrift _type;       /* Drift function type */
  int _rankFex;       /* Rank of the external drift */
};
