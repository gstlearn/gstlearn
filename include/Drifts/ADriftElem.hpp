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

#include "Drifts/ADrift.hpp"
#include "Drifts/EDrift.hpp"
#include "Basic/Vector.hpp"
#include "Basic/AStringable.hpp"
#include "Basic/IClonable.hpp"
#include "Covariances/CovContext.hpp"


/* Elementary Drift function
 * */

class Db;

class ADriftElem : public ADrift, public IClonable
{
public:
  ADriftElem(const EDrift& type, const CovContext& ctxt, int rankFex = 0);
  ADriftElem(const ADriftElem &r);
  ADriftElem& operator= (const ADriftElem &r);
  virtual ~ADriftElem();

  virtual IClonable* clone() const override = 0;

  ///////////////////////////////////////////////////
  /// ASpaceObject AStringable
  virtual std::string toString(int level = 0) const override;
  ///////////////////////////////////////////////////

  ///////////////////////////////////////////////////
  /// ASpaceObject Interface
  virtual bool isConsistent(const ASpace* space) const override;
  ///////////////////////////////////////////////////

  ///////////////////////////////////////////////////
  /// ADrift Interface
  virtual int getNVariables() const override
  {
    return _ctxt.getNVar();
   }

  virtual String getDriftSymbol() const = 0;
  virtual String getDriftName() const = 0;
  virtual int    getOrderIRF() const = 0;
  virtual double eval(const Db* db,int iech) const override = 0;

  void setContext(const CovContext& ctxt);

  int getOrderIrf() const { return _orderIRF; }
  void setOrderIrf(int orderIrf) { _orderIRF = orderIrf; }
  int getRankFex() const { return _rankFex; }
  void setRankFex(int rankFex) { _rankFex = rankFex; }
  const EDrift& getType() const { return _type; }
  void setType(const EDrift& type) { _type = type; }

private:
  CovContext  _ctxt;  /* Context (space, irfDegree, field, ...) */
  EDrift _type;       /* Drift function type */
  int _rankFex;       /* Rank of the external drift */
  int _orderIRF;      /* Rank of the IRF induced by presence of Drift Elementary */
};
