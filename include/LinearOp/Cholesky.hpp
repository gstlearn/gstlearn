/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* Created on: 9 avr. 2019 by N. Desassis                                     */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"

#include "LinearOp/ALinearOp.hpp"
#include "Basic/VectorNumT.hpp"

class GSTLEARN_EXPORT Cholesky: public ALinearOp
{
public:
  Cholesky(const cs* mat = nullptr);
  Cholesky(const Cholesky &m);
  Cholesky& operator=(const Cholesky &m);
  virtual ~Cholesky();

  int reset(const cs* mat = nullptr);
  int getSize() const override;
  void evalInverse(const VectorDouble& inv, VectorDouble& outv) const override;
  void simulate(VectorDouble& inv, VectorDouble& outv);
  void stdev(VectorDouble& vcur, bool flagStDev = false);
  void printout(const char *title, bool verbose = false) const;

protected:
  void _evalDirect(const VectorDouble& inv, VectorDouble& outv) const override;

private:
  bool _isDefined() const { return _mat != nullptr; }
  bool _isDecomposed() const { return _matS != nullptr && _matN != nullptr; }
  void _clean();
  void _decompose(bool verbose = false);

private:
  cs* _mat;
  css *_matS;
  csn *_matN;
  mutable VectorDouble _work;
};
