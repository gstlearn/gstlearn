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

#include "Enum/EPowerPT.hpp"
#include "Basic/ICloneable.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/VectorT.hpp"
#include "Mesh/AMesh.hpp"
#include <memory>

#ifndef SWIG

#include <Eigen/Core>
#include <Eigen/Dense>
#endif

#ifndef SWIG
#  include "LinearOp/ALinearOpEigenCG.hpp"
#else
#  include "LinearOp/ALinearOp.hpp"
#endif
namespace gstlrn {
class CovAniso;
class EConsElem;
class ICloneable;

}
/**
 * \brief Shift Operator for performing the basic tasks of SPDE
 */
#ifndef SWIG
DECLARE_EIGEN_TRAITS(AShiftOp)
#endif

namespace gstlrn {

class GSTLEARN_EXPORT AShiftOp: public ICloneable,
#ifndef SWIG
  public ALinearOpEigenCG<AShiftOp>
#else
  public ALinearOp
#endif
{
public:
  AShiftOp(CovAniso* cova = nullptr, Id napices = 0);
  AShiftOp(const AShiftOp& shift);
  AShiftOp& operator=(const AShiftOp& shift);
  virtual void prodLambda(const VectorDouble& x,
                          VectorDouble& y,
                          const EPowerPT& power) const;
  virtual ~AShiftOp();
  virtual double getMaxEigenValue() const;

  virtual void normalizeLambdaBySills(const AMesh*) = 0;
  const VectorDouble& getLambdas() const { return _Lambda; }
  virtual double getLambda(Id iapex) const { return _Lambda[iapex]; }
  virtual double logDetLambda() const;
  static std::shared_ptr<CovAniso> cloneAndCast(const CovAniso* cova);
  static std::shared_ptr<CovAniso> cloneAndCast(const std::shared_ptr<CovAniso> &cova);
  Id getSize() const override { return _napices; }

#ifndef SWIG
    virtual void addProdLambda(const constvect x, vect y, const EPowerPT& power) const;
    void prodLambda(const constvect x, vect y, const EPowerPT& power) const;
    void prodLambda(const VectorDouble& x, vect y, const EPowerPT& power) const;
    void prodLambda(const constvect x, VectorDouble& y, const EPowerPT& power) const;
#endif
#ifndef SWIG
    Id _addToDest(const constvect inv, vect outv) const override = 0;
#endif

private:
    virtual double _getMaxEigenValue() const = 0;
protected:
    std::shared_ptr<CovAniso>& _getCovAniso();
    void _setCovAniso(const CovAniso* cova);
    bool _isNoStat();
    bool _isGlobalHH();

protected:
    VectorDouble _Lambda;
    Id _napices;
    // Following list of members are there to ease the manipulation and reduce
    // argument list
    std::shared_ptr<CovAniso> _cova;
};

} // namespace gstlrn

#ifndef SWIG
  DECLARE_EIGEN_PRODUCT(AShiftOp)
#endif
