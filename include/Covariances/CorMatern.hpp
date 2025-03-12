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

#include "Basic/Tensor.hpp"
#include "Basic/VectorNumT.hpp"
#include "Covariances/CorAniso.hpp"
#include "Matrix/MatrixSquareSymmetric.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"
#include "Basic/ICloneable.hpp"
#include "Space/SpacePoint.hpp"
#include <vector>

class ACov;
class CorAniso;
/**
 * \brief
 * This class describes the Gneiting correlation function.
 *
 */
class GSTLEARN_EXPORT CorMatern: public ACov
{
public:
  CorMatern(const VectorDouble &ranges = VectorDouble(),
            const VectorDouble &angle = VectorDouble(),
            const VectorDouble& coeffScales = VectorDouble(), 
            const VectorDouble& params = VectorDouble(),
            bool flagRange = true);
  CorMatern(const CorMatern& r);
  CorMatern& operator=(const CorMatern& r);
  virtual ~CorMatern();
  IMPLEMENT_CLONING(CorMatern)

  bool isConsistent(const ASpace* space) const override
  {
    DECLARE_UNUSED(space)
    return true;
  }
  /// ACov Interface
 

  virtual int getNVar() const override { return _nVar; }
  double getCorMax(int ivar, int jvar) const { return _corMax.getValue(ivar, jvar); }
  double computeScale(int ivar, int jvar) const; 
  double computeParam(int ivar, int jvar) const;   
  protected:
  virtual double _eval(const SpacePoint& p1,
                     const SpacePoint& p2,
                     int ivar                = 0,
                     int jvar                = 0,
                     const CovCalcMode* mode = nullptr) const override;
  void _optimizationSetTarget(SpacePoint& pt) const override;

 
private:
  void _optimizationPreProcess(int mode, const std::vector<SpacePoint>& ps) const override;
  void _optimizationPostProcess() const override;
  
private:
    int _nVar;
    const CorAniso* _corRef;
    mutable CorAniso _corMatern;
    VectorDouble _coeffScales; //scale factor for each variable 
                               //starting from the second one 
                               //(first one is guided by the tensor of 
                               //_corMatern)
    VectorDouble _params; //parameters of the Matern correlation function
    
    MatrixSquareSymmetric _corMax;
    VectorDouble _angles; 
  

};

