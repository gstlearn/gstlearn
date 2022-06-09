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

#include "gstlearn_export.hpp"
#include "geoslib_define.h"

#include "Neigh/ENeigh.hpp"
#include "Neigh/NeighMoving.hpp"

class Db;

class GSTLEARN_EXPORT NeighMovingWithScreens: public NeighMoving
{
public:
  using VectorScreens = std::vector<std::pair<VectorDouble, VectorDouble>>;

  NeighMovingWithScreens(
    int ndim = 2,
    const VectorDouble &segments = VectorDouble(),
    bool flag_xvalid = false
  );
  NeighMovingWithScreens(const NeighMovingWithScreens& r) = default;
  NeighMovingWithScreens& operator=(const NeighMovingWithScreens& r) = default;
  virtual ~NeighMovingWithScreens() = default;

  /// Interface for ANeighParam
  virtual ENeigh getType() const override { return ENeigh::MOVING_WITH_SCREENS; }

  int reset(
    int ndim,
    bool flag_xvalid,
    const VectorDouble &screens,
    int nmaxi,
    double radius,
    int nmini = 1,
    int nsect = 1,
    int nsmax = ITEST,
    VectorDouble coeffs = VectorDouble(),
    VectorDouble angles = VectorDouble(),
    double distcont = TEST
  );

  static NeighMovingWithScreens* create(
    int ndim,
    bool flag_xvalid,
    const VectorDouble &screens,
    int nmaxi,
    double radius = TEST,
    int nmini = 1,
    int nsect = 1,
    int nsmax = ITEST,
    VectorDouble coeffs = VectorDouble(),
    VectorDouble angles = VectorDouble(),
    double distcont = TEST
  );

  // TODO Override to add _screens
  // String toString(const AStringFormat* strfmt = nullptr) const override;
  // int dumpToNF(const String& neutralFilename, bool verbose = false) const;
  // static NeighMovingWithScreens* createFromNF(const String& neutralFilename, bool verbose = false);
  
  inline const VectorScreens& getScreens() const { return _screens; }
  void setScreens(const VectorDouble &screens);

protected:
  // TODO Override to add _screens
  // int _deserialize(std::istream& is, bool verbose = false) override;
  // int _serialize(std::ostream& os, bool verbose = false) const override;

private:
  VectorScreens _screens;  /* Screens separating points */
};
