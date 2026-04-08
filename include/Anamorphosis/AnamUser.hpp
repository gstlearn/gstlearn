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

#include "geoslib_define.h"
#include "gstlearn_export.hpp"

#include "Anamorphosis/AnamContinuous.hpp"

namespace gstlrn
{

  class GSTLEARN_EXPORT AnamUser: public AnamContinuous
  {
  private:
    double (*_y2z_function)(double);
    double (*_z2y_function)(double);

  public:
    AnamUser();
    AnamUser(const AnamUser& m);
    AnamUser& operator=(const AnamUser& m);
    virtual ~AnamUser();

    /// ICloneable Interface
    IMPLEMENT_CLONING(AnamUser)

    /// AStringable Interface
    String toString(const AStringFormat* strfmt = nullptr) const override;

    /// ASerializable Interface
    String getNFName() const override { return "AnamUser"; }
#ifdef HDF5
    bool deserializeH5(H5::Group& grp) override;
    bool serializeH5(H5::Group& grp) const override;
#endif

    /// AAnam Interface
    const EAnam& getType() const override { return EAnam::fromKey("EXTERNAL"); }

    bool isChangeSupportDefined() const override { return false; }

    /// AnamContinuous Interface
    void calculateMeanAndVariance() override;
    double transformToRawValue(double h) const override;
    double rawToTransformValue(double h) const override;

    void setY2zFunction(double (*y2z_function)(double))
    {
      _y2z_function = y2z_function;
    }

    void setZ2yFunction(double (*z2y_function)(double))
    {
      _z2y_function = z2y_function;
    }

  protected:
    bool _deserializeAscii(std::istream& is) override;
    bool _serializeAscii(std::ostream& os) const override;
  };
} // namespace gstlrn
