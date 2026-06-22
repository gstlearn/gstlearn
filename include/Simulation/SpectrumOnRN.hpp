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

#include "Basic/VectorNumT.hpp"
#include "Basic/VectorT.hpp"
#include "Matrix/MatrixDense.hpp"
#include "Simulation/BooleanObject.hpp"
#include "geoslib_define.h"
#include "gstlearn_export.hpp"
#include <vector>

namespace gstlrn
{

  class Db;

  /////////////////////////////////////////////////////////////////////////////////
  // abstract class defining the interface of the spectrum on RN
  /////////////////////////////////////////////////////////////////////////////////

  class GSTLEARN_EXPORT SpectrumOnRN
  {
  public:
    SpectrumOnRN(Id nv = 1, Id nd = 2, Id ns = 1000);
    SpectrumOnRN(const SpectrumOnRN& r);
    SpectrumOnRN& operator=(const SpectrumOnRN& r);
    virtual ~SpectrumOnRN();

    Id getNVar() const { return _nv; };

    Id getNDim() const { return _nd; };

    virtual Id getNs() const { return _ns; };

    virtual bool setGamma(const MatrixDense& gamma);

    // virtual functions
#ifndef SWIG
    virtual std::unique_ptr<SpectrumOnRN> clone() const = 0;
#endif
    virtual void
      _computeValues(const VectorDouble& coor, VectorDouble& values) = 0;
    virtual MatrixDense getOmega(Id ifac = 0, Id is = 0) const = 0;
    virtual VectorDouble getPhi(Id ifac = 0, Id is = 0) const = 0;
    virtual MatrixDense getProjection(Id ifac = 0, Id is = 0) const = 0;
    // to simulate Gneiting model
    virtual MatrixDense getOmega0(Id ifac = 0, Id is = 0) const = 0;
    virtual VectorDouble getXi0(Id ifac = 0, Id is = 0) const = 0;

    virtual MatrixDense getGamma(Id is = 0) const
    {
      DECLARE_UNUSED(is)
      return _gamma;
    };

    // interface defined for SpectrumOnRNFactorized
    virtual bool isFactorized() const = 0;
    virtual Id getNFac() const = 0;
    virtual bool addFactor(
      const MatrixDense& omega,
      const VectorDouble& phi,
      const MatrixDense& proj = MatrixDense(),
      const MatrixDense& omega0 = MatrixDense(),
      const VectorDouble& xi0 = VectorDouble()) = 0;

    // interface defined for SpectrumOnRNList
    virtual bool isList() const = 0;
    virtual Id getNSpectrum() const = 0;
#ifndef SWIG
    virtual bool addSpectrum(std::unique_ptr<SpectrumOnRN> sp)
    {
      DECLARE_UNUSED(sp)
      return false;
    };
#endif

    void compute(
      Db* dbout,
      const VectorBool& activeArray,
      VectorVectorDouble& tab);
    MatrixDense computeToMatrix(Db* dbout);

    bool _isValidNs(Id i) const { return (i >= 0) && (i < getNs()); };

    bool _isValidNv(Id i) const { return (i >= 0) && (i < getNVar()); };

    bool _isValidNd(Id i) const { return (i >= 0) && (i < getNDim()); };

    bool _isValidNf(Id i) const { return (i >= 0) && (i < getNFac()); };

    bool _isValidNSpectrum(Id i) const
    {
      return (i >= 0) && (i < getNSpectrum());
    };

    virtual void _initWorkingVariables() = 0;

  protected:
    // spectrum for R^n
    Id _nv;
    Id _nd;
    Id _ns;
    MatrixDense _gamma; // matrix [ns, nv] of the weigths per variable

    // working variables
    VectorDouble
      _u; // working vector of length _ns used by the method _computeValues
    VectorDouble
      _v; // working vector of length _ns used by the method _computeValues
    double _SQRT2;
  };

  /////////////////////////////////////////////////////////////////////////////////
  // concrete class defining a simple spectrum on RN
  /////////////////////////////////////////////////////////////////////////////////

  class GSTLEARN_EXPORT SpectrumOnRNSimple: public SpectrumOnRN
  {
  public:
    SpectrumOnRNSimple(Id nvar = 1, Id ndim = 2, Id ns = 1000);
    SpectrumOnRNSimple(const SpectrumOnRNSimple& r);
    SpectrumOnRNSimple& operator=(const SpectrumOnRNSimple& r);
    virtual ~SpectrumOnRNSimple();

    // virtual functions
#ifndef SWIG
    std::unique_ptr<SpectrumOnRN> clone() const override
    {
      return std::make_unique<SpectrumOnRNSimple>(*this);
    };
#endif
    void
      _computeValues(const VectorDouble& coor, VectorDouble& values) override;
    MatrixDense getOmega(Id ifac = 0, Id is = 0) const override;
    VectorDouble getPhi(Id ifac = 0, Id is = 0) const override;
    // to simulate Gneiting model
    MatrixDense getOmega0(Id ifac = 0, Id is = 0) const override;
    VectorDouble getXi0(Id ifac = 0, Id is = 0) const override;

    MatrixDense getProjection(Id ifac = 0, Id is = 0) const override;

    // interface defined for SpectrumOnRNFactorized
    bool isFactorized() const override { return false; };

    Id getNFac() const override { return 1; };

    bool addFactor(
      const MatrixDense& omega,
      const VectorDouble& phi,
      const MatrixDense& proj = MatrixDense(),
      const MatrixDense& omega0 = MatrixDense(),
      const VectorDouble& xi0 = VectorDouble()) override;

    // interface defined for SpectrumOnRNList
    bool isList() const override { return false; };

    Id getNSpectrum() const override { return 1; };
#ifndef SWIG
    bool addSpectrum(std::unique_ptr<SpectrumOnRN> sp) override
    {
      DECLARE_UNUSED(sp)
      return false;
    };
#endif
    void _initWorkingVariables() override;

  protected:
    // spectrum for R^n
    MatrixDense _omega;
    VectorDouble _phi;
    MatrixDense _omega0;
    VectorDouble _xi0;
  };

  /////////////////////////////////////////////////////////////////////////////////
  // concrete class defining a factorized spectrum on RN
  /////////////////////////////////////////////////////////////////////////////////

  class GSTLEARN_EXPORT SpectrumOnRNFactorized: public SpectrumOnRN
  {
  public:
    SpectrumOnRNFactorized(Id nvar = 1, Id ndim = 2, Id ns = 1000);
    SpectrumOnRNFactorized(const SpectrumOnRNFactorized& r);
    SpectrumOnRNFactorized& operator=(const SpectrumOnRNFactorized& r);
    virtual ~SpectrumOnRNFactorized();

    // virtual functions
#ifndef SWIG
    std::unique_ptr<SpectrumOnRN> clone() const override
    {
      return std::make_unique<SpectrumOnRNFactorized>(*this);
    };
#endif
    void
      _computeValues(const VectorDouble& coor, VectorDouble& values) override;
    MatrixDense getOmega(Id ifac = 0, Id is = 0) const override;
    VectorDouble getPhi(Id ifac = 0, Id is = 0) const override;
    // to simulate Gneiting model
    MatrixDense getOmega0(Id ifac = 0, Id is = 0) const override;
    VectorDouble getXi0(Id ifac = 0, Id is = 0) const override;

    MatrixDense getProjection(Id ifac = 0, Id is = 0) const override;

    // interface defined for SpectrumOnRNFactorized
    bool isFactorized() const override { return true; };

    Id getNFac() const override;
    bool addFactor(
      const MatrixDense& omega,
      const VectorDouble& phi,
      const MatrixDense& proj = MatrixDense(),
      const MatrixDense& omega0 = MatrixDense(),
      const VectorDouble& xi0 = VectorDouble()) override;

    // interface defined for SpectrumOnRNList
    bool isList() const override { return false; };

    Id getNSpectrum() const override { return 1; };
#ifndef SWIG
    bool addSpectrum(std::unique_ptr<SpectrumOnRN> sp) override
    {
      DECLARE_UNUSED(sp)
      return false;
    };
#endif
    void _initWorkingVariables() override;

  protected:
    // spectrum for R^n
    Id _nf;
    std::vector<MatrixDense> _omegas;
    std::vector<VectorDouble> _phis;
    std::vector<MatrixDense> _projections;
    // working variables
    std::vector<VectorDouble> _xproj;
  };

  /////////////////////////////////////////////////////////////////////////////////
  // concrete class defining a spectrum on RN for a CovList: a list of covaraince
  /////////////////////////////////////////////////////////////////////////////////
  class GSTLEARN_EXPORT SpectrumOnRNList: public SpectrumOnRN
  {
  public:
    SpectrumOnRNList(Id nvar = 1, Id ndim = 2, Id ns = 1000);
    SpectrumOnRNList(const SpectrumOnRNList& r);
    SpectrumOnRNList& operator=(const SpectrumOnRNList& r);
    virtual ~SpectrumOnRNList();

    // virtual functions
#ifndef SWIG
    std::unique_ptr<SpectrumOnRN> clone() const override
    {
      return std::make_unique<SpectrumOnRNList>(*this);
    };
#endif
    void
      _computeValues(const VectorDouble& coor, VectorDouble& values) override;
    MatrixDense getOmega(Id ifac = 0, Id is = 0) const override;
    VectorDouble getPhi(Id ifac = 0, Id is = 0) const override;
    MatrixDense getProjection(Id ifac = 0, Id is = 0) const override;
    // to simulate Gneiting model
    MatrixDense getOmega0(Id ifac = 0, Id is = 0) const override;
    VectorDouble getXi0(Id ifac = 0, Id is = 0) const override;

    MatrixDense getGamma(Id is = 0) const override;
    bool setGamma(const MatrixDense& gamma) override;

    // interface defined for SpectrumOnRNFactorized
    bool isFactorized() const override { return false; };

    Id getNFac() const override { return -1; };

    bool addFactor(
      const MatrixDense& omega,
      const VectorDouble& phi,
      const MatrixDense& proj = MatrixDense(),
      const MatrixDense& omega0 = MatrixDense(),
      const VectorDouble& xi0 = VectorDouble()) override;

    // interface defined for SpectrumOnRNList
    bool isList() const override { return true; };

    Id getNSpectrum() const override { return _sps.size(); };
#ifndef SWIG
    bool addSpectrum(std::unique_ptr<SpectrumOnRN> sp) override;
#endif
    void _initWorkingVariables() override;

  protected:
#ifndef SWIG
    std::vector<std::unique_ptr<SpectrumOnRN>>
      _sps; // vector of elementary spectrums
#endif

    // temporary variables
  };

} // namespace gstlrn
