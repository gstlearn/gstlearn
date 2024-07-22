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
#include <functional>

#include "Model/Model.hpp"
#include "LinearOp/ALinearOpMulti.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "Basic/VectorNumT.hpp"
#include "Basic/AStringable.hpp"

class Model;

/**
 * Class to store objects for SPDE
 */
class GSTLEARN_EXPORT PrecisionOpMulti : public ALinearOpMulti, public AStringable
{
public:
  PrecisionOpMulti(Model* model = nullptr, 
                   const std::vector<AMesh*>& meshes = std::vector<AMesh*>());
  PrecisionOpMulti(const PrecisionOpMulti &m)= delete;
  PrecisionOpMulti& operator= (const PrecisionOpMulti &m)= delete;
  virtual ~PrecisionOpMulti();

  /// AStringable Interface
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  /// ALinearOpMulti Interface
  virtual int sizes() const override;
  virtual int size(int imesh) const override;
  virtual void _evalDirect(const VectorVectorDouble &inv,
                           VectorVectorDouble &outv) const override { ; }

  int  setModel(Model* model);
  int  setMeshes(const std::vector<AMesh*>& meshes);
  void clearMeshes();
  void addMesh(AMesh* mesh);

  int evalSimulateInPlace(const VectorDouble& vecin,
                                VectorDouble& vecout);
  
  int evalDirectInPlace(const VectorDouble& vecin,
                              VectorDouble& vecout);
  VectorDouble evalDirect(const VectorDouble& vecin);
  VectorDouble evalSimulate(const VectorDouble& vecin);

  private: 
  int _evalOperator(const VectorDouble& vecin,
                          VectorDouble& vecout,
                    const std::function<int()> &func,
                    const std::function<int(int)> & jvarlimit,
                    const std::function<double(int,int,int)> & getTerm,
                    const std::function<void(int,VectorDouble&,VectorDouble&)> &opFunc) const;
  bool _isValidModel(Model* model);
  bool _isValidMeshes(const std::vector<AMesh*>& meshes);
  bool _matchModelAndMeshes();
  int  _getNVar() const;
  int  _getNCov() const;
  int  _getNMesh() const;
  int  _buildInvSills();
  int  _buildCholSills();
  void _popsClear();

private:
  bool _isValid;
  VectorInt _covList;
  VectorInt _nmeshList;
  std::vector<PrecisionOp*> _pops;

  std::vector<MatrixSquareSymmetric> _invSills; // Inverse of the Sills
  std::vector<MatrixSquareSymmetric> _cholSills; // Cholesky of the Sills

  Model* _model; // Not to be deleted
  std::vector<AMesh*> _meshes; // Not to be deleted
};
