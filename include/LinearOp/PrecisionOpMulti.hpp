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

#include "Model/Model.hpp"
#include "LinearOp/ALinearOpMulti.hpp"
#include "LinearOp/PrecisionOp.hpp"
#include "Basic/VectorNumT.hpp"

class Model;

/**
 * Class to store objects for SPDE
 */
class GSTLEARN_EXPORT PrecisionOpMulti : public ALinearOpMulti {

public:
  PrecisionOpMulti(Model* model = nullptr, 
                   const std::vector<AMesh*>& meshes = std::vector<AMesh*>());
  PrecisionOpMulti(const PrecisionOpMulti &m)= delete;
  PrecisionOpMulti& operator= (const PrecisionOpMulti &m)= delete;
  virtual ~PrecisionOpMulti();

  int setModel(Model* model);
  int setMeshes(const std::vector<AMesh*>& meshes);
  void evalDirect(const VectorDouble &vecin, VectorDouble &vecout);

private:
  bool _isValidModel(Model* model);
  bool _isValidMeshes(const std::vector<AMesh*>& meshes);
  bool _matchModelAndMeshes();
  int  _getNVar() const;
  int  _getNCov() const;
  int  _getNMesh() const;
  int  _getSize() const;
  int  _buildInvSills();

private:
  bool _isValid;
  VectorInt _covList;
  VectorInt _nmeshList;
  std::vector<MatrixSquareSymmetric> _invSills;
  std::vector<PrecisionOp> _pops;

  Model* _model; // Not to be deleted
  std::vector<AMesh*> _meshes; // Not to be deleted
};
