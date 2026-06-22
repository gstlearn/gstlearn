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
#include "Basic/VectorHelper.hpp"
#include "Covariances/CovContext.hpp"
#include "Drifts/DriftFactory.hpp"
#include "Drifts/DriftList.hpp"
#include "Drifts/DriftM.hpp"
#include "utils.hpp"

using namespace gstlrn;

/****************************************************************************/
/*!
** Main Program for testing generic DriftIRF
**
*****************************************************************************/
int main(int argc, char* argv[])
{
  std::stringstream sfn;
  sfn << gslBaseName(__FILE__) << ".out";
  StdoutRedirect sr(sfn.str(), argc, argv);

  // Test case 1: order=0, ndim=1
  message("Test case 1: order=0, ndim=1\n");
  CovContext ctxt1(1, 1);
  DriftList* drifts1 = DriftFactory::createDriftListFromIRF(0, 0, ctxt1);
  message("Number of drifts: %d\n", drifts1->getNDrift());
  for (Id i = 0; i < drifts1->getNDrift(); i++)
  {
    message("  %s\n", drifts1->getDrift(i)->getDriftName().c_str());
  }
  delete drifts1;

  // Test case 2: order=1, ndim=2
  message("\nTest case 2: order=1, ndim=2\n");
  CovContext ctxt2(1, 2);
  DriftList* drifts2 = DriftFactory::createDriftListFromIRF(1, 0, ctxt2);
  message("Number of drifts: %d\n", drifts2->getNDrift());
  for (Id i = 0; i < drifts2->getNDrift(); i++)
  {
    message("  %s\n", drifts2->getDrift(i)->getDriftName().c_str());
  }
  delete drifts2;

  // Test case 3: order=2, ndim=2
  message("\nTest case 3: order=2, ndim=2\n");
  CovContext ctxt3(1, 2);
  DriftList* drifts3 = DriftFactory::createDriftListFromIRF(2, 0, ctxt3);
  message("Number of drifts: %d\n", drifts3->getNDrift());
  for (Id i = 0; i < drifts3->getNDrift(); i++)
  {
    message("  %s\n", drifts3->getDrift(i)->getDriftName().c_str());
  }
  delete drifts3;

  // Test case 4: order=3, ndim=2 (now should work!)
  message("\nTest case 4: order=3, ndim=2\n");
  CovContext ctxt4(1, 2);
  DriftList* drifts4 = DriftFactory::createDriftListFromIRF(3, 0, ctxt4);
  message("Number of drifts: %d\n", drifts4->getNDrift());
  for (Id i = 0; i < drifts4->getNDrift(); i++)
  {
    message("  %s\n", drifts4->getDrift(i)->getDriftName().c_str());
  }
  delete drifts4;

  // Test case 5: order=3, ndim=3
  message("\nTest case 5: order=3, ndim=3\n");
  CovContext ctxt5(1, 3);
  DriftList* drifts5 = DriftFactory::createDriftListFromIRF(3, 0, ctxt5);
  message("Number of drifts: %d\n", drifts5->getNDrift());
  for (Id i = 0; i < drifts5->getNDrift(); i++)
  {
    message("  %s\n", drifts5->getDrift(i)->getDriftName().c_str());
  }
  delete drifts5;

  // Test case 6: order=4, ndim=1
  message("\nTest case 6: order=4, ndim=1\n");
  CovContext ctxt6(1, 1);
  DriftList* drifts6 = DriftFactory::createDriftListFromIRF(4, 0, ctxt6);
  message("Number of drifts: %d\n", drifts6->getNDrift());
  for (Id i = 0; i < drifts6->getNDrift(); i++)
  {
    message("  %s\n", drifts6->getDrift(i)->getDriftName().c_str());
  }
  delete drifts6;

  // Test case 7: order=2, ndim=3 (verify compatibility with old implementation)
  message("\nTest case 7: order=2, ndim=3 (backward compatibility test)\n");
  CovContext ctxt7(1, 3);
  DriftList* drifts7 = DriftFactory::createDriftListFromIRF(2, 0, ctxt7);
  message("Number of drifts: %d\n", drifts7->getNDrift());
  for (Id i = 0; i < drifts7->getNDrift(); i++)
  {
    message("  %s\n", drifts7->getDrift(i)->getDriftName().c_str());
  }
  delete drifts7;

  // Test case 8: order=1, ndim=3, nfex=2 (with external drifts)
  message("\nTest case 8: order=1, ndim=3, nfex=2\n");
  CovContext ctxt8(1, 3);
  DriftList* drifts8 = DriftFactory::createDriftListFromIRF(1, 2, ctxt8);
  message("Number of drifts: %d\n", drifts8->getNDrift());
  for (Id i = 0; i < drifts8->getNDrift(); i++)
  {
    message("  %s\n", drifts8->getDrift(i)->getDriftName().c_str());
  }
  delete drifts8;

  return 0;
}
