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

GSTLEARN_EXPORT void pile_reset(int type);
GSTLEARN_EXPORT void piles_reset(void);
GSTLEARN_EXPORT int pile_next(int type);
GSTLEARN_EXPORT void pile_manage(int type, int rank, int mode, char* ptr);
GSTLEARN_EXPORT int pile_correct(int type, int rank, int mode);
GSTLEARN_EXPORT char* pile_get(int type, int rank);
GSTLEARN_EXPORT void piles_dump(void);
