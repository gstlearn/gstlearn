/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://github.com/gstlearn                                       */
/* License: BSD 3 clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "Enum/AEnum.hpp"

#define ENUM_CONVTYPE EConvType, UNKNOWN,\
                 UNKNOWN,    -1, "Unknown",\
                 UNIFORM,     0, "Uniform",\
                 EXPONENTIAL, 1, "Exponential",\
                 GAUSSIAN,    2, "Gaussian",\
                 SINCARD,     3, "Cardinal Sine"

ENUM_DECLARE(ENUM_CONVTYPE)
