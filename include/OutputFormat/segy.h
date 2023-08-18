/******************************************************************************/
/*                                                                            */
/*                            gstlearn C++ Library                            */
/*                                                                            */
/* Copyright (c) (2023) MINES PARIS / ARMINES                                 */
/* Authors: gstlearn Team                                                     */
/* Website: https://gstlearn.org                                              */
/* License: BSD 3-clause                                                      */
/*                                                                            */
/******************************************************************************/
#pragma once

#include "gstlearn_export.hpp"
#include "Basic/NamingConvention.hpp"

#include "geoslib_d.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <iostream>

class Grid;

#ifndef SWIG
typedef enum
{
  SEGY_NUM    = 0,
  SEGY_ILINE  = 1,
  SEGY_XLINE  = 2,
  SEGY_XTRACE = 3,
  SEGY_YTRACE = 4,
  SEGY_ZMIN   = 5,
  SEGY_ZMAX   = 6,
  SEGY_VMIN   = 7,
  SEGY_VMAX   = 8,
  SEGY_THICK  = 9,
  SEGY_NB     = 10,
  SEGY_AUXTOP = 11,
  SEGY_AUXBOT = 12,
  SEGY_COUNT  = 13,
} ENUM_SEGY;

struct binaryFileHeader
{
  int JOB_ID;
  int LINE_NUM;
  int REEL_NUM;
  short int NUM_OF_TRACE;
  short int NUM_OF_AUX;
  short int INTERVAL_MS;
  short int INTERVAL_MS_ORI;
  unsigned short int NUM_OF_SAMPLES;
  unsigned short int NUM_OF_SAMPLES_ORI;
  short int SAMPLE_FORMAT;
  short int ENSEMBLE;
  short int TRACE_SORT;
  short int VERT_SUM;
  short int SWEEP_FREQ_START;
  short int SWEEP_FREQ_END;
  short int SWEEP_LENGTH;
  short int SWEEP_TYPE;
  short int SWEEP_NUM_CHANNEL;
  short int SWEEP_TAPER_LEN_START;
  short int SWEEP_TAPER_LEN_END;
  short int TAPER_TYPE;
  short int CORRELATED;
  short int BINARY_GAIN;
  short int AMP_RECOR;
  short int MEASURE_SYSTEM;
  short int IMPULSE_POLAR;
  short int POLAR_CODE;
  char UNNASSIGNED1[240];
  short int SEGY_REV_NUM;
  short int FIXED_LEN;
  short int NUM_EXT_HEAD;
  char UNNASSIGNED2[94];
};

struct traceHead
{
  int TRACE_SEQ_GLOBAL;
  int TRACE_SEQ_LOCAL;
  int ORI_RECORD_NUM;
  int TRACE_NUM_FIELD;
  int SOURCE_POINT;
  int ENSEMBLE_NUM;
  int ENS_TRACE_NUM;
  short int TRACE_CODE;
  short int NUM_VERT_SUM;
  short int NUM_HORZ_SUM;
  short int DATA_USE;
  int DIST_CENT_RECV;
  int RECV_GRP_ELEV;
  int SURF_ELEV_SRC;
  int SOURCE_DEPTH;
  int DATUM_ELEV_RECV;
  int DATUM_ELAV_SRC;
  int WATER_DEPTH_SRC;
  int WATER_DEPTH_GRP;
  short int SCALE_DEPTH;
  short int SCALE_COOR;
  int SRC_COOR_X;
  int SRC_COOR_Y;
  int GRP_COOR_X;
  int GRP_COOR_Y;
  short int COOR_UNIT;
  short int WEATHER_VEL;
  short int SWEATHER_VEL;
  short int UPHOLE_T_SRC;
  short int UPHOLE_T_GRP;
  short int SRC_STA_CORRC;
  short int GRP_STA_CORRC;
  short int TOTAL_STA;
  short int LAG_TIME_A;
  short int LAG_TIME_B;
  short int DELAY_T;
  short int MUTE_T_STRT;
  short int MUTE_T_END;
  unsigned short int NUM_OF_SAMPL;
  unsigned short int SAMPLE_INTRVL;
  short int GAIN_TYPE;
  short int GAIN_CONST;
  short int GAIN_INIT;
  short int CORRLTD;
  short int SWEEP_FREQ_START;
  short int SWEEP_FREQ_END;
  short int SWEEP_LENGTH;
  short int SWEEP_TYPE;
  short int SWEEP_TAPER_LEN_START;
  short int SWEEP_TAPER_LEN_END;
  short int TAPER_TYPE;
  short int ALIAS_FREQ;
  short int ALIAS_SLOPE;
  short int NOTCH_FREQ;
  short int NOTCH_SLOPE;
  short int LOWCUT_FREQ;
  short int HIGHCUT_FREQ;
  short int LOWCUT_SLOPE;
  short int HIGHCUT_SLOPE;
  short int YEAR;
  short int DAY;
  short int HOUR;
  short int MINUTE;
  short int SECOND;
  short int TIME_CODE;
  short int WEIGHT_FACT;
  short int GEOPHNE_ROLL;
  short int GEOPHNE_TRACE;
  short int GEOPHNE_LAST;
  short int GAP_SIZE;
  short int OVER_TRAVEL;
  int ENS_COOR_X;
  int ENS_COOR_Y;
  int INLINE;
  int CROSS;
  int SHOOTPOINT;
  short int SHOOTPOINT_SCALE;
  short int TRACE_UNIT;
  char TRANSD_CONST[6];
  short int TRANSD_UNIT;
  short int TRACE_IDENT;
  short int SCALE_TIME;
  short int SRC_ORIENT;
  char SRC_DIRECTION[6];
  char SRC_MEASUREMT[6];
  short int SRC_UNIT;
  char UNNASSIGNED1[6];
};
#endif

struct SegYArg
{
  int error;
  int ndescr;
  int npertrace;
  int ntraces;
  VectorVectorDouble tab;
  VectorVectorDouble descr;
  VectorVectorDouble cotes;
};

/***************************************/
/* Prototyping the functions in segy.c */
/***************************************/

GSTLEARN_EXPORT Grid segy_summary(const char *filesegy,
                                   DbGrid *surf2D = nullptr,
                                   const String &name_top = String(),
                                   const String &name_bot = String(),
                                   double thickmin = TEST,
                                   int option = 0,
                                   int nz_ss = ITEST,
                                   int verbose = 0,
                                   int iline_min = ITEST,
                                   int iline_max = ITEST,
                                   int xline_min = ITEST,
                                   int xline_max = ITEST,
                                   double modif_high = TEST,
                                   double modif_low = TEST,
                                   double modif_scale = TEST);
GSTLEARN_EXPORT SegYArg segy_array(const char *filesegy,
                                   DbGrid *surf2D = nullptr,
                                   const String& name_top = String(),
                                   const String& name_bot = String(),
                                   const String& top_aux  = String(),
                                   const String& bot_aux  = String(),
                                   double thickmin = TEST,
                                   int option = 0,
                                   int nz_ss = ITEST,
                                   int verbose = 0,
                                   int iline_min = ITEST,
                                   int iline_max = ITEST,
                                   int xline_min = ITEST,
                                   int xline_max = ITEST,
                                   double modif_high = TEST,
                                   double modif_low = TEST,
                                   double modif_scale = TEST);
GSTLEARN_EXPORT int db_segy(const char *filesegy,
                            DbGrid *grid3D,
                            DbGrid *surf2D = nullptr,
                            const String &name_top = String(),
                            const String &name_bot = String(),
                            double thickmin = TEST,
                            int option = 0,
                            int nz_ss = ITEST,
                            int verbose = 0,
                            int iline_min = ITEST,
                            int iline_max = ITEST,
                            int xline_min = ITEST,
                            int xline_max = ITEST,
                            double modif_high = TEST,
                            double modif_low = TEST,
                            double modif_scale = TEST,
                            const NamingConvention& namconv = NamingConvention("SEGY"));

