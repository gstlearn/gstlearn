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
#include "Basic/NamingConvention.hpp"
#include "Db/DbGrid.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>

namespace gstlrn
{
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
  Id JOB_ID;
  Id LINE_NUM;
  Id REEL_NUM;
  short NUM_OF_TRACE;
  short NUM_OF_AUX;
  short INTERVAL_MS;
  short INTERVAL_MS_ORI;
  unsigned short NUM_OF_SAMPLES;
  unsigned short NUM_OF_SAMPLES_ORI;
  short SAMPLE_FORMAT;
  short ENSEMBLE;
  short TRACE_SORT;
  short VERT_SUM;
  short SWEEP_FREQ_START;
  short SWEEP_FREQ_END;
  short SWEEP_LENGTH;
  short SWEEP_TYPE;
  short SWEEP_NUM_CHANNEL;
  short SWEEP_TAPER_LEN_START;
  short SWEEP_TAPER_LEN_END;
  short TAPER_TYPE;
  short CORRELATED;
  short BINARY_GAIN;
  short AMP_RECOR;
  short MEASURE_SYSTEM;
  short IMPULSE_POLAR;
  short POLAR_CODE;
  char UNNASSIGNED1[240];
  short SEGY_REV_NUM;
  short FIXED_LEN;
  short NUM_EXT_HEAD;
  char UNNASSIGNED2[94];
};

struct traceHead
{
  Id TRACE_SEQ_GLOBAL;
  Id TRACE_SEQ_LOCAL;
  Id ORI_RECORD_NUM;
  Id TRACE_NUM_FIELD;
  Id SOURCE_POINT;
  Id ENSEMBLE_NUM;
  Id ENS_TRACE_NUM;
  short TRACE_CODE;
  short NUM_VERT_SUM;
  short NUM_HORZ_SUM;
  short DATA_USE;
  Id DIST_CENT_RECV;
  Id RECV_GRP_ELEV;
  Id SURF_ELEV_SRC;
  Id SOURCE_DEPTH;
  Id DATUM_ELEV_RECV;
  Id DATUM_ELAV_SRC;
  Id WATER_DEPTH_SRC;
  Id WATER_DEPTH_GRP;
  short SCALE_DEPTH;
  short SCALE_COOR;
  Id SRC_COOR_X;
  Id SRC_COOR_Y;
  Id GRP_COOR_X;
  Id GRP_COOR_Y;
  short COOR_UNIT;
  short WEATHER_VEL;
  short SWEATHER_VEL;
  short UPHOLE_T_SRC;
  short UPHOLE_T_GRP;
  short SRC_STA_CORRC;
  short GRP_STA_CORRC;
  short TOTAL_STA;
  short LAG_TIME_A;
  short LAG_TIME_B;
  short DELAY_T;
  short MUTE_T_STRT;
  short MUTE_T_END;
  unsigned short NUM_OF_SAMPL;
  unsigned short SAMPLE_INTRVL;
  short GAIN_TYPE;
  short GAIN_CONST;
  short GAIN_INIT;
  short CORRLTD;
  short SWEEP_FREQ_START;
  short SWEEP_FREQ_END;
  short SWEEP_LENGTH;
  short SWEEP_TYPE;
  short SWEEP_TAPER_LEN_START;
  short SWEEP_TAPER_LEN_END;
  short TAPER_TYPE;
  short ALIAS_FREQ;
  short ALIAS_SLOPE;
  short NOTCH_FREQ;
  short NOTCH_SLOPE;
  short LOWCUT_FREQ;
  short HIGHCUT_FREQ;
  short LOWCUT_SLOPE;
  short HIGHCUT_SLOPE;
  short YEAR;
  short DAY;
  short HOUR;
  short MINUTE;
  short SECOND;
  short TIME_CODE;
  short WEIGHT_FACT;
  short GEOPHNE_ROLL;
  short GEOPHNE_TRACE;
  short GEOPHNE_LAST;
  short GAP_SIZE;
  short OVER_TRAVEL;
  Id ENS_COOR_X;
  Id ENS_COOR_Y;
  Id INLINE;
  Id CROSS;
  Id SHOOTPOINT;
  short SHOOTPOINT_SCALE;
  short TRACE_UNIT;
  char TRANSD_CONST[6];
  short TRANSD_UNIT;
  short TRACE_IDENT;
  short SCALE_TIME;
  short SRC_ORIENT;
  char SRC_DIRECTION[6];
  char SRC_MEASUREMT[6];
  short SRC_UNIT;
  char UNNASSIGNED1[6];
};
#endif

struct SegYArg
{
  Id error;
  Id ndescr;
  Id npertrace;
  Id ntraces;
  VectorVectorDouble tab;
  VectorVectorDouble descr;
  VectorVectorDouble cotes;
};

/***************************************/
/* Prototyping the functions in segy.c */
/***************************************/

GSTLEARN_EXPORT Grid segy_summary(const char *filesegy,
                                  DbGrid *surf2D = nullptr,
                                  const String &name_top = "",
                                  const String &name_bot = "",
                                  double thickmin = TEST,
                                  Id option = 0,
                                  Id nz_ss = ITEST,
                                  Id verbOption = 1,
                                  Id iline_min = ITEST,
                                  Id iline_max = ITEST,
                                  Id xline_min = ITEST,
                                  Id xline_max = ITEST,
                                  double modif_high = TEST,
                                  double modif_low = TEST,
                                  double modif_scale = TEST,
                                  Id codefmt = 1);
GSTLEARN_EXPORT SegYArg segy_array(const char *filesegy,
                                   DbGrid *surf2D = nullptr,
                                   const String& top_name = "",
                                   const String& bot_name = "",
                                   const String& top_aux  = "",
                                   const String& bot_aux  = "",
                                   double thickmin = TEST,
                                   Id option = 0,
                                   Id nz_ss = ITEST,
                                   Id verbOption = 0,
                                   Id iline_min = ITEST,
                                   Id iline_max = ITEST,
                                   Id xline_min = ITEST,
                                   Id xline_max = ITEST,
                                   double modif_high = TEST,
                                   double modif_low = TEST,
                                   double modif_scale = TEST,
                                   Id codefmt = 1);
GSTLEARN_EXPORT Id db_segy(const char *filesegy,
                            DbGrid *grid3D,
                            DbGrid *surf2D = nullptr,
                            const String &name_top = "",
                            const String &name_bot = "",
                            double thickmin = TEST,
                            Id option = 0,
                            Id nz_ss = ITEST,
                            Id verbOption = 0,
                            Id iline_min = ITEST,
                            Id iline_max = ITEST,
                            Id xline_min = ITEST,
                            Id xline_max = ITEST,
                            double modif_high = TEST,
                            double modif_low = TEST,
                            double modif_scale = TEST,
                            Id codefmt = 1,
                            const NamingConvention& namconv = NamingConvention("SEGY"));

}
