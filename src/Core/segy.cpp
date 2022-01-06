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
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "geoslib_enum.h"
#include "Basic/Utilities.hpp"
#include "Basic/File.hpp"
#include "Db/Db.hpp"
#include "segy.h"

struct RefPt
{
  int iline;
  int xline;
  double xtrace;
  double ytrace;
};

struct RefStats
{
  int nbtrace;
  int nbtrace_in;
  int nbtrace_def;
  int nbvalue_in;

  double zminl;
  double zmaxl;
  double vminl;
  double vmaxl;
  double thicl;

  int ilming;
  int ilmaxg;
  int xlming;
  int xlmaxg;

  double xming;
  double xmaxg;
  double yming;
  double ymaxg;
  double zming;
  double zmaxg;
  double vming;
  double vmaxg;
  double thicg;

  double modif_low;
  double modif_high;
  double modif_scale;
};

#define MAT(i,j)            (mat[2 * (i) + (j)])
#define SWAP_INT16(x)	(x) = ((0x00ff & ((x))>>8) | (0xff00 & ((x))<<8))
#define SWAP_INT32(x)	(x) = ((0x000000ff & ((x))>>24) | \
                             (0x0000ff00 & ((x))>>8) |  \
                             (0x00ff0000 & ((x))<<8) |  \
                             (0xff000000 & ((x))<<24))

/****************************************************************************/
/*!
 ** Routines for reading / decoding information stored in SEGY file
 **
 *****************************************************************************/
static int st_to_f(int x)
{
  int y;
  y = SWAP_INT32(x);
  return y;
}

static short int st_to_s(short int x)
{
  short int y;
  y = SWAP_INT16(x);
  return y;
}

static unsigned short int st_to_u(unsigned short int a)
{
  unsigned short int tmp = a >> 8;
  unsigned short int b;

  b = (a << 8) | (tmp);
  return b;
}

static int st_to_i(int a)
{
  unsigned short int tmp1 = (a >> 16);
  unsigned short int tmp2 = (a * 0x0000FFFF);
  tmp2 = st_to_u(tmp2);
  tmp1 = st_to_u(tmp1);

  int b = (int) tmp2;
  b = b << 16;
  b = b | (int) tmp1;

  return b;
}

static float st_ibm2ieee(const float ibm)
{

  // DEFINE AND SETUP UNION VARIABLES
  union value
  {
    float f;
    unsigned char c[4];
  };
  value src;

  // ASSIGN PARAMETER TO UNION
  src.f = ibm;

  // CONVERT TO FLOAT
  long IntMantissa = ((long) src.c[1] << 16) + ((long) src.c[2] << 8)
                     + src.c[3];

  float Mantissa = float(IntMantissa) / float(0x1000000);
  float PosResult = Mantissa
      * (float) pow(16.0, double((src.c[0] & 0x7F) - 64));

  if (src.c[0] & 0x80)
    return -PosResult;
  else
    return PosResult;
}

/****************************************************************************/
/*!
 ** Internal function to scale the coordinate by the scaling factor
 **
 ** \return Scaled value
 **
 ** \param[in]  coor      Integer value containing the input coordinate
 ** \param[in]  scale     Integer value giving the scale factor
 **
 *****************************************************************************/
static double st_scaling(int coor, int scale)
{
  double value = coor;
  if (scale < 0)
    value /= (double) ABS(scale);
  else
    value *= (double) ABS(scale);
  return value;
}

/****************************************************************************/
/*!
 ** Internal function for printing the Trace Header
 **
 ** \param[in]  Theader   Pointer to the Trace Header contents
 ** \param[in]  numTrace  Rank of the trace
 **
 *****************************************************************************/
static void st_print_traceHead(traceHead *Theader, int numTrace)
{
  message("\nTrace Header #%d\n", numTrace);
  message("TRACE_SEQ_GLOBAL    : %i \n", st_to_i(Theader->TRACE_SEQ_GLOBAL));
  message("TRACE_SEQ_LOCAL     : %i \n", st_to_i(Theader->TRACE_SEQ_LOCAL));
  message("ORI_RECORD_NUM      : %i \n", st_to_i(Theader->ORI_RECORD_NUM));
  message("TRACE_NUM_FIELD     : %i \n", st_to_i(Theader->TRACE_NUM_FIELD));
  message("SOURCE_POINT        : %i \n", st_to_i(Theader->SOURCE_POINT));
  message("ENSEMBLE_NUM        : %i \n", st_to_i(Theader->ENSEMBLE_NUM));
  message("ENS_TRACE_NUM       : %i \n", st_to_i(Theader->ENS_TRACE_NUM));
  message("TRACE_CODE          : %hi\n", st_to_s(Theader->TRACE_CODE));
  message("NUM_VERT_SUM        : %hi\n", st_to_s(Theader->NUM_VERT_SUM));
  message("NUM_HORZ_SUM        : %hi\n", st_to_s(Theader->NUM_HORZ_SUM));
  message("DATA_USE            : %hi\n", st_to_s(Theader->DATA_USE));
  message("DIST_CENT_RECV      : %i \n", st_to_i(Theader->DIST_CENT_RECV));
  message("RECV_GRP_ELEV       : %i \n", st_to_i(Theader->RECV_GRP_ELEV));
  message("SURF_ELEV_SRC       : %i \n", st_to_i(Theader->SURF_ELEV_SRC));
  message("SOURCE_DEPTH        : %i \n", st_to_i(Theader->SOURCE_DEPTH));
  message("DATUM_ELEV_RECV     : %i \n", st_to_i(Theader->DATUM_ELEV_RECV));
  message("DATUM_ELAV_SRC      : %i \n", st_to_i(Theader->DATUM_ELAV_SRC));
  message("WATER_DEPTH_SRC     : %i \n", st_to_i(Theader->WATER_DEPTH_SRC));
  message("WATER_DEPTH_GRP     : %i \n", st_to_i(Theader->WATER_DEPTH_GRP));
  message("SCALE_DEPTH         : %hi\n", st_to_s(Theader->SCALE_DEPTH));
  message("SCALE_COOR          : %hi\n", st_to_s(Theader->SCALE_COOR));
  message("SRC_COOR_X          : %i \n", st_to_i(Theader->SRC_COOR_X));
  message("SRC_COOR_Y          : %i \n", st_to_i(Theader->SRC_COOR_Y));
  message("GRP_COOR_X          : %i \n", st_to_i(Theader->GRP_COOR_X));
  message("GRP_COOR_Y          : %i \n", st_to_i(Theader->GRP_COOR_Y));
  message("COOR_UNIT           : %hi\n", st_to_s(Theader->COOR_UNIT));
  message("WEATHER_VEL         : %hi\n", st_to_s(Theader->WEATHER_VEL));
  message("SWEATHER_VEL        : %hi\n", st_to_s(Theader->SWEATHER_VEL));
  message("UPHOLE_T_SRC        : %hi\n", st_to_s(Theader->UPHOLE_T_SRC));
  message("UPHOLE_T_GRP        : %hi\n", st_to_s(Theader->UPHOLE_T_GRP));
  message("SRC_STA_CORRC       : %hi\n", st_to_s(Theader->SRC_STA_CORRC));
  message("GRP_STA_CORRC       : %hi\n", st_to_s(Theader->GRP_STA_CORRC));
  message("TOTAL_STA           : %hi\n", st_to_s(Theader->TOTAL_STA));
  message("LAG_TIME_A          : %hi\n", st_to_s(Theader->LAG_TIME_A));
  message("LAG_TIME_B          : %hi\n", st_to_s(Theader->LAG_TIME_B));
  message("DELAY_T             : %hi\n", st_to_s(Theader->DELAY_T));
  message("MUTE_T_STRT         : %hi\n", st_to_s(Theader->MUTE_T_STRT));
  message("MUTE_T_END          : %hi\n", st_to_s(Theader->MUTE_T_END));
  message("NUM_OF_SAMPL        : %hi\n", st_to_s(Theader->NUM_OF_SAMPL));
  message("SAMPLE_INTRVL       : %hi\n", st_to_s(Theader->SAMPLE_INTRVL));
  message("GAIN_TYPE           : %hi\n", st_to_s(Theader->GAIN_TYPE));
  message("GAIN_CONST          : %hi\n", st_to_s(Theader->GAIN_CONST));
  message("GAIN_INIT           : %hi\n", st_to_s(Theader->GAIN_INIT));
  message("CORRLTD             : %hi\n", st_to_s(Theader->CORRLTD));
  message("SWEEP_FREQ_START    : %hi\n", st_to_s(Theader->SWEEP_FREQ_START));
  message("SWEEP_FREQ_END      : %hi\n", st_to_s(Theader->SWEEP_FREQ_END));
  message("SWEEP_LENGTH        : %hi\n", st_to_s(Theader->SWEEP_LENGTH));
  message("SWEEP_TYPE          : %hi\n", st_to_s(Theader->SWEEP_TYPE));
  message("SWEEP_TAPER_LEN_STT : %hi\n",
          st_to_s(Theader->SWEEP_TAPER_LEN_START));
  message("SWEEP_TAPER_LEN_END : %hi\n", st_to_s(Theader->SWEEP_TAPER_LEN_END));
  message("TAPER_TYPE          : %hi\n", st_to_s(Theader->TAPER_TYPE));
  message("ALIAS_FREQ          : %hi\n", st_to_s(Theader->ALIAS_FREQ));
  message("ALIAS_SLOPE         : %hi\n", st_to_s(Theader->ALIAS_SLOPE));
  message("NOTCH_FREQ          : %hi\n", st_to_s(Theader->NOTCH_FREQ));
  message("NOTCH_SLOPE         : %hi\n", st_to_s(Theader->NOTCH_SLOPE));
  message("LOWCUT_FREQ         : %hi\n", st_to_s(Theader->LOWCUT_FREQ));
  message("HIGHCUT_FREQ        : %hi\n", st_to_s(Theader->HIGHCUT_FREQ));
  message("LOWCUT_SLOPE        : %hi\n", st_to_s(Theader->LOWCUT_SLOPE));
  message("HIGHCUT_SLOPE       : %hi\n", st_to_s(Theader->HIGHCUT_SLOPE));
  message("YEAR                : %hi\n", st_to_s(Theader->YEAR));
  message("DAY                 : %hi\n", st_to_s(Theader->DAY));
  message("HOUR                : %hi\n", st_to_s(Theader->HOUR));
  message("MINUTE              : %hi\n", st_to_s(Theader->MINUTE));
  message("SECOND              : %hi\n", st_to_s(Theader->SECOND));
  message("TIME_CODE           : %hi\n", st_to_s(Theader->TIME_CODE));
  message("WEIGHT_FACT         : %hi\n", st_to_s(Theader->WEIGHT_FACT));
  message("GEOPHONE_ROLL       : %hi\n", st_to_s(Theader->GEOPHNE_ROLL));
  message("GEOPHONE_TRACE      : %hi\n", st_to_s(Theader->GEOPHNE_TRACE));
  message("GEOPHONE_LAST       : %hi\n", st_to_s(Theader->GEOPHNE_LAST));
  message("GAP_SIZE            : %hi\n", st_to_s(Theader->GAP_SIZE));
  message("OVER_TRAVEL         : %hi\n", st_to_s(Theader->OVER_TRAVEL));
  message("ENS_COOR_X          : %i \n", st_to_f(Theader->ENS_COOR_X));
  message("ENS_COOR_Y          : %i \n", st_to_f(Theader->ENS_COOR_Y));
  message("INLINE              : %i \n", st_to_i(Theader->INLINE));
  message("CROSS               : %i \n", st_to_i(Theader->CROSS));
  message("SHOOTPOINT          : %i \n", st_to_i(Theader->SHOOTPOINT));
  message("SHOOTPOINT_SCALE    : %hi\n", st_to_s(Theader->SHOOTPOINT_SCALE));
  message("TRACE_UNIT          : %hi\n", st_to_s(Theader->TRACE_UNIT));
  message("TRANSD_UNIT         : %hi\n", st_to_s(Theader->TRANSD_UNIT));
  message("TRACE_IDENT         : %hi\n", st_to_s(Theader->TRACE_IDENT));
  message("SCALE_TIME          : %hi\n", st_to_s(Theader->SCALE_TIME));
  message("SRC_ORIENT          : %hi\n", st_to_s(Theader->SRC_ORIENT));
  message("SRC_UNIT            : %hi\n", st_to_s(Theader->SRC_UNIT));
}

/****************************************************************************/
/*!
 ** Internal function for printing the Binary File Header
 **
 ** \param[in]  Bheader   Pointer to the Binary File header
 **
 *****************************************************************************/
static void st_print_BFileHead(binaryFileHeader *Bheader)
{
  message("\nBinary File Header\n");
  message("JOB_ID               :%i  \n", st_to_i(Bheader->JOB_ID));
  message("LINE_NUM             :%i  \n", st_to_i(Bheader->LINE_NUM));
  message("REEL_NUM             :%i  \n", st_to_i(Bheader->REEL_NUM));
  message("NUM_OF_TRACE         :%hi \n", st_to_s(Bheader->NUM_OF_TRACE));
  message("NUM_OF_AUX           :%hi \n", st_to_s(Bheader->NUM_OF_AUX));
  message("INTERVAL_MS          :%hi \n", st_to_s(Bheader->INTERVAL_MS));
  message("INTERVAL_MS_ORI      :%hi \n", st_to_s(Bheader->INTERVAL_MS_ORI));
  message("NUM_OF_SAMPLES       :%hi \n", st_to_s(Bheader->NUM_OF_SAMPLES));
  message("NUM_OF_SAMPLES_ORI   :%hi \n", st_to_s(Bheader->NUM_OF_SAMPLES_ORI));
  message("SAMPLE_FORMAT        :%hi \n", st_to_s(Bheader->SAMPLE_FORMAT));
  message("ENSEMBLE             :%hi \n", st_to_s(Bheader->ENSEMBLE));
  message("TRACE_SORT           :%hi \n", st_to_s(Bheader->TRACE_SORT));
  message("VERT_SUM             :%hi \n", st_to_s(Bheader->VERT_SUM));
  message("SWEEP_FREQ_START     :%hi \n", st_to_s(Bheader->SWEEP_FREQ_START));
  message("SWEEP_FREQ_END       :%hi \n", st_to_s(Bheader->SWEEP_FREQ_END));
  message("SWEEP_LENGTH         :%hi \n", st_to_s(Bheader->SWEEP_LENGTH));
  message("SWEEP_TYPE           :%hi \n", st_to_s(Bheader->SWEEP_TYPE));
  message("SWEEP_NUM_CHANNEL    :%hi \n", st_to_s(Bheader->SWEEP_NUM_CHANNEL));
  message("SWEEP_TAPER_LEN_START:%hi \n",
          st_to_s(Bheader->SWEEP_TAPER_LEN_START));
  message("SWEEP_TAPER_LEN_END  :%hi \n",
          st_to_s(Bheader->SWEEP_TAPER_LEN_END));
  message("TAPER_TYPE           :%hi \n", st_to_s(Bheader->TAPER_TYPE));
  message("CORRELATED           :%hi \n", st_to_s(Bheader->CORRELATED));
  message("BINARY_GAIN          :%hi \n", st_to_s(Bheader->BINARY_GAIN));
  message("AMP_RECOR            :%hi \n", st_to_s(Bheader->AMP_RECOR));
  message("MEASURE_SYSTEM       :%hi \n", st_to_s(Bheader->MEASURE_SYSTEM));
  message("IMPULSE_POLAR        :%hi \n", st_to_s(Bheader->IMPULSE_POLAR));
  message("POLAR_CODE           :%hi \n", st_to_s(Bheader->POLAR_CODE));
  message("SEGY_REV_NUM         :%hi \n", st_to_s(Bheader->SEGY_REV_NUM));
  message("FIXED_LEN            :%hi \n", st_to_s(Bheader->FIXED_LEN));
  message("NUM_EXT_HEAD         :%hi \n", st_to_s(Bheader->NUM_EXT_HEAD));
}

/****************************************************************************/
/*!
 ** Internal function for reading the File header
 **
 ** \param[in]  file       Pointer to the File stream
 ** \param[in]  verbOption Verbose option
 **
 ** \param[out] NPerTrace   Number of samples per trace
 ** \param[out] delta    Distance between two consecutive vertical samples
 **
 *****************************************************************************/
static int st_readFileHeader(FILE *file,
                             int verbOption,
                             int *NPerTrace,
                             double *delta)
{
  unsigned char TFileHead_[3200];
  binaryFileHeader BFileHead_;

  if (fread(&TFileHead_, 1, sizeof(TFileHead_), file) == 0) return (1);
  if (fread(&BFileHead_, 1, sizeof(BFileHead_), file) == 0) return (1);
  if (verbOption >= 2) st_print_BFileHead(&BFileHead_);
  *NPerTrace = st_to_s(BFileHead_.NUM_OF_SAMPLES);
  *delta = (double) st_to_s(BFileHead_.INTERVAL_MS) / 1000;
  return (0);
}

static int st_read_trace(FILE *file,
                         int code,
                         int numtrace,
                         int nPerTrace,
                         double delta,
                         int verbOption,
                         VectorDouble &values,
                         VectorDouble &cotes)
{
  short svalue;
  int minsamp, maxsamp, ivalue, imaxvalue, iminvalue;
  float fvalue, fmaxvalue, fminvalue, cote;

  cote       = 0;
  ivalue     = 0;
  imaxvalue  = 0;
  iminvalue  = 0;
  fminvalue  = 0.0;
  fmaxvalue  = 0.0;
  minsamp    = -1;
  maxsamp    = -1;

  if (verbOption >= 2) message("Trace %5d: ", numtrace);
  if (code == 1 || code == 5)
  {
    for (int i = 0; i < nPerTrace; i++)
    {
      if (fread(&fvalue, 4, 1, file) == 0) return (1);

      if (code == 1) fvalue = st_ibm2ieee(fvalue);

      if (fvalue > fmaxvalue)
      {
        fmaxvalue = fvalue;
        maxsamp = i + 1;
      }
      if (fvalue < fminvalue)
      {
        fminvalue = fvalue;
        minsamp = i + 1;
      }
      values[i] = (double) fvalue;
      if (ABS(values[i]) < 0.01) values[i] = TEST;
      cotes[i] = (double) cote;
      cote -= (float) delta;
    }
    if (verbOption >= 2)
      message("Min(%6d) = %11.0f - Max(%6d) = %11.0f\n", minsamp, fminvalue,
              maxsamp, fmaxvalue);
  }

  if (code == 2 || code == 3)
  {
    for (int i = 0; i < nPerTrace; i++)
    {
      if (code == 2)
      {
        if (fread(&ivalue, 4, 1, file) == 0) return (1);
        ivalue = st_to_f(ivalue);
      }
      if (code == 3)
      {
        if (fread(&svalue, 2, 1, file) == 0) return (1);
        ivalue = (int) svalue;
      }

      if (ivalue > imaxvalue)
      {
        imaxvalue = ivalue;
        maxsamp = i + 1;
      }
      if (ivalue < iminvalue)
      {
        iminvalue = ivalue;
        minsamp = i + 1;
      }
      values[i] = (double) ivalue;
      if (ABS(values[i]) < 0.01) values[i] = TEST;
      cotes[i] = (double) cote;
      cote -= (float) delta;
    }
    if (verbOption >= 2)
      message("Min(%6d) = %11d - Max(%6d) = %11d\n", minsamp, iminvalue,
              maxsamp, imaxvalue);
  }

  return (minsamp < 0 || maxsamp < 0);
}

/****************************************************************************/
/*!
 ** Identify Bottom and Top surfaces (if present)
 **
 ** \return Error returned code
 **
 ** \param[in] verbOption Verbose Option
 ** \param[in]  surfaces  Db containing the top, Bottom and Reference surfaces
 **                       This file is optional
 ** \param[in]  name_bot  Name of variable containing the Bottom Surface (or empty)
 ** \param[in]  flag_bot  Flag for defining a Bottom surface
 ** \param[in]  name_top  Name of variable containing the Top Surface (or empty)
 ** \param[in]  flag_top  Flag for defining a Top surface
 **
 ** \param[out] iatt_top  Attribute index for the top surface
 ** \param[out] iatt_bot  Attribute index for the Bottom surface
 **
 *****************************************************************************/
static int st_surface_identify(int verbOption,
                               Db *surfaces,
                               const String &name_bot,
                               int flag_bot,
                               int *iatt_bot,
                               const String &name_top,
                               int flag_top,
                               int *iatt_top)
{
  *iatt_bot = -1;

  if (flag_bot)
  {
    if (name_bot.empty())
    {
      messerr("When flattening using Bottom surface");
      messerr("you must provide a Surface file and a valid Bottom variable");
      return 1;
    }
    *iatt_bot = surfaces->getAttribute(name_bot);
    if (*iatt_bot < 0) return 1;
  }

  *iatt_top = -1;
  if (flag_top)
  {
    if (name_top.empty())
    {
      messerr("When flattening using Top surface");
      messerr("you must provide a Surface file and a valid Top variable");
      return 1;
    }
    *iatt_top = surfaces->getAttribute(name_top);
    if (*iatt_top < 0) return 1;
  }

  if (verbOption && (flag_bot || flag_top))
  {
    mestitle(2, "Horizontalization:");
    if (flag_top) message("- Top surface: %s\n", name_top.c_str());
    if (flag_bot) message("- Bottom surface: %s\n", name_bot.c_str());
  }
  return 0;
}

/****************************************************************************/
/*!
 ** Define the cutoffs along a vertical for a given trace
 **
 ** \return Error returned code
 **
 ** \param[in]  surfaces  Db containing the top, Bottom and Reference surfaces
 **                       This file is optional
 ** \param[in]  iatt_top  Rank of attribute containing the Top Surface (or 0)
 ** \param[in]  iatt_bot  Arnk of the attribute containing the Bottom Surface (or 0)
 ** \param[in]  xtrace    Coordinate of the trace along X
 ** \param[in]  ytrace    Coordinate of the trace along Y
 ** \param[in]  thickmin  Minimum thickness (if defined)
 **
 ** \param[out] cztop     Coordinate for the Top along Z (or TEST)
 ** \param[out] czbot     Coordinate for the Bottom along Z (or TEST)
 **
 *****************************************************************************/
static int st_get_cuts(Db *surfaces,
                       int iatt_top,
                       int iatt_bot,
                       double xtrace,
                       double ytrace,
                       double thickmin,
                       double *cztop,
                       double *czbot)
{
  int rank;

  // Initializations

  *cztop = *czbot = TEST;

  // Check if the surface file has been defined

  if (surfaces == nullptr) return (0);
  if (iatt_bot < 0 && iatt_top < 0) return (0);

  // The surface is valid, return the vertical bounds along trace

  VectorDouble coor(2);
  coor[0] = xtrace;
  coor[1] = ytrace;

  // Find the index of grid node

  rank = surfaces->coordinateToRank(coor);
  if (rank < 0) return (1);

  // Check if the top has been defined

  if (iatt_top >= 0) *cztop = surfaces->getArray(rank, iatt_top);

  // Check if the top has been defined

  if (iatt_bot >= 0) *czbot = surfaces->getArray(rank, iatt_bot);

  if (FFFF(*cztop) && FFFF(*czbot)) return (1);
  if (!FFFF(*cztop) && !FFFF(*czbot) && (*cztop) < (*czbot)) return (1);

  if (thickmin > 0.)
  {
    double thick = (*cztop) - (*czbot);
    if (thick < thickmin) return 1;
  }
  return (0);
}

/****************************************************************************/
/*!
 ** Check that calculations are correct
 **
 ** \param[in]  refpt     Array of RefPt structure pointers
 ** \param[in]  refstats  Structure for Statistics
 ** \param[in]  nbrefpt   Number of reference points
 ** \param[in]  x0        Origin of the grid along X
 ** \param[in]  y0        Origin of the grid along Y
 ** \param[in]  dx        Mesh of the grid along X
 ** \param[in]  dy        Mesh of the grid along Y
 ** \param[in]  sint      Sine of the rotation angle
 ** \param[in]  cost      Cosine of the rotation angle
 **
 ** \remarks The internal flag debug must be set to 1 to perform this task
 **
 *****************************************************************************/
static void st_verify_refpt(RefPt refpt[3],
                            RefStats &refstats,
                            int nbrefpt,
                            double x0,
                            double y0,
                            double dx,
                            double dy,
                            double sint,
                            double cost)
{
  int di, dj;
  double xn, yn;
  static int debug = 0;

  if (!debug) return;
  for (int i = 0; i < nbrefpt; i++)
  {
    di = refpt[i].iline - refstats.ilming;
    dj = refpt[i].xline - refstats.xlming;
    xn = x0 + di * dx * cost - dj * dy * sint;
    yn = y0 + di * dx * sint + dj * dy * cost;
    message("- Reference Point #%d / %d\n", i + 1, nbrefpt);
    message("  . Inline        = %d\n", refpt[i].iline);
    message("  . Crossline     = %d\n", refpt[i].xline);
    message("  . Trace along X = %12.4lf (Pred. = %12.4lf)\n", refpt[i].xtrace,
            xn);
    message("  . Trace along Y = %12.4lf (Pred. = %12.4lf)\n", refpt[i].ytrace,
            yn);
  }
}

/****************************************************************************/
/*!
 ** Derive grid information from the three Reference Points
 **
 ** \param[in]  refpt     Array of RefPt structure pointers
 ** \param[in]  refstats  Structure for Statistics
 ** \param[in]  dz        Vertical mesh
 **
 ** \param[out] def_grid  Grid output structure
 **
 *****************************************************************************/
static void st_grid_from_3refpt(RefPt refpt[3],
                                RefStats &refstats,
                                double dz,
                                Grid &def_grid)
{
  double dx12, dy12, di12, dj12, dx13, dy13, di13, dj13, mat[4], a[2], di10,
      dj10;
  double dxs, dxc, dys, dyc, dx, dy, theta, det, x0, y0, cost, sint;
  int nz;

  // Comparing points 1 and 2

  dx12 = refpt[1].xtrace - refpt[0].xtrace;
  dy12 = refpt[1].ytrace - refpt[0].ytrace;
  di12 = refpt[1].iline - refpt[0].iline;
  dj12 = refpt[1].xline - refpt[0].xline;

  // Comparing points 1 and 3

  dx13 = refpt[2].xtrace - refpt[0].xtrace;
  dy13 = refpt[2].ytrace - refpt[0].ytrace;
  di13 = refpt[2].iline - refpt[0].iline;
  dj13 = refpt[2].xline - refpt[0].xline;

  // Fill the matrix and vector

  MAT(0,0) = +di12;
  MAT(0,1) = -dj12;
  MAT(1,0) = +di13;
  MAT(1,1) = -dj13;
  a[0] = dx12;
  a[1] = dx13;

  det = (MAT(0,0) * MAT(1, 1) - MAT(1,0) * MAT(0, 1));
  dxc = (a[0] * MAT(1, 1) - a[1] * MAT(0, 1)) / det;
  dys = (MAT(0,0) * a[1] - MAT(1,0) * a[0]) / det;

  MAT(0,0) = +di12;
  MAT(0,1) = +dj12;
  MAT(1,0) = +di13;
  MAT(1,1) = +dj13;
  a[0] = dy12;
  a[1] = dy13;

  det = (MAT(0,0) * MAT(1, 1) - MAT(1,0) * MAT(0, 1));
  dxs = (a[0] * MAT(1, 1) - a[1] * MAT(0, 1)) / det;
  dyc = (MAT(0,0) * a[1] - MAT(1,0) * a[0]) / det;

  dx = sqrt(dxs * dxs + dxc * dxc);
  dy = sqrt(dys * dys + dyc * dyc);
  theta = 0.5 * (atan2(dxs, dxc) + atan2(dys, dyc)) * 180. / GV_PI;
  cost = 0.5 * (dxc / dx + dyc / dy);
  sint = 0.5 * (dxs / dx + dys / dy);

  // Origin of the grid

  di10 = refpt[0].iline - refstats.ilming;
  dj10 = refpt[0].xline - refstats.xlming;
  x0 = refpt[0].xtrace - di10 * dx * cost + dj10 * dy * sint;
  y0 = refpt[0].ytrace - di10 * dx * sint - dj10 * dy * cost;
  nz = static_cast<int>((refstats.zmaxl - refstats.zminl) / dz);

  // Fill the Grid structure

  def_grid = Grid(3);
  def_grid.setNX(0, refstats.ilmaxg - refstats.ilming + 1);
  def_grid.setNX(1, refstats.xlmaxg - refstats.xlming + 1);
  def_grid.setNX(2, nz);
  def_grid.setX0(0, x0);
  def_grid.setX0(1, y0);
  def_grid.setX0(2, refstats.zminl);
  def_grid.setDX(0, dx);
  def_grid.setDX(1, dy);
  def_grid.setDX(2, dz);
  def_grid.setRotationByAngle(theta);

  // Check that reference traces are well honored

  st_verify_refpt(refpt, refstats, 3, x0, y0, dx, dy, sint, cost);
}

/****************************************************************************/
/*!
 ** Derive grid information from the two Reference Points
 **
 ** \param[in]  refpt     Array of RefPt structure pointers
 ** \param[in]  refstats  Structure for Statistics
 ** \param[in]  dz        Vertical mesh
 **
 ** \param[out] def_grid  Grid output structure
 **
 *****************************************************************************/
static void st_grid_from_2refpt(RefPt refpt[3],
                                RefStats &refstats,
                                double dz,
                                Grid &def_grid)
{
  double dx12, dy12, di12, dj12, di10, dj10;
  double dx, dy, theta, x0, y0, sint, cost;
  int nz;

  // Initializations

  cost = 0.;
  sint = 0.;
  dx = 0.;
  dy = 0.;
  theta = 0.;

  // Comparing points 1 and 2

  dx12 = refpt[1].xtrace - refpt[0].xtrace;
  dy12 = refpt[1].ytrace - refpt[0].ytrace;
  di12 = refpt[1].iline - refpt[0].iline;
  dj12 = refpt[1].xline - refpt[0].xline;

  if (dj12 == 0)
  {
    theta = atan2(dy12, dx12) * 180. / GV_PI;
    dx = sqrt((dx12 * dx12 + dy12 * dy12) / (di12 * di12));
    dy = 0.;
    cost = dx12 / (dx * di12);
    sint = dy12 / (dx * di12);
  }
  else if (di12 == 0)
  {
    theta = atan2(dy12, dx12) * 180. / GV_PI - 90.;
    dy = sqrt((dx12 * dx12 + dy12 * dy12) / (dj12 * dj12));
    dx = 0.;
    sint = -dx12 / (dy * dj12);
    cost = dy12 / (dy * dj12);
  }

  // Origin of the grid

  di10 = refpt[0].iline - refstats.ilming;
  dj10 = refpt[0].xline - refstats.xlming;
  x0 = refpt[0].xtrace - di10 * dx * cost + dj10 * dy * sint;
  y0 = refpt[0].ytrace - di10 * dx * sint - dj10 * dy * cost;
  nz = static_cast<int>((refstats.zmaxl - refstats.zminl) / dz);

  // Fill the Grid structure

  def_grid = Grid(3);
  def_grid.setNX(0, refstats.ilmaxg - refstats.ilming + 1);
  def_grid.setNX(1, refstats.xlmaxg - refstats.xlming + 1);
  def_grid.setNX(2, nz);
  def_grid.setX0(0, x0);
  def_grid.setX0(1, y0);
  def_grid.setX0(2, refstats.zminl);
  def_grid.setDX(0, MAX(dy / 2., dx));
  def_grid.setDX(1, MAX(dx / 2., dy));
  def_grid.setDX(2, dz);
  def_grid.setRotationByAngle(theta);

  // Check that reference traces are well honored

  st_verify_refpt(refpt, refstats, 2, x0, y0, dx, dy, sint, cost);
}

/****************************************************************************/
/*!
 ** Store an additional reference point
 **
 ** \return Number of RefPt already stored
 **
 ** \param[in]  nbrefpt   Number of RefPt currently allocated
 ** \param[in]  refpt     Array of RefPt structure pointers
 ** \param[in]  iline     Inline number
 ** \param[in]  xline     Crossline number
 ** \param[in]  xtrace    Coordinate of the trace along X
 ** \param[in]  ytrace    Coordinate of the trace along y
 **
 *****************************************************************************/
static int st_store_refpt(int nbrefpt,
                          RefPt refpt[3],
                          int iline,
                          int xline,
                          double xtrace,
                          double ytrace)
{
  RefPt *ref0, *ref1, *ref2;

  // Check if the newly defined trace should be referenced or not

  if (nbrefpt >= 3)
    return (nbrefpt);
  else if (nbrefpt == 1)
  {
    ref1 = &refpt[0];
    if (ref1->iline == iline && ref1->xline == xline) return (nbrefpt);
  }
  else if (nbrefpt >= 2)
  {
    ref1 = &refpt[0];
    ref2 = &refpt[1];
    if (ref1->iline == ref2->iline && ref1->iline == iline) return (nbrefpt);
    if (ref1->xline == ref2->xline && ref1->xline == xline) return (nbrefpt);
  }

  // Store the new Reference point

  ref0 = &refpt[nbrefpt];
  ref0->iline = iline;
  ref0->xline = xline;
  ref0->xtrace = xtrace;
  ref0->ytrace = ytrace;
  nbrefpt++;
  return (nbrefpt);
}

/****************************************************************************/
/*!
 ** Print the characteristics of the resulting grid
 **
 ** \param[in]  def_grid  Pointer to the Grid structure
 **
 *****************************************************************************/
static void st_print_grid(const Grid &def_grid)
{
  message("\n");
  message("- Resulting Grid of SEGY traces\n");
  def_grid.display();
}

/****************************************************************************/
/*!
 ** Complete the information in the Squeeze and Stretch feature
 **
 ** \return 1 if the vector only contains TEST values
 **
 ** \param[in]     ntab   Dimension of the vector
 ** \param[in,out] tab    Vector
 **
 *****************************************************************************/
static int st_complete_squeeze_and_stretch(int ntab, VectorDouble &tab)
{
  double v1;

  // Eliminate the NA values downwards

  v1 = TEST;
  for (int i = 0; i < ntab; i++)
  {
    if (!FFFF(tab[i]))
      v1 = tab[i];
    else
      tab[i] = v1;
  }
  if (FFFF(v1)) return 1;

  // Eliminate the NA values upwards

  v1 = TEST;
  for (int j = 0; j < ntab; j++)
  {
    int i = ntab - 1 - j;
    if (!FFFF(tab[i]))
      v1 = tab[i];
    else
      tab[i] = v1;
  }
  return 0;
}

/****************************************************************************/
/*!
 ** Get the vertical limits in indices
 **
 *****************************************************************************/
static bool st_vertical_limits(double z0,
                               double delta,
                               double cztop,
                               double czbot,
                               double cote,
                               int option,
                               int nz,
                               int *iz1_ret,
                               int *iz2_ret)
{
  int iz1, iz2;

  iz1 = iz2 = -1;
  switch (option)
  {
    case 0: // No flattening
      iz1 = iz2 = static_cast<int>((cote - z0) / delta);
      break;

    case -1: // Flattening from Bottom surface
      if (FFFF(czbot)) return true;
      iz1 = iz2 = static_cast<int>((cote - czbot) / delta);
      break;

    case 1: // Flattening from Top surface
      if (FFFF(cztop)) return true;
      iz1 = iz2 = static_cast<int>((nz - 1) - (cztop - cote) / delta);
      break;

    case 2: // Vertical averaging (iz1 & iz2 are not used)
      iz1 = iz2 = 0;
      break;

    case -2: // Squeeze and stretch
      if (FFFF(czbot)) return true;
      if (FFFF(cztop)) return true;
      double dz = (cztop - czbot) / (double) (nz - 1);
      iz1 = static_cast<int>(floor((cote - czbot) / dz));
      iz2 = static_cast<int>(ceil((cote - czbot + delta) / dz));
      break;
  }

  // Calculate the vertical index

  iz1 = MAX(0, MIN(iz1, nz-1));
  iz2 = MAX(0, MIN(iz2, nz-1));

  // Returning arguments

  *iz1_ret = iz1;
  *iz2_ret = iz2;
  return (iz1 > iz2);
}

/****************************************************************************/
/*!
 ** Calculate the average
 **
 *****************************************************************************/
static double st_get_average(int nz, const VectorDouble &writes)
{
  double dmean = 0.;
  int nmean = 0;
  for (int iz = 0; iz < nz; iz++)
  {
    if (FFFF(writes[iz])) continue;
    dmean += writes[iz];
    nmean++;
  }
  dmean = (nmean > 0) ? dmean / (double) nmean :
                        TEST;
  return dmean;
}

/****************************************************************************/
/*!
 ** Processing a Trace
 **
 *****************************************************************************/
static int st_load_trace(int nPerTrace,
                         int nz,
                         int option,
                         double z0,
                         double delta,
                         double czbot,
                         double cztop,
                         const VectorDouble &cotes,
                         const VectorDouble &values,
                         VectorDouble &writes,
                         int *nbvalues,
                         RefStats &refstats)
{
  int iz1, iz2;
  *nbvalues = 0;
  for (int i = 0; i < nPerTrace; i++)
    writes[i] = TEST;

  refstats.zminl = TEST;
  refstats.zmaxl = TEST;
  refstats.vminl = TEST;
  refstats.vmaxl = TEST;

  for (int i = 0; i < nPerTrace; i++)
  {

    // Check that the elevation must be kept
    double cote = cotes[i];
    double dval = values[i];
    if (!FFFF(czbot) && cote < czbot) continue;
    if (!FFFF(cztop) && cote > cztop) continue;
    if (!FFFF(dval))
    {
      (*nbvalues)++;
      refstats.nbvalue_in++;
      if (!FFFF(refstats.modif_high) && dval > refstats.modif_high)
        dval = refstats.modif_high;
      if (!FFFF(refstats.modif_low) && dval < refstats.modif_low)
        dval = refstats.modif_low;
      if (!FFFF(refstats.modif_scale)) dval /= refstats.modif_scale;
    }

    // Dispatch

    if (st_vertical_limits(z0, delta, cztop, czbot, cote, option, nz, &iz1,
                           &iz2)) continue;

    // Update statistics

    if (cote < refstats.zming || FFFF(refstats.zming)) refstats.zming = cote;
    if (cote > refstats.zmaxg || FFFF(refstats.zmaxg)) refstats.zmaxg = cote;
    if (cote < refstats.zminl || FFFF(refstats.zminl)) refstats.zminl = cote;
    if (cote > refstats.zmaxl || FFFF(refstats.zmaxl)) refstats.zmaxl = cote;

    if (!FFFF(dval))
    {
      if (dval < refstats.vming || FFFF(refstats.vming)) refstats.vming = dval;
      if (dval > refstats.vmaxg || FFFF(refstats.vmaxg)) refstats.vmaxg = dval;
      if (dval < refstats.vminl || FFFF(refstats.vminl)) refstats.vminl = dval;
      if (dval > refstats.vmaxl || FFFF(refstats.vmaxl)) refstats.vmaxl = dval;
    }

    // Writing in temporary array

    for (int iz = iz1; iz <= iz2; iz++)
      writes[iz] = dval;
  }

  if ((*nbvalues) > 0)
  {
    refstats.nbvalue_in += (*nbvalues);
    refstats.nbtrace_def++;
    refstats.thicl = refstats.zmaxl - refstats.zminl;
    if (refstats.thicl > refstats.thicg || FFFF(refstats.thicg))
      refstats.thicg = refstats.thicl;
  }

  return ((*nbvalues) <= 0);
}

/****************************************************************************/
/*!
 ** Printout
 **
 *****************************************************************************/
static void st_print_results(int nPerTrace,
                             int flag_surf,
                             double delta,
                             RefStats &refstats)
{
  mestitle(1, "Extracting information from the SEGY file");
  if (flag_surf)
    message(
        "(Statistics are calculated taking Surface information into account)\n");
  message("- Number of samples per Trace    = %d \n", nPerTrace);
  message("- Interval between samples       = %lf\n", delta);
  message("- Minimum Inline number          = %d \n", refstats.ilming);
  message("- Maximum Inline number          = %d \n", refstats.ilmaxg);
  message("- Minimum Xline number           = %d \n", refstats.xlming);
  message("- Maximum Xline number           = %d \n", refstats.xlmaxg);
  message("- Minimum coordinate along X     = %12.4lf\n", refstats.xming);
  message("- Maximum coordinate along X     = %12.4lf\n", refstats.xmaxg);
  message("- Minimum coordinate along Y     = %12.4lf\n", refstats.yming);
  message("- Maximum coordinate along Y     = %12.4lf\n", refstats.ymaxg);
  message("- Minimum coordinate along Z     = %12.4lf\n", refstats.zming);
  message("- Maximum coordinate along Z     = %12.4lf\n", refstats.zmaxg);
  message("- Maximum thickness along Z      = %12.4lf\n", refstats.thicg);
  message("\n");

  message("- Number of traces read\n");
  message("  . Total                        = %d\n", refstats.nbtrace);
  message("  . Within Unit                  = %d\n", refstats.nbtrace_in);
  message("  . With information             = %d\n", refstats.nbtrace_def);
  message("- Number of valid values read    = %d\n", refstats.nbvalue_in);
  if (!FFFF(refstats.modif_high))
    message("  . Upper Truncation bound       = %lf\n", refstats.modif_high);
  if (!FFFF(refstats.modif_low))
    message("  . Lower Truncation bound       = %lf\n", refstats.modif_low);
  if (!FFFF(refstats.modif_scale))
    message("  . Scaling value                = %lf\n", refstats.modif_scale);
  message("  . Minimum value                = %lf\n", refstats.vming);
  message("  . Maximum value                = %lf\n", refstats.vmaxg);
}

/****************************************************************************/
/*!
 ** Get Trace characteristics
 **
 *****************************************************************************/
void st_get_trace_params(traceHead *Theader,
                         int *iline,
                         int *xline,
                         double *delta,
                         double *xtrace,
                         double *ytrace)
{
  int scacsv = st_to_s(Theader->SCALE_COOR);
  *delta = (double) st_to_s(Theader->SAMPLE_INTRVL) / 1000;
  *xtrace = st_scaling(st_to_f(Theader->ENS_COOR_X), scacsv);
  *ytrace = st_scaling(st_to_f(Theader->ENS_COOR_Y), scacsv);
  *iline = st_to_i(Theader->INLINE);
  *xline = st_to_i(Theader->CROSS);
}

/****************************************************************************/
/*!
 ** Reject trace due to trace boundary specifications
 **
 *****************************************************************************/
static int st_reject_trace(int iline,
                           int xline,
                           int iline_min,
                           int iline_max,
                           int xline_min,
                           int xline_max)
{
  if (!IFFFF(iline_min) && iline < iline_min) return (1);
  if (!IFFFF(iline_max) && iline > iline_max) return (1);
  if (!IFFFF(xline_min) && xline < xline_min) return (1);
  if (!IFFFF(xline_max) && xline > xline_max) return (1);
  return (0);
}

/****************************************************************************/
/*!
 ** Update statistics
 **
 *****************************************************************************/
static void st_refstats_update(int iline,
                               int xline,
                               double xtrace,
                               double ytrace,
                               RefStats &refstats)
{
  if (xtrace < refstats.xming || FFFF(refstats.xming)) refstats.xming = xtrace;
  if (xtrace > refstats.xmaxg || FFFF(refstats.xmaxg)) refstats.xmaxg = xtrace;
  if (ytrace < refstats.yming || FFFF(refstats.yming)) refstats.yming = ytrace;
  if (ytrace > refstats.ymaxg || FFFF(refstats.ymaxg)) refstats.ymaxg = ytrace;
  if (iline < refstats.ilming || IFFFF(refstats.ilming))
    refstats.ilming = iline;
  if (iline > refstats.ilmaxg || IFFFF(refstats.ilmaxg))
    refstats.ilmaxg = iline;
  if (xline < refstats.xlming || IFFFF(refstats.xlming))
    refstats.xlming = xline;
  if (xline > refstats.xlmaxg || IFFFF(refstats.xlmaxg))
    refstats.xlmaxg = xline;
  refstats.nbtrace_in++;
}

/****************************************************************************/
/*!
 ** Initialize statistics
 **
 *****************************************************************************/
static void st_refstats_init(RefStats &refstats,
                             double modif_high,
                             double modif_low,
                             double modif_scale)
{
  refstats.nbtrace = 0;
  refstats.nbtrace_in = 0;
  refstats.nbtrace_def = 0;
  refstats.nbvalue_in = 0;

  refstats.ilming = ITEST;
  refstats.ilmaxg = ITEST;
  refstats.xlming = ITEST;
  refstats.xlmaxg = ITEST;

  refstats.zminl = TEST;
  refstats.zmaxl = TEST;
  refstats.vminl = TEST;
  refstats.vmaxl = TEST;
  refstats.thicl = TEST;

  refstats.xming = TEST;
  refstats.xmaxg = TEST;
  refstats.yming = TEST;
  refstats.ymaxg = TEST;
  refstats.zming = TEST;
  refstats.zmaxg = TEST;
  refstats.vming = TEST;
  refstats.vmaxg = TEST;
  refstats.thicg = TEST;

  refstats.modif_low = modif_low;
  refstats.modif_high = modif_high;
  refstats.modif_scale = modif_scale;
}
;

/****************************************************************************/
/*!
 ** Read the contents of a SEGY file
 **
 ** \returns A RetArg structure which contains:
 ** \returns - a vector of trace vectors
 ** \returns - a vector of trace descriptors
 ** \returns The Descriptor for each trace contains:
 ** \returns 0: Absolute rank for the trace number
 ** \returns 1: Cross-Line number
 ** \returns 2: In-Line number
 ** \returns 3: Coordinate along X
 ** \returns 4: Coordinate along Y
 ** \returns 5: Minimum Elevation
 ** \returns 6: Maximum Elevation
 ** \returns 7: Minimum Value
 ** \returns 8: Maximum value
 **
 ** \param[in]  filesegy    Name of the SEGY file
 ** \param[in]  surf2D      Db containing the top, Bottom and Reference surfaces
 **                         This file is optional
 ** \param[in]  name_top    Rank of variable containing the Top Surface (or 0)
 ** \param[in]  name_bot    Rank of variable containing the Bottom Surface (or 0)
 ** \param[in]  thickmin    Minimum thickness (or 0)
 ** \param[in]  option      Flattening option:
 **                          0 no flattening;
 **                          1 flattening from top;
 **                         -1 flattening from bottom
 **                         -2 squeeze and stretch option
 **                          2 averaging from 3-D to 2-D
 ** \param[in]  verbOption  Verbose option
 ** \param[in]  iline_min   Minimum Inline number included (if defined)
 ** \param[in]  iline_max   Maximum Inline number included (if defined)
 ** \param[in]  xline_min   Minimum Xline number included (if defined)
 ** \param[in]  xline_max   Maximum Xline number included (if defined)
 ** \param[in]  nz_ss       Number of layers for different options (see details)
 ** \param[in]  modif_high  Upper truncation (when defined)
 ** \param[in]  modif_low   Lower truncation (when defined)
 ** \param[in]  modif_scale Scaling value (when defined)
 **
 ** \details In the case of Squeeze and Stretch (S&S), the number of layers
 ** \details is meaningless. It is fixed by the user.
 **
 *****************************************************************************/
SegYArg segy_array(const char *filesegy,
                                   Db *surf2D,
                                   const String &name_top,
                                   const String &name_bot,
                                   double thickmin,
                                   int option,
                                   int verbOption,
                                   int iline_min,
                                   int iline_max,
                                   int xline_min,
                                   int xline_max,
                                   int nz_ss,
                                   double modif_high,
                                   double modif_low,
                                   double modif_scale)
{
  traceHead traceHead_;
  FILE *file;
  double xtrace, ytrace, z0, cztop, czbot, delta;
  int nPerTrace, code, nz, flag_surf, flag_top, flag_bot, iatt_bot, iatt_top;
  int nbrefpt, iline, xline, nbvalues;
  RefPt refpt[3];
  RefStats refstats;
  VectorDouble values, cotes, writes;
  SegYArg segyarg;

  // Initializations

  code = 1;
  file = nullptr;
  nbrefpt = nz = 0;
  delta = z0 = 0.;
  st_refstats_init(refstats, modif_high, modif_low, modif_scale);

  // Preliminary checks

  flag_surf = (surf2D != nullptr);
  flag_top = flag_surf && (option == 1 || option == -2);
  flag_bot = flag_surf && (option == -1 || option == -2);
  if (st_surface_identify(verbOption, surf2D, name_bot, flag_bot, &iatt_bot,
                          name_top, flag_top, &iatt_top)) return segyarg;

  // Open Input SEGY file

  if ((file = gslFopen(filesegy, "r")) == NULL)
  {
    messerr("ERROR:  cannot find input file %s", filesegy);
    return segyarg;
  }

  // GET HEADER INFO

  if (st_readFileHeader(file, verbOption, &nPerTrace, &delta)) return segyarg;
  nz = nPerTrace;
  z0 = -delta * nPerTrace;
  if (option == -2 && !IFFFF(nz_ss)) nz = nz_ss;

  // Loop on the traces

  values.resize(nPerTrace);
  cotes.resize(nPerTrace);
  writes.resize(nPerTrace);
  while (1)
  {

    // Read the Trace Header

    if (fread(&traceHead_, 1, sizeof(traceHead_), file) == 0) break;
    refstats.nbtrace++;
    if (verbOption >= 3) st_print_traceHead(&traceHead_, refstats.nbtrace);
    st_get_trace_params(&traceHead_, &iline, &xline, &delta, &xtrace, &ytrace);

    // Read the Trace contents

    if (st_read_trace(file, code, refstats.nbtrace, nPerTrace, delta,
                      verbOption, values, cotes)) continue;

    // Reject trace due to trace boundary specifications

    if (st_reject_trace(iline, xline, iline_min, iline_max, xline_min,
                        xline_max)) continue;

    // Store the reference point

    nbrefpt = st_store_refpt(nbrefpt, refpt, iline, xline, xtrace, ytrace);

    // Compare to the 2-D Surface information
    // If 'surfaces' are not provided, bounds are set to TEST.

    if (st_get_cuts(surf2D, iatt_top, iatt_bot, xtrace, ytrace, thickmin,
                    &cztop, &czbot)) continue;

    // Update statistics

    st_refstats_update(iline, xline, xtrace, ytrace, refstats);

    // Looking at samples along the trace

    if (st_load_trace(nPerTrace, nz, option, z0, delta, czbot, cztop, cotes,
                      values, writes, &nbvalues, refstats)) continue;

    // Particular case of squeeze and stretch

    if (option == -2)
    {
      if (st_complete_squeeze_and_stretch(nz, writes)) continue;
    }

    // Store the results

    segyarg.ndescr = SEGY_COUNT;
    VectorDouble descr(SEGY_COUNT);
    descr[SEGY_NUM] = refstats.nbtrace;
    descr[SEGY_ILINE] = iline;
    descr[SEGY_XLINE] = xline;
    descr[SEGY_XTRACE] = xtrace;
    descr[SEGY_YTRACE] = ytrace;
    descr[SEGY_ZMIN] = refstats.zminl;
    descr[SEGY_ZMAX] = refstats.zmaxl;
    descr[SEGY_VMIN] = refstats.vminl;
    descr[SEGY_VMAX] = refstats.vmaxl;
    descr[SEGY_THICK] = refstats.thicl;
    descr[SEGY_NB] = nbvalues;

    segyarg.descr.push_back(descr);

    segyarg.npertrace = nz;
    VectorDouble local(nz);
    for (int iz = 0; iz < nz; iz++)
      local[iz] = writes[iz];
    segyarg.tab.push_back(local);
  }

  // Optional printout

  if (verbOption)
  {
    st_print_results(nPerTrace, flag_surf, delta, refstats);
  }

  if (file != NULL) fclose(file);

  return segyarg;
}

/****************************************************************************/
/*!
 ** Read the contents of a SEGY file and returns the structure SegyRead
 ** which captures the main characteristics of the SEGY grid
 **
 ** \param[in]  filesegy    Name of the SEGY file
 ** \param[in]  surf2D      Db containing the top, Bottom and Reference surfaces
 **                         This file is optional
 ** \param[in]  name_top    Rank of variable containing the Top Surface (or 0)
 ** \param[in]  name_bot    Rank of variable containing the Bottom Surface (or 0)
 ** \param[in]  thickmin    Minimum thickness (or 0)
 ** \param[in]  option      Flattening option:
 **                          0 no flattening;
 **                          1 flattening from top;
 **                         -1 flattening from bottom
 **                         -2 squeeze and stretch option
 **                          2 averaging from 3-D to 2-D
 ** \param[in]  verbOption  Verbose option
 ** \param[in]  iline_min   Minimum Inline number included (if defined)
 ** \param[in]  iline_max   Maximum Inline number included (if defined)
 ** \param[in]  xline_min   Minimum Xline number included (if defined)
 ** \param[in]  xline_max   Maximum Xline number included (if defined)
 ** \param[in]  nz_ss       Number of layers for different options (see details)
 ** \param[in]  modif_high  Upper truncation (when defined)
 ** \param[in]  modif_low   Lower truncation (when defined)
 ** \param[in]  modif_scale Scaling value (when defined)
 **
 ** \details: In the case of Squeeze and Stretch (S&S), the number of layers
 ** \details: is meaningless. It is fixed by the user, unless defined
 ** \details: by the output grid (if flag_store == 1)
 **
 *****************************************************************************/
Grid segy_summary(const char *filesegy,
                                   Db *surf2D,
                                   const String &name_top,
                                   const String &name_bot,
                                   double thickmin,
                                   int option,
                                   int verbOption,
                                   int iline_min,
                                   int iline_max,
                                   int xline_min,
                                   int xline_max,
                                   int nz_ss,
                                   double modif_high,
                                   double modif_low,
                                   double modif_scale)
{
  traceHead traceHead_;
  FILE    *file;
  double   xtrace,ytrace,z0,cztop,czbot,delta;
  int      iline,xline,nbvalues;
  int      nPerTrace,code,nz,flag_surf,flag_top,flag_bot,iatt_top,iatt_bot,nbrefpt;
  RefPt    refpt[3];
  RefStats refstats;
  Grid def_grid;
  VectorDouble values, cotes, writes;

  // Initializations

  code = 1;
  file = nullptr;
  nbrefpt = nz = 0;
  delta = z0 = 0.;
  st_refstats_init(refstats, modif_high, modif_low, modif_scale);

  // Preliminary checks

  flag_surf = (surf2D != nullptr);
  flag_top = flag_surf && (option == 1 || option == -2);
  flag_bot = flag_surf && (option == -1 || option == -2);
  if (st_surface_identify(verbOption, surf2D, name_bot, flag_bot, &iatt_bot,
                          name_top, flag_top, &iatt_top)) return def_grid;

  // Open Input SEGY file

  if ((file = gslFopen(filesegy, "r")) == NULL)
  {
    messerr("ERROR:  cannot find input file %s", filesegy);
    return def_grid;
  }

  // GET HEADER INFO

  if (st_readFileHeader(file, verbOption, &nPerTrace, &delta)) return def_grid;
  nz = nPerTrace;
  z0 = -delta * nPerTrace;
  if (option == -2 && !IFFFF(nz_ss)) nz = nz_ss;

  // Working arrays

  values.resize(nPerTrace);
  cotes.resize(nPerTrace);
  writes.resize(nPerTrace);

  // Loop on the traces

  while (1)
  {

    // Read the Trace Header

    if (fread(&traceHead_, 1, sizeof(traceHead_), file) == 0) break;
    refstats.nbtrace++;
    if (verbOption >= 3) st_print_traceHead(&traceHead_, refstats.nbtrace);
    st_get_trace_params(&traceHead_, &iline, &xline, &delta, &xtrace, &ytrace);

    // Read the Trace contents

    if (st_read_trace(file, code, refstats.nbtrace, nPerTrace, delta,
                      verbOption, values, cotes)) continue;

    // Reject trace due to trace boundary specifications

    if (st_reject_trace(iline, xline, iline_min, iline_max, xline_min,
                        xline_max)) continue;

    // Store the reference point

    nbrefpt = st_store_refpt(nbrefpt, refpt, iline, xline, xtrace, ytrace);

    // Compare to the 2-D Surface information
    // If 'surfaces' is not provided, bounds are set to TEST.

    if (st_get_cuts(surf2D, iatt_top, iatt_bot, xtrace, ytrace, thickmin,
                    &cztop, &czbot)) continue;

    // Update statistics

    st_refstats_update(iline, xline, xtrace, ytrace, refstats);

    // Looking at samples along the trace

    if (st_load_trace(nPerTrace, nz, option, z0, delta, czbot, cztop, cotes,
                      values, writes, &nbvalues, refstats)) continue;
  }

  // Extract the grid characteristics

  if (nbrefpt == 3)
    st_grid_from_3refpt(refpt, refstats, delta, def_grid);
  else if (nbrefpt == 2) st_grid_from_2refpt(refpt, refstats, delta, def_grid);
  if (option == 2)
  {
    def_grid.setX0(2, 0.);
    def_grid.setDX(2, 1.);
    def_grid.setNX(2, 1);
  }

  // Optional printout

  if (verbOption)
  {
    st_print_results(nPerTrace, flag_surf, delta, refstats);
    st_print_grid(def_grid);
  }

  if (file != NULL) fclose(file);
  return def_grid;
}

/****************************************************************************/
/*!
 ** Read the contents of a SEGY file
 **
 ** \return Error return code
 **
 ** \param[in]  filesegy    Name of the SEGY file
 ** \param[in]  grid3D      Db containing the resulting 3-D grid
 ** \param[in]  surf2D      Db containing the top, Bottom and Reference surfaces
 **                         This file is optional
 ** \param[in]  name_top    Rank of variable containing the Top Surface (or 0)
 ** \param[in]  name_bot    Rank of variable containing the Bottom Surface (or 0)
 ** \param[in]  thickmin    Minimum thickness (or 0)
 ** \param[in]  option      Flattening option:
 **                          0 no flattening;
 **                          1 flattening from top;
 **                         -1 flattening from bottom
 **                         -2 squeeze and stretch option
 **                          2 averaging from 3-D to 2-D
 ** \param[in]  verbOption  Verbose option
 ** \param[in]  iline_min   Minimum Inline number included (if defined)
 ** \param[in]  iline_max   Maximum Inline number included (if defined)
 ** \param[in]  xline_min   Minimum Xline number included (if defined)
 ** \param[in]  xline_max   Maximum Xline number included (if defined)
 ** \param[in]  nz_ss       Number of layers for different options (see details)
 ** \param[in]  modif_high  Upper truncation (when defined)
 ** \param[in]  modif_low   Lower truncation (when defined)
 ** \param[in]  modif_scale Scaling value (when defined)
 ** \param[in]  namconv     Naming convention
 **
 ** \details: In the case of Squeeze and Stretch (S&S), the number of layers
 ** \details: is meaningless. It is fixed by the user, unless defined
 ** \details: by the output grid (if flag_store == 1)
 **
 ** \remarks For filling the 3-D grid
 ** \remarks - 2-D characteristics of the grid are taken into account
 ** \remarks - no attention is paid to the vertical mesh value.
 **
 *****************************************************************************/
int db_segy(const char *filesegy,
                            Db *grid3D,
                            Db *surf2D,
                            const String &name_top,
                            const String &name_bot,
                            double thickmin,
                            int option,
                            int verbOption,
                            int iline_min,
                            int iline_max,
                            int xline_min,
                            int xline_max,
                            int nz_ss,
                            double modif_high,
                            double modif_low,
                            double modif_scale,
                            const NamingConvention& namconv)
{
  traceHead traceHead_;
  FILE *file;
  double xtrace, ytrace, z0, coor[3], cztop, czbot, delta;
  int iline, xline, nbvalues, iatt;
  int indg[3], nPerTrace, code, nz, rank;
  int flag_surf, flag_top, flag_bot, iatt_top, iatt_bot;
  int nbrefpt;
  RefPt refpt[3];
  RefStats refstats;
  VectorDouble values, cotes, writes;

  // Initializations

  code = 1;
  file = nullptr;
  nbrefpt = nz = 0;
  delta = z0 = 0.;
  for (int i = 0; i < 3; i++)
    indg[i] = 0;
  iatt = -1;
  st_refstats_init(refstats, modif_high, modif_low, modif_scale);

  // Preliminary checks

  flag_surf = (surf2D != nullptr);
  flag_top = flag_surf && (option == 1 || option == -2);
  flag_bot = flag_surf && (option == -1 || option == -2);
  if (st_surface_identify(verbOption, surf2D, name_bot, flag_bot, &iatt_bot,
                          name_top, flag_top, &iatt_top)) return 1;
  z0 = grid3D->getX0(2);
  nz = grid3D->getNX(2);

  // Open Input SEGY file

  file = gslFopen(filesegy, "r");
  if (file == NULL)
  {
    messerr("ERROR:  cannot find input file %s", filesegy);
    return 1;
  }

  // GET HEADER INFO

  if (st_readFileHeader(file, verbOption, &nPerTrace, &delta)) return 1;

  // Allocate the new vector in the output file

  iatt = grid3D->addFieldsByConstant(1, TEST);
  if (iatt < 0) return 1;

  // Working arrays

  values.resize(nPerTrace);
  cotes.resize(nPerTrace);
  writes.resize(nPerTrace);

  // Loop on the traces

  while (1)
  {

    // Read the Trace Header

    if (fread(&traceHead_, 1, sizeof(traceHead_), file) == 0) break;
    refstats.nbtrace++;
    if (verbOption >= 3) st_print_traceHead(&traceHead_, refstats.nbtrace);
    st_get_trace_params(&traceHead_, &iline, &xline, &delta, &xtrace, &ytrace);

    // Read the Trace contents

    if (st_read_trace(file, code, refstats.nbtrace, nPerTrace, delta,
                      verbOption, values, cotes)) continue;

    // Reject trace due to trace boundary specifications

    if (st_reject_trace(iline, xline, iline_min, iline_max, xline_min,
                        xline_max)) continue;

    // Store the reference point

    nbrefpt = st_store_refpt(nbrefpt, refpt, iline, xline, xtrace, ytrace);

    // Find the position of the trace within the 3-D grid

    coor[0] = xtrace;
    coor[1] = ytrace;
    coor[2] = z0 + delta;
    if (option == 2) coor[2] = 0.;
    if (point_to_grid(grid3D, coor, 0, indg) != 0) continue;

    // Compare to the 2-D Surface information
    // If 'surfaces' is not provided, bounds are set to TEST.

    if (st_get_cuts(surf2D, iatt_top, iatt_bot, xtrace, ytrace, thickmin,
                    &cztop, &czbot)) continue;

    // Update statistics

    st_refstats_update(iline, xline, xtrace, ytrace, refstats);

    // Looking at samples along the trace

    if (st_load_trace(nPerTrace, nz, option, z0, delta, czbot, cztop, cotes,
                      values, writes, &nbvalues, refstats)) continue;

    // Particular case of squeeze and stretch

    if (option == -2)
    {
      if (st_complete_squeeze_and_stretch(nz, writes)) continue;
    }

    // Storing

    if (option == 2)
    {
      double dmean = st_get_average(nz, writes);
      indg[2] = 0;
      rank = db_index_grid_to_sample(grid3D, indg);
      if (rank >= 0) grid3D->setArray(rank, iatt, dmean);
    }
    else
    {
      for (int iz = 0; iz < nz; iz++)
      {
        indg[2] = iz;
        rank = db_index_grid_to_sample(grid3D, indg);
        if (rank < 0) continue;
        grid3D->setArray(rank, iatt, writes[iz]);
      }
    }
  }

  // Optional printout

  if (verbOption)
  {
    st_print_results(nPerTrace, flag_surf, delta, refstats);
  }

  // Assign Name and Locator to the newly created variable

  namconv.setNamesAndLocators(grid3D, iatt, String());

  if (file != NULL) fclose(file);
  return 0;
}

