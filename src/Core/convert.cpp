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
#include "Basic/AException.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/File.hpp"
#include "vtk.h"
#include "geoslib_e.h"
#include "geoslib_old_f.h"

#include <sstream>

/*! \cond */
#define COLOR_MASK   -1
#define COLOR_FFFF   -2
#define COLOR_LOWER  -3
#define COLOR_UPPER  -4
#define N_SAMPLE(nx,nsample) ((int) ((nx-1) / nsample) + 1)
#define ZYCOR_NULL_CH "  0.1000000E+31"
#define F2G(ix,iy,iz)   (tab[(ix) + nx[0] * ((iy) + nx[1] * (iz))]) 
#define BF_TYPE 0x4D42             /* "MB" */

struct CSV_Encoding {
  FILE* file;                // Stream used for writing into CSV file
  int nitem;                 // Number of items per line
  int current;               // Rank of the current item
  int nlines;                // Number of lines printed
  bool flag_integer;         // true for Integer encoding
  const char *char_sep;      // Separator between consecutive fields
  const char *na_string;     // Substitute for NA
};

static CSV_Encoding *CSV_ENCODE = NULL;

/*! \endcond */

/****************************************************************************/
/*!
**   Color rank for the sample of the Db
**
** \return  Color rank
** \return  COLOR_MASK  : The sample is masked off
** \return  COLOR_FFFF  : The value if FFFF
** \return  COLOR_LOWER : The value if below vmin and must be acknowledged
** \return  COLOR_UPPER : The value if above vmax and must be acknowledged
**
** \param[in]  db         Db structure
** \param[in]  iech       Rank of the sample
** \param[in]  icol       Rank of the attribute
** \param[in]  ncolor     Number of colors
** \param[in]  flag_low   Returns the first color for value below vmin if 0
**                        or returns -3 otherwise
** \param[in]  flag_high  Returns the last color for value above vmax if 0
**                        or returns -4 otherwise
** \param[in]  vmin       Minimum value to be represented
** \param[in]  vmax       Maximum value to be represented
**
*****************************************************************************/
static int st_color_rank(Db    *db,
                         int    iech,
                         int    icol,
                         int    ncolor,
                         int    flag_low,
                         int    flag_high,
                         double vmin,
                         double vmax)
{
  double value;
  int    ival;

  /* Check if the sample is masked off */
  if (! db->getSelection(iech)) return(COLOR_MASK);

  /* Read the value */
  value = db->getArray(iech,icol);

  /* Check if the value is defined */
  if (FFFF(value)) return(COLOR_FFFF);

  /* Find the color */
  ival = (int) (ncolor * (value - vmin) / (vmax - vmin));

  /* Value lower then vmin */
  if (ival < 0)
  {
    if (flag_low) 
      return(COLOR_LOWER);
    else
      return(ival);
  }

  /* Value larger than vmax */
  if (ival >= ncolor)
  {
    if (flag_high)
      return(COLOR_UPPER);
    else
      return(ncolor-1);
  }

  /* Return the rank of the color */
  return(ival);
}

/****************************************************************************/
/*!
** Convert a color rank into the Red, Green, Blue color decompoisition
**
** \param[in]  rank        Rank of the color
** \param[in]  flag_color_scale 1 if the color scale must be used
**                              0 use the grey scale instead
** \param[in]  red         Array of Red intensity for color scale
** \param[in]  green       Array of Green intensity for color scale
** \param[in]  blue        Array of Blue intensity for color scale
** \param[in]  mask_red    Red intensity for masked value
** \param[in]  mask_green  Green intensity for masked value
** \param[in]  mask_blue   Blue intensity for masked value
** \param[in]  ffff_red    Red intensity for FFFF value
** \param[in]  ffff_green  Green intensity for FFFF value
** \param[in]  ffff_blue   Blue intensity for FFFF value
** \param[in]  low_red     Red intensity for lower value
** \param[in]  low_green   Green intensity for lower value
** \param[in]  low_blue    Blue intensity for lower value
** \param[in]  high_red    Red intensity for higher value
** \param[in]  high_green  Green intensity for higher value
** \param[in]  high_blue   Blue intensity for higher value
**  
** \param[out] ired        Value for the red beam
** \param[out] igreen      Value for the green beam
** \param[out] iblue       Value for the blue beam
**
*****************************************************************************/
static void st_color_in_rgb(int     rank,
                            int     flag_color_scale,
                            int    *red,
                            int    *green,
                            int    *blue,
                            int     mask_red,
                            int     mask_green,
                            int     mask_blue,
                            int     ffff_red,
                            int     ffff_green,
                            int     ffff_blue,
                            int     low_red,
                            int     low_green,
                            int     low_blue,
                            int     high_red,
                            int     high_green,
                            int     high_blue,
                            unsigned char *ired,
                            unsigned char *igreen,
                            unsigned char *iblue)
{
  switch (rank)
  {
  case COLOR_MASK:
    *ired   = mask_red;
    *igreen = mask_green;
    *iblue  = mask_blue;
    break;
    
  case COLOR_FFFF:
    *ired   = ffff_red;
    *igreen = ffff_green;
    *iblue  = ffff_blue;
    break;
    
  case COLOR_LOWER:
    *ired   = low_red;
    *igreen = low_green;
    *iblue  = low_blue;
    break;
    
  case COLOR_UPPER:
    *ired   = high_red;
    *igreen = high_green;
    *iblue  = high_blue;
    break;
    
  default:
    if (flag_color_scale)
    {
      *ired   = red  [rank];
      *igreen = green[rank];
      *iblue  = blue [rank];
    }
    else
    {
      *ired   = rank;
      *igreen = rank;
      *iblue  = rank;
    }
  }
}

/****************************************************************************/
/*!
**   Encode a line for IFPEN file
**
** \param[in]  file       FILE descriptor
** \param[in]  mode       Type of encoding
** \li                     0 : Comment
** \li                     1 : Integer value
** \li                     2 : Real value
** \param[in]  comment    Comment string (or NULL)
** \param[in]  valint     Integer value
** \param[in]  valrel     Float value
** \param[in]  combis     Second comment (or NULL)
**
*****************************************************************************/
static void st_ifpen_write(FILE       *file,
                           int         mode,
                           const char *comment,
                           int         valint,
                           double      valrel,
                           const char *combis)
{
  char line[100]; 

  /* Initialize the string */

  (void) gslStrcpy(line,"");

  /* Comment */

  if (comment != NULL)
    (void) gslSPrintf(&line[strlen(line)],"%s",comment);

  /* Encoding the value */

  if (mode == 1)
  {
    (void) gslSPrintf(&line[strlen(line)]," %d",valint);
  }
  else if (mode == 2)
  {
    (void) gslSPrintf(&line[strlen(line)]," %lf",valrel);
  }

  /* Secondary comment */

  if (combis != NULL)
    (void) gslSPrintf(&line[strlen(line)]," %s",combis);

  /* Print the line */

  fprintf(file,"%s\n",line);
}

/****************************************************************************/
/*!
**   Decode a line for IFPEN file
**
** \param[in]  file       FILE descriptor
** \param[in]  mode       Type of encoding
** \li                     0 : Comment
** \li                     1 : Integer value
** \li                     2 : Real value
** \param[in]  comment    Comment string (or NULL)
**
** \param[out]  valint     Integer value
** \param[out]  valrel     Float value
**
*****************************************************************************/
static int st_ifpen_read(FILE       *file,
                         int         mode,
                         const char *comment,
                         int        *valint,
                         double     *valrel)
{
  char line[100]; 
  int  start;

  /* Reading the line */

  if (fgets(line,100,file) == NULL) return(1);
  line[strlen(line)-1] = '\0';

  /* Check the comment */

  start = 0;
  if (comment != NULL)
  {
    if (strcmp(line,comment) < 0) return(1);
    start = static_cast<int> (strlen(comment));
  }

  /* Decoding the value */

  if (mode == 1)
  {
    if (sscanf(&line[start],"%d",valint) != 1) return(1);
  }
  else if (mode == 2)
  {
    if (sscanf(&line[start],"%lf",valrel) != 1) return(1);
  }
  return(0);
}

/****************************************************************************/
/*!
**   Print the grid characteristics (verbose option)
**
** \param[in]  ndim      Space dimension
** \param[in]  nx        Array of number of grid nodes
** \param[in]  dx        Array of grid mesh
** \param[in]  x0        Array of grid origin coordinates
**
*****************************************************************************/
static void st_print_verbose(int    ndim,
                             int    nx[3],
                             double dx[3],
                             double x0[3])
{
  if (ndim >= 1) message ("NX = %d\n" ,nx[0]);
  if (ndim >= 2) message ("NY = %d\n" ,nx[1]);
  if (ndim >= 3) message ("NZ = %d\n" ,nx[2]);
  if (ndim >= 1) message ("DX = %lf\n",dx[0]);
  if (ndim >= 2) message ("DY = %lf\n",dx[1]);
  if (ndim >= 3) message ("DZ = %lf\n",dx[2]);
  if (ndim >= 1) message ("X0 = %lf\n",x0[0]);
  if (ndim >= 2) message ("Y0 = %lf\n",x0[1]);
  if (ndim >= 3) message ("Z0 = %lf\n",x0[2]);
}

/****************************************************************************/
/*!
**   Read the Grid characteristics from a ZYCOR file
**
** \return  Error return code
**
** \param[in]  file       Pointer to the ZYCOR file
** \param[in]  verbose    1 for a verbose output; 0 otherwise
**
** \param[out]  nx        Array of number of grid nodes
** \param[out]  x0        Array of grid origin coordinates
** \param[out]  dx        Array of grid mesh
** \param[out]  test      Value for non-information
**
*****************************************************************************/
static int st_grid_read_zycor_header(FILE   *file,
                                     int     verbose,
                                     int    *nx,
                                     double *x0,
                                     double *dx,
                                     double *test)
{
  char   string[100];
  double xf[2],rbid1,rbid2,rbid3;
  int    nval,ibid1,ibid2,ibid3,error;

  /* Initializations */

  error = 1;

  /* Define the delimitors */

  _file_delimitors('!',',','_');

  /* Read the lines */

  if (_record_read(file,"%s",string)) goto label_end;
  if (string[0] != '@')
  {
    messerr("Missing string starting with (@). Instead: '%s'",string);
    goto label_end;
  }
  if (_record_read(file,"%s",string)) goto label_end;
  if (strcmp(string,"GRID"))
  {
    messerr("Missing string (GRID). Instead: '%s'",string);
    goto label_end;
  }
  if (_record_read(file,"%d" ,&nval))  goto label_end;
  if (_record_read(file,"%d" ,&ibid1)) goto label_end;
  if (_record_read(file,"%lg",test))   goto label_end;
  if (_record_read(file,"%s" ,string)) goto label_end;
  if (_record_read(file,"%d" ,&ibid2)) goto label_end;
  if (_record_read(file,"%d" ,&ibid3)) goto label_end;
  if (_record_read(file,"%d" ,&nx[1])) goto label_end;
  if (_record_read(file,"%d" ,&nx[0])) goto label_end;
  if (_record_read(file,"%lf",&x0[0])) goto label_end;
  if (_record_read(file,"%lf",&xf[0])) goto label_end;
  if (_record_read(file,"%lf",&x0[1])) goto label_end;
  if (_record_read(file,"%lf",&xf[1])) goto label_end;
  if (_record_read(file,"%lf",&rbid1)) goto label_end;
  if (_record_read(file,"%lf",&rbid2)) goto label_end;
  if (_record_read(file,"%lf",&rbid3)) goto label_end;

  if (_record_read(file,"%s" ,string)) goto label_end;
  if (strcmp(string,"@"))
  {
    messerr("Missing string (@). Instead: %s",string);
    goto label_end;
  }

  /* Final calculations */

  dx[0] = (xf[0] - x0[0]) / (nx[0] - 1);
  dx[1] = (xf[1] - x0[1]) / (nx[1] - 1);

  /* Verbose optional printout */

  if (verbose) 
  {
    mestitle(0,"ZYCOR File Characteristics :");
    st_print_verbose(2,nx,dx,x0);
  }

  /* Set the error return code */

  error = 0;

label_end:

  /* Reset the delimitors */

  _file_delimitors('#',' ',' ');
  return(error);
}

/****************************************************************************/
/*!
**   Read the Grid characteristics from a ZYCOR file
**
** \return  Error return code
**
** \param[in]  filename  Name of the ZYCOR file
** \param[in]  verbose   Verbose flag
**
** \param[out]  nx        Array of number of grid nodes
** \param[out]  x0        Array of grid origin coordinates
** \param[out]  dx        Array of grid mesh
**
*****************************************************************************/
GEOSLIB_API int db_grid_read_zycor1(const char   *filename,
                                    int     verbose,
                                    int    *nx,
                                    double *x0,
                                    double *dx)
{
  FILE  *file;
  double test;
  int    error;

  /* Initializations */

  error = 1;
  file  = nullptr;

  /* Open the file */

  file = gslFopen(filename,"r");
  if (file == nullptr)
  {
    messerr("Error when opening the ZYCOR file %s for reading",filename);
    goto label_end;
  }

  /* Read the grid characteristics */

  if (st_grid_read_zycor_header(file,verbose,nx,x0,dx,&test)) goto label_end;

  /* Set the error return code */

  error = 0;

label_end:
  if (file != nullptr) fclose(file);
  return(error);
}

/****************************************************************************/
/*!
**   Read an array of values from a ZYCOR file
**
** \return  Error return code
**
** \param[in]  filename  Name of the ZYCOR file
** \param[in]  nx_r      Array of number of grid nodes
** \param[in]  x0_r      Array of grid origin coordinates
** \param[in]  dx_r      Array of grid mesh
**
** \param[out] tab       Array of values
**
*****************************************************************************/
GEOSLIB_API int db_grid_read_zycor2(const char   *filename,
                                    int    *nx_r,
                                    double *x0_r,
                                    double *dx_r,
                                    double *tab)
{
  FILE  *file;
  int    nx[2],nech,lec,error,ix,iy;
  double x0[2],dx[2],test,value;

  /* Initializations */

  error = 1;
  file  = nullptr;

  /* Open the file */

  file = gslFopen(filename,"r");
  if (file == nullptr)
  {
    messerr("Error when opening the ZYCOR file %s for reading",filename);
    goto label_end;
  }

  /* Read the grid characteristics */

  if (st_grid_read_zycor_header(file,0,nx,x0,dx,&test)) goto label_end;

  /* Check the compatibility of the input arguments with the current grid */

  if (nx[0] != nx_r[0] || nx[1] != nx_r[1] || 
      x0[0] != x0_r[0] || x0[1] != x0_r[1] ||
      dx[0] != dx_r[0] || dx[1] != dx_r[1])
  {
    messerr("The grid parameters are incompatible");
    goto label_end;
  }

  /* Read the array of real values */

  lec  = 0;
  nech = nx[0] * nx[1];
  for (ix=0; ix<nx[0]; ix++)
    for (iy=0; iy<nx[1]; iy++)
    {
      if (_record_read(file,"%lf",&value)) break;
      if (value == test) value = TEST;
      
      if (lec >= nech) 
      {
        messerr("The number of values in the file exceeds the grid dimension");
        goto label_end;
      }
      tab[(nx[1] - iy - 1) * nx[0] + ix] = value;
      lec++;
    }

  /* Final check */

  if (lec != nech)
  {
    messerr("The number of decoded values (%d) is not equal to the expected number of grid nodes(%d)",
            lec,nech);
    goto label_end;
  }

  /* Set the error return code */

  error = 0;

label_end:
  record_close();
  if (file != nullptr) fclose(file);
  return(error);
}

/****************************************************************************/
/*!
**   Read a byte from the binary file
**
** \return  Returned integer value of the byte
**
** \param[in]  file  File pointer
**
*****************************************************************************/
static int st_in(FILE *file)

{
  int c;

  if ((c = fgetc(file)) != EOF) return( c );
  
  if (feof(file))
    message(" End-of-file reached\n" );
  if (ferror(file))
    message(" A READ error occured (%d)\n",ferror(file));

  return(c);
}

/****************************************************************************/
/*!
**   Compose several bytes into an integer
**
** \return  Returned integer
**
** \param[in]  file  File pointer
** \param[in]  nb    Number of bytes to be considered
**
*****************************************************************************/
static int st_compose(FILE *file,
                      int nb)

{
  int i,value,factor;
  unsigned char c;

  /* Initializations */

  value  = 0;
  factor = 1;

  for (i=0; i<nb; i++)
  {
    c = st_in(file);
    value  += c * factor;
    factor *= 0x100;
  }
  return(value);
}

/****************************************************************************/
/*!
**   Read the Grid characteristics from a BMP file
**
** \return  Error return code
**
** \param[in]  file       Pointer to the BMP file
** \param[in]  verbose    1 for a verbose output
**
** \param[out]  nx        Array of number of grid nodes
** \param[out]  x0        Array of grid origin coordinates
** \param[out]  dx        Array of grid mesh
** \param[out]  ir        Array of Red indices
** \param[out]  ig        Array of Green indices
** \param[out]  ib        Array of Blue indices
** \param[out]  nbits     Number of bits per pixel
**
*****************************************************************************/
static int st_grid_read_bmp_header(FILE   *file,
                                   int     verbose,
                                   int    *nx,
                                   double *x0,
                                   double *dx,
                                   int    *ir,
                                   int    *ig,
                                   int    *ib,
                                   int    *nbits)
{
  int icol,offset,compress,ival,ncol,ndx,ndy;

  /* Initializations */

  x0[0] = 0.;
  x0[1] = 0.;
  dx[0] = 1.;
  dx[1] = 1.;

  /* Title */

  if (verbose) mestitle(0,"BMP File Characteristics :");

  /* Reading the file header */

  ival = st_compose(file,2);
  if (verbose) message("Signature = %d\n",ival);
  ival = st_compose(file,4);
  if (verbose) message("File Size = %d\n",ival);
  ival = st_compose(file,4);
  offset = st_compose(file,4);
  if (verbose) message("Offset = %d\n",offset);

  /* Reading the bitmap information header */

  ival = st_compose(file,4);
  if (verbose) message("Header size = %d\n",ival);
  nx[0] = st_compose(file,4);
  if (verbose) message("Width (pixels) = %d\n",nx[0]);
  nx[1] = st_compose(file,4);
  if (verbose) message("Height (pixels) = %d\n",nx[1]);
  ival = st_compose(file,2);
  if (verbose) message("Number of planes = %d\n",ival);
  (*nbits) = st_compose(file,2);
  if (verbose) message("Number of bits per pixel = %d\n",(*nbits));
  compress = st_compose(file,4);
  if (verbose) message("Compression = %d\n",compress);
  ival = st_compose(file,4);
  if (verbose) message("Image Size = %d\n",ival);
  ndx = st_compose(file,4);
  if (verbose) message("Pixel per meter along X = %d\n",ndx);
  ndy = st_compose(file,4);
  if (verbose) message("Pixel per meter along Y = %d\n",ndy);
  ncol = st_compose(file,4);
  if (verbose) message("Number of colors used = %d\n",ncol);
  ival = st_compose(file,4);
  if (verbose) message("Number of important colors = %d\n",ival);
  if (ncol > 256) 
  {
    messerr("Your file seems to contain more than 256 colors");
    messerr("It is probably not a valid BMP file");
    return(1);
  }

  /* Read the palette (optional) */

  for (icol=0; icol<ncol; icol++) 
  {
    ir[icol] = st_in(file);
    ig[icol] = st_in(file);
    ib[icol] = st_in(file);
    ival     = st_in(file);
    if (verbose) message("RGB (%3d) : R(%3d) G(%3d) B(%3d)\n",
                         icol+1,ir[icol],ig[icol],ib[icol]);
  }    

  /* Final checks */

  if (compress != 0)
  {
    messerr("Error : only uncompressed images are treated here");
    return(1);
  }

  /* Final results */

  if (ndx > 0) dx[0] = 100. / (double) ndx;
  if (verbose) message("Mesh along X = %lf\n",dx[0]);
  if (ndy > 0) dx[1] = 100. / (double) ndy;
  if (verbose) message("Mesh along Y = %lf\n",dx[1]);

  /* Verbose optional printout */

  if (verbose) 
  {
    mestitle(0,"BMP File Characteristics :");
    st_print_verbose(2,nx,dx,x0);
  }

  return(0);
}

/****************************************************************************/
/*!
**   Read the Grid characteristics from a BMP file
**
** \return  Error return code
**
** \param[in]  filename Name of the BMP file
** \param[in]  verbose  Verbose option
**
** \param[out]  nx        Array of number of grid nodes
** \param[out]  x0        Array of grid origin coordinates
** \param[out]  dx        Array of grid mesh
**
*****************************************************************************/
GEOSLIB_API int db_grid_read_bmp1(const char   *filename,
                                  int     verbose,
                                  int    *nx,
                                  double *x0,
                                  double *dx)
{
  FILE *file;
  int   error,nbits,ir[256],ig[256],ib[256];

  /* Initializations */

  error = 1;

 /* Open the file */

  file = gslFopen(filename,"rb");
  if (file == nullptr)
  {
    messerr("Error when opening the BMP file %s for reading",filename);
    return(1);
  }

  /* Read the grid characteristics */

  if (st_grid_read_bmp_header(file,verbose,nx,x0,dx,ir,ig,ib,
                              &nbits)) goto label_end;

  /* Set the error return code */

  error = 0;

label_end:
  if (file != nullptr) fclose(file);
  return(error);
}

/****************************************************************************/
/*!
** db_grid_read_bmp2
**
**   Read an array of values from a BMP file
**
** \return  Error return code
**
** \param[in]  filename  Name of the BMP file
** \param[in]  nx_r      Array of number of grid nodes
** \param[in]  x0_r      Array of grid origin coordinates
** \param[in]  dx_r      Array of grid mesh
**
** \param[out] tab       Array of values read
**
*****************************************************************************/
GEOSLIB_API int db_grid_read_bmp2(const char   *filename,
                                  int    *nx_r,
                                  double *x0_r,
                                  double *dx_r,
                                  double *tab)
{
  FILE  *file;
  int    nx[2],error,ix,jy,nbits,npad,noct,ecr,ir[256],ig[256],ib[256],c;
  double x0[2],dx[2];
  unsigned char ctab;

  /* Initializations */

  error = 1;

  /* Open the file */

  file = gslFopen(filename,"rb");
  if (file == nullptr)
  {
    messerr("Error when opening the BMP file %s for reading",filename);
    return(1);
  }

  /* Read the grid characteristics */

  if (st_grid_read_bmp_header(file,0,nx,x0,dx,ir,ig,ib,&nbits)) goto label_end;

  /* Check the compatibility of the input arguments with the current grid */

  if (nx[0] != nx_r[0] || nx[1] != nx_r[1] || 
      x0[0] != x0_r[0] || x0[1] != x0_r[1] ||
      dx[0] != dx_r[0] || dx[1] != dx_r[1])
  {
    messerr("The grid parameters are incompatible");
    goto label_end;
  }

  /* Reading the image (from bottom to up) */

  noct = nx[0] * nbits / 8;
  npad = 4 - noct % 4;
  if (npad == 4) npad = 0;

  ecr = 0;
  for (jy=0; jy<nx[1]; jy++)
  {
    for (ix=0; ix<nx[0]; ix++,ecr++)
    {
      ctab = 0;
      if (nbits ==  8)
      {
        c = st_in(file);
        rgb2num(ir[c],ig[c],ib[c],0,&ctab);
      }
      else if (nbits == 24)
        rgb2num(st_in(file),st_in(file),st_in(file),0,&ctab);
      else if (nbits == 32)
        rgb2num(st_in(file),st_in(file),st_in(file),st_in(file),&ctab);
      else
      {
        messerr("Error: Number of bits/pixel (%d) is not treated",nbits);
        goto label_end;
      }
      tab[ecr] = ctab;
    }

    /* Reading the padding */

    for (ix=0; ix<npad; ix++) (void) st_in(file);
  }

  /* Set the error return code */

  error = 0;

label_end:
  if (file != nullptr) fclose(file);
  return(error);
}

/****************************************************************************/
/*!
**   Write the Grid in a ZYCOR file
**
** \return  Error return code
**
** \param[in]  filename  Name of the ZYCOR file
** \param[in]  db        Db structure to be written
** \param[in]  icol      Rank of the attribute
**
*****************************************************************************/
GEOSLIB_API int db_grid_write_zycor(const char *filename,
                                    Db     *db,
                                    int     icol)
{
  FILE  *file;
  int    i,nx[2],jj,ii,kk,yy,ind,loop;
  double rbid,x0[2],xf[2],dx[2];
  double buff[5];               /* Size = nbyline */
  char   card[100];             /* Size = nbyline * 20 */
  static int    nbyline = 5;
  static double testval = 1.e30;

  /* Preliminary checks */

  if (! is_grid(db) || db->getNDim() != 2)
  {
    messerr("The Db structure should correspond to a 2-D Grid");
    return(1);
  }

  /* Open the file */

  file = gslFopen(filename,"w");
  if (file == nullptr)
  {
    messerr("Error when opening the ZYCOR file %s for writing",filename);
    return(1);
  }

  /* Write a comment */

  fprintf(file,"!\n");
  fprintf(file,"!  File created by RGeostats package\n");
  fprintf(file,"!\n");

  /* Title line */

  fprintf(file,"@GRID ZYCOR FILE    ,   GRID,  %d\n",nbyline);
  fprintf(file,"     15, %13lg,    ,    0,     1\n",testval);

  /* Grid description */

  for (i=0; i<2; i++)
  {
    nx[i] = db->getNX(i);
    x0[i] = db->getX0(i);
    dx[i] = db->getDX(i);
    xf[i] = x0[i] + (nx[i] - 1) * dx[i];
  }
  
  rbid = 0.;
  fprintf(file,"%6d, %6d, %13lf, %13lf, %13lf, %13lf\n",
          nx[1],nx[0],x0[0],xf[0],x0[1],xf[1]);
  fprintf(file," %15lf, %15lf, %15lf\n",rbid,rbid,rbid);
  fprintf(file,"@\n");

  /* The set of values */

  for ( jj = nx[0]-1; jj>=0; jj-- )
  {
     kk = 0;
     ii = ( ( nx[1] * nx[0] ) - ( jj+1 ) );
     for ( loop=1; loop <= nx[1]; loop++ )
     {
        buff[kk++] = db->getArray(ii,icol);
        ii -= nx[0];
        if ( kk == nbyline )
        {
           for ( yy=0; yy<nbyline; yy++ )
           {
              ind = yy * 15;
              if (! FFFF(buff[yy]))
              {
                 sprintf ( &card[ind], "%15g", buff[yy] );
              }
              else
              {
                 memcpy ( &card[ind], (char *)ZYCOR_NULL_CH, 15 );
              }
           }
           sprintf ( &card[15*nbyline], "\n" );
           fprintf ( file, "%s", card );
           kk = 0;
        }
     }

     if ( kk > 0 ) 
     {
        for ( yy=0; yy<kk; yy++ )
        {
           ind = yy * 15;
           if (! FFFF(buff[yy]))
           {
              sprintf ( &card[ind], "%15g", buff[yy] );
           }
           else
           {
              memcpy ( &card[ind], (char *)ZYCOR_NULL_CH, 15 );
           }
        }
        sprintf ( &card[15*kk], "\n" );
        fprintf ( file, "%s", card );
     }
  }

  if (file != nullptr) fclose(file);
  return(0);
}

/****************************************************************************/
/*!
**   Print an integer 
**
** \param[in]  file  File pointer
** \param[in]  mode  Type of writing
** \li                0 : 16-bit unsigned integer
** \li                1 : 32-bit unsigned integer
** \li                2 : 32-bit signed integer
** \param[in]  ival  Integer value to be written
**
*****************************************************************************/
static void st_out(FILE *file, int mode, unsigned int ival)
{
  int jval;

  switch (mode)
    {
    case 0:			/* Unsigned 16-bit */
      putc(ival, file);
      putc(ival >> 8, file);
      break;

    case 1:			/* Unsigned 32-bit */
      putc(ival, file);
      putc(ival >> 8, file);
      putc(ival >> 16, file);
      putc(ival >> 24, file);
      break;

    case 2:			/* Signed 32-bit */
      jval = (int) ival;
      putc(jval, file);
      putc(jval >> 8, file);
      putc(jval >> 16, file);
      putc(jval >> 24, file);
      break;
    }
}

/****************************************************************************/
/*!
**   Write the contents of a Grid into a BMP file
**
** \return  Error return code
**
** \param[in]  filename    Name of the BMP file
** \param[in]  db          Db structure to be written
** \param[in]  icol        Rank of the attribute
** \param[in]  nsamplex    Sampling ratio along X (1 = no sampling)
** \param[in]  nsampley    Sampling ratio along Y (1 = no sampling)
** \param[in]  nmult       Multiplication factor
** \param[in]  ncolor      Number of colors
** \param[in]  flag_low    Returns the first color for value below vmin if 0
**                         or returns -3 otherwise
** \param[in]  flag_high   Returns the last color for value above vmax if 0
**                         or returns -4 otherwise
** \param[in]  valmin      Minimum value represented (or FFFF)
** \param[in]  valmax      Maximum value represented (or FFFF)
** \param[in]  red         Array of Red intensity for color scale
** \param[in]  green       Array of Green intensity for color scale
** \param[in]  blue        Array of Blue intensity for color scale
** \param[in]  mask_red    Red intensity for masked value
** \param[in]  mask_green  Green intensity for masked value
** \param[in]  mask_blue   Blue intensity for masked value
** \param[in]  ffff_red    Red intensity for FFFF value
** \param[in]  ffff_green  Green intensity for FFFF value
** \param[in]  ffff_blue   Blue intensity for FFFF value
** \param[in]  low_red     Red intensity for lower value
** \param[in]  low_green   Green intensity for lower value
** \param[in]  low_blue    Blue intensity for lower value
** \param[in]  high_red    Red intensity for higher value
** \param[in]  high_green  Green intensity for higher value
** \param[in]  high_blue   Blue intensity for higher value
**  
** \remark  If ncolor=0, the colors scale is generated as the grey scale
**
*****************************************************************************/
GEOSLIB_API int db_grid_write_bmp(const char *filename,
                                  Db     *db,
                                  int     icol,
                                  int     nsamplex,
                                  int     nsampley,
                                  int     nmult,
                                  int     ncolor,
                                  int     flag_low,
                                  int     flag_high,
                                  double  valmin,
                                  double  valmax,
                                  int    *red,
                                  int    *green,
                                  int    *blue,
                                  int     mask_red,
                                  int     mask_green,
                                  int     mask_blue,
                                  int     ffff_red,
                                  int     ffff_green,
                                  int     ffff_blue,
                                  int     low_red,
                                  int     low_green,
                                  int     low_blue,
                                  int     high_red,
                                  int     high_green,
                                  int     high_blue)
{
  FILE    *file;
  int     *indg,infosize,headersize,imagesize,ix,iy,i,nx,ny,iech,rank;
  int      imult,jmult,ipad,flag_color_scale,idim,number,color,width,height;
  unsigned char ired,igreen,iblue;
  double   value,vmin,vmax;
  
  /* Preliminary checks */

  indg = nullptr;
  file = nullptr;
  flag_color_scale = (ncolor > 0 && 
                      red   != nullptr &&
                      green != nullptr &&
                      blue  != nullptr);
  if (! flag_color_scale) ncolor = 256;
  if (! is_grid(db))
  {
    messerr("The Db structure should correspond to a Grid");
    return(1);
  }
  if (db->getNDim() > 2)
  {
    number = 1;
    for (idim=2; idim<db->getNDim(); idim++) number *= db->getNX(idim);
    if (number > 1)
    {
      messerr("The Db structure corresponds to a 3-D Grid");
      messerr("Only the first XOY plane is written");
      return(1);
    }
  }

  /* Core allocation */

  indg = db_indg_alloc(db);
  if (indg == nullptr) goto label_end;

  /* Initializations */

  nx = db->getNX(0);
  ny = db->getNX(1);

  /* Calculate the statistics */

  vmin =  1.e30;
  vmax = -1.e30;
  for (i=0; i<nx * ny; i++)
  {
    if (! db->getSelection(i)) continue;
    value = db->getArray(i,icol);
    if (FFFF(value)) continue;
    if (value < vmin) vmin = value;
    if (value > vmax) vmax = value;
  }
  if (! FFFF(valmin)) vmin = valmin;
  if (! FFFF(valmax)) vmax = valmax;

 /* Open the file */

  file = gslFopen(filename,"wb");
  if (file == nullptr)
  {
    messerr("Error when opening the BMP file %s for writing",filename);
    return(1);
  }

  /* Figure out the constants */
  infosize   = 40;
  headersize = 14;
  width      = nmult * N_SAMPLE(nx,nsamplex);
  height     = nmult * N_SAMPLE(ny,nsampley);
  imagesize  = 3 * width * height;
  
  /* Write the file header, bitmap information, and bitmap pixel data... */
  st_out(file,0, BF_TYPE);
  st_out(file,1, headersize);	/* Size of File Header */
  st_out(file,0, 0);		/* Reserved */
  st_out(file,0, 0);		/* Reserved */
  st_out(file,1, headersize + infosize); /* Offset */
  
  st_out(file,1, infosize);	/* Size of Information Block */
  st_out(file,2, width);	/* Width */
  st_out(file,2, height);	/* Height */
  st_out(file,0, 1);		/* Number of planes */
  st_out(file,0, 24);		/* 24-bits per pixel */
  st_out(file,1, 0);		/* No compression */
  st_out(file,1, imagesize);	/* Image size */
  st_out(file,2, 0);
  st_out(file,2, 0);
  st_out(file,1, 0);
  st_out(file,1, 0);

  /* Writing the pixels */

  indg[0] = 0;
  indg[1] = 0;
  indg[2] = 0;
  ipad = nx * nmult;
  ipad = ipad - 4 * ((int) (ipad/4));
  for (iy=0; iy<ny; iy++)
  {
    if (iy % nsampley != 0) continue;
    for (jmult=0; jmult<nmult; jmult++)
    {
      for (ix=0; ix<nx; ix++)
      {
	if (ix % nsamplex != 0) continue;
	indg[0] = ix;
	indg[1] = iy;
	iech  = db_index_grid_to_sample(db,indg);
	rank  = st_color_rank(db,iech,icol,ncolor,flag_low,flag_high,vmin,vmax);
	st_color_in_rgb(rank,flag_color_scale,red,green,blue,
			mask_red,mask_green,mask_blue,
			ffff_red,ffff_green,ffff_blue,
			low_red , low_green, low_blue,
			high_red,high_green,high_blue,
			&ired,&igreen,&iblue);

	for (imult=0; imult<nmult; imult++)
        {
          (void) fwrite(&iblue ,1,1,file);
          (void) fwrite(&igreen,1,1,file);
          (void) fwrite(&ired  ,1,1,file);
        }
      }
      
      /* Write the padding */

      for (i=0; i<ipad; i++)
      {
        color = 0;
        (void) fwrite(&color,1,1,file);
      }
    }
  }
  
 label_end:
  indg = db_indg_free(indg);
  if (file != nullptr) fclose(file);
  return(0);
}

/****************************************************************************/
/*!
**   Dump an attribute from a Db into an ASCII file using the
**   IRAP format
**
** \param[in]  filename  Name of the output ASCII file
** \param[in]  db        Target Db
** \param[in]  icol      Rank of the target attribute
** \param[in]  nsamplex  Sampling ratio along X (1 = no sampling)
** \param[in]  nsampley  Sampling ratio along Y (1 = no sampling)
**
** \remark  If filename is not specified, the print is routed towards
** \remark  the standard output
**
*****************************************************************************/
GEOSLIB_API int db_grid_write_irap(const char *filename,
                                   Db   *db,
                                   int   icol,
                                   int   nsamplex,
                                   int   nsampley)
{
  FILE  *file;
  double xmin,xmax,ymin,ymax,value,dx,dy;
  int   *indg,iech,ix,iy,nx,ny,necr,error;

  /* Initializations */
  error = 1;
  file  = nullptr;
  indg  = nullptr;

  /* Check that the file is a grid */
  if (db == NULL) return(1);
  if (! is_grid(db))
  {
    messerr("The IRAP dump is only designed for a Grid DB");
    goto label_end;
  }
  if (db->getNDim() != 2)
  {
    messerr("The IRAP dump is only designed for 2-D Grid");
    goto label_end;
  }

  /* Core allocation */

  indg = db_indg_alloc(db);
  if (indg == nullptr) goto label_end;

  /* Open the output file */

  file = gslFopen(filename,"w");
  if (file == nullptr) goto label_end;

  /* Preliminary calculations */
  nx   = N_SAMPLE(db->getNX(0),nsamplex);
  ny   = N_SAMPLE(db->getNX(1),nsampley);
  dx   = db->getDX(0) * nsamplex;
  dy   = db->getDX(1) * nsampley;
  xmin = db->getX0(0);
  ymin = db->getX0(1);
  xmax = xmin + dx * (nx - 1);
  ymax = ymin + dy * (ny - 1);

  /* Write the header */
  fprintf(file,"%d %d %lf %lf\n",nx,ny,dx,dy);
  fprintf(file,"%lf %lf %lf %lf\n",xmin,xmax,ymin,ymax);

  necr = 0;
  for (iy=0; iy<ny; iy++)
  {
    if (iy % nsampley != 0) continue;
    for (ix=0; ix<nx; ix++)
    {
      if (ix % nsamplex != 0) continue;
      indg[0] = ix;
      indg[1] = iy;
      iech    = db_index_grid_to_sample(db,indg);
      value   = db->getArray(iech,icol);
      if (FFFF(value)) value = 9999990.;
      fprintf(file,"%10.3lf ",value);
      necr++;
      if (necr == 6)
      {
        fprintf(file,"\n");
        necr = 0;
      }
    }
  }
  if (necr > 0) fprintf(file,"\n");

  /* Set the error return code */

  error = 0;

 label_end:
  indg = db_indg_free(indg);
  if (file != nullptr) fclose(file);
  return(error);
}

/****************************************************************************/
/*!
**   Write the Grid in a IFPEN .PROP file
**
** \return  Error return code
**
** \param[in]  filename  Name of the IFPEN .PROP file
** \param[in]  db        Db structure to be written
** \param[in]  ncol      Number of attributes
** \param[in]  icols     Rank(s) of the attribute
**
*****************************************************************************/
GEOSLIB_API int db_grid_write_prop(const char *filename,
                                   Db     *db,
                                   int     ncol,
                                   int    *icols)
{
  FILE   *file;
  int     i,j,idim,ndim,flag_not_rotz,ntot,nx[3];
  double  value;
  VectorDouble angles;
  static  double valnull = 3.0;

  /* Preliminary checks */

  if (! is_grid(db))
  {
    messerr("The Db structure should correspond to a Grid");
    return(1);
  }
  ndim = db->getNDim();
  ntot = 1;
  for (idim=0; idim<3; idim++) 
  {
    nx[idim] = (idim < ndim) ? db->getNX(idim) : 1;
    ntot *= nx[idim];
  }

  if (db->isGridRotated())
  {
    angles = db->getGrid().getRotAngles();
    flag_not_rotz = 0;
    for (idim=1; idim<ndim; idim++)
      if (angles[idim] != 0.) flag_not_rotz = 1;
    if (flag_not_rotz)
    {
      messerr("The Grid rotation may not involve Oy nor Ox angles");
     return(1);
    }
  }

  /* Open the file */

  file = gslFopen(filename,"w");
  if (file == nullptr)
  {
    messerr("Error when opening the IFPEN .PROP file %s for writing",filename);
    return(1);
  }

  /* Write the header */

  st_ifpen_write(file,0,"##########################",0,0.,NULL);
  st_ifpen_write(file,0,"FILE_DESCRIPTION         # PROP",0,0.,NULL);
  st_ifpen_write(file,0,"APPLICATION              #",0,0.,"# CobraFlow");
  st_ifpen_write(file,0,"SURVEY_NAME              #",0,0.,NULL);
  st_ifpen_write(file,0,"MATRIX_NAME              # VPCMatrix_test_export",
		 0,0.,NULL);
  st_ifpen_write(file,0,"METHOD                   # BY_CPV",0,0.,NULL);
  st_ifpen_write(file,2,"FLOAT_NULL_VALUE         #",0,valnull,NULL);
  st_ifpen_write(file,0,"ROW_COLUMN_ORIENTATION   # ROW",0,0.,NULL);
  st_ifpen_write(file,0,"REPRESENTATION_CODE      # ASCII",0,0.,NULL);
  st_ifpen_write(file,0,"##########################",0,0.,NULL);
  st_ifpen_write(file,2,"ANGLE                    #",0,angles[0],"# DEG");
  st_ifpen_write(file,1,"ROW_COUNT                #",nx[1],0.,NULL);
  st_ifpen_write(file,1,"COLUMN_COUNT             #",nx[0],0.,NULL);
  st_ifpen_write(file,2,"ROW_DISTANCE             #",0,db->getDX(1),"# m");
  st_ifpen_write(file,2,"COLUMN_DISTANCE          #",0,db->getDX(0),"# m");
  st_ifpen_write(file,1,"LAYER_COUNT              #",nx[2],0.,NULL);
  st_ifpen_write(file,2,"X_ORIGIN                 #",0,db->getX0(0),"# m");
  st_ifpen_write(file,2,"Y_ORIGIN                 #",0,db->getX0(1),"# m");
  st_ifpen_write(file,1,"FACIES_COUNT             #",ncol,0.,NULL);
  st_ifpen_write(file,0,"DATA_PROP                # CHANNEL1",0,0.,
		 "# Facies proportion");
  st_ifpen_write(file,0,"##########################",0,0.,NULL);

  /* Grid description */

  for (j=0; j<ncol; j++)
    for (i=0; i<ntot; i++)
    {
      value = db->getArray(i,icols[j]);
      st_ifpen_write(file,2,NULL,0,value,NULL);
    }

  if (file != nullptr) fclose(file);
  return(0);
}

/****************************************************************************/
/*!
**   Write the Grid in a ECLIPSE file
**
** \return  Error return code
**
** \param[in]  filename  Name of the ECLIPSE file
** \param[in]  db        Db structure to be written
** \param[in]  icol      Rank of the attribute
**
*****************************************************************************/
GEOSLIB_API int db_grid_write_eclipse(const char   *filename,
                                      Db     *db,
                                      int     icol)
{
  FILE  *file;
  int    i,idim,nxyz,ninline;
  double value,valprt;
  static int nbyline = 6;
  static double valtest = -9999.;

  /* Preliminary checks */

  if (! is_grid(db))
  {
    messerr("The Db structure should correspond to a Grid");
    return(1);
  }
  nxyz = 1;
  for (idim=0; idim<db->getNDim(); idim++) nxyz *= db->getNX(idim);

  /* Open the file */

  file = gslFopen(filename,"w");
  if (file == nullptr)
  {
    messerr("Error when opening the ECLIPSE file %s for writing",filename);
    return(1);
  }

  /* Write a comment */

  fprintf(file,"Facies\n");

  /* Write the set of values */
  
  ninline = 0;
  for (i=0; i<nxyz; i++)
  {
    valprt = valtest;
    if (db->getSelection(i))
    {
      value = db->getArray(i,icol);
      if (! FFFF(value)) valprt = value;
    }
    fprintf(file,"%lf ",valprt);
    ninline++;
    if (ninline == nbyline)
    {
      fprintf(file,"\n");
      ninline = 0;
    }
  }
  if (ninline > 0) fprintf(file,"\n");

  if (file != nullptr) fclose(file);
  return(0);
}

/****************************************************************************/
/*!
**   Read the header from a IFPEN .PROP file
**
** \return  Error return code
**
** \param[in]  file      Pointer to the IFPEN .PROP file
** \param[in]  verbose   1 for a verbose option; 0 otherwise
**
** \param[out]  ncol      Number of attributes
** \param[out]  nx        Array of number of grid nodes
** \param[out]  x0        Array of grid origin coordinates
** \param[out]  dx        Array of grid mesh
** \param[out]  test      Value for non-information
**
*****************************************************************************/
static int st_grid_read_prop_header(FILE   *file,
                                    int     verbose,
                                    int    *ncol,
                                    int    *nx,
                                    double *x0,
                                    double *dx,
                                    double *test)
{
  int    i,error,dumint;
  double dumrel,anglez;

  /* Initializations */

  error = 1;
  for (i=0; i<3; i++)
  {
    nx[i] = 1;
    dx[i] = 1.;
    x0[i] = 0.;
  }

  /* Read the header */

  if (st_ifpen_read(file,0,"##########################",
		    &dumint,&dumrel)) goto label_end;
  if (st_ifpen_read(file,0,"FILE_DESCRIPTION         #",
		    &dumint,&dumrel)) goto label_end;
  if (st_ifpen_read(file,0,"APPLICATION              #",
		    &dumint,&dumrel)) goto label_end;
  if (st_ifpen_read(file,0,"SURVEY_NAME              #",
		    &dumint,&dumrel)) goto label_end;
  if (st_ifpen_read(file,0,"MATRIX_NAME              #",
		    &dumint,&dumrel)) goto label_end;
  if (st_ifpen_read(file,0,"METHOD                   #",
		    &dumint,&dumrel)) goto label_end;
  if (st_ifpen_read(file,2,"FLOAT_NULL_VALUE         #",
		    &dumint,test)) goto label_end;
  if (st_ifpen_read(file,0,"ROW_COLUMN_ORIENTATION   #",
		    &dumint,&dumrel)) goto label_end;
  if (st_ifpen_read(file,0,"REPRESENTATION_CODE      #",
		    &dumint,&dumrel)) goto label_end;
  if (st_ifpen_read(file,0,"##########################",
		    &dumint,&dumrel)) goto label_end;
  if (st_ifpen_read(file,2,"ANGLE                    #",
		    &dumint,&anglez)) goto label_end;
  if (st_ifpen_read(file,1,"ROW_COUNT                #",
		    &nx[1],&dumrel)) goto label_end;
  if (st_ifpen_read(file,1,"COLUMN_COUNT             #",
		    &nx[0],&dumrel)) goto label_end;
  if (st_ifpen_read(file,2,"ROW_DISTANCE             #",
		    &dumint,&dx[1])) goto label_end;
  if (st_ifpen_read(file,2,"COLUMN_DISTANCE          #",
		    &dumint,&dx[0])) goto label_end;
  if (st_ifpen_read(file,1,"LAYER_COUNT              #",
		    &nx[2],&dumrel)) goto label_end;
  if (st_ifpen_read(file,2,"X_ORIGIN                 #",
		    &dumint,&x0[0])) goto label_end;
  if (st_ifpen_read(file,2,"Y_ORIGIN                 #",
		    &dumint,&x0[1])) goto label_end;
  if (st_ifpen_read(file,1,"FACIES_COUNT             #",
		    ncol,&dumrel)) goto label_end;
  if (st_ifpen_read(file,0,"DATA_PROP                #",
		    &dumint,&dumrel)) goto label_end;
  if (st_ifpen_read(file,0,"##########################",
		    &dumint,&dumrel)) goto label_end;

  /* The rotation angle is not valid */

  if (anglez != 0.)
  {
    messerr("The rotation angle of the IFPEN .PROP file is not valid here");
    goto label_end;
  }

  /* Verbose optional printout */

  if (verbose) 
  {
    mestitle(0,"IFPEN .PROP File Characteristics :");
    st_print_verbose(3,nx,dx,x0);
  }

  error = 0;

 label_end:
  return(error);
}

/****************************************************************************/
/*!
**   Read the Grid characteristics from a IFPEN .PROP file
**
** \return  Error return code
**
** \param[in]  filename  Name of the IFPEN .PROP file
** \param[in]  verbose   Verbose flag
**
** \param[out]  ncol      Number of attributes
** \param[out]  nx        Array of number of grid nodes
** \param[out]  x0        Array of grid origin coordinates
** \param[out]  dx        Array of grid mesh
**
*****************************************************************************/
GEOSLIB_API int db_grid_read_prop1(const char *filename,
                                   int     verbose,
                                   int    *ncol,
                                   int    *nx,
                                   double *x0,
                                   double *dx)
{
  FILE  *file;
  double test;
  int    error;

  /* Initializations */

  error = 1;

 /* Open the file */

  file = gslFopen(filename,"r");
  if (file == nullptr)
  {
    messerr("Error when opening the IFPEN .PROP file %s for reading",filename);
    return(1);
  }

  /* Read the grid characteristics */

  if (st_grid_read_prop_header(file,verbose,ncol,nx,x0,dx,
                               &test)) goto label_end;

  /* Set the error return code */

  error = 0;

label_end:
  fclose(file);
  return(error);
}

/****************************************************************************/
/*!
**   Read an array of values from a IFPEN .PROP file
**
** \return  Error return code
**
** \param[in]  filename  Name of the IFPEN .PROP file
** \param[in]  ncol_r    Number of attributes to be read
** \param[in]  nx_r      Array of number of grid nodes
** \param[in]  x0_r      Array of grid origin coordinates
** \param[in]  dx_r      Array of grid mesh
**
** \param[out] tab       Array of values
**
*****************************************************************************/
GEOSLIB_API int db_grid_read_prop2(const char   *filename,
                                   int     ncol_r,
                                   int    *nx_r,
                                   double *x0_r,
                                   double *dx_r,
                                   double *tab)
{
  FILE  *file;
  int    nx[3],nech,lec,error,ix,iy,iz,ncol,icol;
  double x0[3],dx[3],test,value;

  /* Initializations */

  error = 1;
  file  = nullptr;

  /* Open the file */

  file = gslFopen(filename,"r");
  if (file == nullptr)
  {
    messerr("Error when opening the IFPEN .PROP file %s for reading",filename);
    goto label_end;
  }

  /* Read the grid characteristics */

  if (st_grid_read_prop_header(file,0,&ncol,nx,x0,dx,
                               &test)) goto label_end;

  /* Check the compatibility of the input arguments with the current grid */

  if (ncol  != ncol_r  ||
      nx[0] != nx_r[0] || nx[1] != nx_r[1] || 
      x0[0] != x0_r[0] || x0[1] != x0_r[1] ||
      dx[0] != dx_r[0] || dx[1] != dx_r[1])
  {
    messerr("The grid parameters are incompatible");
    goto label_end;
  }

  /* Read the array of real values */

  lec  = 0;
  nech = nx[0] * nx[1] * nx[2];

  for (icol=0; icol<ncol; icol++)
    for (ix=0; ix<nx[0]; ix++)
      for (iy=0; iy<nx[1]; iy++)
        for (iz=0; iz<nx[2]; iz++)
        {
          if (_record_read(file,"%lf",&value)) break;
          if (value == test) value = TEST;
          tab[lec] = value;
          lec++;
        }

  if (lec != nech * ncol)
  {
    messerr("The number of decoded values (%d) is not equal to the expected number of grid nodes (%d) x number of attributes (%d)",
            lec,nech,ncol);
    goto label_end;
  }

  /* Set the error return code */

  error = 0;

label_end:
  record_close();
  if (file != nullptr) fclose(file);
  return(error);
}

/****************************************************************************/
/*!
**   Write the Db in a VTK file
**
** \return  Error return code
**
** \param[in]  filename  Name of the VTK file
** \param[in]  db        Db structure to be written
** \param[in]  cols      Rank(s) of the attribute
** \param[in]  names     Array of "selected" variable names
**
*****************************************************************************/
GEOSLIB_API int db_write_vtk(const char *filename,
                             Db     *db,
                             const VectorInt& cols,
                             const VectorString& names)
{
  int    *vardim,*center;
  int     dims[3],error,nech,ndim,flag_grid,useBinary,ecr,iad,nactive,ncol;
  float  *points,*xcoor,*ycoor,*zcoor,factx,facty,factz,fact,factvar,value;
  float **tab;
  std::vector<char *> vc;

  /* Preliminary checks */

  ndim = db->getNDim();
  ncol = static_cast<int> (cols.size());
  if (ndim > 3)
  {
    messerr("VTK files are limited to 3D");
    return(1);
  }
  if (ncol <= 0)
  {
    messerr("You must define a set of columns to be written");
    return(1);
  }

  /* Initializations */

  error     = 1;
  vardim    = center = nullptr;
  tab       = nullptr;
  points    = xcoor = ycoor = zcoor = nullptr;
  nech      = db->getSampleNumber();
  nactive   = db->getActiveSampleNumber();
  flag_grid = is_grid(db);
  useBinary = (int) get_keypone("VTK_Use_Binary",1.);
  factx     = get_keypone("VTK_Fact_X",1.);
  facty     = get_keypone("VTK_Fact_Y",1.);
  factz     = get_keypone("VTK_Fact_Z",1.);
  factvar   = get_keypone("VTK_Fact_Var",1.);

  /* Define the reading parameters */

  for (int idim=0; idim<3; idim++)
    dims[idim] = (idim < ndim) ? db->getNX(idim) : 1;

  /* Core allocation */

  vardim = (int *)    mem_alloc(sizeof(int)   * ncol,0);
  if (vardim == nullptr) goto label_end;
  center = (int *)    mem_alloc(sizeof(int)   * ncol,0);
  if (center == nullptr) goto label_end;
  for (int icol=0; icol<ncol; icol++)
  {
    vardim[icol] = 1;
    center[icol] = 1;
  }
  tab    = (float **) mem_alloc(sizeof(float *) * ncol,0);
  if (tab == nullptr) goto label_end;
  for (int icol=0; icol<ncol; icol++)
  {
    if (flag_grid)
      tab[icol] = (float *) mem_alloc(sizeof(float) * nech,0);
    else
      tab[icol] = (float *) mem_alloc(sizeof(float) * nactive,0);
    if (tab[icol] == nullptr) goto label_end;
  }
  if (flag_grid)
  {
    xcoor = (float *) mem_alloc(sizeof(float) * dims[0],0);
    if (xcoor == nullptr) goto label_end;
    for (int i=0; i<dims[0]; i++) 
      xcoor[i] = (float) (factx * (db->getX0(0) + i * db->getDX(0)));
    ycoor = (float *) mem_alloc(sizeof(float) * dims[1],0);
    if (ycoor == nullptr) goto label_end;
    for (int i=0; i<dims[1]; i++) 
      ycoor[i] = (float) (facty * (db->getX0(1) + i * db->getDX(1)));
    zcoor = (float *) mem_alloc(sizeof(float) * dims[2],0);
    if (zcoor == nullptr) goto label_end;
    for (int i=0; i<dims[2]; i++) 
      zcoor[i] = (float) (factz * (db->getX0(2) + i * db->getDX(2)));
  }
  else
  {
    points = (float *) mem_alloc(sizeof(float) * 3 * nactive,0);
    if (points == nullptr) goto label_end;
  }

  /* Read the coordinates (for points only) */

  if (! flag_grid)
  {
    ecr = 0;
    for (int iech=0; iech<nech; iech++)
    {
      if (! db->isActive(iech)) continue;
      for (int idim=0; idim<3; idim++)
      {
        fact = 1.;
        if (idim == 0) fact = factx;
        if (idim == 1) fact = facty;
        if (idim == 2) fact = factz;
        points[ecr++] = (idim < ndim) ? (float) (fact * db->getCoordinate(iech,idim)) : 0.;
      }
    }
  }

  /* Load the array */

  for (int icol=0; icol<ncol; icol++)
  {
    if (! flag_grid)
    {
      ecr = 0;
      for (int iech=0; iech<nech; iech++)
        if (db->isActive(iech))
        {
          value = (float) (db->getArray(iech,cols[icol]));
          if (FFFF(value))
            tab[icol][ecr] = (float) (TEST);
          else
            tab[icol][ecr] = factvar * value;
          ecr++;
        }
    }
    else
    {
      ecr = 0;
      for (int iz=0; iz<dims[2]; iz++)
        for (int iy=0; iy<dims[1]; iy++)
          for (int ix=0; ix<dims[0]; ix++)
          {
            iad = ix + dims[0] * (iy + dims[1] * iz);
            if (db->isActive(iad))
            {
              value = (float) (db->getByColumn(iad,cols[icol]));
              if (FFFF(value)) 
                tab[icol][ecr] = (float) (TEST);
              else
                tab[icol][ecr] = factvar * value;
            }
            else
              tab[icol][ecr] = (float) (TEST);
            ecr++;
          }
    }
  }

  vc = util_vs_to_vs(names);

  /* Write the file */

  if (flag_grid)
    write_rectilinear_mesh(filename,useBinary,dims,xcoor,ycoor,zcoor,
                           ncol,vardim,center,vc.data(),tab);
  else
    write_point_mesh(filename,useBinary,nactive,points,
                     ncol,vardim,vc.data(),tab);

  /* Set the error return code */

  error = 0;

label_end:
  vardim = (int    *) mem_free((char *) vardim);
  center = (int    *) mem_free((char *) center);
  if (tab != nullptr)
  {
    for (int icol=0; icol<ncol; icol++)
      tab[icol] = (float *) mem_free((char *) tab[icol]);
    tab = (float **) mem_free((char *) tab);
  }
  points = (float *) mem_free((char *) points);
  xcoor  = (float *) mem_free((char *) xcoor);
  ycoor  = (float *) mem_free((char *) ycoor);
  zcoor  = (float *) mem_free((char *) zcoor);
  return(error);
}

/****************************************************************************/
/*!
**   Read the next line from a LAS file
**
** \return  Error return code
**
** \param[in]      s_length   Length of the string array
** \param[in]      file       File structure
** \param[in]      flag_up    Convert to uppercase
** \param[in,out]  numline    Rank of the line
** \param[out] string         New line
**
*****************************************************************************/
static int st_read_next(int   s_length,
                        FILE *file,
                        int   flag_up,
                        int  *numline,
                        char *string)
{
  int size;

  (*numline)++;
  if (fgets(string,s_length,file) == NULL) return(1);
  size = static_cast<int> (strlen(string));

  // Suppress the trailing newline
  if (string[size-1] == '\n') string[size-1] = '\0';
  
  // Convert to uppercase (optional)
  if (flag_up) string_to_uppercase(string);

  return(0);
}

/****************************************************************************/
/*!
**  Read the next lines from a LAS file until the search file is found
**
** \return  0 if the string is found; 1 for end-of-file
**
** \param[in]      s_length   Length of the string array
** \param[in]      file       File structure
** \param[in]      target     Target line to be found
** \param[in,out]  numline    Rank of the line
** \param[out]     string     New line
**
*****************************************************************************/
static int st_read_find(int   s_length,
                        FILE *file,
                        const char *target,
                        int  *numline,
                        char *string)
{
  char big_target[1000];
  
  (void) gslStrcpy(big_target,target);
  string_to_uppercase(big_target);

  /* Check the current line */
  if (strstr(string,big_target)) return(0);
  
  while (1)
  {
    if (st_read_next(s_length,file,1,numline,string)) return(1);
    if (strstr(string,big_target)) return(0);
  }
  return(1);
}

/****************************************************************************/
/*!
**   Read the characteristics from a LAS file
**
** \return  Error return code
**
** \param[in]  filename   Pointer to the LAS file
** \param[in]  verbose    1 for a verbose output; 0 otherwise
** \param[in]  xwell      Coordinate of the well (along X)
** \param[in]  ywell      Coordinate of the well (along Y)
** \param[in]  cwell      Code of the well
**
** \param[out] nvarout    Number of variables
** \param[out] nechout    Number of samples
** \param[out] var_names  Array of variable names
** \param[out] tab        Array containing variables (sorted by sample)
**
** \remarks The arrays 'var_names' and 'tab' must be freed by calling functions
**
*****************************************************************************/
GEOSLIB_API int db_well_read_las(const char   *filename,
                                 int      verbose,
                                 double   xwell,
                                 double   ywell,
                                 double   cwell,
                                 int     *nvarout,
                                 int     *nechout,
                                 char  ***var_names,
                                 double **tab)
{
  FILE *file;
  char string[1000],*lcur,sep_blank[2],sep_point[2],*token,**varloc;
  int  error,numline,nvar,nech,ecr,nquant,nvarlu;
  double *tabloc,test,value;
  static int s_length = 1000;
  static int sizemax = 10;
  static int nvar_special = 3;
  static int quantum = 10000;
  const  char *name_special[] = {"X", "Y", "CODE"};

  // Initializations */

  error = 1;
  test = TEST;
  sep_blank[0] = ' ';
  sep_point[0] = '.';
  *nvarout = *nechout = 0;
  file   = nullptr;
  varloc = nullptr;
  tabloc = nullptr;
  
  // Open the file

  file = gslFopen(filename,"r");
  if (file == nullptr)
  {
    messerr("Error when opening the LAS file %s for reading",filename);
    goto label_end;
  }
  
  // Add the names of the special variables
  if (verbose) message("Adding Special variables (coordinates)\n");
  nvar = 0;
  varloc = (char **)
    mem_realloc((char *) varloc,sizeof(char *) * nvar_special,1);
  for (int ivar=0; ivar<nvar_special; ivar++)
  {
    varloc[ivar] = (char *) mem_alloc(sizeof(char) * sizemax,1);
    (void) gslStrcpy(varloc[ivar], name_special[ivar]);
    nvar++;
  }

  numline = 0;
  (void) gslStrcpy(string,"");

  // Decode the header
  if (verbose) message("Looking for ~Well\n");
  if (st_read_find(s_length,file,"~Well",&numline,string)) goto label_end;
  while (1)
  {
    if (st_read_next(s_length,file,1,&numline,string)) goto label_end;

    // Looking for the TEST value (keyword "NULL")
    if (strstr(string,"NULL"))
      sscanf(&string[8],"%lf",&test);

    // Looking for the next delimitor
    if (strstr(string,"~")) break;
  }

  // Decode the variable names
  if (verbose) message("Decoding the List of Variables\n");
  if (st_read_find(s_length,file,"~Curve",&numline,string)) goto label_end;
  while (1)
  {
    if (st_read_next(s_length,file,1,&numline,string)) goto label_end;

    // Skipping the comments
    if (strstr(string,"#")) continue;

    // Looking for the variable header
    if (strstr(string,"~")) break;

    // Reading the variable name
    token = gslStrtok(string, sep_point);
    if (token == NULL) break;

    // Add the variable to the list
    string[strlen(token)] = '\0';
    varloc = (char **) mem_realloc((char *) varloc,sizeof(char *) * (nvar+1),1);
    varloc[nvar] = (char *) mem_alloc(sizeof(char) * sizemax,1);
    (void) gslStrncpy(varloc[nvar], token, sizemax);
    string_strip_blanks(varloc[nvar],0);
    if (verbose) message("Variable: %s\n",varloc[nvar]);
    nvar++;
  }
  
  /* Decoding the array of data */

  if (verbose) message("Reading the Data\n");
  if (st_read_find(s_length,file,"~A",&numline,string)) goto label_end;
  nech = ecr = nquant = 0;
  while (1)
  {
    if (st_read_next(s_length,file,1,&numline,string)) break;

    if (nech >= nquant * quantum) 
    {
      nquant++;
      tabloc = (double *)
        mem_realloc((char *) tabloc,sizeof(double) * nvar * nquant * quantum,1);
    }

    // Add the special variables
    tabloc[ecr++] = xwell;
    tabloc[ecr++] = ywell;
    tabloc[ecr++] = cwell;
    nvarlu = 3;
    
    // Add other variables

    lcur = string;
    while (nvarlu < nvar)
    {
      token = gslStrtok(lcur, sep_blank);
      lcur  = NULL;
      if (token == NULL || sscanf(token,"%lf",&value) == EOF) break;
      if (value == test) value = TEST;
      tabloc[ecr++] = value;
      nvarlu++;
    }

    if (nvarlu < nvar)
    {
        if (fgets(string,s_length,file) == NULL) break;
    }
    nech++;
  }
  
  /* Final resize */

  tabloc = (double *)
    mem_realloc((char *) tabloc,sizeof(double) * nvar * nech,1);
    
  /* Printout optional */
  
  if (verbose)
  {
    message("Number of variables = %d\n",nvar);
    message("Number of samples   = %d\n",nech);
  }

  /* Returning arguments */

  *nvarout   = nvar;
  *nechout   = nech;
  *var_names = varloc;
  *tab       = tabloc;
  
  /* Set the error return code */

  error = 0;

label_end:
  if (error)
  {
    varloc = (char  **) mem_free((char *) varloc);
    tabloc = (double *) mem_free((char *) tabloc);
  }
  if (file != nullptr) fclose(file);
  return(error);
}

/****************************************************************************/
/*!
**   Read the Grid from a F2G file (FLUMY format)
**
** \return  Error return code
**
** \param[in]  filename   Name of the F2G file
** \param[in]  verbose    1 for a verbose output; 0 otherwise
**
** \param[out]  nx        Array of number of grid nodes
** \param[out]  x0        Array of grid origin coordinates
** \param[out]  dx        Array of grid mesh
** \param[out]  angle     Rotation angle around Z-axis
** \param[out]  ncol      Number of fields that have been read (1 to 3)
** \param[out]  tab_arg   A pointer to the newly read data
**
** \remarks The returned array 'tab_arg' must be freed by the calling function
**
*****************************************************************************/
GEOSLIB_API int db_grid_read_f2g(const char *filename,
                                 int      verbose,
                                 int      nx[3],
                                 double   x0[3],
                                 double   dx[3],
                                 double  *angle,
                                 int     *ncol,
                                 double **tab_arg)
{ 
  FILE   *file;
  char    string[100],refchar[100],valtest[10],valread[10];
  double *tab,dum,value;
  int     error,size,ndim,version,nused;

  /* Initializations */

  error = 1;
  file  = nullptr;
  tab   = nullptr;
  nx[0] = nx[1] = nx[2] = 1;
  x0[0] = x0[1] = x0[2] = 0.;
  dx[0] = dx[1] = dx[2] = 1.;

  /* Open the file */

  file = gslFopen(filename,"r");
  if (file == nullptr)
  {
    messerr("Error when opening the F2G file %s for reading",filename);
    goto label_end;
  }

  /* Read the header */

  if (_record_read(file,"%s",string)) goto label_end;
  (void) gslStrcpy(refchar,"F2G_DIM");
  if (strcmp(string,refchar)) goto label_key;
  if (_record_read(file,"%d",&ndim)) goto label_contents;

  if (_record_read(file,"%s",string)) goto label_end;
  (void) gslStrcpy(refchar,"F2G_VERSION");
  if (strcmp(string,refchar)) goto label_key;
  if (_record_read(file,"%d",&version)) goto label_contents;

  if (_record_read(file,"%s",string)) goto label_end;
  (void) gslStrcpy(refchar,"F2G_LOCATION");
  if (strcmp(string,refchar)) goto label_key;
  for (int idim=0; idim<3; idim++) // Always three parameters
    if (_record_read(file,"%lf",&x0[idim])) goto label_contents;

  if (_record_read(file,"%s",string)) goto label_end;
  (void) gslStrcpy(refchar,"F2G_ROTATION");
  if (strcmp(string,refchar)) goto label_key;
  if (_record_read(file,"%lf",angle)) goto label_contents;

  if (_record_read(file,"%s",string)) goto label_end;
  (void) gslStrcpy(refchar,"F2G_ORIGIN");
  if (strcmp(string,refchar)) goto label_key;
  for (int idim=0; idim<ndim; idim++)
    if (_record_read(file,"%lf",&dum)) goto label_contents;

  if (_record_read(file,"%s",string)) goto label_end;
  (void) gslStrcpy(refchar,"F2G_NB_NODES");
  if (strcmp(string,refchar)) goto label_key;
  for (int idim=0; idim<ndim; idim++)
    if (_record_read(file,"%d",&nx[idim])) goto label_contents;

  if (_record_read(file,"%s",string)) goto label_end;
  (void) gslStrcpy(refchar,"F2G_LAGS");
  if (strcmp(string,refchar)) goto label_key;
  for (int idim=0; idim<ndim; idim++)
    if (_record_read(file,"%lf",&dx[idim])) goto label_contents;

  if (_record_read(file,"%s",string)) goto label_end;
  (void) gslStrcpy(refchar,"F2G_ORDER");
  if (strcmp(string,refchar)) goto label_key; 
  // We need to read the three next strings (orders) 
  // Only the order +Y +X +Z is interfaced
  if (_record_read(file,"%s",string)) goto label_end; 
  if (strcmp(string,"+Y")) goto label_end;
  if (_record_read(file,"%s",string)) goto label_end; 
  if (strcmp(string,"+X")) goto label_end;
  if (_record_read(file,"%s",string)) goto label_end; 
  if (strcmp(string,"+Z")) goto label_end;

  if (_record_read(file,"%s",string)) goto label_end;
  (void) gslStrcpy(refchar,"F2G_NB_VARIABLES");
  if (strcmp(string,refchar)) goto label_key;
  if (_record_read(file,"%d",ncol)) goto label_contents;

  for (int i=0; i<(*ncol); i++)
  {
    if (_record_read(file,"%s",string)) goto label_end; // Variable Name
    (void) gslSPrintf(refchar,"F2G_VARIABLE_%d",i+1);
    if (strcmp(string,refchar)) goto label_key;
    // We need to read the Name even if ignored
    if (_record_read(file,"%s",string)) goto label_end; 

    if (_record_read(file,"%s",string)) goto label_end; // NA value
    (void) gslSPrintf(refchar,"F2G_UNDEFINED_%d",i+1);
    if (strcmp(string,refchar)) goto label_key;
    if (_record_read(file,"%s",valtest)) goto label_contents;
  }

  if (_record_read(file,"%s",string)) goto label_end;
  (void) gslStrcpy(refchar,"F2G_VALUES");
  if (strcmp(string,refchar)) goto label_key;

  size = nx[0] * nx[1] * nx[2];
  tab  = (double *) mem_alloc(sizeof(double) * size,0);
  if (tab == nullptr) 
  {
    messerr("Core problem when reading the values of the F2G file");
    goto label_end;
  }
  for (int i=0; i<size; i++) tab[i] = TEST;

  nused = 0;
  for (int iz=0; iz<nx[2]; iz++)
    for (int ix=0; ix<nx[0]; ix++)
      for (int iy=0; iy<nx[1]; iy++)
        for (int icol=0; icol<(*ncol); icol++)
        {
          if (_record_read(file,"%s",valread)) goto label_short;
          if (! strcmp(valread,valtest))
            value = TEST;
          else
          {
            value = atof(valread);
            nused++;
          }
          F2G(ix,iy,iz) = value;
        }

  /* Patch the origin of the grid along vertical (to match RGeostats) */

  x0[2] = x0[2] - dx[2] / 2.;
    
  /* Verbose optional printout */

label_short:
  if (verbose) 
  {
    mestitle(0,"F2G File Characteristics :");
    st_print_verbose(3,nx,dx,x0);
    message("Number of valid values read = %d (out of %d)\n",nused,size);
  }

  /* Set the error return code */

  error = 0;
  goto label_end;

label_key:
  messerr("Error when looking for keyword %s",refchar);
  goto label_end;

label_contents:
  messerr("Error when looking at contents of keyword %s",refchar);
  goto label_end;

label_end:
  if (file != nullptr) fclose(file);
  *tab_arg = tab;
  return(error);
}

/****************************************************************************/
/*!
**   Read the Data frame from a CSV file. Reserved for numerical data frame.
**
** \return  Error return code
**
** \param[in]  filename    Name of the CSV file
** \param[in]  verbose     1 for a verbose output; 0 otherwise
** \param[in]  flag_header 1 if the first line of the file contains the
**                         variable names
** \param[in]  nskip       Number of lines to skip
** \param[in]  char_sep    Character used as a column separator
** \param[in]  char_dec    Character used as a decimal
** \param[in]  na_string   String used for absent information
** \param[in]  ncol_max    Maximum number of columns (or -1)
** \param[in]  nrow_max    Maximum number of rows (or -1)
**
** \param[out]  ncol_arg   Number of columns
** \param[out]  nrow_arg   Number of rows
** \param[out]  names      Array containing the variable names
** \param[out]  tab        Array of values
**
** \remarks The returned array 'tab' is organized by sample
**
*****************************************************************************/
GEOSLIB_API int csv_table_read(const String& filename,
                               int           verbose,
                               int           flag_header,
                               int           nskip,
                               char          char_sep,
                               char          char_dec,
                               const String& na_string,
                               int           ncol_max,
                               int           nrow_max,
                               int*          ncol_arg,
                               int*          nrow_arg,
                               VectorString& names,
                               VectorDouble& tab)
{
  std::string line;
  std::ifstream file(filename.c_str());
  if (!file.is_open())
  {
    messerr("Error when opening the CSV file %s for reading",filename.c_str());
    return 1;
  }

  // Remove windows stuff at the file beginning
  skipBOM(file);

  // Initialization
  names.clear();
  tab.clear();
  int ncol = 0;

  // Define the variable names
  if (flag_header)
  {
    std::getline(file, line);
    if (!line.empty())
    {
      line = trimRight(line);
      std::istringstream iss(line);
      std::string word;
      while (std::getline(iss, word, char_sep))
      {
        word = trim(word, "\"\'");
        word = trim(word);
        names.push_back(word);
        if (verbose) message("Column Name (%d): %s\n",ncol+1,word.c_str());
        ncol++;
        if (ncol_max > 0 && ncol >= ncol_max) break;
      }
    }

    if (verbose) message("Number of columns = %d\n",ncol);
  }

  // Skip some lines (optional)
  if (nskip > 0)
  {
    int iskip = 0;
    while (iskip<nskip && !file.eof())
    {
      std::getline(file, line);
      iskip++;
    }
  }

  // Read the values:
  int ncol2 = 0;
  int nrow = 0;
  while (!file.eof())
  {
    std::getline(file, line);
    if (!line.empty())
    {
      ncol2 = 0;
      std::istringstream iss(line);
      std::string word;
      while (std::getline(iss, word, char_sep))
      {
        if (word == na_string)
          tab.push_back(TEST);
        else
          tab.push_back(toDouble(word, char_dec));

        ncol2++;
        if (ncol_max > 0 && ncol2 >= ncol_max) break;
        if (ncol     > 0 && ncol2 >= ncol) break;
      }
      if (ncol <= 0) ncol = ncol2;
      nrow++;
    }
    if (nrow_max > 0 && nrow >= nrow_max) break;
  }

  // Optional printout
  if (verbose)
  {
    message("Data table read (%s) successfully\n",filename.c_str());
    message("- Number of columns = %d\n",ncol);
    message("- Number of rows    = %d\n",nrow);
  }

  *ncol_arg = ncol;
  *nrow_arg = nrow;

  return 0;
}

/****************************************************************************/
/*!
**   Write the Grid in a ENI_XYZ ASCII file
**
** \return  Error return code
**
** \param[in]  filename  Name of the ECLIPSE file
** \param[in]  db        Db structure to be written
** \param[in]  icol      Rank of the attribute
 **
 *****************************************************************************/
GEOSLIB_API int db_grid_write_XYZ(const char *filename, Db *db, int icol)
{
  FILE  *file;
  int    lec;

  /* Preliminary checks */

  if (! is_grid(db))
  {
    messerr("The Db structure should correspond to a Grid");
    return(1);
  }
  if (db->getNDim() != 2)
  {
    messerr("This FORMAT is limited to the 2-D case");
    return(1);
  }

  /* Open the file */

  file = gslFopen(filename,"w");
  if (file == nullptr)
  {
    messerr("Error when opening the XYZ file %s for writing",filename);
    return(1);
  }

  /* Write a comment */

  fprintf(file,"FDASCII 0 0 0 0 1E30\n");
  fprintf(file,"->\n");

  /* Write the set of values */

  lec = 0;
  for (int ix=0; ix<db->getNX(0); ix++)
    for (int iy=0; iy<db->getNX(1); iy++)
    {
      for (int i = 0; i < db->getNDim(); i++)
        fprintf(file, "%lf,", db->getCoordinate(lec, i));
      double value = db->getArray(lec, icol);
      if (FFFF(value))
        fprintf(file,"1E+30\n");
      else
        fprintf(file,"%lf\n",value);
      lec++;
    }

  if (file != nullptr) fclose(file);
  return(0);
}

/****************************************************************************/
/*!
**   Write a STRING element into the (opened) CSV file
**
** \param[in]  string       String to be written
**
** \remark: This function uses CSV_ENCODING static structure
** \remark: which must have been initiated beforehand
**
*****************************************************************************/
GEOSLIB_API void csv_print_string(const char *string)
{
  if (CSV_ENCODE == NULL)
    my_throw("You must initiate CSV_ENCODING first");

  (void) fprintf(CSV_ENCODE->file, "%s", string);
  if (CSV_ENCODE->current < CSV_ENCODE->nitem - 1)
  {
    (void) fprintf(CSV_ENCODE->file, "%s", CSV_ENCODE->char_sep);
    CSV_ENCODE->current++;
  }
  else
  {
    (void) fprintf(CSV_ENCODE->file, "\n");
    CSV_ENCODE->nlines++;
    CSV_ENCODE->current = 0;
  }
}

/****************************************************************************/
/*!
**   Write a DOUBLE element into the (opened) CSV file
**
** \param[in]  value        Real value to be written
**
** \remark: This function uses CSV_ENCODING static structure
** \remark: which must have been initiated beforehand
**
*****************************************************************************/
GEOSLIB_API void csv_print_double(double value)
{
  if (CSV_ENCODE == NULL)
    my_throw("You must initiate CSV_ENCODING first");

  if (FFFF(value))
    (void) fprintf(CSV_ENCODE->file, "%s", CSV_ENCODE->na_string);
  else
  {
    if (CSV_ENCODE->flag_integer)
      (void) fprintf(CSV_ENCODE->file, "%d", (int) value);
    else
      (void) fprintf(CSV_ENCODE->file, "%lf", value);
  }
  if (CSV_ENCODE->current < CSV_ENCODE->nitem - 1)
  {
    (void) fprintf(CSV_ENCODE->file, "%s", CSV_ENCODE->char_sep);
    CSV_ENCODE->current++;
  }
  else
  {
    (void) fprintf(CSV_ENCODE->file, "\n");
    CSV_ENCODE->nlines++;
    CSV_ENCODE->current = 0;
  }
}

/****************************************************************************/
/*!
**   Force the printing of End-Of-Line into the (opened) CSV file
**
** \remark: This function uses CSV_ENCODING static structure
** \remark: which must have been initiated beforehand
**
*****************************************************************************/
GEOSLIB_API void csv_print_eol(void)
{
  if (CSV_ENCODE->current <= 0) return;

  (void) fprintf(CSV_ENCODE->file, "\n");
  CSV_ENCODE->current = 0;
  CSV_ENCODE->nlines++;
}

/****************************************************************************/
/*!
**   Manage the Utility to write into a CSV file
**
** \return  Error return code
**
** \param[in]  filename     Name of the CSV file
** \param[in]  mode         1 for opening File; -1 for closing File
** \param[in]  nitem        Number of items per line
** \param[in]  flag_integer true if the numerical values must be printed as integer
** \param[in]  char_sep     Character used as a column separator
** \param[in]  na_string    String used for absent information
** \param[in]  verbose      Verbose flag
**
** \remark: This procedure manages an internal structure (declared as static)
** \remark: When opened, you can use csv_print_string() or csv_print_double()
** \remark: in order to store items in the file
** \remark: Do not forget to use csv_manage(-1,...) to close the file
**
*****************************************************************************/
GEOSLIB_API int csv_manage(const char *filename,
                           int mode,
                           int nitem,
                           bool flag_integer,
                           const char *char_sep,
                           const char *na_string,
                           bool verbose)
{
  // Dispatch

  if (mode > 0)
  {
    // Initiate the CSV_ENCODE structure

    if (CSV_ENCODE != NULL)
      CSV_ENCODE = (CSV_Encoding *) mem_free((char * ) CSV_ENCODE);
    CSV_ENCODE = (CSV_Encoding *) mem_alloc(sizeof(CSV_Encoding), 1);
    CSV_ENCODE->file = gslFopen(filename,"w");
    if (CSV_ENCODE->file == nullptr)
    {
      messerr("Error when opening the CSV file %s for writing",filename);
      (void) csv_manage(filename, -1, nitem, flag_integer, char_sep, na_string);
      return 1;
    }
    CSV_ENCODE->nitem        = nitem;
    CSV_ENCODE->current      = 0;
    CSV_ENCODE->nlines       = 0;
    CSV_ENCODE->flag_integer = flag_integer;
    CSV_ENCODE->char_sep     = char_sep;
    CSV_ENCODE->na_string    = na_string;

    // Optional printout

    if (verbose)
    {
      if (CSV_ENCODE->flag_integer)
        mestitle(1,"CSV Integer Encoding");
      else
        mestitle(1,"CSV Float Encoding\n");
      message("File Name                      = %s\n",filename);
      message("Number of items per line       = %d\n",CSV_ENCODE->nitem);
      message("Separator between items        = %s\n",CSV_ENCODE->char_sep);
      message("String for missing information = %s\n",CSV_ENCODE->na_string);
    }
  }
  else
  {
    // Write the last record (if necessary)
    csv_print_eol();

    if (CSV_ENCODE->file != NULL) fclose(CSV_ENCODE->file);

    // Option printout
    if (verbose)
    {
      if (CSV_ENCODE->flag_integer)
        message("CSV Integer Encoding : Summary\n");
      else
        message("CSV Float Encoding : Summary\n");
      message("Number of lines successfully written = %d\n",CSV_ENCODE->nlines);
    }

    if (CSV_ENCODE != NULL)
      CSV_ENCODE = (CSV_Encoding *) mem_free((char * ) CSV_ENCODE);
  }
  return 0;
}

/****************************************************************************/
/*!
**   Write the Data frame into a CSV file. Reserved for numerical data frame.
**
** \return  Error return code
**
** \param[in]  db           Name of the Db
** \param[in]  filename     Name of the CSV file
** \param[in]  flag_header  1 if the variable names must be dumed out
** \param[in]  flag_allcol  1 if all the columns available must be dumped out
** \param[in]  flag_coor    1 if the coordinates must be dumped out
** \param[in]  flag_integer true if the numerical values must be printed as integer
** \param[in]  char_sep     Character used as a column separator
** \param[in]  na_string    String used for absent information
**
** \remarks: This procedure dumps the Z-variables and optionally the X-variables
**
*****************************************************************************/
GEOSLIB_API int db_write_csv(Db *db,
                             const char *filename,
                             int flag_header,
                             int flag_allcol,
                             int flag_coor,
                             bool flag_integer,
                             const char *char_sep,
                             const char *na_string)
{
  if (db == nullptr) return 1;
  int ncol = db->getFieldNumber();
  int ndim = db->getNDim();
  int nech = db->getSampleNumber();
  int nvar = db->getVariableNumber();

  // Count the number of items per line

  int nitem = 0;
  if (flag_allcol)
    nitem = ncol;
  else
  {
    nitem = nvar;
    if (flag_coor)
      nitem += ndim;
  }

  // Initiate the CSV_Encoding structure

  if (csv_manage(filename, 1, nitem, flag_integer, char_sep, na_string))
    return 1;

  /* Dump the header */

  if (flag_header)
  {
    // Case where all columns are dumped out

     if (flag_allcol)
     {
       for (int rank = 0; rank < ncol; rank++)
       {
         csv_print_string(db_name_get_by_att(db,rank).c_str());
       }
     }
     else
     {
       int rank = 0;
       if (flag_coor)
         for (int idim = 0; idim < ndim; idim++)
         {
           int iatt = db_attribute_identify(db,ELoc::X,idim);
           csv_print_string(db_name_get_by_att(db,iatt).c_str());
           rank++;
         }
       for (int ivar = 0; ivar < nvar; ivar++)
       {
         int iatt = db_attribute_identify(db,ELoc::Z,ivar);
         csv_print_string(db_name_get_by_att(db,iatt).c_str());
         rank++;
       }
     }
  }

  // Dump the samples (one sample per line)

  for (int iech = 0; iech < nech; iech++)
  {
    if (!db->isActive(iech)) continue;

    if (flag_allcol)
    {
      for (int rank = 0; rank < ncol; rank++)
        csv_print_double(db->getByColumn(iech, rank));
    }
    else
    {
      int rank = 0;
      if (flag_coor)
        for (int idim = 0; idim < ndim; idim++)
        {
          int iatt = db_attribute_identify(db, ELoc::X, idim);
          csv_print_double(db->getCoordinate(iech, iatt));
          rank++;
        }
      for (int ivar = 0; ivar < nvar; ivar++)
      {
        int iatt = db_attribute_identify(db, ELoc::Z, ivar);
        csv_print_double(db->getVariable(iech, iatt));
        rank++;
      }
    }
  }

  // Close the file
  (void) csv_manage(filename, -1, nitem, flag_integer, char_sep, na_string);

  return 0;
}
