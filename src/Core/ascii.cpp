/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_e.h"
#include "geoslib_enum.h"
#include "Anamorphosis/Anam.hpp"
#include "Anamorphosis/AnamDiscreteDD.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Anamorphosis/AnamEmpirical.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Basic/Utilities.hpp"
#include "Covariances/CovAniso.hpp"

/*! \cond */
#define OLD 0
#define NEW 1

#define NODES(inode,i)   (nodes[6 * (inode) + (i)])
#define FROM_TYPE(inode) (nodes[6 * (inode) + 0])
#define FROM_RANK(inode) (nodes[6 * (inode) + 1])
#define FROM_VERS(inode) (nodes[6 * (inode) + 2])
#define NODE_TYPE(inode) (nodes[6 * (inode) + 3])
#define NODE_RANK(inode) (nodes[6 * (inode) + 4])
#define FACIES(inode)    (nodes[6 * (inode) + 5])

static int ASCII_BUFFER_LENGTH = 0;
static int ASCII_BUFFER_QUANT = 1000;
static char *ASCII_BUFFER = NULL;
static FILE *FILE_MEM = NULL;
static char FILE_NAME_MEM[BUFFER_LENGTH];

/*! \endcond */

static char STUDY[BUFFER_LENGTH] = "./";
static char EXT_DAT[] = "dat";
static char EXT_OUT[] = "out";
static char Fichier_environ[] = "Environ";
static char Fichier_donnees[] = "Data";
static char Fichier_grid[] = "Grid";
static char Fichier_vario[] = "Vario";
static char Fichier_model[] = "Model";
static char Fichier_neigh[] = "Neigh";
static char Fichier_polygon[] = "Polygon";
static char Fichier_option[] = "Option";
static char Fichier_rule[] = "Rule";
static char Fichier_simu[] = "Simu";
static char Fichier_frac[] = "Frac";

/****************************************************************************/
/*!
 **  Allocate the buffer
 **
 *****************************************************************************/
static void st_buffer_initiate(void)

{
  ASCII_BUFFER_LENGTH = ASCII_BUFFER_QUANT;
  ASCII_BUFFER = mem_realloc(ASCII_BUFFER, ASCII_BUFFER_QUANT, 1);
  ASCII_BUFFER[0] = '\0';
  _erase_current_string();
}

/****************************************************************************/
/*!
 **  Close the buffer
 **
 ** \param[out]  buffer     Output buffer
 ** \param[out]  buf_length Buffer length
 **
 *****************************************************************************/
static void st_buffer_close(char **buffer, int *buf_length)
{
  int long1;

  long1 = 0;
  if (ASCII_BUFFER != NULL)
  {
    long1 = strlen(ASCII_BUFFER);
    if (long1 <= 0)
      ASCII_BUFFER = mem_free(ASCII_BUFFER);
    else
      ASCII_BUFFER = mem_realloc(ASCII_BUFFER, long1, 1);
  }
  *buf_length = long1;
  *buffer = ASCII_BUFFER;
}

/****************************************************************************/
/*!
 **  Read the next record
 **
 ** \return Error return code
 **
 ** \param[in]  title      Name of the quantity to be read
 ** \param[in]  format     Encoding format
 ** \param[in]  ...        Value to be written
 **
 *****************************************************************************/
static int st_record_read(const char *title, const char *format, ...)
{
  va_list ap;
  int error;

  error = 0;
  va_start(ap, format);

  if (FILE_MEM != (FILE *) NULL)
  {
    error = _file_read(FILE_MEM, format, ap);
  }
  else
  {
    error = _buffer_read(&ASCII_BUFFER, format, ap);
  }

  if (error > 0)
  {
    messerr("Error when reading '%s' from %s", title, FILE_NAME_MEM);
    print_current_line();
  }

  va_end(ap);
  return (error);
}

/****************************************************************************/
/*!
 **  Write the next record
 **
 ** \param[in]  format     Encoding format
 ** \param[in]  ...        Value to be written
 **
 *****************************************************************************/
static void st_record_write(const char *format, ...)
{
  va_list ap;
  char buf[1000];
  int long1, long2;

  va_start(ap, format);
  if (FILE_MEM != (FILE *) NULL)
  {
    _file_write(FILE_MEM, format, ap);
  }
  else
  {
    _buffer_write(buf, format, ap);
    long1 = strlen(buf);
    long2 = (ASCII_BUFFER != NULL) ? strlen(ASCII_BUFFER) :
                                     0;
    while (long1 + long2 > ASCII_BUFFER_LENGTH)
    {
      ASCII_BUFFER_LENGTH += ASCII_BUFFER_QUANT;
      ASCII_BUFFER = mem_realloc(ASCII_BUFFER, ASCII_BUFFER_LENGTH, 1);
    }
    (void) strcat(ASCII_BUFFER, buf);
  }

  va_end(ap);
}

/****************************************************************************/
/*!
 **   Create the File Name by patching the generic name
 **
 ** \param[in]  ref_name Reference Name
 ** \param[in]  rank     Rank of the name
 ** \param[in]  mode     Choice of the added extension
 ** \li                  0 for reading - extension ".dat"
 ** \li                  1 for writing - extension ".out"
 ** \li                 -1 no extension
 **
 ** \param[out] file_name Output filename
 **
 ** \remark  When the rank is 0, the generic name is returned
 ** \remark  Otherwise the rank is combined in the name
 **
 *****************************************************************************/
static void st_filename_patch(const char *ref_name,
                              int rank,
                              int mode,
                              char *file_name)
{
  if (rank == 0)
  {
    switch (mode)
    {
      case 0:
        (void) sprintf(file_name, "%s/%s.%s", STUDY, ref_name, EXT_DAT);
        break;

      case 1:
        (void) sprintf(file_name, "%s/%s.%s", STUDY, ref_name, EXT_OUT);
        break;

      case -1:
        (void) sprintf(file_name, "%s/%s", STUDY, ref_name);
        break;
    }
  }
  else
  {
    switch (mode)
    {
      case 0:
        (void) sprintf(file_name, "%s/%s%1d.%s", STUDY, ref_name, rank,
                       EXT_DAT);
        break;

      case 1:
        (void) sprintf(file_name, "%s/%s%1d.%s", STUDY, ref_name, rank,
                       EXT_OUT);
        break;

      case -1:
        (void) sprintf(file_name, "%s/%s%1d", STUDY, ref_name, rank);
        break;
    }
  }
  return;
}

/****************************************************************************/
/*!
 **   Returns the name of the external file
 **
 ** \param[in]  filein    External filename
 ** \param[in]  mode      0 for read; 1 for write
 **
 ** \param[out] filename  Output filename
 **
 *****************************************************************************/
GEOSLIB_API void ascii_external_filename(const char *filein,
                                         int mode,
                                         char *filename)
{
  st_filename_patch(filein, 0, mode, filename);
}

/****************************************************************************/
/*!
 **   Returns the name of the file
 **
 ** \param[in]  type      Type of the file to be named
 ** \param[in]  rank      Rank of the file (optional)
 ** \param[in]  mode      0 for read; 1 for write
 **
 ** \param[out] filename  Output filename
 **
 *****************************************************************************/
GEOSLIB_API void ascii_filename(const char *type,
                                int rank,
                                int mode,
                                char *filename)
{
  if (!strcmp(type, "Environ"))
    st_filename_patch(Fichier_environ, rank, mode, filename);
  else if (!strcmp(type, "Data"))
    st_filename_patch(Fichier_donnees, rank, mode, filename);
  else if (!strcmp(type, "Grid"))
    st_filename_patch(Fichier_grid, rank, mode, filename);
  else if (!strcmp(type, "Vario"))
    st_filename_patch(Fichier_vario, rank, mode, filename);
  else if (!strcmp(type, "Model"))
    st_filename_patch(Fichier_model, rank, mode, filename);
  else if (!strcmp(type, "Neigh"))
    st_filename_patch(Fichier_neigh, rank, mode, filename);
  else if (!strcmp(type, "Rule"))
    st_filename_patch(Fichier_rule, rank, mode, filename);
  else if (!strcmp(type, "Simu"))
    st_filename_patch(Fichier_simu, rank, mode, filename);
  else if (!strcmp(type, "Polygon"))
    st_filename_patch(Fichier_polygon, rank, mode, filename);
  else if (!strcmp(type, "Option"))
    st_filename_patch(Fichier_option, rank, mode, filename);
  else if (!strcmp(type, "Frac"))
    st_filename_patch(Fichier_frac, rank, mode, filename);
  else
  {
    messageAbort("The file type %s is not referenced", type);
  }
}

/****************************************************************************/
/*!
 **   Set the Study for the test data
 **
 ** \param[in]  study Local name of the study
 **
 *****************************************************************************/
GEOSLIB_API void ascii_study_define(const char *study)

{
  (void) strcpy(STUDY, study);
  return;
}

/****************************************************************************/
/*!
 **   Close the ASCII file
 **
 ** \param[in]  file       FILE structure to be close
 **
 *****************************************************************************/
static void st_file_close(FILE *file)
{
  FILE_MEM = NULL;
  fclose(file);
}

/****************************************************************************/
/*!
 **   Open an ASCII file
 **
 ** \return  FILE returned pointer
 **
 ** \param[in]  filename Local file name
 ** \param[in]  filetype Type of the file (optional [NULL] when NEW)
 ** \param[in]  mode     type of file (OLD or NEW)
 ** \param[in]  verbose  Verbose option if the file cannot be opened
 **
 *****************************************************************************/
static FILE *st_file_open(const char *filename,
                          const char *filetype,
                          int mode,
                          int verbose)
{
  FILE *file;
  char idtype[LONG_SIZE];

  /* Open the file */

  file = FILE_MEM = _file_open(filename, mode);
  (void) strcpy(FILE_NAME_MEM, filename);

  if (file == (FILE *) NULL)
  {
    if (verbose) messerr("Error when opening the file %s", filename);
    FILE_MEM = NULL;
    return (file);
  }

  if (debug_query("interface")) message("Opening the File = %s\n", filename);

  /* Check against the file type */

  if (mode == OLD)
  {
    if (st_record_read("File Type", "%s", idtype))
    {
      FILE_MEM = NULL;
      return (NULL);
    }
    if (strcmp(idtype, filetype))
    {
      messerr(
          "Error: in the File (%s), its Type (%s) does not match the requested one (%s)",
          filename, idtype, filetype);
      FILE_MEM = NULL;
      return (NULL);
    }
  }
  else
  {
    if (filetype != NULL)
    {
      st_record_write("%s", filetype);
      st_record_write("\n");
    }
  }

  return (file);
}

/****************************************************************************/
/*!
 **   Deliver a message for successfull creation
 **
 ** \param[in]  filename Local file name
 ** \param[in]  filetype Type of the file
 ** \param[in]  verbose  Verbose option if the file cannot be opened
 **
 *****************************************************************************/
static void st_create_message(const char *filename,
                              const char *filetype,
                              int verbose)
{
  if (!verbose) return;

  if (filename != NULL)
    message("File %s (%s) created successfully\n", filename, filetype);
  else
    message("Buffer created successfully\n", filetype);

  return;
}

/****************************************************************************/
/*!
 **  Write the table (real values)
 **
 ** \param[in]  string     String of the table records (or NULL)
 ** \param[in]  ntab       Number of values in the table
 ** \param[in]  tab        Array of real values to be written
 **
 *****************************************************************************/
static void st_table_write(const char *string, int ntab, const double *tab)
{
  int i;
  char local[LONG_SIZE];

  for (i = 0; i < ntab; i++)
  {
    st_record_write("%lf", tab[i]);
    if (string != NULL)
    {
      (void) sprintf(local, "%s (%d)", string, i + 1);
      st_record_write("#", local);
    }
    else
    {
      st_record_write("\n");
    }
  }
}

/****************************************************************************/
/*!
 **  Write the table (integer values)
 **
 ** \param[in]  string     String of the table records (or NULL)
 ** \param[in]  ntab       Number of values in the table
 ** \param[in]  itab       Array of integer values to be written
 **
 *****************************************************************************/
static void st_tablei_write(const char *string, int ntab, int *itab)
{
  int i;
  char local[LONG_SIZE];

  for (i = 0; i < ntab; i++)
  {
    st_record_write("%d", itab[i]);
    if (string != NULL)
    {
      (void) sprintf(local, "%s (%d)", string, i + 1);
      st_record_write("#", local);
    }
    else
    {
      st_record_write("\n");
    }
  }
}

/****************************************************************************/
/*!
 **  Read the table
 **
 ** \param[in]  ntab       Number of values in the table
 **
 ** \param[out]  tab        Array to be written
 **
 ** \remark The receiving array 'tab' must be dimensioned beforehand
 **
 *****************************************************************************/
static int st_table_read(int ntab, double *tab)
{
  int i;

  for (i = 0; i < ntab; i++)
    if (st_record_read("Reading Table", "%lf", &tab[i])) return (1);
  return (0);
}

/****************************************************************************/
/*!
 **   Read variables from an ASCII file
 **
 ** \return  Error return code
 **
 ** \param[in]  file       FILE structure
 **
 ** \param[out]  natt_r    Number of attributes
 ** \param[out]  ndim_r    Dimension of the space
 ** \param[out]  nech_r    Number of samples
 ** \param[out]  tabatt    Array containing the attribute names
 ** \param[out]  tabnum    Array containing the attribute ranks
 ** \param[out]  tabnam    Array containing the attribute names
 ** \param[out]  tab       Array of real values
 **
 *****************************************************************************/
static void st_variables_read(FILE *file,
                              int *natt_r,
                              int *ndim_r,
                              int *nech_r,
                              std::vector<ENUM_LOCS>& tabatt,
                              VectorInt& tabnum,
                              VectorString& tabnam,
                              VectorDouble& tab)
{
  char line[LONG_SIZE];
  int  inum, natt, ndim, nval, ecr, mult;
  ENUM_LOCS iatt;
  double value;

  /* Initializations */

  natt = nval = ndim = 0;

  /* Read the number of variables */

  if (st_record_read("Number of Variables", "%d", &natt)) goto label_end;

  /* Decoding the locators */

  ecr = 0;
  while (1)
  {
    if (ecr >= natt) break;
    if (st_record_read("Locator Name", "%s", line)) goto label_end;
    if (locatorIdentify(line, &iatt, &inum, &mult)) break;
    tabatt.push_back(iatt);
    tabnum.push_back(inum);
    if (iatt == LOC_X) ndim++;
    ecr++;
  }

  /* Decoding the names */

  ecr = 0;
  while (1)
  {
    if (ecr >= natt) break;
    if (st_record_read("Variable Name", "%s", line)) goto label_end;
    tabnam.push_back(line);
    ecr++;
  }

  /* Read the numeric values */

  while (1)
  {
    if (st_record_read("Numerical value", "%lf", &value)) goto label_end;
    tab.push_back(value);
    nval++;
  }

  label_end:

  /* Returning arguments */

  *natt_r = natt;
  *nech_r = (natt > 0) ? nval / natt :
                         0;
  *ndim_r = ndim;
  return;
}

/****************************************************************************/
/*!
 **   Write variables in an ASCII file
 **
 ** \return  Error return code
 **
 ** \param[in]  file         FILE structure
 ** \param[in]  db           Pointer to the Db structure
 ** \param[in]  flag_grid    1 if the Db should be printed as a grid file
 ** \param[in]  flag_calcul  1 if the data must be printed
 **
 *****************************************************************************/
static int st_variables_write(FILE *file,
                              Db *db,
                              int flag_grid,
                              int flag_calcul)
{
  int icol, iech, item, ecr, ncol;
  ENUM_LOCS locatorType;

  /* Preliminary check */

  if (db->getFieldNumber() <= 0 || get_NECH(db) <= 0) return (0);
  if (flag_grid && !flag_calcul) return (0);

  /* Count the number of variables to be written */

  ncol = 0;
  for (icol = 0; icol < db->getFieldNumber(); icol++)
  {
    if (! db->getLocatorByColumn(icol, &locatorType, &item)) continue;
    if (flag_grid && locatorType == LOC_X) continue;
    ncol++;
  }
  st_record_write("%d", ncol);
  st_record_write("#", "Number of variables");

  /* Print the locators */

  st_record_write("#", "Locators");
  for (icol = ecr = 0; icol < db->getFieldNumber(); icol++)
  {
    if (! db->getLocatorByColumn(icol, &locatorType, &item)) continue;
    if (flag_grid && locatorType == LOC_X) continue;
    if (ecr >= ncol) break;
    String string = getLocatorName(locatorType, item);
    st_record_write("%s", string.c_str());
    ecr++;
  }
  st_record_write("\n");

  /* Print the variable names */

  st_record_write("#", "Names");
  for (icol = ecr = 0; icol < db->getFieldNumber(); icol++)
  {
    if (! db->getLocatorByColumn(icol, &locatorType, &item)) continue;
    if (flag_grid && locatorType == LOC_X) continue;
    if (ecr >= ncol) break;
    st_record_write("%s", db_name_get_by_att(db, icol).data());
  }
  st_record_write("\n");

  /* Print the array of values */

  st_record_write("#", "Array of values");
  for (iech = 0; iech < get_NECH(db); iech++)
  {
    if (!flag_grid && !db->getSelection(iech)) continue;
    for (icol = 0; icol < db->getFieldNumber(); icol++)
    {
      if (! db->getLocatorByColumn(icol, &locatorType, &item)) continue;
      if (flag_grid && locatorType == LOC_X) continue;
      st_record_write("%lf", get_ARRAY(db, iech, icol));
    }
    st_record_write("\n");
  }
  return (0);
}

/****************************************************************************/
/*!
 **   Read the Environment definition file
 **
 ** \param[in] file_name  Name of the ASCII file
 ** \param[in] verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
GEOSLIB_API void ascii_environ_read(char *file_name, int verbose)

{
  FILE *file;
  char name[10];
  int debug;

  /* Opening the Data file */

  file = st_file_open(file_name, "Environ", OLD, verbose);
  if (file == (FILE *) NULL) return;

  /* Reading the environment */

  while (1)
  {
    if (st_record_read("Debug Keyword", "%s", name)) goto label_end;
    if (st_record_read("Debug Value", "%d", &debug)) goto label_end;
    debug_define(name, debug);
  }

  label_end: st_file_close(file);
  return;
}

/****************************************************************************/
/*!
 **   Read a Db
 **
 ** \return  Pointer to the Db descriptor
 **
 ** \param[in]  file_name  Name of the ASCII file
 ** \param[in]  must_grid  Desired grid type
 ** \li                     0 : Set of isolated points
 ** \li                     1 : Regular grid
 ** \li                    -1 : Any organization
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
GEOSLIB_API Db *ascii_db_read(const char *file_name, int must_grid, int verbose)
{
  Db *db;
  FILE *file;
  int ndim, ndim2, idim, ntot, natt, nech, i, flag_grid;
  VectorInt tabnum;
  std::vector<ENUM_LOCS> tabatt;
  VectorInt nx;
  VectorString tabnam;
  VectorDouble x0;
  VectorDouble dx;
  VectorDouble tab;
  static int flag_add_rank = 1;

  /* Initializations */

  db = (Db *) NULL;
  natt = ndim = nech = ntot = 0;

  /* Opening the Data file */

  file = st_file_open(file_name, "Db", OLD, verbose);
  if (file == (FILE *) NULL) return (db);

  /* Check the grid organization */

  if (st_record_read("Grid Flag", "%d", &flag_grid)) goto label_end;
  if (must_grid >= 0 && must_grid != flag_grid)
  {
    messerr(
        "The ASCII file does not have the requested organization (Point or Grid)");
    goto label_end;
  }

  /* Grid case; read the grid header */

  if (flag_grid)
  {

    /* Decoding the header */

    if (st_record_read("Space Dimension", "%d", &ndim)) goto label_end;

    /* Core allocation */

    nx.resize(ndim);
    dx.resize(ndim);
    x0.resize(ndim);

    /* Read the grid characteristics */

    ntot = 1;
    for (idim = 0; idim < ndim; idim++)
    {
      if (st_record_read("Grid Number of Nodes", "%d", &nx[idim]))
        goto label_end;
      if (st_record_read("Grid Origin", "%lf", &x0[idim])) goto label_end;
      if (st_record_read("Grid Mesh", "%lf", &dx[idim])) goto label_end;
      ntot *= nx[idim];
    }
  }

  /* Reading the tail of the file */

  st_variables_read(file, &natt, &ndim2, &nech, tabatt, tabnum, tabnam, tab);

  /* Creating the Db */

  if (flag_grid)
  {
    if (natt > 0 && nech != ntot)
    {
      messerr("The number of lines read from the Grid file (%d)", nech);
      messerr("is not a multiple of the number of samples (%d)", ntot);
      messerr("The Grid Db is created with no sample attached");
      natt = 0;
    }
    db = db_create_grid(0, ndim, natt, 1, flag_add_rank,
                        nx, x0, dx, VectorDouble(), tab);
    if (db == (Db *) NULL) goto label_end;
  }
  else
  {
    db = db_create_point(nech, natt, 1, flag_add_rank, tab);
    if (db == (Db *) NULL) goto label_end;
  }

  /* Loading the names */

  if (natt > 0)
    for (i = 0; i < natt; i++)
  {
    if (db_name_set(db, i+flag_add_rank, tabnam[i]))
      messerr("Error in db_name_set");
  }

  /* Create the locators */

  if (natt > 0)
    for (i = 0; i < natt; i++)
      db->setLocatorByAttribute(i+flag_add_rank, tabatt[i], tabnum[i]);

  if (debug_query("interface") && debug_query("db")) db_print(db, 1);

  /* Core deallocation */

  label_end: st_file_close(file);
  return (db);
}

/****************************************************************************/
/*!
 **   Read the variogram definition file
 **
 ** \return  Pointer to the Vario structure
 **
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
static Vario *st_ascii_vario_read(int verbose)
{
  Vario *vario;
  int idim, ndim, ivar, jvar, nvar, ndir, idir, npas, opt_code;
  int flag_regular, flag_calcul, i, ecr;
  double dpas, tolang, scale, tolcode, toldis;
  VectorDouble codir, grincr, vars;

  /* Initializations */

  vario = (Vario *) NULL;

  /* Create the Vario structure */

  if (st_record_read("Space Dimension", "%d", &ndim)) goto label_end;
  if (st_record_read("Number of Variables", "%d", &nvar)) goto label_end;
  if (st_record_read("Number of Variogram Directions", "%d", &ndir))
    goto label_end;
  if (st_record_read("Scale", "%lf", &scale)) goto label_end;

  /* Read the variances (optional) */

  if (st_record_read("Variogram calculation Option", "%d", &flag_calcul))
    goto label_end;
  vars.resize(nvar * nvar);
  if (flag_calcul)
    for (ivar = ecr = 0; ivar < nvar; ivar++)
      for (jvar = 0; jvar < nvar; jvar++, ecr++)
        if (st_record_read("Experimental Variance term", "%lf",
                           &vars[ecr])) goto label_end;

  /* Initialize the variogram structure */

  vario = variogram_init("vg", scale, VectorDouble());
  if (vario == (Vario *) NULL) return (vario);
  vario->resize(ndim, nvar);
  vario->setVars(vars);

  /* Reading the variogram calculation directions */

  for (idir = 0; idir < ndir; idir++)
  {
    if (st_record_read("Regular Variogram Calculation", "%d", &flag_regular))
      goto label_end;
    if (st_record_read("Number of Variogram Lags", "%d", &npas)) goto label_end;
    if (st_record_read("Variogram Code Option", "%d", &opt_code))
      goto label_end;
    if (st_record_read("Tolerance on Code", "%lf", &tolcode)) goto label_end;
    if (st_record_read("Lag Value", "%lf", &dpas)) goto label_end;
    if (st_record_read("Tolerance on Distance", "%lf", &toldis)) goto label_end;
    if (st_record_read("Tolerance on Direction", "%lf", &tolang))
      goto label_end;
    codir.resize(ndim);
    grincr.resize(ndim);
    for (idim = 0; idim < ndim; idim++)
      if (st_record_read("Direction vector", "%lf", &codir[idim]))
        goto label_end;
    for (idim = 0; idim < ndim; idim++)
      if (st_record_read("Grid Increment", "%lf", &grincr[idim]))
        goto label_end;

    if (variogram_direction_add(vario, npas, opt_code, 0, dpas,
                                toldis, tolang, TEST, TEST, tolcode,
                                VectorDouble(), codir, grincr))
      goto label_end;

    /* Read the arrays of results (optional) */

    if (!flag_calcul) continue;
    Dir& dir = vario->getDirs(idir);
    dir.resize(nvar, vario->getFlagAsym());
    for (i = 0; i < dir.getSize(); i++)
    {
      double sw, hh, gg;
      if (st_record_read("Experimental Variogram Weight", "%lf", &sw))
        goto label_end;
      dir.setSw(i,sw);
      if (st_record_read("Experimental Variogram Distance", "%lf", &hh))
        goto label_end;
      dir.setHh(i,hh);
      if (st_record_read("Experimental Variogram Value", "%lf", &gg))
        goto label_end;
      dir.setGg(i,gg);
    }
  }

  label_end: if (debug_query("interface") && debug_query("model"))
    variogram_print(vario, 1);
  return (vario);
}

/****************************************************************************/
/*!
 **   Read the variogram definition file
 **
 ** \return  Pointer to the Vario structure
 **
 ** \param[in]  file_name  Name of the ASCII file
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
GEOSLIB_API Vario *ascii_vario_read(const char *file_name, int verbose)
{
  Vario *vario;
  FILE *file;

  /* Opening the Data file */

  vario = (Vario *) NULL;
  file = st_file_open(file_name, "Vario", OLD, verbose);
  if (file == (FILE *) NULL) return (vario);

  /* Read the variogram */

  vario = st_ascii_vario_read(verbose);

  /* Close the file */

  st_file_close(file);

  return (vario);
}

/****************************************************************************/
/*!
 **   Read the variogram from a buffer
 **
 ** \return  Pointer to the Vario structure
 **
 ** \param[in]  buffer     Buffer to be read
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
GEOSLIB_API Vario *ascii_vario_read_buffer(char *buffer, int verbose)
{
  Vario *vario;

  /* Attaching the buffer */

  FILE_MEM = NULL;
  ASCII_BUFFER = buffer;
  _erase_current_string();

  /* Read the variogram */

  vario = st_ascii_vario_read(verbose);

  /* Detach the buffer */

  ASCII_BUFFER = NULL;

  return (vario);
}

/****************************************************************************/
/*!
 **   Read the model definition file
 **
 ** \return  Pointer to the Model structure
 **
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
static Model *st_ascii_model_read(int verbose)

{
  Model *model;
  int ndim, nvar, type, ncova, nbfl, icova, jvar, ibfl, ivar, idim, jdim;
  int flag_aniso, flag_rotation, lec;
  double field, range, value, param, radius;
  VectorDouble aniso_ranges, aniso_rotmat;

  /* Initializations */

  model = (Model *) NULL;
  ndim = nvar = ncova = nbfl = 0;
  field = radius = 0.;

  /* Create the Model structure */

  if (st_record_read("Space Dimension", "%d", &ndim)) goto label_end;
  if (st_record_read("Number of Variables", "%d", &nvar)) goto label_end;
  if (st_record_read("Field dimension", "%lf", &field)) goto label_end;
  if (st_record_read("Radius for Model", "%lf", &radius)) goto label_end;
  if (st_record_read("Number of Basic Structures", "%d", &ncova))
    goto label_end;
  if (st_record_read("Number of Basic Drift Functions", "%d", &nbfl))
    goto label_end;

  //TODO Move this horrible setting which is used only because the Space
  // is not known beforehand
  if (! ASpaceObject::hasGlobalSpace())
    ASpaceObject::createGlobalSpace(SPACE_RN,ndim);
  model = model_init(ndim, nvar, field, 0, radius, false);
  if (model == (Model *) NULL) return (model);

  /* Reading the covariance part */

  for (icova = 0; icova < ncova; icova++)
  {
    flag_aniso = flag_rotation = 0;

    if (st_record_read("Covariance Type", "%d", &type)) goto label_end;
    if (st_record_read("Isotropic Range", "%lf", &range)) goto label_end;
    if (st_record_read("Model third Parameter", "%lf", &param)) goto label_end;

    if (st_record_read("Flag for Anisotropy", "%d", &flag_aniso))
      goto label_end;
    if (flag_aniso)
    {
      aniso_ranges.resize(ndim);
      // In fact, the file contains the anisotropy coefficients
      // After reading, we must turn them into anisotropic ranges
      for (idim = 0; idim < ndim; idim++)
        if (st_record_read("Anisotropy coefficient", "%lf",
                           &aniso_ranges[idim])) goto label_end;
      for (idim = 0; idim < ndim; idim++) aniso_ranges[idim] *= range;

      if (st_record_read("Flag for Anisotropy Rotation", "%d", &flag_rotation))
        goto label_end;
      if (flag_rotation)
      {
        // Warning: the storage in the File is performed by Column
        // whereas the internal storage (Cova) is by column
        aniso_rotmat.resize(ndim * ndim);
        lec = 0;
        for (idim = 0; idim < ndim; idim++)
          for (jdim = 0; jdim < ndim; jdim++)
            if (st_record_read("Anisotropy Rotation Matrix", "%lf",
                               &aniso_rotmat[lec++])) goto label_end;
      }
    }
    if (model_add_cova(model, type, flag_aniso, flag_rotation, range, param,
                       aniso_ranges, aniso_rotmat, VectorDouble()))
      goto label_end;
  }

  /* Reading the drift part */

  for (ibfl = 0; ibfl < nbfl; ibfl++)
  {
    if (st_record_read("Drift Function", "%d", &type)) goto label_end;
    if (model_add_drift(model, type, 0)) goto label_end;
  }

  /* Reading the matrix of means (only if nbfl <= 0) */

  if (nbfl <= 0)
    for (ivar = 0; ivar < nvar; ivar++)
    {
      double mean;
      if (st_record_read("Mean of Variable", "%lf", &mean))
        goto label_end;
      model->setMean(ivar,mean);
    }

  /* Reading the matrices of sills (optional) */

  for (icova = 0; icova < ncova; icova++)
  {
    for (ivar = 0; ivar < nvar; ivar++)
      for (jvar = 0; jvar < nvar; jvar++)
      {
        if (st_record_read("Matrix of Sills", "%lf", &value)) goto label_end;
        model->setSill(icova,ivar,jvar,value);
      }
  }

  /* Reading the variance-covariance at the origin (optional) */

  for (ivar = 0; ivar < nvar; ivar++)
    for (jvar = 0; jvar < nvar; jvar++)
    {
      if (st_record_read("Variance-covariance at Origin", "%lf", &value))
        goto label_end;
      model->setCovar0(ivar, jvar, value);
    }

  label_end: (void) model_setup(model);
  if (debug_query("interface") && debug_query("model")) model->display();
  return (model);
}

/****************************************************************************/
/*!
 **   Read the model definition file from an ASCII file
 **
 ** \return  Pointer to the Model structure
 **
 ** \param[in]  file_name  Name of the ASCII file
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
GEOSLIB_API Model *ascii_model_read(const char *file_name, int verbose)

{
  Model *model;
  FILE *file;

  /* Opening the Data file */

  model = (Model *) NULL;
  file = st_file_open(file_name, "Model", OLD, verbose);
  if (file == (FILE *) NULL) return (model);

  /* Read the Model structure */

  model = st_ascii_model_read(verbose);

  /* Close the file */

  st_file_close(file);

  return (model);
}

/****************************************************************************/
/*!
 **   Read the model from a buffer
 **
 ** \return  Pointer to the Model structure
 **
 ** \param[in]  buffer     Buffer to be read
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
GEOSLIB_API Model *ascii_model_read_buffer(char *buffer, int verbose)
{
  Model *model;

  /* Attaching the buffer */

  FILE_MEM = NULL;
  ASCII_BUFFER = buffer;
  _erase_current_string();

  /* Read the Model structure */

  model = st_ascii_model_read(verbose);

  /* Detach the buffer */

  ASCII_BUFFER = NULL;

  return (model);
}

/****************************************************************************/
/*!
 **   Read the neighborhood definition file
 **
 ** \return  Pointer to the Neigh structure
 **
 ** \param[in]  file_name  Name of the ASCII file
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 ** \remark  For MOVING neighborhood, only isotropic case is considered
 **
 *****************************************************************************/
GEOSLIB_API Neigh *ascii_neigh_read(const char *file_name, int verbose)

{
  Neigh *neigh;
  FILE *file;
  int type, idim, ndim, flag_sector, flag_xvalid, nmini, nmaxi, nsect, nsmax,
      skip;
  int flag_aniso, flag_rotation, lec, jdim;
  double width, dmax;
  VectorDouble radius;
  VectorDouble nbgh_coeffs;
  VectorDouble nbgh_rotmat;

  /* Initializations */

  neigh = (Neigh *) NULL;

  /* Opening the Data file */

  file = st_file_open(file_name, "Neigh", OLD, verbose);
  if (file == (FILE *) NULL) return (neigh);

  /* Create the Model structure */

  if (st_record_read("Space Dimension", "%d", &ndim)) goto label_end;
  if (st_record_read("Neighborhood type", "%d", &type)) goto label_end;

  /* Core allocation */

  radius.resize(ndim);
  nbgh_coeffs.resize(ndim);
  nbgh_rotmat.resize(ndim * ndim);

  switch (type)
  {
    case NEIGH_UNIQUE:
      neigh = neigh_init_unique(ndim);
      break;

    case NEIGH_BENCH:
      if (st_record_read("Flag for Cross-validation", "%d", &flag_xvalid))
        goto label_end;
      if (st_record_read("Bench Width", "%lf", &width)) goto label_end;
      neigh = neigh_init_bench(ndim, flag_xvalid, width);
      break;

    case NEIGH_MOVING:
      flag_aniso = flag_rotation = 0;
      if (st_record_read("Flag for Cross-validaiton", "%d", &flag_xvalid))
        goto label_end;
      if (st_record_read("Neighborhood sector search", "%d", &flag_sector))
        goto label_end;
      if (st_record_read("Neighborhood Width", "%lf", &width)) goto label_end;
      if (st_record_read("Minimum Number of samples", "%d", &nmini))
        goto label_end;
      if (st_record_read("Maximum Number of samples", "%d", &nmaxi))
        goto label_end;
      if (st_record_read("Optimum Number of samples per sector", "%d", &nsect))
        goto label_end;
      if (st_record_read("Maximum Number of samples per sector", "%d", &nsmax))
        goto label_end;
      if (st_record_read("Maximum Isotropic Radius", "%lf", &dmax))
        goto label_end;
      if (st_record_read("Flag for Anisotropy", "%d", &flag_aniso))
        goto label_end;
      if (flag_aniso)
      {
        for (idim = 0; idim < ndim; idim++)
          if (st_record_read("Anisotropy Coefficient", "%lf",
                             &nbgh_coeffs[idim])) goto label_end;

        if (st_record_read("Flag for Anisotropy Rotation", "%d",
                           &flag_rotation)) goto label_end;
        if (flag_rotation)
        {
          for (idim = lec = 0; idim < ndim; idim++)
            for (jdim = 0; jdim < ndim; jdim++, lec++)
              if (st_record_read("Anisotropy Rotation Matrix", "%lf",
                                 &nbgh_rotmat[lec])) goto label_end;
        }
      }
      if (!nbgh_coeffs.empty()) for (idim = 0; idim < ndim; idim++)
        nbgh_coeffs[idim] *= dmax;

      neigh = neigh_init(ndim, NEIGH_MOVING, flag_xvalid, flag_sector,
                         flag_aniso, flag_rotation, 0, nmini, nmaxi, nsect,
                         nsmax, 0, 0., dmax, 0., nbgh_coeffs, nbgh_rotmat,
                         VectorDouble());
      break;

    case NEIGH_IMAGE:
      if (st_record_read("Flag for Cross-Validation", "%d", &flag_xvalid))
        goto label_end;
      if (st_record_read("Skipping factor", "%d", &skip)) goto label_end;
      for (idim = 0; idim < ndim; idim++)
        if (st_record_read("Image Neighborhood Radius", "%lf", &radius[idim]))
          goto label_end;
      neigh = neigh_init_image(ndim, flag_xvalid, skip, radius);
  }

  label_end: if (debug_query("interface") || debug_query("nbgh"))
    neigh_print(neigh);
  st_file_close(file);
  return (neigh);
}

/****************************************************************************/
/*!
 **   Read the Rule definition file
 **
 ** \return  Pointer to the Rule structure
 **
 ** \param[in]  file_name  Name of the ASCII file
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
GEOSLIB_API Rule *ascii_rule_read(const char *file_name, int verbose)
{
  Rule *rule;
  FILE *file;
  int mode_rule, i, inode, nb_node, lec, nfacies, ny1, ny2, ngrf;
  double rho, slope, sh_down, sh_up;
  VectorInt nodes;
  VectorDouble shift(3);

  /* Initializations */

  rule = (Rule *) NULL;

  /* Opening the Data file */

  file = st_file_open(file_name, "Rule", OLD, verbose);
  if (file == (FILE *) NULL) return (rule);

  /* Create the Rule structure */

  if (st_record_read("Rule definition", "%d", &mode_rule)) goto label_end;
  if (st_record_read("Correlation Coefficient of GRFs", "%lf", &rho))
    goto label_end;
  if (st_record_read("Slope for Shadow Rule", "%lf", &slope)) goto label_end;
  if (st_record_read("Lower Threshold for Shadow Rule", "%lf", &sh_down))
    goto label_end;
  if (st_record_read("Upper Threshold for Shadow Rule", "%lf", &sh_up))
    goto label_end;
  if (st_record_read("Shift along first direction", "%lf", &shift[0]))
    goto label_end;
  if (st_record_read("Shift along second direction", "%lf", &shift[1]))
    goto label_end;
  if (st_record_read("Shift along third direction", "%lf", &shift[2]))
    goto label_end;

  /* Read the number of nodes */

  if (st_record_read("Number of Rule Nodes", "%d", &nb_node)) goto label_end;
  nodes.resize(6 * nb_node);

  /* Loop on the nodes for reading: */
  /* - from_type: Type of the parent */
  /* - from_rank: Rank of the parent */
  /* - from_vers: Orientation of the parent */
  /* - node_type: 0 (idle) - 1 (Thresh along Y1) - 2 (Thresh along Y2) */
  /* - node_rank: Rank of the node (starting from 1) */
  /* - facies   : Rank of the facies */
  for (inode = lec = 0; inode < nb_node; inode++)
    for (i = 0; i < 6; i++)
      if (st_record_read("Rule Node Definition", "%d", &nodes[lec++]))
        goto label_end;

  /* Create the structure */

  rule = rule_init(mode_rule, rho, slope, sh_down, sh_up, shift, nodes,
                   &nfacies, &ngrf, &ny1, &ny2);

  label_end: st_file_close(file);
  return (rule);
}

/****************************************************************************/
/*!
 **   Read the Simulation Characteristics
 **
 ** \param[in]  file_name  Name of the ASCII file
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 ** \param[out]  nbsimu    Number of simulations
 ** \param[out]  nbtuba    Number of turning bands
 ** \param[out]  seed      Seed for the random number generator
 **
 *****************************************************************************/
GEOSLIB_API void ascii_simu_read(char *file_name,
                                 int verbose,
                                 int *nbsimu,
                                 int *nbtuba,
                                 int *seed)
{
  FILE *file;

  /* Initializations */

  (*nbsimu) = 0;
  (*nbtuba) = 100;
  (*seed) = 0;

  /* Opening the Simulation Definition file */

  file = st_file_open(file_name, "Simu", OLD, verbose);
  if (file == (FILE *) NULL) return;

  /* Read the parameters */

  if (st_record_read("Number of simulations", "%d", nbsimu)) return;
  if (st_record_read("Number of Turning Bands", "%d", nbtuba)) return;
  if (st_record_read("Random Seed", "%d", seed)) return;

  st_file_close(file);
  return;
}

/****************************************************************************/
/*!
 **   Write a Db
 **
 ** \return  Error returned code
 **
 ** \param[in]  file_name    Name of the ASCII file
 ** \param[in]  db           Pointer to the Db structure to be written
 ** \param[in]  no_grid      1 to forbid the Grid organization in printout
 ** \param[in]  flag_calcul  1 if the data must be printed (for grid only)
 ** \param[in]  verbose      Verbose option if the file cannot be opened
 **
 *****************************************************************************/
GEOSLIB_API int ascii_db_write(const char *file_name,
                               Db *db,
                               int no_grid,
                               int flag_calcul,
                               int verbose)
{
  FILE *file;
  int error, idim, flag_grid;

  /* Initializations */

  error = 1;
  flag_grid = is_grid(db);
  if (no_grid) flag_grid = 0;

  /* Opening the Data file */

  file = st_file_open(file_name, "Db", NEW, verbose);
  if (file == (FILE *) NULL) return (1);

  /* Writing the file organization */

  st_record_write("%d", flag_grid);
  st_record_write("#", "File organization (0:Points; 1:Grid)");

  if (flag_grid)
  {

    /* Writing the header */

    st_record_write("%d", db->getNDim());
    st_record_write("#", "Space Dimension");

    /* Writing the grid characteristics */

    st_record_write("#", "Grid characteristics (NX,X0,DX)");
    for (idim = 0; idim < db->getNDim(); idim++)
    {
      st_record_write("%d",  db->getNX(idim));
      st_record_write("%lf", db->getX0(idim));
      st_record_write("%lf", db->getDX(idim));
      st_record_write("\n");
    }
  }

  /* Writing the tail of the file */

  if (st_variables_write(file, db, flag_grid, flag_calcul)) goto label_end;

  /* Set the error return code */

  st_create_message(file_name, "Db", verbose);
  error = 0;

  label_end: st_file_close(file);
  return (error);
}

/****************************************************************************/
/*!
 **   Write a Model
 **
 ** \param[in]  model      Pointer to the Model structure to be written
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
static void st_ascii_model_write(Model *model, int verbose)
{
  ADriftElem *drift;
  int icova, ibfl, ivar, jvar;

  /* Write the Model structure */
  st_record_write("%d",  model->getDimensionNumber());
  st_record_write("%d",  model->getVariableNumber());
  st_record_write("%lf", model->getField());
  st_record_write("%lf", model->getContext().getBallRadius());
  st_record_write("#", "General parameters");
  st_record_write("%d", model->getCovaNumber());
  st_record_write("#", "Number of basic covariance terms");
  st_record_write("%d", model->getDriftNumber());
  st_record_write("#", "Number of drift terms");

  /* Writing the covariance part */

  for (icova = 0; icova < model->getCovaNumber(); icova++)
  {
    CovAniso* cova = model->getCova(icova);
    st_record_write("%d",  cova->getType());
    st_record_write("%lf", cova->getRange());
    st_record_write("%lf", cova->getParam());
    st_record_write("#", "Covariance characteristics");

    // Writing the Anisotropy information

    st_record_write("%d", cova->getFlagAniso());
    st_record_write("#", "Anisotropy Flag");

    if (! cova->getFlagAniso()) continue;
    for (int idim = 0; idim < model->getDimensionNumber(); idim++)
      st_record_write("%lf", cova->getAnisoCoeffs(idim));
    st_record_write("#", "Anisotropy Coefficients");
    st_record_write("%d", cova->getFlagRotation());
    st_record_write("#", "Anisotropy Rotation Flag");
    if (! cova->getFlagRotation()) continue;
    // Storing the rotation matrix by Column (compatibility)
    for (int idim = 0; idim < model->getDimensionNumber(); idim++)
      for (int jdim = 0; jdim < model->getDimensionNumber(); jdim++)
        st_record_write("%lf", cova->getAnisoRotMat(jdim,idim));
    st_record_write("#", "Anisotropy Rotation Matrix");
  }

  /* Writing the drift part */

  for (ibfl = 0; ibfl < model->getDriftNumber(); ibfl++)
  {
    drift = model->getDrift(ibfl);
    st_record_write("%d", drift->getType());
    st_record_write("#", "Drift characteristics");
  }

  /* Writing the matrix of means (if nbfl <= 0) */

  if (model->getDriftNumber() <= 0)
    for (ivar = 0; ivar < model->getVariableNumber(); ivar++)
  {
    st_record_write("%lf", model->getContext().getMean(ivar));
    st_record_write("#", "Mean of Variables");
  }

  /* Writing the matrices of sills (optional) */

  for (icova = 0; icova < model->getCovaNumber(); icova++)
  {
    for (ivar = 0; ivar < model->getVariableNumber(); ivar++)
      for (jvar = 0; jvar < model->getVariableNumber(); jvar++)
        st_record_write("%lf",
                        model->getSill(icova,ivar,jvar));
    st_record_write("#", "Matrix of sills");
  }

  /* Writing the variance-covariance at the origin (optional) */

  for (ivar = 0; ivar < model->getVariableNumber(); ivar++)
    for (jvar = 0; jvar < model->getVariableNumber(); jvar++)
      st_record_write("%lf",
                      model->getContext().getCovar0(ivar, jvar));
  st_record_write("#", "Var-Covar at origin");

  return;
}

/****************************************************************************/
/*!
 **   Write a Model
 **
 ** \return  Error returned code
 **
 ** \param[in]  file_name  Name of the ASCII file
 ** \param[in]  model      Pointer to the Model structure to be written
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
GEOSLIB_API int ascii_model_write(const char *file_name,
                                  Model *model,
                                  int verbose)
{
  FILE *file;

  /* Opening the Data file */

  file = st_file_open(file_name, "Model", NEW, verbose);
  if (file == (FILE *) NULL) return (1);

  /* Write the Model structure */

  st_ascii_model_write(model, verbose);

  /* Set the error return code */

  st_create_message(file_name, "Model", verbose);

  /* Close the file */

  st_file_close(file);

  return (0);
}

/****************************************************************************/
/*!
 **   Write a Model into a buffer
 **
 ** \return  Error returned code
 **
 ** \param[in]  model      Pointer to the Model structure to be written
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 ** \param[out] buffer       Buffer containing the contents of Model
 ** \param[out] buf_length   Length of the output buffer
 **
 *****************************************************************************/
GEOSLIB_API int ascii_model_write_buffer(Model *model,
                                         int verbose,
                                         char **buffer,
                                         int *buf_length)
{
  /* Initiate the buffer */

  st_buffer_initiate();

  /* Write the Model structure */

  st_ascii_model_write(model, verbose);

  /* Set the error return code */

  st_create_message(NULL, "Model", verbose);

  /* Close the buffer */

  st_buffer_close(buffer, buf_length);

  return (0);
}

/****************************************************************************/
/*!
 **   Write a Neigh
 **
 ** \return  Error returned code
 **
 ** \param[in]  file_name  Name of the ASCII file
 ** \param[in]  neigh      Pointer to the Neigh structure to be written
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
GEOSLIB_API int ascii_neigh_write(const char *file_name,
                                  Neigh *neigh,
                                  int verbose)
{
  FILE *file;
  int error, idim, jdim, ecr;

  /* Opening the Data file */

  file = st_file_open(file_name, "Neigh", NEW, verbose);
  if (file == (FILE *) NULL) return (1);

  /* Create the Model structure */

  st_record_write("%d", neigh->getNDim());
  st_record_write("#", "Space Dimension");
  st_record_write("%d", neigh->getType());
  st_record_write("#", "Neighborhood type");

  switch (neigh->getType())
  {
    case NEIGH_UNIQUE:
      break;

    case NEIGH_BENCH:
      st_record_write("%d", neigh->getFlagXvalid());
      st_record_write("#", "Cross-Validation flag");
      st_record_write("%lf", neigh->getWidth());
      st_record_write("#", "Bench Width");
      break;

    case NEIGH_MOVING:
      st_record_write("%d", neigh->getFlagXvalid());
      st_record_write("#", "Cross-Validation flag");
      st_record_write("%d", neigh->getFlagSector());
      st_record_write("#", "Use angular sectors");
      st_record_write("%lf", neigh->getWidth());
      st_record_write("#", "Bench Width");
      st_record_write("%d", neigh->getNMini());
      st_record_write("%d", neigh->getNMaxi());
      st_record_write("%d", neigh->getNSect());
      st_record_write("%d", neigh->getNSMax());
      st_record_write("#", "Parameters (nmini,nmaxi,nsect,nsmax)");
      st_record_write("%lf", neigh->getRadius());
      st_record_write("#", "Maximum distance radius");
      st_record_write("%d", neigh->getFlagAniso());
      st_record_write("#", "Anisotropy Flag");
      if (!neigh->getFlagAniso()) break;
      for (idim = 0; idim < neigh->getNDim(); idim++)
        st_record_write("%lf", neigh->getAnisoCoeff(idim));
      st_record_write("#", "Anisotropy Coefficients");
      st_record_write("%d", neigh->getFlagRotation());
      st_record_write("#", "Anisotropy Rotation Flag");
      if (!neigh->getFlagRotation()) break;
      for (idim = ecr = 0; idim < neigh->getNDim(); idim++)
        for (jdim = 0; jdim < neigh->getNDim(); jdim++, ecr++)
          st_record_write("%lf", neigh->getAnisoRotMat(ecr));
      st_record_write("#", "Anisotropy Rotation Matrix");
      break;

    case NEIGH_IMAGE:
      st_record_write("%d", neigh->getFlagXvalid());
      st_record_write("#", "Cross-Validation flag");
      st_record_write("%d", neigh->getSkip());
      for (idim = 0; idim < neigh->getNDim(); idim++)
        st_record_write("%lf", neigh->getImageRadius(idim));
      st_record_write("#", "Image neighborhood parameters");
  }

  /* Set the error return code */

  st_create_message(file_name, "Neigh", verbose);
  error = 0;

  st_file_close(file);
  return (error);
}

/****************************************************************************/
/*!
 **   Write a Vario (internal function)
 **
 ** \param[in]  vario        Pointer to the Vario structure to be written
 ** \param[in]  flag_calcul  1 if variogram calculations are stored
 **
 *****************************************************************************/
static void st_ascii_vario_write(Vario *vario, int flag_calcul)
{
  int idim, i, ivar, jvar, idir;
  double value;

  /* Write the Vario structure */

  st_record_write("%d", vario->getDimensionNumber());
  st_record_write("#", "Space Dimension");
  st_record_write("%d", vario->getVariableNumber());
  st_record_write("#", "Number of variables");
  st_record_write("%d", vario->getDirectionNumber());
  st_record_write("#", "Number of directions");
  st_record_write("%lf", vario->getScale());
  st_record_write("#", "Scale");
  st_record_write("%d", flag_calcul);
  st_record_write("#", "Calculation Flag");

  /* Dumping the Variances */

  if (flag_calcul)
  {
    st_record_write("#", "Variance");
    for (ivar = 0; ivar < vario->getVariableNumber(); ivar++)
    {
      for (jvar = 0; jvar < vario->getVariableNumber(); jvar++)
        st_record_write("%lf", vario->getVars(ivar,jvar));
      st_record_write("\n");
    }
  }

  /* Loop on the directions */

  for (idir = 0; idir < vario->getDirectionNumber(); idir++)
  {
    st_record_write("#", "Direction characteristics");
    Dir& dir = vario->getDirs(idir);
    st_record_write("%d", dir.getLagRegular());
    st_record_write("#", "Regular lags");
    st_record_write("%d", dir.getNPas());
    st_record_write("#", "Number of lags");
    st_record_write("%d", dir.getOptionCode());
    st_record_write("%lf", dir.getTolCode());
    st_record_write("#", "Code selection: Option - Tolerance");
    st_record_write("%lf", dir.getDPas());
    st_record_write("#", "Lag value");
    st_record_write("%lf", dir.getTolDist());
    st_record_write("#", "Tolerance on distance");
    st_record_write("%lf", dir.getTolAngle());
    st_record_write("#", "Tolerance on angle");

    for (idim = 0; idim < vario->getDimensionNumber(); idim++)
      st_record_write("%lf", dir.getCodir(idim));
    st_record_write("#", "Direction coefficients");

    for (idim = 0; idim < vario->getDimensionNumber(); idim++)
      st_record_write("%lf", dir.getGrincr(idim));
    st_record_write("#", "Direction increments on grid");

    if (!flag_calcul) continue;
    st_record_write("#", "Variogram results (Weight, Distance, Variogram)");
    for (i = 0; i < dir.getSize(); i++)
    {
      value = FFFF(dir.getSw(i)) ? 0. : dir.getSw(i);
      st_record_write("%lf", value);
      value = FFFF(dir.getHh(i)) ? 0. : dir.getHh(i);
      st_record_write("%lf", value);
      value = FFFF(dir.getGg(i)) ? 0. : dir.getGg(i);
      st_record_write("%lf", value);
      st_record_write("\n");
    }
  }

  return;
}

/****************************************************************************/
/*!
 **   Write a Vario
 **
 ** \return  Error returned code
 **
 ** \param[in]  file_name    Name of the ASCII file
 ** \param[in]  vario        Pointer to the Vario structure to be written
 ** \param[in]  verbose      Verbose option if the file cannot be opened
 ** \param[in]  flag_calcul  1 if variogram calculations are stored
 **
 *****************************************************************************/
GEOSLIB_API int ascii_vario_write(const char *file_name,
                                  Vario *vario,
                                  int verbose,
                                  int flag_calcul)
{
  FILE *file;

  /* Opening the Data file */

  file = st_file_open(file_name, "Vario", NEW, verbose);
  if (file == (FILE *) NULL) return (1);

  /* Write the Vario structure */

  st_ascii_vario_write(vario, flag_calcul);

  /* Set the error return code */

  st_create_message(file_name, "Vario", verbose);

  /* Close the file */

  st_file_close(file);

  return (0);
}

/****************************************************************************/
/*!
 **   Write a Vario into a buffer
 **
 ** \return  Error returned code
 **
 ** \param[in]  vario        Pointer to the Vario structure to be written
 ** \param[in]  verbose      Verbose option if the file cannot be opened
 ** \param[in]  flag_calcul  1 if variogram calculations are stored
 **
 ** \param[out] buffer       Buffer containing the contents of Vario
 ** \param[out] buf_length   Length of the output buffer
 **
 *****************************************************************************/
GEOSLIB_API int ascii_vario_write_buffer(Vario *vario,
                                         int verbose,
                                         int flag_calcul,
                                         char **buffer,
                                         int *buf_length)
{
  /* Initiate the buffer */

  st_buffer_initiate();

  /* Write the Vario structure */

  st_ascii_vario_write(vario, flag_calcul);

  /* Set the error return code */

  st_create_message(NULL, "Vario", verbose);

  /* Close the buffer */

  st_buffer_close(buffer, buf_length);

  return (0);
}

/****************************************************************************/
/*!
 **   Fill the Rule characteristics recursively
 **
 ** \param[in]  file       FILE structure
 ** \param[in]  node       Node structure
 ** \param[in]  from_type  Type of the calling node
 ** \param[in]  from_rank  Rank of the calling node
 ** \param[in]  from_vers  Orientation of the calling node
 ** \param[in]  rank       Current rank
 **
 ** \param[out]  rank      Current rank
 **
 *****************************************************************************/
static void st_rule_define(FILE *file,
                           Node *node,
                           int from_type,
                           int from_rank,
                           int from_vers,
                           int *rank)
{
  int cur_rank;

  /* Calling node */

  st_record_write("%d", from_type);
  st_record_write("%d", from_rank);
  st_record_write("%d", from_vers);

  /* Current node */

  st_record_write("%d", node->getOrient());
  if (IFFFF(node->getFacies()))
  {
    cur_rank = *rank = (*rank) + 1;
    st_record_write("%d", cur_rank);
    st_record_write("%d", 0);
  }
  else
  {
    cur_rank = *rank;
    st_record_write("%d", cur_rank);
    st_record_write("%d", node->getFacies());
  }

  /* Comment */

  st_record_write("#", "Node characteristics");

  if (node->getR1() != (Node *) NULL)
    st_rule_define(file, node->getR1(), node->getOrient(), cur_rank, 1, rank);
  if (node->getR2() != (Node *) NULL)
    st_rule_define(file, node->getR2(), node->getOrient(), cur_rank, 2, rank);

  return;
}

/****************************************************************************/
/*!
 **   Write a Rule
 **
 ** \return  Error returned code
 **
 ** \param[in]  file_name  Name of the ASCII file
 ** \param[in]  rule       Pointer to the Rule structure to be written
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
GEOSLIB_API int ascii_rule_write(const char *file_name, Rule *rule, int verbose)
{
  FILE *file;
  int error, nb_node, nfacies, nmax_tot, ny1_tot, ny2_tot, rank;
  double prop_tot;

  /* Opening the Data file */

  error = 1;
  file = st_file_open(file_name, "Rule", NEW, verbose);
  if (file == (FILE *) NULL) return (1);

  /* Create the Rule structure */

  st_record_write("%d", rule->getModeRule());
  st_record_write("#", "Type of Rule");
  st_record_write("%lf", rule->getRho());
  st_record_write("#", "Correlation coefficient between GRFs");
  st_record_write("%lf", rule->getSlope());
  st_record_write("%lf", rule->getShDown());
  st_record_write("%lf", rule->getShDsup());
  st_record_write("#", "Parameters for Shadow option");
  st_record_write("%lf", rule->getShift(0));
  st_record_write("%lf", rule->getShift(1));
  st_record_write("%lf", rule->getShift(2));
  st_record_write("#", "Parameters for Shift option");

  /* Count the number of nodes */

  rule->statistics(0,&nb_node,&nfacies,&nmax_tot,&ny1_tot,&ny2_tot,&prop_tot);
  st_record_write("%d", nb_node);
  st_record_write("#", "Number of nodes");

  /* Fill the nodes characteristics recursively */

  rank = 0;
  st_rule_define(file, rule->getMainNode(), 0, 0, 0, &rank);

  /* Set the error return code */

  st_create_message(file_name, "Rule", verbose);
  error = 0;

  st_file_close(file);
  return (error);
}

/****************************************************************************/
/*!
 **   Write Polygons
 **
 ** \return  Error returned code
 **
 ** \param[in]  file_name  Name of the ASCII file
 ** \param[in]  polygon    Pointer to the Polygons structure to be written
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
GEOSLIB_API int ascii_polygon_write(const char *file_name,
                                    Polygons *polygon,
                                    int verbose)
{
  FILE *file;
  int i, error, ipol;

  /* Opening the Data file */

  file = st_file_open(file_name, "Polygon", NEW, verbose);
  if (file == (FILE *) NULL) return (1);

  /* Create the Model structure */

  st_record_write("%d", polygon->getPolySetNumber());
  st_record_write("#", "Number of Polygons");

  /* Writing the covariance part */

  for (ipol = 0; ipol < polygon->getPolySetNumber(); ipol++)
  {
    const PolySet& polyset = polygon->getPolySet(ipol);
    st_record_write("%d", polyset.getNVertices());
    st_record_write("#", "Number of Vertices");

    for (i = 0; i < polyset.getNVertices(); i++)
    {
      st_record_write("%lf", polyset.getX(i));
      st_record_write("%lf", polyset.getY(i));
      st_record_write("\n");
    }
  }

  /* Set the error return code */

  st_create_message(file_name, "Polygon", verbose);
  error = 0;

  st_file_close(file);
  return (error);
}

/****************************************************************************/
/*!
 **   Read the Polygons definition file from an ASCII file
 **
 ** \return  Pointer to the Polygons structure
 **
 ** \param[in]  file_name  Name of the ASCII file
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
GEOSLIB_API Polygons *ascii_polygon_read(const char *file_name, int verbose)
{
  Polygons *polygon;
  FILE *file;
  int error, npol, ipol, nvert, i;
  double zmin, zmax;
  VectorDouble x, y;

  /* Initializations */

  error = 1;
  polygon = (Polygons *) NULL;
  zmin = zmax = TEST;

  /* Opening the Data file */

  file = st_file_open(file_name, "Polygon", OLD, verbose);
  if (file == (FILE *) NULL) return (polygon);

  /* Create the Model structure */

  if (st_record_read("Number of Polygons", "%d", &npol)) goto label_end;
  polygon = polygon_create();
  if (polygon == (Polygons *) NULL) goto label_end;

  /* Loop on the PolySets */

  for (ipol = 0; ipol < npol; ipol++)
  {
    if (st_record_read("Number of Vertices", "%d", &nvert)) goto label_end;
    x.resize(nvert);
    y.resize(nvert);

    /* Loop on the Vertices */

    for (i = 0; i < nvert; i++)
    {
      if (st_record_read("X-Coordinate of a Polyset", "%lf", &x[i]))
        goto label_end;
      if (st_record_read("Y-Coordinate of a Polyset", "%lf", &y[i]))
        goto label_end;
    }

    /* Add the polyset */

    polygon = polygon_add(polygon, x, y, zmin, zmax);
    if (polygon == (Polygons *) NULL) goto label_end;
  }

  /* Set the error return code */

  error = 0;

  label_end: if (error) polygon = polygon_free(polygon);
  if (debug_query("interface")) polygon_print(polygon, 1);
  st_file_close(file);
  return (polygon);
}

/****************************************************************************/
/*!
 **   Write an Anam structure
 **
 ** \return  Error returned code
 **
 ** \param[in]  file_name    Name of the ASCII file
 ** \param[in]  anam         Pointer to the Anam structure to be written
 ** \param[in]  verbose      Verbose option if the file cannot be opened
 ** \param[in]  flag_calcul  1 if anamorphosis calculations are stored
 **
 *****************************************************************************/
GEOSLIB_API int ascii_anam_write(const char *file_name,
                                 const Anam* anam,
                                 int verbose,
                                 int flag_calcul)
{
  FILE *file;
  int error;

  /* Opening the Data file */

  file = st_file_open(file_name, "Anam", NEW, verbose);
  if (file == (FILE *) NULL) return (1);

  /* Write the Anam structure */

  st_record_write("%d", anam->getType());
  st_record_write("#", "Type of Anamorphosis");

  /* Hermitian case */

  if (anam->getType() == ANAM_HERMITIAN)
  {
    const AnamHermite* anam_hermite = dynamic_cast<const AnamHermite*>(anam);
    st_record_write("%d", anam_hermite->getNbPoly());
    st_record_write("#", "Number of Hermite Polynomials");
    st_record_write("%lf", anam_hermite->getRCoef());
    st_record_write("#", "Change of support coefficient");
    st_record_write("%lf", anam_hermite->getAzmin());
    st_record_write("%lf", anam_hermite->getAzmax());
    st_record_write("#", "Absolute Values for Z");
    st_record_write("%lf", anam_hermite->getAymin());
    st_record_write("%lf", anam_hermite->getAymax());
    st_record_write("#", "Absolute Values for Y");
    st_record_write("%lf", anam_hermite->getPzmin());
    st_record_write("%lf", anam_hermite->getPzmax());
    st_record_write("#", "Practical Values for Z");
    st_record_write("%lf", anam_hermite->getPymin());
    st_record_write("%lf", anam_hermite->getPymax());
    st_record_write("#", "Practical Values for Y");
    st_record_write("%d", flag_calcul);
    st_record_write("#", "Storage of Calculation results");

    if (flag_calcul)
    {
      st_record_write("%lf", anam_hermite->getVariance());
      st_record_write("#", "Calculated variance");
      st_table_write("Hermite Polynomial", anam_hermite->getNbPoly(),
                     anam_hermite->getPsiHn().data());
    }
  }
  else if (anam->getType() == ANAM_EMPIRICAL)
  {
    const AnamEmpirical* anam_empirical = dynamic_cast<const AnamEmpirical*>(anam);
    st_record_write("%d", anam_empirical->getNDisc());
    st_record_write("#", "Number of Discretization lags");
    st_record_write("%lf", anam_empirical->getAzmin());
    st_record_write("%lf", anam_empirical->getAzmax());
    st_record_write("#", "Absolute Values for Z");
    st_record_write("%lf", anam_empirical->getAymin());
    st_record_write("%lf", anam_empirical->getAymax());
    st_record_write("#", "Absolute Values for Y");
    st_record_write("%lf", anam_empirical->getPzmin());
    st_record_write("%lf", anam_empirical->getPzmax());
    st_record_write("#", "Practical Values for Z");
    st_record_write("%lf", anam_empirical->getPymin());
    st_record_write("%lf", anam_empirical->getPymax());
    st_record_write("#", "Practical Values for Y");
    st_record_write("%ld", anam_empirical->getSigma2e());
    st_record_write("#", "additional variance");
    st_record_write("%d", flag_calcul);
    st_record_write("#", "Storage of Calculation results");

    if (flag_calcul)
      st_table_write("Coefficients", 2 * anam_empirical->getNDisc(),
                     anam_empirical->getTDisc().data());
  }
  else if (anam->getType() == ANAM_DISCRETE_DD)
  {
    const AnamDiscreteDD* anam_discrete_DD = dynamic_cast<const AnamDiscreteDD*>(anam);
    st_record_write("%d", anam_discrete_DD->getNCut());
    st_record_write("#", "Number of cutoffs");
    st_record_write("%d", anam_discrete_DD->getNClass());
    st_record_write("#", "Number of classes");
    st_record_write("%d", anam_discrete_DD->getNElem());
    st_record_write("#", "Number of elements");
    st_record_write("%lf", anam_discrete_DD->getSCoef());
    st_record_write("#", "Change of support coefficient");
    st_record_write("%lf", anam_discrete_DD->getMu());
    st_record_write("#", "Additional Mu coefficient");
    st_table_write("Cutoff value", anam_discrete_DD->getNCut(),
                   anam_discrete_DD->getZCut().data());
    st_table_write("PCA Z2Y", anam_discrete_DD->getNCut() * anam_discrete_DD->getNCut(),
                   anam_discrete_DD->getPcaZ2F().data());
    st_table_write("PCA Y2Z", anam_discrete_DD->getNCut() * anam_discrete_DD->getNCut(),
                   anam_discrete_DD->getPcaF2Z().data());
    st_record_write("%d", flag_calcul);
    st_record_write("#", "Storage of Calculation results");
    if (flag_calcul)
      st_table_write("DD Stats",
                     anam_discrete_DD->getNClass() * anam_discrete_DD->getNElem(),
                     anam_discrete_DD->getStats().getValues().data());
  }
  else if (anam->getType() == ANAM_DISCRETE_IR)
  {
    const AnamDiscreteIR* anam_discrete_IR = dynamic_cast<const AnamDiscreteIR*>(anam);
    st_record_write("%d", anam_discrete_IR->getNCut());
    st_record_write("#", "Number of cutoffs");
    st_record_write("%d", anam_discrete_IR->getNClass());
    st_record_write("#", "Number of classes");
    st_record_write("%d", anam_discrete_IR->getNElem());
    st_record_write("#", "Number of elements");
    st_record_write("%lf", anam_discrete_IR->getRCoef());
    st_record_write("#", "Change of support coefficient");
    st_table_write("Cutoff value", anam_discrete_IR->getNCut(),
                   anam_discrete_IR->getZCut().data());
    st_record_write("%d", flag_calcul);
    st_record_write("#", "Storage of Calculation results");
    if (flag_calcul)
      st_table_write("IR Stats",
                     anam_discrete_IR->getNClass() * anam_discrete_IR->getNElem(),
                     anam_discrete_IR->getStats().getValues().data());
  }

  /* Set the error return code */

  st_create_message(file_name, "Anam", verbose);
  error = 0;

  st_file_close(file);
  return (error);
}

/****************************************************************************/
/*!
 **   Read the Anam definition file
 **
 ** \return  Pointer to the Anam structure
 **
 ** \param[in]  file_name  Name of the ASCII file
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
GEOSLIB_API Anam *ascii_anam_read(const char *file_name, int verbose)
{
  FILE *file;
  double azmin, azmax, aymin, aymax, pzmin, pzmax, pymin, pymax;
  double variance, sigma2e, r, s, mu;
  int nbpoly, flag_calcul, type, ndisc, nCut, nClass, nElem;
  VectorDouble hermite, tdisc, zCut, pcaf2z, pcaz2f, stats;

  /* Initializations */

  variance = r = s = mu = TEST;

  /* Opening the Data file */

  file = st_file_open(file_name, "Anam", OLD, verbose);
  if (file == (FILE *) NULL) return (nullptr);

  /* Read the general anamorphosis parameters */

  if (st_record_read("Anamorphosis Type", "%d", &type)) goto label_end;

  /* Read the remaining parameters */

  if (type == ANAM_HERMITIAN)
  {
    AnamHermite* anam_hermite = new AnamHermite;

    if (st_record_read("Number of Hermite Polynomials", "%d", &nbpoly))
      goto label_end;
    if (st_record_read("Change of Support Coefficient", "%lf", &r))
      goto label_end;
    if (st_record_read("Minimum absolute Z-value", "%lf", &azmin))
      goto label_end;
    if (st_record_read("Maximum absolute Z-value", "%lf", &azmax))
      goto label_end;
    if (st_record_read("Minimum absolute Y-value", "%lf", &aymin))
      goto label_end;
    if (st_record_read("Maximum absolute Y-value", "%lf", &aymax))
      goto label_end;
    if (st_record_read("Minimum Experimental Z-value", "%lf", &pzmin))
      goto label_end;
    if (st_record_read("Maximum Experimental Z-value", "%lf", &pzmax))
      goto label_end;
    if (st_record_read("Minimum Experimental Y-value", "%lf", &pymin))
      goto label_end;
    if (st_record_read("Maximum Experimental Y-value", "%lf", &pymax))
      goto label_end;
    if (st_record_read("Calculation Flag", "%d", &flag_calcul)) goto label_end;

    if (flag_calcul)
    {
      hermite.resize(nbpoly);
      if (st_record_read("Experimental Variance", "%lf", &variance))
        goto label_end;
      if (st_table_read(nbpoly, hermite.data())) goto label_end;
    }

    anam_update_hermitian(anam_hermite, nbpoly, pymin, pzmin, pymax, pzmax, aymin,
                          azmin, aymax, azmax, r, variance, hermite);
    if (debug_query("interface")) anam_hermite->display();
    st_file_close(file);
    return (anam_hermite);
  }
  else if (type == ANAM_EMPIRICAL)
  {
    AnamEmpirical* anam_empirical = new(AnamEmpirical);

    if (st_record_read("Number of Discretization classes", "%d", &ndisc))
      goto label_end;
    if (st_record_read("Minimum absolute Z-value", "%lf", &azmin))
      goto label_end;
    if (st_record_read("Maximum absolute Z-value", "%lf", &azmax))
      goto label_end;
    if (st_record_read("Minimum absolute Y-value", "%lf", &aymin))
      goto label_end;
    if (st_record_read("Maximum absolute Y-value", "%lf", &aymax))
      goto label_end;
    if (st_record_read("Minimum Experimental Z-value", "%lf", &pzmin))
      goto label_end;
    if (st_record_read("Maximum Experimental Z-value", "%lf", &pzmax))
      goto label_end;
    if (st_record_read("Minimum Experimental Y-value", "%lf", &pymin))
      goto label_end;
    if (st_record_read("Maximum Experimental Y-value", "%lf", &pymax))
      goto label_end;
    if (st_record_read("Experimental Error Variance", "%lf", &sigma2e))
      goto label_end;
    if (st_record_read("Calculation Flag", "%d", &flag_calcul)) goto label_end;

    if (flag_calcul)
    {
      tdisc.resize(2 * ndisc);
      if (st_table_read(2 * ndisc, tdisc.data())) goto label_end;
    }
    anam_update_empirical(anam_empirical, ndisc, pymin, pzmin, pymax, pzmax, aymin,
                              azmin, aymax, azmax, sigma2e, tdisc);
    if (debug_query("interface")) anam_empirical->display();
    st_file_close(file);
    return (anam_empirical);
  }
  else if (type == ANAM_DISCRETE_DD)
  {
    AnamDiscreteDD* anam_discrete_DD = new(AnamDiscreteDD);

    if (st_record_read("Number of Cutoffs", "%d", &nCut)) goto label_end;
    if (st_record_read("Number of Classes", "%d", &nClass)) goto label_end;
    if (st_record_read("Number of Statistic Columns", "%d", &nElem))
      goto label_end;
    if (st_record_read("Anamorphosis 's' coefficient", "%lf", &s))
      goto label_end;
    if (st_record_read("Anamorphosis 'mu' coefficient", "%lf", &mu))
      goto label_end;
    zCut.resize(nCut);
    if (st_table_read(nCut, zCut.data())) goto label_end;
    pcaz2f.resize(nCut * nCut);
    pcaf2z.resize(nCut * nCut);
    if (st_table_read(nCut * nCut, pcaz2f.data())) goto label_end;
    if (st_table_read(nCut * nCut, pcaf2z.data())) goto label_end;
    if (st_record_read("Calculation Flag", "%d", &flag_calcul)) goto label_end;

    if (flag_calcul)
    {
      stats.resize(nClass * nElem);
      if (st_table_read(nClass * nElem, stats.data())) goto label_end;
    }
    anam_update_discrete_DD(anam_discrete_DD, nCut, s, mu, zCut, pcaz2f, pcaf2z, stats);
    if (debug_query("interface")) anam_discrete_DD->display();
    st_file_close(file);
    return (anam_discrete_DD);
  }
  else if (type == ANAM_DISCRETE_IR)
  {
    AnamDiscreteIR* anam_discrete_IR = new(AnamDiscreteIR);

    if (st_record_read("Number of Cutoffs", "%d", &nCut)) goto label_end;
    if (st_record_read("Number of Classes", "%d", &nClass)) goto label_end;
    if (st_record_read("Number of Statistic Columns", "%d", &nElem))
      goto label_end;
    if (st_record_read("Anamorphosis 'r' coefficient", "%lf", &r))
      goto label_end;
    zCut.resize(nCut);
    if (st_table_read(nCut, zCut.data())) goto label_end;
    if (st_record_read("Calculation Flag", "%d", &flag_calcul)) goto label_end;

    if (flag_calcul)
    {
      stats.resize(nClass * nElem);
      if (st_table_read(nClass * nElem, stats.data())) goto label_end;
    }
    anam_update_discrete_IR(anam_discrete_IR, nCut, r, zCut, stats);
    if (debug_query("interface")) anam_discrete_IR->display();
    st_file_close(file);
    return (anam_discrete_IR);
  }

  label_end:
  return (nullptr);
}

/****************************************************************************/
/*!
 **   Check if an option is defined in the Options ASCII file
 **
 ** \return  1 if the option is defined and 0 otherwise
 **
 ** \param[in]  file_name    Name of the ASCII file
 ** \param[in]  verbose      Verbose option if the file cannot be opened
 ** \param[in]  option_name  Keyword for the requested option
 ** \param[in]  type         Answer type
 ** \li                      0 : Logical (returned as 0 or 1)
 ** \li                      1 : integer
 ** \li                      2 : real (returned as a double)
 **
 ** \param[out]  answer      Answer
 **
 *****************************************************************************/
GEOSLIB_API int ascii_option_defined(const char *file_name,
                                     int verbose,
                                     const char *option_name,
                                     int type,
                                     void *answer)
{
  FILE *file;
  char keyword[100], keyval[100];
  double rval;
  int lrep, ival;

  /* Initializations */

  lrep = 0;

  /* Opening the Data file */

  file = st_file_open(file_name, "Option", OLD, verbose);
  if (file == (FILE *) NULL) return (lrep);

  /* Implicit loop on the lines of the file */

  while (1)
  {
    if (st_record_read("Option Keyword", "%s", keyword)) goto label_end;
    if (st_record_read("Option Key-value", "%s", keyval)) goto label_end;
    if (strcmp(keyword, option_name)) continue;

    /* The keyword matches the option name */
    switch (type)
    {
      case 0:
        ival = 0;
        if (!strcmp(keyval, "Y") || !strcmp(keyval, "YES")
            || !strcmp(keyval, "y") || !strcmp(keyval, "yes")
            || atoi(keyval) == 1) ival = 1;
        *((int *) answer) = ival;
        break;

      case 1:
        ival = atoi(keyval);
        *((int *) answer) = ival;
        break;

      case 2:
        rval = atof(keyval);
        *((double *) answer) = rval;
        break;
    }
    lrep = 1;
    goto label_end;
  }

  label_end: st_file_close(file);
  return (lrep);
}

/****************************************************************************/
/*!
 **   Write a set of Fracture
 **
 ** \return  Error returned code
 **
 ** \param[in]  file_name  Name of the ASCII file
 ** \param[in]  frac       Pointer to the Frac_Environ structure to be written
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
GEOSLIB_API int ascii_frac_write(const char *file_name,
                                 Frac_Environ *frac,
                                 int verbose)
{
  FILE *file;
  int family, ifault, error;

  /* Opening the Data file */

  file = st_file_open(file_name, "Frac", NEW, verbose);
  if (file == (FILE *) NULL) return (1);

  /* Create the Frac_Environ structure */

  st_record_write("%d", frac->nfamilies);
  st_record_write("#", "Number of families");
  st_record_write("%d", frac->nfaults);
  st_record_write("#", "Number of main faults");
  st_record_write("%lf", frac->xmax);
  st_record_write("#", "Maximum horizontal distance");
  st_record_write("%lf", frac->ymax);
  st_record_write("#", "Maximum vertical distance");
  st_record_write("%lf", frac->deltax);
  st_record_write("#", "Dilation along the horizontal axis");
  st_record_write("%lf", frac->deltay);
  st_record_write("#", "Dilation along the vertical axis");
  st_record_write("%lf", frac->mean);
  st_record_write("#", "Mean of thickness distribution");
  st_record_write("%lf", frac->stdev);
  st_record_write("#", "Stdev of thickness distribution");

  /* Loop on the families */

  for (family = 0; family < frac->nfamilies; family++)
  {
    st_record_write("#", "Characteristics of family");
    const Frac_Fam& frac_fam = frac->frac_fams[family];
    st_record_write("%lf", frac_fam.orient);
    st_record_write("#", "Mean orientation");
    st_record_write("%lf", frac_fam.dorient);
    st_record_write("#", "Tolerance for orientation");
    st_record_write("%lf", frac_fam.theta0);
    st_record_write("#", "Reference Poisson intensity");
    st_record_write("%lf", frac_fam.alpha);
    st_record_write("#", "Power dependency between layer and intensity");
    st_record_write("%lf", frac_fam.ratcst);
    st_record_write("#", "Ratio of cste vs. shaped intensity");
    st_record_write("%lf", frac_fam.prop1);
    st_record_write("#", "Survival probability (constant term)");
    st_record_write("%lf", frac_fam.prop2);
    st_record_write("#", "Survival probability (length dependent term)");
    st_record_write("%lf", frac_fam.aterm);
    st_record_write("#", "Survival probability (cumulative length exponent)");
    st_record_write("%lf", frac_fam.bterm);
    st_record_write("#", "Survival probability (layer thickness exponent)");
    st_record_write("%lf", frac_fam.range);
    st_record_write("#", "Fracture repulsion area Range");
  }

  /* Loop on the main faults */

  for (ifault = 0; ifault < frac->nfaults; ifault++)
  {
    st_record_write("#", "Characteristics of main fault");
    const Frac_Fault& frac_fault = frac->frac_faults[ifault];
    st_record_write("%lf", frac_fault.coord);
    st_record_write("#", "Abscissae of the first Fault point");
    st_record_write("%lf", frac_fault.orient);
    st_record_write("#", "Fault orientation");

    /* Loop on the families */

    for (family = 0; family < frac->nfamilies; family++)
    {
      st_record_write("%lf", frac_fault.thetal[family]);
      st_record_write("%lf", frac_fault.rangel[family]);
      st_record_write("%lf", frac_fault.thetar[family]);
      st_record_write("%lf", frac_fault.ranger[family]);
      st_record_write("#", "Max-left, Range-Left, Max-Right, Range-Right");
    }
  }

  /* Set the error return code */

  st_create_message(file_name, "Frac", verbose);
  error = 0;

  st_file_close(file);
  return (error);
}

/****************************************************************************/
/*!
 **   Read the fracture definition file
 **
 ** \return  Pointer to the Frac_Environ structure
 **
 ** \param[in]  file_name  Name of the ASCII file
 ** \param[in]  verbose    Verbose option if the file cannot be opened
 **
 *****************************************************************************/
GEOSLIB_API Frac_Environ *ascii_frac_read(const char *file_name, int verbose)
{
  Frac_Environ *frac;
  FILE *file;
  int family, ifault, nfamilies, nfaults;
  double thetal, rangel, thetar, ranger, xmax, ymax, deltax, deltay, mean,
      stdev;
  double orient, dorient, coord, theta0, alpha, prop1, prop2, aterm, bterm,
      ratcst, range;

  /* Initializations */

  frac = (Frac_Environ *) NULL;

  /* Opening the Data file */

  file = st_file_open(file_name, "Frac", OLD, verbose);
  if (file == (FILE *) NULL) return (frac);

  /* Create the structure */

  if (st_record_read("Number of Families", "%d", &nfamilies)) goto label_end;
  if (st_record_read("Number of Faults", "%d", &nfaults)) goto label_end;
  if (st_record_read("Maximum X-value", "%lf", &xmax)) goto label_end;
  if (st_record_read("Maximum Y-value", "%lf", &ymax)) goto label_end;
  if (st_record_read("Increment along X", "%lf", &deltax)) goto label_end;
  if (st_record_read("Increment along Y", "%lf", &deltay)) goto label_end;
  if (st_record_read("Mean value", "%lf", &mean)) goto label_end;
  if (st_record_read("Standard Deviation value", "%lf", &stdev)) goto label_end;

  frac = fracture_alloc_environ(nfamilies, xmax, ymax, deltax, deltay, mean,
                                stdev);

  /* Loop on the families */

  for (family = 0; family < nfamilies; family++)
  {
    if (st_record_read("Fracture Orientation", "%lf", &orient)) goto label_end;
    if (st_record_read("Tolerance on the Orientation", "%lf", &dorient))
      goto label_end;
    if (st_record_read("Reference Angle", "%lf", &theta0)) goto label_end;
    if (st_record_read("Fracture 'alpha' coefficient", "%lf", &alpha))
      goto label_end;
    if (st_record_read("Constant Ratio", "%lf", &ratcst)) goto label_end;
    if (st_record_read("First Proportion", "%lf", &prop1)) goto label_end;
    if (st_record_read("Second Proportion", "%lf", &prop2)) goto label_end;
    if (st_record_read("Coefficient 'aterm'", "%lf", &aterm)) goto label_end;
    if (st_record_read("Coefficient 'bterm'", "%lf", &bterm)) goto label_end;
    if (st_record_read("Range value", "%lf", &range)) goto label_end;

    fracture_update_family(frac, family, orient, dorient, theta0, alpha, ratcst,
                           prop1, prop2, aterm, bterm, range);
  }

  /* Loop on the faults */

  for (ifault = 0; ifault < nfaults; ifault++)
  {
    if (st_record_read("Fault coordinate", "%lf", &coord)) goto label_end;
    if (st_record_read("Fault orientation", "%lf", &orient)) goto label_end;
    (void) fracture_add_fault(frac, coord, orient);

    for (family = 0; family < nfamilies; family++)
    {
      if (st_record_read("Coefficient 'theta' on the left", "%lf", &thetal))
        goto label_end;
      if (st_record_read("Coefficient 'range' on the left", "%lf", &rangel))
        goto label_end;
      if (st_record_read("Coefficient 'theta' on the right", "%lf", &thetar))
        goto label_end;
      if (st_record_read("Coefficient 'range' on the right", "%lf", &ranger))
        goto label_end;

      fracture_update_fault(frac, ifault, family, thetal, thetar, rangel,
                            ranger);
    }
  }

  label_end: if (debug_query("interface")) fracture_print(frac);
  st_file_close(file);
  return (frac);
}

/****************************************************************************/
/*!
 **   Write a Table (of real values)
 **
 ** \return  Error returned code
 **
 ** \param[in]  file_name    Name of the ASCII file
 ** \param[in]  verbose      Verbose option if the file cannot be opened
 ** \param[in]  ntab         Number of samples
 ** \param[in]  tab          Array of real values to be written
 **
 *****************************************************************************/
GEOSLIB_API int ascii_table_write(const char *file_name,
                                  int verbose,
                                  int ntab,
                                  double *tab)
{
  FILE *file;

  /* Opening the Data file */

  file = st_file_open(file_name, NULL, NEW, verbose);
  if (file == (FILE *) NULL) return (1);
  st_record_write("%d", ntab);
  st_record_write("\n");
  st_table_write(NULL, ntab, tab);

  /* Set the error return code */

  st_create_message(file_name, "Table", verbose);

  st_file_close(file);
  return (0);
}

/****************************************************************************/
/*!
 **   Write a Table (of integer values)
 **
 ** \return  Error returned code
 **
 ** \param[in]  file_name    Name of the ASCII file
 ** \param[in]  verbose      Verbose option if the file cannot be opened
 ** \param[in]  ntab         Number of samples
 ** \param[in]  itab         Array of integer values to be written
 **
 *****************************************************************************/
GEOSLIB_API int ascii_tablei_write(const char *file_name,
                                   int verbose,
                                   int ntab,
                                   int *itab)
{
  FILE *file;

  /* Opening the Data file */

  file = st_file_open(file_name, NULL, NEW, verbose);
  if (file == (FILE *) NULL) return (1);

  st_record_write("%d", ntab);
  st_record_write("\n");
  st_tablei_write(NULL, ntab, itab);

  /* Set the error return code */

  st_create_message(file_name, "Table", verbose);

  st_file_close(file);
  return (0);
}

/****************************************************************************/
/*!
 **   Read a CSV file and load the results into a Db
 **
 ** \return  Pointer to the Db descriptor
 **
 ** \param[in]  file_name     Name of the ASCII file
 ** \param[in]  verbose       Verbose option if the file cannot be opened
 ** \param[in]  flag_header   1 if the first line of the file contains the
 **                           variable names
 ** \param[in]  nskip         Number of lines to skip ** \param[in]  order  manner in which values in tab are ordered
 ** \param[in]  order         manner in which values in tab are ordered
 ** \param[in]  char_sep      Character used as a column separator
 ** \param[in]  char_dec      Character used as a decimal
 ** \param[in]  na_string     String used for absent information
 ** \param[in]  ncol_max      Maximum number of columns (or -1)
 ** \param[in]  nrow_max      Maximum number of rows (or -1)
 ** \param[in]  flag_add_rank 1 To add the rank number
 **
 *****************************************************************************/
GEOSLIB_API Db *db_read_csv(const char *file_name,
                            int verbose,
                            int flag_header,
                            int nskip,
                            ENUM_LOAD_DATA order,
                            const char *char_sep,
                            const char *char_dec,
                            const char *na_string,
                            int ncol_max,
                            int nrow_max,
                            int flag_add_rank)
{
  Db *db;
  int ncol, nrow;
  VectorString names;
  VectorDouble tab;

  /* Initializations */

  db = (Db *) NULL;

  /* Reading the CSV file */

  if (csv_table_read(file_name, verbose, flag_header, nskip, char_sep, char_dec,
                     na_string, ncol_max, nrow_max, &ncol, &nrow, names, tab))
    goto label_end;

  /* Creating the Db */

  db = db_create_point(nrow, ncol, order, flag_add_rank, tab);
  if (db == (Db *) NULL) goto label_end;

  /* Loading the names */

  for (int i = 0; i < ncol; i++)
  {
    int j = (flag_add_rank) ? i+1 : i;
    if (db_name_set(db, j, names[i])) messerr("Error in db_name_set");
  }

  /* Core deallocation */

  label_end: return (db);
}
