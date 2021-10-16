/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                     */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "Variogram/Vario.hpp"
#include "Anamorphosis/Anam.hpp"
#include "Anamorphosis/AnamDiscreteDD.hpp"
#include "Anamorphosis/AnamDiscreteIR.hpp"
#include "Anamorphosis/AnamEmpirical.hpp"
#include "Anamorphosis/AnamHermite.hpp"
#include "Basic/Utilities.hpp"
#include "Basic/File.hpp"
#include "Covariances/CovAniso.hpp"
#include "geoslib_e.h"
#include "geoslib_enum.h"
#include "geoslib_old_f.h"

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

static bool st_file_exists(const char* file_name)
{
  return gslFileExist(file_name, "r");
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
    long1 = static_cast<int> (strlen(buf));
    long2 = (ASCII_BUFFER != NULL) ? static_cast<int> (strlen(ASCII_BUFFER)) : 0;
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
  (void) gslStrcpy(STUDY, gslArraySize(STUDY), study);
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
  (void) gslStrcpy(FILE_NAME_MEM, gslArraySize(FILE_NAME_MEM), filename);

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
      (void) gslSPrintf(local, LONG_SIZE, "%s (%d)", string, i + 1);
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
      (void) gslSPrintf(local, LONG_SIZE, "%s (%d)", string, i + 1);
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
  if (! st_file_exists(file_name)) return nullptr;
  Db* db = new Db(String(file_name), verbose);
  return (db);
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
GEOSLIB_API Vario *ascii_vario_read(const char *file_name, bool verbose)
{
  if (! st_file_exists(file_name)) return nullptr;
  Vario* vario = new Vario(file_name,verbose);
  return (vario);
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
  if (! st_file_exists(file_name)) return nullptr;
  Model *model = new Model(file_name,verbose);
  model_setup(model);
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
  if (! st_file_exists(file_name)) return nullptr;
  Neigh* neigh = new Neigh(file_name,verbose);
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
  if (! st_file_exists(file_name)) return nullptr;
  Rule* rule = new Rule(file_name,verbose);
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

  st_record_write("%d", anam->getType().getValue());
  st_record_write("#", "Type of Anamorphosis");

  /* Hermitian case */

  if (anam->getType() == EAnam::HERMITIAN)
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
  else if (anam->getType() == EAnam::EMPIRICAL)
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
  else if (anam->getType() == EAnam::DISCRETE_DD)
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
  else if (anam->getType() == EAnam::DISCRETE_IR)
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
  int nbpoly, flag_calcul, atype, ndisc, nCut, nClass, nElem;
  VectorDouble hermite, tdisc, zCut, pcaf2z, pcaz2f, stats;
  EAnam type;

  /* Initializations */

  variance = r = s = mu = TEST;

  /* Opening the Data file */

  file = st_file_open(file_name, "Anam", OLD, verbose);
  if (file == (FILE *) NULL) return (nullptr);

  /* Read the general anamorphosis parameters */

  if (st_record_read("Anamorphosis Type", "%d", &atype)) goto label_end;
  type = EAnam::fromValue(atype);
  if (type == EAnam::UNDEFINED) goto label_end;

  /* Read the remaining parameters */

  if (type == EAnam::HERMITIAN)
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
  else if (type == EAnam::EMPIRICAL)
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
  else if (type == EAnam::DISCRETE_DD)
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
  else if (type == EAnam::DISCRETE_IR)
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
 ** \param[in]  nskip         Number of lines to skip
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

  db = db_create_point(nrow, ncol, ELoadBy::SAMPLE, flag_add_rank, tab);
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
