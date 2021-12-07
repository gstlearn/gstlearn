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
#include "geoslib_old_f.h"
#include "geoslib_f_private.h"
#include "Basic/Utilities.hpp"
#include "Basic/Law.hpp"
#include "Db/Db.hpp"

#include <math.h>


/*! \cond */
// get_keypone default values
#define NDISC_DEF        1000
#define LOW0_DEF        1.e-8
#define LOW1_DEF        1.e-6
#define FLAG_CHECK_DEF      1

#define NPART 5
#define XXF(i)              (desc.xy[2 * (i)])
#define YYF(i)              (desc.xy[2 * (i) + 1])
#define INSIDE(x,y)         ((y) > -EPS && (x) > -EPS && (x) < (xmax + EPS))
#define MEM_LAYER(i)        (locinfo[(i) * ninfos])
#define MEM_THETA1(ifam,i)  (locinfo[(i) * ninfos + 1 + (ifam) * NPART + 0])
#define MEM_THETA2(ifam,i)  (locinfo[(i) * ninfos + 1 + (ifam) * NPART + 1])
#define MEM_PROPSUR(ifam,i) (locinfo[(i) * ninfos + 1 + (ifam) * NPART + 2])
#define MEM_FRAC(ifam,i)    (locinfo[(i) * ninfos + 1 + (ifam) * NPART + 3])
#define MEM_TOTAL(ifam,i)   (locinfo[(i) * ninfos + 1 + (ifam) * NPART + 4])
#define VALID_DISC(idisc)   ((idisc) >= 0 && (idisc) < NDISC)
#define FRACSEG(ifrac,i)    (frac_segs[nbyfrac * (ifrac) + (i)])
#define WELL(iw,i)          (well[2*(iw) + (i)])
#define WELLOUT(is,i)       (wellout[nbyout * (is) + (i)])
#define TRAJ(ip,i)          (traj[2*(ip)+(i)])

static double XORIGIN, STEP, LOW0, LOW1;

static double EPS = 1.e-3;
static int VERBOSE = 0;
static int NDISC = 0;
static int FLAG_CHECK = 0;

/*! \endcond */

/****************************************************************************/
/*!
 **  Count the number of end points
 **
 ** \return  Number of end points
 **
 ** \param[in]  frac_list    Pointer to the Frac_List structure
 **
 *****************************************************************************/
static int st_end_point_count(Frac_List *frac_list)

{
  int number;

  number = 0;
  for (int i = 0; i < frac_list->nfracs; i++)
    number += frac_list->frac_descs[i].npoint;
  return (number);
}

/****************************************************************************/
/*!
 **  Manage the Frac_Desc structure
 **
 ** \param[in]  frac_desc    Pointer to the Frac_Desc structure
 **
 *****************************************************************************/
static void st_frac_desc_init(Frac_Desc &frac_desc)
{
  frac_desc.family = 0;
  frac_desc.npoint = 0;
  frac_desc.orient = 0.;
}

/****************************************************************************/
/*!
 **  Manage the Frac_List structure
 **
 ** \return  Pointer to the Frac_List structure
 **
 ** \param[in]  mode         0 initialization; 1 allocation; -1 deallocation
 ** \param[in]  frac_list    Pointer to the Frac_List structure
 **
 *****************************************************************************/
Frac_List* fracture_manage_list(int mode, Frac_List *frac_list)
{
  int nfracs;

  /* Dispatch */

  switch (mode)
  {
    case 0:
      frac_list = new (Frac_List);
      frac_list->nfracs = 0;
      break;

    case 1:
      if (frac_list == nullptr) frac_list = fracture_manage_list(0, NULL);
      nfracs = frac_list->nfracs;
      frac_list->frac_descs.resize(nfracs + 1);
      st_frac_desc_init(frac_list->frac_descs[nfracs]);
      frac_list->nfracs++;
      break;

    case -1:
      if (frac_list == nullptr) return (frac_list);
      delete frac_list;
      frac_list = (Frac_List*) nullptr;
      break;
  }
  return (frac_list);
}

/****************************************************************************/
/*!
 **  Print a list of simulated fractures
 **
 ** \param[in]  title     Title to be printed
 ** \param[in]  frac_list Frac_List structure
 ** \param[in]  level     0 summary; 1 fracture; 2 end-points
 **
 *****************************************************************************/
void fracture_list_print(const char *title,
                                         Frac_List *frac_list,
                                         int level)
{
  message("%s\n", title);
  message("Current number of simulated fractures = %d\n", frac_list->nfracs);
  if (level <= 0) return;

  /* Loop on the fractures */

  for (int ifrac = 0; ifrac < frac_list->nfracs; ifrac++)
  {
    const Frac_Desc &desc = frac_list->frac_descs[ifrac];
    message("Fracture #%3d (family=%d): %d segment(s) - start at level #%lf\n",
            ifrac + 1, desc.family, desc.npoint - 1, YYF(0));

    if (level <= 1) continue;

    message("     X           Y \n");
    for (int j = 0; j < desc.npoint; j++)
      message(" %10.3lf %10.3lf\n", XXF(j), YYF(j));
  }
  return;
}

/****************************************************************************/
/*!
 **  Allocate the Frac_Environ structure
 **
 ** \return  Pointer to the Frac_Environ structure
 **
 ** \param[in]  nfamilies    Number of families
 ** \param[in]  xmax         Maximum horizontal distance
 ** \param[in]  ymax         Maximum vertical distance
 ** \param[in]  deltax       Horizontal dilation
 ** \param[in]  deltay       Vertical dilation
 ** \param[in]  mean         Mean of thickness law
 ** \param[in]  stdev        Standard deviation of thickness law
 **
 *****************************************************************************/
Frac_Environ* fracture_alloc_environ(int nfamilies,
                                                     double xmax,
                                                     double ymax,
                                                     double deltax,
                                                     double deltay,
                                                     double mean,
                                                     double stdev)
{
  Frac_Environ *frac_environ;
  int family;

  /* Allocation */

  frac_environ = new Frac_Environ;
  frac_environ->nfamilies = nfamilies;
  frac_environ->nfaults = 0;
  frac_environ->xmax = xmax;
  frac_environ->ymax = ymax;
  frac_environ->deltax = deltax;
  frac_environ->deltay = deltay;
  frac_environ->xextend = xmax + 2 * deltax;
  frac_environ->mean = mean;
  frac_environ->stdev = stdev;

  /* Family structures */

  frac_environ->frac_fams.resize(nfamilies);
  for (family = 0; family < nfamilies; family++)
  {
    Frac_Fam &frac_family = frac_environ->frac_fams[family];
    frac_family.orient = 0.;
    frac_family.dorient = 0.;
    frac_family.theta0 = 0.;
    frac_family.alpha = 0.;
    frac_family.ratcst = 0.;
    frac_family.prop1 = 0.;
    frac_family.prop2 = 0.;
    frac_family.aterm = 0.;
    frac_family.bterm = 0.;
    frac_family.range = 0.;
  }

  return (frac_environ);
}

/****************************************************************************/
/*!
 **  Deallocate the Frac_Environ structure
 **
 ** \return  Pointer to the Frac_Environ structure
 **
 ** \param[in]  frac_environ Pointer to the Frac_Environ structure
 **
 *****************************************************************************/
Frac_Environ* fracture_dealloc_environ(Frac_Environ *frac_environ)

{
  if (frac_environ == nullptr) return (frac_environ);
  delete frac_environ;
  frac_environ = (Frac_Environ*) nullptr;
  return (frac_environ);
}

/****************************************************************************/
/*!
 **  Add a major Fault to the Frac_Environ structure
 **
 ** \return Rank of the added fault
 **
 ** \param[in]  frac_environ Pointer to the Frac_Environ structure
 ** \param[in]  fault_coord  Abscissae of the first fault point
 ** \param[in]  fault_orient Orientation of the fault
 **
 *****************************************************************************/
int fracture_add_fault(Frac_Environ *frac_environ,
                                       double fault_coord,
                                       double fault_orient)
{
  int nfaults;

  /* Initializations */

  if (frac_environ == nullptr) messageAbort("Frac_Environ must already exist");
  nfaults = frac_environ->nfaults;

  /* Resize the structure */

  frac_environ->frac_faults.resize(nfaults + 1);
  Frac_Fault &fault = frac_environ->frac_faults[nfaults];

  fault.thetal.resize(frac_environ->nfamilies);
  fault.thetar.resize(frac_environ->nfamilies);
  fault.rangel.resize(frac_environ->nfamilies);
  fault.ranger.resize(frac_environ->nfamilies);

  /* Initialize the new Frac_Fault structure */

  fault.coord = fault_coord;
  fault.orient = fault_orient;

  /* Increase the number of faults */

  frac_environ->nfaults++;
  return (nfaults);
}

/****************************************************************************/
/*!
 **  Update Frac_Environ for a family
 **
 ** \param[in]  frac_environ Pointer to the Frac_Environ structure
 ** \param[in]  family       Rank of family
 ** \param[in]  orient       Average fault orientation
 ** \param[in]  dorient      Tolerance for Orientation
 ** \param[in]  theta0       Reference Poisson Intensity
 ** \param[in]  alpha        Intensity law exponent
 ** \param[in]  ratcst       Ratio of constant vs shaped intensity
 ** \param[in]  prop1        Survival constant probability
 ** \param[in]  prop2        Survival length-dependent probability
 ** \param[in]  aterm        Survival cumulated length exponent
 ** \param[in]  bterm        Survival thickness exponent
 ** \param[in]  range        Range of fracture repulsion area
 **
 *****************************************************************************/
void fracture_update_family(Frac_Environ *frac_environ,
                                            int family,
                                            double orient,
                                            double dorient,
                                            double theta0,
                                            double alpha,
                                            double ratcst,
                                            double prop1,
                                            double prop2,
                                            double aterm,
                                            double bterm,
                                            double range)
{
  if (family < 0 || family >= frac_environ->nfamilies) return;
  Frac_Fam &frac_family = frac_environ->frac_fams[family];
  frac_family.orient = orient;
  frac_family.dorient = dorient;
  frac_family.theta0 = theta0;
  frac_family.alpha = alpha;
  frac_family.ratcst = ratcst;
  frac_family.prop1 = prop1;
  frac_family.prop2 = prop2;
  frac_family.aterm = aterm;
  frac_family.bterm = bterm;
  frac_family.range = range;
  return;
}

/****************************************************************************/
/*!
 **  Update Frac_Fault for a family
 **
 ** \param[in]  frac_environ Pointer to the Frac_Environ structure
 ** \param[in]  ifault       Rank of the fault
 ** \param[in]  family       Rank of family
 ** \param[in]  fault_thetal Maximum intensity on the left
 ** \param[in]  fault_thetar Maximum intensity on the right
 ** \param[in]  fault_rangel Range of density decrease on the left
 ** \param[in]  fault_ranger Range of density decrease on the right
 **
 *****************************************************************************/
void fracture_update_fault(Frac_Environ *frac_environ,
                                           int ifault,
                                           int family,
                                           double fault_thetal,
                                           double fault_thetar,
                                           double fault_rangel,
                                           double fault_ranger)
{
  if (ifault < 0 || ifault >= frac_environ->nfaults) return;
  Frac_Fault &fault = frac_environ->frac_faults[ifault];
  fault.thetal[family] = fault_thetal;
  fault.thetar[family] = fault_thetar;
  fault.rangel[family] = fault_rangel;
  fault.ranger[family] = fault_ranger;
}

/****************************************************************************/
/*!
 **  Print a Fracture
 **
 ** \param[in]  frac_environ Frac_Environ structure
 **
 *****************************************************************************/
void fracture_print(Frac_Environ *frac_environ)

{
  if (frac_environ == nullptr) return;

  /* General characteristics */

  mestitle(0, "Geometry");
  message("Field extension (horizontal)    = %lf\n", frac_environ->xmax);
  message("Field extension (vertical)      = %lf\n", frac_environ->ymax);
  message("Field dilation (horizontal)     = %lf\n", frac_environ->deltax);
  message("Field dilation (vertical)       = %lf\n", frac_environ->deltay);
  message("Mean of thickness law           = %lf\n", frac_environ->mean);
  message("St. dev. of thickness law       = %lf\n", frac_environ->stdev);
  message("Number of families              = %d\n", frac_environ->nfamilies);
  message("Number of faults                = %d\n", frac_environ->nfaults);

  /* Loop on the families */

  for (int j = 0; j < frac_environ->nfamilies; j++)
  {
    mestitle(2, "Family #%d/%d", j + 1, frac_environ->nfamilies);
    const Frac_Fam &family = frac_environ->frac_fams[j];
    message("Average Fault Orientation       = %lf deg\n", family.orient);
    message("Tolerance for Orientation       = %lf deg\n", family.dorient);
    message("Reference Poisson Intensity     = %lf\n", family.theta0);
    message("Intensity from thick. exponent  = %lf\n", family.alpha);
    message("Intensity Constant/Shaped ratio = %lf\n", family.ratcst);
    message("Survival constant probability   = %lf\n", family.prop1);
    message("Survival length-dependent proba = %lf\n", family.prop2);
    message("Survival cumul. length exponent = %lf\n", family.aterm);
    message("Survival thickness exponent     = %lf\n", family.bterm);
    message("Fracture repulsion area Range   = %lf\n", family.range);
  }

  /* Loop on the faults */

  for (int i = 0; i < frac_environ->nfaults; i++)
  {
    mestitle(1, "Fault #%d/%d", i + 1, frac_environ->nfaults);
    const Frac_Fault &fault = frac_environ->frac_faults[i];
    message("Location of the Fault           = %lf\n", fault.coord);
    message("Fault orientation               = %lf deg\n", fault.orient);

    for (int j = 0; j < frac_environ->nfamilies; j++)
    {
      mestitle(2, "Fault #%d - Family #%d/%d", i + 1, j + 1,
               frac_environ->nfamilies);
      message("Intensity maximum value (left)  = %lf\n", fault.thetal[j]);
      message("Intensity range (left)          = %lf\n", fault.rangel[j]);
      message("Intensity maximum value (right) = %lf\n", fault.thetar[j]);
      message("Intensity range (right)         = %lf\n", fault.ranger[j]);
    }
  }
  return;
}

/****************************************************************************/
/*!
 **  Calculate the abscissae of a fault at a given elevation
 **
 ** \return The fault abscissae
 **
 ** \param[in]  fault  Frac_Fault structure
 ** \param[in]  cote   Ordinate of the fracture starting point
 **
 *****************************************************************************/
static double st_fault_abscissae(const Frac_Fault &fault, double cote)
{
  return (fault.coord + cote * tan(ut_deg2rad(fault.orient)));
}

/****************************************************************************/
/*!
 **  Calculate the attenuation according to a cubic covariance shape
 **
 ** \return Returned value
 **
 ** \param[in]  h      Normalized distance
 **
 *****************************************************************************/
static double st_cubic(double h)
{
  double value, h2;

  if (h >= 1) return (0.);

  h2 = h * h;
  value = 1. - h2 * (7. - h * (35. / 4. - h2 * (7. / 2. - 3. / 4. * h2)));

  return (value);
}

/****************************************************************************/
/*!
 **  Project the contribution of a density attached to a fault
 **
 ** \return Density shape due to the Fault
 **
 ** \param[in]  fault  Frac_Fault structure
 ** \param[in]  side   -1 for left and +1 for right
 ** \param[in]  family Rank of the family
 ** \param[in]  cote   Ordinate of the fracture starting point
 ** \param[in]  xx     Discretization abscissae
 **
 ** \remarks The density shape function is analoguous to the Cubic Covariance
 **
 *****************************************************************************/
static double st_density_update(const Frac_Fault &fault,
                                int side,
                                int family,
                                double cote,
                                double xx)
{
  double x0, dist, range, theta, h, value;

  x0 = st_fault_abscissae(fault, cote);
  if (side < 0)
  {
    if (xx > x0) return (0.);
    range = fault.rangel[family];
    theta = fault.thetal[family];
  }
  else
  {
    if (xx < x0) return (0.);
    range = fault.ranger[family];
    theta = fault.thetar[family];
  }

  dist = ABS(x0 - xx);
  h = dist / range;
  value = st_cubic(h);
  return (value * theta);
}

/****************************************************************************/
/*!
 **  Check that the current discretization point is not separated from the
 **  current fault by other main faults
 **
 ** \return The discretized point is not separated by another fault
 **
 ** \param[in]  frac_environ Frac_Environ structure
 ** \param[in]  ifault0      Rank of the target fault
 ** \param[in]  x0           Location of the discretized point
 **
 *****************************************************************************/
static int st_same_fault_side(Frac_Environ *frac_environ,
                              int ifault0,
                              double x0)
{
  double x1, x2;
  int ifault;

  x1 = frac_environ->frac_faults[ifault0].coord;
  for (ifault = 0; ifault < frac_environ->nfaults; ifault++)
  {
    if (ifault == ifault0) continue;
    x2 = frac_environ->frac_faults[ifault].coord;

    if ((x0 - x2) * (x2 - x1) >= 0) return (0);
  }
  return (1);
}

/****************************************************************************/
/*!
 **  Calculate the cumulated density
 **
 ** \returns 1 if there is no more room for fractures
 **
 ** \param[in]  title        Title of the printout
 ** \param[in]  flag_print   1 if the message must be printed
 ** \param[in]  denstab      Discretized density array
 **
 ** \param[out] totdens      Cumulated density
 **
 *****************************************************************************/
static int st_density_cumulate(const char *title,
                               int flag_print,
                               double *denstab,
                               double *totdens)
{
  int flag_stop, idisc;
  double total;

  total = 0.;
  for (idisc = 0; idisc < NDISC; idisc++)
    total += denstab[idisc];

  (*totdens) = total;
  flag_stop = (total <= LOW1 * NDISC);

  if (VERBOSE)
  {
    if (flag_stop)
      message("Fracture simulation interrupted: no more room available\n");
    else
    {
      if (flag_print)
        message("- Cumulated Distribution: %s = %lf\n", title, total);
    }
  }
  return (flag_stop);
}

/****************************************************************************/
/*!
 **  Generate the density along the scene
 **
 ** \param[in]  frac_environ Frac_Environ structure
 ** \param[in]  family       Rank of the family
 ** \param[in]  cote         Ordinate of the fracture starting point
 **
 ** \param[out] denstab      Discretization density array
 **
 *****************************************************************************/
static void st_generate_density(Frac_Environ *frac_environ,
                                int family,
                                double cote,
                                double *denstab)
{
  int idisc, ifault;
  double x0, value, ratcst, total;

  /* Initialize the discretization grid */

  ratcst = frac_environ->frac_fams[family].ratcst;
  for (idisc = 0; idisc < NDISC; idisc++)
    denstab[idisc] = 0.;

  /* Create the shaped density */

  if (ratcst < 1.)
  {

    /* Loop on the discretization steps */

    for (idisc = 0; idisc < NDISC; idisc++)
    {
      x0 = XORIGIN + STEP * (idisc + 0.5);

      /* Loop on the faults */

      value = 0.;
      for (ifault = 0; ifault < frac_environ->nfaults; ifault++)
      {
        Frac_Fault &fault = frac_environ->frac_faults[ifault];
        if (!st_same_fault_side(frac_environ, ifault, x0)) continue;

        /* Left side */
        value = MAX(value, st_density_update(fault, -1, family, cote, x0));

        /* Right side */
        value = MAX(value, st_density_update(fault, +1, family, cote, x0));
      }
      denstab[idisc] = value;
    }
  }

  /* Add the constant and the shaped factors */

  for (idisc = 0; idisc < NDISC; idisc++)
    denstab[idisc] = denstab[idisc] * (1. - ratcst) + (1. / NDISC) * ratcst;

  if (VERBOSE)
    (void) st_density_cumulate("Main Fault         ", 1, denstab, &total);
  return;
}

/****************************************************************************/
/*!
 **  Get the rank of the discretized cell
 **
 ** \return Rank of the discretized cell where the cumulated density reaches
 ** \return the target density
 **
 ** \param[in]  cumdens      Target cumulated density
 ** \param[in]  denstab      Discretized density array
 **
 *****************************************************************************/
static int st_get_discretized_rank(double cumdens, double *denstab)
{
  double local;
  int idisc;

  local = 0.;
  for (idisc = 0; idisc < NDISC; idisc++)
  {
    local += denstab[idisc];
    if (local > cumdens) return (idisc);
  }
  return (NDISC - 1);
}

/****************************************************************************/
/*!
 **  Search for the segment of a fracture belonging to the current layer
 **
 ** \return 1 if the segment belongs to the layer; 0 otherwise
 **
 ** \param[in]  desc         Frac_Desc structure
 ** \param[in]  cote         Ordinate of the fracture starting point
 **
 ** \param[in]  xd           Abscissae of the first end-point
 ** \param[in]  yd           Ordinate of the first end-point
 ** \param[in]  xe           Abscissae of the second end-point
 ** \param[in]  ye           Ordinate of the second end-point
 **
 *****************************************************************************/
static int st_belong_to_layer(Frac_Desc &desc,
                              double cote,
                              double *xd,
                              double *yd,
                              double *xe,
                              double *ye)
{
  int i;

  for (i = 0; i < desc.npoint - 1; i++)
  {
    if (ABS(YYF(i) - cote) > EPS) continue;
    *xd = XXF(i);
    *yd = YYF(i);
    *xe = XXF(i + 1);
    *ye = YYF(i + 1);
    return (1);
  }
  return (0);
}

/****************************************************************************/
/*!
 **  Check fracture against the other ones and possibly intersect the current one
 **
 ** \param[in]  frac_list    Frac_List structure
 ** \param[in]  cote         Ordinate of the fracture starting point
 ** \param[in]  ifrac0       Rank of the current fracture
 **
 *****************************************************************************/
static void st_check_fracture_intersect(Frac_List *frac_list,
                                        double cote,
                                        int ifrac0)
{
  double xd1, xe1, yd1, ye1, xd2, xe2, yd2, ye2, x, y;

  /* Initializations */

  if (!FLAG_CHECK) return;
  Frac_Desc &desc = frac_list->frac_descs[ifrac0];
  xd1 = XXF(desc.npoint - 2);
  yd1 = YYF(desc.npoint - 2);
  xe1 = XXF(desc.npoint - 1);
  ye1 = YYF(desc.npoint - 1);

  /* Check against the previous fractures */

  for (int ifrac = 0; ifrac < frac_list->nfracs; ifrac++)
  {
    if (ifrac == ifrac0) continue;
    Frac_Desc &desc2 = frac_list->frac_descs[ifrac];

    /* Look for the segment belonging to the current layer */

    if (!st_belong_to_layer(desc2, cote, &xd2, &yd2, &xe2, &ye2)) continue;

    /* Perform the intersection between two segments */

    if (segment_intersect(xd1, yd1, xe1, ye1, xd2, yd2, xe2, ye2, &x, &y))
      continue;

    /* Update the endpoint in case of intersection */

    XXF(desc.npoint - 1) = x;
    YYF(desc.npoint - 1) = y;
  }
}

/****************************************************************************/
/*!
 **  Add a fracture
 **
 ** \return Rank of the fracture
 **
 ** \param[in]  ifrac        Rank of the fracture (if old) or -1 for a new one
 ** \param[in]  family       Rank of the family
 ** \param[in]  xx           Abscissae of the starting point
 ** \param[in]  cote         Ordinate of the fracture starting point
 ** \param[in]  thick        Thickness of the layer
 ** \param[in]  orient       Orientation
 ** \param[in,out] frac_list Frac_List structure
 **
 ** \param[out] xp           Abscissae of the ending point
 **
 *****************************************************************************/
static int st_frac_add(int ifrac,
                       int family,
                       double xx,
                       double cote,
                       double thick,
                       double orient,
                       Frac_List *frac_list,
                       double *xp)
{
  int number;

  /* Create a new Fracture */

  if (ifrac < 0)
  {
    frac_list = fracture_manage_list(1, frac_list);
    ifrac = frac_list->nfracs - 1;
  }
  Frac_Desc &desc = frac_list->frac_descs[ifrac];

  /* Set the attributes of the new segment */

  number = (desc.npoint == 0) ? 2 :
                                desc.npoint + 1;
  desc.xy.resize(2 * number);
  if (desc.npoint == 0)
  {
    desc.family = family;
    desc.orient = orient;
    XXF(desc.npoint) = xx;
    YYF(desc.npoint) = cote;
    desc.npoint++;
  }
  desc.family = family;
  desc.orient = orient;
  XXF(desc.npoint) = *xp = thick * tan(ut_deg2rad(orient)) + xx;
  YYF(desc.npoint) = thick + cote;
  desc.npoint++;

  st_check_fracture_intersect(frac_list, cote, ifrac);

  return (ifrac);
}

/****************************************************************************/
/*!
 **  Update the discretized intensity array to account for repulsion
 **
 ** \param[in]  x0           Abscissae of the fracture
 ** \param[in]  range        Range for the repulsion
 ** \param[in,out] denstab   Discretized density array
 **
 ** \remarks  We first designate the central cell of the discretized density
 ** \remarks  array where the new fault is located
 ** \remarks  The density of the central cell is set to LOW0
 ** \remarks  The density of the peripheral cells is set to LOW1
 ** \remarks  A cell is peripheral if located in the 'range' of the central one
 **
 *****************************************************************************/
static void st_update_repulsion(double x0, double range, double *denstab)
{
  int idisc, i, idisc0, nrange;

  /* Initializations */

  idisc0 = (int) ((x0 - XORIGIN) / STEP);
  nrange = MAX(1, (int ) (range / STEP));

  /* Central cell */
  if (VALID_DISC(idisc0)) denstab[idisc0] = LOW0;

  /* Peripheral cell on the left */
  for (i = 0; i < nrange; i++)
  {
    idisc = idisc0 - i - 1;
    if (VALID_DISC(idisc)) denstab[idisc] = LOW1;
  }

  /* Peripheral cell on the right */
  for (i = 0; i < nrange; i++)
  {
    idisc = idisc0 + i + 1;
    if (VALID_DISC(idisc)) denstab[idisc] = LOW1;
  }
}

/****************************************************************************/
/*!
 **  Simulate the fractures
 **
 ** \return  Number of fracture created in this layer
 **
 ** \param[in]  frac_environ Frac_Environ structure
 ** \param[in]  family       Rank of the family
 ** \param[in]  cote         Ordinate of the fracture starting point
 ** \param[in]  thick        Thickness of the layer
 ** \param[in]  theta        Intensity of the layer
 ** \param[in]  denstab      Discretized density array
 ** \param[in,out] frac_list Frac_List structure
 **
 *******************************************************************F**********/
static int st_simulate_fractures(Frac_Environ *frac_environ,
                                 int family,
                                 double cote,
                                 double thick,
                                 double theta,
                                 double *denstab,
                                 Frac_List *frac_list)
{
  int ifrac, idisc, nfracs, neff;
  double orient, dorient, cumdens, xdeb, xfin, xx, xp, angle, total;

  /* Initializations */

  Frac_Fam &frac_fam = frac_environ->frac_fams[family];
  orient = frac_fam.orient;
  dorient = frac_fam.dorient;
  nfracs = law_poisson(theta * frac_environ->xextend);

  /* Loop on the fractures to be simulated */

  for (ifrac = neff = 0; ifrac < nfracs; ifrac++)
  {

    /* Calculate cumulative density */

    if (st_density_cumulate("Fracture generation", 0, denstab, &total)) break;

    /* Simulate the location along to the regionalized density */

    cumdens = total * law_uniform(0., 1.);

    /* Get the rank of the discretized cell */

    idisc = st_get_discretized_rank(cumdens, denstab);
    xdeb = XORIGIN + STEP * (idisc);
    xfin = XORIGIN + STEP * (idisc + 1);
    xx = law_uniform(xdeb, xfin);
    angle = law_uniform(MAX(-90., orient - dorient),
                        MIN(+90., orient + dorient));

    /* Add a new fracture */

    (void) st_frac_add(-1, family + 1, xx, cote, thick, angle, frac_list, &xp);

    /* Update the density to account for repulsion */

    st_update_repulsion(xx, frac_fam.range, denstab);
    neff++;
  }

  if (VERBOSE)
  {
    message("New Simulated fractures in Bench  = %d \n", neff);
    message("Total count of fractures          = %d \n", frac_list->nfracs);
  }
  return (nfracs);
}

/****************************************************************************/
/*!
 **  Calculate the fracture extension
 **
 ** \return  The fracture extension
 **
 ** \param[in]  desc         Current fracture description
 ** \param[in]  cote         Selected layer or TEST for all layers
 ** \param[in]  dcote        Tolerance on the layer elevation
 **
 *****************************************************************************/
static double st_fracture_extension(Frac_Desc &desc, double cote, double dcote)
{
  double dist, distx, disty;
  int i;

  dist = 0.;
  for (i = 0; i < desc.npoint - 1; i++)
  {
    distx = XXF(i+1) - XXF(i);
    disty = YYF(i+1) - YYF(i);
    if (!FFFF(cote) && (YYF(i) < cote - dcote || YYF(i+1) < cote - dcote))
      continue;
    dist += sqrt(distx * distx + disty * disty);
  }
  return (dist);
}

/******************************************************************F**********/
/*!
 **  Check if the fracture must be interrupted
 **
 ** \return  1 if the fracture must be interrupted
 **
 ** \param[in]  frac_fam     Frac_Fam structure
 ** \param[in]  desc         Current fracture description
 ** \param[in]  thick        Thickness of the current layer
 **
 *****************************************************************************/
static int st_fracture_interrupt(Frac_Fam &frac_fam,
                                 Frac_Desc &desc,
                                 double thick)
{
  double rndval, dist, prop1, prop2, expa, expb, proba;

  /* Constant probability */
  prop1 = frac_fam.prop1;

  /* Length dependent probability */
  prop2 = frac_fam.prop2;

  /* Cumulated fracture length dependent exponent */
  dist = st_fracture_extension(desc, TEST, TEST);
  expa = exp(-dist / frac_fam.aterm);

  /* Thickness dependent exponent */
  expb = exp(-thick / frac_fam.bterm);

  /* Final probability */
  proba = (prop1 + prop2 * (1. - expa)) * expb;

  rndval = law_uniform(0., 1);
  return (rndval > proba);
}

/****************************************************************************/
/*!
 **  Correct density due to fractures of previous families (if any)
 **
 ** \param[in]  frac_environ Frac_Environ structure
 ** \param[in]  family       Rank of the family (starting from 0)
 ** \param[in]  cote         Ordinate of the fracture starting point
 ** \param[in]  denstab      Discretization density array
 ** \param[in,out]  frac_list Frac_List structure
 **
 *****************************************************************************/
static void st_correct_density(Frac_Environ *frac_environ,
                               int family,
                               double cote,
                               double *denstab,
                               Frac_List *frac_list)
{
  int i, ifrac;
  double total;

  /* Initializations */

  Frac_Fam &frac_fam = frac_environ->frac_fams[family];

  /* Loop on the fractures */

  for (ifrac = 0; ifrac < frac_list->nfracs; ifrac++)
  {
    Frac_Desc &desc = frac_list->frac_descs[ifrac];
    if (desc.family >= family + 1) continue;

    /* Loop on the segments of the fault */

    for (i = 0; i < desc.npoint - 1; i++)
    {
      if (ABS(YYF(i) - cote) <= EPS)
        st_update_repulsion(XXF(i), frac_fam.range, denstab);
      if (ABS(YYF(i+1) - cote) <= EPS)
        st_update_repulsion(XXF(i + 1), frac_fam.range, denstab);
    }
  }

  if (VERBOSE)
    (void) st_density_cumulate("Previous families  ", 1, denstab, &total);
}

/****************************************************************************/
/*!
 **  Extend already existing fractures
 **
 ** \return The proportion of survival
 **
 ** \param[in]  frac_environ  Frac_Environ structure
 ** \param[in]  family        Rank of the family
 ** \param[in]  cote          Ordinate of the fracture starting point
 ** \param[in]  thick         Thickness of the layer
 ** \param[in]  denstab       Discretization density array
 ** \param[in,out]  frac_list Frac_List structure
 **
 *****************************************************************************/
static double st_extend_fractures(Frac_Environ *frac_environ,
                                  int family,
                                  double cote,
                                  double thick,
                                  double *denstab,
                                  Frac_List *frac_list)
{
  int ifrac, next, ntotal;
  double xx, xp, angle, propsur, total;

  /* Initializations */

  Frac_Fam &frac_fam = frac_environ->frac_fams[family];

  /* Loop on the already existing fractures */

  next = ntotal = 0;
  for (ifrac = 0; ifrac < frac_list->nfracs; ifrac++)
  {
    Frac_Desc &desc = frac_list->frac_descs[ifrac];
    if (desc.family != family + 1) continue;
    if (ABS(YYF(desc.npoint-1) - cote) > EPS) continue;
    angle = desc.orient;

    /* Evaluate the persistence probability */

    ntotal++;
    if (st_fracture_interrupt(frac_fam, desc, thick)) continue;

    /* The fracture is extended */

    next++;
    xx = XXF(desc.npoint - 1);
    (void) st_frac_add(ifrac, desc.family, xx, cote, thick, angle, frac_list,
                       &xp);

    /* Update the density to account for repulsion */

    st_update_repulsion(xx, frac_fam.range, denstab);
  }

  propsur = (ntotal > 0) ? (double) next / (double) ntotal :
                           0.;
  if (VERBOSE && ntotal > 0)
  {
    message("Survival proportion               = %lf (%d from %d)\n", propsur,
            next, ntotal);
    (void) st_density_cumulate("After survival     ", 1, denstab, &total);
  }
  return (propsur);
}

/****************************************************************************/
/*!
 **  Derive the intensity of the current layer
 **
 ** \return Layer intensity
 **
 ** \param[in]  frac_environ Pointer to the Frac_Environ structure
 ** \param[in]  family       Rank of family
 ** \param[in]  thick        Layer thickness
 **
 *****************************************************************************/
static double st_layer_intensity(Frac_Environ *frac_environ,
                                 int family,
                                 double thick)
{
  Frac_Fam &frac_fam = frac_environ->frac_fams[family];
  double theta1 = frac_fam.theta0 / pow(thick, frac_fam.alpha);

  if (VERBOSE) message("Initial Intensity                 = %lf\n", theta1);

  return (theta1);
}

/****************************************************************************/
/*!
 **  Derive the intensity of the current layer given the survival
 **
 ** \return Layer actual intensity
 **
 ** \param[in]  theta1       Apparent layer intensity
 ** \param[in]  thetap       Intensity of the previous layer
 ** \param[in]  propsur      Actual survival proportion
 **
 *****************************************************************************/
static double st_derive_intensity(double theta1, double thetap, double propsur)
{
  double theta2 = theta1;

  if (thetap > 0) theta2 = MAX(0., theta1 - propsur * thetap);

  if (VERBOSE) message("Intensity corrected from survival = %lf\n", theta2);

  return (theta2);
}

/****************************************************************************/
/*!
 **  Generate the series of layers
 **
 ** \param[in,out] frac_environ Pointer to the Frac_Environ structure
 ** \param[in]  nlayers      Number of layers
 ** \param[in]  y0           Ordinate of the origin
 **
 ** \param[out] thicks       Array of layer thicknesses
 **
 ** \remarks The origin is calculated so that a layer edge is located at 0.
 **
 *****************************************************************************/
static void st_layers_manage(Frac_Environ *frac_environ,
                             int *nlayers,
                             double *y0,
                             double **thicks)
{
  double *thicks_loc;
  double mean, stdev, var, total, thick, cote_up, cote_down, thick_min, theta,
      k;
  int number, nquant, flag_rnd;
  int quantum = 10;

  /* Initializations */

  thicks_loc = *thicks;
  thick_min = frac_environ->mean / 10.;
  mean = frac_environ->mean;
  stdev = frac_environ->stdev;
  var = stdev * stdev;
  theta = var / mean;
  flag_rnd = var > 0.;
  k = (flag_rnd) ? (mean * mean) / var :
                   0.;

  /* Allocation */

  number = nquant = 0;

  /* Downwards */

  total = 0.;
  while (total < frac_environ->deltay)
  {
    thick = (flag_rnd) ? theta * law_gamma(k) :
                         mean;
    if (thick < thick_min) continue;
    total += thick;
    if (number >= nquant * quantum)
    {
      nquant++;
      thicks_loc = (double*) mem_realloc((char* ) thicks_loc,
                                         sizeof(double) * nquant * quantum, 1);
    }
    thicks_loc[number++] = thick;
  }
  cote_down = -total;

  /* Upwards */

  total = 0.;
  while (total < frac_environ->ymax)
  {
    thick = (flag_rnd) ? theta * law_gamma(k) :
                         mean;
    if (thick < thick_min) continue;
    total += thick;
    if (number >= nquant * quantum)
    {
      nquant++;
      thicks_loc = (double*) mem_realloc((char* ) thicks_loc,
                                         sizeof(double) * nquant * quantum, 1);
    }
    thicks_loc[number++] = thick;
  }
  cote_up = total;

  /* Resize */

  thicks_loc = (double*) mem_realloc((char* ) thicks_loc,
                                     sizeof(double) * number, 1);
  *y0 = cote_down;
  *nlayers = number;

  /* Optional information */

  if (VERBOSE)
  {
    mestitle(0, "Layer generation");
    message("Thickness law - Mean               = %lf\n", mean);
    message("Thickness law - St. Dev.           = %lf\n", stdev);
    message("Minimum simulated level            = %lf\n", cote_down);
    message("Maximum simulated level            = %lf\n", cote_up);
    message("Number of layers                   = %d \n", number);
  }

  *thicks = thicks_loc;
  return;
}

/****************************************************************************/
/*!
 **  Read the series of layers
 **
 ** \param[in]  nlayers_in   Number of elevations to be read
 ** \param[in]  elevations   Array of elevations to be read
 **
 ** \param[out] nlayers      Number of layers
 ** \param[out] y0           Ordinate of the origin
 ** \param[out] thicks       Array of layer thicknesses
 **
 *****************************************************************************/
static void st_layers_read(int nlayers_in,
                           double *elevations,
                           int *nlayers,
                           double *y0,
                           double **thicks)
{
  double *thicks_loc;
  double cote, cote_down, cote_up;

  /* Initializations */

  thicks_loc = (double*) mem_alloc(sizeof(double) * (nlayers_in - 1), 1);

  cote = elevations[0];
  for (int i = 1; i < nlayers_in; i++)
  {
    thicks_loc[i - 1] = elevations[i] - cote;
    cote = elevations[i];
  }
  cote_down = (*y0) = elevations[0];
  cote_up = elevations[nlayers_in - 1];
  *nlayers = nlayers_in - 1;

  /* Optional information */

  if (VERBOSE)
  {
    mestitle(0, "Layer (read)");
    message("Minimum simulated level            = %lf\n", cote_down);
    message("Maximum simulated level            = %lf\n", cote_up);
    message("Number of layers                   = %d \n", nlayers_in);
  }

  *thicks = thicks_loc;
  return;
}

/****************************************************************************/
/*!
 **  Simulate a set of fractures
 **
 ** \param[in]  frac_environ   Pointer to the Frac_Environ structure
 ** \param[in]  flag_sim_layer TRUE for simulating layers
 **                            FALSE if they are read
 ** \param[in]  flag_sim_fract TRUE for simulating the fractures
 **                            FALSE establish the layers and the main faults
 ** \param[in]  nlayers_in     Number of input layers (used if flag_sim_layer=F)
 ** \param[in]  elevations     Array of elevations (used if flag_sim_layer=F)
 ** \param[in]  seed           Seed for the random number generator
 ** \param[in]  verbose        Verbose option
 ** \param[in,out] frac_list   List of the previously defined fractures
 **                            (NULL for the first usage)
 **                            At end, list of the finally defined fractures
 **
 ** \param[out]  nlayers       Number of layers
 ** \param[out]  ninfo_arg     Number of information per layer
 ** \param[out]  layinfo       Array of layer information
 **                            Dimension: ninfos * (nlayers+1)
 **
 ** \remark The number of discretization steps used to establish the fracture
 ** \remark density can be defined using:
 ** \remark set_keypair("Fracture_Discretization_Count",newval)
 **
 ** \remark The option for checking the fracture intersect or not (0 or 1)
 ** \remark is defined using set_keypair("Fracture_Check_Intersect",newval)
 **
 ** \remark The lower value used for checking the fracture repulsion
 ** \remark for the central cell is defined using:
 ** \remark set_keypair("Fracture_Repulsion_Low0",newval)
 **
 ** \remark The lower value used for checking the fracture repulsion
 ** \remark for the peripheral cells is defined using:
 ** \remark set_keypair("Fracture_Repulsion_Low1",newval)
 **
 *****************************************************************************/
int fracture_simulate(Frac_Environ *frac_environ,
                                      int flag_sim_layer,
                                      int flag_sim_fract,
                                      int seed,
                                      int verbose,
                                      int nlayers_in,
                                      double *elevations,
                                      int *nlayers,
                                      int *ninfo_arg,
                                      double **layinfo,
                                      Frac_List *frac_list)
{
  double *denstab, *thicks, *locinfo;
  double cote, thick, thickp, theta1, theta2, thetap, xx, y0, angle, propsur;
  int ilayer, ifault, family, ifrac, nfracs, i, ninfos, np1;

  /* Initializations */

  NDISC = (int) get_keypone("Fracture_Discretization_Count", NDISC_DEF);
  FLAG_CHECK = (int) get_keypone("Fracture_Check_Intersect", FLAG_CHECK_DEF);
  LOW0 = get_keypone("Fracture_Repulsion_Low0", LOW0_DEF);
  LOW1 = get_keypone("Fracture_Repulsion_Low1", LOW1_DEF);

  STEP = frac_environ->xextend / NDISC;
  XORIGIN = -frac_environ->deltax;
  VERBOSE = verbose;
  y0 = xx = 0.;
  thicks = denstab = nullptr;
  ninfos = *ninfo_arg = 1 + frac_environ->nfamilies * NPART;
  law_set_random_seed(seed);

  /* Keypair options echoed in the verbose case */

  if (VERBOSE)
  {
    message("Options set by the keypair mechanism:\n");
    message("Fracture_Discretization_Count = %d \n", NDISC);
    message("Fracture_Check_Intersect      = %d \n", FLAG_CHECK);
    message("Fracture_Repulsion_Low0       = %lg\n", LOW0);
    message("Fracture_Repulsion_Low1       = %lg\n", LOW1);
  }

  /* Determine the layers */

  if (flag_sim_layer)
    st_layers_manage(frac_environ, nlayers, &y0, &thicks);
  else
    st_layers_read(nlayers_in, elevations, nlayers, &y0, &thicks);

  /* Core allocation */

  np1 = (*nlayers) + 1;
  denstab = (double*) mem_alloc(sizeof(double) * NDISC, 1);
  locinfo = (double*) mem_alloc(sizeof(double) * np1 * ninfos, 1);
  for (i = 0; i < ninfos * np1; i++)
    locinfo[i] = 0.;

  /* Define the layers */

  cote = y0;
  for (ilayer = 0; ilayer < (*nlayers); ilayer++)
  {
    thick = thicks[ilayer];
    MEM_LAYER(ilayer) = cote;
    cote += thick;
  }
  MEM_LAYER(*nlayers) = cote;

  /* Define the main faults */

  for (ifault = 0; ifault < frac_environ->nfaults; ifault++)
  {
    Frac_Fault &fault = frac_environ->frac_faults[ifault];
    angle = fault.orient;

    /* Loop on the layers */

    ifrac = -1;
    cote = y0;
    xx = st_fault_abscissae(fault, cote);
    for (ilayer = 0; ilayer < (*nlayers); ilayer++)
    {
      thick = thicks[ilayer];
      ifrac = st_frac_add(ifrac, 0, xx, cote, thick, angle, frac_list, &xx);
      cote += thick;
    }
  }
  if (VERBOSE)
  {
    message("Total number of main faults        = %d \n", frac_list->nfracs);
  }

  if (!flag_sim_fract) goto label_suite;

  /* Loop on the different fault famillies */

  for (family = 0; family < frac_environ->nfamilies; family++)
  {
    if (VERBOSE)
      mestitle(0, "Processing Family #%d/%d", family + 1,
               frac_environ->nfamilies);

    /* Loop on the layers */

    cote = y0;
    thetap = thickp = 0.;
    for (ilayer = 0; ilayer < (*nlayers); ilayer++)
    {
      thick = thicks[ilayer];
      if (VERBOSE)
      {
        mestitle(1, "Processing Layer #%d/%d", ilayer + 1, (*nlayers));
        message("Elevation of the layer bottom     = %lf\n", cote);
        message("Thickness of the layer            = %lf\n", thick);
      }

      /* Derive the layer intensity */

      theta1 = st_layer_intensity(frac_environ, family, thick);

      /* Generate the density function */

      st_generate_density(frac_environ, family, cote, denstab);

      /* Correct the density due to already existing faults */

      st_correct_density(frac_environ, family, cote, denstab, frac_list);

      /* Extend previously existing fractures */

      propsur = st_extend_fractures(frac_environ, family, cote, thick, denstab,
                                    frac_list);

      /* Derive the layer intensity, given the survival from previous layers */

      theta2 = st_derive_intensity(theta1, thetap, propsur);

      /* Simulate the fractures abscissa */

      nfracs = st_simulate_fractures(frac_environ, family, cote, thick, theta2,
                                     denstab, frac_list);

      /* Shift the ordinate */

      cote += thick;
      thetap = theta1;
      thickp = thick;
      MEM_THETA1(family,ilayer) = theta1;
      MEM_THETA2(family,ilayer) = theta2;
      MEM_PROPSUR(family,ilayer) = propsur;
      MEM_FRAC(family,ilayer) = (double) nfracs;
      MEM_TOTAL(family,ilayer) = (double) frac_list->nfracs;
    }
  }

  label_suite: denstab = (double*) mem_free((char* ) denstab);
  thicks = (double*) mem_free((char* ) thicks);
  *nlayers = (*nlayers) + 1;
  *layinfo = locinfo;
  return (0);
}

/****************************************************************************/
/*!
 **  Fracture input function
 **
 ** \return  Pointer to the newly created Fracture
 **
 ** \param[in] frac_def : Frac_Environ default structure
 **
 *****************************************************************************/
Frac_Environ* fracture_input(Frac_Environ *frac_def)

{
  Frac_Environ *frac_environ;
  int flag_new, flag_def, ifault, nfamilies, family, flag_def_new;
  double coord, thetal, thetar, rangel, ranger, mean, stdev, deltax, deltay,
      range;
  double orient, dorient, xmax, ymax, theta0, prop1, prop2, aterm, bterm, alpha,
      ratcst;

  /* Initializations */

  frac_environ = nullptr;
  flag_def = (frac_def != nullptr);

  /* Ask for the main questions */

  nfamilies = (flag_def) ? frac_def->nfamilies :
                           0;
  nfamilies = _lire_int("Number of fracture families", flag_def, nfamilies, 1,
                        ITEST);

  mestitle(1, "Scene Geometry");
  xmax = (flag_def) ? frac_def->xmax :
                      0.;
  xmax = _lire_double("Field Extension along X-axis", flag_def, xmax, 0., TEST);
  ymax = (flag_def) ? frac_def->ymax :
                      0.;
  ymax = _lire_double("Field Extension along Y-axis", flag_def, ymax, 0., TEST);
  deltax = (flag_def) ? frac_def->deltax :
                        0.;
  deltax = _lire_double("Field Dilation along X-axis ", flag_def, deltax, 0.,
                        TEST);
  deltay = (flag_def) ? frac_def->deltay :
                        0.;
  deltay = _lire_double("Field Dilation along Y-axis ", flag_def, deltay, 0.,
                        TEST);
  mean = (flag_def) ? frac_def->mean :
                      1.;
  mean = _lire_double("Mean of thickness law       ", 1, mean, 0., TEST);
  stdev = (flag_def) ? frac_def->stdev :
                       1.;
  stdev = _lire_double("St. dev. of thickness law   ", 1, stdev, 0., TEST);

  frac_environ = fracture_alloc_environ(nfamilies, xmax, ymax, deltax, deltay,
                                        mean, stdev);

  /* Loop on fracture families */

  mestitle(1, "Family linked parameters");
  for (family = 0; family < frac_environ->nfamilies; family++)
  {
    message("General parameters for family %d\n", family + 1);
    Frac_Fam &fam_def = frac_def->frac_fams[family];

    flag_def_new = (flag_def && family < frac_def->nfamilies);
    orient = (flag_def_new) ? fam_def.orient :
                              0.;
    orient = _lire_double("Mean Fracture Orientation      ", 1, orient, -89.,
                          89.);
    dorient = (flag_def_new) ? fam_def.dorient :
                               0.;
    dorient = _lire_double("Tolerance for Orientation      ", 1, dorient, 0.,
                           45.);
    theta0 = (flag_def_new) ? fam_def.theta0 :
                              1.;
    theta0 = _lire_double("Reference Poisson intensity    ", 1, theta0, 0.,
                          TEST);
    alpha = (flag_def_new) ? fam_def.alpha :
                             1.;
    alpha = _lire_double("Intensity from thick. exponent ", 1, alpha, 0., 1.);
    ratcst = (flag_def_new) ? fam_def.ratcst :
                              0.;
    ratcst = _lire_double("Intensity Constant/Shaped ratio", flag_def, ratcst,
                          0., 1.);
    prop1 = (flag_def_new) ? fam_def.prop1 :
                             0.;
    prop1 = _lire_double("Survival constant probability  ", flag_def, prop1, 0.,
                         1.);
    prop2 = (flag_def_new) ? fam_def.prop2 :
                             0.;
    prop2 = _lire_double("Survival length-dependent proba", flag_def, prop2, 0.,
                         1.);
    aterm = (flag_def_new) ? fam_def.aterm :
                             0.;
    aterm = _lire_double("Survival cumul. length exponent", flag_def, aterm, 0.,
                         TEST);
    bterm = (flag_def_new) ? fam_def.bterm :
                             0.;
    bterm = _lire_double("Survival thickness exponent    ", flag_def, bterm, 0.,
                         TEST);
    range = (flag_def_new) ? fam_def.range :
                             0.;
    range = _lire_double("Fracture repulsion area Range  ", flag_def, range, 0.,
                         TEST);

    fracture_update_family(frac_environ, family, orient, dorient, theta0, alpha,
                           ratcst, prop1, prop2, aterm, bterm, range);
  }

  /* Implicit loop on the faults */

  flag_new = 1;
  while (flag_new)
  {
    message("\n");
    flag_def_new = 0;
    if (flag_def) flag_def_new = frac_environ->nfaults < frac_def->nfaults;
    flag_new = _lire_logical("Add a new Deterministic Fault", flag_def,
                             flag_def_new);
    if (!flag_new) break;

    coord =
        (flag_def_new) ? frac_def->frac_faults[frac_environ->nfaults].coord :
                         0.;
    coord = _lire_double("Fault Location   ", flag_def_new, coord, 0., xmax);
    orient =
        (flag_def_new) ? frac_def->frac_faults[frac_environ->nfaults].orient :
                         0.;
    orient = _lire_double("Fault orientation", flag_def_new, orient, -89., 89.);
    ifault = fracture_add_fault(frac_environ, coord, orient);

    Frac_Fault &fault_def = frac_def->frac_faults[ifault];

    for (family = 0; family < frac_environ->nfamilies; family++)
    {
      message("Fault parameters for family %d\n", family + 1);
      thetal =
          (flag_def_new && family < frac_def->nfamilies) ? fault_def.thetal[family] :
                                                           0.;
      thetal = _lire_double("Intensity max. value (left)    ", flag_def_new,
                            thetal, 0., TEST);
      rangel =
          (flag_def_new && family < frac_def->nfamilies) ? fault_def.rangel[family] :
                                                           0.;
      rangel = _lire_double("Intensity max. distance (left) ", flag_def_new,
                            rangel, 0., TEST);
      thetar =
          (flag_def_new && family < frac_def->nfamilies) ? fault_def.thetar[family] :
                                                           0.;
      thetar = _lire_double("Intensity max. value (right)   ", flag_def_new,
                            thetar, 0., TEST);
      ranger =
          (flag_def_new && family < frac_def->nfamilies) ? fault_def.ranger[family] :
                                                           0.;
      ranger = _lire_double("Intensity max. distance (right)", flag_def_new,
                            ranger, 0., TEST);
      fracture_update_fault(frac_environ, ifault, family, thetal, thetar,
                            rangel, ranger);
    }
  }

  return (frac_environ);
}

/****************************************************************************/
/*!
 **  Return the number of values per Fracture Segment
 **
 *****************************************************************************/
static int st_get_nbyfrac(void)
{
  static int nbyfrac = 7;
  return (nbyfrac);
}

/****************************************************************************/
/*!
 **  Return the number of values per Wellout information
 **
 *****************************************************************************/
static int st_get_nbyout(void)
{
  static int nbyout = 8;
  return (nbyout);
}

/****************************************************************************/
/*!
 **  Export the Fractures
 **
 ** \param[in]  frac_list    : Pointer to the Frac_List structure
 **
 ** \param[out] nfracs_arg  : Number of fractures exported
 ** \param[out] nbyfrac_arg : Number of elements per fracture
 ** \param[out] fracs_arg   : Array of fracture segments
 **
 ** \remarks The allocated array must be freed
 **
 *****************************************************************************/
void fracture_export(Frac_List *frac_list,
                                     int *nfracs_arg,
                                     int *nbyfrac_arg,
                                     double **fracs_arg)
{
  int nfracs, nbyfrac, ifrac, ip, ntotal;
  double *frac_segs;

  /* Initializations */

  frac_segs = nullptr;
  nfracs = 0;
  nbyfrac = st_get_nbyfrac();
  ntotal = st_end_point_count(frac_list);

  /* Core allocation */

  frac_segs = (double*) mem_alloc(sizeof(double) * nbyfrac * ntotal, 1);

  /* Loading the fractures */

  for (ifrac = 0; ifrac < frac_list->nfracs; ifrac++)
  {
    Frac_Desc &desc = frac_list->frac_descs[ifrac];

    for (ip = 0; ip < desc.npoint; ip++)
    {
      FRACSEG(nfracs,0) = ifrac + 1;
      FRACSEG(nfracs,1) = ip + 1;
      FRACSEG(nfracs,2) = desc.family;
      FRACSEG(nfracs,3) = XXF(ip);
      FRACSEG(nfracs,4) = YYF(ip);
      FRACSEG(nfracs,5) = desc.orient;
      FRACSEG(nfracs,6) = (ip == 0);
      nfracs++;
    }
  }

  /* Resizing */

  frac_segs = (double*) mem_realloc((char* ) frac_segs,
                                    sizeof(double) * nbyfrac * nfracs, 1);

  /* Returning arguments */

  *nfracs_arg = nfracs;
  *nbyfrac_arg = nbyfrac;
  *fracs_arg = frac_segs;
  return;
}

/****************************************************************************/
/*!
 **  Import the Fractures
 **
 ** \return Pointer to the Frac_List structure
 **
 ** \param[in]  nval         : Total number of values
 ** \param[in]  frac_segs    : Array of fracture segments
 **
 *****************************************************************************/
Frac_List* fracture_import(int nval, double *frac_segs)
{
  int i, icur, ifam, nbyfrac, flag_new, nseg;
  double xx, yy, orient;
  Frac_List *frac_list;

  /* Preliminary checks */

  frac_list = nullptr;
  nbyfrac = st_get_nbyfrac();
  if (!isMultiple(nval, nbyfrac))
  {
    messerr("The number of values read (%d) should be a multiple of %d", nval,
            nbyfrac);
    return (frac_list);
  }
  nseg = nval / nbyfrac;

  /* Initializations */

  frac_list = fracture_manage_list(0, NULL);

  /* Loop on the end points */

  icur = -1;
  for (i = 0; i < nseg; i++)
  {
    ifam = (int) FRACSEG(i, 2);
    xx = FRACSEG(i, 3);
    yy = FRACSEG(i, 4);
    orient = FRACSEG(i, 5);
    flag_new = (int) FRACSEG(i, 6);

    if (flag_new || icur < 0)
    {
      frac_list = fracture_manage_list(1, frac_list);
      icur++;
    }
    Frac_Desc &desc = frac_list->frac_descs[icur];

    desc.family = ifam;
    desc.orient = orient;
    desc.xy.resize(2 * (desc.npoint + 1));
    desc.xy[2 * desc.npoint] = xx;
    desc.xy[2 * desc.npoint + 1] = yy;
    desc.npoint++;
  }

  return (frac_list);
}

/****************************************************************************/
/*!
 **  Extract the array of fracture lengthes
 **
 ** \return The returned array
 **
 ** \param[in]  frac_list    Frac_List structure
 ** \param[in]  family       Rank of the family or ITEST for all
 ** \param[in]  cote         Selected layer or TEST for all layers
 ** \param[in]  dcote        Tolerance on the layer elevation
 **
 ** \param[out] ntab         Number of values in the returned array
 **
 ** \remark The returned array must be freed by the calling routine
 **
 *****************************************************************************/
double* fracture_extract_length(Frac_List *frac_list,
                                                int family,
                                                double cote,
                                                double dcote,
                                                int *ntab)
{
  int i, ecr;
  double *tab, value;

  /* Core allocation */

  (*ntab) = 0;
  tab = (double*) mem_alloc(sizeof(double) * frac_list->nfracs, 0);
  if (tab == nullptr) return (tab);

  /* Loading the fractures */

  ecr = 0;
  for (i = 0; i < frac_list->nfracs; i++)
  {
    Frac_Desc &desc = frac_list->frac_descs[i];

    /* Selection according to the family criterion */

    if (!IFFFF(family) && desc.family != family) continue;

    /* Calculate the fracture length */

    value = st_fracture_extension(desc, cote, dcote);
    if (value <= 0.) continue;
    tab[ecr++] = value;
  }

  /* Sort the intersections */

  ut_sort_double(0, ecr, NULL, tab);

  /* Rescale memory */

  tab = (double*) mem_realloc((char* ) tab, sizeof(double) * ecr, 0);
  if (tab == nullptr) return (tab);

  (*ntab) = ecr;
  return (tab);
}

/****************************************************************************/
/*!
 **  Extract the fracture interdistances
 **
 ** \return The returned array
 **
 ** \param[in]  frac_list    Frac_List structure
 ** \param[in]  family       Rank of the family or ITEST for all
 ** \param[in]  cote         Selected layer or TEST for all layers
 ** \param[in]  dcote        Tolerance on the layer elevation
 **
 ** \param[out] ntab         Number of values in the returned array
 **
 ** \remark The returned array must be freed by the calling routine
 **
 *****************************************************************************/
double* fracture_extract_dist(Frac_List *frac_list,
                                              int family,
                                              double cote,
                                              double dcote,
                                              int *ntab)
{
  int i, ip, number, ecr, ndist;
  double *tab, value, valref;

  /* Core allocation */

  (*ntab) = 0;
  number = st_end_point_count(frac_list);
  tab = (double*) mem_alloc(sizeof(double) * number, 0);
  if (tab == nullptr) return (tab);

  /* Loading the abcissae of the fracture (with the target layer) */

  ecr = 0;
  for (i = 0; i < frac_list->nfracs; i++)
  {
    Frac_Desc &desc = frac_list->frac_descs[i];

    /* Loop on the end points */

    for (ip = 0; ip < desc.npoint; ip++)
    {

      /* Selection according to the family criterion */

      if (!IFFFF(family) && desc.family != family) continue;

      /* Selection according to the layer */

      if (!FFFF(cote) && !FFFF(dcote) &&
      ABS(cote - YYF(ip)) > dcote) continue;

      tab[ecr++] = XXF(ip);
    }
  }
  ndist = ecr;
  if (ndist <= 0) goto label_end;

  /* Sort the intersections */

  ut_sort_double(0, ndist, NULL, tab);

  /* Calculate the inter-fracture distances */

  valref = tab[0];
  for (i = 1; i < ndist; i++)
  {
    value = tab[i];
    tab[i - 1] = value - valref;
    valref = value;
  }

  /* Sort the fracture inter-distances */

  ndist--;
  ut_sort_double(0, ndist, NULL, tab);

  /* Rescale memory */

  label_end: tab = (double*) mem_realloc((char* ) tab, sizeof(double) * ndist,
                                         0);
  if (tab == nullptr) return (tab);
  (*ntab) = ndist;
  return (tab);
}

/****************************************************************************/
/*!
 **  Plunge a segment in a block Db
 **
 ** \param[in]  dbgrid    Db structure
 ** \param[in]  iptr      Pointer to the Perm variable
 ** \param[in]  indg      Array for grid indices
 ** \param[in]  delta     Increment
 ** \param[in]  value     Value assigned to the fracture segment
 ** \param[in]  x1        First coordinate of the first point
 ** \param[in]  y1        Second coordinate of the first point
 ** \param[in]  x2        First coordinate of the second point
 ** \param[in]  y2        Second coordinate of the second point
 **
 *****************************************************************************/
static void st_plunge_segment(Db *dbgrid,
                              int iptr,
                              int *indg,
                              double delta,
                              double value,
                              double x1,
                              double y1,
                              double x2,
                              double y2)
{
  double deltax, deltay, dist, coor[2];
  int number, i, iech;

  deltax = x2 - x1;
  deltay = y2 - y1;
  dist = sqrt(deltax * deltax + deltay * deltay);
  number = (int) floor(dist / delta);

  // Loop on the discretization lags

  for (i = 0; i <= number; i++)
  {
    coor[0] = x1 + deltax * i / number;
    coor[1] = y1 + deltay * i / number;

    if (point_to_grid(dbgrid, coor, 0, indg) < 0) continue;
    iech = db_index_grid_to_sample(dbgrid, indg);
    dbgrid->updArray(iech, iptr, 5, value);
  }
}

/****************************************************************************/
/*!
 **  Plunge a segment painted with gradually changing permeability in a block Db
 **
 ** \param[in]  dbgrid    Db structure
 ** \param[in]  iptr      Pointer to the Perm variable
 ** \param[in]  indg      Array for grid indices
 ** \param[in]  delta     Increment
 ** \param[in]  traj      Vector describing the trajectory
 ** \param[in]  ntraj     Number of end points of the trajectory
 ** \param[in]  perm1     Permeability at origin of curvilinear distance
 ** \param[in]  perm2     Permeability at range of curvilinear distance
 ** \param[in]  range     Range of the permeability change
 **
 *****************************************************************************/
static void st_plunge_segment_gradual(Db *dbgrid,
                                      int iptr,
                                      int *indg,
                                      double delta,
                                      double *traj,
                                      int ntraj,
                                      double perm1,
                                      double perm2,
                                      double range)
{
  double deltax, deltay, dist, incr, incx, incy, total, x1, y1, x2, y2, perm, h,
      red;
  double coor[2];
  int number, i, iech;

  /* Loop on the segments */

  total = 0.;
  for (int ip = 0; ip < ntraj - 1; ip++)
  {
    x1 = TRAJ(ip, 0);
    y1 = TRAJ(ip, 1);
    x2 = TRAJ(ip + 1, 0);
    y2 = TRAJ(ip + 1, 1);
    deltax = x2 - x1;
    deltay = y2 - y1;
    dist = sqrt(deltax * deltax + deltay * deltay);
    number = (int) floor(dist / delta);
    incx = deltax / number;
    incy = deltay / number;
    incr = dist / number;

    /* Loop on the discretization lags */

    for (i = 0; i <= number; i++)
    {
      coor[0] = x1 + i * incx;
      coor[1] = y1 + i * incy;
      total += incr;

      if (point_to_grid(dbgrid, coor, 0, indg) < 0) continue;
      iech = db_index_grid_to_sample(dbgrid, indg);
      if (iech < 0) continue;

      if (FFFF(range))
        perm = perm2;
      else
      {
        h = total / range;
        red = st_cubic(h);
        perm = perm2 + (perm1 - perm2) * red;
      }
      dbgrid->setArray(iech, iptr, perm);
    }
  }
}

/****************************************************************************/
/*!
 **  Plunge a subset of the simulated fractures into a block
 **
 **  \return Error return code
 **
 ** \param[in,out] dbgrid    Db structure
 ** \param[in]  frac_list    Pointer to the Frac_List structure
 ** \param[in]  locinfo      Array of fracture environment description
 ** \param[in]  n_layers     Number of layers in 'locinfo'
 ** \param[in]  nfamilies    Number of families
 ** \param[in]  xmax         Maximum extension along horizontal axis
 ** \param[in]  permtab      Permabilities per family (starting from 0)
 ** \param[in]  perm_mat     Permability for the matrix
 ** \param[in]  perm_bench   Permability along the bench edge
 ** \param[in]  ndisc        Number of discretization steps
 **
 *****************************************************************************/
int fracture_to_block(Db *dbgrid,
                                      Frac_List *frac_list,
                                      double *locinfo,
                                      int n_layers,
                                      int nfamilies,
                                      double xmax,
                                      double *permtab,
                                      double perm_mat,
                                      double perm_bench,
                                      int ndisc)
{
  int *indg, iptr, error, idim, ilayer, ninfos;
  double dmin, delta, cote;

  /* Initializations */

  error = 1;
  indg = nullptr;
  ninfos = 1 + nfamilies * NPART;

  /* Preliminary checks */

  if (!is_grid(dbgrid))
  {
    messerr("The Db file must be organized as a regular Grid");
    goto label_end;
  }
  if (dbgrid->getNDim() != 2)
  {
    messerr("This application is limited to 2-D grid");
    goto label_end;
  }
  dmin = 1.e30;
  for (idim = 0; idim < dbgrid->getNDim(); idim++)
  {
    if (dbgrid->getDX(idim) < dmin) dmin = dbgrid->getDX(idim);
  }
  delta = dmin / (double) ndisc;

  // Allocate the new variable 

  indg = db_indg_alloc(dbgrid);
  if (indg == nullptr) goto label_end;
  iptr = dbgrid->addFields(1, perm_mat);

  // Plunge the environment

  if (perm_bench > 0)
  {
    for (ilayer = 0; ilayer < n_layers; ilayer++)
    {
      cote = MEM_LAYER(ilayer);
      st_plunge_segment(dbgrid, iptr, indg, delta, perm_bench, 0., cote, xmax,
                        cote);
    }
  }

  /* Loop on the fractures */

  for (int ifrac = 0; ifrac < frac_list->nfracs; ifrac++)
  {
    Frac_Desc &desc = frac_list->frac_descs[ifrac];

    /* Loop on the segments */

    for (int ip = 0; ip < desc.npoint - 1; ip++)
    {
      st_plunge_segment(dbgrid, iptr, indg, delta, permtab[desc.family],
                        XXF(ip), YYF(ip), XXF(ip + 1), YYF(ip + 1));
    }
  }

  /* Set the error return code */

  error = 0;

  label_end: indg = db_indg_free(indg);
  return (error);
}

/****************************************************************************/
/*!
 **  Add an end point to the wellout array
 **
 **  \return New allocated wellout array
 **
 ** \param[in]  wellout      Array provided as input
 ** \param[in]  nbyout       Number of attribute by end point
 ** \param[in]  x            Ascissae of the intersection
 ** \param[in]  y            Ordinate of the intersection
 ** \param[in]  ifrac        Rank of the fracture (starting from 1)
 ** \param[in]  ip           Rank of the segment (starting from 1)
 ** \param[in]  family       Family to which the fracture belongs
 ** \param[in]  perm         Assigned permeability or TEST
 ** \param[in,out] nout      Initial/Final number of end points
 **
 *****************************************************************************/
static double* st_wellout_add(double *wellout,
                              int nbyout,
                              double x,
                              double y,
                              int ifrac,
                              int ip,
                              int family,
                              double perm,
                              int *nout)
{
  int nloc;

  nloc = *nout;
  wellout = (double*) mem_realloc((char* ) wellout,
                                  sizeof(double) * nbyout * (nloc + 1), 0);
  if (wellout == nullptr) return (wellout);

  WELLOUT(nloc,0) = x; /* First coordinate */
  WELLOUT(nloc,1) = y; /* Second coordinate */
  WELLOUT(nloc,2) = (double) ifrac + 1; /* Rank of the fracture */
  WELLOUT(nloc,3) = (double) ip + 1; /* Rank of segment in frature */
  WELLOUT(nloc,4) = (double) family; /* Rank of family */
  WELLOUT(nloc,5) = perm; /* Starting permeability */
  WELLOUT(nloc,6) = perm; /* Ending permeability */
  WELLOUT(nloc,7) = TEST; /* Range of permeability change */

  *nout = nloc + 1;
  return (wellout);
}

/****************************************************************************/
/*!
 **  Plunge a well line in a set of fractures
 **
 **  \return Array of the intersections (should be checked against NULL)
 **
 ** \param[in]  nval         Number of well information
 ** \param[in]  well         Array giving the well trajectory
 ** \param[in]  frac_list    Pointer to the Frac_List structure
 ** \param[in]  xmax         Maximum extension along horizontal axis
 ** \param[in]  permtab      Permabilities per family (starting from 0) Optional
 **
 ** \param[out] nint         Number of intersections
 ** \param[out] ncol         Number of attributes per intersection
 **
 ** \remark Output array must be freed by the calling function
 **
 *****************************************************************************/
double* fracture_to_well(int nval,
                                         double *well,
                                         Frac_List *frac_list,
                                         double xmax,
                                         double *permtab,
                                         int *nint,
                                         int *ncol)
{
  int nout, nbyout, nw_xy, family;
  double *wellout, x1, y1, x2, y2, x, y, perm;

  /* Initializations */

  nout = 0;
  nbyout = st_get_nbyout();
  wellout = nullptr;

  /* Preliminary checks */

  if (!isMultiple(nval, 2))
  {
    messerr("The number of values read from 'well' should be a multiple of 2");
    goto label_end;
  }
  nw_xy = nval / 2;
  if (nw_xy < 2)
  {
    messerr("Number of end points for the well line must not be less than 2");
    goto label_end;
  }

  /* Loop on the line segments */

  for (int iw = 0; iw < nw_xy - 1; iw++)
  {
    x1 = WELL(iw, 0);
    y1 = WELL(iw, 1);
    x2 = WELL(iw + 1, 0);
    y2 = WELL(iw + 1, 1);

    /* Store the starting point */

    wellout = st_wellout_add(wellout, nbyout, x1, y1, -1, -1, 0, 0., &nout);
    if (wellout == nullptr) goto label_end;

    /* Loop on the fractures */

    for (int ifrac = 0; ifrac < frac_list->nfracs; ifrac++)
    {
      Frac_Desc &desc = frac_list->frac_descs[ifrac];
      family = desc.family;

      /* Loop on the segments */

      for (int ip = 0; ip < desc.npoint - 1; ip++)
      {
        if (segment_intersect(x1, y1, x2, y2, XXF(ip), YYF(ip), XXF(ip + 1),
                              YYF(ip + 1), &x, &y)) continue;
        if (!INSIDE(x, y)) continue;

        /* Store the new intersection */

        perm = (permtab == nullptr) ? 0. :
                                      permtab[family];
        wellout = st_wellout_add(wellout, nbyout, x, y, ifrac, ip, family, perm,
                                 &nout);
        if (wellout == nullptr) goto label_end;
      }
    }

    /* Store the ending point */

    wellout = st_wellout_add(wellout, nbyout, x2, y2, -1, -1, 0, 0., &nout);
    if (wellout == nullptr) goto label_end;
  }

  label_end: (*nint) = nout;
  (*ncol) = nbyout;
  return (wellout);
}

/****************************************************************************/
/*!
 **  Add a sample to the trajectory
 **
 ** \param[in]  traj         Array for storing the trajectory
 ** \param[in]  x            First coordinate
 ** \param[in]  y            Second coordinate
 ** \param[in,out] ntraj     Number of samples in the trajectory
 **
 *****************************************************************************/
static void st_traj_add(double *traj, double x, double y, int *ntraj)
{
  int nloc;

  nloc = *ntraj;
  TRAJ(nloc,0) = x;
  TRAJ(nloc,1) = y;
  *ntraj = nloc + 1;
}

/****************************************************************************/
/*!
 **  Plunge a line trajectory and the modified permeabilities within an
 **  existing block
 **
 **  \return Error return code
 **
 ** \param[in,out] dbgrid    Db structure
 ** \param[in]  frac_list    Pointer to the Frac_List structure
 ** \param[in]  col_perm     Existing attribute for permeability (or ITEST)
 ** \param[in]  col_fluid    Existing attribute for fluid (or ITEST)
 ** \param[in]  flag_fluid   1 for performing the Fluid filling
 ** \param[in]  val_fluid    Value assigned to the fluid
 ** \param[in]  wellout      Pointer to the new wellout information
 ** \param[in]  nval         Number of values for wellout informations
 ** \param[in]  ndisc        Number of discretization steps
 ** \param[in]  verbose      Verbose flag
 **
 *****************************************************************************/
int fracture_well_to_block(Db *dbgrid,
                                           Frac_List *frac_list,
                                           int col_perm,
                                           int col_fluid,
                                           int flag_fluid,
                                           double val_fluid,
                                           double *wellout,
                                           int nval,
                                           int ndisc,
                                           int verbose)
{
  int *indg;
  int iptr_perm, iptr_fluid, nbyout, ifrac, ip, error, idim, npoint, ntraj,
      nw_xy;
  double *traj, dmin, delta, xx, yy, x1, y1, x2, y2, perm1, perm2, range;

  /* Initializations */

  error = 1;
  nbyout = st_get_nbyout();
  indg = nullptr;
  traj = nullptr;
  iptr_perm = 0;
  iptr_fluid = 0;

  /* Preliminary checks */

  if (!isMultiple(nval, nbyout))
  {
    messerr("The number of values read (%d) in well should be a multiple of %d",
            nval, nbyout);
    goto label_end;
  }
  nw_xy = nval / nbyout;
  if (!is_grid(dbgrid))
  {
    messerr("The Db file must be organized as a regular Grid");
    goto label_end;
  }
  if (dbgrid->getNDim() != 2)
  {
    messerr("This application is limited to 2-D grid");
    goto label_end;
  }
  dmin = 1.e30;
  for (idim = 0; idim < dbgrid->getNDim(); idim++)
  {
    if (dbgrid->getDX(idim) < dmin) dmin = dbgrid->getDX(idim);
  }
  delta = dmin / (double) ndisc;

  /* Allocate the new variable */

  indg = db_indg_alloc(dbgrid);
  if (indg == nullptr) goto label_end;

  iptr_perm = dbgrid->addFields(1, 0);
  if (!IFFFF(col_perm)) db_attribute_copy(dbgrid, col_perm, iptr_perm);

  if (flag_fluid)
  {
    iptr_fluid = dbgrid->addFields(1, 0.);
    if (!IFFFF(col_fluid)) db_attribute_copy(dbgrid, col_fluid, iptr_fluid);
  }

  /* Verbose option */

  if (verbose)
    print_matrix("Well information", 0, 0, nbyout, nw_xy, NULL, wellout);

  /* Paint the fluid (optional) */

  if (flag_fluid)
  {
    for (int iw = 0; iw < nw_xy - 1; iw++)
    {
      x1 = WELLOUT(iw, 0);
      y1 = WELLOUT(iw, 1);
      x2 = WELLOUT(iw + 1, 0);
      y2 = WELLOUT(iw + 1, 1);
      st_plunge_segment(dbgrid, iptr_fluid, indg, delta, val_fluid, x1, y1, x2,
                        y2);
    }
  }

  /* Loop on the end points of the wellout */

  for (int iw = 0; iw < nw_xy; iw++)
  {
    xx = WELLOUT(iw, 0);
    yy = WELLOUT(iw, 1);
    ifrac = (int) WELLOUT(iw, 2) - 1;
    ip = (int) WELLOUT(iw, 3) - 1;
    perm1 = WELLOUT(iw, 5);
    perm2 = WELLOUT(iw, 6);
    range = WELLOUT(iw, 7);
    if (ifrac < 0) continue;

    if (ifrac >= frac_list->nfracs)
    {
      messerr("The well information (line %d/%d) is invalid:", iw + 1, nw_xy);
      messerr("- Fracture number (%d) should lie within [1,%d]", ifrac + 1,
              frac_list->nfracs);
      continue;
    }

    Frac_Desc &desc = frac_list->frac_descs[ifrac];
    if (ip < 0 || ip >= desc.npoint)
    {
      messerr("The well information (line %d/%d) is invalid:", iw + 1, nw_xy);
      messerr("- For fracture (%d, Segment number (%d) is not within [1,%d]",
              ifrac + 1, ip + 1, desc.npoint);
      continue;
    }
    npoint = desc.npoint;

    /* Paint the blocks with the stationary permeability */

    for (int i = 0; i < npoint - 1; i++)
      st_plunge_segment(dbgrid, iptr_perm, indg, delta, perm2, XXF(i), YYF(i),
                        XXF(i + 1), YYF(i + 1));

    /* Allocate the trajectory */

    traj = (double*) mem_alloc(sizeof(double) * (npoint + 1) * 2, 0);
    if (traj == nullptr) goto label_end;

    /* Trajectory upwards */

    ntraj = 0;
    st_traj_add(traj, xx, yy, &ntraj);
    for (int i = ip; i >= 0; i--)
      st_traj_add(traj, XXF(ip), YYF(ip), &ntraj);
    st_plunge_segment_gradual(dbgrid, iptr_perm, indg, delta, traj, ntraj,
                              perm1, perm2, range);

    /* Trajectory downwards */

    ntraj = 0;
    st_traj_add(traj, xx, yy, &ntraj);
    for (int i = ip; i < npoint; i++)
      st_traj_add(traj, XXF(ip), YYF(ip), &ntraj);
    st_plunge_segment_gradual(dbgrid, iptr_perm, indg, delta, traj, ntraj,
                              perm1, perm2, range);

    /* Core deallocation */

    traj = (double*) mem_free((char* ) traj);
  }

  /* Set the error return code */

  error = 0;

  label_end: traj = (double*) mem_free((char* ) traj);
  indg = db_indg_free(indg);
  return (error);
}

