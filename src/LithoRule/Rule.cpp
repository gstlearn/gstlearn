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
#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Basic/AException.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/Node.hpp"
#include "geoslib_f.h"
#include "geoslib_enum.h"

#define THRESH_IDLE 0
#define THRESH_Y1   1
#define THRESH_Y2   2

#define NODES(inode,i)       (nodes[6 * (inode) + (i)])
#define FROM_TYPE(inode)     (nodes[6 * (inode) + 0])
#define FROM_RANK(inode)     (nodes[6 * (inode) + 1])
#define FROM_VERS(inode)     (nodes[6 * (inode) + 2])
#define NODE_TYPE(inode)     (nodes[6 * (inode) + 3])
#define NODE_RANK(inode)     (nodes[6 * (inode) + 4])
#define FACIES(inode)        (nodes[6 * (inode) + 5])

static const VectorString symbol = {"F","S","T"};
static int GAUSS_MODE = 1;

/****************************************************************************/
/*!
**  Set the Rule Mode
**
** \param[in]  rule_mode   1 for Gaussian; 0 for absence of conversion
**
** \remarks The absence of conversion is used in order to evaluate the
** \remarks real thresholds along the rule for representing the rule
** \remarks by proportions rather than in gaussian scale
**
*****************************************************************************/
void set_rule_mode(int rule_mode)
{
  GAUSS_MODE = rule_mode;
}

/****************************************************************************/
/*!
**  Get the current Rule Mode
**
** \return  Returns the current mode (1 for Gaussian; 0 for Raw)
**
*****************************************************************************/
int get_rule_mode(void)
{
  return(GAUSS_MODE);
}

/****************************************************************************/
/*!
**  Get the lower or upper bound
**
** \return  Lower of Upper Bound
**
** \param[in]  mode    <0 for lower bound; >0 for upper bound
**
*****************************************************************************/
double get_rule_extreme(int mode)
{
  if (mode < 0)
  {
    if (GAUSS_MODE)
      return(THRESH_INF);
    else
      return(0);
  }
  else
  {
    if (GAUSS_MODE)
      return(THRESH_SUP);
    else
      return(1.);
  }
}

Rule::Rule(int mode_rule,
           double rho,
           double slope,
           double sh_down,
           double sh_dsup,
           const VectorDouble& shift)
    : AStringable(),
      ASerializable(),
      _modeRule(mode_rule),
      _flagProp(0),
      _rho(rho),
      _shDsup(sh_dsup),
      _shDown(sh_down),
      _slope(slope),
      _dMax(TEST),
      _tgte(TEST),
      _incr(TEST),
      _shift(shift),
      _xyz(),
      _ind1(),
      _ind2(),
      _mainNode(nullptr)
{
}

Rule::Rule(const VectorInt& n_type, const VectorInt& n_facs, double rho)
    : AStringable(),
      ASerializable(),
      _modeRule(RULE_STD),
      _flagProp(0),
      _rho(rho),
      _shDsup(TEST),
      _shDown(TEST),
      _slope(TEST),
      _dMax(TEST),
      _tgte(TEST),
      _incr(TEST),
      _shift(),
      _xyz(),
      _ind1(),
      _ind2(),
      _mainNode(nullptr)
{
  int ipos = 0;
  int n_fac = 0;
  int n_y1 = 0;
  int n_y2 = 0;
  _mainNode = new Node("main", n_type, n_facs, &ipos, &n_fac, &n_y1, &n_y2);
}

Rule::Rule(const VectorString& nodnames, double rho)
    : AStringable(),
      ASerializable(),
      _modeRule(RULE_STD),
      _flagProp(0),
      _rho(rho),
      _shDsup(TEST),
      _shDown(TEST),
      _slope(TEST),
      _dMax(TEST),
      _tgte(TEST),
      _incr(TEST),
      _shift(),
      _xyz(),
      _ind1(),
      _ind2(),
      _mainNode(nullptr)
{
  VectorInt n_type;
  VectorInt n_facs;
  _nodNamesToIds(nodnames, n_type, n_facs);
  int ipos = 0;
  int n_fac = 0;
  int n_y1 = 0;
  int n_y2 = 0;
  _mainNode = new Node("main", n_type, n_facs, &ipos, &n_fac, &n_y1, &n_y2);
}

/**
 * Defining the Lithotype Rule in the case of Shadow
 * @param slope   Slope definition
 * @param sh_dsup Maximum threshold
 * @param sh_down Minimum threshold
 * @param shift   Shift orientation
 */
Rule::Rule(double slope,
           double sh_dsup,
           double sh_down,
           const VectorDouble& shift)
    : AStringable(),
      ASerializable(),
      _modeRule(RULE_SHADOW),
      _flagProp(0),
      _rho(0.),
      _shDsup(sh_dsup),
      _shDown(sh_down),
      _slope(slope),
      _dMax(TEST),
      _tgte(TEST),
      _incr(TEST),
      _shift(shift),
      _xyz(),
      _ind1(),
      _ind2(),
      _mainNode(nullptr)
{
//  rule->main         = st_node_alloc("S1",THRESH_Y1,0);
//  rule->main->r1     = st_node_alloc("T1",THRESH_Y2,0);
//  rule->main->r1->r1 = st_node_alloc("F3",THRESH_IDLE,SHADOW_SHADOW);
//  rule->main->r1->r2 = st_node_alloc("F2",THRESH_IDLE,SHADOW_WATER);
//  rule->main->r2     = st_node_alloc("F1",THRESH_IDLE,SHADOW_ISLAND);
}

/**
 * Definition of the Lithotype Rule in the case of Shift
 * @param shift Vector defining the shift
 */
Rule::Rule(const VectorDouble& shift)
    : AStringable(),
      ASerializable(),
      _modeRule(RULE_SHIFT),
      _flagProp(0),
      _rho(0.),
      _shDsup(0.),
      _shDown(0.),
      _slope(0.),
      _dMax(TEST),
      _tgte(TEST),
      _incr(TEST),
      _shift(shift),
      _xyz(),
      _ind1(),
      _ind2(),
      _mainNode(nullptr)
{
}

Rule::Rule(const String& neutralFileName, bool verbose)
    : AStringable(),
      ASerializable(),
      _modeRule(RULE_STD),
      _flagProp(0),
      _rho(0.),
      _shDsup(0.),
      _shDown(0.),
      _slope(0.),
      _dMax(TEST),
      _tgte(TEST),
      _incr(TEST),
      _shift(),
      _xyz(),
      _ind1(),
      _ind2(),
      _mainNode(nullptr)
{
  if (deSerialize(neutralFileName, verbose))
    my_throw("Problem reading the Neutral File");
}

Rule::Rule(const Rule& m)
    : _modeRule(m._modeRule),
      _flagProp(m._flagProp),
      _rho(m._rho),
      _shDsup(m._shDsup),
      _shDown(m._shDown),
      _slope(m._slope),
      _dMax(m._dMax),
      _tgte(m._tgte),
      _incr(m._incr),
      _shift(m._shift),
      _mainNode(new Node(*m._mainNode))
{
}

Rule& Rule::operator=(const Rule& m)
{
  if (this != &m)
  {
    _modeRule = m._modeRule;
    _flagProp = m._flagProp;
    _rho = m._rho;
    _shDsup = m._shDsup;
    _shDown = m._shDown;
    _slope = m._slope;
    _dMax = m._dMax;
    _tgte = m._tgte;
    _incr = m._incr;
    _shift = m._shift;
    _mainNode = new Node(*m._mainNode);
  }
  return *this;
}

Rule::~Rule()
{
  if (_mainNode != nullptr)
    delete _mainNode;
}

void Rule::init(int nfacies)
{
  if (nfacies <= 1)
    my_throw("This class requires at least TWO facies");
  VectorString nodnames = _buildNodNames(nfacies);
  VectorInt n_type;
  VectorInt n_facs;
  _nodNamesToIds(nodnames, n_type, n_facs);
  int ipos = 0;
  int n_fac = 0;
  int n_y1 = 0;
  int n_y2 = 0;
  _mainNode = new Node("main", n_type, n_facs, &ipos, &n_fac, &n_y1, &n_y2);
}

/**
 * Initialization of the nodes of the Rule (from the ASCII file).
 * @param nodes Node description (6 * number of nodes)
 * 0 : Type of the parent
 * 1 : Rank of the parent
 * 2 : Orientation of the parent
 * 3 : Information type: 0 (idle) - 1 (Threshold along Y1) - 2 (Threshold along Y2)
 * 4 : Rank of the node (starting from 1)
 * 5 : Rank of the facies
 * @return Newly created Rule structure
 */
int Rule::init(const VectorInt& nodes)
{
  int nb_node = static_cast<int> (nodes.size()) / 6;
  std::vector<Node *> n1tab(nb_node, nullptr);
  std::vector<Node *> n2tab(nb_node, nullptr);

  // Loop on the nodes

  for (int inode = 0; inode < nb_node; inode++)
  {

    /* Check the validity of the current node information */

    if (NODE_TYPE(inode) != THRESH_IDLE &&
        NODE_TYPE(inode) != THRESH_Y1   &&
        NODE_TYPE(inode) != THRESH_Y2)
    {
      messerr("Error in the type of the node #%d (%d)",
              inode + 1,NODE_TYPE(inode));
      return 1;
    }
    if (NODE_RANK(inode) < 1 || NODE_RANK(inode) > nb_node)
    {
      messerr("Error: the rank of the node #%d (%d) must lie within [1;%d]",
              inode + 1, NODE_RANK(inode), nb_node);
      return 1;
    }
    if ((NODE_TYPE(inode) == THRESH_Y1 && n1tab[NODE_RANK(inode) - 1] != (Node *) NULL) ||
        (NODE_TYPE(inode) == THRESH_Y2 && n2tab[NODE_RANK(inode) - 1] != (Node *) NULL))
    {
      messerr("Error: Node #%d (%s%d) has already been created", inode + 1,
              symbol[NODE_TYPE(inode)].c_str(), NODE_RANK(inode));
      return 1;
    }

    /* Check the validity of the parent node */

    if (inode > 0)
    {
      int found = -1;
      for (int jnode = 0; jnode < inode && found < 0; jnode++)
      {
        if (FROM_TYPE(inode) == NODE_TYPE(jnode) &&
            FROM_RANK(inode) == NODE_RANK(jnode)) found = jnode;
      }
      if (found < 0)
      {
        messerr("Node #%d (%s%d) refers to unknown parent (%s%d)", inode + 1,
                symbol[NODE_TYPE(inode)].c_str(), NODE_RANK(inode),
                symbol[FROM_TYPE(inode)].c_str(), FROM_RANK(inode));
        return 1;
      }
    }

    /* Create the nodes */

    int facies = (NODE_TYPE(inode) == THRESH_IDLE) ? FACIES(inode) : 0;

    std::stringstream name;
    if (NODE_TYPE(inode) == THRESH_IDLE)
      name << symbol[NODE_TYPE(inode)] << FACIES(inode);
    else
      name << symbol[NODE_TYPE(inode)];

    // Allocate the new node

    Node* node_loc = new Node(name.str(), NODE_TYPE(inode), facies);
    if (inode == 0) _mainNode = node_loc;

    /* Link to the previous pointer */

    if (FROM_TYPE(inode) == THRESH_Y1)
    {
      if (FROM_VERS(inode) == 1)
        n1tab[FROM_RANK(inode) - 1]->setR1(node_loc);
      else
        n1tab[FROM_RANK(inode) - 1]->setR2(node_loc);
    }
    if (FROM_TYPE(inode) == THRESH_Y2)
    {
      if (FROM_VERS(inode) == 1)
        n2tab[FROM_RANK(inode) - 1]->setR1(node_loc);
      else
        n2tab[FROM_RANK(inode) - 1]->setR2(node_loc);
    }

    /* Store the pointer */

    switch (NODE_TYPE(inode))
    {
      case THRESH_IDLE:
        break;

      case THRESH_Y1:
        n1tab[NODE_RANK(inode) - 1] = node_loc;
        break;

      case THRESH_Y2:
        n2tab[NODE_RANK(inode) - 1] = node_loc;
        break;
    }
  }
  return 0;
}

std::string Rule::toString(int level) const
{
  std::stringstream sstr;
  sstr << _display(false, false);
  return sstr.str();
}

void Rule::display(bool flagProp, bool flagThresh) const
{
  messageFlush(_display(flagProp, flagThresh));
}

String Rule::_display(bool flagProp, bool flagThresh) const
{
  std::stringstream sstr;
  int node_tot,nfac_tot,nmax_tot,ny1_tot,ny2_tot;
  double prop_tot;

  sstr << toTitle(0,"Lithotype Rule");

  if (statistics(0, &node_tot, &nfac_tot, &nmax_tot, &ny1_tot, &ny2_tot,
                 &prop_tot)) return sstr.str();
  if (prop_tot <= 0.) flagProp = flagThresh = false;

  sstr << "- Number of nodes               = " << node_tot << std::endl;
  sstr << "- Number of facies              = " << nfac_tot << std::endl;
  sstr << "- Number of thresholds along G1 = " << ny1_tot  << std::endl;
  sstr << "- Number of thresholds along G2 = " << ny2_tot  << std::endl;
  sstr << std::endl;

  switch (_modeRule)
  {
    case RULE_STD:
      if (ABS(_rho) > 0.)
        sstr << "Correlation between the two GRFs = " << _rho << std::endl;
      if (_rho == 1.)
        sstr << "(With the current option, only the first GRF is used)" << std::endl;
      sstr << _mainNode->nodePrint(flagProp, flagThresh);
      break;

    case RULE_SHIFT:
      sstr << toVector("Translation Vector",_shift) << std::endl;
      sstr << "(With the current option, only the first GRF is used)" << std::endl;
      sstr << _mainNode->nodePrintShadow(flagProp, flagThresh);
      break;

    case RULE_SHADOW:
      sstr << toVector("Normalized Translation Vector = ",_shift) << std::endl;
      sstr << "Slope for shadow                  = " << _slope << "(degrees)" << std::endl;
      sstr << "Upwards shift for the threshold   = " << _shDsup << std::endl;
      sstr << "Downwards shift for the threshold = " << _shDown << std::endl;
      sstr << std::endl;
      sstr << "Note for non-stationary case:" << std::endl;
      sstr << "- P1 gives the proportion of Island" << std::endl;
      sstr << "- P2 gives the value of Upwards shift" << std::endl;
      sstr << "- P3 gives the value of Downwards shift" << std::endl;
      sstr << "(With the current option, only the first GRF is used)" << std::endl;
      sstr << _mainNode->nodePrintShadow(flagProp, flagThresh);
      break;
  }
  return sstr.str();
}

int Rule::getFaciesNumber() const
{
  int node_tot, nfac_tot, nmax_tot, ny1_tot, ny2_tot;
  double prop_tot;
  if (statistics(0, &node_tot, &nfac_tot, &nmax_tot, &ny1_tot, &ny2_tot,
                 &prop_tot)) return 0;
  return nfac_tot;
}

int Rule::getGRFNumber() const
{
  int ny2 = getY2Number();
  return (ny2 > 0 ? 2 : 1);
}

int Rule::getY1Number() const
{
  int node_tot, nfac, nmax, ny1, ny2;
  double prop;
  if (statistics(0, &node_tot, &nfac, &nmax, &ny1, &ny2, &prop)) return 0;
  return ny1;
}

int Rule::getY2Number() const
{
  int node_tot, nfac, nmax, ny1, ny2;
  double prop;
  if (statistics(0, &node_tot, &nfac, &nmax, &ny1, &ny2, &prop)) return 0;
  if (getModeRule() == RULE_SHADOW || getModeRule() == RULE_SHIFT) ny2 = 0;
  if (_rho == 1.) ny2 = 0;
  return ny2;
}

bool Rule::isYUsed(int igrf) const
{
  if (igrf == 0)
    return getY1Number() > 0;
  else
    return getY2Number() > 0;
}

VectorInt Rule::whichGRFUsed() const
{
  VectorInt flag(2);
  for (int igrf = 0; igrf < 2; igrf++)
    flag[igrf] = isYUsed(igrf);
  return flag;
}

/****************************************************************************/
/*!
**  Calculates the statistics from the Lithotype Rule
**
** \return  Error return code
**
** \param[in]  verbose  1 for a verbose output; 0 otherwise
**
** \param[out]  node_tot Number of nodes
** \param[out]  nfac_tot Number of facies
** \param[out]  nmax_tot Number of (different) facies
** \param[out]  ny1_tot  Number of thresholds for Y1
** \param[out]  ny2_tot  Number of thresholds for Y2
** \param[out]  prop_tot Total proportion
**
** \remark The tolerance on the sum of the proportions can be defined
** \remark using set_keypair("TolSumProportions",newtol)
**
*****************************************************************************/
int Rule::statistics(int  verbose,
                     int *node_tot,
                     int *nfac_tot,
                     int *nmax_tot,
                     int *ny1_tot,
                     int *ny2_tot,
                     double *prop_tot) const
{
  int nfac,ifac,ntot;

  /* Establish the statistics on the Lithotype Rule */

  _mainNode->getStatistics(node_tot,nfac_tot,ny1_tot,ny2_tot,prop_tot);

  /* Check that the facies are defined */

  nfac = (*nfac_tot);
  VectorInt facies = VectorInt(nfac);
  for (ifac=0; ifac<nfac; ifac++) facies[ifac] = 0;
  if (_mainNode->isValid(facies)) return 1;

  /* Check that the first consecutive facies are defined */

  ntot = 0;
  for (ifac=0; ifac<nfac; ifac++)
    if (facies[ifac] > 0) ntot = ifac + 1;
  for (ifac=0; ifac<nfac; ifac++)
  {
    if (facies[ifac] <= 0)
    {
      messerr("The facies (%d) is not defined",ifac+1);
      return(1);
    }
  }
  (*nmax_tot) = ntot;

  /* If proportions are defined, check the sum of the proportions */

  if (getFlagProp())
  {
    if (ABS((*prop_tot) - 1.) > EPSILON2)
    {
      messerr("Error: Cumulated proportions not equal to 1 (%lf)",(*prop_tot));
      messerr("Tolerance                          = %lf",EPSILON2);
      messerr("Number of nodes                    = %d",(*node_tot));
      messerr("Number of facies                   = %d",(*nfac_tot));
      messerr("Number of different facies numbers = %d",(*nmax_tot));
      messerr("Number of thresholds along Y1      = %d",(*ny1_tot));
      messerr("Number of thresholds along Y2      = %d",(*ny2_tot));
    }
    else
    {
      _mainNode->scaleProp(*prop_tot);
      *prop_tot = 1;
    }
  }

  /* Optional printout */

  if (verbose)
  {
    mestitle(1,"Lithotype Rule");
    message("Number of nodes      = %d\n" ,(*node_tot));
    message("Number of facies     = %d\n" ,(*nfac_tot));
    message("Maximum facies rank  = %d\n" ,(*nmax_tot));
    message("Cumulated proportion = %lf\n",(*prop_tot));
  }

  return(0);
}

/****************************************************************************/
/*!
**  Define the particularities of the PGS model
**
** \return  Error return code
**
** \param[in]  db              Db structure
** \param[in]  dbprop          Db structure used for proportions
** \param[in]  model           Model structure (only used for shift option)
** \param[in]  flag_grid_check 1 if grid is compulsory; 0 otherwise
**                             (only for SHIFT)
** \param[in]  flag_stat       1 for stationary; 0 otherwise
**
*****************************************************************************/
int Rule::particularities(Db *db,
                          const Db *dbprop,
                          Model *model,
                          int flag_grid_check,
                          int flag_stat)
{
  int ndim = (model != (Model *) NULL) ? model->getDimensionNumber() : 0;
  VectorDouble wxyz(ndim);
  double rhoval;

  /* Dispatch */

  _xyz.resize(ndim);
  switch(getModeRule())
  {
    case RULE_STD:
      return(0);

    case RULE_SHIFT:
      double hval = 0.;
      for (int idim=0; idim<ndim; idim++)
      {
        _xyz[idim] = _shift[idim];
        hval += _xyz[idim] * _xyz[idim];
      }
      hval = sqrt(hval);

      /* Calculate the covariance between the two GRF */

      for (int idim = 0; idim < ndim; idim++) wxyz[idim] = _xyz[idim];
      if (model->getVariableNumber() == 1)
        model_evaluate(model, 0, 0, -1, 0, 1, 0, 0, 0, MEMBER_LHS, 1, wxyz,
                       &hval, &rhoval);
      else
        model_evaluate(model, 0, 1, -1, 0, 1, 0, 0, 0, MEMBER_LHS, 1, wxyz,
                       &hval, &rhoval);
      setRho(rhoval);

      /* Translate the shift into grid increments */

      if (_st_shift_on_grid(db,ndim,flag_grid_check)) return(1);
      break;
  }
  return(0);
}

/****************************************************************************/
/*!
**  Define the particularities of the PGS model (for Shadow)
**
** \return  Error return code
**
** \param[in]  db              Db structure
** \param[in]  dbprop          Db structure used for proportions
** \param[in]  model           Model structure (only used for shift option)
** \param[in]  flag_grid_check 1 if grid is compulsory; 0 otherwise
**                             (only for SHIFT)
** \param[in]  flag_stat       1 for stationary; 0 otherwise
**
*****************************************************************************/
int Rule::particularities_shadow(Db *db,
                                 const Db *dbprop,
                                 Model *model,
                                 int flag_grid_check,
                                 int flag_stat)
{
  double sh_dsup_max,sh_down_max;
  int ndim = (model != (Model *) NULL) ? model->getDimensionNumber() : 0;

  double norme = 0.;
  for (int idim=0; idim<ndim; idim++)
    norme += _shift[idim] * _shift[idim];
  norme = sqrt(norme);
  if (norme <= 0) return(1);
  for (int idim=0; idim<ndim; idim++) _shift[idim] /= norme;

  _incr = 1.e30;
  for (int idim=0; idim<ndim; idim++)
    if (_shift[idim] != 0) _incr = MIN(_incr, db->getDX(idim));

  /* Calculate the maximum distance */

  _tgte = tan(ut_deg2rad(_slope));
  _st_shadow_max(dbprop,flag_stat,&sh_dsup_max,&sh_down_max);
  _dMax = (_tgte > 0) ? (sh_dsup_max + sh_down_max) / _tgte : 0.;

  return(0);
}

void Rule::_st_shadow_max(const Db *dbprop,
                          int flag_stat,
                          double *sh_dsup_max,
                          double *sh_down_max)
{
  int iech;
  double val2,val3;

  if (flag_stat || dbprop == (Db *) NULL)
  {
    /* Stationary case */

    *sh_dsup_max = getShDsup();
    *sh_down_max = getShDown();
  }
  else
  {
    *sh_dsup_max = *sh_down_max = 0.;
    for (iech=0; iech<dbprop->getSampleNumber(); iech++)
    {
      val2 = dbprop->getProportion(iech,1);
      if (val2 > (*sh_dsup_max)) (*sh_dsup_max) = val2;
      val3 = dbprop->getProportion(iech,2);
      if (val3 > (*sh_down_max)) (*sh_down_max) = val3;
    }
  }
  return;
}

int Rule::_st_shift_on_grid(Db     *db,
                            int     ndim,
                            int     flag_grid_check)
{
  VectorDouble xyz(ndim);
  _ind1.resize(ndim);

  if (db == (Db *) NULL || ! is_grid(db))
  {
    if (! flag_grid_check) return(0);
    messerr("The shift Rule requires a Grid Db");
    return(1);
  }

  for (int idim=0; idim<ndim; idim++)
    _xyz[idim] = _shift[idim] + db->getX0(idim);

  (void) point_to_grid(db,_xyz.data(),-1,_ind1.data());

  /* Check that the translation is significant */

  int ntot = 0;
  for (int idim=0; idim<ndim; idim++)
    ntot += ABS(_ind1[idim]);
  if (ntot <= 0)
  {
    messerr("The shift of the Lithotype Rule cannot be rendered");
    messerr("using the Output Grid characteristics");
    return(1);
  }
  return(0);
}

/****************************************************************************/
/*!
**  Evaluation of the value from a grid by inverse square distance
**  interpolation from the 4 surrounding grid nodes
**
** \return  Interpolated value or FFFF if out of grid
**
** \param[in]  dbgrid Db structure
** \param[in]  isimu  Rank of the simulation
** \param[in]  icase  Rank of the Simulation storage
** \param[in]  nbsimu Number of simulations
**
** \param[out]  xyz0  Working array
**
*****************************************************************************/
double Rule::st_grid_eval(Db *dbgrid,
                          int isimu,
                          int icase,
                          int nbsimu,
                          VectorDouble& xyz0)
{
  double top = 0.;
  double bot = 0.;
  int ndim = dbgrid->getNDim();

  /* First point */
  int iech = db_index_grid_to_sample(dbgrid,_ind2.data());
  double z = dbgrid->getSimvar(LOC_SIMU,iech,isimu,0,icase,nbsimu,1);
  if (! FFFF(z))
  {
    double d2 = 0.;
    grid_to_point(dbgrid,_ind2.data(),NULL,xyz0.data());
    for (int idim=0; idim<ndim; idim++)
    {
      double delta = _xyz[idim] - xyz0[idim];
      d2 += delta * delta;
    }
    if (d2 <= 0.) return(z);
    top += z / d2;
    bot += 1./ d2;
  }

  /* Second point */
  _ind2[0] += 1;
  iech = db_index_grid_to_sample(dbgrid,_ind2.data());
  z = dbgrid->getSimvar(LOC_SIMU,iech,isimu,0,icase,nbsimu,1);
  if (! FFFF(z))
  {
    double d2 = 0.;
    grid_to_point(dbgrid,_ind2.data(),NULL,xyz0.data());
    for (int idim=0; idim<ndim; idim++)
    {
      double delta = _xyz[idim] - xyz0[idim];
      d2 += delta * delta;
    }
    if (d2 <= 0.) return(z);
    top += z / d2;
    bot += 1./ d2;
  }

  /* Third point */
  _ind2[1] += 1;
  iech = db_index_grid_to_sample(dbgrid,_ind2.data());
  z = dbgrid->getSimvar(LOC_SIMU,iech,isimu,0,icase,nbsimu,1);
  if (! FFFF(z))
  {
    double d2 = 0.;
    grid_to_point(dbgrid,_ind2.data(),NULL,xyz0.data());
    for (int idim=0; idim<ndim; idim++)
    {
      double delta = _xyz[idim] - xyz0[idim];
      d2 += delta * delta;
    }
    if (d2 <= 0.) return(z);
    top += z / d2;
    bot += 1./ d2;
  }

  /* Fourth point */
  _ind2[0] -= 1;
  iech = db_index_grid_to_sample(dbgrid,_ind2.data());
  z = dbgrid->getSimvar(LOC_SIMU,iech,isimu,0,icase,nbsimu,1);
  if (! FFFF(z))
  {
    double d2 = 0.;
    grid_to_point(dbgrid,_ind2.data(),NULL,xyz0.data());
    for (int idim=0; idim<ndim; idim++)
    {
      double delta = _xyz[idim] - xyz0[idim];
      d2 += delta * delta;
    }
    if (d2 <= 0.) return(z);
    top += z / d2;
    bot += 1./ d2;
  }

  /* Final interpolation */
  _ind2[1] -= 1;
  z = (bot != 0) ? top / bot : TEST;
  return(z);
}

void Rule::updateShift()
{
  Node* node;
  node  = _mainNode->getR2();
  double seuil = node->getT1min();
  node  = _mainNode->getR1()->getR1();
  node->setT2max(seuil);
  node  = _mainNode->getR1()->getR2();
  node->setT2min(seuil);
}

void Rule::_nodNamesToIds(const VectorString& nodes,
                          VectorInt& n_type,
                          VectorInt& n_facs)
{
  int nb_node = static_cast<int> (nodes.size());
  n_type.resize(nb_node,0);
  n_facs.resize(nb_node,0);

  for (int i = 0; i < nb_node; i++)
    decodeInList(symbol,nodes[i],&n_type[i],&n_facs[i]);
}

/**
 * Returns the proportion of a given facies
 * @param facies Facies rank (starting from 1)
 * @return Proportion of the given Facies
 */
double Rule::getProportion(int facies)
{
  double prop;
  if (_mainNode->getProportion(facies,&prop))
    return prop;
  else
    return TEST;
}

/**
 * Return the vector of bounds for a given facies
 * @param facies Rank of the target facies (starting from 1)
 * @return The vector of bounds organized as [t1min, t1max, t2min, t2max]
 */
VectorDouble Rule::getThresh(int facies) const
{
  int fac_ret;
  int rank = 0;
  double t1min, t1max, t2min, t2max;

  if (!_mainNode->getThresh(1, facies, &rank, &fac_ret, &t1min, &t1max, &t2min,
                            &t2max)) return VectorDouble();
  VectorDouble bounds(4);
  bounds[0] = t1min;
  bounds[1] = t1max;
  bounds[2] = t2min;
  bounds[3] = t2max;
  return bounds;
}

VectorDouble Rule::getThreshFromRectangle(int rect, int *facies)
{
  VectorDouble bounds;

  int rank = 0;
  double t1min, t1max, t2min, t2max;

  if (!_mainNode->getThresh(2, rect, &rank, facies, &t1min, &t1max, &t2min,
                            &t2max)) return bounds;
  bounds.resize(4);
  bounds[0] = t1min;
  bounds[1] = t1max;
  bounds[2] = t2min;
  bounds[3] = t2max;
  return bounds;
}

/****************************************************************************/
/*!
**  Convert the two underlying GRFs into facies
**
** \return  The facies rank or 0 (facies not found)
**
** \param[in]  y1     Value of the first underlying GRF
** \param[in]  y2     Value of the second underlying GRF
**
** \remark  If one of the two GRF is undefined, the resulting facies is 0
**
*****************************************************************************/
int Rule::getFaciesFromGaussian(double y1, double y2) const
{
  double facies;

  if (FFFF(y1) || FFFF(y2)) return(0);
  if (! _mainNode->gaussianToFacies(y1,y2,&facies)) return(0);
  return((int) facies);
}

/**
 * Define constant proportions
 * @param proportions The vector of constant proportions.
 * It should be dimensioned to the number of facies.
 * If absent, all proportions are considered equal.
 * @return
 */
int Rule::setProportions(const VectorDouble& proportions)
{
  int    node_tot,nfac_tot,nmax_tot,ny1_tot,ny2_tot;
  double prop_tot;

  // Set the proportions when the input argument is left empty

  VectorDouble props = proportions;
  if (props.empty())
  {
    int nfacies = getFaciesNumber();
    props = VectorDouble(nfacies, 1. / (double) nfacies);
  }

  /* Set the proportions */

  if (_mainNode->proportionDefine(props)) return 1;
  _flagProp = 1;

  /* Calculate the cumulative proportions */

  statistics(0, &node_tot, &nfac_tot, &nmax_tot, &ny1_tot, &ny2_tot, &prop_tot);

  /* Convert from proportions to thresholds */

  _mainNode->proportionToThresh(_rho,
                                get_rule_extreme(-1), get_rule_extreme(+1),
                                get_rule_extreme(-1), get_rule_extreme(+1));

  /* Debug printout (optional) */

  if (debug_query("props")) display(true,true);

  return(0);
}

int Rule::deSerialize(const String& filename, bool verbose)
{
  int nb_node;

  if (_fileOpen(filename, "Rule", "r", verbose)) return 1;

  /* Create the Rule structure */

  _shift.resize(3);
  if (_recordRead("Rule definition", "%d", &_modeRule)) return 1;
  if (_recordRead("Correlation Coefficient of GRFs", "%lf", &_rho)) return 1;
  if (_recordRead("Slope for Shadow Rule", "%lf", &_slope)) return 1;
  if (_recordRead("Lower Threshold for Shadow Rule", "%lf", &_shDown)) return 1;
  if (_recordRead("Upper Threshold for Shadow Rule", "%lf", &_shDsup)) return 1;
  if (_recordRead("Shift along first direction", "%lf", &_shift[0]))  return 1;
  if (_recordRead("Shift along second direction", "%lf", &_shift[1])) return 1;
  if (_recordRead("Shift along third direction", "%lf", &_shift[2]))  return 1;

  /* Read the number of nodes */

  if (_recordRead("Number of Rule Nodes", "%d", &nb_node)) return 1;
  VectorInt nodes(6 * nb_node);

  /* Loop on the nodes for reading: */
  /* - from_type: Type of the parent */
  /* - from_rank: Rank of the parent */
  /* - from_vers: Orientation of the parent */
  /* - node_type: 0 (idle) - 1 (Thresh along Y1) - 2 (Thresh along Y2) */
  /* - node_rank: Rank of the node (starting from 1) */
  /* - facies   : Rank of the facies */
  int lec = 0;
  for (int inode =  0; inode < nb_node; inode++)
    for (int i = 0; i < 6; i++)
      if (_recordRead("Rule Node Definition", "%d", &nodes[lec++])) return 1;
  init(nodes);

  _fileClose(verbose);

  return 0;
}

int Rule::serialize(const String& filename, bool verbose) const
{
  int nb_node, nfacies, nmax_tot, ny1_tot, ny2_tot, rank;
  double prop_tot;

  if (_fileOpen(filename, "Rule", "w", verbose)) return 1;

  /* Create the Rule structure */

  _recordWrite("%d", getModeRule());
  _recordWrite("#", "Type of Rule");
  _recordWrite("%lf", getRho());
  _recordWrite("#", "Correlation coefficient between GRFs");
  double slope = (FFFF(getSlope())) ? 0. : getSlope();
  _recordWrite("%lf", slope);
  double shdown = (FFFF(getShDown())) ? 0. : getShDown();
  _recordWrite("%lf", shdown);
  double shdsup = (FFFF(getShDsup())) ? 0. : getShDsup();
  _recordWrite("%lf", shdsup);
  _recordWrite("#", "Parameters for Shadow option");
  _recordWrite("%lf", getShift(0));
  _recordWrite("%lf", getShift(1));
  _recordWrite("%lf", getShift(2));
  _recordWrite("#", "Parameters for Shift option");

  /* Count the number of nodes */

  statistics(0,&nb_node,&nfacies,&nmax_tot,&ny1_tot,&ny2_tot,&prop_tot);
  _recordWrite("%d", nb_node);
  _recordWrite("#", "Number of nodes");

  /* Fill the nodes characteristics recursively */

  rank = 0;
  _ruleDefine(getMainNode(), 0, 0, 0, &rank);

  _fileClose(verbose);

  return 0;
}

void Rule::_ruleDefine(const Node *node,
                       int from_type,
                       int from_rank,
                       int from_vers,
                       int *rank) const
{
  int cur_rank;

  /* Calling node */

  _recordWrite("%d", from_type);
  _recordWrite("%d", from_rank);
  _recordWrite("%d", from_vers);

  /* Current node */

  _recordWrite("%d", node->getOrient());
  if (node->getFacies() <= 0)
  {
    cur_rank = *rank = (*rank) + 1;
    _recordWrite("%d", cur_rank);
    _recordWrite("%d", 0);
  }
  else
  {
    cur_rank = *rank;
    _recordWrite("%d", cur_rank);
    _recordWrite("%d", node->getFacies());
  }

  /* Comment */

  _recordWrite("#", "Node characteristics");

  if (node->getR1() != nullptr)
    _ruleDefine(node->getR1(), node->getOrient(), cur_rank, 1, rank);
  if (node->getR2() != nullptr)
    _ruleDefine(node->getR2(), node->getOrient(), cur_rank, 2, rank);
}

VectorString Rule::_buildNodNames(int nfacies)
{
  VectorString nodnames;

  for (int i = 1; i < nfacies; i++)
    nodnames.push_back("S");

  for (int i = 0; i < nfacies; i++)
    nodnames.push_back(incrementStringVersion("F",i+1,""));

  return nodnames;
}
