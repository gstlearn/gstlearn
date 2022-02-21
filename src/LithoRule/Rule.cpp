/******************************************************************************/
/* COPYRIGHT ARMINES, ALL RIGHTS RESERVED                                     */
/*                                                                            */
/* THE CONTENT OF THIS WORK CONTAINS CONFIDENTIAL AND PROPRIETARY             */
/* INFORMATION OF ARMINES. ANY DUPLICATION, MODIFICATION,                  prop   */
/* DISTRIBUTION, OR DISCLOSURE IN ANY FORM, IN WHOLE, OR IN PART, IS STRICTLY */
/* PROHIBITED WITHOUT THE PRIOR EXPRESS WRITTEN PERMISSION OF ARMINES         */
/*                                                                            */
/* TAG_SOURCE_CG                                                              */
/******************************************************************************/
#include "geoslib_f.h"
#include "geoslib_old_f.h"
#include "geoslib_enum.h"

#include "Basic/Utilities.hpp"
#include "Basic/String.hpp"
#include "Basic/AException.hpp"
#include "Basic/OptDbg.hpp"
#include "Model/Model.hpp"
#include "LithoRule/Rule.hpp"
#include "LithoRule/RuleStringFormat.hpp"
#include "LithoRule/Node.hpp"
#include "Db/Db.hpp"

#include <sstream>

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

Rule::Rule(double rho)
    : AStringable(),
      ASerializable(),
      _modeRule(ERule::STD),
      _flagProp(0),
      _rho(rho),
      _mainNode(nullptr)
{
}

Rule::Rule(const Rule& m)
    : AStringable(m),
      ASerializable(m),
      _modeRule(m._modeRule),
      _flagProp(m._flagProp),
      _rho(m._rho),
      _mainNode(new Node(*m._mainNode))
{
}

Rule& Rule::operator=(const Rule& m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    ASerializable::operator=(m);
    _modeRule = m._modeRule;
    _flagProp = m._flagProp;
    _rho = m._rho;
    _mainNode = new Node(*m._mainNode);
  }
  return *this;
}

Rule::~Rule()
{
  _clear();
}

void Rule::_clear()
{
  if (_mainNode != nullptr)
    delete _mainNode;
}

/**
 * This entry is specific for the Rule inference procedure
 * @param n_type Vector of items: 0 for Facies, 1 for G1 threshold; 2 for G2 threshold
 * @param n_facs Vector of facies: XXX for Facies, 0 for thresholds
 * @param rho GRFs correlation coefficient
 * @return
 */
int Rule::resetFromNumericalCoding(const VectorInt& n_type, const VectorInt& n_facs, double rho)
{
  _clear();
  _modeRule = ERule::STD;
  _rho = rho;
  setMainNodeFromNodNames(n_type, n_facs);
  return 0;
}

int Rule::resetFromFaciesCount(int nfacies, double rho)
{
  _clear();
  _modeRule = ERule::STD;
  _rho = rho;
  VectorString nodnames = buildNodNames(nfacies);
  setMainNodeFromNodNames(nodnames);
  return 0;
}

int Rule::resetFromNames(const VectorString& nodnames, double rho)
{
  _clear();
  _modeRule = ERule::STD;
  _rho = rho;
  setMainNodeFromNodNames(nodnames);
  return 0;
}

int Rule::resetFromCodes(const VectorInt& nodes, double rho)
{
  _clear();
  _modeRule = ERule::STD;
  _rho = rho;
  setMainNodeFromNodNames(nodes);
  return 0;
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
int Rule::setMainNodeFromNodNames(const VectorInt& nodes)
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

String Rule::toString(const AStringFormat* strfmt) const
{
  std::stringstream sstr;
  int node_tot,nfac_tot,nmax_tot,ny1_tot,ny2_tot;
  double prop_tot;

  const RuleStringFormat* rulefmt = dynamic_cast<const RuleStringFormat*>(strfmt);
  RuleStringFormat dsf;
  if (rulefmt != nullptr) dsf = *rulefmt;

  sstr << toTitle(0,"Lithotype Rule");

  if (statistics(0, &node_tot, &nfac_tot, &nmax_tot, &ny1_tot, &ny2_tot,
                 &prop_tot)) return sstr.str();
  if (prop_tot <= 0.)
  {
    dsf.setFlagProp(0);
    dsf.setFlagThresh(0);
  }

  sstr << "- Number of nodes               = " << node_tot << std::endl;
  sstr << "- Number of facies              = " << nfac_tot << std::endl;
  sstr << "- Number of thresholds along G1 = " << ny1_tot  << std::endl;
  sstr << "- Number of thresholds along G2 = " << ny2_tot  << std::endl;

  sstr << displaySpecific();

  sstr << std::endl;
  sstr << _mainNode->nodePrint(dsf.getFlagProp(), dsf.getFlagThresh());

  return sstr.str();
}

String Rule::displaySpecific() const
{
  std::stringstream sstr;
  if (ABS(_rho) > 0.)
    sstr << "Correlation between the two GRFs = " << _rho << std::endl;
  if (_rho == 1.)
    sstr << "(As the correlation is set to 1, only the first GRF is used)" << std::endl;
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
  if (getModeRule() == ERule::SHADOW || getModeRule() == ERule::SHIFT) ny2 = 0;
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
int Rule::particularities(Db* /*db*/,
                          const Db* /*dbprop*/,
                          Model* /*model*/,
                          int /*flag_grid_check*/,
                          int /*flag_stat*/) const
{
  return(0);
}

bool Rule::checkModel(const Model* /*model*/, int /*nvar*/) const
{
  return true;
}

void Rule::updateShift() const
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
  {
    decodeInList(symbol,nodes[i],&n_type[i],&n_facs[i]);

    // Check that the Facies rank is defined
    if (n_type[i] == 0)
    {
      if (n_facs[i] <= 0)
      {
        messerr("The Rule definition using 'nodnames' is incorrect");
        messerr("The element (%d) refers to a Facies with no Number",i+1);
      }
    }
  }
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
int Rule::setProportions(const VectorDouble& proportions) const
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

  if (OptDbg::query(EDbg::PROPS))
  {
    RuleStringFormat rulefmt(1);
    display(&rulefmt);
  }

  return(0);
}

int Rule::_deserialize(std::istream& is, bool /*verbose*/)
{
  int nb_node;

  /* Create the Rule structure */

  int mrule = 0;
  bool ret = _recordRead<int>(is, "Rule definition", mrule);
  ret = ret && _recordRead<double>(is, "Correlation Coefficient of GRFs", _rho);
  if (! ret) return 1;
  _modeRule = ERule::fromValue(mrule);

  // Specific case

  if (_deserializeSpecific(is)) return 1;

  /* Read the number of nodes */

  ret = ret && _recordRead<int>(is, "Number of Rule Nodes", nb_node);
  if (! ret) return 1;
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
      ret = ret && _recordRead<int>(is, "Rule Node Definition", nodes[lec++]);
  setMainNodeFromNodNames(nodes);

  return 0;
}

int Rule::_serialize(std::ostream& os, bool /*verbose*/) const
{
  int nb_node, nfacies, nmax_tot, ny1_tot, ny2_tot, rank;
  double prop_tot;

  /* Create the Rule structure */

  bool ret = _recordWrite<int>(os, "Type of Rule", getModeRule().getValue());
  ret = ret && _recordWrite<double>(os, "Correlation coefficient between GRFs", getRho());

  // Specific parameters

  _serializeSpecific(os);

  /* Count the number of nodes */

  statistics(0,&nb_node,&nfacies,&nmax_tot,&ny1_tot,&ny2_tot,&prop_tot);
  ret = ret && _recordWrite<int>(os, "Number of nodes", nb_node);

  /* Fill the nodes characteristics recursively */

  rank = 0;
  _ruleDefine(os, getMainNode(), 0, 0, 0, &rank);

  return ret ? 0 : 1;
}

void Rule::_ruleDefine(std::ostream& os,
                       const Node *node,
                       int from_type,
                       int from_rank,
                       int from_vers,
                       int *rank) const
{
  int cur_rank;

  /* Calling node */

  bool ret = _recordWrite<int>(os, "", from_type);
  ret = ret && _recordWrite<int>(os, "", from_rank);
  ret = ret && _recordWrite<int>(os, "", from_vers);

  /* Current node */

  ret = ret && _recordWrite<int>(os, "", node->getOrient());
  if (node->getFacies() <= 0)
  {
    cur_rank = *rank = (*rank) + 1;
    ret = ret && _recordWrite<int>(os, "", cur_rank);
    ret = ret && _recordWrite<int>(os, "", 0);
  }
  else
  {
    cur_rank = *rank;
    ret = ret && _recordWrite(os, "", cur_rank);
    ret = ret && _recordWrite(os, "", node->getFacies());
  }

  /* Comment */

  _commentWrite(os, "Node characteristics");

  if (node->getR1() != nullptr)
    _ruleDefine(os, node->getR1(), node->getOrient(), cur_rank, 1, rank);
  if (node->getR2() != nullptr)
    _ruleDefine(os, node->getR2(), node->getOrient(), cur_rank, 2, rank);
}

VectorString Rule::buildNodNames(int nfacies)
{
  VectorString nodnames;

  for (int i = 1; i < nfacies; i++)
    nodnames.push_back("S");

  for (int i = 0; i < nfacies; i++)
    nodnames.push_back(incrementStringVersion("F",i+1,""));

  return nodnames;
}

void Rule::setMainNodeFromNodNames(const VectorInt& n_type,
                                   const VectorInt& n_facs)
{
  int ipos = 0;
  int n_fac = 0;
  int n_y1 = 0;
  int n_y2 = 0;
  _mainNode = new Node("main", n_type, n_facs, &ipos, &n_fac, &n_y1, &n_y2);
}

void Rule::setMainNodeFromNodNames(const VectorString& nodnames)
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

/****************************************************************************/
/*!
**  Combine the underlying GRF into a facies value at data points
**
** \return  Error return code
**
** \param[in]  propdef    Props structure
** \param[in]  dbin       Db input structure
** \param[in]  dbout      Db output structure
** \param[in]  flag_used  1 if the gaussian is used; 0 otherwise
** \param[in]  ipgs       Rank of the PGS
** \param[in]  isimu      Rank of the simulation
** \param[in]  nbsimu     Number of simulations
**
** \remark Attributes ELoc::GAUSFAC are mandatory
** \remark Attributes ELoc::FACIES are mandatory
**
*****************************************************************************/
int Rule::gaus2facData(PropDef* propdef,
                       Db* dbin,
                       Db* /*dbout*/,
                       int *flag_used,
                       int ipgs,
                       int isimu,
                       int nbsimu)
{
  double y[2],facies,t1min,t1max,t2min,t2max;

  /* Initializations */

  check_mandatory_attribute("rule_gaus2fac_data",dbin,ELoc::GAUSFAC);

  /* Processing the translation */

  for (int iech=0; iech<dbin->getSampleNumber(); iech++)
  {
    if (! dbin->isActive(iech)) continue;

    /* Initializations */

    facies = TEST;
    for (int igrf=0; igrf<2; igrf++) y[igrf] = TEST;
    if (rule_thresh_define(propdef,dbin,this,ITEST,
                           iech,isimu,nbsimu,1,
                           &t1min,&t1max,&t2min,&t2max)) return 1;

    for (int igrf=0; igrf<2; igrf++)
    {
      int icase = get_rank_from_propdef(propdef,ipgs,igrf);
      y[igrf] = (flag_used[igrf]) ?
        dbin->getSimvar(ELoc::GAUSFAC,iech,isimu,0,icase,nbsimu,1) : 0.;
    }
    facies = getFaciesFromGaussian(y[0],y[1]);

    /* Combine the underlying GRFs to derive Facies */

    dbin->setSimvar(ELoc::FACIES,iech,isimu,0,ipgs,nbsimu,1,facies);
  }
  return 0;
}

/****************************************************************************/
/*!
**  Combine the underlying GRF into a facies value
**
** \return  Error return code
**
** \param[in]  propdef    Props structure
** \param[in]  dbout      Db output structure
** \param[in]  flag_used  1 if the gaussian is used; 0 otherwise
** \param[in]  ipgs       Rank of the PGS
** \param[in]  isimu      Rank of the simulation
** \param[in]  nbsimu     Number of simulations
**
** \remark Attributes ELoc::FACIES and ELoc::SIMU are mandatory
**
*****************************************************************************/
int Rule::gaus2facResult(PropDef  *propdef,
                         Db *dbout,
                         int *flag_used,
                         int ipgs,
                         int isimu,
                         int nbsimu) const
{
  int    ndim,iech,igrf,icase;
  double t1min,t1max,t2min,t2max,facies,y[2];

  /* Initializations */

  check_mandatory_attribute("rule_gaus2fac_result",dbout,ELoc::FACIES);
  check_mandatory_attribute("rule_gaus2fac_result",dbout,ELoc::SIMU);
  ndim   = dbout->getNDim();
  VectorDouble xyz(ndim);

  /* Processing the translation */

  for (iech=0; iech<dbout->getSampleNumber(); iech++)
  {
    if (! dbout->isActive(iech)) continue;

    /* Initializations */

    facies = TEST;
    for (igrf=0; igrf<2; igrf++) y[igrf] = TEST;

    if (rule_thresh_define(propdef,dbout,this,ITEST,
                           iech,isimu,nbsimu,1,
                           &t1min,&t1max,&t2min,&t2max)) return 1;
    for (igrf=0; igrf<2; igrf++)
    {
      icase = get_rank_from_propdef(propdef,ipgs,igrf);
      y[igrf] = (flag_used[igrf]) ?
          dbout->getSimvar(ELoc::SIMU,iech,isimu,0,icase,nbsimu,1) : 0.;
    }
    facies = getFaciesFromGaussian(y[0],y[1]);

    /* Combine the underlying GRFs to derive Facies */

    dbout->setSimvar(ELoc::FACIES,iech,isimu,0,ipgs,nbsimu,1,facies);
  }
  return 0;
}

/****************************************************************************/
/*!
**  Check if the current replicate can be added
**
** \return  1 if the point is a duplicate; 0 otherwise
**
** \param[in]  dbin       Db structure
** \param[in]  dbout      Db output structure
** \param[in]  jech       Rank of the replicate
**
*****************************************************************************/
int Rule::replicateInvalid(Db *dbin, Db *dbout, int jech) const
{
  for (int iech=0; iech<jech; iech++)
  {
    bool similar = false;
    for (int idim=0; idim<dbin->getNDim() && ! similar; idim++)
    {
      double delta = ABS(dbin->getCoordinate(iech, idim) -
                         dbin->getCoordinate(jech, idim));
      if (delta >= dbout->getUnit(idim)) similar = true;
    }
    if (similar)
    {
      message("Replicate invalid\n");
      return(1);
    }
  }
  return(0);
}

/****************************************************************************/
/*!
**  Set the bounds and possibly add replicates
**
** \return  Error return code
**
** \param[in]  propdef    PropDef structure
** \param[in]  dbin       Db structure
** \param[in]  dbout      Db grid structure
** \param[in]  isimu      Rank of the simulation (if EProcessOper::CONDITIONAL)
** \param[in]  igrf       Rank of the GRF
** \param[in]  ipgs       Rank of the GS
** \param[in]  nbsimu     Number of simulations (if EProcessOper::CONDITIONAL)
**
*****************************************************************************/
int Rule::evaluateBounds(PropDef* propdef,
                         Db* dbin,
                         Db* /*dbout*/,
                         int isimu,
                         int igrf,
                         int ipgs,
                         int nbsimu) const
{
  int    iech,nadd,nech,facies,nstep;
  double t1min,t1max,t2min,t2max;

  /* Initializations */

  if (dbin == nullptr) return(0);
  nadd = nstep = 0;
  nech = dbin->getSampleNumber();

  /* Dispatch */

  for (iech = 0; iech < nech; iech++)
  {
    if (!dbin->isActive(iech)) continue;
    facies = (int) dbin->getVariable(iech, 0);
    if (rule_thresh_define(propdef, dbin, this, facies, iech, isimu, nbsimu, 1,
                           &t1min, &t1max, &t2min, &t2max)) return (1);
    if (igrf == 0)
    {
      dbin->setLowerBound(iech, get_rank_from_propdef(propdef, ipgs, igrf),
                          t1min);
      dbin->setUpperBound(iech, get_rank_from_propdef(propdef, ipgs, igrf),
                          t1max);
    }
    else
    {
      dbin->setLowerBound(iech, get_rank_from_propdef(propdef, ipgs, igrf),
                          t2min);
      dbin->setUpperBound(iech, get_rank_from_propdef(propdef, ipgs, igrf),
                          t2max);
    }
  }

  if (igrf == 0 && nadd > 0)
  {
    message("Initial count of data = %d\n",nech);
    message("Number of replicates  = %d\n",nadd);
  }
  return(0);
}

int Rule::dumpToNF(const String& neutralFilename, bool verbose) const
{
  std::ofstream os;
  int ret = 1;

  if (_fileOpenWrite(neutralFilename, "Rule", os, verbose))
  {
    ret = _serialize(os, verbose);
    if (ret && verbose) messerr("Problem writing in the Neutral File.");
    os.close();
  }
  return ret;
}

Rule* Rule::create(double rho)
{
  return new Rule(rho);
}

Rule* Rule::createFromNF(const String& neutralFilename, bool verbose)
{
  Rule* rule = nullptr;
  std::ifstream is;
  if (_fileOpenRead(neutralFilename, "Rule", is, verbose))
  {
    rule = new Rule();
    rule->setModeRule(ERule::STD);
    if (rule->_deserialize(is, verbose))
    {
      if (verbose) messerr("Problem reading the Neutral File.");
      delete rule;
      rule = nullptr;
    }
    is.close();
  }
  return rule;
}

Rule* Rule::createFromNames(const VectorString& nodnames,double rho)
{
  Rule* rule = new Rule();
  if (rule->resetFromNames(nodnames, rho))
  {
    messerr("Problem when creating Rule from a vector of Names");
    delete rule;
    return nullptr;
  }
  return rule;
}
Rule* Rule::createFromCodes(const VectorInt& nodes,double rho)
{
  Rule* rule = new Rule();
  if (rule->resetFromCodes(nodes, rho))
  {
    messerr("Problem when creating Rule from a vector of Codes");
    delete rule;
    return nullptr;
  }
  return rule;
}
Rule* Rule::createFromNumericalCoding(const VectorInt& n_type,
                                      const VectorInt& n_facs,
                                      double rho)
{
  Rule* rule = new Rule();
  if (rule->resetFromNumericalCoding(n_type, n_facs, rho))
  {
    messerr("Problem when creating Rule from Numerical Coding");
    delete rule;
    return nullptr;
  }
  return rule;
}
Rule* Rule::createFromFaciesCount(int nfacies, double rho)
{
  Rule* rule = new Rule();
  if (rule->resetFromFaciesCount(nfacies, rho))
  {
    messerr("Problem when creating Rule from a number of Facies");
    delete rule;
    return nullptr;
  }
  return rule;
}
