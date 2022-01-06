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
#include "Basic/Law.hpp"
#include "LithoRule/Node.hpp"
#include "Basic/MathFunc.hpp"
#include "LithoRule/Rule.hpp"
#include "geoslib_f.h"
#include "geoslib_enum.h"

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

Node::Node(const String& nodnam, int orient, int facies)
    : AStringable(),
      _nodnam(nodnam),
      _r1(nullptr),
      _r2(nullptr),
      _orient(orient),
      _facies((orient == THRESH_IDLE) ? facies : ITEST),
      _prop(0.),
      _thresh(0.),
      _p1(0.),
      _p2(0.),
      _t1min(0.),
      _t1max(0.),
      _t2min(0.),
      _t2max(0.),
      _cdf1min(0.),
      _cdf1max(0.),
      _cdf2min(0.),
      _cdf2max(0.)
{

}

Node::Node(const String& nodnam, const VectorInt& n_type, const VectorInt& n_facs,
           int *ipos, int *n_fac, int *n_y1, int *n_y2)
    : AStringable(),
      _nodnam(nodnam),
      _r1(nullptr),
      _r2(nullptr),
      _orient(0),
      _facies(0),
      _prop(0.),
      _thresh(0.),
      _p1(0.),
      _p2(0.),
      _t1min(0.),
      _t1max(0.),
      _t2min(0.),
      _t2max(0.),
      _cdf1min(0.),
      _cdf1max(0.),
      _cdf2min(0.),
      _cdf2max(0.)
{
  std::stringstream sstr;

  /* Decode the next operation */

  int rank = 0;
  int jpos = *ipos;
  int type = n_type[jpos];
  (*ipos)++;

  /* Compose the name of the node */

  switch (type)
  {
    case THRESH_IDLE:
      rank = ++(*n_fac);
      _facies = n_facs[jpos];
      break;

    case THRESH_Y1:
      rank = ++(*n_y1);
      _facies = 0;
      break;

    case THRESH_Y2:
      rank = ++(*n_y2);
      _facies = 0;
      break;
  }
  _orient = type;
  sstr << symbol[type] << rank;
  _nodnam = sstr.str();

  /* Create the node */

  switch (type)
  {
    case THRESH_Y1:
    case THRESH_Y2:
      sstr << _nodnam << "[Low]";
      _r1 = new Node(sstr.str(),n_type,n_facs,ipos,n_fac,n_y1,n_y2);
      if (_r1 == nullptr) return;
      sstr << _nodnam << "[Sup]";
      _r2 = new Node(sstr.str(),n_type,n_facs,ipos,n_fac,n_y1,n_y2);
      if (_r2 == nullptr) return;
      break;
  }
  return;
}

Node::Node(bool /*flagShadow*/)
    : AStringable(),
      _nodnam("S1"),
      _r1(nullptr),
      _r2(nullptr),
      _orient(THRESH_Y1),
      _facies(0),
      _prop(0.),
      _thresh(0.),
      _p1(0.),
      _p2(0.),
      _t1min(0.),
      _t1max(0.),
      _t2min(0.),
      _t2max(0.),
      _cdf1min(0.),
      _cdf1max(0.),
      _cdf2min(0.),
      _cdf2max(0.)
{
  _r1 = new Node("T1",THRESH_Y2,0);
  _r1->_r1 = new Node("F3",THRESH_IDLE,SHADOW_SHADOW);
  _r1->_r2 = new Node("F2",THRESH_IDLE,SHADOW_WATER);
  _r2 =  new Node("F1",THRESH_IDLE,SHADOW_ISLAND);
}

Node::Node(const Node& m)
    : AStringable(m),
      _nodnam(m._nodnam),
      _r1(nullptr),
      _r2(nullptr),
      _orient(m._orient),
      _facies(m._facies),
      _prop(m._prop),
      _thresh(m._thresh),
      _p1(m._p1),
      _p2(m._p2),
      _t1min(m._t1min),
      _t1max(m._t1max),
      _t2min(m._t2min),
      _t2max(m._t2max),
      _cdf1min(m._cdf1min),
      _cdf1max(m._cdf1max),
      _cdf2min(m._cdf2min),
      _cdf2max(m._cdf2max)
{
  if (m._r1 != nullptr)
    _r1 = new Node(*m._r1);
  if (m._r2 != nullptr)
    _r2 = new Node(*m._r2);
}

Node& Node::operator=(const Node& m)
{
  if (this != &m)
  {
    AStringable::operator=(m);
    _nodnam = m._nodnam;
    if (m._r1 != nullptr)
      _r1 = new Node(*m._r1);
    else
      _r1 = nullptr;
    if (m._r2 != nullptr)
      _r2 = new Node(*m._r2);
    else
      _r2 = nullptr;
    _orient = m._orient;
    _facies = m._facies;
    _prop = m._prop;
    _thresh = m._thresh;
    _p1 = m._p1;
    _p2 = m._p2;
    _t1min = m._t1min;
    _t1max = m._t1max;
    _t2min = m._t2min;
    _t2max = m._t2max;
    _cdf1min = m._cdf1min;
    _cdf1max = m._cdf1max;
    _cdf2min = m._cdf2min;
    _cdf2max = m._cdf2max;
  }
  return *this;
}

Node::~Node()
{
  if (_r1 != (Node *) NULL)
  {
    delete _r1;
    _r1 = nullptr;
  }
  if (_r2 != (Node *) NULL)
  {
    delete _r2;
    _r2 = nullptr;
  }
}

/****************************************************************************/
/*!
**  Recursively print the Node
**
** \param[in]  flagProp   true if the Proportions must be printed
** \param[in]  flagThresh true if the Threshold must be printed
**
*****************************************************************************/
String Node::nodePrint(bool flagProp, bool flagThresh) const

{
  std::stringstream sstr;
  if (_r1 != (Node *) NULL)
    sstr << _r1->nodePrint(flagProp, flagThresh);
  if (_r2 != (Node *) NULL)
    sstr << _r2->nodePrint(flagProp, flagThresh);

  switch (_orient)
  {
    case THRESH_IDLE:
      sstr << "Node " << _nodnam << " - Facies " << _facies;

      if (flagProp)
      {
        sstr << " - Proportion = " << _prop << std::endl;
        if (flagThresh)
        {
          sstr << "          Y1 in [" << _t1min << " ; " << _t1max
               << "]" << std::endl;
          sstr << "          Y2 in [" << _t2min << " ; " << _t2max
               << "]" << std::endl;
        }
      }
      else
        sstr << std::endl;
      break;
  }
  return sstr.str();
}

/****************************************************************************/
/*!
**  Recursively print the Node
**
** \param[in]  flagProp   true if the Proportions must be printed
** \param[in]  flagThresh true if the Threshold must be printed
**
*****************************************************************************/
String Node::nodePrintShadow(bool flagProp, bool flagThresh) const
{
  std::stringstream sstr;
  if (_r1 != (Node *) NULL)
    sstr << _r1->nodePrintShadow(flagProp, flagThresh);
  if (_r2 != (Node *) NULL)
    sstr << _r2->nodePrintShadow(flagProp, flagThresh);

  switch (_orient)
  {
    case THRESH_IDLE:
      if (_facies == SHADOW_ISLAND)
        sstr << "Node " << _nodnam << " - Island";
      if (_facies == SHADOW_WATER)
        sstr << "Node " << _nodnam << " - Water";
      if (_facies == SHADOW_SHADOW)
        sstr << "Node " << _nodnam << " - Shadow";

      if (flagProp)
      {
        sstr << " - Proportion = " << _prop << std::endl;
        if (flagThresh)
        {
          sstr << "            Y1 in [" << _t1min << " ; " << _t1max
               << "]" << std::endl;
          sstr << "            Y@ in [" << _t2min << " ; " << _t2max
               << "]" << std::endl;
        }
      }
      else
        sstr << std::endl;
      break;
  }
  return sstr.str();
}

/****************************************************************************/
/*!
**  Calculates the statistics from a given node recursively
**
** \param[in,out]  node_tot Number of nodes
** \param[in,out]  nfac_tot Number of facies
** \param[in,out]  ny1_tot  Number of thresholds for Y1
** \param[in,out]  ny2_tot  Number of thresholds for Y2
** \param[in,out]  prop_tot Total proportion
**
** \remark  The variables 'nfac_tot' and 'prop_tot' must set to zero
** \remark  as arguments of the call
**
*****************************************************************************/
void Node::getStatistics(int *node_tot,
                         int *nfac_tot,
                         int *ny1_tot,
                         int *ny2_tot,
                         double *prop_tot)
{
  *node_tot = 0;
  *nfac_tot = 0;
  *ny1_tot  = 0;
  *ny2_tot  = 0;
  *prop_tot = 0.;
  _getStatistics(node_tot,nfac_tot,ny1_tot,ny2_tot,prop_tot);
}

void Node::_getStatistics(int *node_tot,
                          int *nfac_tot,
                          int *ny1_tot,
                          int *ny2_tot,
                          double *prop_tot)
{
  double p1_loc,p2_loc;

  p1_loc = 0.;
  if (_r1 != (Node *) NULL)
    _r1->_getStatistics(node_tot,nfac_tot,ny1_tot,ny2_tot,&p1_loc);

  p2_loc = 0.;
  if (_r2 != (Node *) NULL)
    _r2->_getStatistics(node_tot,nfac_tot,ny1_tot,ny2_tot,&p2_loc);

  _p1 = p1_loc;
  _p2 = p2_loc;
  (*prop_tot) += p1_loc + p2_loc;
  (*node_tot) += 1;

  switch (_orient)
  {
    case THRESH_IDLE:
      (*nfac_tot) += 1;
      (*prop_tot) += _prop;
      break;

    case THRESH_Y1:
      (*ny1_tot) += 1;
      break;

    case THRESH_Y2:
      (*ny2_tot) += 1;
      break;
  }
  return;
}

/****************************************************************************/
/*!
**  Count the presence of each facies
**
** \return  Error return code
**
** \param[in]  facies Array for counting the presence of facies
**
*****************************************************************************/
int Node::isValid(VectorInt& facies)
{
  if (_r1 != (Node *) NULL)
  {
    if (_r1->isValid(facies)) return(1);
  }
  if (_r2 != (Node *) NULL)
  {
    if (_r2->isValid(facies)) return(1);
  }

  if (_orient != THRESH_IDLE) return(0);
  int nfac = static_cast<int> (facies.size());
  if (IFFFF(_facies))
  {
    messerr("The facies of node %s has not been defined",_nodnam.c_str());
    return(1);
  }
  if (_facies < 1 || _facies > nfac)
  {
    messerr("Error in the facies rank (%d) at node %s: it should lie within [1,%d]",
            _facies,_nodnam.c_str(),nfac);
    return(1);
  }
  facies[_facies - 1]++;
  return(0);
}

/****************************************************************************/
/*!
**  Scale the proportions of all facies
**  new_prop = old_prop / scale
**
** \param[in]  scale Scaling factor
**
*****************************************************************************/
void Node::scaleProp(double scale)
{
  if (_r1 != (Node *) NULL) _r1->scaleProp(scale);
  if (_r2 != (Node *) NULL) _r2->scaleProp(scale);

  if (_orient == THRESH_IDLE) _prop /= scale;
  return;
}

String Node::toString(const AStringFormat* /*strfmt*/) const
{
  std::stringstream sstr;
  sstr << nodePrint(false, false);
  return sstr.str();
}

int Node::proportionDefine(const VectorDouble& props)
{
  static double eps = 1.e-3;

  if (_r1 != (Node *) NULL)
  {
    if (_r1->proportionDefine(props)) return(1);
  }
  if (_r2 != (Node *) NULL)
  {
    if (_r2->proportionDefine(props)) return(1);
  }

  /* Assign the proportion of the current facies */

  if (_orient == THRESH_IDLE) // [FO] 21/09/13 => nothing to do in splitting nodes (facies=0)
  {
    int facies = _facies;
    if (IFFFF(facies)) return(0);

    double propval = props[facies - 1];
    if (! FFFF(propval))
    {
      if (propval < (0. - eps) || propval > (1. + eps))
      {
        messerr("Wrong proportion for facies %d (%lf): it should lie in [0,1]",
                facies,propval);
        return(1);
      }
      if (propval < 0.) propval = 0.;
      if (propval > 1.) propval = 1.;
    }
    _prop = propval;
  }
  return(0);
}

/**
 * For a given facies, return the Proportion
 * @param facies Target facies (starting from 1)
 * @param prop   Returned proportion
 * @return 1 if the facies is found; 0 otherwise
 */
int Node::getProportion(int facies, double *prop)
{
  if (_r1 != (Node *) NULL)
  {
    if (_r1->getProportion(facies,prop)) return(1);
  }
  if (_r2 != (Node *) NULL)
  {
    if (_r2->getProportion(facies,prop)) return(1);
  }
  if (_facies != facies) return(0);
  *prop = _prop;
  return(1);
}

/****************************************************************************/
/*!
**  Get the threshold for a facies
**
** \return  1 if the facies is found; 0 otherwise
**
** \param[in]  mode   Stopping criterion
** \li                 1 : when the requested facies is found
** \li                 2 : when the rectangle rank is reached
** \param[in]  istop  Value of the stopping criterion
** \param[in]  rank   Index of the rectangle (input)
**
** \param[out]  facies Value of the facies
** \param[out]  t1min  Minimum threshold along Y1
** \param[out]  t1max  Maximum threshold along Y1
** \param[out]  t2min  Minimum threshold along Y2
** \param[out]  t2max  Maximum threshold along Y2
**
** \remark Argument 'rank' must be initialized to 0 in the calling function
**
*****************************************************************************/
int Node::getThresh(int  mode,
                    int  istop,
                    int *rank,
                    int *facies,
                    double *t1min,
                    double *t1max,
                    double *t2min,
                    double *t2max)
{
  if (_r1 != (Node *) NULL)
  {
    if (_r1->getThresh(mode, istop, rank, facies, t1min, t1max, t2min, t2max))
      return (1);
  }
  if (_r2 != (Node *) NULL)
  {
    if (_r2->getThresh(mode, istop, rank, facies, t1min, t1max, t2min, t2max))
      return (1);
  }

  /* Stopping criterion */

  if (mode == 1)
  {
    if (_facies != istop) return(0);
  }
  else
  {
    if (! (IFFFF(_facies))) (*rank)++;
    if ((*rank) != istop) return(0);
  }

  *facies = _facies;
  *t1min  = _t1min;
  *t1max  = _t1max;
  *t2min  = _t2min;
  *t2max  = _t2max;
  return(1);
}

/****************************************************************************/
/*!
**  Recursively deduce the threshold from cumulative proportions
**
** \param[in]  rho   Correlation between the GRFs
** \param[in]  t1min Minimum bound for the first GRF
** \param[in]  t1max Maximum bound for the first GRF
** \param[in]  t2min Minimum bound for the second GRF
** \param[in]  t2max Maximum bound for the second GRF
**
*****************************************************************************/
void Node::proportionToThresh(double rho,
                              double t1min,
                              double t1max,
                              double t2min,
                              double t2max)
{
  _t1min   = t1min;
  _t1max   = t1max;
  _t2min   = t2min;
  _t2max   = t2max;
  _cdf1min = _transform(-1,t1min);
  _cdf1max = _transform(-1,t1max);
  _cdf2min = _transform(-1,t2min);
  _cdf2max = _transform(-1,t2max);

  _thresh = _threshFromPropcum(rho);

  if (_orient == THRESH_Y1)
  {
    if (_r1 != (Node *) NULL)
      _r1->proportionToThresh(rho, t1min, _thresh, t2min, t2max);
    if (_r2 != (Node *) NULL)
      _r2->proportionToThresh(rho, _thresh, t1max, t2min, t2max);
  }
  else
  {
    if (_r1 != (Node *) NULL)
      _r1->proportionToThresh(rho, t1min, t1max, t2min, _thresh);
    if (_r2 != (Node *) NULL)
      _r2->proportionToThresh(rho, t1min, t1max, _thresh, t2max);
  }
  return;
}

/****************************************************************************/
/*!
**  Perform the transform between real and gaussian scales
**
** \return  Returned value
**
** \param[in]  mode    <0 from gaussian to real; >0 from real to gaussian
** \param[in]  value   Input value
**
*****************************************************************************/
double Node::_transform(int mode, double value)
{
  if (mode < 0)
  {
    if (get_rule_mode())
      return(law_cdf_gaussian(value));
    else
      return(value);
  }
  else
  {
    if (get_rule_mode())
      return(law_invcdf_gaussian(value));
    else
      return(value);
  }
}

/****************************************************************************/
/*!
**  Derive a threshold from the bounds and the lower and upper
**  cumulative proportions
**
** return  Threshold value or TEST if an error occurs
**
** \param[in]  rho  Correlation between the GRFs
**
*****************************************************************************/
double Node::_threshFromPropcum(double rho)
{
  double gval,sump;
  static double eps_small = 1.e-04;

  /* Initializations */

  if (_orient == THRESH_IDLE)
  {
    return TEST;
  }

  if (rho == 0.)
  {
    sump = _p1 + _p2;

    /* Case of two uncorrelated GRFs */

    if (_orient == THRESH_Y1)
    {
      if (ABS(sump) > eps_small)
        gval = (_cdf1min * _p2 + _cdf1max * _p1) / sump;
      else
        gval = _cdf1min;
    }
    else
    {
      if (ABS(sump) > eps_small)
        gval = (_cdf2min * _p2 + _cdf2max * _p1) / sump;
      else
        gval = _cdf2min;
    }

    if (gval <      eps_small) gval = 0.;
    if (gval > 1. - eps_small) gval = 1.;
    return _transform(1,gval);
  }
  else
  {

    /* Case of two correlated GRFs */

    return _threshDichotomy(rho);
  }
}

/****************************************************************************/
/*!
**  Find a threshold by dichotomy method
**  This method is restricted to the case of correlated GRFs
**
** \return  The requested threshold
**
** \param[in]  rho  Correlation between the two GRFs (in ]0,1[)
**
*****************************************************************************/
double Node::_threshDichotomy(double rho)
{
  double low[2],sup[2],mini[2],maxi[2],prop,error;
  int    i,infin[2],ier,ind;
  static double abseps = 1.e-08;
  static double releps = 0.;
  static int maxpts    = 40000;

  /* Initializations */

  low[0] = mini[0] = _t1min;
  low[1] = mini[1] = _t2min;
  sup[0] = maxi[0] = _t1max;
  sup[1] = maxi[1] = _t2max;
  for ( i = 0 ; i < 2 ; i++ ) infin[i] = mvndst_infin(low[i],sup[i]);
  ind = (_orient == THRESH_Y1) ? 0 : 1;

  /* Convergence loop */

  do
  {
    sup[ind]   = (mini[ind] + maxi[ind]) / 2.;
    if (sup[ind] >= get_rule_extreme(+1)) break;
    infin[ind] = mvndst_infin(low[ind],sup[ind]);
    mvndst(2,low,sup,infin,&rho,maxpts,abseps,releps,&error,&prop,&ier);
    if (ier) messageAbort("Fatal error in mvndst");
    if (prop < _p1)
      mini[ind] = sup[ind];
    else
      maxi[ind] = sup[ind];
  } while (ABS(_p1 - prop) > abseps);

  return(sup[ind]);
}

/****************************************************************************/
/*!
**  Convert the two underlying GRFs into facies
**
** \return  1 if the facies is found; 0 otherwise
**
** \param[in]  y1     Value of the first underlying GRF
** \param[in]  y2     Value of the second underlying GRF
**
** \param[out]  facies Facies value
**
*****************************************************************************/
int Node::gaussianToFacies(double y1, double y2, double *facies)
{
  if (_r1 != (Node *) NULL)
  {
    if (_r1->gaussianToFacies(y1,y2,facies)) return(1);
  }
  if (_r2 != (Node *) NULL)
  {
    if (_r2->gaussianToFacies(y1,y2,facies)) return(1);
  }

  if (_orient != THRESH_IDLE) return(0);
  if (_t1min > get_rule_extreme(-1) && y1 < _t1min) return(0);
  if (_t1max < get_rule_extreme(+1) && y1 > _t1max) return(0);
  if (_t2min > get_rule_extreme(-1) && y2 < _t2min) return(0);
  if (_t2max < get_rule_extreme(+1) && y2 > _t2max) return(0);
  *facies = (double) _facies;
  return(1);
}

/****************************************************************************/
/*!
**  Recursive call for encoding the node characteristics
**
** \param[in]  nodes       Array for node characteristics
** \param[in]  parent_type Type of the parent (THRESH_IDLE, THRESH_Y1 or Y2)
** \param[in]  parent_rank Rank of the parent
** \param[in]  parent_vers Orientation for the parent
** \param[in]  rank        Node number in the node pile
** \param[in]  n_fac       Number of the facies
** \param[in]  n_y1        Number of the threshold along Y1
** \param[in]  n_y2        Number of the threshold along Y2
**
*****************************************************************************/
void Node::_getInfo(int *nodes,
                    int parent_type,
                    int parent_rank,
                    int parent_vers,
                    int *rank,
                    int *n_fac,
                    int *n_y1,
                    int *n_y2) const
{
  int type,number;

  /* Load the parameters of the parent node */

  NODES(*rank,0) = parent_type;
  NODES(*rank,1) = parent_rank;
  NODES(*rank,2) = parent_vers;

  /* Load the parameters of the current node */

  NODES(*rank,3) = _orient;
  switch (_orient)
  {
    case THRESH_IDLE:
      NODES(*rank,4) = ++(*n_fac);
      NODES(*rank,5) = _facies;
      break;

    case THRESH_Y1:
      NODES(*rank,4) = ++(*n_y1);
      NODES(*rank,5) = 0;
      break;

    case THRESH_Y2:
      NODES(*rank,4) = ++(*n_y2);
      NODES(*rank,5) = 0;
      break;
  }
  type   = NODES(*rank,3);
  number = NODES(*rank,4);

  /* Process the subsequent tree */

  if (_r1 != (Node *) NULL)
  {
    (*rank)++;
    _r1->_getInfo(nodes,type,number,1,rank,n_fac,n_y1,n_y2);
  }
  if (_r2 != (Node *) NULL)
  {
    (*rank)++;
    _r2->_getInfo(nodes,type,number,2,rank,n_fac,n_y1,n_y2);
  }
  return;
}

void Node::getInfo(int *nodes) const
{
  int rank = 0;
  int n_fac = 0;
  int n_y1 = 0;
  int n_y2 = 0;
  _getInfo(nodes, 0, 0, 0, &rank, &n_fac, &n_y1, &n_y2);
}
