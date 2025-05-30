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

#include "Basic/AStringable.hpp"
#include "Basic/NamingConvention.hpp"
#include "Matrix/MatrixDense.hpp"

#include "Fractures/FracDesc.hpp"
#include "Fractures/FracEnviron.hpp"

#define NPART 5
#define NBYFRAC 7
#define NBYWOUT 8

class DbGrid;

class GSTLEARN_EXPORT FracList: public AStringable
{
public:
  FracList(int ndisc = 1000,
           bool flag_check = true,
           double low0 = EPSILON8,
           double low1 = EPSILON6,
           double eps = EPSILON3);
  FracList(const FracList& r);
  FracList& operator=(const FracList& r);
  virtual ~FracList();

  /// Interface for AStringable
  virtual String toString(const AStringFormat* strfmt = nullptr) const override;

  int getNFracs() const { return (int) _descs.size(); }

  int simulate(const FracEnviron& envir,
               bool flag_sim_layer,
               bool flag_sim_fract,
               int seed,
               bool verbose,
               const VectorDouble& elevations = VectorDouble());
  void addDescription(const FracDesc& description = FracDesc());

  MatrixDense fractureExport() const;
  MatrixDense layinfoExport() const { return _layinfo; };
  static FracList* fractureImport(const VectorDouble& frac_segs,
                                  const VectorDouble& layinfo = VectorDouble(),
                                  int nfamilies = 0);
  int fractureToBlock(DbGrid *dbgrid,
                      double xmax,
                      VectorDouble& permtab,
                      double perm_mat,
                      double perm_bench,
                      int ndisc = 1000.,
                      bool verbose = false,
                      const NamingConvention& namconv = NamingConvention("Fractures"));
  VectorDouble fractureToWell(int nval,
                              const VectorDouble& well,
                              double xmax,
                              const VectorDouble& permtab,
                              int *nint,
                              int *ncol);
  int fractureWellToBlock(DbGrid *dbgrid,
                          int col_perm,
                          int col_fluid,
                          int flag_fluid,
                          double val_fluid,
                          const VectorDouble& wellout,
                          int nval,
                          int ndisc = 1000.,
                          bool verbose = false);
  VectorDouble fractureExtractLength(int ifam, double cote, double dcote);
  VectorDouble fractureExtractDist(int ifam, double cote, double dcote);

  // Pipe for Class Description
  void setFamily(int i, int ifam) { _descs[i].setFamily(ifam); }
  void setOrient(int i, double orient) { _descs[i].setOrient(orient); }
  void addPoint(int i, double xx, double yy) { _descs[i].addPoint(xx,yy); }

private:
  static int _getRank(int ifam, int shift) { return (1 + ifam * NPART + shift); }
  void _setMemLayer(int i, double value)             { _layinfo.setValue(i,0,value); }
  void _setMemTheta1(int i, int ifam, double value)  { _layinfo.setValue(i,_getRank(ifam,0),value); }
  void _setMemTheta2(int i, int ifam, double value)  { _layinfo.setValue(i,_getRank(ifam,1),value); }
  void _setMemPropsur(int i, int ifam, double value) { _layinfo.setValue(i,_getRank(ifam,2),value); }
  void _setMemFrac(int i, int ifam, double value)    { _layinfo.setValue(i,_getRank(ifam,3),value); }
  void _setMemTotal(int i, int ifam, double value)   { _layinfo.setValue(i,_getRank(ifam,4),value); }
  double _getMemLayer(int i)                         { return _layinfo.getValue(i,0); }

  VectorDouble _layersManage(const FracEnviron& envir, double *y0) const;
  VectorDouble _layersRead(const VectorDouble& elevations, double *y0) const;
  int _fracAdd(int ifrac,
               int ifam,
               double xx,
               double cote,
               double thick,
               double orient,
               double* xp);
  void _checkFractureIntersect(double cote, int ifrac0);
  bool _belongToLayer(const FracDesc& desc,
                      double cote,
                      double *xd,
                      double *yd,
                      double *xe,
                      double *ye) const;
  double _layerIntensity(const FracFamily& family,
                         double thick) const;
  void _generateDensity(const FracEnviron& envir,
                        const FracFamily& family,
                        int ifam,
                        double cote,
                        VectorDouble& denstab);
  void _correctDensity(const FracFamily& family,
                       int ifam,
                       double cote,
                       VectorDouble& denstab);
  double _deriveIntensity(double theta1,
                          double thetap,
                          double propsur) const;
  double _extendFractures(const FracFamily& family,
                          int ifam,
                          double cote,
                          double thick,
                          VectorDouble& denstab);
  static bool _sameFaultSide(const FracEnviron& envir, int ifault0, double x0);
  static double _densityUpdate(const FracFault& fault, int side, int ifam, double cote, double xx);
  double _densityCumulate(const VectorDouble& denstab) const;
  bool _noRoomForMoreFracture(const VectorDouble& denstab) const;
  void _updateRepulsion(double x0, double range, VectorDouble& denstab);
  static bool _fractureInterrupt(const FracFamily& family,
                                 const FracDesc& desc,
                                 double thick);
  static double _faultAbscissae(const FracFault& fault, double cote);
  static double _cubic(double h);
  static double _fractureExtension(const FracDesc& desc, double cote, double dcote);
  int _simulateFractures(const FracEnviron& envir,
                         const FracFamily& family,
                         int ifam,
                         double cote,
                         double thick,
                         double theta,
                         VectorDouble& denstab);
  int _getDiscretizedRank(double cumdens, const VectorDouble& denstab) const;
  int _getEndPointCount() const;
  bool _isValidDisc(int idisc) const;

  static void _plungeSegment(DbGrid* dbgrid,
                             int iptr,
                             double delta,
                             double value,
                             double x1,
                             double y1,
                             double x2,
                             double y2);
  static void _welloutAdd(VectorDouble& wellout,
                          double x,
                          double y,
                          int ifrac,
                          int ip,
                          int family,
                          double perm);
  static void _trajAdd(VectorDouble& traj, double x, double y);
  static void _plungeSegmentGradual(DbGrid* dbgrid,
                                    int iptr,
                                    double delta,
                                    VectorDouble& traj,
                                    double perm1,
                                    double perm2,
                                    double range);

private:
  // Array of fracture descriptions
  std::vector<FracDesc> _descs;
  MatrixDense _layinfo;
  int _nlayers;
  // The number of discretization steps used to establish the fracture density
  int _ndisc;
  // The option for checking the fracture intersect or not
  bool _flagCheck;
  // The lower value used for checking the fracture repulsion for the central cell
  double _low0;
  // The lower value used for checking the fracture repulsion for the peripheral cells
  double _low1;
  double _xorigin;
  double _step;
  double _eps;
  bool _verbose;
};
