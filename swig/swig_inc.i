%feature(director) IProj;
%feature(director) ACov;
%feature(director) AFunction;
%feature(director) AFunctional;
%feature(director) ABiTargetCheck;

%{
  #include "gstlearn_export.hpp"
  #include "geoslib_define.h"
  #include "geoslib_enum.h"
  #include "geoslib_d.h"
  #include "geoslib_f.h"
  #include "geoslib_old_f.h"

  #include "Transform/ATransform.hpp"
  #include "Transform/ATransformWithAutoDiff.hpp"
  #include "Transform/TuckeyGH.hpp"
  #include "Transform/YeoJohnson.hpp"

  #include "Enum/AEnum.hpp"
  #include "Enum/EKrigOpt.hpp"
  #include "Enum/ESPDECalcMode.hpp"
  #include "Enum/EAnam.hpp"
  #include "Enum/ECst.hpp"
  #include "Enum/EDbg.hpp"
  #include "Enum/ELaw.hpp"
  #include "Enum/EShape.hpp"
  #include "Enum/EConvDir.hpp"
  #include "Enum/ECalcVario.hpp"
  #include "Enum/EConvType.hpp"
  #include "Enum/ECov.hpp"
  #include "Enum/ECSV.hpp"
  #include "Enum/ETape.hpp"
  #include "Enum/ELoadBy.hpp"
  #include "Enum/ELoc.hpp"
  #include "Enum/EOperator.hpp"
  #include "Enum/EPowerPT.hpp"
  #include "Enum/ERule.hpp"
  #include "Enum/EConsElem.hpp"
  #include "Enum/EConsType.hpp"
  #include "Enum/EModelProperty.hpp"
  #include "Enum/EMorpho.hpp"
  #include "Enum/ENeigh.hpp"
  #include "Enum/ESpaceType.hpp"
  #include "Enum/ESelectivity.hpp"
  #include "Enum/EStatOption.hpp"
  #include "Enum/EDirGen.hpp"
  #include "Enum/EGaussInv.hpp"
  #include "Enum/ECalcMember.hpp"
  #include "Enum/EPostUpscale.hpp"
  #include "Enum/EPostStat.hpp"
  #include "Enum/ESimuType.hpp"
  #include "Enum/ERole.hpp"

  #include "Basic/VectorT.hpp"
  #include "Basic/VectorNumT.hpp"
  #include "Basic/ICloneable.hpp"
  #include "Basic/VectorHelper.hpp"
  #include "Basic/AFunctional.hpp"
  #include "Basic/AFunction.hpp"
  #include "Basic/ArgumentTest.hpp"
  #include "Basic/AStringable.hpp"
  #include "Basic/AStringFormat.hpp"
  #include "Basic/ASerializable.hpp"
  #include "Basic/Tensor.hpp"
  #include "Basic/Grid.hpp"
  #include "Basic/String.hpp"
  #include "Basic/Interval.hpp"
  #include "Basic/Limits.hpp"
  #include "Basic/Utilities.hpp"
  #include "Basic/CSVformat.hpp"
  #include "Basic/FunctionalSpirale.hpp"
  #include "Basic/NamingConvention.hpp"
  #include "Basic/OptDbg.hpp"
  #include "Basic/OptCst.hpp"
  #include "Basic/OptCustom.hpp"
  #include "Basic/File.hpp"
  #include "Basic/Limits.hpp"
  #include "Basic/Plane.hpp"
  #include "Basic/FFT.hpp"
  #include "Basic/PolyLine2D.hpp"
  #include "Basic/Law.hpp"
  #include "Basic/LawStable.hpp"
  #include "Basic/MathFunc.hpp"
  #include "Basic/Indirection.hpp"
  #include "Basic/Message.hpp"
  #include "Basic/WarningMacro.hpp"

  #include "Geometry/GeometryHelper.hpp"
  #include "Geometry/Rotation.hpp"
  #include "Geometry/ABiTargetCheck.hpp"
  #include "Geometry/BiTargetCheckBench.hpp"
  #include "Geometry/BiTargetCheckCell.hpp"
  #include "Geometry/BiTargetCheckDistance.hpp"
  #include "Geometry/BiTargetCheckFaults.hpp"
  #include "Geometry/BiTargetCheckCode.hpp"
  #include "Geometry/BiTargetCheckDate.hpp"
  #include "Geometry/BiTargetCheckGeometry.hpp"

  #include "Arrays/AArray.hpp"
  #include "Arrays/Array.hpp"
  #include "Arrays/BImage.hpp"
  #include "Arrays/BImageStringFormat.hpp"

  #include "Faults/Faults.hpp"

  #include "Boolean/ShapeParameter.hpp"
  #include "Boolean/AShape.hpp"
  #include "Boolean/ShapeParallelepiped.hpp"
  #include "Boolean/ShapeEllipsoid.hpp"
  #include "Boolean/ShapeParaboloid.hpp"
  #include "Boolean/ShapeHalfEllipsoid.hpp"
  #include "Boolean/ShapeHalfParaboloid.hpp"
  #include "Boolean/ShapeHalfSinusoid.hpp"
  #include "Boolean/ModelBoolean.hpp"

  #include "Space/ASpace.hpp"
  #include "Space/SpaceComposite.hpp"
  #include "Space/ASpaceObject.hpp"
  #include "Space/SpacePoint.hpp"
  #include "Space/SpaceTarget.hpp"
  #include "Space/SpaceRN.hpp"
  #include "Space/SpaceSN.hpp"
  #include "Space/SpaceShape.hpp"

  #include "Skin/ISkinFunctions.hpp"
  #include "Skin/Skin.hpp"

  #include "Calculators/ACalculator.hpp"
  #include "Calculators/ACalcDbVarCreator.hpp"
  #include "Calculators/ACalcDbToDb.hpp"
  #include "Calculators/CalcMigrate.hpp"
  #include "Calculators/ACalcInterpolator.hpp"
  #include "Calculators/CalcStatistics.hpp"
  #include "Calculators/CalcGridToGrid.hpp"
  #include "Calculators/CalcSimuPost.hpp"
  #include "Calculators/CalcSimuPostDemo.hpp"
  #include "Calculators/CalcSimuPostPropByLayer.hpp"

  #include "Mesh/AMesh.hpp"
  #include "Mesh/MeshEFaulted.hpp"
  #include "Mesh/MeshEStandard.hpp"
  #include "Mesh/MeshETurbo.hpp"
  #include "Mesh/MeshSpherical.hpp"
  #include "Mesh/MeshSphericalExt.hpp"
  #include "Mesh/VectorMeshes.hpp"

  #include "Polynomials/APolynomial.hpp"
  #include "Polynomials/ClassicalPolynomial.hpp"
  #include "Polynomials/Hermite.hpp"
  #include "Polynomials/MonteCarlo.hpp"

  #include "LinearOp/CGParam.hpp"
  #include "LinearOp/LogStats.hpp"
  #include "LinearOp/ALinearOp.hpp"
  #include "LinearOp/APreconditioner.hpp"
  #include "LinearOp/MultiGridSolver.hpp"
  #include "LinearOp/MultiGridSPDE.hpp"
  #include "LinearOp/ASimulable.hpp"
  #include "LinearOp/ASimulableMatrix.hpp"
  #include "LinearOp/LinearOpCGSolver.hpp"
  #include "LinearOp/ALinearOpMulti.hpp"
  #include "LinearOp/AShiftOp.hpp"
  #include "LinearOp/ShiftOpStencil.hpp"
  #include "LinearOp/ShiftOpMatrix.hpp"
  #include "LinearOp/IPrecisionOp.hpp"
  #include "LinearOp/PrecisionOp.hpp"
  #include "LinearOp/PrecisionOpMatrix.hpp"
  #include "LinearOp/TurboOptimizer.hpp"
  #include "LinearOp/IProj.hpp"
  #include "LinearOp/ScaleOp.hpp"
  #include "LinearOp/ProjMatrix.hpp"
  #include "LinearOp/ProjMulti.hpp"
  #include "LinearOp/ProjMultiMatrix.hpp"
  #include "LinearOp/ProjZero.hpp"
  #include "LinearOp/ProjComposition.hpp"
  #include "LinearOp/PrecisionOpMulti.hpp"
  #include "LinearOp/PrecisionOpMultiMatrix.hpp"
  #include "LinearOp/IOptimCost.hpp"
  #include "LinearOp/OptimCostBinary.hpp"
  #include "LinearOp/OptimCostColored.hpp"
  #include "LinearOp/ProjConvolution.hpp"
  #include "LinearOp/SPDEOp.hpp"
  #include "LinearOp/SPDEOpMatrix.hpp"
  #include "LinearOp/MatrixSymmetricSim.hpp"
  #include "LinearOp/ACholesky.hpp"
  #include "LinearOp/CholeskyDense.hpp"
  #include "LinearOp/CholeskySparse.hpp"
  #include "LinearOp/LinearOpHelper.hpp"

  #include "Neigh/ANeigh.hpp"
  #include "Neigh/NeighUnique.hpp"
  #include "Neigh/NeighImage.hpp"
  #include "Neigh/NeighMoving.hpp"
  #include "Neigh/NeighBench.hpp"
  #include "Neigh/NeighCell.hpp"

  #include "Variogram/Variograms.hpp"
  #include "Variogram/AVario.hpp"
  #include "Variogram/VarioParam.hpp"
  #include "Variogram/Vario.hpp"
  #include "Variogram/DirParam.hpp"
  #include "Variogram/VMap.hpp"
  #include "Variogram/VCloud.hpp"
  #include "Variogram/VarioOrder.hpp"

  #include "Basic/ParamInfo.hpp"
  #include "Basic/ListParams.hpp"
  #include "Model/ModelGeneric.hpp"
  #include "Model/GaussianProcess.hpp"
  #include "Model/ModelCovList.hpp"
  #include "Model/Model.hpp"
  #include "Model/ModelOptimParam.hpp"
  #include "Model/ElemNostat.hpp"
  #include "Model/Option_AutoFit.hpp"
  #include "Model/Option_VarioFit.hpp"
  #include "Model/Constraints.hpp"
  #include "Model/ConsItem.hpp"
  #include "Model/CovParamId.hpp"

  #include "Covariances/ParamId.hpp"
  #include "Covariances/TabNoStat.hpp"
  #include "Covariances/TabNoStatCovAniso.hpp"
  #include "Covariances/TabNoStatSills.hpp"
  #include "Covariances/ANoStat.hpp"
  #include "Covariances/NoStatArray.hpp"
  #include "Covariances/NoStatOnMesh.hpp"
  #include "Covariances/NoStatFunctional.hpp"
  #include "Covariances/ACov.hpp"
  #include "Covariances/CovBase.hpp"
  #include "Covariances/CovProportional.hpp"
  #include "Covariances/AKernel.hpp"
  #include "Covariances/CovAnisoList.hpp"
  #include "Covariances/CovAnisoList.hpp"
  #include "Covariances/CovAniso.hpp"
  #include "Covariances/CovGradientGeneric.hpp"
  #include "Covariances/CovGradientAnalytic.hpp"
  #include "Covariances/CorAniso.hpp"
  #include "Covariances/CorFactorized.hpp"
  #include "Covariances/CorGaussianMixture.hpp"
  #include "Covariances/CorGneiting.hpp"
  #include "Covariances/CorMatern.hpp"
  #include "Covariances/CovLMCTapering.hpp"
  #include "Covariances/CovLMCConvolution.hpp"
  #include "Covariances/CovLMCAnamorphosis.hpp"
  #include "Covariances/CovContext.hpp"
  #include "Covariances/CovCalcMode.hpp"
  #include "Covariances/KernelBesselJ.hpp"
  #include "Covariances/KernelMatern.hpp"
  #include "Covariances/KernelCauchy.hpp"
  #include "Covariances/KernelCauchyGen.hpp"
  #include "Covariances/KernelCosExp.hpp"
  #include "Covariances/KernelCosinus.hpp"
  #include "Covariances/KernelCubic.hpp"
  #include "Covariances/KernelExponential.hpp"
  #include "Covariances/KernelGamma.hpp"
  #include "Covariances/KernelGaussian.hpp"
  #include "Covariances/KernelGC1.hpp"
  #include "Covariances/KernelGC3.hpp"
  #include "Covariances/KernelGC5.hpp"
  #include "Covariances/KernelGCspline2.hpp"
  #include "Covariances/KernelGCspline.hpp"
  #include "Covariances/KernelLinear.hpp"
  #include "Covariances/KernelNugget.hpp"
  #include "Covariances/KernelMarkov.hpp"
  #include "Covariances/KernelPenta.hpp"
  #include "Covariances/KernelPower.hpp"
  #include "Covariances/KernelReg1D.hpp"
  #include "Covariances/KernelSincard.hpp"
  #include "Covariances/KernelSpherical.hpp"
  #include "Covariances/KernelStable.hpp"
  #include "Covariances/KernelStorkey.hpp"
  #include "Covariances/KernelTriangle.hpp"
  #include "Covariances/KernelWendland0.hpp"
  #include "Covariances/KernelWendland1.hpp"
  #include "Covariances/KernelWendland2.hpp"
  #include "Covariances/KernelGeometric.hpp"
  #include "Covariances/KernelPoisson.hpp"
  #include "Covariances/KernelLinearSph.hpp"
  #include "Covariances/CovDiffusionAdvection.hpp"
  #include "Covariances/CovHelper.hpp"

  #include "Drifts/ADrift.hpp"
  #include "Drifts/DriftList.hpp"
  #include "Drifts/DriftM.hpp"
  #include "Drifts/DriftF.hpp"
  #include "Drifts/DriftFactory.hpp"

  #include "Matrix/AMatrix.hpp"
  #include "Matrix/MatrixDense.hpp"
  #include "Matrix/MatrixSparse.hpp"
  #include "Matrix/MatrixSquare.hpp"
  #include "Matrix/NF_Triplet.hpp"
  #include "Matrix/MatrixSymmetric.hpp"
  #include "Matrix/MatrixFactory.hpp"
  #include "Matrix/MatrixInt.hpp"
  #include "Matrix/Table.hpp"
  #include "MLayers/MLayers.hpp"
  #include "LinearOp/InvNuggetOp.hpp"

  #include "API/SPDE.hpp"
  #include "API/TestInheritance.hpp"
  #include "API/Style.hpp"
  #include "API/SPDEParam.hpp"
  #include "API/Potential.hpp"

  #include "Db/Db.hpp"
  #include "Db/DbGrid.hpp"
  #include "Db/DbLine.hpp"
  #include "Db/DbGraphO.hpp"
  #include "Db/DbMeshTurbo.hpp"
  #include "Db/DbMeshStandard.hpp"
  #include "Db/DbStringFormat.hpp"
  #include "Db/DbHelper.hpp"
  #include "Db/RankHandler.hpp"

  #include "DataBase/ColID.hpp"
  #include "DataBase/RoleID.hpp"
  #include "DataBase/DbCol.hpp"
  #include "DataBase/DbData.hpp"
  #include "DataBase/Category.hpp"
  #include "DataBase/Dictionary.hpp"
  #include "DataBase/VectorCategory.hpp"

  #include "Anamorphosis/CalcAnamTransform.hpp"
  #include "Anamorphosis/AAnam.hpp"
  #include "Anamorphosis/AnamContinuous.hpp"
  #include "Anamorphosis/AnamDiscrete.hpp"
  #include "Anamorphosis/AnamUser.hpp"
  #include "Anamorphosis/AnamHermite.hpp"
  #include "Anamorphosis/AnamEmpirical.hpp"
  #include "Anamorphosis/AnamDiscreteDD.hpp"
  #include "Anamorphosis/AnamDiscreteIR.hpp"
  #include "Anamorphosis/PPMT.hpp"

  #include "Gibbs/GibbsMMulti.hpp"
  #include "Gibbs/GibbsUMulti.hpp"
  #include "Gibbs/GibbsMultiMono.hpp"
  #include "Gibbs/GibbsUMultiMono.hpp"
  #include "Gibbs/GibbsUPropMono.hpp"

  #include "Morpho/Morpho.hpp"

  #include "Polygon/Polygons.hpp"
  #include "Polygon/PolyElem.hpp"

  #include "Stats/Classical.hpp"
  #include "Stats/PCA.hpp"
  #include "Stats/PCAStringFormat.hpp"
  #include "Stats/Selectivity.hpp"
  #include "Stats/Regression.hpp"

  #include "LithoRule/Node.hpp"
  #include "LithoRule/Rule.hpp"
  #include "LithoRule/RuleStringFormat.hpp"
  #include "LithoRule/RuleProp.hpp"
  #include "LithoRule/RuleShift.hpp"
  #include "LithoRule/RuleShadow.hpp"
  #include "LithoRule/PropDef.hpp"

  #include "Estimation/Estimations.hpp"
  #include "Estimation/KrigingSystem.hpp"
  #include "Estimation/KrigingAlgebra.hpp"
  #include "Estimation/CalcKriging.hpp"
  #include "Estimation/CalcKrigingFactors.hpp"
  #include "Estimation/CalcKrigingGradient.hpp"
  #include "Estimation/CalcSimpleInterpolation.hpp"
  #include "Estimation/CalcImage.hpp"
  #include "Estimation/CalcGlobal.hpp"
  #include "Estimation/KrigOpt.hpp"
  #include "Estimation/AModelOptim.hpp"
  #include "Estimation/AModelOptimFactory.hpp"
  #include "Estimation/ALikelihood.hpp"
  #include "Estimation/Vecchia.hpp"
  #include "Estimation/Likelihood.hpp"

  #include "OutputFormat/AOF.hpp"
  #include "OutputFormat/FileLAS.hpp"
  #include "OutputFormat/FileVTK.hpp"
  #include "OutputFormat/GridArcGis.hpp"
  #include "OutputFormat/GridBmp.hpp"
  #include "OutputFormat/GridEclipse.hpp"
  #include "OutputFormat/GridF2G.hpp"
  #include "OutputFormat/GridIfpEn.hpp"
  #include "OutputFormat/GridIrap.hpp"
  #include "OutputFormat/GridXYZ.hpp"
  #include "OutputFormat/GridZycor.hpp"
  #include "OutputFormat/segy.h"

  #include "Polynomials/Chebychev.hpp"

  #include "Simulation/Simulations.hpp"
  #include "Simulation/ACalcSimulation.hpp"
  #include "Simulation/ACalcSimuGaussian.hpp"
  #include "Simulation/CalcSimuTurningBands.hpp"
  #include "Simulation/CalcSimuPGS.hpp"
  #include "Simulation/TurningBandDirection.hpp"
  #include "Simulation/TurningBandOperate.hpp"
  #include "Simulation/CalcSimuSpectral.hpp"
  #include "Simulation/CalcSimuBoolean.hpp"
  #include "Simulation/SimuSpectralRN.hpp"
  #include "Simulation/SimuSpectralS2.hpp"
  #include "Simulation/SpectrumOnRN.hpp"
  #include "Simulation/BooleanObject.hpp"
  #include "Simulation/SimuBooleanParam.hpp"
  #include "Simulation/SimuSpherical.hpp"
  #include "Simulation/SimuSphericalParam.hpp"
  #include "Simulation/CalcSimuSubstitution.hpp"
  #include "Simulation/SimuSubstitutionParam.hpp"
  #include "Simulation/CalcSimuPartition.hpp"
  #include "Simulation/SimuPartitionParam.hpp"
  #include "Simulation/SimuFFTParam.hpp"
  #include "Simulation/CalcSimuFFT.hpp"
  #include "Simulation/SimuRefineParam.hpp"
  #include "Simulation/CalcSimuRefine.hpp"
  #include "Simulation/CalcSimuEden.hpp"

  #include "Fractures/FracEnviron.hpp"
  #include "Fractures/FracFamily.hpp"
  #include "Fractures/FracFault.hpp"
  #include "Fractures/FracDesc.hpp"
  #include "Fractures/FracList.hpp"

  #include "Tree/Ball.hpp"
  #include "Tree/BallFaulted.hpp"
  #include "Tree/KNN.hpp"

  #include "Spatial/Projection.hpp"
  #include "Spatial/SpatialIndices.hpp"

  #include "PluriGaussian/TracePGS.hpp"
  #include "PluriGaussian/CorPGS.hpp"
  #include "PluriGaussian/DiscretePGS.hpp"
  #include "PluriGaussian/CalcModelPGS.hpp"

  #include "Core/Acknowledge.hpp"
  #include "Core/Seismic.hpp"

  #include "API/newAPIs.hpp"
  // Mask some warning generated by SWIG:
  //DISABLE_WARNING_DECLARATION_MASKED
  //DISABLE_WARNING_EXPR_COND_ASSIGNMENT
  //DISABLE_WARNING_UNREFERENCED_FORMAL_PARAMETER
   using namespace gstlrn;
%}

////////////////////////////
//        Typemaps        //
////////////////////////////

// Mandatory for using swig::asptr and swig::from for std::vectors
%include std_vector.i
%include std_string.i
%template(DoNotUseVectorIntStd)     std::vector< int >;
%template(DoNotUseVectorLongStd)    std::vector< long >;
%template(DoNotUseVectorLLongStd)   std::vector< long long >;
%template(DoNotUseVectorSizeT)      std::vector< size_t >; // Keep size_t here otherwise asptr fails!
%template(DoNotUseVectorDoubleStd)  std::vector< double >;
%template(DoNotUseVectorStringStd)  std::vector< std::string >; // Keep std::string here otherwise asptr fails!
%template(DoNotUseVectorFloatStd)   std::vector< float >;
%template(DoNotUseVectorUCharStd)   std::vector< unsigned char >; // Keep unsigned char here
%template(DoNotUseVectorBoolStd)    std::vector< bool >;
%template(DoNotUseVVectorIntStd)    std::vector< std::vector< int > >;
%template(DoNotUseVVectorLongStd)   std::vector< std::vector< long > >;
%template(DoNotUseVVectorLLongStd)  std::vector< std::vector< long long > >;
%template(DoNotUseVVectorDoubleStd) std::vector< std::vector< double > >;
%template(DoNotUseVVectorFloatStd)  std::vector< std::vector< float > >;

%template(VectorECov)              std::vector< gstlrn::ECov >;
%template(VectorEStatOption)       std::vector< gstlrn::EStatOption >;
%template(VectorESelectivity)      std::vector< gstlrn::ESelectivity >;
%template(VectorDirParam)          std::vector< gstlrn::DirParam >;
%template(VectorPolyElem)          std::vector< gstlrn::PolyElem >;
%template(VectorInterval)          std::vector< gstlrn::Interval >;
%template(VectorEPostStat)         std::vector< gstlrn::EPostStat >;
%template(VectorSpacePoint)        std::vector< gstlrn::SpacePoint >;
%template(VectorABiTargetCheck)    std::vector< gstlrn::ABiTargetCheck* >;
%template(VectorProjMatrix)        std::vector< gstlrn::ProjMatrix* >;
%template(VectorConstProjMatrix)   std::vector< const gstlrn::ProjMatrix*>;
%template(VectorConstIProj)        std::vector< const gstlrn::IProj*>;
%template(VVectorConstProjMatrix)  std::vector< std::vector< const gstlrn::ProjMatrix*> >;
%template(VVectorConstIProj)       std::vector< std::vector< const gstlrn::IProj*> >;
%template(VecMeshes)               std::vector< const gstlrn::AMesh*>;
%template(VectorMatrixSquare)      std::vector< gstlrn::MatrixSquare >;
%template(VectorRule)              std::vector< gstlrn::Rule>;

///////////////////////////////////////
// Conversion Target language => C++ //
///////////////////////////////////////
namespace gstlrn {

// Note : Before including this file :
//        - vectorToCpp, vectorVectorToCpp, matrixDenseToCpp, matrixSparseToCpp and convertToCpp
//          functions must be defined in ToCpp fragment

// Convert scalar arguments by value
%typemap(in, fragment="ToCpp") Id,
                               double,
                               String,
                               float,
                               UChar,
                               bool
{
  try
  {
    int errcode = convertToCpp($input, $1);
    if (!SWIG_IsOK(errcode))
      %argument_fail(errcode, "$type", $symname, $argnum);
  }
  catch(...)
  {
    messerr("Error while converting argument #$argnum of type '$type' in '$symname' function");
  }
}
%typemap(in, fragment="ToCpp") std::string_view (String tmp)
{
  try
  {
    int errcode = convertToCpp($input, tmp);
    $1 = tmp;
    if (!SWIG_IsOK(errcode))
      %argument_fail(errcode, "$type", $symname, $argnum);
  }
  catch(...)
  {
    messerr("Error while converting argument #$argnum of type '$type' in '$symname' function");
  }
}

// Convert scalar argument by reference
// Don't add String or char here otherwise "res2 not declared" / "alloc1 not declared"
%typemap(in, fragment="ToCpp") Id*     (Id val),     const Id*     (Id val),
                               Id&     (Id val),     const Id&     (Id val),
                               double* (double val), const double* (double val),
                               double& (double val), const double& (double val),
                               float*   (float val), const float*   (float val),
                               float&   (float val), const float&   (float val),
                               UChar*   (UChar val), const UChar*   (UChar val),
                               UChar&   (UChar val), const UChar&   (UChar val),
                               bool*     (bool val), const bool*     (bool val),
                               bool&     (bool val), const bool&     (bool val)
{
  try
  {
    int errcode = convertToCpp($input, val);
    if (!SWIG_IsOK(errcode))
      %argument_fail(errcode, "$type", $symname, $argnum);
    $1 = &val;
  }
  catch(...)
  {
    messerr("Error while converting argument #$argnum of type '$type' in '$symname' function");
  }
}

%typemap(in, fragment="ToCpp") VectorInt    (void *argp),
                               VectorDouble (void *argp),
                               VectorString (void *argp),
                               VectorFloat  (void *argp),
                               VectorUChar  (void *argp),
                               VectorBool   (void *argp)
{
  // Try to convert from any target language vector
  int errcode = vectorToCpp($input, $1);
  if (!SWIG_IsOK(errcode))
  {
    // Try direct conversion of Vectors by value (see swigtypes.swg)
    errcode = SWIG_ConvertPtr($input, &argp, $&descriptor, %convertptr_flags);
    if (SWIG_IsOK(errcode))
    {
      if (!argp) {
        %argument_nullref("$type", $symname, $argnum);
      }
      else {
        $&ltype temp = %reinterpret_cast(argp, $&ltype);
        $1 = *temp;
        if (SWIG_IsNewObj(errcode)) %delete(temp);
      }
    }
    else
      %argument_fail(errcode, "$type", $symname, $argnum);
  }
}

%typemap(in, fragment="ToCpp") VectorVectorInt    (void *argp),
                               VectorVectorDouble (void *argp),
                               VectorVectorFloat  (void *argp)
{
  // Try to convert from any target language vector
  int errcode = vectorVectorToCpp($input, $1);
  if (!SWIG_IsOK(errcode))
  {
    // Try direct conversion of VectorVectors by value (see swigtypes.swg)
    errcode = SWIG_ConvertPtr($input, &argp, $&descriptor, %convertptr_flags);
    if (SWIG_IsOK(errcode))
    {
      if (!argp) {
        %argument_nullref("$type", $symname, $argnum);
      }
      else {
        $&ltype temp = %reinterpret_cast(argp, $&ltype);
        $1 = *temp;
        if (SWIG_IsNewObj(errcode)) %delete(temp);
      }
    }
    else {
      %argument_fail(errcode, "$type", $symname, $argnum);
    }
  }
}

// Typemap for missing pointer types are treated separately below.
%typemap(in, fragment="ToCpp") const VectorString* (void *argp, VectorString vec),
                               const VectorFloat*  (void *argp, VectorFloat vec),
                               const VectorUChar*  (void *argp, VectorUChar vec),
                               const VectorBool*   (void *argp, VectorBool vec)
{
  // Try to convert from any target language vector
  int errcode = vectorToCpp($input, vec);
  if (errcode == SWIG_NullReferenceError)
  {
    $1 = nullptr;
  }
  else if (!SWIG_IsOK(errcode))
  {
    try
    {
      // Try direct conversion of Vectors by reference/pointer (see swigtypes.swg)
      errcode = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);
      if (SWIG_IsOK(errcode))
      {
        if (!argp) {
          %argument_nullref("$type", $symname, $argnum);
        }
        $1 = %reinterpret_cast(argp, $ltype);
      }
      else {
        %argument_fail(errcode, "$type", $symname, $argnum);
      }
    }
    catch(...)
    {
      %argument_fail(errcode, "$type", $symname, $argnum);
    }
  }
  else {
    $1 = &vec;
  }
}

// Special case for pointer to VectorDouble
// TODO: All the special cases for pointers are based upon a volunteer memory leak
// (in order to ensure that the local pointer is not immediately deleted).
// It could be solved by using a share_memory type of pointer as ...
// auto vec = std::make_shared<VectorDouble>(); A ameliorer
//
%typemap(in, fragment="ToCpp") const VectorDouble* (void *argp)
{
  // Try to convert from any target language vector
  auto* vec = new VectorDouble();
  int errcode = vectorToCpp($input, *vec);
  if (errcode == SWIG_NullReferenceError)
  {
    $1 = nullptr;
  }
  else if (!SWIG_IsOK(errcode))
  {
    try
    {
      // Try direct conversion of Vectors by reference/pointer (see swigtypes.swg)
      errcode = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);
      if (SWIG_IsOK(errcode))
      {
        if (!argp) {
          %argument_nullref("$type", $symname, $argnum);
        }
        $1 = %reinterpret_cast(argp, $ltype);
      }
      else {
        %argument_fail(errcode, "$type", $symname, $argnum);
      }
    }
    catch(...)
    {
      %argument_fail(errcode, "$type", $symname, $argnum);
    }
  }
  else {
    $1 = vec;
  }
}

// Special case for pointer to VectorInt
%typemap(in, fragment="ToCpp") const VectorInt*    (void *argp, VectorInt vec)
{
  // Try to convert from any target language vector

  auto* vec = new VectorInt();
  int errcode = vectorToCpp($input, *vec);
  if (errcode == SWIG_NullReferenceError)
  {
    $1 = nullptr;
  }
  else if (!SWIG_IsOK(errcode))
  {
    try
    {
      // Try direct conversion of VectorVectors by reference/pointer (see swigtypes.swg)
      errcode = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);
      if (SWIG_IsOK(errcode))
      {
        if (!argp) {
          %argument_nullref("$type", $symname, $argnum);
        }
        $1 = %reinterpret_cast(argp, $ltype);
      }
      else {
        %argument_fail(errcode, "$type", $symname, $argnum);
      }
    }
    catch(...)
    {
      %argument_fail(errcode, "$type", $symname, $argnum);
    }
  }
  else
  {
    $1 = vec;
  }
}

// Special case for pointer to VectorVectorInt
%typemap(in, fragment="ToCpp") const VectorVectorInt*    (void *argp, VectorVectorInt vec)
{
  // Try to convert from any target language vector
  auto* vec = new VectorVectorInt();
  int errcode = vectorVectorToCpp($input, *vec);
  if (errcode == SWIG_NullReferenceError)
  {
    $1 = nullptr;
  }
  else if (!SWIG_IsOK(errcode))
  {
    try
    {
      // Try direct conversion of VectorVectors by reference/pointer (see swigtypes.swg)
      errcode = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);
      if (SWIG_IsOK(errcode))
      {
        if (!argp) {
          %argument_nullref("$type", $symname, $argnum);
        }
        $1 = %reinterpret_cast(argp, $ltype);
      }
      else {
        %argument_fail(errcode, "$type", $symname, $argnum);
      }
    }
    catch(...)
    {
      %argument_fail(errcode, "$type", $symname, $argnum);
    }
  }
  else
  {
    $1 = vec;
  }
}

%typemap(in, fragment="ToCpp") const VectorInt&    (void *argp, VectorInt vec),
                               const VectorDouble& (void *argp, VectorDouble vec),
                               const VectorString& (void *argp, VectorString vec),
                               const VectorFloat&  (void *argp, VectorFloat vec),
                               const VectorUChar&  (void *argp, VectorUChar vec),
                               const VectorBool&   (void *argp, VectorBool vec)
{
  // Try to convert from any target language vector
  int errcode = vectorToCpp($input, vec);
  if (errcode == SWIG_NullReferenceError)
  {
    $1 = &vec;
  }
  else if (!SWIG_IsOK(errcode))
  {
    try
    {
      // Try direct conversion of Vectors by reference/pointer (see swigtypes.swg)
      errcode = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);
      if (SWIG_IsOK(errcode))
      {
        if (!argp) {
          %argument_nullref("$type", $symname, $argnum);
        }
        $1 = %reinterpret_cast(argp, $ltype);
      }
      else {
        %argument_fail(errcode, "$type", $symname, $argnum);
      }
    }
    catch(...)
    {
      %argument_fail(errcode, "$type", $symname, $argnum);
    }
  }
  else {
    $1 = &vec;
  }
}

%typemap(in, fragment="ToCpp") const VectorVectorInt&    (void *argp, VectorVectorInt vec),
                               const VectorVectorDouble& (void *argp, VectorVectorDouble vec),
                               const VectorVectorDouble* (void *argp, VectorVectorDouble vec),
                               const VectorVectorFloat&  (void *argp, VectorVectorFloat vec),
                               const VectorVectorFloat*  (void *argp, VectorVectorFloat vec)
{
  // Try to convert from any target language vector
  int errcode = vectorVectorToCpp($input, vec);
  if (errcode == SWIG_NullReferenceError)
  {
    $1 = nullptr;
  }
  else if (!SWIG_IsOK(errcode))
  {
    try
    {
      // Try direct conversion of VectorVectors by reference/pointer (see swigtypes.swg)
      errcode = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);
      if (SWIG_IsOK(errcode))
      {
        if (!argp) {
          %argument_nullref("$type", $symname, $argnum);
        }
        $1 = %reinterpret_cast(argp, $ltype);
      }
      else {
        %argument_fail(errcode, "$type", $symname, $argnum);
      }
    }
    catch(...)
    {
      %argument_fail(errcode, "$type", $symname, $argnum);
    }
  }
  else
  {
    $1 = &vec;
  }
}

%typemap(in, fragment="ToCpp")
VectorNumT<double>&& (void *argp, VectorNumT<double> vec),
VectorNumT<float>&& (void *argp, VectorNumT<float> vec),
VectorNumT<long long>&& (void *argp, VectorNumT<long long> vec),
VectorNumT<UChar>&& (void *argp, VectorNumT<UChar> vec),
VectorNumT<bool>&& (void *argp, VectorNumT<bool> vec),
VectorT<UChar>&& (void *argp, VectorT<UChar> vec)
{
  int errcode = vectorToCpp($input, vec);

  if (errcode == SWIG_NullReferenceError)
  {
    $1 = &vec;
  }
  else if (!SWIG_IsOK(errcode))
  {
    try
    {
      errcode = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);
      if (SWIG_IsOK(errcode))
      {
        if (!argp)
          %argument_nullref("$type", $symname, $argnum);

        $1 = %reinterpret_cast(argp, $ltype);
      }
      else
      {
        %argument_fail(errcode, "$type", $symname, $argnum);
      }
    }
    catch(...)
    {
      %argument_fail(errcode, "$type", $symname, $argnum);
    }
  }
  else
  {
    $1 = &vec;
  }
}

%typemap(in, fragment="ToCpp")
VectorT<std::string>&& (void *argp, VectorT<std::string> vec)
{
  int errcode = vectorToCpp($input, vec);

  if (errcode == SWIG_NullReferenceError)
  {
    $1 = &vec;
  }
  else if (!SWIG_IsOK(errcode))
  {
    try
    {
      errcode = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);

      if (SWIG_IsOK(errcode))
      {
        if (!argp)
          %argument_nullref("$type", $symname, $argnum);

        $1 = %reinterpret_cast(argp, $ltype);
      }
      else
      {
        %argument_fail(errcode, "$type", $symname, $argnum);
      }
    }
    catch(...)
    {
      %argument_fail(errcode, "$type", $symname, $argnum);
    }
  }
  else
  {
    $1 = &vec;
  }
}

%typemap(in, fragment="ToCpp")
VectorT<String> (VectorT<String> vec)
{
  int errcode = vectorToCpp($input, vec);

  if (!SWIG_IsOK(errcode))
    %argument_fail(errcode, "$type", $symname, $argnum);

  $1 = vec;
}

%typemap(in, fragment="ToCpp") const MatrixDense&     (void *argp, MatrixDense mat),
                               const MatrixDense*     (void *argp, MatrixDense mat),
                               const MatrixSquare&    (void *argp, MatrixSquare mat),
                               const MatrixSquare*    (void *argp, MatrixSquare mat),
                               const MatrixSymmetric& (void *argp, MatrixSymmetric mat),
                               const MatrixSymmetric* (void *argp, MatrixSymmetric mat)
{
  // Try to convert from any target language vector
  int errcode = matrixDenseToCpp($input, mat);
  if (errcode == SWIG_NullReferenceError)
  {
    $1 = nullptr;
  }
  else if (!SWIG_IsOK(errcode))
  {
    try
    {
      // Try direct conversion of Matrices by reference/pointer (see swigtypes.swg)
      errcode = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);
      if (SWIG_IsOK(errcode))
      {
        if (!argp) {
          %argument_nullref("$type", $symname, $argnum);
        }
        else
          $1 = %reinterpret_cast(argp, $ltype);
      }
      else {
        %argument_fail(errcode, "$type", $symname, $argnum);
      }
    }
    catch(...)
    {
      %argument_fail(errcode, "$type", $symname, $argnum);
    }
  }
  else
  {
    $1 = &mat;
  }
}

%typemap(in, fragment="ToCpp") const MatrixSparse&     (void *argp, MatrixSparse mat),
                               const MatrixSparse*     (void *argp, MatrixSparse mat)
{
  // Try to convert from any target language vector
  int errcode = matrixSparseToCpp($input, mat);
  if (errcode == SWIG_NullReferenceError)
  {
    $1 = nullptr;
  }
  else if (!SWIG_IsOK(errcode))
  {
    try
    {
      // Try direct conversion of Matrices by reference/pointer (see swigtypes.swg)
      errcode = SWIG_ConvertPtr($input, &argp, $descriptor, %convertptr_flags);
      if (SWIG_IsOK(errcode))
      {
        if (!argp) {
          %argument_nullref("$type", $symname, $argnum);
        }
        else
          $1 = %reinterpret_cast(argp, $ltype);
      }
      else {
        %argument_fail(errcode, "$type", $symname, $argnum);
      }
    }
    catch(...)
    {
      %argument_fail(errcode, "$type", $symname, $argnum);
    }
  }
  else
  {
    $1 = &mat;
  }
}

////////////////////////////////////////////////
// Conversion C++ => Target language
//
// Note : Before including this file :
//        - vectorFromCpp, vectorVectorFromCpp,
//        - matrixDenseFromCpp, matrixSparseFromCpp, objectFromCpp
//          functions must be defined in FromCpp fragment

%typemap(out, fragment="FromCpp") Id,
                                  double,
                                  String,
                                  float,
                                  UChar,
                                  bool
{
  $result = objectFromCpp($1);
}
%typemap(out, fragment="FromCpp") std::string_view
{
  String tmp{$1};
  $result = objectFromCpp(tmp);
}

%typemap(out, fragment="FromCpp") Id*,     const Id*,     Id&,     const Id&,
                                  double*, const double*, double&, const double&,
                                  String*, const String*, String&, const String&,
                                  float*,  const float*,  float&,  const float&,
                                  UChar*,  const UChar*,  UChar&,  const UChar&,
                                  bool*,   const bool*,   bool&,   const bool&
{
  $result = objectFromCpp(*$1);
}

%typemap(out, fragment="FromCpp") VectorBool
{
  int errcode = vectorFromCpp(&($result), $1);
  if (!SWIG_IsOK(errcode))
    SWIG_exception_fail(SWIG_ArgError(errcode),
                        "in method $symname, wrong return value: $type");
}

%typemap(out, fragment="FromCpp") VectorInt,
                                  VectorDouble,
                                  VectorString,
                                  VectorFloat,
                                  VectorUChar,
                                  VectorBool
{
  int errcode = vectorFromCpp(&($result), $1);
  if (!SWIG_IsOK(errcode))
    SWIG_exception_fail(SWIG_ArgError(errcode), "in method $symname, wrong return value: $type");
}

%typemap(out, fragment="FromCpp") VectorInt*,    VectorInt&,
                                  VectorDouble*, VectorDouble&,
                                  VectorString*, VectorString&,
                                  VectorFloat*,  VectorFloat&,
                                  VectorUChar*,  VectorUChar&,
                                  VectorBool*,   VectorBool&
{
  int errcode = vectorFromCpp(&($result), *$1);
  if (!SWIG_IsOK(errcode))
    SWIG_exception_fail(SWIG_ArgError(errcode), "in method $symname, wrong return value: $type");
}

%typemap(out, fragment="FromCpp") VectorVectorInt,
                                  VectorVectorDouble,
                                  VectorVectorFloat
{
  int errcode = vectorVectorFromCpp(&($result), $1);
  if (!SWIG_IsOK(errcode))
    SWIG_exception_fail(SWIG_ArgError(errcode), "in method $symname, wrong return value: $type");
}

%typemap(out, fragment="FromCpp") VectorVectorInt*,    VectorVectorInt&,
                                  VectorVectorDouble*, VectorVectorDouble&,
                                  VectorVectorFloat*,  VectorVectorFloat&
{
  int errcode = vectorVectorFromCpp(&($result), *$1);
  if (!SWIG_IsOK(errcode))
    SWIG_exception_fail(SWIG_ArgError(errcode), "in method $symname, wrong return value: $type");
}

//%typemap(out, fragment="FromCpp") MatrixDense,
//                                  MatrixSquare,
//                                  MatrixSymmetric
//{
//  int errcode = matrixDenseFromCpp(&($result), $1);
//  if (!SWIG_IsOK(errcode))
//    SWIG_exception_fail(SWIG_ArgError(errcode), "in method $symname, wrong return value: $type");
//}

//%typemap(out, fragment="FromCpp") MatrixDense* MatrixDense::create
//{
//  int errcode = matrixDenseFromCppCreate(&($result), *$1);
//  if (!SWIG_IsOK(errcode))
//    SWIG_exception_fail(SWIG_ArgError(errcode), "in method $symname, wrong return value: $type");
//}

//%typemap(out, fragment="FromCpp") MatrixDense& MatrixDense::create
//{
//  int errcode = matrixDenseFromCppCreate(&($result), *$1);
//  if (!SWIG_IsOK(errcode))
//    SWIG_exception_fail(SWIG_ArgError(errcode), "in method $symname, wrong return value: $type");
//}

//%typemap(out, fragment="FromCpp") MatrixDense*,     MatrixDense&,
//                                  MatrixSquare*,   MatrixSquare&,
//                                  MatrixSymmetric*, MatrixSymmetric&
//{
//  int errcode = matrixDenseFromCpp(&($result), *$1);
//  if (!SWIG_IsOK(errcode))
//    SWIG_exception_fail(SWIG_ArgError(errcode), "in method $symname, wrong return value: $type");
//}

//%typemap(out, fragment="FromCpp") MatrixSparse
//{
//  int errcode = matrixSparseFromCpp(&($result), $1);
//  if (!SWIG_IsOK(errcode))
//    SWIG_exception_fail(SWIG_ArgError(errcode), "in method $symname, wrong return value: $type");
//}

//%typemap(out, fragment="FromCpp") MatrixSparse*,     MatrixSparse&
//{
//  int errcode = matrixSparseFromCpp(&($result), *$1);
//  if (!SWIG_IsOK(errcode))
//    SWIG_exception_fail(SWIG_ArgError(errcode), "in method $symname, wrong return value: $type");
//}

//%typemap(out, fragment="FromCpp") MatrixSparse* MatrixSparse::create
//{
//  int errcode = matrixSparseFromCppCreate(&($result), *$1);
//  if (!SWIG_IsOK(errcode))
//    SWIG_exception_fail(SWIG_ArgError(errcode), "in method $symname, wrong return value: $type");
//}

//%typemap(out, fragment="FromCpp") MatrixSparse& MatrixSparse::create
//{
//  int errcode = matrixSparseFromCppCreate(&($result), *$1);
//  if (!SWIG_IsOK(errcode))
//    SWIG_exception_fail(SWIG_ArgError(errcode), "in method $symname, wrong return value: $type");
//}

} // namespace gstlrn

%extend gstlrn::Grid {
  double indiceToCoordinate(int idim0, const VectorInt& indice,
                            const VectorDouble& percent = {},
                            bool flag_rotate            = true) const
  {
    return $self->indiceToCoordinate(idim0, indice.getVector(), percent.getVector(), flag_rotate);
  }

  Id indiceToRank(const VectorInt &indice) const
  {
    return $self->indiceToRank(indice.getVector());
  }

  void rankToIndice(Id rank, VectorInt &indices, bool minusOne = false) const
  {
    return $self->rankToIndice(rank, indices.getVector(), minusOne);
  }

  void indicesToCoordinateInPlace(const VectorInt& indice,
                                  VectorDouble& coor,
                                  const VectorDouble& percent = VectorDouble(),
                                  bool flag_rotate=true) const
  {
    return $self->indicesToCoordinateInPlace(indice, coor, percent, flag_rotate);
  }
};

%extend gstlrn::Rule {
  // Don't return std::array to wrapping languages
  VectorDouble getThresh(Id facies) const
  {
    const auto thresh = $self->getThresh(facies);
    return VectorDouble{thresh.begin(), thresh.end()};
  }
};

%extend gstlrn::KNN {
  VectorInt getIndices(int rank = 0) const {
    const auto view = $self->getIndices(rank);
    return {view.begin(), view.end()};
  }
}

%extend gstlrn::MatrixDense {
/**
  * @brief List of methods from class MatrixDense exported for Target Language
  */
  MatrixDense linearCombination(double addition,
                         double val1,
                         const MatrixDense& other1,
                         double val2     = 0.,
                         const MatrixDense& other2 = MatrixDense(),
                         double val3     = 0.,
                         const MatrixDense& other3 = MatrixDense())
  {
    return AMatrix::linearCombination(addition, val1, other1, val2, other2, val3, other3);
  }
  VectorDouble prodMatVec(const VectorDouble& x, bool transpose = false)
  {
    return AMatrix::product(*$self, x, transpose, true);
  }
  void prodMatVecInPlace(const VectorDouble& x, VectorDouble& y, bool transpose = false) const
  {
    AMatrix::productInPlace(y, *$self, x, transpose, true);
  }
  VectorDouble prodVecMat(const VectorDouble& x, bool transpose = false) const
  {
    return AMatrix::product(x, *$self, transpose, true);
  }
  void prodVecMatInPlace(const VectorDouble& x, VectorDouble& y, bool transpose = false) const
  {
    AMatrix::productInPlace(y, x, *$self, transpose, true);
  }
  void AMatrix::prodMat(const MatrixDense* matY, bool transposeY)
  {
    AMatrix::prodMatMatInPlace(*self, *self, *matY, false, transposeY);
  }
  void AMatrix::prodMatMatInPlace(const MatrixDense* x,
                                  const MatrixDense* y,
                                  bool transposeX = false,
                                  bool transposeY = false)
  {
    AMatrix::prodMatMatInPlace(*self, *x, *y, transposeX, transposeY);
  }
}

%extend gstlrn::MatrixSparse {
/**
  * @brief List of methods from class MatrixSparse exported for Target Language
  */
  MatrixSparse linearCombination(double addition,
                         double val1,
                         const MatrixSparse& other1,
                         double val2     = 0.,
                         const MatrixSparse& other2 = MatrixSparse(),
                         double val3     = 0.,
                         const MatrixSparse& other3 = MatrixSparse())
  {
    return AMatrix::linearCombination(addition, val1, other1, val2, other2, val3, other3);
  }
  VectorDouble prodMatVec(const VectorDouble& x, bool transpose = false)
  {
    return AMatrix::product(*$self, x, transpose, true);
  }
  void prodMatVecInPlace(const VectorDouble& x, VectorDouble& y, bool transpose = false) const
  {
    AMatrix::productInPlace(y, *$self, x, transpose, true);
  }
  VectorDouble prodVecMat(const VectorDouble& x, bool transpose = false) const
  {
     return AMatrix::product(x, *$self, transpose, true);
  }
  void prodVecMatInPlace(const VectorDouble& x, VectorDouble& y, bool transpose = false) const
  {
    AMatrix::productInPlace(y, x, *$self, transpose, true);
  }
  void AMatrix::prodMat(const MatrixSparse* matY, bool transposeY)
  {
    AMatrix::prodMatMatInPlace(*self, *self, *matY, false, transposeY);
  }
  void AMatrix::prodMatMatInPlace(const MatrixSparse* x,
                                  const MatrixSparse* y,
                                  bool transposeX = false,
                                  bool transposeY = false)
  {
    AMatrix::prodMatMatInPlace(*self, *x, *y, transposeX, transposeY);
  }
}

%extend gstlrn::VectorHelper {
/**
  * @brief List of operators from class VectorHelper exported for Target Language
  */
  static VectorDouble addVD(const VectorDouble& v1, const VectorDouble& v2)
  {
    return VectorHelper::add(v1, v2);
  }
  static VectorDouble addVDCst(const VectorDouble& v1, double v2)
  {
    return VectorHelper::addCst(v1, v2);
  }
  static VectorDouble subtractVD(const VectorDouble& v1, const VectorDouble& v2)
  {
    return VectorHelper::subtract(v1, v2);
  }
  static VectorDouble subtractVDCst(const VectorDouble& v1, double v2, bool flagOpposite = false)
  {
    return VectorHelper::subtractCst(v1, v2, flagOpposite);
  }
  static VectorDouble multiplyVD(const VectorDouble& v1, const VectorDouble& v2)
  {
    return VectorHelper::multiply(v1, v2);
  }
  static VectorDouble multiplyVDCst(const VectorDouble& v1, double v2)
  {
    return VectorHelper::multiplyCst(v1, v2);
  }
  static VectorDouble divideVD(const VectorDouble& v1, const VectorDouble& v2)
  {
    return VectorHelper::divide(v1, v2);
  }
  static VectorDouble divideVDCst(const VectorDouble& v1, double v2, bool flagOpposite = false)
  {
    return VectorHelper::divideCst(v1, v2, flagOpposite);
  }
  static VectorInt addVI(const VectorInt& v1, const VectorInt& v2)
  {
    return VectorHelper::add(v1, v2);
  }
  static VectorInt addVICst(const VectorInt& v1, Id v2)
  {
    return VectorHelper::addCst(v1, v2);
  }
  static VectorInt subtractVI(const VectorInt& v1, const VectorInt& v2)
  {
    return VectorHelper::subtract(v1, v2);
  }
  static VectorInt subtractVICst(const VectorInt& v1, Id v2, bool flagOpposite = false)
  {
    return VectorHelper::subtractCst(v1, v2, flagOpposite);
  }
  static VectorInt multiplyVI(const VectorInt& v1, const VectorInt& v2)
  {
    return VectorHelper::multiply(v1, v2);
  }
  static VectorInt multiplyVICst(const VectorInt& v1, Id v2)
  {
    return VectorHelper::multiplyCst(v1, v2);
  }
  static VectorInt divideVI(const VectorInt& v1, const VectorInt& v2)
  {
    return VectorHelper::divide(v1, v2);
  }
  static VectorInt divideVICst(const VectorInt& v1, Id v2, bool flagOpposite = false)
  {
    return VectorHelper::divideCst(v1, v2, flagOpposite);
  }

}

%extend gstlrn::DbData {

//----> In DbData: getValue()

double getValueD(ColID&& colid, Id isample, double def = -1.)
{
  return $self->getValue<double>(std::move(colid), isample).value_or(def);
}

float getValueF(ColID&& colid, Id isample, float def = -1.)
{
  return $self->getValue<float>(std::move(colid), isample).value_or(def);
}

Id getValueI(ColID&& colid, Id isample, Id def = -1)
{
  return $self->getValue<Id>(std::move(colid), isample).value_or(def);
}

UChar getValueU(ColID&& colid, Id isample, UChar def = 0)
{
  return $self->getValue<UChar>(std::move(colid), isample).value_or(def);
}

String getValueS(ColID&& colid, Id isample)
{
  return $self->getValue<String>(
    std::move(colid), isample).value_or(String("failed"));
}

bool getValueB(ColID&& colid, Id isample, bool def = false)
{
  return $self->getValue<bool>(std::move(colid), isample).value_or(def);
}

Category getValueC(ColID&& colid, Id isample,
          const Category& def = Category())
{
    return $self->getValue<Category>(std::move(colid), isample).value_or(def);
}

//----> In DbData: setValue()

void setValueD(ColID&& colid, Id isample, double value)
{
  $self->setValue<double>(std::move(colid), isample, value);
}

void setValueF(ColID&& colid, Id isample, float value)
{
  $self->setValue<float>(std::move(colid), isample, value);
}

void setValueI(ColID&& colid, Id isample, Id value)
{
  $self->setValue<Id>(std::move(colid), isample, value);
}

void setValueU(ColID&& colid, Id isample, UChar value)
{
  $self->setValue<UChar>(std::move(colid), isample, value);
}

void setValueS(ColID&& colid, Id isample, const String& value)
{
  $self->setValue<String>(std::move(colid), isample, value);
}

void setValueB(ColID&& colid, Id isample, bool value)
{
  $self->setValue<bool>(std::move(colid), isample, value);
}

void setValueC(ColID&& colid, Id isample, const Category& value)
{
    $self->setValue<Category>(std::move(colid), isample, value);
}

//----> In DbData: getColumn()

VectorDouble getColumnD(ColID&& colid)
{
  return $self->getColumn<VectorDouble>(std::move(colid));
}

VectorFloat getColumnF(ColID&& colid)
{
  return $self->getColumn<VectorFloat>(std::move(colid));
}

VectorInt getColumnI(ColID&& colid)
{
  return $self->getColumn<VectorInt>(std::move(colid));
}

VectorUChar getColumnU(ColID&& colid)
{
  return $self->getColumn<VectorUChar>(std::move(colid));
}

VectorString getColumnS(ColID&& colid)
{
  return $self->getColumn<VectorString>(std::move(colid));
}

VectorBool getColumnB(ColID&& colid)
{
  return $self->getColumn<VectorBool>(std::move(colid));
}

VectorCategory getColumnC(ColID&& colid)
{
  return $self->getColumn<VectorCategory>(std::move(colid));
}

}

// Prevent memory leaks from 'create*' and 'clone' methods

// The following file should contain all 'createFrom*' methods
%include swig/newobject.i
// This is for all 'create' methods
%newobject *::create;
// So bad that the following syntax doesn't work:
// %newobject *::create*;
// This is for all 'clone' methods
%newobject *::clone;

%{
  #include <memory>
%}

//quick and dirty way to handle shared_ptrs for one python test
// This is not a good practice, but it works for now
// better solution would be to use a shared_ptr
%typemap(in) const std::shared_ptr<const gstlrn::ASimulable> & {
  gstlrn::ASimulable* ptr = nullptr;
  int res = SWIG_ConvertPtr($input, (void**)&ptr, SWIGTYPE_p_gstlrn__ASimulable, 0);

  if (SWIG_IsOK(res) && ptr != nullptr) {
    *$1 = std::shared_ptr<const gstlrn::ASimulable>(ptr);
  } else {
    SWIG_exception_fail(SWIG_TypeError, "Expected ASimulable-derived object");
  }
}

%typemap(in) const std::shared_ptr<const gstlrn::MatrixSparse> & {
  gstlrn::MatrixSparse* ptr = nullptr;
  int res = SWIG_ConvertPtr($input, (void**)&ptr, SWIGTYPE_p_gstlrn__MatrixSparse, 0);

  if (SWIG_IsOK(res) && ptr != nullptr) {
    *$1 = std::shared_ptr<const gstlrn::MatrixSparse>(ptr);
  } else {
    SWIG_exception_fail(SWIG_TypeError, "Expected MatrixSparse object");
  }
}

%include <std_shared_ptr.i>
%template(ASpaceSharedPtr)    std::shared_ptr<const gstlrn::ASpace>;
%template(ASpaceSharedPtrVector)   std::vector< gstlrn::ASpaceSharedPtr>;
%template(MatrixSparseSh) std::shared_ptr<const gstlrn::MatrixSparse >;
