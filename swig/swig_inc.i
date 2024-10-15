%feature(director) IProjMatrix;
%feature(director) AFunctional;

//%feature(director) ICLoneable;
%feature(director) ABiTargetCheck;

%{
  #include "gstlearn_export.hpp"
  #include "geoslib_define.h"
  #include "geoslib_enum.h"
  #include "geoslib_d.h"
  #include "geoslib_f.h"
  #include "geoslib_old_f.h"
  
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
  #include "Basic/MathFunc.hpp"
  #include "Basic/Indirection.hpp"
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
  #include "Mesh/MeshEStandard.hpp"
  #include "Mesh/MeshETurbo.hpp"
  #include "Mesh/MeshSpherical.hpp"
  #include "Mesh/MeshSphericalExt.hpp"
  
  #include "Polynomials/APolynomial.hpp"
  #include "Polynomials/ClassicalPolynomial.hpp"
  #include "Polynomials/Hermite.hpp"
  #include "Polynomials/MonteCarlo.hpp"
  
  #include "LinearOp/CGParam.hpp"
  #include "LinearOp/LogStats.hpp"
  #include "LinearOp/ALinearOp.hpp"
  #include "LinearOp/ASimulable.hpp"
  #include "LinearOp/LinearOpCGSolver.hpp"
  #include "LinearOp/ALinearOpMulti.hpp"
  #include "LinearOp/ShiftOpCs.hpp"
  #include "LinearOp/PrecisionOp.hpp"
  #include "LinearOp/PrecisionOpCs.hpp"
  #include "LinearOp/TurboOptimizer.hpp"
  #include "LinearOp/IProjMatrix.hpp"
  #include "LinearOp/ScaleOp.hpp"
  #include "LinearOp/ProjMatrix.hpp"
  #include "LinearOp/ProjMulti.hpp"
  #include "LinearOp/ProjMultiMatrix.hpp"
  #include "LinearOp/PrecisionOpMulti.hpp"
  #include "LinearOp/PrecisionOpMultiMatrix.hpp"
  #include "LinearOp/PrecisionOpMultiConditional.hpp"
  #include "LinearOp/IOptimCost.hpp"
  #include "LinearOp/OptimCostBinary.hpp"
  #include "LinearOp/OptimCostColored.hpp"
  #include "LinearOp/ProjConvolution.hpp"
  #include "LinearOp/SPDEOp.hpp"
  #include "LinearOp/SPDEOpMatrix.hpp"
  #include "LinearOp/MatrixSquareSymmetricSim.hpp"

  #include "Neigh/ANeigh.hpp"
  #include "Neigh/NeighUnique.hpp"
  #include "Neigh/NeighImage.hpp"
  #include "Neigh/NeighMoving.hpp"
  #include "Neigh/NeighBench.hpp"
  #include "Neigh/NeighCell.hpp"
  
  #include "Variogram/AVario.hpp"
  #include "Variogram/VarioParam.hpp"
  #include "Variogram/Vario.hpp"
  #include "Variogram/DirParam.hpp"
  #include "Variogram/VMap.hpp"
  #include "Variogram/VCloud.hpp"
  
  #include "Model/Model.hpp"
  #include "Model/Option_AutoFit.hpp"
  #include "Model/Option_VarioFit.hpp"
  #include "Model/Constraints.hpp"
  #include "Model/ConsItem.hpp"
  #include "Model/CovParamId.hpp"
  #include "Model/CovParamId.hpp"
  
  #include "Covariances/ParamId.hpp"
  #include "Covariances/TabNoStat.hpp"
  #include "Covariances/TabNoStatCovAniso.hpp"
  #include "Covariances/ANoStat.hpp"
  #include "Covariances/NoStatArray.hpp"
  #include "Covariances/NoStatFunctional.hpp"
  #include "Covariances/ACov.hpp"
  #include "Covariances/ACovFunc.hpp"
  #include "Covariances/ACovAnisoList.hpp"
  #include "Covariances/CovAniso.hpp"
  #include "Covariances/ACovGradient.hpp"
  #include "Covariances/CovGneiting.hpp"
  #include "Covariances/CovLMCTapering.hpp"
  #include "Covariances/CovLMCConvolution.hpp"
  #include "Covariances/CovLMCAnamorphosis.hpp"
  #include "Covariances/CovLMGradient.hpp"
  #include "Covariances/CovContext.hpp"
  #include "Covariances/CovCalcMode.hpp"
  #include "Covariances/CovBesselJ.hpp"
  #include "Covariances/CovMatern.hpp"
  #include "Covariances/CovCauchy.hpp"
  #include "Covariances/CovCosExp.hpp"
  #include "Covariances/CovCosinus.hpp"
  #include "Covariances/CovCubic.hpp"
  #include "Covariances/CovExponential.hpp"
  #include "Covariances/CovGamma.hpp"
  #include "Covariances/CovGaussian.hpp"
  #include "Covariances/CovGC1.hpp"
  #include "Covariances/CovGC3.hpp"
  #include "Covariances/CovGC5.hpp"
  #include "Covariances/CovGCspline2.hpp"
  #include "Covariances/CovGCspline.hpp"
  #include "Covariances/CovLinear.hpp"
  #include "Covariances/CovNugget.hpp"
  #include "Covariances/CovMarkov.hpp"
  #include "Covariances/CovPenta.hpp"
  #include "Covariances/CovPower.hpp"
  #include "Covariances/CovReg1D.hpp"
  #include "Covariances/CovSincard.hpp"
  #include "Covariances/CovSpherical.hpp"
  #include "Covariances/CovStable.hpp"
  #include "Covariances/CovStorkey.hpp"
  #include "Covariances/CovTriangle.hpp"
  #include "Covariances/CovWendland0.hpp"
  #include "Covariances/CovWendland1.hpp"
  #include "Covariances/CovWendland2.hpp"
  #include "Covariances/CovGeometric.hpp"
  #include "Covariances/CovPoisson.hpp"
  #include "Covariances/CovLinearSph.hpp"
  #include "Covariances/CovDiffusionAdvection.hpp"
  #include "Covariances/CovHelper.hpp"
  
  #include "Drifts/ADrift.hpp"
  #include "Drifts/DriftList.hpp"
  #include "Drifts/DriftM.hpp"
  #include "Drifts/DriftF.hpp"
  #include "Drifts/DriftFactory.hpp"
  
  #include "Matrix/AMatrix.hpp"
  #include "Matrix/AMatrixDense.hpp"
  #include "Matrix/MatrixSparse.hpp"
  #include "Matrix/LinkMatrixSparse.hpp"
  #include "Matrix/AMatrixSquare.hpp"
  #include "Matrix/NF_Triplet.hpp"
  #include "Matrix/MatrixRectangular.hpp"
  #include "Matrix/MatrixSquareGeneral.hpp"
  #include "Matrix/MatrixSquareSymmetric.hpp"
  #include "Matrix/MatrixFactory.hpp"
  #include "Matrix/MatrixInt.hpp"
  #include "Matrix/Table.hpp"
  
  #include "API/SPDE.hpp"
  #include "API/PGSSPDE.hpp"
  #include "API/TestInheritance.hpp"
  #include "API/Style.hpp"
  #include "API/SPDEParam.hpp"
  
  #include "Db/Db.hpp"
  #include "Db/DbGrid.hpp"
  #include "Db/DbLine.hpp"
  #include "Db/DbGraphO.hpp"
  #include "Db/DbMeshTurbo.hpp"
  #include "Db/DbMeshStandard.hpp"
  #include "Db/DbStringFormat.hpp"
  #include "Db/DbHelper.hpp"
  
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
  
  #include "Morpho/Morpho.hpp"
  
  #include "Polygon/Polygons.hpp"
  #include "Polygon/PolyElem.hpp"
  
  #include "Stats/Classical.hpp"
  #include "Stats/PCA.hpp"
  #include "Stats/PCAStringFormat.hpp"
  #include "Stats/Selectivity.hpp"
  #include "Stats/Regression.hpp"
  
  #include "LithoRule/Rule.hpp"
  #include "LithoRule/RuleStringFormat.hpp"
  #include "LithoRule/RuleProp.hpp"
  
  #include "Estimation/KrigingSystem.hpp"
  #include "Estimation/KrigingCalcul.hpp"
  #include "Estimation/CalcKriging.hpp"
  #include "Estimation/CalcKrigingFactors.hpp"
  #include "Estimation/CalcSimpleInterpolation.hpp"
  #include "Estimation/CalcImage.hpp"
  #include "Estimation/CalcGlobal.hpp"
  
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
  #include "Simulation/ACalcSimulation.hpp"
  #include "Simulation/CalcSimuTurningBands.hpp"
  #include "Simulation/TurningBandDirection.hpp"
  #include "Simulation/TurningBandOperate.hpp"
  #include "Simulation/SimuSpectral.hpp"
  #include "Simulation/BooleanObject.hpp"
  #include "Simulation/SimuBoolean.hpp"
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
  #include "Tree/KNN.hpp"
    
  #include "Spatial/Projection.hpp"
  #include "Spatial/SpatialIndices.hpp"

  #include "Core/Potential.hpp"
      
  // Mask some warning generated by SWIG:
  //DISABLE_WARNING_DECLARATION_MASKED
  //DISABLE_WARNING_EXPR_COND_ASSIGNMENT
  //DISABLE_WARNING_UNREFERENCED_FORMAL_PARAMETER
%}

////////////////////////////
//        Typemaps        //
////////////////////////////

// Mandatory for using swig::asptr and swig::from for std::vectors
%include std_vector.i
%include std_string.i
%template(DoNotUseVectorIntStd)     std::vector< int >;
%template(DoNotUseVectorDoubleStd)  std::vector< double >;
%template(DoNotUseVectorStringStd)  std::vector< std::string >; // Keep std::string here otherwise asptr fails!
%template(DoNotUseVectorFloatStd)   std::vector< float >;
%template(DoNotUseVectorUCharStd)   std::vector< unsigned char >; // Keep unsigned char here
%template(DoNotUseVectorBoolStd)    std::vector< bool >;
%template(DoNotUseVVectorIntStd)    std::vector< std::vector< int > >;
%template(DoNotUseVVectorDoubleStd) std::vector< std::vector< double > >;
%template(DoNotUseVVectorFloatStd)  std::vector< std::vector< float > >; 

%template(VectorECov)              std::vector< ECov >;
%template(VectorEStatOption)       std::vector< EStatOption >;
%template(VectorESelectivity)      std::vector< ESelectivity >;
%template(VectorDirParam)          std::vector< DirParam >;
%template(VectorPolyElem)          std::vector< PolyElem >;
%template(VectorInterval)          std::vector< Interval >; 
%template(VectorEPostStat)         std::vector< EPostStat >;
%template(VectorSpacePoint)        std::vector< SpacePoint >;
%template(VectorABiTargetCheck)    std::vector< ABiTargetCheck* >;
%template(VectorProjMatrix)        std::vector< ProjMatrix* >;
%template(VectorConstProjMatrix)   std::vector< const ProjMatrix*>;
%template(VectorConstIProjMatrix)  std::vector< const IProjMatrix*>;
%template(VVectorConstProjMatrix)  std::vector< std::vector< const ProjMatrix*> >;
%template(VVectorConstIProjMatrix) std::vector< std::vector< const IProjMatrix*> >;
%template(VectorMeshes)            std::vector< const AMesh*>;
////////////////////////////////////////////////
// Conversion Target language => C++

// Note : Before including this file :
//        - vectorToCpp, vectorVectorToCpp, matrixDenseToCpp, matrixSparseToCpp and convertToCpp 
//          functions must be defined in ToCpp fragment

// Convert scalar arguments by value
%typemap(in, fragment="ToCpp") int,
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
%typemap(in, fragment="ToCpp") std::string_view
{
  try
  {
    static String tmp;
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
%typemap(in, fragment="ToCpp") int*       (int val), const int*       (int val),
                               int&       (int val), const int&       (int val),
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

%typemap(in, fragment="ToCpp") const VectorInt&    (void *argp, VectorInt vec),
                               const VectorInt*    (void *argp, VectorInt vec),
                               const VectorDouble& (void *argp, VectorDouble vec),
                               const VectorDouble* (void *argp, VectorDouble vec),
                               const VectorString& (void *argp, VectorString vec),
                               const VectorString* (void *argp, VectorString vec),
                               const VectorFloat&  (void *argp, VectorFloat vec),
                               const VectorFloat*  (void *argp, VectorFloat vec),
                               const VectorUChar&  (void *argp, VectorUChar vec),
                               const VectorUChar*  (void *argp, VectorUChar vec),
                               const VectorBool&   (void *argp, VectorBool vec),
                               const VectorBool*   (void *argp, VectorBool vec)
{
  // Try to convert from any target language vector
  int errcode = vectorToCpp($input, vec);
  if (!SWIG_IsOK(errcode))
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
                               const VectorVectorInt*    (void *argp, VectorVectorInt vec),
                               const VectorVectorDouble& (void *argp, VectorVectorDouble vec),
                               const VectorVectorDouble* (void *argp, VectorVectorDouble vec),
                               const VectorVectorFloat&  (void *argp, VectorVectorFloat vec),
                               const VectorVectorFloat*  (void *argp, VectorVectorFloat vec)
{
  // Try to convert from any target language vector
  int errcode = vectorVectorToCpp($input, vec);
  if (!SWIG_IsOK(errcode))
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

%typemap(in, fragment="ToCpp") const MatrixRectangular&     (void *argp, MatrixRectangular mat),
                               const MatrixRectangular*     (void *argp, MatrixRectangular mat),
                               const MatrixSquareGeneral&   (void *argp, MatrixSquareGeneral mat),
                               const MatrixSquareGeneral*   (void *argp, MatrixSquareGeneral mat),
                               const MatrixSquareSymmetric& (void *argp, MatrixSquareSymmetric mat),
                               const MatrixSquareSymmetric* (void *argp, MatrixSquareSymmetric mat)
{
  // Try to convert from any target language vector
  int errcode = matrixDenseToCpp($input, mat);
  if (!SWIG_IsOK(errcode))
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
  if (!SWIG_IsOK(errcode))
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

%typemap(out, fragment="FromCpp") int,
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

%typemap(out, fragment="FromCpp") int*,    const int*,    int&,    const int&,
                                  double*, const double*, double&, const double&,
                                  String*, const String*, String&, const String&,
                                  std::string_view*, const std::string_view*,
                                  std::string_view&, const std::string_view&,
                                  float*,  const float*,  float&,  const float&,
                                  UChar*,  const UChar*,  UChar&,  const UChar&,
                                  bool*,   const bool*,   bool&,   const bool&
{
  $result = objectFromCpp(*$1);
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

//%typemap(out, fragment="FromCpp") MatrixRectangular, 
//                                  MatrixSquareGeneral, 
//                                  MatrixSquareSymmetric
//{
//  int errcode = matrixDenseFromCpp(&($result), $1);
//  if (!SWIG_IsOK(errcode))
//    SWIG_exception_fail(SWIG_ArgError(errcode), "in method $symname, wrong return value: $type");
//}

//%typemap(out, fragment="FromCpp") MatrixRectangular* MatrixRectangular::create
//{
//  int errcode = matrixDenseFromCppCreate(&($result), *$1);
//  if (!SWIG_IsOK(errcode))
//    SWIG_exception_fail(SWIG_ArgError(errcode), "in method $symname, wrong return value: $type");
//}

//%typemap(out, fragment="FromCpp") MatrixRectangular& MatrixRectangular::create
//{
//  int errcode = matrixDenseFromCppCreate(&($result), *$1);
//  if (!SWIG_IsOK(errcode))
//    SWIG_exception_fail(SWIG_ArgError(errcode), "in method $symname, wrong return value: $type");
//}

//%typemap(out, fragment="FromCpp") MatrixRectangular*,     MatrixRectangular&,
//                                  MatrixSquareGeneral*,   MatrixSquareGeneral&,
//                                  MatrixSquareSymmetric*, MatrixSquareSymmetric&
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

