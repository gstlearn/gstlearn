
%{
  #include "gstlearn_export.hpp"
  #include "geoslib_define.h"
  #include "geoslib_enum.h"
  #include "geoslib_d.h"
  #include "geoslib_f.h"
  #include "geoslib_old_f.h"
  
  #include "csparse_d.h"
  #include "csparse_f.h"
  
  #include "Basic/VectorT.hpp"
  #include "Basic/VectorNumT.hpp"
  #include "Basic/ICloneable.hpp"
  #include "Basic/Vector.hpp"
  #include "Basic/AFunction.hpp"
  #include "Basic/ArgumentTest.hpp"
  #include "Basic/AStringable.hpp"
  #include "Basic/AStringFormat.hpp"
  #include "Basic/ASerializable.hpp"
  #include "Basic/Rotation.hpp"
  #include "Basic/Tensor.hpp"
  #include "Basic/Grid.hpp"
  #include "Basic/String.hpp"
  #include "Basic/Interval.hpp"
  #include "Basic/Limits.hpp"
  #include "Basic/Utilities.hpp"
  #include "Basic/CSVformat.hpp"
  #include "Basic/AFunctional.hpp"
  #include "Basic/FunctionalSpirale.hpp"
  #include "Basic/Law.hpp"
  #include "Basic/Table.hpp"
  #include "Basic/Utilities.hpp"
  #include "Basic/NamingConvention.hpp"
  #include "Basic/OptDbg.hpp"
  #include "Basic/OptCst.hpp"
  #include "Basic/OptCustom.hpp"
  #include "Basic/EDbg.hpp"
  #include "Basic/ECst.hpp"
  #include "Basic/File.hpp"
  #include "Basic/Limits.hpp"
  #include "Basic/Plane.hpp"
  #include "Basic/FFT.hpp"
  #include "Basic/Array.hpp"
  #include "Basic/Geometry.hpp"
  #include "Basic/PolyLine2D.hpp"
  
  #include "Faults/Faults.hpp"
  
  #include "Boolean/ETLaw.hpp"
  #include "Boolean/ETShape.hpp"
  #include "Boolean/ShapeParameter.hpp"
  #include "Boolean/AShape.hpp"
  #include "Boolean/ShapeParallelepiped.hpp"
  #include "Boolean/ShapeEllipsoid.hpp"
  #include "Boolean/ShapeParaboloid.hpp"
  #include "Boolean/ShapeHalfEllipsoid.hpp"
  #include "Boolean/ShapeHalfParaboloid.hpp"
  #include "Boolean/ShapeHalfSinusoid.hpp"
  #include "Boolean/ModelBoolean.hpp"
  
  #include "Enum/EKrigOpt.hpp"
  
  #include "Space/Space.hpp"
  #include "Space/ASpace.hpp"
  #include "Space/ASpaceObject.hpp"
  #include "Space/SpacePoint.hpp"
  #include "Space/SpaceRN.hpp"
  #include "Space/SpaceShape.hpp"
  
  #include "Skin/ISkinFunctions.hpp"
  #include "Skin/Skin.hpp"
  
  #include "Calculators/ACalculator.hpp"
  #include "Calculators/ACalcInterpolator.hpp"
  
  /*
  #include "Interfaces/geoslib_f_swig.h"
  #include "Interfaces/AParam.hpp"
  #include "Interfaces/AVariable.hpp"
  #include "Interfaces/AVariableTemplate.hpp"
  #include "Interfaces/Category.hpp"
  #include "Interfaces/Database.hpp"
  #include "Interfaces/Dictionary.hpp"
  #include "Interfaces/interface_d.hpp"
  #include "Interfaces/ParamCSV.hpp"
  #include "Interfaces/ParamGrid.hpp"
  #include "Interfaces/Param.hpp"
  #include "Interfaces/VariableBool.hpp"
  #include "Interfaces/VariableCategorical.hpp"
  #include "Interfaces/VariableDouble.hpp"
  #include "Interfaces/VariableInt.hpp"
  #include "Interfaces/VariableString.hpp"
  */
  
  #include "Mesh/AMesh.hpp"
  #include "Mesh/MeshFactory.hpp"
  #include "Mesh/MeshEStandard.hpp"
  #include "Mesh/MeshETurbo.hpp"
  #include "Mesh/MeshSpherical.hpp"
  
  #include "Enum/AEnum.hpp"
  
  #include "Polynomials/APolynomial.hpp"
  #include "Polynomials/ClassicalPolynomial.hpp"
  #include "Polynomials/Hermite.hpp"
  #include "Polynomials/MonteCarlo.hpp"
  
  #include "LinearOp/ALinearOp.hpp"
  #include "LinearOp/ALinearOpMulti.hpp"
  #include "LinearOp/ShiftOpCs.hpp"
  #include "LinearOp/PrecisionOp.hpp"
  #include "LinearOp/PrecisionOpCs.hpp"
  #include "LinearOp/TurboOptimizer.hpp"
  #include "LinearOp/IProjMatrix.hpp"
  #include "LinearOp/ProjMatrix.hpp"
  #include "LinearOp/PrecisionOpMultiConditional.hpp"
  #include "LinearOp/IOptimCost.hpp"
  #include "LinearOp/OptimCostBinary.hpp"
  #include "LinearOp/OptimCostColored.hpp"
  #include "LinearOp/EPowerPT.hpp"
  
  #include "Neigh/ANeighParam.hpp"
  #include "Neigh/NeighUnique.hpp"
  #include "Neigh/NeighImage.hpp"
  #include "Neigh/NeighMoving.hpp"
  #include "Neigh/NeighBench.hpp"
  #include "Neigh/ENeigh.hpp"
  #include "Neigh/NeighWork.hpp"
  
  #include "Variogram/VarioParam.hpp"
  #include "Variogram/Vario.hpp"
  #include "Variogram/DirParam.hpp"
  #include "Variogram/ECalcVario.hpp"
  
  #include "Model/Model.hpp"
  #include "Model/ANoStat.hpp"
  #include "Model/NoStatArray.hpp"
  #include "Model/NoStatFunctional.hpp"
  #include "Model/Option_AutoFit.hpp"
  #include "Model/Option_VarioFit.hpp"
  #include "Model/Constraints.hpp"
  #include "Model/ConsItem.hpp"
  #include "Model/CovParamId.hpp"
  #include "Model/EModelProperty.hpp"
  #include "Model/EConsElem.hpp"
  #include "Model/EConsType.hpp"
  #include "Model/CovParamId.hpp"
  
  #include "Covariances/ACov.hpp"
  #include "Covariances/ACovFunc.hpp"
  #include "Covariances/ACovAnisoList.hpp"
  #include "Covariances/CovAniso.hpp"
  #include "Covariances/ACovGradient.hpp"
  #include "Covariances/CovLMC.hpp"
  #include "Covariances/CovLMCTapering.hpp"
  #include "Covariances/CovLMCConvolution.hpp"
  #include "Covariances/CovLMCAnamorphosis.hpp"
  #include "Covariances/CovLMGradient.hpp"
  #include "Covariances/CovContext.hpp"
  #include "Covariances/CovCalcMode.hpp"
  #include "Covariances/CovBesselJ.hpp"
  #include "Covariances/CovBesselK.hpp"
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
  #include "Covariances/CovDiffusionAdvection.hpp"
  #include "Covariances/ECov.hpp"
  #include "Covariances/ETape.hpp"
  #include "Covariances/EConvType.hpp"
  #include "Covariances/EConvDir.hpp"
  
  #include "Drifts/ADrift.hpp"
  #include "Drifts/ADriftElem.hpp"
  #include "Drifts/DriftList.hpp"
  #include "Drifts/Drift1.hpp"
  #include "Drifts/DriftF.hpp"
  #include "Drifts/DriftFactory.hpp"
  #include "Drifts/DriftX.hpp"
  #include "Drifts/DriftX2.hpp"
  #include "Drifts/DriftX2Y.hpp"
  #include "Drifts/DriftX3.hpp"
  #include "Drifts/DriftXY.hpp"
  #include "Drifts/DriftXY2.hpp"
  #include "Drifts/DriftXZ.hpp"
  #include "Drifts/DriftY.hpp"
  #include "Drifts/DriftY2.hpp"
  #include "Drifts/DriftY3.hpp"
  #include "Drifts/DriftYZ.hpp"
  #include "Drifts/DriftZ.hpp"
  #include "Drifts/DriftZ2.hpp"
  #include "Drifts/DriftZ3.hpp"
  #include "Drifts/EDrift.hpp"
  
  #include "Matrix/AMatrix.hpp"
  #include "Matrix/AMatrixSquare.hpp"
  #include "Matrix/MatrixRectangular.hpp"
  #include "Matrix/MatrixSquareDiagonal.hpp"
  #include "Matrix/MatrixSquareDiagonalCst.hpp"
  #include "Matrix/MatrixSquareGeneral.hpp"
  #include "Matrix/MatrixSquareSymmetric.hpp"
  #include "Matrix/MatrixInt.hpp"
  
  #include "API/SPDE.hpp"
  #include "API/PGSSPDE.hpp"
  #include "API/ESPDECalcMode.hpp"
  
  #include "Db/Db.hpp"
  #include "Db/DbGrid.hpp"
  #include "Db/DbStringFormat.hpp"
  #include "Db/ELoadBy.hpp"
  #include "Db/ELoc.hpp"
  
  #include "Anamorphosis/AAnam.hpp"
  #include "Anamorphosis/AnamContinuous.hpp"
  #include "Anamorphosis/AnamDiscrete.hpp"
  #include "Anamorphosis/AnamUser.hpp"
  #include "Anamorphosis/AnamHermite.hpp"
  #include "Anamorphosis/AnamEmpirical.hpp"
  #include "Anamorphosis/AnamDiscreteDD.hpp"
  #include "Anamorphosis/AnamDiscreteIR.hpp"
  #include "Anamorphosis/EAnam.hpp"
  
  #include "Gibbs/GibbsMMulti.hpp"
  #include "Gibbs/GibbsUMulti.hpp"
  
  #include "Morpho/Morpho.hpp"
  #include "Polygon/Polygons.hpp"
  #include "Polygon/PolySet.hpp"
  
  #include "Stats/Classical.hpp"
  #include "Stats/PCA.hpp"
  #include "Stats/PCAStringFormat.hpp"
  #include "Stats/Selectivity.hpp"
  #include "Stats/EStatOption.hpp"
  #include "Stats/ESelectivity.hpp"
  
  #include "LithoRule/Rule.hpp"
  #include "LithoRule/RuleStringFormat.hpp"
  #include "LithoRule/RuleProp.hpp"
  #include "LithoRule/ERule.hpp"
  
  #include "Estimation/KrigingSystem.hpp"
  #include "Estimation/CalcKriging.hpp"
  
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
  
  #include "Polynomials/Chebychev.hpp"
  #include "Simulation/ACalcSimulation.hpp"
  #include "Simulation/CalcSimuTurningBands.hpp"
  #include "Simulation/TurningDirection.hpp"
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
  #include "Simulation/SimuFFT.hpp"
  #include "Simulation/SimuRefineParam.hpp"
  #include "Simulation/SimuRefine.hpp"
  #include "Simulation/CalcSimuEden.hpp"
  
  #include "Fractures/FracEnviron.hpp"
  #include "Fractures/FracFamily.hpp"
  #include "Fractures/FracFault.hpp"
  #include "Fractures/FracDesc.hpp"
  #include "Fractures/FracList.hpp"
  
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


////////////////////////////////////////////////
// Conversion Target language => C++

// Note : Before including this file :
//        - vectorToCpp, vectorVectorToCpp and convertToCpp 
//          functions must be defined in ToCpp fragment

// Convert scalar arguments by value
%typemap(in, fragment="ToCpp") int,
                               double,
                               String,
                               float,
                               UChar,
                               bool
{
  const int errcode = convertToCpp($input, $1);
  if (!SWIG_IsOK(errcode))
    %argument_fail(errcode, "$type", $symname, $argnum);
}

// Convert scalar argument by reference
%typemap(in, fragment="ToCpp") int*       (int val), const int*       (int val),
                               int&       (int val), const int&       (int val),
                               double* (double val), const double* (double val),
                               double& (double val), const double& (double val), // Don't add String here otherwise "res2 not declared"
                               float*   (float val), const float*   (float val),
                               float&   (float val), const float&   (float val),
                               UChar*   (UChar val), const UChar*   (UChar val),
                               UChar&   (UChar val), const UChar&   (UChar val),
                               bool*     (bool val), const bool*    (bool val),
                               bool&     (bool val), const bool&    (bool val)
{
  const int errcode = convertToCpp($input, val);
  if (!SWIG_IsOK(errcode))
    %argument_fail(errcode, "$type", $symname, $argnum);
  $1 = &val;
}

%typemap(in, fragment="ToCpp") VectorInt,
                               VectorDouble,
                               VectorString,
                               VectorFloat,
                               VectorUChar,
                               VectorBool
{ 
  const int errcode = vectorToCpp($input, $1);
  if (!SWIG_IsOK(errcode))
    %argument_fail(errcode, "$type", $symname, $argnum);
}

%typemap(in, fragment="ToCpp") VectorVectorInt,
                               VectorVectorDouble,
                               VectorVectorFloat
{ 
  const int errcode = vectorVectorToCpp($input, $1);
  if (!SWIG_IsOK(errcode))
    %argument_fail(errcode, "$type", $symname, $argnum);
}

%typemap(in, fragment="ToCpp") const VectorInt&    (VectorInt vec),
                               const VectorInt*    (VectorInt vec),
                               const VectorDouble& (VectorDouble vec),
                               const VectorDouble* (VectorDouble vec),
                               const VectorString& (VectorString vec),
                               const VectorString* (VectorString vec),
                               const VectorFloat&  (VectorFloat vec),
                               const VectorFloat*  (VectorFloat vec),
                               const VectorUChar&  (VectorUChar vec),
                               const VectorUChar*  (VectorUChar vec),
                               const VectorBool&   (VectorBool vec),
                               const VectorBool*   (VectorBool vec)
{
  const int errcode = vectorToCpp($input, vec);
  if (!SWIG_IsOK(errcode))
    %argument_fail(errcode, "$type", $symname, $argnum);
  $1 = &vec;
}

%typemap(in, fragment="ToCpp") const VectorVectorInt&    (VectorVectorInt vec),
                               const VectorVectorInt*    (VectorVectorInt vec),
                               const VectorVectorDouble& (VectorVectorDouble vec),
                               const VectorVectorDouble* (VectorVectorDouble vec),
                               const VectorVectorFloat&  (VectorVectorFloat vec),
                               const VectorVectorFloat*  (VectorVectorFloat vec)
{
  const int errcode = vectorVectorToCpp($input, vec);
  if (!SWIG_IsOK(errcode))
    %argument_fail(errcode, "$type", $symname, $argnum);
  $1 = &vec;
}

////////////////////////////////////////////////
// Conversion C++ => Target language

// Note : Before including this file :
//        - vectorFromCpp, vectorVectorFromCpp and objectFromCpp 
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

%typemap(out, fragment="FromCpp") int*,    const int*,    int&,    const int&,
                                  double*, const double*, double&, const double&,
                                  String*, const String*, String&, const String&,
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
  const int errcode = vectorFromCpp(&($result), $1);
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
  const int errcode = vectorFromCpp(&($result), *$1);
  if (!SWIG_IsOK(errcode))
    SWIG_exception_fail(SWIG_ArgError(errcode), "in method $symname, wrong return value: $type");
}

%typemap(out, fragment="FromCpp") VectorVectorInt,
                                  VectorVectorDouble,
                                  VectorVectorFloat
{
  const int errcode = vectorVectorFromCpp(&($result), $1);
  if (!SWIG_IsOK(errcode))
    SWIG_exception_fail(SWIG_ArgError(errcode), "in method $symname, wrong return value: $type");
}

%typemap(out, fragment="FromCpp") VectorVectorInt*,    VectorVectorInt&,
                                  VectorVectorDouble*, VectorVectorDouble&,
                                  VectorVectorFloat*, VectorVectorFloat&
{
  const int errcode = vectorVectorFromCpp(&($result), *$1);
  if (!SWIG_IsOK(errcode))
    SWIG_exception_fail(SWIG_ArgError(errcode), "in method $symname, wrong return value: $type");
}

