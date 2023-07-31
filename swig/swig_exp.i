// Sensitive file ! Keep Order

// Global files to be wrapped from C++ library
// Remind that swig %include doesn't follow #include inclusion.
// You must cite below each single header file that you want to export!
// Put low level headers in first positions (otherwise Syntax error in input(1).)
%include gstlearn_export.hpp // Do not forget this file in priority (for SWIG preprocessor)

// Export VectorXXX classes
%include Basic/VectorT.hpp
%include Basic/VectorNumT.hpp
%template(VectorTInt)         VectorT< int >;
%template(VectorTDouble)      VectorT< double >;
%template(VectorString)       VectorT< String >;
%template(VectorTFloat)       VectorT< float >;
%template(VectorTUChar)       VectorT< UChar >;
%template(VectorBool)         VectorT< UChar >; // See VectorT.hpp
%template(VectorInt)          VectorNumT< int >;
%template(VectorDouble)       VectorNumT< double >;
%template(VectorFloat)        VectorNumT< float >;
%template(VectorUChar)        VectorNumT< UChar >;
%template(VectorVectorInt)    VectorT< VectorNumT< int > >;
%template(VectorVectorDouble) VectorT< VectorNumT< double > >;
%template(VectorVectorFloat)  VectorT< VectorNumT< float > >;

%include Basic/ICloneable.hpp
%include Basic/RepeatMacro.hpp
%include Basic/RepeatMacroSwig.hpp

%include geoslib_define.h
%include geoslib_enum.h
%include geoslib_d.h
%include geoslib_f.h

%include Enum/AEnum.hpp
%include Enum/EKrigOpt.hpp
%include Enum/ESPDECalcMode.hpp
%include Enum/EAnam.hpp
%include Enum/ECst.hpp
%include Enum/EDbg.hpp
%include Enum/ELaw.hpp
%include Enum/EShape.hpp
%include Enum/EConvDir.hpp
%include Enum/ECalcVario.hpp
%include Enum/EConvType.hpp
%include Enum/ECov.hpp
%include Enum/ETape.hpp
%include Enum/ELoadBy.hpp
%include Enum/ELoc.hpp
%include Enum/EDrift.hpp
%include Enum/EPowerPT.hpp
%include Enum/ERule.hpp
%include Enum/EConsElem.hpp
%include Enum/EConsType.hpp
%include Enum/EModelProperty.hpp
%include Enum/EMorpho.hpp
%include Enum/ENeigh.hpp
%include Enum/ESpaceType.hpp
%include Enum/ESelectivity.hpp
%include Enum/EStatOption.hpp
%include Enum/EDirGen.hpp
%include Enum/EGaussInv.hpp
%include Enum/ECalcMember.hpp
%include Enum/EPostUpscale.hpp
%include Enum/EPostStat.hpp

%include Basic/ArgumentTest.hpp
%include Basic/AStringable.hpp
%include Basic/AStringFormat.hpp
%include Basic/ASerializable.hpp

%include Basic/NamingConvention.hpp

%include Calculators/ACalculator.hpp
%include Calculators/ACalcDbVarCreator.hpp
%include Calculators/ACalcDbToDb.hpp
%include Calculators/CalcMigrate.hpp
%include Calculators/ACalcInterpolator.hpp
%include Calculators/CalcStatistics.hpp
%include Calculators/CalcGridToGrid.hpp
%include Calculators/CalcSimuPost.hpp
%include Calculators/CalcSimuPostDemo.hpp
%include Calculators/CalcSimuPostPropByLayer.hpp

%include Basic/Tensor.hpp
%include Basic/Grid.hpp
%include Basic/String.hpp
%include Basic/Interval.hpp
%include Basic/Limits.hpp
%include Basic/Utilities.hpp
%include Basic/CSVformat.hpp
%include Basic/AFunction.hpp
%include Basic/AFunctional.hpp
%include Basic/FunctionalSpirale.hpp
%include Basic/RepeatMacro.hpp
%include Basic/RepeatMacroSwig.hpp
%include Basic/OptDbg.hpp
%include Basic/OptCst.hpp
%include Basic/OptCustom.hpp
%include Basic/File.hpp
%include Basic/VectorHelper.hpp
%include Basic/Plane.hpp
%include Basic/FFT.hpp
%include Basic/PolyLine2D.hpp
%include Basic/Law.hpp
%include Basic/MathFunc.hpp
%include Basic/Indirection.hpp

%include Geometry/GeometryHelper.hpp
%include Geometry/Rotation.hpp
%include Geometry/ABiTargetCheck.hpp
%include Geometry/BiTargetCheckBench.hpp
%include Geometry/BiTargetCheckCell.hpp
%include Geometry/BiTargetCheckDistance.hpp
%include Geometry/BiTargetCheckFaults.hpp
%include Geometry/BiTargetCheckCode.hpp
%include Geometry/BiTargetCheckDate.hpp

%include Arrays/AArray.hpp
%include Arrays/Array.hpp
%include Arrays/BImage.hpp
%include Arrays/BImageStringFormat.hpp

%include Faults/Faults.hpp

%include Boolean/ShapeParameter.hpp
%include Boolean/AShape.hpp
%include Boolean/ShapeParallelepiped.hpp
%include Boolean/ShapeEllipsoid.hpp
%include Boolean/ShapeParaboloid.hpp
%include Boolean/ShapeHalfEllipsoid.hpp
%include Boolean/ShapeHalfParaboloid.hpp
%include Boolean/ShapeHalfSinusoid.hpp
%include Boolean/ModelBoolean.hpp

%include Space/ASpace.hpp
%include Space/ASpaceObject.hpp
%include Space/SpacePoint.hpp
%include Space/SpaceTarget.hpp
%include Space/SpaceRN.hpp
%include Space/SpaceShape.hpp

%include Skin/ISkinFunctions.hpp
%include Skin/Skin.hpp

%include Mesh/AMesh.hpp
%include Mesh/MeshEStandard.hpp
%include Mesh/MeshETurbo.hpp
%include Mesh/MeshSpherical.hpp
%include Mesh/MeshSphericalExt.hpp

%include Polynomials/APolynomial.hpp
%include Polynomials/ClassicalPolynomial.hpp
%include Polynomials/Hermite.hpp
%include Polynomials/MonteCarlo.hpp
%include Polynomials/Chebychev.hpp

%include LinearOp/ALinearOp.hpp
%include LinearOp/ALinearOpMulti.hpp
%include LinearOp/ShiftOpCs.hpp
%include LinearOp/PrecisionOp.hpp
%include LinearOp/PrecisionOpCs.hpp
%include LinearOp/TurboOptimizer.hpp
%include LinearOp/IProjMatrix.hpp
%include LinearOp/ProjMatrix.hpp
%include LinearOp/PrecisionOpMultiConditional.hpp
%include LinearOp/ProjConvolution.hpp
%include LinearOp/IOptimCost.hpp
%include LinearOp/OptimCostBinary.hpp
%include LinearOp/OptimCostColored.hpp
%include LinearOp/Cholesky.hpp

%include Model/ANoStat.hpp
%include Model/NoStatArray.hpp
%include Model/NoStatFunctional.hpp

%include Neigh/ANeigh.hpp
%include Neigh/NeighUnique.hpp
%include Neigh/NeighImage.hpp
%include Neigh/NeighMoving.hpp
%include Neigh/NeighBench.hpp
%include Neigh/NeighCell.hpp

%include Variogram/VarioParam.hpp
%include Variogram/Vario.hpp
%include Variogram/VarioParam.hpp
%include Variogram/DirParam.hpp

%include Model/Model.hpp
%include Model/Option_AutoFit.hpp
%include Model/Option_VarioFit.hpp
%include Model/Constraints.hpp
%include Model/ConsItem.hpp
%include Model/CovParamId.hpp
%include Model/CovParamId.hpp

%include Covariances/ACov.hpp
%include Covariances/ACovFunc.hpp
%include Covariances/ACovAnisoList.hpp
%include Covariances/CovAniso.hpp
%include Covariances/ACovGradient.hpp
%include Covariances/CovLMC.hpp
%include Covariances/CovLMCTapering.hpp
%include Covariances/CovLMCConvolution.hpp
%include Covariances/CovLMCAnamorphosis.hpp
%include Covariances/CovLMGradient.hpp
%include Covariances/CovContext.hpp
%include Covariances/CovCalcMode.hpp
%include Covariances/CovBesselJ.hpp
%include Covariances/CovBesselK.hpp
%include Covariances/CovCauchy.hpp
%include Covariances/CovCosExp.hpp
%include Covariances/CovCosinus.hpp
%include Covariances/CovCubic.hpp
%include Covariances/CovExponential.hpp
%include Covariances/CovGamma.hpp
%include Covariances/CovGaussian.hpp
%include Covariances/CovGC1.hpp
%include Covariances/CovGC3.hpp
%include Covariances/CovGC5.hpp
%include Covariances/CovGCspline2.hpp
%include Covariances/CovGCspline.hpp
%include Covariances/CovLinear.hpp
%include Covariances/CovNugget.hpp
%include Covariances/CovPenta.hpp
%include Covariances/CovPower.hpp
%include Covariances/CovReg1D.hpp
%include Covariances/CovSincard.hpp
%include Covariances/CovSpherical.hpp
%include Covariances/CovStable.hpp
%include Covariances/CovStorkey.hpp
%include Covariances/CovTriangle.hpp
%include Covariances/CovWendland0.hpp
%include Covariances/CovWendland1.hpp
%include Covariances/CovWendland2.hpp
%include Covariances/CovMarkov.hpp
%include Covariances/CovDiffusionAdvection.hpp
%include Covariances/CovHelper.hpp

%include Drifts/ADrift.hpp
%include Drifts/ADriftElem.hpp
%include Drifts/DriftList.hpp
%include Drifts/Drift1.hpp
%include Drifts/DriftF.hpp
%include Drifts/DriftFactory.hpp
%include Drifts/DriftX.hpp
%include Drifts/DriftX2.hpp
%include Drifts/DriftX2Y.hpp
%include Drifts/DriftX3.hpp
%include Drifts/DriftXY.hpp
%include Drifts/DriftXY2.hpp
%include Drifts/DriftXZ.hpp
%include Drifts/DriftY.hpp
%include Drifts/DriftY2.hpp
%include Drifts/DriftY3.hpp
%include Drifts/DriftYZ.hpp
%include Drifts/DriftZ.hpp
%include Drifts/DriftZ2.hpp
%include Drifts/DriftZ3.hpp

%include Matrix/AMatrix.hpp
%include Matrix/AMatrixSquare.hpp
%include Matrix/LinkMatrixSparse.hpp
%include Matrix/MatrixSparse.hpp
%include Matrix/MatrixRectangular.hpp
%include Matrix/MatrixSquareGeneral.hpp
%include Matrix/MatrixSquareSymmetric.hpp
%include Matrix/MatrixSquareDiagonal.hpp
%include Matrix/MatrixSquareDiagonalCst.hpp
%include Matrix/MatrixInt.hpp
%include Matrix/Table.hpp

%include API/SPDE.hpp
%include API/PGSSPDE.hpp
%include API/TestInheritance.hpp
%include API/Style.hpp

%include Db/Db.hpp
%include Db/DbGrid.hpp
%include Db/DbStringFormat.hpp

%include Anamorphosis/CalcAnamTransform.hpp
%include Anamorphosis/AAnam.hpp
%include Anamorphosis/AnamContinuous.hpp
%include Anamorphosis/AnamDiscrete.hpp
%include Anamorphosis/AnamUser.hpp
%include Anamorphosis/AnamHermite.hpp
%include Anamorphosis/AnamEmpirical.hpp
%include Anamorphosis/AnamDiscreteDD.hpp
%include Anamorphosis/AnamDiscreteIR.hpp
%include Anamorphosis/PPMT.hpp

%include Gibbs/AGibbs.hpp
%include Gibbs/GibbsMulti.hpp
%include Gibbs/GibbsMMulti.hpp
%include Gibbs/GibbsUMulti.hpp

%include Morpho/Morpho.hpp

%include Polygon/Polygons.hpp
%include Polygon/PolyElem.hpp

%include Stats/Classical.hpp
%include Stats/PCA.hpp
%include Stats/PCAStringFormat.hpp
%include Stats/Selectivity.hpp

%include LithoRule/Rule.hpp
%include LithoRule/RuleStringFormat.hpp
%include LithoRule/RuleProp.hpp

%include Estimation/KrigingSystem.hpp
%include Estimation/CalcKriging.hpp
%include Estimation/CalcKrigingFactors.hpp
%include Estimation/CalcSimpleInterpolation.hpp
%include Estimation/CalcImage.hpp
%include Estimation/CalcGlobal.hpp

%include OutputFormat/AOF.hpp
%include OutputFormat/FileLAS.hpp
%include OutputFormat/FileVTK.hpp
%include OutputFormat/GridArcGis.hpp
%include OutputFormat/GridBmp.hpp
%include OutputFormat/GridEclipse.hpp
%include OutputFormat/GridF2G.hpp
%include OutputFormat/GridIfpEn.hpp
%include OutputFormat/GridIrap.hpp
%include OutputFormat/GridXYZ.hpp
%include OutputFormat/GridZycor.hpp
%include OutputFormat/segy.h

%include Simulation/ACalcSimulation.hpp
%include Simulation/CalcSimuTurningBands.hpp
%include Simulation/TurningDirection.hpp
%include Simulation/BooleanObject.hpp
%include Simulation/SimuBoolean.hpp
%include Simulation/SimuBooleanParam.hpp
%include Simulation/SimuSpherical.hpp
%include Simulation/SimuSphericalParam.hpp
%include Simulation/CalcSimuSubstitution.hpp
%include Simulation/SimuSubstitutionParam.hpp
%include Simulation/CalcSimuPartition.hpp
%include Simulation/SimuPartitionParam.hpp
%include Simulation/SimuFFTParam.hpp
%include Simulation/CalcSimuFFT.hpp
%include Simulation/SimuRefineParam.hpp
%include Simulation/SimuRefine.hpp
%include Simulation/CalcSimuEden.hpp

%include Fractures/FracEnviron.hpp
%include Fractures/FracFamily.hpp
%include Fractures/FracFault.hpp
%include Fractures/FracDesc.hpp
%include Fractures/FracList.hpp

%include Skin/Skin.hpp

%include Tree/Ball.hpp

// For suppressing SWIG warning due to -keyword option (if used)
#pragma SWIG nowarn=511
#pragma SWIG nowarn=506
#pragma SWIG nowarn=509
