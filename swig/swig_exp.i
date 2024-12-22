// Sensitive file ! Keep Order

// Global files to be wrapped from C++ library
// Remind that swig %include doesn't follow #include inclusion.
// Put low level headers in first positions (otherwise Syntax error in input(1).)
%include gstlearn_export.hpp // Do not forget this file in priority (for SWIG preprocessor)

// Export VectorXXX classes
%include Basic/VectorT.hpp
%include Basic/VectorNumT.hpp
%template(VectorTInt)         VectorT< int >;
%template(VectorTDouble)      VectorT< double >;
%template(VectorTFloat)       VectorT< float >;
%template(VectorBool)         VectorT< UChar >; // See VectorT.hpp
%template(VectorString)       VectorT< String >;
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
%include Enum/EOperator.hpp
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
%include Geometry/BiTargetCheckGeometry.hpp

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
%include Space/SpaceComposite.hpp
%include Space/ASpaceObject.hpp
%include Space/SpacePoint.hpp
%include Space/SpaceTarget.hpp
%include Space/SpaceRN.hpp
%include Space/SpaceShape.hpp

%include LinearOp/ALinearOp.hpp
%include LinearOp/ASimulable.hpp
%include Matrix/AMatrix.hpp
%include Matrix/AMatrixDense.hpp
%include Matrix/MatrixSparse.hpp
%include Matrix/MatrixRectangular.hpp
%include Matrix/AMatrixSquare.hpp
%include Matrix/NF_Triplet.hpp
%include Matrix/MatrixSquareGeneral.hpp
%include Matrix/MatrixSquareSymmetric.hpp
%include Matrix/MatrixFactory.hpp
%include Matrix/MatrixInt.hpp
%include Matrix/Table.hpp

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

%include LinearOp/CGParam.hpp
%include LinearOp/LogStats.hpp
%include LinearOp/LinearOpCGSolver.hpp
%include LinearOp/ALinearOpMulti.hpp
%include LinearOp/ScaleOp.hpp
%include LinearOp/AShiftOp.hpp
%include LinearOp/ShiftOpStencil.hpp
%include LinearOp/ShiftOpMatrix.hpp
%include LinearOp/PrecisionOp.hpp
%include LinearOp/PrecisionOpMatrix.hpp
%include LinearOp/SPDEOp.hpp
%include LinearOp/SPDEOpMatrix.hpp
%include LinearOp/TurboOptimizer.hpp
%include LinearOp/IProj.hpp
%include LinearOp/ProjMatrix.hpp
%include LinearOp/ProjMulti.hpp
%include LinearOp/ProjMultiMatrix.hpp
%include LinearOp/PrecisionOpMulti.hpp
%include LinearOp/PrecisionOpMultiMatrix.hpp
%include LinearOp/PrecisionOpMultiConditional.hpp
%include LinearOp/ProjConvolution.hpp
%include LinearOp/IOptimCost.hpp
%include LinearOp/OptimCostBinary.hpp
%include LinearOp/OptimCostColored.hpp
%include LinearOp/MatrixSquareSymmetricSim.hpp

%include Neigh/ANeigh.hpp
%include Neigh/NeighUnique.hpp
%include Neigh/NeighImage.hpp
%include Neigh/NeighMoving.hpp
%include Neigh/NeighBench.hpp
%include Neigh/NeighCell.hpp

%include Variogram/AVario.hpp
%include Variogram/VarioParam.hpp
%include Variogram/Vario.hpp
%include Variogram/VarioParam.hpp
%include Variogram/DirParam.hpp
%include Variogram/VMap.hpp
%include Variogram/VCloud.hpp

%include Model/Model.hpp
%include Model/Option_AutoFit.hpp
%include Model/Option_VarioFit.hpp
%include Model/Constraints.hpp
%include Model/ConsItem.hpp
%include Model/CovParamId.hpp
%include Model/CovParamId.hpp

%include Covariances/ParamId.hpp
%include Covariances/TabNoStat.hpp
%include Covariances/TabNoStatCovAniso.hpp
%include Covariances/ANoStat.hpp
%include Covariances/NoStatArray.hpp
%include Covariances/NoStatFunctional.hpp
%include Covariances/ACov.hpp
%include Covariances/CovBase.hpp
%include Covariances/ACor.hpp
%include Covariances/CorAniso.hpp
%include Covariances/ACovFunc.hpp
%include Covariances/ACovAnisoList.hpp
%include Covariances/CovAniso.hpp
%include Covariances/ACovGradient.hpp
%include Covariances/CovGneiting.hpp
%include Covariances/CovLMCTapering.hpp
%include Covariances/CovLMCConvolution.hpp
%include Covariances/CovLMCAnamorphosis.hpp
%include Covariances/CovLMGradient.hpp
%include Covariances/CovContext.hpp
%include Covariances/CovCalcMode.hpp
%include Covariances/CovBesselJ.hpp
%include Covariances/CovMatern.hpp
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
%include Covariances/CovGeometric.hpp
%include Covariances/CovPoisson.hpp
%include Covariances/CovLinearSph.hpp
%include Covariances/CovDiffusionAdvection.hpp
%include Covariances/CovHelper.hpp

%include Drifts/ADrift.hpp
%include Drifts/DriftList.hpp
%include Drifts/DriftM.hpp
%include Drifts/DriftF.hpp
%include Drifts/DriftFactory.hpp

%include API/SPDE.hpp
%include API/PGSSPDE.hpp
%include API/TestInheritance.hpp
%include API/Style.hpp
%include API/SPDEParam.hpp

%include Db/Db.hpp
%include Db/DbGrid.hpp
%include Db/DbLine.hpp
%include Db/DbGraphO.hpp
%include Db/DbMeshTurbo.hpp
%include Db/DbMeshStandard.hpp
%include Db/DbStringFormat.hpp
%include Db/DbHelper.hpp

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
%include Stats/Regression.hpp

%include LithoRule/Rule.hpp
%include LithoRule/RuleStringFormat.hpp
%include LithoRule/RuleProp.hpp

%include Estimation/KrigingSystem.hpp
%include Estimation/KrigingCalcul.hpp
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
%include Simulation/TurningBandDirection.hpp
%include Simulation/TurningBandOperate.hpp
%include Simulation/SimuSpectral.hpp
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
%include Simulation/CalcSimuRefine.hpp
%include Simulation/CalcSimuEden.hpp

%include Fractures/FracEnviron.hpp
%include Fractures/FracFamily.hpp
%include Fractures/FracFault.hpp
%include Fractures/FracDesc.hpp
%include Fractures/FracList.hpp

%include Tree/Ball.hpp
%include Tree/KNN.hpp

%include Spatial/Projection.hpp
%include Spatial/SpatialIndices.hpp

%include Core/Acknowledge.hpp
%include Core/Potential.hpp
%include Core/Seismic.hpp

// For suppressing SWIG warning due to -keyword option (if used)
#pragma SWIG nowarn=511
#pragma SWIG nowarn=506
#pragma SWIG nowarn=509


%template(LinearOpCGSolver) LinearOpCGSolver< ScaleOp >;
%template(LinearSPDEOpCGSolver) LinearOpCGSolver< SPDEOp >;
