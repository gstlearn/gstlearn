// Sensitive file ! Keep Order


// Global files to be wrapped from C++ library
// Remind that swig %include doesn't follow #include inclusion.
// You must cite below each single header file that you want to export!
// Put low level headers in first positions (otherwise Syntax error in input(1).)

%include gstlearn_export.hpp // Do not forget this file in priority (for SWIG preprocessor)

// Export VectorXXX classes
%include Basic/VectorT.hpp
typedef unsigned char UChar; // Add to allow aliasing Boolean and UChar
typedef std::string String;  // Add to allow aliasing String and std::string
%template(VectorTInt)         gstlrn::VectorT< long long >;
%template(VectorTDouble)      gstlrn::VectorT< double >;
%template(VectorTFloat)       gstlrn::VectorT< float >;
%template(VectorBool)         gstlrn::VectorT< UChar >; // See VectorT.hpp
%template(VectorString)       gstlrn::VectorT< String >;

%include Basic/VectorNumT.hpp
%template(VectorInt)          gstlrn::VectorNumT< long long >;
%template(VectorDouble)       gstlrn::VectorNumT< double >;
%template(VectorFloat)        gstlrn::VectorNumT< float >;
%template(VectorUChar)        gstlrn::VectorNumT< UChar >;

%template(VectorTVectorInt)    gstlrn::VectorT< VectorNumT< long long > >;
%template(VectorTVectorDouble) gstlrn::VectorT< VectorNumT< double > >;
%template(VectorTVectorFloat)  gstlrn::VectorT< VectorNumT< float > >;

%template(VectorVectorInt)    gstlrn::VectorNumT< VectorNumT< long long > >;
%template(VectorVectorDouble) gstlrn::VectorNumT< VectorNumT< double > >;
%template(VectorVectorFloat)  gstlrn::VectorNumT< VectorNumT< float > >;

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
%include Enum/ECSV.hpp
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
%include Enum/EFormatNF.hpp
%include Enum/ESimuType.hpp
%include Enum/ERole.hpp

%include Basic/ArgumentTest.hpp
%include Basic/AStringable.hpp
%include Basic/AStringFormat.hpp
%include Basic/ASerializable.hpp
%include Basic/Message.hpp
%include Basic/NamingConvention.hpp

%include Transform/ATransform.hpp
%include Transform/ATransformWithAutoDiff.hpp
%include Transform/TuckeyGH.hpp
%include Transform/YeoJohnson.hpp

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
%include Basic/LawStable.hpp
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
%include Space/SpaceSN.hpp

%include LinearOp/ALinearOp.hpp
%include LinearOp/APreconditioner.hpp
%include LinearOp/MultiGridSolver.hpp
%include LinearOp/MultiGridSPDE.hpp
%include LinearOp/ASimulable.hpp
%include LinearOp/ASimulableMatrix.hpp
%include Matrix/AMatrix.hpp
%include Matrix/MatrixDense.hpp
%include Matrix/MatrixSparse.hpp
%include LinearOp/InvNuggetOp.hpp
%include Matrix/MatrixSquare.hpp
%include Matrix/NF_Triplet.hpp
%include Matrix/MatrixSymmetric.hpp
%include Matrix/MatrixFactory.hpp
%include Matrix/MatrixInt.hpp
%include Matrix/Table.hpp
%include LinearOp/LinearOpHelper.hpp

%include MLayers/MLayers.hpp

%include Skin/ISkinFunctions.hpp
%include Skin/Skin.hpp

%include Mesh/AMesh.hpp
%include Mesh/MeshEStandard.hpp
%include Mesh/MeshEFaulted.hpp
%include Mesh/MeshETurbo.hpp
%include Mesh/MeshSpherical.hpp
%include Mesh/MeshSphericalExt.hpp
%include Mesh/VectorMeshes.hpp

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
%include LinearOp/IPrecisionOp.hpp
%include LinearOp/PrecisionOp.hpp
%include LinearOp/PrecisionOpMatrix.hpp
%include LinearOp/SPDEOp.hpp
%include LinearOp/SPDEOpMatrix.hpp
%include LinearOp/TurboOptimizer.hpp
%include LinearOp/IProj.hpp
%include LinearOp/ProjZero.hpp
%include LinearOp/ProjComposition.hpp
%include LinearOp/ProjMatrix.hpp
%include LinearOp/ProjMulti.hpp
%include LinearOp/ProjMultiMatrix.hpp
%include LinearOp/PrecisionOpMulti.hpp
%include LinearOp/PrecisionOpMultiMatrix.hpp
%include LinearOp/ProjConvolution.hpp
%include LinearOp/IOptimCost.hpp
%include LinearOp/OptimCostBinary.hpp
%include LinearOp/OptimCostColored.hpp
%include LinearOp/MatrixSymmetricSim.hpp
%include LinearOp/ACholesky.hpp
%include LinearOp/CholeskyDense.hpp
%include LinearOp/CholeskySparse.hpp

%include Neigh/ANeigh.hpp
%include Neigh/NeighUnique.hpp
%include Neigh/NeighImage.hpp
%include Neigh/NeighMoving.hpp
%include Neigh/NeighBench.hpp
%include Neigh/NeighCell.hpp

%include Variogram/Variograms.hpp
%include Variogram/AVario.hpp
%include Variogram/VarioParam.hpp
%include Variogram/Vario.hpp
%include Variogram/VarioParam.hpp
%include Variogram/DirParam.hpp
%include Variogram/VMap.hpp
%include Variogram/VCloud.hpp
%include Variogram/VarioOrder.hpp

%include Basic/ParamInfo.hpp
%include Basic/ListParams.hpp
%include Model/ModelGeneric.hpp
%include Model/GaussianProcess.hpp
%include Model/ModelCovList.hpp
%include Model/Model.hpp
%include Model/ModelOptimParam.hpp
%include Model/ElemNostat.hpp
%include Model/Option_AutoFit.hpp
%include Model/Option_VarioFit.hpp
%include Model/Constraints.hpp
%include Model/ConsItem.hpp
%include Model/CovParamId.hpp
%include Model/CovParamId.hpp

%include Covariances/ParamId.hpp
%include Covariances/TabNoStat.hpp
%include Covariances/TabNoStatCovAniso.hpp
%include Covariances/TabNoStatSills.hpp
%include Covariances/ANoStat.hpp
%include Covariances/NoStatArray.hpp
%include Covariances/NoStatOnMesh.hpp
%include Covariances/NoStatFunctional.hpp
%include Covariances/ACov.hpp
%include Covariances/CovBase.hpp
%include Covariances/CovProportional.hpp
%include Covariances/AKernel.hpp
%include Covariances/CovList.hpp
%include Covariances/CovAnisoList.hpp
%include Covariances/CovAniso.hpp
%include Covariances/CovGradientGeneric.hpp
%include Covariances/CovGradientAnalytic.hpp
%include Covariances/CorAniso.hpp
%include Covariances/CorFactorized.hpp
%include Covariances/CorGaussianMixture.hpp
%include Covariances/CorGneiting.hpp
%include Covariances/CorMatern.hpp
%include Covariances/CovLMCTapering.hpp
%include Covariances/CovLMCConvolution.hpp
%include Covariances/CovLMCAnamorphosis.hpp
%include Covariances/CovContext.hpp
%include Covariances/CovCalcMode.hpp
%include Covariances/KernelBesselJ.hpp
%include Covariances/KernelMatern.hpp
%include Covariances/KernelCauchy.hpp
%include Covariances/KernelCauchyGen.hpp
%include Covariances/KernelCosExp.hpp
%include Covariances/KernelCosinus.hpp
%include Covariances/KernelCubic.hpp
%include Covariances/KernelExponential.hpp
%include Covariances/KernelGamma.hpp
%include Covariances/KernelGaussian.hpp
%include Covariances/KernelGC1.hpp
%include Covariances/KernelGC3.hpp
%include Covariances/KernelGC5.hpp
%include Covariances/KernelGCspline2.hpp
%include Covariances/KernelGCspline.hpp
%include Covariances/KernelLinear.hpp
%include Covariances/KernelNugget.hpp
%include Covariances/KernelPenta.hpp
%include Covariances/KernelPower.hpp
%include Covariances/KernelReg1D.hpp
%include Covariances/KernelSincard.hpp
%include Covariances/KernelSpherical.hpp
%include Covariances/KernelStable.hpp
%include Covariances/KernelStorkey.hpp
%include Covariances/KernelTriangle.hpp
%include Covariances/KernelWendland0.hpp
%include Covariances/KernelWendland1.hpp
%include Covariances/KernelWendland2.hpp
%include Covariances/KernelMarkov.hpp
%include Covariances/KernelGeometric.hpp
%include Covariances/KernelPoisson.hpp
%include Covariances/KernelLinearSph.hpp
%include Covariances/CovDiffusionAdvection.hpp
%include Covariances/CovHelper.hpp

%include Drifts/ADrift.hpp
%include Drifts/DriftList.hpp
%include Drifts/DriftM.hpp
%include Drifts/DriftF.hpp
%include Drifts/DriftFactory.hpp

%include API/SPDE.hpp
%include API/TestInheritance.hpp
%include API/Style.hpp
%include API/SPDEParam.hpp
%include API/Potential.hpp

// 1. Load typemaps BEFOREHAND so that SWIG knows them
#ifdef SWIGPYTHON
%include "typemaps_optional_python.i"
%include "typemaps_colID_python.i"
#endif

#ifdef SWIGR
%include "typemaps_optional_r.i"
%include "typemaps_colID_r.i"
#endif

// 2. Define types
%include "DataBase/RoleID.hpp"
%include "DataBase/ColID.hpp"

// 3. Mpad the headers that use ColID
%include "DataBase/DbData.hpp"
%include "DataBase/DbCol.hpp"
//%include DataBase/Dictionary.hpp
//%include DataBase/VectorCategory.hpp

%include Db/Db.hpp
%include Db/DbGrid.hpp
%include Db/DbLine.hpp
%include Db/DbGraphO.hpp
%include Db/DbMeshTurbo.hpp
%include Db/DbMeshStandard.hpp
%include Db/DbStringFormat.hpp
%include Db/DbHelper.hpp
%include Db/RankHandler.hpp

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
%include Gibbs/GibbsMultiMono.hpp
%include Gibbs/GibbsUMultiMono.hpp
%include Gibbs/GibbsUPropMono.hpp

%include Morpho/Morpho.hpp

%include Polygon/Polygons.hpp
%include Polygon/PolyElem.hpp

%include Stats/Classical.hpp
%include Stats/PCA.hpp
%include Stats/PCAStringFormat.hpp
%include Stats/Selectivity.hpp
%include Stats/Regression.hpp

%include LithoRule/Node.hpp
%include LithoRule/Rule.hpp
%include LithoRule/RuleShadow.hpp
%include LithoRule/RuleShift.hpp
%include LithoRule/RuleStringFormat.hpp
%include LithoRule/RuleProp.hpp
%include LithoRule/PropDef.hpp

%include Estimation/Estimations.hpp
%include Estimation/KrigingSystem.hpp
%include Estimation/KrigingAlgebra.hpp
%include Estimation/CalcKriging.hpp
%include Estimation/CalcKrigingFactors.hpp
%include Estimation/CalcKrigingGradient.hpp
%include Estimation/CalcSimpleInterpolation.hpp
%include Estimation/CalcImage.hpp
%include Estimation/CalcGlobal.hpp
%include Estimation/KrigOpt.hpp
%include Estimation/AModelOptim.hpp
%include Estimation/AModelOptimFactory.hpp
%include Estimation/ALikelihood.hpp
%include Estimation/Vecchia.hpp
%include Estimation/Likelihood.hpp


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

%include Simulation/Simulations.hpp
%include Simulation/ACalcSimulation.hpp
%include Simulation/ACalcSimuGaussian.hpp
%include Simulation/CalcSimuTurningBands.hpp
%include Simulation/CalcSimuPGS.hpp
%include Simulation/TurningBandDirection.hpp
%include Simulation/TurningBandOperate.hpp
%include Simulation/CalcSimuSpectral.hpp
%include Simulation/CalcSimuBoolean.hpp
%include Simulation/SimuSpectralRN.hpp
%include Simulation/SimuSpectralS2.hpp
%include Simulation/SpectrumOnRN.hpp
%include Simulation/BooleanObject.hpp
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
%include Tree/BallFaulted.hpp
%include Tree/KNN.hpp

%include Spatial/Projection.hpp
%include Spatial/SpatialIndices.hpp

%include PluriGaussian/TracePGS.hpp
%include PluriGaussian/CorPGS.hpp
%include PluriGaussian/DiscretePGS.hpp
%include PluriGaussian/CalcModelPGS.hpp

%include Core/Acknowledge.hpp
%include Core/Seismic.hpp

%include/API/newAPIs.hpp

// For suppressing SWIG warning due to -keyword option (if used)
#pragma SWIG nowarn=511
#pragma SWIG nowarn=506
#pragma SWIG nowarn=509

%template(LinearOpCGSolver) LinearOpCGSolver< ScaleOp >;
%template(LinearSPDEOpCGSolver) LinearOpCGSolver< SPDEOp >;

// In DbData: addColumn()
%template(addColumnEmptyD)  gstlrn::DbData::addColumnEmpty<VectorDouble>;
%template(addColumnEmptyF)  gstlrn::DbData::addColumnEmpty<VectorFloat>;
%template(addColumnEmptyI)  gstlrn::DbData::addColumnEmpty<VectorInt>;
%template(addColumnEmptyU)  gstlrn::DbData::addColumnEmpty<VectorUChar>;
%template(addColumnEmptyS)  gstlrn::DbData::addColumnEmpty<VectorString>;
%template(addColumnEmptyB)  gstlrn::DbData::addColumnEmpty<VectorBool>;
//%template(addColumnEmptyC)  gstlrn::DbData::addColumnEmpty<VectorCategory>;

%template(addColumnD)  gstlrn::DbData::addColumn<VectorDouble>;
%template(addColumnF)  gstlrn::DbData::addColumn<VectorFloat>;
%template(addColumnI)  gstlrn::DbData::addColumn<VectorInt>;
%template(addColumnU)  gstlrn::DbData::addColumn<VectorUChar>;
%template(addColumnS)  gstlrn::DbData::addColumn<VectorString>;
%template(addColumnB)  gstlrn::DbData::addColumn<VectorBool>;
//%template(addColumnC)  gstlrn::DbData::addColumn<VectorCategory>;
