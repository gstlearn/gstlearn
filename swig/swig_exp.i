%ignore *::operator=;


%{
#define SWIG_FILE_WITH_INIT

%}


// TODO: How to mask IClonable and clone method?
class IClonable{};
//%ignore *::clone;

%include stl.i
// Cast strings into native type of the target language
%include std_string.i
// Cast vectors of integers into native type of the target language
%include std_vector.i

// Keep order in the file
%template(VectorDouble)         std::vector< double >;
%template(VectorInt)            std::vector< int >;
%template(VectorString)         std::vector< std::string >;
%template(VectorBool)           std::vector< bool >;
%template(VectorUChar)          std::vector< unsigned char >;
%template(VectorVectorInt)      std::vector< std::vector< int > >;
%template(VectorVectorDouble)   std::vector< std::vector< double > >;

%template(VectorEnumCovs)       std::vector< ECov >;    // Not a pointers list

%template(VectorCTable)         std::vector<CTable*>;
%template(VectorDir)            std::vector<DirParam>;  // Not a pointers list
%template(VectorDirection)      std::vector<Direction*>;
%template(VectorDrft)           std::vector<Drift*>;
%template(VectorFrac_Desc)      std::vector<Frac_Desc*>;
%template(VectorFrac_Fam)       std::vector<Frac_Fam*>;
%template(VectorFrac_Fault)     std::vector<Frac_Fault*>;
%template(VectorLocal_Split)    std::vector<Local_Split*>;
%template(VectorPolySet)        std::vector<PolySet*>;
%template(VectorQChol)          std::vector<QChol*>;
%template(VectorSPDE_SS_Option) std::vector<SPDE_SS_Option*>;
%template(VectorSubPlan)        std::vector<SubPlan*>;
%template(VectorToken_Def)      std::vector<Token_Def*>;
%template(VectorToken_Par)      std::vector<Token_Par*>;
%template(VectorIntervals)      std::vector<Interval*>;


//%include "numpy.i"

//%init %{
//import_array();
//%}


//%apply (double IN_ARRAY1[ANY]){(VectorDouble(ANY))};
//%apply (double IN_ARRAY2[ANY][ANY){(VectorVectorDouble(ANY))};


//%apply (double IN_ARRAY2[ANY][ANY]){(AMatrix(ANY,ANY))};
//%apply (double IN_ARRAY2[ANY][ANY]){(MatrixSquareGeneral(ANY,ANY))};


// Remind that swig %include doesn't follow #include inclusion.
// You must cite below each single header file that you want to export!
// Put low level headers in first positions (otherwise Syntax error in input(1).)
%include gstlearn_export.hpp
%include Basic/Vector.hpp
%include csparse_d.h
%include csparse_f.h
%include geoslib_define.h
%include geoslib_enum.h
%include geoslib_d.h
%include geoslib_f.h

%include Basic/RepeatMacro.hpp
%include Basic/RepeatMacroSwig.hpp

%include Enum/AEnum.hpp
%include Enum/EKrigOpt.hpp

%include Basic/ArgumentTest.hpp
%include Basic/AStringable.hpp
%include Basic/AStringFormat.hpp
%include Basic/ASerializable.hpp
%include Basic/NamingConvention.hpp
%include Basic/Rotation.hpp
%include Basic/Tensor.hpp
%include Basic/Grid.hpp
%include Basic/String.hpp
%include Basic/Interval.hpp
%include Basic/Limits.hpp
%include Basic/Utilities.hpp
%include Basic/CSVformat.hpp
%include Basic/AFunctional.hpp
%include Basic/FunctionalSpirale.hpp
%include Basic/RepeatMacro.hpp
%include Basic/RepeatMacroSwig.hpp
%include Basic/Table.hpp
%include Basic/Utilities.hpp
%include Basic/OptDbg.hpp
%include Basic/OptCst.hpp
%include Basic/OptCustom.hpp
%include Basic/EDbg.hpp
%include Basic/ECst.hpp
%include Basic/File.hpp

%include Space/Space.hpp
%include Space/ASpace.hpp
%include Space/ASpaceObject.hpp
%include Space/SpacePoint.hpp
%include Space/SpaceRN.hpp
%include Space/SpaceShape.hpp
/*
%include Interfaces/geoslib_f_swig.h
%include Interfaces/ACalculator.hpp
%include Interfaces/AParam.hpp
%include Interfaces/AVariable.hpp
%include Interfaces/AVariableTemplate.hpp
%include Interfaces/Category.hpp
%include Interfaces/Database.hpp
%include Interfaces/Dictionary.hpp
%include Interfaces/interface_d.hpp
%include Interfaces/ParamCSV.hpp
%include Interfaces/ParamGrid.hpp
%include Interfaces/Param.hpp
*/
%include Mesh/AMesh.hpp
%include Mesh/MeshFactory.hpp
%include Mesh/MeshEStandard.hpp
%include Mesh/MeshETurbo.hpp

%include Polynomials/Hermite.hpp
%include Polynomials/MonteCarlo.hpp

%include LinearOp/ALinearOp.hpp
%include LinearOp/ALinearOpMulti.hpp
%include LinearOp/ShiftOpCs.hpp
%include LinearOp/PrecisionOp.hpp
%include LinearOp/PrecisionOpCs.hpp
%include LinearOp/TurboOptimizer.hpp
%include LinearOp/IProjMatrix.hpp
%include LinearOp/ProjMatrix.hpp
%include LinearOp/PrecisionOpMultiConditional.hpp
%include LinearOp/IOptimCost.hpp
%include LinearOp/OptimCostBinary.hpp
%include LinearOp/OptimCostColored.hpp
%include LinearOp/EPowerPT.hpp

%include Model/ANoStat.hpp
%include Model/NoStatArray.hpp
%include Model/NoStatFunctional.hpp

%include Neigh/ANeighParam.hpp
%include Neigh/NeighUnique.hpp
%include Neigh/NeighImage.hpp
%include Neigh/NeighMoving.hpp
%include Neigh/NeighBench.hpp
%include Neigh/ENeigh.hpp
%include Neigh/NeighWork.hpp

%include Variogram/VarioParam.hpp
%include Variogram/Vario.hpp
%include Variogram/VarioParam.hpp
%include Variogram/DirParam.hpp
%include Variogram/ECalcVario.hpp

%include Model/Model.hpp
%include Model/Option_AutoFit.hpp
%include Model/Option_VarioFit.hpp
%include Model/Constraints.hpp
%include Model/ConsItem.hpp
%include Model/CovParamId.hpp
%include Model/EModelProperty.hpp
%include Model/EConsElem.hpp
%include Model/EConsType.hpp
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
%include Covariances/CovP8.hpp
%include Covariances/CovPenta.hpp
%include Covariances/CovPower.hpp
%include Covariances/CovReg1D.hpp
%include Covariances/CovSincard.hpp
%include Covariances/CovSpherical.hpp
%include Covariances/CovStable.hpp
%include Covariances/CovStorkey.hpp
%include Covariances/CovTriangle.hpp
%include Covariances/CovWendland1.hpp
%include Covariances/CovWendland2.hpp
%include Covariances/ECov.hpp
%include Covariances/ETape.hpp
%include Covariances/EConvType.hpp
%include Covariances/EConvDir.hpp

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
%include Drifts/EDrift.hpp

%include Matrix/AMatrix.hpp
%include Matrix/AMatrixSquare.hpp
%include Matrix/MatrixRectangular.hpp
%include Matrix/MatrixSquareDiagonal.hpp
%include Matrix/MatrixSquareDiagonalCst.hpp
%include Matrix/MatrixSquareGeneral.hpp
%include Matrix/MatrixSquareSymmetric.hpp

%include API/SPDE.hpp
%include API/PGSSPDE.hpp
%include API/ESPDECalcMode.hpp

%include Db/Db.hpp
%include Db/DbGrid.hpp
%include Db/DbStringFormat.hpp
%include Db/ELoadBy.hpp
%include Db/ELoc.hpp

%include Anamorphosis/AAnam.hpp
%include Anamorphosis/AnamContinuous.hpp
%include Anamorphosis/AnamDiscrete.hpp
%include Anamorphosis/AnamUser.hpp
%include Anamorphosis/AnamHermite.hpp
%include Anamorphosis/AnamEmpirical.hpp
%include Anamorphosis/AnamDiscreteDD.hpp
%include Anamorphosis/AnamDiscreteIR.hpp
%include Anamorphosis/EAnam.hpp

%include Gibbs/AGibbs.hpp
%include Gibbs/GibbsMulti.hpp
%include Gibbs/GibbsMMulti.hpp
%include Gibbs/GibbsUMulti.hpp

%include Morpho/Morpho.hpp
%include Polygon/Polygons.hpp
%include Polygon/PolySet.hpp

%include Stats/Classical.hpp
%include Stats/PCA.hpp
%include Stats/PCAStringFormat.hpp
%include Stats/Selectivity.hpp

%include LithoRule/Rule.hpp
%include LithoRule/RuleStringFormat.hpp
%include LithoRule/RuleProp.hpp
%include LithoRule/ERule.hpp

%include Estimation/KrigingSystem.hpp

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

%include Simulation/ASimulation.hpp
%include Simulation/TurningBands.hpp
%include Simulation/TurningDirection.hpp

%include segy.h

/*
// Definition of AVariableTemplate for useful type
%template(AVariableInt) AVariableTemplate<int>;
%template(AVariableDouble) AVariableTemplate<double>;
%template(AVariableBool) AVariableTemplate<bool>;
%template(AVariableString) AVariableTemplate<String>;

// THEN include our class that inherited AVariableTemplate<T>
%include Interfaces/VariableInt.hpp
%include Interfaces/VariableDouble.hpp
%include Interfaces/VariableBool.hpp
%include Interfaces/VariableString.hpp
*/

/// https://blog.mbedded.ninja/programming/languages/python/python-swig-bindings-from-cplusplus/
%feature("director");

// For suppressing SWIG warning for overloaded methods
#pragma SWIG nowarn=509
