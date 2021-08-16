%ignore *::operator=;

class IClonable{};

%include stl.i
%include std_vector.i
%include std_string.i

// Keep order in the file
namespace std {
%template(VectorDouble)       vector< double >;
%template(VectorInt)          vector< int >;
%template(VectorString)       vector< string >;
%template(VectorBool)         vector< bool >;
%template(VectorUChar)        vector< unsigned char >;
%template(VectorVectorInt)    vector< vector< int > >;
%template(VectorVectorDouble) vector< vector< double > >;
%template(VectorEnumCovs)     vector< ENUM_COVS >;
};

%template(VectorCova) std::vector<Cova*>;
%template(VectorCTable) std::vector<CTable*>;
%template(VectorDir) std::vector<Dir>;  // Not a pointers list
%template(VectorDirection) std::vector<Direction*>;
%template(VectorDrft) std::vector<Drift*>;
%template(VectorFrac_Desc) std::vector<Frac_Desc*>;
%template(VectorFrac_Fam) std::vector<Frac_Fam*>;
%template(VectorFrac_Fault) std::vector<Frac_Fault*>;
%template(VectorLocal_Split) std::vector<Local_Split*>;
%template(VectorPolySet) std::vector<PolySet*>;
%template(VectorQChol) std::vector<QChol*>;
%template(VectorSPDE_SS_Option) std::vector<SPDE_SS_Option*>;
%template(VectorSubPlan) std::vector<SubPlan*>;
%template(VectorToken_Def) std::vector<Token_Def*>;
%template(VectorToken_Par) std::vector<Token_Par*>;
%template(VectorIntervals) std::vector<Interval*>;

%include Basic/Vector.hpp
%include csparse_d.h
%include csparse_f.h
%include geoslib_define.h
%include geoslib_enum.h
%include geoslib_d.h
%include geoslib_f.h
%include Basic/ArgumentTest.hpp
%include Basic/AStringable.hpp
%include Basic/ASerializable.hpp
%include Basic/Rotation.hpp
%include Basic/Tensor.hpp
%include Basic/GridC.hpp
%include Basic/String.hpp
%include Basic/Interval.hpp
%include Basic/Limits.hpp
%include Basic/Utilities.hpp
%include Basic/CSVformat.hpp
%include Space/Space.hpp
%include Space/ASpace.hpp
%include Space/ASpaceObject.hpp
%include Space/SpacePoint.hpp
%include Space/SpaceRN.hpp
%include Space/SpaceShape.hpp
%include Interfaces/geoslib_f_swig.h
%include Interfaces/ACalculator.hpp
%include Interfaces/AParam.hpp
%include Interfaces/AVariable.hpp
%include Interfaces/AVariableTemplate.hpp
%include Interfaces/CalcVarioExp.hpp
%include Interfaces/Category.hpp
%include Interfaces/Database.hpp
%include Interfaces/Dictionary.hpp
%include Interfaces/interface_d.hpp
%include Interfaces/ParamCSV.hpp
%include Interfaces/ParamGrid.hpp
%include Interfaces/Param.hpp
%include Interfaces/VarioDir.hpp
%include Interfaces/VarioExp.hpp
%include Interfaces/VarioValue.hpp
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
%include Model/ANoStat.hpp
%include Model/NoStatArray.hpp
%include Neigh/Neigh.hpp
%include Variogram/Vario.hpp
%include Variogram/Dir.hpp
%include Model/Model.hpp
%include Model/Cova.hpp
%include Model/Option_AutoFit.hpp
%include Model/Option_VarioFit.hpp
//%include Model/ConsItem.hpp
%include Model/Constraints.hpp
%include Model/ConsItem.hpp
%include Covariances/ACov.hpp
%include Covariances/ACovFunc.hpp
%include Covariances/ACovAnisoList.hpp
%include Covariances/CovAniso.hpp
%include Covariances/ACovGradient.hpp
%include Covariances/CovLMC.hpp
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
%include MatrixC/AMatrixC.hpp
%include MatrixC/AMatrixCSquare.hpp
%include MatrixC/MatrixCRectangular.hpp
%include MatrixC/MatrixCSDiag.hpp
%include MatrixC/MatrixCSDiagCst.hpp
%include MatrixC/MatrixCSGeneral.hpp
%include MatrixC/MatrixCSSym.hpp
%include Db/Db.hpp
%include Anamorphosis/Anam.hpp
%include Anamorphosis/AnamContinuous.hpp
%include Anamorphosis/AnamDiscrete.hpp
%include Anamorphosis/AnamUser.hpp
%include Anamorphosis/AnamHermite.hpp
%include Anamorphosis/AnamEmpirical.hpp
%include Anamorphosis/AnamDiscreteDD.hpp
%include Anamorphosis/AnamDiscreteIR.hpp
%include Morpho/Morpho.hpp
%include Polygon/Polygons.hpp
%include Polygon/PolySet.hpp
%include Stats/Classical.hpp
%include LithoRule/Rule.hpp
%include segy.h

/*Definition of AVariableTemplate for useful type*/
%template(AVariableInt) AVariableTemplate<int>;
%template(AVariableDouble) AVariableTemplate<double>;
%template(AVariableBool) AVariableTemplate<bool>;
%template(AVariableString) AVariableTemplate<String>;

/*THEN include our class that inherited AVariableTemplate<T>*/
%include Interfaces/VariableInt.hpp
%include Interfaces/VariableDouble.hpp
%include Interfaces/VariableBool.hpp
%include Interfaces/VariableString.hpp

/// https://blog.mbedded.ninja/programming/languages/python/python-swig-bindings-from-cplusplus/
%feature("director");

// For suppressing SWIG warning for overloaded methods
#pragma SWIG nowarn=509
