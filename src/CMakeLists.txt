# Why is cmake file GLOB evil?
# https://stackoverflow.com/questions/32411963/why-is-cmake-file-glob-evil
set(SOURCES
  License/MD5Utility.cpp
  License/RegistryUtility.cpp
  License/MACAddressUtility.cpp
  License/TimeUtility.cpp
  License/LicenseUtility.cpp
  License/LicenseKey.cpp
  Matrix/MatrixRectangular.cpp
  Matrix/AMatrix.cpp
  Matrix/MatrixSquareSymmetric.cpp
  Matrix/MatrixFactory.cpp
  Matrix/MatrixSquareDiagonalCst.cpp
  Matrix/MatrixSquareDiagonal.cpp
  Matrix/AMatrixSquare.cpp
  Matrix/MatrixSquareGeneral.cpp
  API/PGSSPDE.cpp
  API/SPDE.cpp
  INIParser.cpp
  Gibbs/GibbsUPropMono.cpp
  Gibbs/GibbsFactory.cpp
  Gibbs/GibbsMultiMono.cpp
  Gibbs/GibbsMulti.cpp
  Gibbs/GibbsMMulti.cpp
  Gibbs/GibbsUMultiMono.cpp
  Gibbs/AGibbs.cpp
  Gibbs/GibbsUMulti.cpp
  LithoRule/RuleProp.cpp
  LithoRule/Rule.cpp
  LithoRule/RuleShadow.cpp
  LithoRule/RuleShift.cpp
  LithoRule/Node.cpp
  Enum/Enums.cpp
  Model/ConsItem.cpp
  Model/NoStatArray.cpp
  Model/ANoStat.cpp
  Model/Constraints.cpp
  Model/ElemNoStat.cpp
  Model/CovInternal.cpp
  Model/Option_AutoFit.cpp
  Model/Tapering.cpp
  Model/ModelNostat.cpp
  Model/Option_VarioFit.cpp
  Model/ModTrans.cpp
  Model/Model.cpp
  Model/NoStatFunctional.cpp
  Model/Convolution.cpp
  Covariances/CovGaussian.cpp
  Covariances/CovLinear.cpp
  Covariances/CovWendland2.cpp
  Covariances/CovGradientFunctional.cpp
  Covariances/CovBesselK.cpp
  Covariances/CovLMGradient.cpp
  Covariances/CovPower.cpp
  Covariances/CovStorkey.cpp
  Covariances/CovGC5.cpp
  Covariances/ACovAnisoList.cpp
  Covariances/CovPenta.cpp
  Covariances/CovGamma.cpp
  Covariances/CovExponential.cpp
  Covariances/CovGC1.cpp
  Covariances/ACov.cpp
  Covariances/CovSincard.cpp
  Covariances/CovLMC.cpp
  Covariances/ACovGradient.cpp
  Covariances/CovGCspline2.cpp
  Covariances/CovGC3.cpp
  Covariances/CovGCspline.cpp
  Covariances/CovTriangle.cpp
  Covariances/CovCalcMode.cpp
  Covariances/CovP8.cpp
  Covariances/CovContext.cpp
  Covariances/CovCosinus.cpp
  Covariances/CovCosExp.cpp
  Covariances/CovReg1D.cpp
  Covariances/ACovFunc.cpp
  Covariances/CovAniso.cpp
  Covariances/CovGradientNumerical.cpp
  Covariances/CovCubic.cpp
  Covariances/CovNugget.cpp
  Covariances/CovWendland1.cpp
  Covariances/CovBesselJ.cpp
  Covariances/CovSpherical.cpp
  Covariances/CovCauchy.cpp
  Covariances/CovFactory.cpp
  Covariances/CovStable.cpp
  Polygon/Polygons.cpp
  Polygon/PolySet.cpp
  Core/surface.cpp
  Core/simbool.cpp
  Core/math.cpp
  Core/stats.cpp
  Core/license.cpp
  Core/spatial.cpp
  Core/csparse.cpp
  Core/sphtriangle.cpp
  Core/vario.cpp
  Core/simsub.cpp
  Core/simpart.cpp
  Core/io.cpp
  Core/db.cpp
  Core/variopgs.cpp
  Core/model_auto.cpp
  Core/cluster.cpp
  Core/fracture.cpp
  Core/simtub.cpp
  Core/matrix.cpp
  Core/spill.cpp
  Core/simsphe.cpp
  Core/acknowledge.cpp
  Core/simfine.cpp
  Core/vtk.cpp
  Core/simfft.cpp
  Core/geophy.cpp
  Core/potential.cpp
  Core/foxleg.cpp
  Core/eden.cpp
  Core/util.cpp
  Core/anam.cpp
  Core/ascii.cpp
  Core/skin.cpp
  Core/pile.cpp
  Core/segy.cpp
  Core/mlayers.cpp
  Core/thresh.cpp
  Core/neigh.cpp
  Core/fft.cpp
  Core/convert.cpp
  Core/memory.cpp
  Core/poly.cpp
  Core/model.cpp
  Core/krige.cpp
  Core/seismic.cpp
  Core/spde.cpp
  Core/dbtools.cpp
  Anamorphosis/AnamContinuous.cpp
  Anamorphosis/AnamDiscrete.cpp
  Anamorphosis/AnamDiscreteDD.cpp
  Anamorphosis/AnamUser.cpp
  Anamorphosis/Anam.cpp
  Anamorphosis/AnamHermite.cpp
  Anamorphosis/AnamDiscreteIR.cpp
  Anamorphosis/AnamEmpirical.cpp
  Interfaces/VariableCategorical.cpp
  Interfaces/VariableInt.cpp
  Interfaces/VariableString.cpp
  Interfaces/Category.cpp
  Interfaces/ParamGrid.cpp
  Interfaces/function.cpp
  Interfaces/ParamCSV.cpp
  Interfaces/Database.cpp
  Interfaces/Dictionary.cpp
  Interfaces/VariableBool.cpp
  Interfaces/AVariable.cpp
  Interfaces/VariableDouble.cpp
  Db/Db.cpp
  Db/PtrGeos.cpp
  LinearOp/PrecisionOp.cpp
  LinearOp/TurboOptimizer.cpp
  LinearOp/ProjMatrix.cpp
  LinearOp/PrecisionOpCs.cpp
  LinearOp/ALinearOpMulti.cpp
  LinearOp/ALinearOp.cpp
  LinearOp/Identity.cpp
  LinearOp/ProdMatVect.cpp
  LinearOp/ShiftOpCs.cpp
  LinearOp/HessianOp.cpp
  LinearOp/OptimCostBinary.cpp
  LinearOp/PrecisionOpMultiConditional.cpp
  LinearOp/OptimCostColored.cpp
  Space/SpaceSN.cpp
  Space/SpaceRN.cpp
  Space/SpacePoint.cpp
  Space/ASpaceObject.cpp
  Space/ASpace.cpp
  Variogram/Vario.cpp
  Variogram/VarioParam.cpp
  Variogram/DirParam.cpp
  Basic/Limits.cpp
  Basic/AStringable.cpp
  Basic/GlobalEnvironment.cpp
  Basic/MathFunc.cpp
  Basic/Utilities.cpp
  Basic/ArgumentTest.cpp
  Basic/Table.cpp
  Basic/File.cpp
  Basic/Vector.cpp
  Basic/FunctionalSpirale.cpp
  Basic/Law.cpp
  Basic/Interval.cpp
  Basic/Tensor.cpp
  Basic/Rotation.cpp
  Basic/AException.cpp
  Basic/ASerializable.cpp
  Basic/String.cpp
  Basic/GridC.cpp
  Basic/AFunctional.cpp
  Basic/CSVformat.cpp
  Basic/NamingConvention.cpp
  Polynomials/Chebychev.cpp
  Polynomials/APolynomial.cpp
  Polynomials/Hermite.cpp
  Polynomials/ClassicalPolynomial.cpp
  Polynomials/MonteCarlo.cpp
  Mesh/delaunay.cpp
  Mesh/MeshETurbo.cpp
  Mesh/predicates.cpp
  Mesh/m_delaunay.cpp
  Mesh/MeshSpherical.cpp
  Mesh/MeshFactory.cpp
  Mesh/MeshEStandard.cpp
  Mesh/AMesh.cpp
  Mesh/tetgen.cpp
  Morpho/Morpho.cpp
  Stats/Classical.cpp
  Stats/PCA.cpp
  Drifts/DriftX.cpp
  Drifts/DriftY.cpp
  Drifts/DriftX2Y.cpp
  Drifts/DriftZ3.cpp
  Drifts/Drift1.cpp
  Drifts/DriftXY.cpp
  Drifts/DriftF.cpp
  Drifts/DriftFactory.cpp
  Drifts/DriftX3.cpp
  Drifts/DriftZ.cpp
  Drifts/DriftX2.cpp
  Drifts/DriftXZ.cpp
  Drifts/DriftY2.cpp
  Drifts/DriftY3.cpp
  Drifts/DriftXY2.cpp
  Drifts/ADrift.cpp
  Drifts/DriftZ2.cpp
  Drifts/ADriftElem.cpp
  Drifts/DriftYZ.cpp
  Drifts/ADriftList.cpp
  Neigh/Neigh.cpp
)
#option(Boost_DEBUG "Set to ON to enable debug output from FindBoost." ON)
#option(Boost_DETAILED_FAILURE_MSG "Set to ON to get detailed information" ON)

# Not for Linux?
find_package(Boost REQUIRED)

# Compile objects only once
#https://stackoverflow.com/a/29824424/3952924
add_library(objlib OBJECT ${SOURCES})

# Shared libraries need -fPIC (will be used also for static)
set_property(TARGET objlib PROPERTY POSITION_INDEPENDENT_CODE 1)

# On Windows, prevent the include sys/time.h (-DNO_TIMER)
target_compile_definitions(objlib PRIVATE NO_TIMER) 

# Include directories (-I)
#target_include_directories(objlib PUBLIC "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>")
target_include_directories(objlib PUBLIC
    "$<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/gstlearn>"
    "$<INSTALL_INTERFACE:${CMAKE_INSTALL_DATAROOTDIR}/gstlearn/swig>"
)
target_include_directories(objlib PUBLIC ${Boost_INCLUDE_DIRS})

# Define objects list for each target library
add_library(gstlearn-shared SHARED $<TARGET_OBJECTS:objlib>)
#target_link_libraries(gstlearn-shared PUBLIC objlib)
add_library(gstlearn-static STATIC $<TARGET_OBJECTS:objlib>)
#target_link_libraries(gstlearn-static PUBLIC objlib)

# Rename the output library names
set_target_properties(gstlearn-shared PROPERTIES OUTPUT_NAME gstlearn) # Will produce libgstlearn.so (on Linux)
set_target_properties(gstlearn-static PROPERTIES OUTPUT_NAME gstlearn) # Will produce libgstlearn.a  (on Linux)

# Rename the export names (used by cmake for importing them like gstlearn::shared)
set_target_properties(gstlearn-shared PROPERTIES EXPORT_NAME shared)
set_target_properties(gstlearn-static PROPERTIES EXPORT_NAME static)

# Link to boost
target_link_libraries(gstlearn-shared PUBLIC Boost::boost)
target_link_libraries(gstlearn-static PUBLIC Boost::boost)

# Install the libraries
# https://stackoverflow.com/questions/14990343/cmake-error-targets-given-no-library-destination-for-shared-library-target
install(TARGETS gstlearn-shared gstlearn-static EXPORT corelibs  # Targets (libraries)
        ARCHIVE DESTINATION lib                                  # Static  library in DESTINATION/lib folder
        LIBRARY DESTINATION lib                                  # Dynamic library in DESTINATION/lib folder
)
# Install the includes  
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include/                   # Includes from source folder
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/gstlearn         # in DESTINATION/include/gstlearn
)
# Install swig interface files
install(DIRECTORY ${CMAKE_SOURCE_DIR}/swig/                      # Swig files from source folder
        DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/gstlearn/swig   # in DESTINATION/share/gstlearn/swig
)
