from os import name

project = "gstlearn"

copyright = ""
author = ""

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.mathjax",
    "sphinx.ext.autosummary",
    "myst_parser",
]

html_theme = "sphinx_rtd_theme"

autodoc_default_options = {
    "members": True,
    "undoc-members": True,
    "inherited-members": True,
}

napoleon_google_docstring = True
napoleon_numpy_docstring = True
napoleon_use_ivar = True
napoleon_use_param = False
autosummary_generate = True

suppress_warnings = ["autodoc", "docutils"]


# fmt: off

def autodoc_skip_member_handler(app, what, name, obj, skip, options):
    short_name = name.split(".")[-1]

    if short_name == "thisown":
        return True

    if short_name.startswith("E_"):
        return True

    functions_to_skip = ['_combineSelection', '_deserializeAscii', '_serializeAscii', 'add', 'addColumns', 'addColumnsByConstant', 'addColumnsByVVD', 'addColumnsRandom', 'addCst', 'addInPlace', 'addMultiplyConstantInPlace', 'addSamples', 'addSelection', 'addSelectionByLimit', 'addSelectionByRanks', 'addSelectionByVariable', 'addSelectionFromDbByConvexHull', 'addSelectionRandom', 'addSquareInPlace', 'areSame', 'arrangeInPlace', 'capInPlace', 'capInPlaceVVD', 'clearLocator', 'clearLocators', 'clearSelection', 'complement', 'compress', 'concatenate', 'concatenateInPlace', 'copy', 'copyByCol', 'copyByName', 'copyByUID', 'count', 'coverSeveralDbs', 'create', 'createEmpty', 'createFillRandom', 'createFromBox', 'createFromCSV', 'createFromDbGrid', 'createFromGridRandomized', 'createFromNF', 'createFromOnePoint', 'createFromSamples', 'createReduce', 'createSamplingDb', 'crossProduct3D', 'crossProduct3DInPlace', 'cumsum', 'cumulIncrement', 'cumulLog', 'cumulate', 'cumulateInPlace', 'deleteColumn', 'deleteColumnByColIdx', 'deleteColumnByUID', 'deleteColumns', 'deleteColumnsByColIdx', 'deleteColumnsByLocator', 'deleteColumnsByUID', 'deleteColumnsByUIDRange', 'deleteSample', 'deleteSamples', 'deleteVersion', 'dumpGeometry', 'dumpNNZ', 'dumpRange', 'dumpStats', 'duplicateColumnByUID', 'expandNameList', 'extensionDiagonal', 'extractAllVersions', 'extractInPlace', 'extractVersion', 'fillUndef', 'filter', 'flatten', 'flattenInPlace', 'fromUID', 'generateRank', 'getActiveAndDefinedArray', 'getActiveArray', 'getAllColumns', 'getAllCoordinates', 'getAllCoordinatesMat', 'getAllNames', 'getArgInt', 'getArgVectorDouble', 'getArgVectorInt', 'getArrayBySample', 'getArrayVec', 'getArrayVecByUID', 'getCenter', 'getCenters', 'getCodeList', 'getColIdx', 'getColIdxByUID', 'getColIdxs', 'getColIdxsByLocator', 'getColIdxsByUID', 'getColumn', 'getColumnByColIdx', 'getColumnByLocator', 'getColumnByUID', 'getColumnPtr', 'getColumns', 'getColumnsActiveAndDefined', 'getColumnsAsMatrix', 'getColumnsAsVVD', 'getColumnsByColIdx', 'getColumnsByColIdxInterval', 'getColumnsByLocator', 'getColumnsByUID', 'getColumnsByUIDInterval', 'getContentsAsTable', 'getCoorMaximum', 'getCoorMinimum', 'getCorrelation', 'getCorrelationAsTable', 'getCosineToDirection', 'getData', 'getDistance', 'getDistance1D', 'getDistanceVecInPlace', 'getExtends', 'getExtension', 'getExtensionDiagonal', 'getExtensionInPlace', 'getExtrema', 'getExtremas', 'getGradient', 'getIncrement1D', 'getIncrements', 'getItem', 'getItemNames', 'getLastColumn', 'getLastName', 'getLastUID', 'getListOfSampleIndicesInPlace', 'getLocVariables', 'getLocator', 'getLocatorByUID', 'getLocators', 'getMaximum', 'getMean', 'getMinimum', 'getMostSignificant', 'getMultipleSelectedRanks', 'getMultipleSelectedVariables', 'getNBounds', 'getNColumn', 'getNFacies', 'getNInterval', 'getNLoc', 'getNOccurence', 'getNSample', 'getNSampleActive', 'getNSampleActiveAndDefined', 'getNVersions', 'getNVersionsByUID', 'getName', 'getNameByUID', 'getNames', 'getNamesByColIdx', 'getNamesByLocator', 'getNamesByUID', 'getRange', 'getRankAbsoluteToRelative', 'getRankAbsoluteToRelativeVec', 'getRankRelativeToAbsolute', 'getRankRelativeToAbsoluteVec', 'getSampleAsSPInPlace', 'getSampleAsSTInPlace', 'getSampleClosestTo', 'getSampleCoordinates', 'getSampleLocators', 'getSampleRanks', 'getSampleRanksInPlace', 'getSampleRanksPerVariable', 'getSamplesAsSP', 'getSamplesFromNbghAsSP', 'getSamplesOneCoordinate', 'getSelection', 'getSelections', 'getSimvar', 'getSortArray', 'getStatsAsTable', 'getStatsByCategoryAsTable', 'getStdv', 'getSumNVersions', 'getTangent', 'getUID', 'getUIDs', 'getUIDsByColIdx', 'getUIDsByLocator', 'getUIDsDefined', 'getValue', 'getValuesByColIdx', 'getValuesByNames', 'getValuesByRanks', 'getValuesByRanksInPlace', 'getVariance', 'getVecCoordinate', 'getWeight', 'getWeights', 'getWithinBounds', 'gridcell_neigh', 'hasLargerDimension', 'hasLocator', 'hasLocatorDefined', 'hasSameDimension', 'identifyNames', 'initThread', 'initVDouble', 'initVInt', 'initVString', 'initVVDouble', 'initVVInt', 'innerProduct', 'inverse', 'isActive', 'isActiveAndDefined', 'isActiveDomain', 'isAllIsotopic', 'isAllUndefined', 'isAllUndefinedByType', 'isColIdxValid', 'isDimensionIndexValid', 'isEqual', 'isEqualExtended', 'isInList', 'isIsotopic', 'isIsotropic', 'isLocatorIndexValid', 'isNVarComparedTo', 'isOldStyle', 'isSampleIndexValid', 'isSampleIndicesValid', 'isSorted', 'linearCombinationInPlace', 'linearCombinationVVDInPlace', 'mean1AndMean2ToStdev', 'mergeInPlace', 'morpho_closing', 'morpho_dilation', 'morpho_double2image', 'morpho_duplicate', 'morpho_image2double', 'morpho_labelsize', 'morpho_negation', 'morpho_opening', 'morpho_union', 'multiplyComplexInPlace', 'normalScore', 'normalize', 'normalizeCodir', 'normalizeFromGaussianDistribution', 'orderRanks', 'pnormVec', 'power', 'printOneSample', 'qnormVec', 'quantiles', 'rangeVals', 'reduce', 'reduceOne', 'reorder', 'resetFromBox', 'resetFromCSV', 'resetFromGridRandomized', 'resetFromOnePoint', 'resetFromSamples', 'resetReduce', 'resetSamplingDb', 'resizeSamples', 'revert', 'sample', 'sampleRanks', 'sequence', 'sequenceInPlace', 'sequenceVD', 'setAllColumns', 'setArgInt', 'setArgVectorDouble', 'setArgVectorInt', 'setArrayBySample', 'setArrayVec', 'setArraysByUID', 'setBound', 'setColumn', 'setColumnByColIdx', 'setColumnByColIdxOldStyle', 'setColumnByUID', 'setColumnByUIDOldStyle', 'setColumnsByColIdx', 'setCoordinate', 'setCoordinates', 'setInterval', 'setItem', 'setLocVariables', 'setLocator', 'setLocatorByUID', 'setLocators', 'setLocatorsByColIdx', 'setLocatorsByUID', 'setName', 'setNameByUID', 'setNames', 'setOldStyle', 'setSampleCoordinates', 'setSimvar', 'setValudByUIDNext', 'setValue', 'setValueByUID', 'setValuesByColIdxInPlace', 'setValuesByNamesInPlace', 'shrinkToValidCols', 'shrinkToValidRows', 'simulateBernoulli', 'simulateBoolean', 'simulateGaussian', 'simulateGaussianInPlace', 'simulateInteger', 'simulateUniform', 'sort', 'sortInPlace', 'sortRanks', 'squeezeAndStretchInPlaceBackward', 'squeezeAndStretchInPlaceForward', 'suppressTest', 'switchLocator', 'transformVD', 'truncateDecimalsInPlace', 'truncateDigitsInPlace', 'unflatten', 'unflattenInPlace', 'unique', 'updSimvar', 'updValue', 'whereElement', 'whereMaximum', 'whereMinimum']
    if name in functions_to_skip:
        return True
    return skip

def setup(app):
    app.connect('autodoc-skip-member', autodoc_skip_member_handler)
