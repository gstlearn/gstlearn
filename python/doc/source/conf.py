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

    functions_to_skip = ['sequence', 'resetSamplingDb', 'initVVDouble', 'concatenateInPlace', 'truncateDigitsInPlace', 'morpho_dilation', 'transformVD', 'normalizeFromGaussianDistribution', 'morpho_opening', 'getExtremas', 'getAllColumns', 'whereMaximum', 'createSamplingDb', 'morpho_double2image', 'gridcell_neigh', 'getName', 'identifyNames', 'simulateUniform', 'getColumnsAsVVD', 'resetFromOnePoint', 'innerProduct', 'getVariance', 'getNamesByColIdx', 'filter', 'sample', 'setArgInt', 'multiplyComplexInPlace', 'addInPlace', 'getColumnsByUIDInterval', 'orderRanks', 'createFillRandom', 'getColumns', 'createReduce', 'crossProduct3DInPlace', 'extractInPlace', 'whereMinimum', 'getAllNames', 'squeezeAndStretchInPlaceForward', 'resetFromBox', 'setLocVariable', 'createFromBox', 'getNamesByLocator', 'unflattenInPlace', 'flatten', 'getColumnsAsMatrix', 'add', 'sampleRanks', 'morpho_image2double', 'createFromDbGrid', 'getMean', 'cumulLog', 'morpho_negation', 'isInList', 'getRange', 'simulateInteger', 'isEqual', 'deleteColumnsByLocator', 'initVVInt', 'getCoorMinimum', 'deleteColumn', 'create', 'updSimvar', 'getNamesByUID', 'cumulateInPlace', 'squeezeAndStretchInPlaceBackward', 'initThread', 'qnormVec', 'reduceOne', 'isDimensionIndexValid', 'dumpRange', 'deleteColumnsByUID', 'isUIDValid', 'cumulIncrement', 'deleteColumnsByColIdx', 'extensionDiagonal', 'isSampleIndexValid', 'simulateBernoulli', 'whereElement', 'getColumnsActiveAndDefined', 'getExtensionDiagonal', 'getMostSignificant', 'setSimvar', 'addMultiplyConstantInPlace', 'reorder', 'getColumnsByUID', 'concatenate', 'getColumnByUID', 'cumsum', 'morpho_union', 'mean1AndMean2ToStdev', 'sortRanks', 'setName', 'createFromOnePoint', 'capInPlace', 'simulateBoolean', 'quantiles', 'copy', 'setArgVectorDouble', 'getCenters', 'getArgVectorDouble', 'setArgVectorInt', 'revert', 'reduce', 'createFromNF', 'capInPlaceVVD', 'createFromGridRandomized', 'initVDouble', 'truncateDecimalsInPlace', 'inverse', 'crossProduct3D', 'addCst', 'isLocatorIndexValid', 'getColumnsByLocator', 'getArgVectorInt', 'isSampleIndicesValid', 'unflatten', 'getExtensionInPlace', 'isIsotropic', 'compress', 'deleteColumnsByUIDRange', 'getNameByLocator', 'fillUndef', 'rangeVals', 'getCenter', 'resetReduce', 'morpho_closing', 'cumulate', 'isSorted', 'resetFromGridRandomized', 'flattenInPlace', 'createFromSamples', 'getSimvar', 'normalScore', 'getColumnByLocator', 'getExtends', 'initVString', 'dumpStats', 'linearCombinationVVDInPlace', 'getExtrema', 'getNames', 'getLocVariable', 'morpho_duplicate', 'initVInt', 'getColumnsByColIdxInterval', 'isEqualExtended', 'getColumnsByColIdx', 'sort', 'createFromCSV', 'count', 'resetFromSamples', 'createEmpty', 'power', 'deleteColumns', 'deleteColumnByUID', 'morpho_labelsize', 'arrangeInPlace', 'dumpNNZ', 'pnormVec', 'getColumn', 'normalizeCodir', 'getNameByUID', 'hasLocVariable', 'getArgInt', 'getCorrelation', 'simulateGaussian', 'expandNameList', 'normalize', 'getNameByColIdx', 'getMaximum', 'sequenceInPlace', 'addSquareInPlace', 'mergeInPlace', 'updLocVariable', 'simulateGaussianInPlace', 'getCoorMaximum', 'sequenceVD', 'resetFromCSV', 'getColumnByColIdx', 'deleteColumnByColIdx', 'getExtension', 'suppressTest', 'setNameByColIdx', 'unique', 'setNameByUID', 'getStdv', 'linearCombinationInPlace', 'getMinimum', 'complement', 'setNameByLocator', 'sortInPlace']
    if name in functions_to_skip:
        return True
    return skip

def setup(app):
    app.connect('autodoc-skip-member', autodoc_skip_member_handler)
