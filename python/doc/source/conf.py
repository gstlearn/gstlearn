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

    functions_to_skip = ['add', 'addCst', 'addInPlace', 'addMultiplyConstantInPlace', 'addSquareInPlace', 'arrangeInPlace', 'capInPlace', 'capInPlaceVVD', 'complement', 'compress', 'concatenate', 'concatenateInPlace', 'copy', 'count', 'create', 'createEmpty', 'createFillRandom', 'createFromBox', 'createFromCSV', 'createFromDbGrid', 'createFromGridRandomized', 'createFromNF', 'createFromOnePoint', 'createFromSamples', 'createReduce', 'createSamplingDb', 'crossProduct3D', 'crossProduct3DInPlace', 'cumsum', 'cumulIncrement', 'cumulLog', 'cumulate', 'cumulateInPlace', 'deleteColumn', 'deleteColumnByColIdx', 'deleteColumnByUID', 'deleteColumns', 'deleteColumnsByColIdx', 'deleteColumnsByLocator', 'deleteColumnsByUID', 'deleteColumnsByUIDRange', 'dumpNNZ', 'dumpRange', 'dumpStats', 'expandNameList', 'extensionDiagonal', 'extractInPlace', 'fillUndef', 'filter', 'flatten', 'flattenInPlace', 'getAllColumns', 'getAllNames', 'getArgInt', 'getArgVectorDouble', 'getArgVectorInt', 'getCenter', 'getCenters', 'getColumn', 'getColumnByColIdx', 'getColumnByLocator', 'getColumnByUID', 'getColumns', 'getColumnsActiveAndDefined', 'getColumnsAsMatrix', 'getColumnsAsVVD', 'getColumnsByColIdx', 'getColumnsByColIdxInterval', 'getColumnsByLocator', 'getColumnsByUID', 'getColumnsByUIDInterval', 'getCoorMaximum', 'getCoorMinimum', 'getCorrelation', 'getExtends', 'getExtension', 'getExtensionDiagonal', 'getExtensionInPlace', 'getExtrema', 'getExtremas', 'getLocVariable', 'getMaximum', 'getMean', 'getMinimum', 'getMostSignificant', 'getName', 'getNameByColIdx', 'getNameByLocator', 'getNameByUID', 'getNames', 'getNamesByColIdx', 'getNamesByLocator', 'getNamesByUID', 'getRange', 'getSimvar', 'getStdv', 'getVariance', 'gridcell_neigh', 'hasLocVariable', 'identifyNames', 'initThread', 'initVDouble', 'initVInt', 'initVString', 'initVVDouble', 'initVVInt', 'innerProduct', 'inverse', 'isDimensionIndexValid', 'isEqual', 'isEqualExtended', 'isInList', 'isIsotropic', 'isLocatorIndexValid', 'isSampleIndexValid', 'isSampleIndicesValid', 'isSorted', 'isUIDValid', 'linearCombinationInPlace', 'linearCombinationVVDInPlace', 'mean1AndMean2ToStdev', 'mergeInPlace', 'morpho_closing', 'morpho_dilation', 'morpho_double2image', 'morpho_duplicate', 'morpho_image2double', 'morpho_labelsize', 'morpho_negation', 'morpho_opening', 'morpho_union', 'multiplyComplexInPlace', 'normalScore', 'normalize', 'normalizeCodir', 'normalizeFromGaussianDistribution', 'orderRanks', 'pnormVec', 'power', 'qnormVec', 'quantiles', 'rangeVals', 'reduce', 'reduceOne', 'reorder', 'resetFromBox', 'resetFromCSV', 'resetFromGridRandomized', 'resetFromOnePoint', 'resetFromSamples', 'resetReduce', 'resetSamplingDb', 'revert', 'sample', 'sampleRanks', 'sequence', 'sequenceInPlace', 'sequenceVD', 'setArgInt', 'setArgVectorDouble', 'setArgVectorInt', 'setLocVariable', 'setName', 'setNameByColIdx', 'setNameByLocator', 'setNameByUID', 'setSimvar', 'simulateBernoulli', 'simulateBoolean', 'simulateGaussian', 'simulateGaussianInPlace', 'simulateInteger', 'simulateUniform', 'sort', 'sortInPlace', 'sortRanks', 'squeezeAndStretchInPlaceBackward', 'squeezeAndStretchInPlaceForward', 'suppressTest', 'transformVD', 'truncateDecimalsInPlace', 'truncateDigitsInPlace', 'unflatten', 'unflattenInPlace', 'unique', 'updLocVariable', 'updSimvar', 'whereElement', 'whereMaximum', 'whereMinimum']
    if name in functions_to_skip:
        return True
    return skip

def setup(app):
    app.connect('autodoc-skip-member', autodoc_skip_member_handler)
