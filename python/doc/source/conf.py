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

    functions_to_skip = ['createFromGridRandomized', 'pnormVec', 'simulateBoolean', 'dumpNNZ', 'getNameByColIdx', 'createFromOnePoint', 'deleteColumn', 'isEqualExtended', 'unflattenInPlace', 'createSamplingDb', 'sequence', 'getColumnsActiveAndDefined', 'getExtension', 'getColumnByUID', 'getColumnsByUIDInterval', 'sampleRanks', 'concatenateInPlace', 'initVInt', 'count', 'crossProduct3D', 'resetFromOnePoint', 'resetFromBox', 'cumulLog', 'morpho_closing', 'reorder', 'quantiles', 'simulateGaussianInPlace', 'morpho_opening', 'simulateInteger', 'createFromSamples', 'getArgInt', 'truncateDigitsInPlace', 'getAllColumns', 'unflatten', 'getMean', 'addSquareInPlace', 'createFromBox', 'reduceOne', 'getNameByLocator', 'sort', 'getLocVariable', 'setNameByUID', 'initVVDouble', 'addCst', 'sequenceInPlace', 'normalScore', 'revert', 'getStdv', 'getCoorMinimum', 'createFromCSV', 'getNames', 'getNamesByColIdx', 'transformVD', 'initVVInt', 'deleteColumnsByLocator', 'sequenceVD', 'cumulIncrement', 'morpho_labelsize', 'morpho_image2double', 'createFillRandom', 'getCoorMaximum', 'cumsum', 'isIsotropic', 'concatenate', 'getCenters', 'gridcell_neigh', 'getExtremas', 'resetSamplingDb', 'orderRanks', 'getExtends', 'initVString', 'resetFromCSV', 'squeezeAndStretchInPlaceForward', 'setSimvar', 'getExtensionDiagonal', 'create', 'capInPlaceVVD', 'getColumn', 'getNamesByLocator', 'getName', 'deleteColumnsByUID', 'whereMaximum', 'identifyNames', 'suppressTest', 'getRange', 'addMultiplyConstantInPlace', 'isInList', 'crossProduct3DInPlace', 'reduce', 'getArgVectorDouble', 'isDimensionIndexValid', 'sample', 'getArgVectorInt', 'setArgVectorDouble', 'arrangeInPlace', 'complement', 'setNameByColIdx', 'getVariance', 'add', 'createEmpty', 'dumpStats', 'mean1AndMean2ToStdev', 'normalizeCodir', 'morpho_double2image', 'morpho_duplicate', 'getColumnsByColIdxInterval', 'inverse', 'linearCombinationInPlace', 'fillUndef', 'getColumnsByColIdx', 'initVDouble', 'normalize', 'setLocVariable', 'isSorted', 'getExtrema', 'getNameByUID', 'createFromNF', 'deleteColumnsByColIdx', 'isLocatorIndexValid', 'rangeVals', 'capInPlace', 'compress', 'simulateUniform', 'setArgVectorInt', 'extractInPlace', 'getColumnByLocator', 'deleteColumns', 'qnormVec', 'simulateBernoulli', 'hasLocVariable', 'simulateGaussian', 'getColumnByColIdx', 'flattenInPlace', 'truncateDecimalsInPlace', 'power', 'deleteColumnsByUIDRange', 'getAllNames', 'unique', 'dumpRange', 'getColumnsByLocator', 'flatten', 'extensionDiagonal', 'getCenter', 'linearCombinationVVDInPlace', 'copy', 'whereElement', 'setArgInt', 'setName', 'getMinimum', 'getExtensionInPlace', 'getMostSignificant', 'cumulateInPlace', 'resetReduce', 'addInPlace', 'createReduce', 'deleteColumnByUID', 'isUIDValid', 'sortRanks', 'morpho_union', 'isEqual', 'resetFromSamples', 'cumulate', 'getColumnsAsMatrix', 'mergeInPlace', 'innerProduct', 'multiplyComplexInPlace', 'normalizeFromGaussianDistribution', 'setNameByLocator', 'deleteColumnByColIdx', 'morpho_negation', 'filter', 'isSampleIndicesValid', 'getNamesByUID', 'whereMinimum', 'getColumnsByUID', 'createFromDbGrid', 'resetFromGridRandomized', 'getCorrelation', 'updSimvar', 'getSimvar', 'sortInPlace', 'getColumnsAsVVD', 'squeezeAndStretchInPlaceBackward', 'updLocVariable', 'initThread', 'isSampleIndexValid', 'expandNameList', 'morpho_dilation', 'getMaximum', 'getColumns']
    if name in functions_to_skip:
        return True
    return skip

def setup(app):
    app.connect('autodoc-skip-member', autodoc_skip_member_handler)
