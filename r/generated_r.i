%insert(s)%{
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalCovMat(...))
}

assign('ModelGeneric_evalCovMat', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalCovMatInPlace(...))
}

assign('ModelGeneric_evalCovMatInPlace', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalCovMatInPlaceFromIdx(...))
}

assign('ModelGeneric_evalCovMatInPlaceFromIdx', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalCovMatSym(...))
}

assign('ModelGeneric_evalCovMatSym', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalCovMatSymInPlace(...))
}

assign('ModelGeneric_evalCovMatSymInPlace', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalCovMatSymInPlaceFromIdx(...))
}

assign('ModelGeneric_evalCovMatSymInPlaceFromIdx', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$eval0Mat(...))
}

assign('ModelGeneric_eval0Mat', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalCovMat0(...))
}

assign('ModelGeneric_evalCovMat0', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalCovMat0InPlace(...))
}

assign('ModelGeneric_evalCovMat0InPlace', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalCovVecRHSInPlace(...))
}

assign('ModelGeneric_evalCovVecRHSInPlace', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalCovMatOptimInPlace(...))
}

assign('ModelGeneric_evalCovMatOptimInPlace', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalCovMatRHSInPlaceFromIdx(...))
}

assign('ModelGeneric_evalCovMatRHSInPlaceFromIdx', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalCovMatSparse(...))
}

assign('ModelGeneric_evalCovMatSparse', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$eval0(...))
}

assign('ModelGeneric_eval0', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalCov(...))
}

assign('ModelGeneric_evalCov', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalNvarIpas(...))
}

assign('ModelGeneric_evalNvarIpas', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalNvarIpasIncr(...))
}

assign('ModelGeneric_evalNvarIpasIncr', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalIvarNlag(...))
}

assign('ModelGeneric_evalIvarNlag', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalIvarIpas(...))
}

assign('ModelGeneric_evalIvarIpas', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalCvv(...))
}

assign('ModelGeneric_evalCvv', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalCvvShift(...))
}

assign('ModelGeneric_evalCvvShift', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalCvvM(...))
}

assign('ModelGeneric_evalCvvM', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalCxv(...))
}

assign('ModelGeneric_evalCxv', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalCxvM(...))
}

assign('ModelGeneric_evalCxvM', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalPointToDb(...))
}

assign('ModelGeneric_evalPointToDb', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalPointToDbAsSP(...))
}

assign('ModelGeneric_evalPointToDbAsSP', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalAverageDbToDb(...))
}

assign('ModelGeneric_evalAverageDbToDb', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalAverageIncrToIncr(...))
}

assign('ModelGeneric_evalAverageIncrToIncr', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evalAveragePointToDb(...))
}

assign('ModelGeneric_evalAveragePointToDb', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$samplingDensityVariance(...))
}

assign('ModelGeneric_samplingDensityVariance', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$specificVolume(...))
}

assign('ModelGeneric_specificVolume', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$coefficientOfVariation(...))
}

assign('ModelGeneric_coefficientOfVariation', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$specificVolumeFromCoV(...))
}

assign('ModelGeneric_specificVolumeFromCoV', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$extensionVariance(...))
}

assign('ModelGeneric_extensionVariance', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$calculateStDev(...))
}

assign('ModelGeneric_calculateStDev', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evaluateMatInPlace(...))
}

assign('ModelGeneric_evaluateMatInPlace', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evaluateOneGeneric(...))
}

assign('ModelGeneric_evaluateOneGeneric', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$evaluateOneIncr(...))
}

assign('ModelGeneric_evaluateOneIncr', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$buildVmapOnDbGrid(...))
}

assign('ModelGeneric_buildVmapOnDbGrid', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$sample(...))
}

assign('ModelGeneric_sample', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$sampleUnitary(...))
}

assign('ModelGeneric_sampleUnitary', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$envelop(...))
}

assign('ModelGeneric_envelop', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$gofToVario(...))
}

assign('ModelGeneric_gofToVario', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$isNoStat(...))
}

assign('ModelGeneric_isNoStat', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$manage(...))
}

assign('ModelGeneric_manage', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$optimizationPreProcessForData(...))
}

assign('ModelGeneric_optimizationPreProcessForData', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$optimizationPostProcess(...))
}

assign('ModelGeneric_optimizationPostProcess', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$getDrift(...))
}

assign('ModelGeneric_getDrift', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$computeDrift(...))
}

assign('ModelGeneric_computeDrift', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$evalDriftValue(...))
}

assign('ModelGeneric_evalDriftValue', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$evalDriftMat(...))
}

assign('ModelGeneric_evalDriftMat', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$evalDriftMatInPlace(...))
}

assign('ModelGeneric_evalDriftMatInPlace', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$evalDriftMatByRanks(...))
}

assign('ModelGeneric_evalDriftMatByRanks', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$evalMeanVecByRanks(...))
}

assign('ModelGeneric_evalMeanVecByRanks', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$evalDriftMatByRanksInPlace(...))
}

assign('ModelGeneric_evalDriftMatByRanksInPlace', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$evalDriftMatByTargetInPlace(...))
}

assign('ModelGeneric_evalDriftMatByTargetInPlace', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$getNDrift(...))
}

assign('ModelGeneric_getNDrift', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$getNDriftEquation(...))
}

assign('ModelGeneric_getNDriftEquation', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$getNExtDrift(...))
}

assign('ModelGeneric_getNExtDrift', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$isFlagLinked(...))
}

assign('ModelGeneric_isFlagLinked', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$getDriftMaxIRFOrder(...))
}

assign('ModelGeneric_getDriftMaxIRFOrder', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$getRankFex(...))
}

assign('ModelGeneric_getRankFex', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$isDriftSampleDefined(...))
}

assign('ModelGeneric_isDriftSampleDefined', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$isDriftFiltered(...))
}

assign('ModelGeneric_isDriftFiltered', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$isDriftDefined(...))
}

assign('ModelGeneric_isDriftDefined', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$isDriftDifferentDefined(...))
}

assign('ModelGeneric_isDriftDifferentDefined', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$getDrifts(...))
}

assign('ModelGeneric_getDrifts', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$evalDrift(...))
}

assign('ModelGeneric_evalDrift', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$evalDriftBySample(...))
}

assign('ModelGeneric_evalDriftBySample', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$evalDriftBySampleInPlace(...))
}

assign('ModelGeneric_evalDriftBySampleInPlace', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$evalDriftCoef(...))
}

assign('ModelGeneric_evalDriftCoef', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$hasDrift(...))
}

assign('ModelGeneric_hasDrift', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$getMean(...))
}

assign('ModelGeneric_getMean', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$getMeans(...))
}

assign('ModelGeneric_getMeans', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$evalDriftVarCoef(...))
}

assign('ModelGeneric_evalDriftVarCoef', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getDriftList(self)$evalDriftVarCoefs(...))
}

assign('ModelGeneric_evalDriftVarCoefs', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getContext(self)$getNVar(...))
}

assign('ModelGeneric_getNVar', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getContext(self)$getNDim(...))
}

assign('ModelGeneric_getNDim', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getContext(self)$getSpace(...))
}

assign('ModelGeneric_getSpace', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getContext(self)$getCovar0(...))
}

assign('ModelGeneric_getCovar0', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getContext(self)$getField(...))
}

assign('ModelGeneric_getField', f , envir = asNamespace('gstlearn'))
%}
%insert(s)%{
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$getActiveFactor(...))
}

assign('Model_getActiveFactor', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$getCovAniso(...))
}

assign('Model_getCovAniso', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$getNCov(...))
}

assign('Model_getNCov', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$getCovType(...))
}

assign('Model_getCovType', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$getRange(...))
}

assign('Model_getRange', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$getRanges(...))
}

assign('Model_getRanges', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$getAngles(...))
}

assign('Model_getAngles', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$getAnam(...))
}

assign('Model_getAnam', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$getParam(...))
}

assign('Model_getParam', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$getCovName(...))
}

assign('Model_getCovName', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$extractCova(...))
}

assign('Model_extractCova', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$getNGradParam(...))
}

assign('Model_getNGradParam', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$getMaximumDistance(...))
}

assign('Model_getMaximumDistance', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$getCovMinIRFOrder(...))
}

assign('Model_getCovMinIRFOrder', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$getAnamNClass(...))
}

assign('Model_getAnamNClass', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$hasAnam(...))
}

assign('Model_hasAnam', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$hasNugget(...))
}

assign('Model_hasNugget', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$getRankNugget(...))
}

assign('Model_getRankNugget', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$getBallRadius(...))
}

assign('Model_getBallRadius', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$hasExternalCov(...))
}

assign('Model_hasExternalCov', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$isChangeSupportDefined(...))
}

assign('Model_isChangeSupportDefined', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$getAnamHermite(...))
}

assign('Model_getAnamHermite', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$evalCovMat(...))
}

assign('Model_evalCovMat', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovAnisoListConst(self)$getCovMode(...))
}

assign('Model_getCovMode', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model_castInCovLMGradientConst(self)$evalZAndGradients(...))
}

assign('Model_evalZAndGradients', f , envir = asNamespace('gstlearn'))
%}
%insert(s)%{
f = function(self,...)
{
   return(ModelCovList_getCovList(self)$getNCov(...))
}

assign('ModelCovList_getNCov', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelCovList_getCovList(self)$getSills(...))
}

assign('ModelCovList_getSills', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelCovList_getCovList(self)$getSill(...))
}

assign('ModelCovList_getSill', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelCovList_getCovList(self)$getTotalSill(...))
}

assign('ModelCovList_getTotalSill', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelCovList_getCovList(self)$getTotalSills(...))
}

assign('ModelCovList_getTotalSills', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelCovList_getCovList(self)$isAllActiveCovList(...))
}

assign('ModelCovList_isAllActiveCovList', f , envir = asNamespace('gstlearn'))
%}
%insert(s)%{
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$setOptimEnabled(...))
}

assign('ModelGeneric_setOptimEnabled', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$attachNoStatDb(...))
}

assign('ModelGeneric_attachNoStatDb', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric_getCov(self)$makeStationary(...))
}

assign('ModelGeneric_makeStationary', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric__getCovModify(self)$setContext(...))
}

assign('ModelGeneric_setContext', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric__getDriftListModify(self)$setFlagLinked(...))
}

assign('ModelGeneric_setFlagLinked', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric__getDriftListModify(self)$setBetaHat(...))
}

assign('ModelGeneric_setBetaHat', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric__getDriftListModify(self)$setFiltered(...))
}

assign('ModelGeneric_setFiltered', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric__getDriftListModify(self)$delDrift(...))
}

assign('ModelGeneric_delDrift', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric__getDriftListModify(self)$delAllDrifts(...))
}

assign('ModelGeneric_delAllDrifts', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric__getDriftListModify(self)$copyCovContext(...))
}

assign('ModelGeneric_copyCovContext', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric__getDriftListModify(self)$setMeans(...))
}

assign('ModelGeneric_setMeans', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric__getDriftListModify(self)$setMean(...))
}

assign('ModelGeneric_setMean', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric__getContextModify(self)$setField(...))
}

assign('ModelGeneric_setField', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric__getContextModify(self)$setCovar0s(...))
}

assign('ModelGeneric_setCovar0s', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelGeneric__getContextModify(self)$setCovar0(...))
}

assign('ModelGeneric_setCovar0', f , envir = asNamespace('gstlearn'))
%}
%insert(s)%{
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$setActiveFactor(...))
}

assign('Model_setActiveFactor', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$getCovAniso(...))
}

assign('Model_getCovAniso', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$setSill(...))
}

assign('Model_setSill', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$setSills(...))
}

assign('Model_setSills', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$setRangeIsotropic(...))
}

assign('Model_setRangeIsotropic', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$setMarkovCoeffs(...))
}

assign('Model_setMarkovCoeffs', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$normalize(...))
}

assign('Model_normalize', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$makeRangeNoStatDb(...))
}

assign('Model_makeRangeNoStatDb', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$makeScaleNoStatDb(...))
}

assign('Model_makeScaleNoStatDb', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$makeAngleNoStatDb(...))
}

assign('Model_makeAngleNoStatDb', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$makeTensorNoStatDb(...))
}

assign('Model_makeTensorNoStatDb', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$makeParamNoStatDb(...))
}

assign('Model_makeParamNoStatDb', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$makeRangeNoStatFunctional(...))
}

assign('Model_makeRangeNoStatFunctional', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$makeScaleNoStatFunctional(...))
}

assign('Model_makeScaleNoStatFunctional', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$makeAngleNoStatFunctional(...))
}

assign('Model_makeAngleNoStatFunctional', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$makeTensorNoStatFunctional(...))
}

assign('Model_makeTensorNoStatFunctional', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$makeParamNoStatFunctional(...))
}

assign('Model_makeParamNoStatFunctional', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$makeRangeStationary(...))
}

assign('Model_makeRangeStationary', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$makeScaleStationary(...))
}

assign('Model_makeScaleStationary', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$makeAngleStationary(...))
}

assign('Model_makeAngleStationary', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$makeTensorStationary(...))
}

assign('Model_makeTensorStationary', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovAnisoList(self)$makeParamStationary(...))
}

assign('Model_makeParamStationary', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(Model__castInCovLMCTapering(self)$setTapeRange(...))
}

assign('Model_setTapeRange', f , envir = asNamespace('gstlearn'))
%}
%insert(s)%{
f = function(self,...)
{
   return(ModelCovList_getCovListModify(self)$delCov(...))
}

assign('ModelCovList_delCov', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelCovList_getCovListModify(self)$delAllCov(...))
}

assign('ModelCovList_delAllCov', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelCovList_getCovListModify(self)$setCovFiltered(...))
}

assign('ModelCovList_setCovFiltered', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelCovList_getCovListModify(self)$makeSillNoStatDb(...))
}

assign('ModelCovList_makeSillNoStatDb', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelCovList_getCovListModify(self)$makeSillStationary(...))
}

assign('ModelCovList_makeSillStationary', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelCovList_getCovListModify(self)$makeSillsStationary(...))
}

assign('ModelCovList_makeSillsStationary', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelCovList_getCovListModify(self)$makeSillNoStatFunctional(...))
}

assign('ModelCovList_makeSillNoStatFunctional', f , envir = asNamespace('gstlearn'))
f = function(self,...)
{
   return(ModelCovList_getCovListModify(self)$setSill(...))
}

assign('ModelCovList_setSill', f , envir = asNamespace('gstlearn'))
%}
