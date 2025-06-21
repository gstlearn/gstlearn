%pythoncode %{
import gstlearn as gl

def evalCovMat(self, *args, **kwargs):
    return self.getCov().evalCovMat(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalCovMat', evalCovMat)

def evalCovMatInPlace(self, *args, **kwargs):
    return self.getCov().evalCovMatInPlace(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalCovMatInPlace', evalCovMatInPlace)

def evalCovMatInPlaceFromIdx(self, *args, **kwargs):
    return self.getCov().evalCovMatInPlaceFromIdx(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalCovMatInPlaceFromIdx', evalCovMatInPlaceFromIdx)

def evalCovMatSym(self, *args, **kwargs):
    return self.getCov().evalCovMatSym(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalCovMatSym', evalCovMatSym)

def evalCovMatSymInPlace(self, *args, **kwargs):
    return self.getCov().evalCovMatSymInPlace(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalCovMatSymInPlace', evalCovMatSymInPlace)

def evalCovMatSymInPlaceFromIdx(self, *args, **kwargs):
    return self.getCov().evalCovMatSymInPlaceFromIdx(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalCovMatSymInPlaceFromIdx', evalCovMatSymInPlaceFromIdx)

def eval0Mat(self, *args, **kwargs):
    return self.getCov().eval0Mat(*args, **kwargs)

setattr(gl.ModelGeneric, 'eval0Mat', eval0Mat)

def evalCovMat0(self, *args, **kwargs):
    return self.getCov().evalCovMat0(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalCovMat0', evalCovMat0)

def evalCovMat0InPlace(self, *args, **kwargs):
    return self.getCov().evalCovMat0InPlace(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalCovMat0InPlace', evalCovMat0InPlace)

def evalCovVecRHSInPlace(self, *args, **kwargs):
    return self.getCov().evalCovVecRHSInPlace(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalCovVecRHSInPlace', evalCovVecRHSInPlace)

def evalCovMatOptimInPlace(self, *args, **kwargs):
    return self.getCov().evalCovMatOptimInPlace(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalCovMatOptimInPlace', evalCovMatOptimInPlace)

def evalCovMatRHSInPlaceFromIdx(self, *args, **kwargs):
    return self.getCov().evalCovMatRHSInPlaceFromIdx(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalCovMatRHSInPlaceFromIdx', evalCovMatRHSInPlaceFromIdx)

def evalCovMatSparse(self, *args, **kwargs):
    return self.getCov().evalCovMatSparse(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalCovMatSparse', evalCovMatSparse)

def eval0(self, *args, **kwargs):
    return self.getCov().eval0(*args, **kwargs)

setattr(gl.ModelGeneric, 'eval0', eval0)

def evalCov(self, *args, **kwargs):
    return self.getCov().evalCov(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalCov', evalCov)

def evalNvarIpas(self, *args, **kwargs):
    return self.getCov().evalNvarIpas(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalNvarIpas', evalNvarIpas)

def evalNvarIpasIncr(self, *args, **kwargs):
    return self.getCov().evalNvarIpasIncr(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalNvarIpasIncr', evalNvarIpasIncr)

def evalIvarNlag(self, *args, **kwargs):
    return self.getCov().evalIvarNlag(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalIvarNlag', evalIvarNlag)

def evalIvarIpas(self, *args, **kwargs):
    return self.getCov().evalIvarIpas(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalIvarIpas', evalIvarIpas)

def evalCvv(self, *args, **kwargs):
    return self.getCov().evalCvv(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalCvv', evalCvv)

def evalCvvShift(self, *args, **kwargs):
    return self.getCov().evalCvvShift(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalCvvShift', evalCvvShift)

def evalCvvM(self, *args, **kwargs):
    return self.getCov().evalCvvM(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalCvvM', evalCvvM)

def evalCxv(self, *args, **kwargs):
    return self.getCov().evalCxv(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalCxv', evalCxv)

def evalCxvM(self, *args, **kwargs):
    return self.getCov().evalCxvM(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalCxvM', evalCxvM)

def evalPointToDb(self, *args, **kwargs):
    return self.getCov().evalPointToDb(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalPointToDb', evalPointToDb)

def evalPointToDbAsSP(self, *args, **kwargs):
    return self.getCov().evalPointToDbAsSP(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalPointToDbAsSP', evalPointToDbAsSP)

def evalAverageDbToDb(self, *args, **kwargs):
    return self.getCov().evalAverageDbToDb(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalAverageDbToDb', evalAverageDbToDb)

def evalAverageIncrToIncr(self, *args, **kwargs):
    return self.getCov().evalAverageIncrToIncr(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalAverageIncrToIncr', evalAverageIncrToIncr)

def evalAveragePointToDb(self, *args, **kwargs):
    return self.getCov().evalAveragePointToDb(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalAveragePointToDb', evalAveragePointToDb)

def samplingDensityVariance(self, *args, **kwargs):
    return self.getCov().samplingDensityVariance(*args, **kwargs)

setattr(gl.ModelGeneric, 'samplingDensityVariance', samplingDensityVariance)

def specificVolume(self, *args, **kwargs):
    return self.getCov().specificVolume(*args, **kwargs)

setattr(gl.ModelGeneric, 'specificVolume', specificVolume)

def coefficientOfVariation(self, *args, **kwargs):
    return self.getCov().coefficientOfVariation(*args, **kwargs)

setattr(gl.ModelGeneric, 'coefficientOfVariation', coefficientOfVariation)

def specificVolumeFromCoV(self, *args, **kwargs):
    return self.getCov().specificVolumeFromCoV(*args, **kwargs)

setattr(gl.ModelGeneric, 'specificVolumeFromCoV', specificVolumeFromCoV)

def extensionVariance(self, *args, **kwargs):
    return self.getCov().extensionVariance(*args, **kwargs)

setattr(gl.ModelGeneric, 'extensionVariance', extensionVariance)

def calculateStDev(self, *args, **kwargs):
    return self.getCov().calculateStDev(*args, **kwargs)

setattr(gl.ModelGeneric, 'calculateStDev', calculateStDev)

def evaluateMatInPlace(self, *args, **kwargs):
    return self.getCov().evaluateMatInPlace(*args, **kwargs)

setattr(gl.ModelGeneric, 'evaluateMatInPlace', evaluateMatInPlace)

def evaluateOneGeneric(self, *args, **kwargs):
    return self.getCov().evaluateOneGeneric(*args, **kwargs)

setattr(gl.ModelGeneric, 'evaluateOneGeneric', evaluateOneGeneric)

def evaluateOneIncr(self, *args, **kwargs):
    return self.getCov().evaluateOneIncr(*args, **kwargs)

setattr(gl.ModelGeneric, 'evaluateOneIncr', evaluateOneIncr)

def buildVmapOnDbGrid(self, *args, **kwargs):
    return self.getCov().buildVmapOnDbGrid(*args, **kwargs)

setattr(gl.ModelGeneric, 'buildVmapOnDbGrid', buildVmapOnDbGrid)

def sample(self, *args, **kwargs):
    return self.getCov().sample(*args, **kwargs)

setattr(gl.ModelGeneric, 'sample', sample)

def sampleUnitary(self, *args, **kwargs):
    return self.getCov().sampleUnitary(*args, **kwargs)

setattr(gl.ModelGeneric, 'sampleUnitary', sampleUnitary)

def envelop(self, *args, **kwargs):
    return self.getCov().envelop(*args, **kwargs)

setattr(gl.ModelGeneric, 'envelop', envelop)

def gofToVario(self, *args, **kwargs):
    return self.getCov().gofToVario(*args, **kwargs)

setattr(gl.ModelGeneric, 'gofToVario', gofToVario)

def isNoStat(self, *args, **kwargs):
    return self.getCov().isNoStat(*args, **kwargs)

setattr(gl.ModelGeneric, 'isNoStat', isNoStat)

def manage(self, *args, **kwargs):
    return self.getCov().manage(*args, **kwargs)

setattr(gl.ModelGeneric, 'manage', manage)

def optimizationPreProcessForData(self, *args, **kwargs):
    return self.getCov().optimizationPreProcessForData(*args, **kwargs)

setattr(gl.ModelGeneric, 'optimizationPreProcessForData', optimizationPreProcessForData)

def optimizationPostProcess(self, *args, **kwargs):
    return self.getCov().optimizationPostProcess(*args, **kwargs)

setattr(gl.ModelGeneric, 'optimizationPostProcess', optimizationPostProcess)

def getDrift(self, *args, **kwargs):
    return self.getDriftList().getDrift(*args, **kwargs)

setattr(gl.ModelGeneric, 'getDrift', getDrift)

def computeDrift(self, *args, **kwargs):
    return self.getDriftList().computeDrift(*args, **kwargs)

setattr(gl.ModelGeneric, 'computeDrift', computeDrift)

def evalDriftValue(self, *args, **kwargs):
    return self.getDriftList().evalDriftValue(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalDriftValue', evalDriftValue)

def evalDriftMat(self, *args, **kwargs):
    return self.getDriftList().evalDriftMat(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalDriftMat', evalDriftMat)

def evalDriftMatInPlace(self, *args, **kwargs):
    return self.getDriftList().evalDriftMatInPlace(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalDriftMatInPlace', evalDriftMatInPlace)

def evalDriftMatByRanks(self, *args, **kwargs):
    return self.getDriftList().evalDriftMatByRanks(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalDriftMatByRanks', evalDriftMatByRanks)

def evalMeanVecByRanks(self, *args, **kwargs):
    return self.getDriftList().evalMeanVecByRanks(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalMeanVecByRanks', evalMeanVecByRanks)

def evalDriftMatByRanksInPlace(self, *args, **kwargs):
    return self.getDriftList().evalDriftMatByRanksInPlace(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalDriftMatByRanksInPlace', evalDriftMatByRanksInPlace)

def evalDriftMatByTargetInPlace(self, *args, **kwargs):
    return self.getDriftList().evalDriftMatByTargetInPlace(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalDriftMatByTargetInPlace', evalDriftMatByTargetInPlace)

def getNDrift(self, *args, **kwargs):
    return self.getDriftList().getNDrift(*args, **kwargs)

setattr(gl.ModelGeneric, 'getNDrift', getNDrift)

def getNDriftEquation(self, *args, **kwargs):
    return self.getDriftList().getNDriftEquation(*args, **kwargs)

setattr(gl.ModelGeneric, 'getNDriftEquation', getNDriftEquation)

def getNExtDrift(self, *args, **kwargs):
    return self.getDriftList().getNExtDrift(*args, **kwargs)

setattr(gl.ModelGeneric, 'getNExtDrift', getNExtDrift)

def isFlagLinked(self, *args, **kwargs):
    return self.getDriftList().isFlagLinked(*args, **kwargs)

setattr(gl.ModelGeneric, 'isFlagLinked', isFlagLinked)

def getDriftMaxIRFOrder(self, *args, **kwargs):
    return self.getDriftList().getDriftMaxIRFOrder(*args, **kwargs)

setattr(gl.ModelGeneric, 'getDriftMaxIRFOrder', getDriftMaxIRFOrder)

def getRankFex(self, *args, **kwargs):
    return self.getDriftList().getRankFex(*args, **kwargs)

setattr(gl.ModelGeneric, 'getRankFex', getRankFex)

def isDriftSampleDefined(self, *args, **kwargs):
    return self.getDriftList().isDriftSampleDefined(*args, **kwargs)

setattr(gl.ModelGeneric, 'isDriftSampleDefined', isDriftSampleDefined)

def isDriftFiltered(self, *args, **kwargs):
    return self.getDriftList().isDriftFiltered(*args, **kwargs)

setattr(gl.ModelGeneric, 'isDriftFiltered', isDriftFiltered)

def isDriftDefined(self, *args, **kwargs):
    return self.getDriftList().isDriftDefined(*args, **kwargs)

setattr(gl.ModelGeneric, 'isDriftDefined', isDriftDefined)

def isDriftDifferentDefined(self, *args, **kwargs):
    return self.getDriftList().isDriftDifferentDefined(*args, **kwargs)

setattr(gl.ModelGeneric, 'isDriftDifferentDefined', isDriftDifferentDefined)

def getDrifts(self, *args, **kwargs):
    return self.getDriftList().getDrifts(*args, **kwargs)

setattr(gl.ModelGeneric, 'getDrifts', getDrifts)

def evalDrift(self, *args, **kwargs):
    return self.getDriftList().evalDrift(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalDrift', evalDrift)

def evalDriftBySample(self, *args, **kwargs):
    return self.getDriftList().evalDriftBySample(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalDriftBySample', evalDriftBySample)

def evalDriftBySampleInPlace(self, *args, **kwargs):
    return self.getDriftList().evalDriftBySampleInPlace(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalDriftBySampleInPlace', evalDriftBySampleInPlace)

def evalDriftCoef(self, *args, **kwargs):
    return self.getDriftList().evalDriftCoef(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalDriftCoef', evalDriftCoef)

def hasDrift(self, *args, **kwargs):
    return self.getDriftList().hasDrift(*args, **kwargs)

setattr(gl.ModelGeneric, 'hasDrift', hasDrift)

def getMean(self, *args, **kwargs):
    return self.getDriftList().getMean(*args, **kwargs)

setattr(gl.ModelGeneric, 'getMean', getMean)

def getMeans(self, *args, **kwargs):
    return self.getDriftList().getMeans(*args, **kwargs)

setattr(gl.ModelGeneric, 'getMeans', getMeans)

def evalDriftVarCoef(self, *args, **kwargs):
    return self.getDriftList().evalDriftVarCoef(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalDriftVarCoef', evalDriftVarCoef)

def evalDriftVarCoefs(self, *args, **kwargs):
    return self.getDriftList().evalDriftVarCoefs(*args, **kwargs)

setattr(gl.ModelGeneric, 'evalDriftVarCoefs', evalDriftVarCoefs)

def getNVar(self, *args, **kwargs):
    return self.getContext().getNVar(*args, **kwargs)

setattr(gl.ModelGeneric, 'getNVar', getNVar)

def getNDim(self, *args, **kwargs):
    return self.getContext().getNDim(*args, **kwargs)

setattr(gl.ModelGeneric, 'getNDim', getNDim)

def getSpace(self, *args, **kwargs):
    return self.getContext().getSpace(*args, **kwargs)

setattr(gl.ModelGeneric, 'getSpace', getSpace)

def getCovar0(self, *args, **kwargs):
    return self.getContext().getCovar0(*args, **kwargs)

setattr(gl.ModelGeneric, 'getCovar0', getCovar0)

def getField(self, *args, **kwargs):
    return self.getContext().getField(*args, **kwargs)

setattr(gl.ModelGeneric, 'getField', getField)
%}
%pythoncode %{
import gstlearn as gl

def getActiveFactor(self, *args, **kwargs):
    return self.castInCovAnisoListConst().getActiveFactor(*args, **kwargs)

setattr(gl.Model, 'getActiveFactor', getActiveFactor)

def getCovAniso(self, *args, **kwargs):
    return self.castInCovAnisoListConst().getCovAniso(*args, **kwargs)

setattr(gl.Model, 'getCovAniso', getCovAniso)

def getNCov(self, *args, **kwargs):
    return self.castInCovAnisoListConst().getNCov(*args, **kwargs)

setattr(gl.Model, 'getNCov', getNCov)

def getCovType(self, *args, **kwargs):
    return self.castInCovAnisoListConst().getCovType(*args, **kwargs)

setattr(gl.Model, 'getCovType', getCovType)

def getRange(self, *args, **kwargs):
    return self.castInCovAnisoListConst().getRange(*args, **kwargs)

setattr(gl.Model, 'getRange', getRange)

def getRanges(self, *args, **kwargs):
    return self.castInCovAnisoListConst().getRanges(*args, **kwargs)

setattr(gl.Model, 'getRanges', getRanges)

def getAngles(self, *args, **kwargs):
    return self.castInCovAnisoListConst().getAngles(*args, **kwargs)

setattr(gl.Model, 'getAngles', getAngles)

def getAnam(self, *args, **kwargs):
    return self.castInCovAnisoListConst().getAnam(*args, **kwargs)

setattr(gl.Model, 'getAnam', getAnam)

def getParam(self, *args, **kwargs):
    return self.castInCovAnisoListConst().getParam(*args, **kwargs)

setattr(gl.Model, 'getParam', getParam)

def getCovName(self, *args, **kwargs):
    return self.castInCovAnisoListConst().getCovName(*args, **kwargs)

setattr(gl.Model, 'getCovName', getCovName)

def extractCova(self, *args, **kwargs):
    return self.castInCovAnisoListConst().extractCova(*args, **kwargs)

setattr(gl.Model, 'extractCova', extractCova)

def getNGradParam(self, *args, **kwargs):
    return self.castInCovAnisoListConst().getNGradParam(*args, **kwargs)

setattr(gl.Model, 'getNGradParam', getNGradParam)

def getMaximumDistance(self, *args, **kwargs):
    return self.castInCovAnisoListConst().getMaximumDistance(*args, **kwargs)

setattr(gl.Model, 'getMaximumDistance', getMaximumDistance)

def getCovMinIRFOrder(self, *args, **kwargs):
    return self.castInCovAnisoListConst().getCovMinIRFOrder(*args, **kwargs)

setattr(gl.Model, 'getCovMinIRFOrder', getCovMinIRFOrder)

def getAnamNClass(self, *args, **kwargs):
    return self.castInCovAnisoListConst().getAnamNClass(*args, **kwargs)

setattr(gl.Model, 'getAnamNClass', getAnamNClass)

def hasAnam(self, *args, **kwargs):
    return self.castInCovAnisoListConst().hasAnam(*args, **kwargs)

setattr(gl.Model, 'hasAnam', hasAnam)

def hasNugget(self, *args, **kwargs):
    return self.castInCovAnisoListConst().hasNugget(*args, **kwargs)

setattr(gl.Model, 'hasNugget', hasNugget)

def getRankNugget(self, *args, **kwargs):
    return self.castInCovAnisoListConst().getRankNugget(*args, **kwargs)

setattr(gl.Model, 'getRankNugget', getRankNugget)

def getBallRadius(self, *args, **kwargs):
    return self.castInCovAnisoListConst().getBallRadius(*args, **kwargs)

setattr(gl.Model, 'getBallRadius', getBallRadius)

def hasExternalCov(self, *args, **kwargs):
    return self.castInCovAnisoListConst().hasExternalCov(*args, **kwargs)

setattr(gl.Model, 'hasExternalCov', hasExternalCov)

def isChangeSupportDefined(self, *args, **kwargs):
    return self.castInCovAnisoListConst().isChangeSupportDefined(*args, **kwargs)

setattr(gl.Model, 'isChangeSupportDefined', isChangeSupportDefined)

def getAnamHermite(self, *args, **kwargs):
    return self.castInCovAnisoListConst().getAnamHermite(*args, **kwargs)

setattr(gl.Model, 'getAnamHermite', getAnamHermite)

def evalCovMat(self, *args, **kwargs):
    return self.castInCovAnisoListConst().evalCovMat(*args, **kwargs)

setattr(gl.Model, 'evalCovMat', evalCovMat)

def getCovMode(self, *args, **kwargs):
    return self.castInCovAnisoListConst().getCovMode(*args, **kwargs)

setattr(gl.Model, 'getCovMode', getCovMode)

def evalZAndGradients(self, *args, **kwargs):
    return self.castInCovLMGradientConst().evalZAndGradients(*args, **kwargs)

setattr(gl.Model, 'evalZAndGradients', evalZAndGradients)
%}
%pythoncode %{
import gstlearn as gl

def getNCov(self, *args, **kwargs):
    return self.getCovList().getNCov(*args, **kwargs)

setattr(gl.ModelCovList, 'getNCov', getNCov)

def getSills(self, *args, **kwargs):
    return self.getCovList().getSills(*args, **kwargs)

setattr(gl.ModelCovList, 'getSills', getSills)

def getSill(self, *args, **kwargs):
    return self.getCovList().getSill(*args, **kwargs)

setattr(gl.ModelCovList, 'getSill', getSill)

def getTotalSill(self, *args, **kwargs):
    return self.getCovList().getTotalSill(*args, **kwargs)

setattr(gl.ModelCovList, 'getTotalSill', getTotalSill)

def getTotalSills(self, *args, **kwargs):
    return self.getCovList().getTotalSills(*args, **kwargs)

setattr(gl.ModelCovList, 'getTotalSills', getTotalSills)

def isAllActiveCovList(self, *args, **kwargs):
    return self.getCovList().isAllActiveCovList(*args, **kwargs)

setattr(gl.ModelCovList, 'isAllActiveCovList', isAllActiveCovList)
%}
%pythoncode %{
import gstlearn as gl

def setOptimEnabled(self, *args, **kwargs):
    return self.getCov().setOptimEnabled(*args, **kwargs)

setattr(gl.ModelGeneric, 'setOptimEnabled', setOptimEnabled)

def attachNoStatDb(self, *args, **kwargs):
    return self.getCov().attachNoStatDb(*args, **kwargs)

setattr(gl.ModelGeneric, 'attachNoStatDb', attachNoStatDb)

def makeStationary(self, *args, **kwargs):
    return self.getCov().makeStationary(*args, **kwargs)

setattr(gl.ModelGeneric, 'makeStationary', makeStationary)

def setContext(self, *args, **kwargs):
    return self._getCovModify().setContext(*args, **kwargs)

setattr(gl.ModelGeneric, 'setContext', setContext)

def setFlagLinked(self, *args, **kwargs):
    return self._getDriftListModify().setFlagLinked(*args, **kwargs)

setattr(gl.ModelGeneric, 'setFlagLinked', setFlagLinked)

def setBetaHat(self, *args, **kwargs):
    return self._getDriftListModify().setBetaHat(*args, **kwargs)

setattr(gl.ModelGeneric, 'setBetaHat', setBetaHat)

def setFiltered(self, *args, **kwargs):
    return self._getDriftListModify().setFiltered(*args, **kwargs)

setattr(gl.ModelGeneric, 'setFiltered', setFiltered)

def delDrift(self, *args, **kwargs):
    return self._getDriftListModify().delDrift(*args, **kwargs)

setattr(gl.ModelGeneric, 'delDrift', delDrift)

def delAllDrifts(self, *args, **kwargs):
    return self._getDriftListModify().delAllDrifts(*args, **kwargs)

setattr(gl.ModelGeneric, 'delAllDrifts', delAllDrifts)

def copyCovContext(self, *args, **kwargs):
    return self._getDriftListModify().copyCovContext(*args, **kwargs)

setattr(gl.ModelGeneric, 'copyCovContext', copyCovContext)

def setMeans(self, *args, **kwargs):
    return self._getDriftListModify().setMeans(*args, **kwargs)

setattr(gl.ModelGeneric, 'setMeans', setMeans)

def setMean(self, *args, **kwargs):
    return self._getDriftListModify().setMean(*args, **kwargs)

setattr(gl.ModelGeneric, 'setMean', setMean)

def setField(self, *args, **kwargs):
    return self._getContextModify().setField(*args, **kwargs)

setattr(gl.ModelGeneric, 'setField', setField)

def setCovar0s(self, *args, **kwargs):
    return self._getContextModify().setCovar0s(*args, **kwargs)

setattr(gl.ModelGeneric, 'setCovar0s', setCovar0s)

def setCovar0(self, *args, **kwargs):
    return self._getContextModify().setCovar0(*args, **kwargs)

setattr(gl.ModelGeneric, 'setCovar0', setCovar0)
%}
%pythoncode %{
import gstlearn as gl

def setActiveFactor(self, *args, **kwargs):
    return self._castInCovAnisoList().setActiveFactor(*args, **kwargs)

setattr(gl.Model, 'setActiveFactor', setActiveFactor)

def getCovAniso(self, *args, **kwargs):
    return self._castInCovAnisoList().getCovAniso(*args, **kwargs)

setattr(gl.Model, 'getCovAniso', getCovAniso)

def setSill(self, *args, **kwargs):
    return self._castInCovAnisoList().setSill(*args, **kwargs)

setattr(gl.Model, 'setSill', setSill)

def setSills(self, *args, **kwargs):
    return self._castInCovAnisoList().setSills(*args, **kwargs)

setattr(gl.Model, 'setSills', setSills)

def setRangeIsotropic(self, *args, **kwargs):
    return self._castInCovAnisoList().setRangeIsotropic(*args, **kwargs)

setattr(gl.Model, 'setRangeIsotropic', setRangeIsotropic)

def setMarkovCoeffs(self, *args, **kwargs):
    return self._castInCovAnisoList().setMarkovCoeffs(*args, **kwargs)

setattr(gl.Model, 'setMarkovCoeffs', setMarkovCoeffs)

def normalize(self, *args, **kwargs):
    return self._castInCovAnisoList().normalize(*args, **kwargs)

setattr(gl.Model, 'normalize', normalize)

def makeRangeNoStatDb(self, *args, **kwargs):
    return self._castInCovAnisoList().makeRangeNoStatDb(*args, **kwargs)

setattr(gl.Model, 'makeRangeNoStatDb', makeRangeNoStatDb)

def makeScaleNoStatDb(self, *args, **kwargs):
    return self._castInCovAnisoList().makeScaleNoStatDb(*args, **kwargs)

setattr(gl.Model, 'makeScaleNoStatDb', makeScaleNoStatDb)

def makeAngleNoStatDb(self, *args, **kwargs):
    return self._castInCovAnisoList().makeAngleNoStatDb(*args, **kwargs)

setattr(gl.Model, 'makeAngleNoStatDb', makeAngleNoStatDb)

def makeTensorNoStatDb(self, *args, **kwargs):
    return self._castInCovAnisoList().makeTensorNoStatDb(*args, **kwargs)

setattr(gl.Model, 'makeTensorNoStatDb', makeTensorNoStatDb)

def makeParamNoStatDb(self, *args, **kwargs):
    return self._castInCovAnisoList().makeParamNoStatDb(*args, **kwargs)

setattr(gl.Model, 'makeParamNoStatDb', makeParamNoStatDb)

def makeRangeNoStatFunctional(self, *args, **kwargs):
    return self._castInCovAnisoList().makeRangeNoStatFunctional(*args, **kwargs)

setattr(gl.Model, 'makeRangeNoStatFunctional', makeRangeNoStatFunctional)

def makeScaleNoStatFunctional(self, *args, **kwargs):
    return self._castInCovAnisoList().makeScaleNoStatFunctional(*args, **kwargs)

setattr(gl.Model, 'makeScaleNoStatFunctional', makeScaleNoStatFunctional)

def makeAngleNoStatFunctional(self, *args, **kwargs):
    return self._castInCovAnisoList().makeAngleNoStatFunctional(*args, **kwargs)

setattr(gl.Model, 'makeAngleNoStatFunctional', makeAngleNoStatFunctional)

def makeTensorNoStatFunctional(self, *args, **kwargs):
    return self._castInCovAnisoList().makeTensorNoStatFunctional(*args, **kwargs)

setattr(gl.Model, 'makeTensorNoStatFunctional', makeTensorNoStatFunctional)

def makeParamNoStatFunctional(self, *args, **kwargs):
    return self._castInCovAnisoList().makeParamNoStatFunctional(*args, **kwargs)

setattr(gl.Model, 'makeParamNoStatFunctional', makeParamNoStatFunctional)

def makeRangeStationary(self, *args, **kwargs):
    return self._castInCovAnisoList().makeRangeStationary(*args, **kwargs)

setattr(gl.Model, 'makeRangeStationary', makeRangeStationary)

def makeScaleStationary(self, *args, **kwargs):
    return self._castInCovAnisoList().makeScaleStationary(*args, **kwargs)

setattr(gl.Model, 'makeScaleStationary', makeScaleStationary)

def makeAngleStationary(self, *args, **kwargs):
    return self._castInCovAnisoList().makeAngleStationary(*args, **kwargs)

setattr(gl.Model, 'makeAngleStationary', makeAngleStationary)

def makeTensorStationary(self, *args, **kwargs):
    return self._castInCovAnisoList().makeTensorStationary(*args, **kwargs)

setattr(gl.Model, 'makeTensorStationary', makeTensorStationary)

def makeParamStationary(self, *args, **kwargs):
    return self._castInCovAnisoList().makeParamStationary(*args, **kwargs)

setattr(gl.Model, 'makeParamStationary', makeParamStationary)

def setTapeRange(self, *args, **kwargs):
    return self._castInCovLMCTapering().setTapeRange(*args, **kwargs)

setattr(gl.Model, 'setTapeRange', setTapeRange)
%}
%pythoncode %{
import gstlearn as gl

def delCov(self, *args, **kwargs):
    return self.getCovListModify().delCov(*args, **kwargs)

setattr(gl.ModelCovList, 'delCov', delCov)

def delAllCov(self, *args, **kwargs):
    return self.getCovListModify().delAllCov(*args, **kwargs)

setattr(gl.ModelCovList, 'delAllCov', delAllCov)

def setCovFiltered(self, *args, **kwargs):
    return self.getCovListModify().setCovFiltered(*args, **kwargs)

setattr(gl.ModelCovList, 'setCovFiltered', setCovFiltered)

def makeSillNoStatDb(self, *args, **kwargs):
    return self.getCovListModify().makeSillNoStatDb(*args, **kwargs)

setattr(gl.ModelCovList, 'makeSillNoStatDb', makeSillNoStatDb)

def makeSillStationary(self, *args, **kwargs):
    return self.getCovListModify().makeSillStationary(*args, **kwargs)

setattr(gl.ModelCovList, 'makeSillStationary', makeSillStationary)

def makeSillsStationary(self, *args, **kwargs):
    return self.getCovListModify().makeSillsStationary(*args, **kwargs)

setattr(gl.ModelCovList, 'makeSillsStationary', makeSillsStationary)

def makeSillNoStatFunctional(self, *args, **kwargs):
    return self.getCovListModify().makeSillNoStatFunctional(*args, **kwargs)

setattr(gl.ModelCovList, 'makeSillNoStatFunctional', makeSillNoStatFunctional)

def setSill(self, *args, **kwargs):
    return self.getCovListModify().setSill(*args, **kwargs)

setattr(gl.ModelCovList, 'setSill', setSill)
%}
