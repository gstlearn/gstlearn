import gstlearn as gl
import gstlearn.test as gt
import numpy as np
from scipy.spatial import distance_matrix

def cova(x,sills=1):
    return np.kron(sills,np.exp(-x/2))

A = np.random.normal(size=(3,3))
sills = gl.VectorDouble((A@A.T).reshape(1,-1)[0])
model = gl.Model.createFromParam(gl.ECov.EXPONENTIAL,range = 2.,flagRange=False,sills=sills)

nx = [10,10]

def matriceReduce(m,ind):
    M = gl.MatrixSquareSymmetric(len(ind))
    for i,j in enumerate(ind):
        for k in range(j,len(ind)):
            M.setValue(i,k,m.getValue(j,ind[k]))
    return M

def modelReduce(model,var):
    ctxt=model.getContext()
    ctxt.setNVar(len(var))
    modeln = gl.Model.create(ctxt)
    covlist = model.getCovAnisoList()
    
    for i in range(covlist.getCovNumber()):
        cova = covlist.getCova(i)
        sills = matriceReduce(cova.getSill(),var)
        covar = gl.CovAniso.createAnisotropicMulti(ctxt,
                                             cova.getType(),
                                             cova.getRanges(),
                                             sills,
                                             cova.getParam(),
                                             cova.getAnisoAngles())
        modeln.addCov(covar)
    return modeln


def createDbIn(ndat,nvar,percent,ndim=2,selDbin=False,measurement_error=False,ndrift = 0,
               flag_isotopic=False,seed=1234):
    db = gl.Db.create()
    np.random.seed(seed)
    for i in range(ndim):
        db["x" +str(i)] = np.random.uniform(size = ndat)
     
    db.setLocators(["x*"],gl.ELoc.X)
        
    indIn = np.arange(ndat)
    if selDbin:
        np.random.shuffle(indIn)
        indIn = np.sort(indIn)
        indIn = indIn[range(int(ndat/2))]
        sel = np.zeros(shape=ndat)
        sel[indIn] = 1
        db["sel"] = sel
        db.setLocator("sel",gl.ELoc.SEL)
      
    #Creation of an heterotopic data set
    indList = [] 
    for i in range(nvar):
        u = np.array([None for i in range(ndat)])
        ind = np.arange(ndat)
        if not flag_isotopic and nvar>1: 
            np.random.shuffle(ind)
            ind = ind[range(int(ndat*percent[i]))]
        ind = np.sort(list(set(ind) & set(indIn)))
        indList += [ind]
        vect = np.array([None for i in range(ndat)])
        vect[ind] = np.random.normal(size = len(ind))
        db["z"+str(i)]=vect
          
    db.setLocators(["z*"],gl.ELoc.Z)
    
    indF = []
    
    for i in range(nvar):
        indF += list(np.array(indList[i]) + ndat * i)
    
    if measurement_error :
        for i in range(nvar):
            db["err"+str(i)] = 0.1 * np.random.uniform(size = ndat)
            
        db.setLocators(["err*"],gl.ELoc.V)
    
    if ndrift>0:
        for i in range(ndrift):
            db["ff" + str(i)] = np.random.normal(size = ndat)
        db.setLocator("ff*",gl.ELoc.F)
    
    return db,indF


def test_kriging(ndat,nx,nvar,percent,model,cova,
                 irf=None,drift=False,measurement_error=False,compute_vars = True,
                 selDbin = True, selDbout = True,flag_isotopic=True,
                 seed=1234,tol=1e-12,eps=1e-3,test = True,verbose=False):
    
    np.random.seed(seed)
    ndrift = 1 if drift else 0
    modeln = modelReduce(model,range(nvar))
    #### Create the description of the case #####
    case = "case:\n"
    
    inter = ""
    if nvar > 1:
        inter = "co-"
    if irf is None and not drift:
        case += "- simple "+ inter+ "kriging\n"
    else :
        if irf is not None :
            case += "- KU with drift of degree " + str(irf) + "\n"
        if drift :
            case +="- with external drift\n"
    if nvar>1:
        case +="- number of covariables for co-kriging " + str(nvar) + "\n"
        if flag_isotopic:
            case += "- isotopic case\n"
        else:
            case += "- heterotopic case\n"
    if measurement_error:
        case += "- with measurement error\n"
    else:
        case += "- without measurement error\n"
    if compute_vars:
        case += "- no dual form\n"
    else:
        case += "- dual\n"
    case += "- selection on Dbin " + str(selDbin) + "\n"
    case += "- selection on Dbout "+ str(selDbout) + "\n"
    case += "- number of data " + str(ndat) + "\n"
    case += "- nx = ["+str(nx[0]) +"," + str(nx[1]) + "]\n"
    
    if verbose:
        print(case)
    ##################################################
    db,indF = createDbIn(ndat,nvar,percent,2,selDbin,measurement_error,ndrift,flag_isotopic,seed)
    
    target = gl.DbGrid.create(nx = nx)
   
    indOut = np.arange(target.getSampleNumber())
    
    if selDbout:
        np.random.shuffle(indOut)
        indOut = indOut[range(int(target.getSampleNumber()/2))]
        indOut = np.sort(indOut)
        sel = np.zeros(shape = target.getSampleNumber())
        sel[indOut] = 1
        target["sel"] = sel
        target.setLocator("sel",gl.ELoc.SEL)
                  
    if irf is not None:
        modeln.setDriftIRF(irf)
    
    if drift :
        target["ff"] = np.random.normal(size = target.getSampleNumber())
        
        target.setLocator("ff",gl.ELoc.F)
        modeln.addDrift(gl.DriftF(0, modeln.getContext()))
      
    v = np.array([db["x0"],db["x1"]]).T
    v0 = np.array([target["x1"][indOut],target["x2"][indOut]]).T
    cov =cova(distance_matrix(v,v),modeln.getSill(0).toTL())[indF,:][:,indF]
    c0 = cova(distance_matrix(v,v0),modeln.getSill(0).toTL())[indF,:]
    
    #Creation of a db2 without selection to build the complete covariance matrix
    db2 = db.clone()
    db2.setLocator("sel")
    vect = gl.VectorDouble(nvar**2 * db2.getSampleNumber()**2)
    
    target2 = target.clone()
    target2.setLocator("sel")
    #modeln.covMatrix(vect,db2,None,-1,-1)
    ##########
    
    covgl = np.array(list(vect)).reshape(nvar * db2.getSampleNumber(),-1)[indF,:][:,indF]
    
    if measurement_error:
        err = db["err*"].T.reshape(-1,)
        np.fill_diagonal(cov,np.diag(cov)+err[indF])
    
    
    #assert(np.linalg.norm(covgl-cov)<tol)   
    
    vect = gl.VectorDouble(nvar**2 * db2.getSampleNumber() * len(indOut))
    #modeln.covMatrix(vect,db2,target,-1,-1)
    #c0gl = np.array(list(vect)).reshape(nvar*db2.getSampleNumber(),nvar*len(indOut))[indF,:]
    
    #assert(np.linalg.norm(c0gl-c0)<tol)
    
    neigh = gl.NeighUnique()
    
    if compute_vars:
        gl.kriging(db,target,modeln,neigh,flag_std = True,flag_varz=True)
        #gl.krigingExperimentalEigen(db,target,modeln,flag_std = True,flag_varz=True)
    else :
        gl.kriging(db,target,modeln,neigh,flag_std = False,flag_varz=False)
        #gl.krigingExperimentalEigen(db,target,modeln,flag_std=False,flag_varz=False)
    
    if irf is not None or drift:
        c0t = np.copy(c0)
        driftd = np.kron(np.eye(nvar),modeln.getDrifts(db2, False))[:,indF]
        driftt = np.kron(np.eye(nvar),modeln.getDrifts(target, True))
        temp = np.linalg.solve(cov,driftd.T)
        term = (driftt-temp.T@c0).T@np.linalg.inv(driftd@temp)@driftd
        c0t= term.T+c0
    else:
        c0t = c0
    
    weights = np.linalg.solve(cov,c0t)
    
    ntarget = target.getSampleNumber(useSel=True)
    krig = (db["z*"].T.reshape(1,-1).T[indF].T@weights).reshape(-1,)
    
    varest = np.sum(c0t*weights,axis=0).T.reshape(-1,)
    
    nt = target.getSampleNumber(useSel=True)
    c0v = np.atleast_2d(np.ones(nt)).T@np.atleast_2d(np.array([modeln.getTotalSill(i,i) for i in range(nvar)]))
    c0v = c0v.T.reshape(-1,)
    var = c0v + varest - 2 * np.sum(weights*c0,axis=0).T.reshape(-1,)
    std = np.sqrt(var)
    
    krigref = target["*estim"][indOut].T.reshape(-1,)
    stdref = target["*stdev"][indOut].T.reshape(-1,)
    varestref = target["*varz"][indOut].T.reshape(-1,)
    if test :
        gt.checkEqualityVector(krigref, krig, tolerance=tol)
        if compute_vars:
            #assert np.linalg.norm((stdref-std)/(eps+std))<tol , "Problem with kriging stdev in " + case
            #assert np.linalg.norm((varestref-varest)/(eps+varest))<tol, "Problem with kriging variance (var_est) in " + case
            gt.checkEqualityVector(varestref, varest, tolerance=tol)
    if verbose:
        print("Test Ok")
    return krig,target,indOut,db , varest


percent = [0.5,0.9,1.]
ndat = 40
nbtests = 0
verbose = False
general = True

if general:
    for irf in [None,0,1]:
        for drift in [False,True]:
            for measurement_error in [True, False]:
                for selDbin in [True, False]:
                    for selDbout in [True, False]:
                        for nx in [[5,5]]:
                            for nvar in [1,2,3]:
                                isolist = [True, False]
                                if nvar >1 :
                                    isolist = [True,False]
                                for iso in isolist:
                                    for cv in [False,True]:
                                        a = test_kriging(ndat,nx,nvar,percent,
                                                         model,cova,compute_vars=cv,
                                                         irf=irf,drift=drift,
                                                         measurement_error=measurement_error,
                                                         selDbin=selDbin,selDbout=selDbout,
                                                         flag_isotopic = iso,
                                                         seed=1234,tol=1e-8,eps=1e-3,verbose=verbose)
                                        nbtests += 1

    print("Number of tests performed =",nbtests)

# Individual test

if not general:
    irf = 0
    drift = False
    measurement_error = False
    selDbin = False
    selDbout = False
    nx = [5,5]
    nvar = 2
    iso = True
    cv = False
                                    
    a = test_kriging(ndat,nx,nvar,percent,
                     model,cova,compute_vars=cv,
                     irf=irf,drift=drift,
                     measurement_error=measurement_error,
                     selDbin=selDbin,selDbout=selDbout,
                     flag_isotopic = iso,
                     seed=1234,tol=1e-8,eps=1e-3,verbose=verbose)
