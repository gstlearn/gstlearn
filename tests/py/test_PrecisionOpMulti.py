# %%
import gstlearn as gl
import numpy as np
import scipy.sparse.linalg
def logit(x,a=1,b=0):
    return 1/(1+np.exp(-a*x+b))


# %%
def create(nvar = 1, multistruct = True, nostatType = "Fake",nx1 = [4,4],nx2 = [4,4]):

    if nostatType == "Fake":
        nostat= True
        p = 0
    if nostatType == "Real":
        nostat = True
        p = 1
    if nostatType == "Stat":
        nostat = False

    meshes = [gl.MeshETurbo(nx1)]
    if multistruct:
        meshes +=[gl.MeshETurbo(nx2)]
    
    eps =0.1

    if nvar == 1:
        sills1 = np.array(2.)
        sills2 = np.array(5.)
    if nvar == 2:
        r1 = 0.
        s11 = 2
        s21 = 1
        s121 = r1 * np.sqrt(s11*s21)
        sills1 = np.array([s11,s121,s121,s21])
      
        r2 = 0.3
        s12 = .01
        s22 = .3
        s122 = r2 * np.sqrt(s12*s22)
        sills2 = np.array([s11,s121,s121,s21])
        
    if nvar == 3:
        np.random.seed(14556)
        sills1 = np.random.normal(size=[3,3])
        sills1 = sills1@sills1.T

        sills2 = np.random.normal(size=[3,3])
        sills2 = sills2@sills2.T
        
    modelMulti = gl.Model.createFromParam(gl.ECov.MATERN,param=1,sills = sills1.reshape(-1),range = 20)

    if multistruct:
        modelMulti2 = gl.Model.createFromParam(gl.ECov.MATERN,param=2,sills = sills2.reshape(-1),range = 10)
        modelMulti.addCov(modelMulti2.getCova(0))

    if nostat:
        grid  = gl.DbGrid.create(nx1)
        
        if nvar == 1:
            if multistruct :
                grid2 = grid.clone()
            rho1 = eps + (grid["x1"]-np.min(grid["x1"])+eps) / (eps + np.max(grid["x1"])-np.min(grid["x1"]))
            rho2 = eps + (grid["x2"]-np.min(grid["x2"])+eps) / (eps + np.max(grid["x2"])-np.min(grid["x2"]))
            grid["rho1"] = p * rho1 + (1-p) * sills1
            modelMulti.getCova(0).attachNoStatDb(grid)
            modelMulti.getCova(0).makeSillNoStatDb("rho1")
            if multistruct :
                grid2["rho2"] = p * rho2 + (1-p) * sills2
                modelMulti.getCova(1).attachNoStatDb(grid2)
                modelMulti.getCova(1).makeSillNoStatDb("rho2")

        if nvar == 2:
            rho = (eps + grid["x1"]-np.min(grid["x1"])) / (eps + np.max(grid["x1"])-np.min(grid["x1"]))
            rho = logit(rho,20,10)
            grid["rho"] = rho * np.sqrt(s11*s21) * p + s121 * (1-p)
            modelMulti.getCova(0).makeSillNoStatDb("rho",0,1,grid)
        if nvar == 3:
            grid["rho"] = sills1[0,1] * np.ones_like(grid["rank"])
            modelMulti.getCova(0).makeSillNoStatDb("rho",0,1,grid)
       
    return modelMulti,meshes


# %%
class PrecisionOpMulti:
    def __init__(self,model,meshes,matrix = False):
        modelC = model.clone()
        if matrix :
            createQ = gl.PrecisionOpMatrix
        else:
            createQ = gl.PrecisionOp
        self.sills = []
        self.invsill = []
        self.cholsill = []
        self.nvertex = []
        self.SigmaMult = []
        self.invSigmaMult = []
        self.Qop = []
        self.temp = []
        self.ncovar = modelC.getCovaNumber()
        
        for i in range(self.ncovar):
            cova = modelC.getCova(i)
 
            self.sills += [cova.getSill().toTL()]
            self.invsill += [np.linalg.inv(self.sills[i])]
            self.cholsill += [np.linalg.cholesky(self.sills[i])]
            
            self.nvertex += [meshes[i].getNApices()]
            self.nvar = self.invsill[0].shape[0]
 
            modelMono = gl.Model.createFromParam(cova.getType(),
                                                 param = cova.getParam(),
                                                 range = cova.getRange(),
                                                 sills = 1)
            covatemp = modelMono.getCova(0)
            self.Qop += [createQ(meshes[i],covatemp)]
            self.temp += [gl.VectorDouble(np.zeros(shape=self.nvertex[i]))]  
        
        
        self.nvertextot = np.sum(self.nvertex)
        self.sizetot = self.nvertextot * self.nvar
        
    def evalDirect(self,inv):
        outv = np.zeros_like(inv)
        s=0
        for ii in range(self.ncovar):
            for i in range(self.nvar):
                self.Qop[ii].evalDirect(inv[(s+self.nvertex[ii]*i):(s+self.nvertex[ii]*(i+1))],self.temp[ii])
                for j in range(self.nvar):
                    outv[(s+self.nvertex[ii]*j):(s+self.nvertex[ii]*(j+1))]+=self.invsill[ii][i,j] * \
                        np.array(self.temp[ii].getVector())
            s+= self.nvertex[ii] * self.nvar
        return outv
    
    def evalSimulate(self,gauss):
        outv = np.zeros_like(gauss)
        s=0
        iadx = 0
        for ii in range(self.ncovar):
            nv = self.nvertex[ii]
            for i in range(self.nvar):
                self.Qop[ii].evalSimulate(
                    gl.VectorDouble(gauss[iadx:(iadx+nv)]),
                    self.temp[ii])
                iady = s
                for j in range(0,self.nvar):
                    outv[iady:(iady+nv)]+=self.cholsill[ii][j,i] *\
                        np.array(self.temp[ii].getVector())
                    iady+=nv
                iadx += nv
            s = iadx

        return outv
    
def computeError(res,ref,ndig =11):
    return round(np.max(np.abs(res-ref)),ndig)

def testCase(result,simulate = True, matrix = True,ndig = 11):
    if simulate:
        op = "simulation"
    else:
        op = "evalDirect"
    if matrix:
        case = "matrix"
    else:
        case = "matrix-free"

    
    if simulate and matrix :
        res = result["ressimumat"]
        ref = result["refsimumat"]
    if simulate and not matrix:
        res = result["ressimu"]
        ref = result["refsimu"]
    if not simulate and matrix:
        res = result["resevalmat"]
        ref = result["refevalmat"]
    if not simulate and not matrix:
        res = result["reseval"]
        ref = result["refeval"]

    vv = computeError(res,ref,ndig)
    message = f" Compare " + op + " with ref in the " + case + " case. Error = " + str(vv)
    print(message)


def computeAll(model,mesh,seed = 12344):

    s = sizetot(model,mesh)
    np.random.seed(seed)
    u = np.random.normal(size=s )
    
    modelCrMatrix = model.clone()
    modelCr       = model.clone()
    modelCsMatrix = model.clone()
    modelCs       = model.clone()

    poprMatrix =          PrecisionOpMulti(modelCrMatrix,mesh,True)
    popr       =          PrecisionOpMulti(modelCr      ,mesh,False)
    popsMatrix = gl.PrecisionOpMultiMatrix(modelCsMatrix,mesh)
    pops       =       gl.PrecisionOpMulti(modelCs      ,mesh)
    result = dict()
    result["ressimu"] = pops.evalSimulate(u)
    result["refsimu"] = popr.evalSimulate(u)

    result["reseval"] = pops.evalDirect(u)
    result["refeval"] = popr.evalDirect(u)

    result["ressimumat"] = popsMatrix.evalSimulate(u)
    result["refsimumat"] = poprMatrix.evalSimulate(u)

    result["resevalmat"] = popsMatrix.evalDirect(u)
    result["refevalmat"] = poprMatrix.evalDirect(u)

    return result

def sizetot(model,mesh):
    s = 0
    for i in mesh:
        s+=i.getNApices()
    return s * model.getVariableNumber()

def testAllSituations(model,mesh,verbose = True,ndig=11):
    
    result = computeAll(model,mesh)
    if verbose:
        for op in [False,True]:
            for mat in [False,True]:
                testCase(result,op,mat,ndig)

    vv = computeError(result["resevalmat"],result["reseval"],ndig)
    print(f" Compare evalDirect between the matrix case and the matrix-free case. Error = " +str(vv))
    print(f"---------------------------")

    return result


def test(nvar = 1, multistruct = True,ndig = 10):
  print(f"-----------------------------------------------")
  print(f"|      Case: nvar = " + str(nvar) + " multistruct = " + str(multistruct) + "      |")
  print(f"-----------------------------------------------")

  resultAll = dict()
  for i in ["Stat","Fake","Real"]:
    if i == "Real" and nvar == 3:
      continue
    print(f"No stationarity case: " + i)
    print(f"---------------------------")
    modelMulti,meshes =create(nvar = nvar,
                              multistruct = multistruct, 
                              nostatType = i,
                              nx1 = [4,4],
                              nx2 = [4,4])
    verbose = True
    if i == "Real":
      verbose = False
    resultAll[i] = testAllSituations(modelMulti,meshes,verbose,ndig)
  print(f"-------------------------------------------------------")
  print(f"Comparisons between fake non stationary and Stationary:") 
  print(f"-------------------------------------------------------")

  vv = computeError(resultAll["Fake"]["reseval"],
                    resultAll["Stat"]["reseval"],ndig)
  print(f" evalDirect in the matrix-free case. Error = " +str(vv))
  vv = computeError(resultAll["Fake"]["resevalmat"],
                    resultAll["Stat"]["resevalmat"],ndig)
  print(f" evalDirect in the matrix case. Error = " +str(vv))
  vv = computeError(resultAll["Fake"]["ressimu"],
                    resultAll["Stat"]["ressimu"],ndig)
  print(f" Simulation in the matrix-free case. Error = " +str(vv))
  vv = computeError(resultAll["Fake"]["ressimumat"],
                    resultAll["Stat"]["ressimumat"],ndig)
  print(f" Simulation in the matrix case. Error = " +str(vv))
  return resultAll

nvar = 1
multistruct =False

# %%
ndig = 8
for nvar in [1,2,3]:
    for multistruct in [False,True]:
        rr = test(nvar,multistruct,ndig)
