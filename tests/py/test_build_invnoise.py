# %%
import gstlearn as gl
import numpy as np

np.random.seed(1312)

# %%
# Nouvelle fonction renvoyant directement l'inverse de la matrice
def computeNew(dat, model, debug=False):
    mat = gl.buildInvNugget(dat, model)
    
    if debug:
        mat.display()
        
    return mat

# %%
#Fonction pour calculer l'inverse de la matrice de covariance correspondant au bruit. 
def computeRef(dat,model,debug=False):
    mat = model.evalCovMatrixSymmetric(db1=dat)

    #matinv = np.linalg.inv(mat.toTL())
    mat.invert()
    
    if debug:
        mat.display()

    return mat

# %%
#Renvoie un vecteur de bool indiquant pour chaque variable d'un échantillon donné
#si elle est observée (1) ou non (0)
def heterobool(dat,iech):
    names = dat.getNamesByLocator(gl.ELoc.Z)
    return ~np.isnan(dat.getValuesByNames(iech,names))

#Fonction qui s'applique au résultat de la fonction précédente et qui renvoit un indice
#d'hétérotopie par décomposition binaire.
#exemples:
# si pour iech, on a la variable 1, pas la variable 2, mais la variable 3
# on renvoit 1 * (2^0) + 0 * (2^1) + 1 * (2^2) = 5
# si toute les variables sont définies, on renvoit 1 + 2 + 4 = 7

def codeHetero(isdefined):
    return np.sum(isdefined *(2**np.arange(len(isdefined))))

# %%
#Vérifie si toutes les erreurs de mesure sont égales pour tous les échantillons d'une même variable.
#ou si il n'y a pas d'erreur de mesure.
def measurementErrorAllEquals(dat):
    ok = True
    nverr = dat.getLocNumber(gl.ELoc.V)
    if nverr==0:
        return ok
    nvar = dat.getLocNumber(gl.ELoc.Z)
    nech = dat.getSampleNumber()
    
    ok = 1
    for i in range(nvar):
        if i<nverr:
            err  = [dat.getLocVariable(gl.ELoc.V,iech,i) for iech in range(nech)]
            ok*=len(np.unique(err))==1
    return ok

# %%
#Renvoit un vecteur contenant l'erreur de mesure d'un échantillon pour
#chaque variable (0 si il n'y a pas d'erreur de mesure pour cette variable)
def getVerr(dat,iech):
    nvar = dat.getLocNumber(gl.ELoc.Z)
    verr = np.zeros(shape=nvar)
    nverr = dat.getLocNumber(gl.ELoc.V)
    for i in range(nvar):
        if i<nverr:
            verr[i] = dat.getLocVariable(gl.ELoc.V,iech,i)
    return verr

#add the measurements error of iech on the diagonal of sillMat
def update(dat,sillMat,iech):
    nvar = dat.getLocNumber(gl.ELoc.Z)
    verr = getVerr(dat,iech)
    for i in range(nvar):
         sillMat[i,i]+=verr[i]

# %%
#Compte le nombre de valeurs présentes (définies)  pour chaque variable
#puis définit l'indice de démarrage de l'indexation pour chaque variable
# (pour le remplissage de l'inverse de la matrice de covariance
def startingIndex(dat):
    names = dat.getNamesByLocator(gl.ELoc.Z)
    nvar = dat.getLocatorNumber(gl.ELoc.Z)
    count = np.array([0 for i in range(nvar)]) #std::vector<int> de zéros de taille nvar
    for iech in np.arange(dat.getSampleNumber(True)): #useSel = True
        count += heterobool(dat,iech)
    starind = 0 * count
    starind[1:] = np.cumsum(count)[:-1] #premier indice à 0 puis deuxième au nombre d'échantillons
                                        #actifs de la première variable, puis troisième à la somme des
                                        #échantillons actifs pour les 2 premières variables ...
    return starind

# %%
#Calcule l'inverse de la matrice de covariance
#Cas Isotopique + Stationnaire + toutes les erreurs de mesure égales
def IsoStat(dat,model):

    #On calcule l'inverse de la matrice de Sills
    sillmat = model.getCova(0).getSill().toTL()
    update(dat,sillmat,0)
    invAic = np.linalg.inv(sillmat) 
    
    #On fournit cette matrice à un modèle pépitique. C'est juste un trick pour ne 
    #pas avoir à recoder la construction de la matrice block avec des diagonales
    #mais ce n'est sans doute pas le plus efficace car ça nécessite de calculer plein de distances inutilement
    #à moins que ce ne soit codé intelligemment dans gstlearn dans le cas de la covariance.
    #D'ailleurs, pour les applications large échelle (matrix free), il me semble qu'on devrait juste construire l'opérateur
    #pour éviter de stocker une immense matrice avec duplication de la même valeur des millions de fois
    #Autre commentaire: on devrait pouvoir passer une matrice de sills plutôt qu'un vecteur double.
    #Pour éviter de tout casser (notamment les typemaps qui permettent de fournir un float dans le cas
    #monovariable), peut-être qu'un createMultivariateFromParam serait le bienvenu. Ou alors voir avec Fabien 
    #si on peut faire un typemap qui convertit les floats en matrice 1x1.
    modeltrick = gl.Model.createFromParam(gl.ECov.NUGGET,sills = invAic.reshape(-1))
    
    #J'ai utilisé une matrice Sparse. Ca fonctionne bien (je ne pensais pas qu'on utiliserait si vite 
    #cette fonction).
    datc = dat.clone()
    datc.clearLocators(gl.ELoc.V)
    result = modeltrick.evalCovMatrixSparse(datc)

    return result

# %%
#Ajoute les triplets pour un iech donné et incrémente le compteur d'indexation
def patch(nvar,isdefined,triplet,starind,count,tempmat):
    indi = 0
    for ivar in range(nvar):
        indj=0
        if isdefined[ivar]:
            for jvar in range(nvar):
                if isdefined[jvar]:
                    triplet.add(starind[ivar]+count[ivar],
                                starind[jvar]+count[jvar],
                                tempmat[indi,:][indj])
                    indj+=1
            indi+=1 
    count += isdefined


# %%
#Calcule l'inverse de la matrice de covariance
#quand on n'est pas dans le cas Isotopique + Stationnaire + toutes les erreurs de mesure égales
def GeneralCase(dat,model):
    nvar = len(dat.getNamesByLocator(gl.ELoc.Z))

    #std::vector<MatrixSymmetric>
    matrices = [None for i in range(2**nvar)]

    starind = startingIndex(dat)

    #Compteur pour chaque variable pour l'hétérotopie
    count = np.array([0 for i in range(nvar)])

    triplet = gl.NF_Triplet()

    sillMat = model.getCova(0).getSill().toTL()

    if measurementErrorAllEquals(dat):
       update(dat,sillMat,0) ##### ATTENTION, je n'ai pas traité le cas où le premier échantillon
                             ##### ne serait pas observé pour une variable et n'aurait potentiellement
                             ##### pas d'erreur de mesure. Quoiqu'il en soit, dans ce cas, je ne sais pas
                             ##### ce que répond la fonction measurementErrorAllEquals
        
    for iech in np.arange(dat.getSampleNumber()):
        isdefined = heterobool(dat,iech) #privilegier un proto heterobool(const &db, int, &vector<int>)
        indexdefined = np.where(isdefined)[0]

        if ~model.getCova(0).isNoStat() and measurementErrorAllEquals(dat):
            hetero = codeHetero(isdefined)
            #Si on n'a pas encore trouvé cette configuration d'hétérotopie,
            #on inverse la sous matrice des AIC correspondante et on la stocke
            if matrices[hetero] is None: #verifier si la matrice est de dimension 0 x 0
                matrices[hetero] = np.linalg.inv(sillMat[indexdefined,:][:,indexdefined])
            tempmat = matrices[hetero] #pointeur
        else:
            sillMatc = sillMat.copy()
            update(dat,sillMatc,iech)
            tempmat= np.linalg.inv(sillMatc[indexdefined,:][:,indexdefined])
            
        patch(nvar,isdefined,triplet,starind,count,tempmat)

    return gl.MatrixSparse.createFromTriplet(triplet)

# %%
#Calcule l'inverse de la matrice de covariance
def computeInvNoise(dat,model,verbose=True, debug=False):
  
    #Si le modèle est stationnaire et qu'il n'y a pas d'hétérotopie
    if not model.getCova(0).isNoStat() and dat.isAllIsotopic() and measurementErrorAllEquals(dat):
        if verbose:
            print("---------------------------------------------------------")
            print(f"Isotopic and Stationary and Measurement error identical:")
            print("---------------------------------------------------------")            
        result = IsoStat(dat,model)
        
    else:
        if verbose:
            print(f"General Case")
        result = GeneralCase(dat,model)
    
    if debug:
        result.display()
    return result

# %%
#Créé une Db pour les tests
def createDb(ndat,hetero,measurement_error,measurement_error_all_equal,nostat=False):
    
    if ndat < 10:
        gl.messerr(f"Data File must contain at least 10 samples")
        return None
    
    dat = gl.Db.create()
    dat["x1"] = np.random.normal(size=ndat)
    dat["x2"] = np.random.normal(size=ndat)
    dat["z1"] = np.random.normal(size=ndat)
    dat["z2"] = np.random.normal(size=ndat)
    dat["z3"] = np.random.normal(size=ndat)

    if measurement_error:
        dat["err1"] = np.random.uniform(size=ndat)
        if measurement_error_all_equal:
            dat["err1"] = 0 * dat["err1"] + 0.5
        dat.setLocator("err1",gl.ELoc.V)

    if nostat:
        dat["v11"] =  np.abs(np.random.normal(size=ndat)) 
        dat["v22"] =  np.abs(np.random.normal(size=ndat)) 
        dat["v12"] =  np.abs(np.random.uniform(size=ndat))*np.sqrt(dat["v11"]*dat["v22"])
        model.getCova(0).makeSillNoStatDb("v11",0,0,dat)
        model.getCova(0).makeSillNoStatDb("v12",0,1,dat)
        model.getCova(0).makeSillNoStatDb("v22",1,1,dat)


    if hetero:
        # Ranks of undefined variables have been reduced to allow smaller data set
        dat.setValue("z2",2,gl.TEST)
        dat.setValue("z2",4,gl.TEST)
        dat.setValue("z2",8,gl.TEST)
        
        dat.setValue("z3",3,gl.TEST)
        dat.setValue("z3",4,gl.TEST)
        dat.setValue("z3",6,gl.TEST)
        dat.setValue("z3",8,gl.TEST)

    dat.setLocators(["x*"],gl.ELoc.X)
    dat.setLocators(["z*"],gl.ELoc.Z)
    return dat

# %%
#Fonction de test
def testInvNoise(dat,model, debug=False):
    if dat == None:
        return
    
    ref = computeRef(dat,model, debug=debug).toTL()
    
    refnew = computeNew(dat, model, debug=debug).toTL()
    
    result = computeInvNoise(dat,model,True, debug=debug).toTL()
    
    #Pour la comparaison, je remets le résultat en dense
    err1 = np.sum(np.abs(ref - refnew))
    if err1 > 1e-13:
        print(f"Pb between ref and new")
        
    error = np.sum(np.abs(refnew.toarray()-ref))
    
    if error < 1e-13:
        print(f"Error = ",1e-13)
    else: 
        print(error)

# %%
ndat= 10

sills = np.array([[1,0.2,0.1], [0.2,2,0.3], [0.1,0.3,3]])
model = gl.Model.createFromParam(gl.ECov.NUGGET,sills = sills.reshape(-1))
debug = False



# %%
print("---------------------------------")
print(f"Hetero without measurement error")
print("---------------------------------")
dat = createDb(ndat,True,False,False)
testInvNoise(dat,model,debug)

# %%
print("-----------------------------------------------------------------------")
print(f"Hetero with measurement errors (not the same for all the observations)")
print("-----------------------------------------------------------------------")
dat = createDb(ndat,True,True,False)
testInvNoise(dat,model,debug)

# %%
print("----------------------------------------------------------------")
print(f"Hetero with the same measurement error for all the observations")
print("----------------------------------------------------------------")
dat = createDb(ndat,True,True,True)
testInvNoise(dat,model,debug)

# %%
print("------------------------------")
print(f"Iso without measurement error")
print("------------------------------")
dat = createDb(ndat,False,False,False)
testInvNoise(dat,model,debug)


# %%
print("----------------------------------------------------------------------")
print(f"Iso,  with measurement errors (not the same for all the observations)")
print("----------------------------------------------------------------------")
dat = createDb(ndat,False,True,False)
testInvNoise(dat,model,debug)

# %%
print("-------------------------------------------------------------")
print(f"Iso,with the same measurement error for all the observations")
print("-------------------------------------------------------------")
dat = createDb(ndat,False,True,True)
testInvNoise(dat,model,debug)

# %%
print("---------------------------------------------------------------------")
print(f"Iso,with the same measurement error for all the observations, nostat")
print("---------------------------------------------------------------------")
dat = createDb(ndat,False,True,True,True)
testInvNoise(dat,model,debug)

# %%



