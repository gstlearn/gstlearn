################################################################################
#                                                                              #
#                         gstlearn Python package                              #
# Copyright (c) (2023) MINES Paris / ARMINES                                   #
# Authors: gstlearn Team                                                       #
# License: BSD 3-clause                                                        #
#                                                                              #
################################################################################
import matplotlib.pyplot     as plt
import matplotlib.patches    as ptc
import matplotlib.transforms as transform
import matplotlib.colors     as mcolors
import numpy                 as np
import numpy.ma              as ma
import gstlearn              as gl
import gstlearn.plot         as gp

from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy                   import shape
from pandas.io               import orc
from plotly.matplotlylib     import mpltools
import math
from plotly.validators.layout.scene import aspectratio
from matplotlib.pyplot import axes

#Set of global values
defaultDims = [[5,5], [8,8]]
defaultXlim = [ None, None ]
defaultYlim = [ None, None ]
defaultAspect = [ 'auto', 1 ]

def setDefaultGeographic(dims=None, xlim=None, ylim=None, aspect=None):
    '''
    Set the default values for all *Geographical* figures in the Graphic Environment.
    dims: Dimensions of each figure
    xlim: Bounds along the horizontal axis (otherwise: no bound is defined)
    ylim: Bounds along the vertical axis (otherwise: no bound is defined)
    aspect: Ratio between dimensions along bith axes (Y/X)
    
    Remark:
    - to reset the default value for 'xlim' (resp. 'ylim'): set xlim='reset'
    '''
    __setDefaultInternal(1, dims=dims, xlim=xlim, ylim=ylim, aspect=aspect)
    
def setDefault(dims=None, xlim=None, ylim=None, aspect=None):
    '''
    Set the default values for all *non-Geographical* figures in the Graphic Environment.
    dims: Dimensions of each figure
    xlim: Bounds along the horizontal axis (otherwise: no bound is defined)
    ylim: Bounds along the vertical axis (otherwise: no bound is defined)
    aspect: Ratio between dimensions along bith axes (Y/X)
    '''
    __setDefaultInternal(0, dims=dims, xlim=xlim, ylim=ylim, aspect=aspect)

def __setDefaultInternal(mode, dims=None, xlim=None, ylim=None, aspect=None):
    global defaultDims
    global defaultXlim
    global defaultYlim
    global defaultAspect
        
    if dims is not None:
        defaultDims[mode] = dims
    if xlim is not None:
        if __isArray(xlim, 2):
            defaultXlim[mode] = xlim
        elif xlim == 'reset':
            defaultXlim[mode] = None
    if ylim is not None:
        if __isArray(ylim, 2):
            defaultYlim[mode] = ylim
        elif ylim == 'reset':
            defaultYlim[mode] = None
    if aspect is not None:
        defaultAspect[mode] = aspect

def printDefault():
    ''' 
    Print the variables defined in the Environment
    - for the non-Geographical figures
    - for the Geographical figures
    '''
    for mode in range(2):
        if mode == 0:
            print("Non geographical defaults:")
        else:
            print("Geographical defaults:")
            
        if defaultDims[mode] is not None:
            print("- Figure dimensions =", defaultDims[mode])
        else:
            print("- Figure dimensions (not defined)")
        if defaultXlim[mode] is not None:
            print("- Limits along X =",defaultXlim[mode])
        else:
            print("- Limits along X (not defined)")
        if defaultYlim[mode] is not None:
            print("- Limits along Y =",defaultYlim[mode])
        else:
            print("- Limits along Y (not defined)")
        if defaultAspect[mode] is not None:
            print("- Aspect =",defaultAspect[mode])
        else:
            print("- Aspect (not defined)")
        
def getColorMap(n, name='gist_rainbow'):
    '''
    Returns a function that maps each index in 0, 1, ..., n-1 to a distinct RGB color
    
    n: requested number of different colors
    name: this argument must be a standard mpl colormap name.
    '''
    return plt.cm.get_cmap(name, n)
    
def __selectItems(nvalues, sitem=-1):
    outs = range(0, nvalues)
    nout = nvalues
    if sitem >= 0:
        outs = range(sitem, sitem+1)
        nout = 1
    return outs, nout

def __isNotCorrect(object, types):
    if object is None:
        print("Argument 'object' must be provided")
        return True
    filetype = type(object).__name__
    if filetype in types:
        return False
    
    print("Argument 'object' (",filetype,") must be a valid type among",types)
    return True

def __defaultVariable(db, name):
    if name is None:
        if db.getLocNumber(gl.ELoc.Z) > 0:
            name = db.getNameByLocator(gl.ELoc.Z,0)
        else : # if no Z locator, choose the last field
            name = db.getLastName()
    else:
        if db.getUID(name) < 0:
            name = db.getLastName()
    return name

def geometry(ax, dims=None, xlim=None, ylim=None, aspect=None):
    '''
    Set the default values for the geometric parameters for one or a set of Axes
    
    ax: matplotlib.Axes or numpy.array of matplotlib.Axes or matplotlib.Figure
    dims: Extension of graphic Axes
    xlim: Range of values along the X-axis
    ylim: Range of values along the Y-axis
    aspect: Y/X ratio
    
    Remark: When 'ax' designates a set of Axes, parameters are applied to all of them.
    '''
    if dims is not None:
        if type(ax) == plt.Figure:
            ax_list = ax.get_axes()
            gspec = ax_list[0].get_gridspec()
            ax_list[0].figure.set_size_inches(dims[0]*gspec.nrows, 
                                              dims[1]*gspec.ncols)
        elif __isArray(dims, 2):
            if type(ax) == np.ndarray:
                ax[0,0].figure.set_size_inches(dims[0]*ax.shape[0], 
                                               dims[1]*ax.shape[1])
            else:
                ax.figure.set_size_inches(dims[0], dims[1])
        else:
            print("'dims' should be [a,b]. Ignored")
        
    if xlim is not None:
        if type(ax) == plt.Figure:
            ax_list = ax.get_axes()
            for ax in ax_list:
                ax.set_xlim(left = xlim[0], right = xlim[1])
        elif __isArray(xlim, 2):
            if type(ax) == np.ndarray:
                ax[ix,iy].set_xlim(left = xlim[0], right = xlim[1])
            else:
                ax.set_xlim(left = xlim[0], right = xlim[1])
        else:
            print("'xlim' should be [a,b] or [None,b] or [a,None]. Ignored")
    
    if ylim is not None:
        if type(ax) == plt.Figure:
            ax_list = ax.get_axes()
            for ax in ax_list:
                ax.set_ylim(left = ylim[0], right = ylim[1])
        if __isArray(ylim, 2):
            if type(ax) == np.ndarray:
                for ix in range(ax.shape[0]):
                    for iy in range(ax.shape[1]):
                        ax[ix,iy].set_ylim(bottom = ylim[0], top = ylim[1])
            else:
                ax.set_ylim(bottom = ylim[0], top = ylim[1])
        else:
           print("'ylim' should be [a,b] or [None,b] or [a,None]. Ignored")
        
    if aspect is not None:
        if type(ax) == plt.Figure:
            ax_list = ax.get_axes()
            for ax in ax_list:
                ax.set_aspect(aspect)
        elif type(ax) == np.ndarray:
            for ix in range(ax.shape[0]):
                for iy in range(ax.shape[1]):
                    ax[ix,iy].set_aspect(aspect)
        else:
            ax.set_aspect(aspect)
            
# This function should return 'ax'. Has been skipped for compatibility with previous version           
#    return ax

def decoration(ax, xlabel=None, ylabel=None, title=None, **kwargs):
    '''
    Add the decoration to a figure.
    
    Parameters
    ----------
    ax: matplotlib.Axes or numpy array of matplotlib.Axes or Matplotlib.Figure
    xlabel: label along the horizontal axis
    ylabel: label along the vertical axis
    title: title contents (for the main for a collection of Axes)
    '''
    if title is not None:
        if type(ax) == plt.Figure:
            ax.suptitle(title, **kwargs)
        elif type(ax) == np.ndarray:
            ax[0,0].figure.suptitle(title, **kwargs)
        else:
            ax.set_title(title, **kwargs)

    if xlabel is not None:
        if type(ax) == plt.Figure:
            print("decoration() for xlabel cannot be used when 'ax' is an array. Ignored")
        elif type(ax) == np.ndarray:
            print("decoration() for xlabel cannot be used when 'ax' is an array. Ignored")
        else:
            ax.set_xlabel(xlabel)
    if ylabel is not None:
        if type(ax) == plt.Figure:
            print("decoration() for ylabel cannot be used when 'ax' is an array. Ignored")
        elif type(ax) == np.ndarray:
            print("decoration() for ylabel cannot be used when 'ax' is an array. Ignored")
        else:
            ax.set_ylabel(ylabel)
 
def __initGeneric(mode=0, nx=1, ny=1, sharex=False, sharey=False, figsize=None):
    ''' Creates a new Geographic figure (possibly containing several subplots)
    
        Parameters
        ----------
        nx, ny:     Number of subplots along X and Y
        sharex, sharey: If the subplots should all share respectively X and Y axis
        figsize: Vector giving the dimension of the picture
        
        Details
        -------
        The use of 'figsize' overwrites the default values stored in geographical
        and non-geographical environments. In the case of multi-figures, it corresponds
        to the overall dimension. The dimension of each figure is derived internally.
        
        Returns
        -------
        ax description
    '''
    if len(plt.get_fignums()) == 0:
        
        # Axes is None and no Figure already exists. Create it
        fig, ax = plt.subplots(nx, ny, squeeze=False, sharex=sharex, sharey=sharey)
        
        if __isArray(ax, 1):
            ax = ax[0,0]

        # Apply the Global Geometry parameters (when defined)
        
        if figsize is None:
            locdims = defaultDims[mode]
        else:
            locdims = [figsize[0] / nx, figsize[1] / ny]
            
        geometry(ax,
                 dims = locdims, 
                 xlim = defaultXlim[mode], 
                 ylim = defaultYlim[mode], 
                 aspect = defaultAspect[mode])
        
    else:
        # Axes is None but a figure already exists, return the (last) Axes of Figure
        fig = plt.gcf()
        ax = plt.gca()
    
    if __isArray(ax, 1):
        ax = ax[0,0]

    return fig, ax
    
def initGeographic(nx=1, ny=1, sharex=False, sharey=False, figsize=None):
    return __initGeneric(1, nx=nx, ny=ny, sharex=sharex, sharey=sharey, figsize=figsize)

def init(nx=1, ny=1, sharex=False, sharey=False, figsize=None):
    return __initGeneric(0, nx=nx, ny=ny, sharex=sharex, sharey=sharey, figsize=figsize)

def __getNewAxes(ax=None, mode=0, nx=1, ny=1, sharex=False, sharey=False):
    ''' Creates a new figure (possibly containing multiple subplots)
    
        Parameters
        ----------
        ax: Axes description (or array of them). See remarks.
        mode: Type of Global value (used for default assignment). 0 for Non-Geographical; 1 for Geographical
        nx, ny:     Number of subplots along X and Y
        sharex, sharey: If the subplots should all share respectively X and Y axis
        
        Returns
        -------
        ax description
        
        Remarks
        -------
        
        If 'ax' does not exist (input argument), a new figure is created. 
        Otherwise, the input argument is simply returned.
    '''
    if ax is None:
        _, ax = __initGeneric(mode=mode, nx=nx, ny=ny, sharex=sharex, sharey=sharey)
    
    if __isArray(ax, 1):
        ax = ax[0,0]
    
    return ax

def __isArray(tab, ndim=None):
    '''
    Check if the input argument is an array (of dimension 'ndim' when defined) or a scalar
    tab:  Argument to be checked
    ndim: Requested dimension. The dimension check is not performed if None 
    '''
    if not hasattr(tab, "__len__"):
        return False
    
    if not isinstance(tab, np.ndarray):
        tab = np.asarray(tab)
        
    if (ndim is not None) and (tab.size is not ndim):
        return False
    
    return True

def __addColorbar(im, ax, legendName = None):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, ax=ax, cax=cax)
    if legendName is not None:
        cbar.ax.set_title(legendName)
    return cbar

def __getDefinedValues(db, name, posX=0, posY=1, corner=None, useSel=True, 
                       compress=False, asGrid=True, flagConvertNanToZero=False):

    if db.isGrid() and asGrid:
        if corner is None:
            corner = np.zeros(db.getNDim())
        
        if db.getNDim() == 1:
            tab = db.getColumn(name, useSel, False)
        else:
            tab = db.getOneSlice(name, posX, posY, corner, useSel)
    else:
        tab = db.getColumn(name, useSel, True)
    tab = np.array(tab).transpose()

    if flagConvertNanToZero:
        tab[np.isnan(tab)] = 0
    
    if compress:
        tab = tab[np.logical_not(np.isnan(tab))]
        
    return tab

def varioElem(vario, ivar=0, jvar=0, *args, **kwargs):
    """
    Plot a single experimental variogram (one direction and fixed pair of variable(s)).
    
    Parameters
    ----------
    ax: matplotlib.Axes (necessary when used as a method of the class)
    vario : experimental variogram to be represented (gstlearn.Vario).
    ivar, jvar : Indices of the variables for the variogram to be represented (the default is 0).
    idir : Index of the direction of the variogram to be represented (the default is 0).
    hmax : Maximum distance to be represented.
    showPairs : Flag for annotating the number of pairs for each lag on the plot (the default is False).
    varColor : color of the horizontal line representing the sill (the default is 'black').
    varLinestyle : linestyle of the horizontal line representing the sill (the default is 'dashed').
    label : Label to be drawn (constructed if not provided)
    flagDrawVariance : Flag to add the variance (default is True)    
    flagLabelDir : Encode the direction in the label (when constructed)
    flagLegend : Flag to display the axes legend.
    **kwargs : arguments passed to matplotlib.pyplot.plot

    Returns
    -------
    ax : axes where the variogram is represented
    """
    ax = __getNewAxes(None, 0)
    return __ax_varioElem(ax, vario, ivar=ivar, jvar=jvar, *args, **kwargs)

def __ax_varioElem(ax, vario, ivar=0, jvar=0, idir=0, hmax=None, showPairs = False,
                   varColor='black', varLinestyle='dashed', 
                   flagDrawVariance = True, flagLabelDir=False, flagLegend=False, 
                   label=None, **kwargs):
    if label is None:
        if flagLabelDir:
            if vario.isDefinedForGrid():
                label = "vario grid={}".format(vario.getGrincrs(idir))
            else:
                angles = gl.GeometryHelper.rotationGetAngles(vario.getCodirs(idir),True)
                label = "vario dir={}".format(np.round(angles,3))
        else:
            label = "vario"
    
    # Plotting the experimental variogram
    gg = vario.getGgVec(idir,ivar,jvar)
    hh = vario.getHhVec(idir,ivar,jvar)
    if len(hh) == 0:
        return None
    
    # Representing the variogram
    res = ax.plot(hh, gg, label = label, **kwargs)
    
    # Maximum distance to be represented
    if hmax is None:
        hmax = np.nanmax(hh)
    
    # Adding the Y=0 axis in multivariate case
    if ivar != jvar:
        ax.hlines(0, 0, hmax, colors="black", linewidth=0.5)

    # Drawing the variance-covariance reference line (optional)
    if flagDrawVariance:
        ax.hlines(vario.getVar(ivar,jvar), 0, hmax, varColor, varLinestyle)
    
    # Representing the number of pairs (optional)
    if showPairs:
        pairs = vario.getSwVec(idir,ivar,jvar)
        for i in range(len(hh)):
            ax.annotate(str(int(pairs[i])), (hh[i],gg[i]), xytext=(0,5), xycoords = 'data',
                        textcoords = 'offset points', ha='center')
    
    # Displaying the Legend (optional)
    if flagLegend:
        ax.legend()
    
    # Hard code the lower bounds
    if vario.drawOnlyPositiveX(ivar, jvar):
        ax.set_xlim(left=0)
    if vario.drawOnlyPositiveY(ivar, jvar):
        ax.set_ylim(bottom=0)
        
    return res

def varmod(vario, model=None, ivar=-1, jvar=-1, axsOld=None, *args, **kwargs):
    """
    Construct a figure for plotting experimental variogram(s) and model.
    
    Parameters
    ----------
    vario : experimental variogram to be represented
    model : optional, variogram model
    ivar, jvar : Indices of the variables for the variogram to be represented. If -1 (default), all 
                 variables are selected and all the simple and crossed variograms are represented.
    idir : Index of the direction of the variogram to be represented. If -1 (default) all available
           directions are selected and multi-directional variograms are represented.
    flagDrawVariance : Flag to add the variance (default is True)  
    varioLinestyle: Linestyle for representing the experimental variogram
    modelLinestyle: Linestyle for representing the Model
    varColor, varLinestyle: parameters for representing variance-covariance line
    envColor, envLinestyle: parameters for representing coregionalization envelop
    nh : number of points between 0 and hmax where the model variogram is calculated (default is 100).
    hmax : Maximum distance to be represented.
    cmap : Optional Color scale
    flagLegend : Flag to display the axes legend.
    axsOld : Reference for the plot(s) within the figure. If None (default),
          it creates a new figure (with multiple axes for multivariate variograms).

    **kwargs : arguments passed to matplotlib.pyplot.plot for all variograms plotted (not models!)
    """
    nvar = vario.getVariableNumber()
    ivarUtil, ivarN = __selectItems(nvar, ivar)
    jvarUtil, jvarN = __selectItems(nvar, jvar)
    axs = __getNewAxes(axsOld, 0, nx=ivarN, ny=jvarN)

    return __ax_varmod(axs, vario=vario, model=model, ivar=ivar, jvar=jvar, 
                       *args, **kwargs)

def __ax_varmod(axs, vario, model=None, ivar=-1, jvar=-1, idir=-1,
                nh = 100, hmax = None, showPairs=False, asCov=False, 
                flagDrawVariance = True,
                varioLinestyle = 'dashed', modelLinestyle = 'solid',
                varColor='black', varLinestyle="dotted",
                envColor='black', envLinestyle="dotted",
                cmap=None, flagLegend=False,
                **kwargs):
    if __isNotCorrect(object=vario, types=["Vario"]):
        return None

    color_in_kwargs = 'color' in kwargs
    
    if hmax is None:
        hmax = vario.getHmax(ivar, jvar, idir)
        
    ndir = vario.getDirectionNumber()
    nvar = vario.getVariableNumber()
    cols = getColorMap(ndir,cmap)
    
    ndirUtil, ivarD = __selectItems(ndir, idir)
    ivarUtil, ivarN = __selectItems(nvar, ivar)
    jvarUtil, jvarN = __selectItems(nvar, jvar)
    
    # Check that the number of subplots in the Figure is correct
    if type(axs) == plt.Figure:
        ax_list = axs.get_axes()
        if len(ax_list) != ivarN * jvarN:
            print("Argument 'axs' does not contain the correct number of subplots")
            exit()
            
    if ndir > 1:
        flagLabelDir = True
    else:
        flagLabelDir = False
        
    # Loop on the variables 
    iaxes = 0 
    for iv in ivarUtil:
        for jv in jvarUtil:
            
            if type(axs) == plt.Figure:
                ax = ax_list[iaxes]
                iaxes = iaxes + 1
            elif type(axs) == np.ndarray:
                ax = axs[iv,jv]
            else:
                ax = axs

            if iv < jv:
                ax.set_visible(False)
                continue;
            
            for idirUtil in ndirUtil:
                if not color_in_kwargs:
                    kwargs.update({'color':cols(idirUtil)})
                if varioLinestyle is not None:
                    kwargs.update({'linestyle': varioLinestyle})
                
                __ax_varioElem(ax, vario, iv, jv, idirUtil, 
                               showPairs=showPairs, hmax=hmax,
                               flagDrawVariance = flagDrawVariance,
                               varColor=varColor, varLinestyle=varLinestyle,  
                               flagLabelDir=flagLabelDir, flagLegend=flagLegend, 
                               **kwargs)

                # Plotting the Model (optional)
                if model is not None:
                    codir = vario.getCodirs(idirUtil)
                    if modelLinestyle is not None:
                        kwargs.update({'linestyle': modelLinestyle})
                    __ax_modelElem(ax, model, ivar=iv, jvar=jv, codir=codir, 
                                   hmax=hmax, nh=nh, asCov=asCov,
                                   envColor=envColor, envLinestyle=envLinestyle, 
                                   flagLabelDir=flagLabelDir, flagLegend=flagLegend, 
                                   **kwargs)

            ax.autoscale(True)
            
            if vario.drawOnlyPositiveX(iv, jv):
                ax.set_xlim(left=0)
            if vario.drawOnlyPositiveY(iv, jv):
                ax.set_ylim(bottom=0)
    
    return axs

def variogram(vario, ivar=0, jvar=0, axsOld=None, *args, **kwargs):
    """
    Plot experimental variogram(s) (can be multidirectional and multivariable or selected ones).
    
    Parameters
    ----------
    axs : matplotlib.Axes or matplotlib.Figure
    vario : experimental variogram to be represented (gstlearn.Vario).
    ivar, jvar : Indices of the variables for the variogram to be represented. If -1 (default), all 
                 variables are selected and all the simple and crossed variograms are represented.
    idir : Index of the direction of the variogram to be represented. If -1 (default) all available
           directions are selected and multidirectional variograms are represented.
    varColor, varLinestyle : parameters for representing variance-covariance line
    hmax : Maximum distance to be represented.
    cmap : Optional Color scale
    flagLegend : Flag to display the axes legend.
    axs : Reference for the plot(s) within the figure. If None (default),
          it creates a new figure (with multiple axes for multivariate variograms).
    **kwargs : arguments passed to matplotlib.pyplot.plot for all variograms plotted

    Returns
    -------
    axs : axes where the variograms are represented
    """
    nvar = vario.getVariableNumber()
    ivarUtil, ivarN = __selectItems(nvar, ivar)
    jvarUtil, jvarN = __selectItems(nvar, jvar)
    axs = __getNewAxes(axsOld, 0, nx=ivarN, ny=jvarN)
    
    return __ax_variogram(axs, vario, ivar=ivar, jvar=jvar, *args, **kwargs)
    
def __ax_variogram(axs, vario, ivar=0, jvar=0, idir=0,
                   varColor='black', varLinestyle='dotted', hmax=None,  
                   cmap = None, flagLegend=False, 
                   *args, **kwargs):
    return __ax_varmod(axs, vario, ivar=ivar, jvar=jvar, idir=idir, 
                       varColor=varColor, varLinestyle=varLinestyle, 
                       hmax=hmax, cmap=cmap, flagLegend=flagLegend, 
                       *args, **kwargs)

def modelElem(modelobj, ivar=0, jvar=0, *args, **kwargs):
    """
    Construct a Layer for plotting a model
    
    Parameters
    ----------
    ax: matplotlib.Axes
    modelobj : variogram model to be represented (gstlearn.Model).
    ivar, jvar : Indices of the variables for the variogram to be represented (the default is 0).
    codir : Vector of the direction of the variogram to be represented. The default is the unit 
            vector in the first space dimension.
    vario, idir: Vario information used to set the direction (when codir is not provided)
    envColor, envLinestyle : color and linestyle for coregionalization envelop 
    nh : number of points between 0 and hmax where the model variogram is calculated (default is 100).
    hmax : Maximum distance to be represented. By default: 3 times the maximum range of the
           basic structures, or 1 if no range is defined.
    asCov : Present the Model as a Covariance (rather than as a Variogram)
    label: Label displayed in the Legend (constructed if not provided)
    flagLabelDir : Add the direction vector to the label (if constructed)
    flagEnvelop: Represent the coregionalization envelop (in multivariate case only)
    flagLegend : Flag to display the axes legend.
    """
    ax = __getNewAxes(None, 0)
    return __ax_modelElem(ax, modelobj, ivar=ivar, jvar=jvar, *args, **kwargs)

def __ax_modelElem(ax, modelobj, ivar=0, jvar=0, codir=None, vario=None, idir=0,
                   nh = 100, hmax = None, asCov=False,
                   envColor='black', envLinestyle='dashed',
                   label=None, flagLabelDir=False, flagEnvelop = True, flagLegend=False, 
                   **kwargs):
    if __isNotCorrect(object=modelobj, types=["Model"]):
        return None

    if codir is None:
        if vario is None:
            codir = [0] * modelobj.getDimensionNumber()
            codir[0] = 1
        else:
            codir = vario.getCodirs(idir)
            
    # if hmax not specified = 3*maximum range of the model's basic structures
    if hmax is None:
        hmax = 0
        for icova in range(modelobj.getCovaNumber()):
            range_max = np.max(modelobj.getCova(icova).getRanges())
            if 3*range_max > hmax:
                hmax = 3*range_max
    if hmax == 0: # if the model has no range defined
        hmax = 1
            
    if label is None:
        if flagLabelDir:
            angles = gl.GeometryHelper.rotationGetAngles(codir,True)
            label = "model dir={}".format(np.round(angles,3))
        else:
            label = "model"

    istart = 0
    for i in range(modelobj.getCovaNumber()):
        if modelobj.getCovName(i) == 'Nugget Effect':
            istart = 1 # do not plot the first lag (h=0) for nugget effect (discontinuity)
     
    # Represent the Model 
    hh = np.linspace(0, hmax, nh+1)
    mode = gl.CovCalcMode()
    mode.setAsVario(not asCov)
    gg = modelobj.sample(hh, codir, ivar, jvar, mode)
    res = ax.plot(hh[istart:], gg[istart:], label=label, **kwargs)
    
    # Represent the coregionalization envelop (optional)
    if ivar != jvar and flagEnvelop:
        ggp = modelobj.envelop(hh, ivar, jvar, +1, codir, mode)
        ax.plot(hh[istart:], ggp[istart:], c = envColor, linestyle = envLinestyle)
        ggm = modelobj.envelop(hh, ivar, jvar, -1, codir,mode)
        ax.plot(hh[istart:], ggm[istart:], c = envColor, linestyle = envLinestyle)
    
    # Draw the Legend (optional)
    if flagLegend:
        ax.legend()
        
    return res

def model(modelobj, *args, **kwargs):
    ax = __getNewAxes(None, 0)
    return __ax_model(ax, modelobj, *args, **kwargs)
    
def __ax_model(ax, modelobj = None, **kwargs):
    __ax_modelElem(ax, modelobj = modelobj, **kwargs)
    
    return ax

def __readCoorPoint(db, nameCoorX=None, nameCoorY=None, useSel=True, posX=0, posY=1):
    
    # Extracting coordinates
    if nameCoorX is not None:
        tabx = db.getColumn(nameCoorX, useSel)
    else:
        if db.getNDim() > 0:
            tabx = db.getCoordinates(posX,useSel)
            
    if nameCoorY is not None:
        taby = db.getColumn(nameCoorY, useSel)
    else:
        if db.getNDim() > 1:
            taby = db.getCoordinates(posY,useSel)
    
    if len(tabx) <= 0 or len(taby) <= 0:
        return None
    
    return tabx, taby
    
def symbol(db, nameColor=None, nameSize=None, *args, **kwargs):
    ax = __getNewAxes(None, 0)
    return __ax_symbol(ax, db, nameColor=nameColor, nameSize=nameSize, 
                       *args, **kwargs)
    
def __ax_symbol(ax, db, nameColor=None, nameSize=None, 
                nameCoorX=None, nameCoorY=None, useSel=True, 
                c='r', s=20, sizmin=10, sizmax=200, flagAbsSize=False, flagCst=False,
                flagLegendColor=False, flagLegendSize=False,
                legendNameColor=None, legendNameSize=None, posX=0, posY=1, **kwargs):
    '''
    Construct a Layer for plotting a point data base, with optional color and size variables
    
    ax: matplotlib.Axes (necessary when used as a method of the class)
    db: Db containing the variable to be plotted
    nameColor: Name of the variable containing the color per sample
    nameSize: Name of the variable containing the size per sample
    nameCoorX: Name of the variable standing for X coordinate 
    nameCoorY: Name of the variable standing for Y coordinate 
    useSel : Boolean to indicate if the selection has to be considered
    c: Constant color (used if 'nameColor' is not defined)
    s: Constant size (used if 'nameSize' is not defined)
    sizmin: Size corresponding to the smallest value (used if 'nameSize' is defined)
    sizmax: Size corresponding to the largest value (used if 'nameSize' is defined)
    flagAbsSize: Represent the Absolute value in Size representation
    flagCst: When True, the size is kept constant (equal to 's')
    flagLegendColor: Flag for representing the Color Legend
    flagLegendSize: Flag for representing the Size Legend
    legendNameColor: Title for the Color Legend (set to 'nameColor' if not defined)
    legendNameSize: Title for the Size Legend (set to 'size_name' if not defined)
    posX: rank of the first coordinate
    posY: rank of the second coordinate
    **kwargs : arguments passed to matplotllib.pyplot.scatter
    '''
    name = ''
    if (nameColor is None) and (nameSize is None):
        nameSize = __defaultVariable(db, None)
    
    # Read the coordinates
    tabx, taby = __readCoorPoint(db, nameCoorX, nameCoorY, useSel, posX, posY)
    nb = len(tabx)
    valid = np.full(nb, True)
    
    # Color of symbol
    colval = c
    if nameColor is not None:
        colval = __getDefinedValues(db, nameColor, 0, 1, None, useSel, 
                                    compress=False, asGrid=False, 
                                    flagConvertNanToZero=False)
        name = name + ' ' + nameColor
        valid = np.logical_and(valid, ~np.isnan(colval))
    else:
        colval = np.full(nb, c)

    # Size of symbol
    sizval = s
    if nameSize is not None:
        sizval = __getDefinedValues(db, nameSize, 0, 1, None, useSel, 
                                    compress=False, asGrid=False, 
                                    flagConvertNanToZero=False)
        m = np.nanmin(sizval)
        M = np.nanmax(sizval)

        if not flagCst:
            if flagAbsSize:
                sizval = np.absolute(sizval)
            m = np.nanmin(sizval)
            M = np.nanmax(sizval)
            if M > m:
                sizval = (sizmax - sizmin) * (sizval - m) / (M - m) + sizmin
                
            name = name + ' ' + nameSize
        else:
            sizval = s * ~np.isnan(sizval)
        valid = np.logical_and(valid, ~np.isnan(sizval))
    else:
        sizval = np.full(nb, s)

    if len(ax.get_title()) <= 0:
        ax.decoration(title = name)
    
    res = ax.scatter(x = tabx[valid], y = taby[valid], 
                     s = sizval[valid], c = colval[valid], **kwargs)

    if flagLegendColor:
        if nameColor is not None:
            if legendNameColor is None:
                legendNameColor = nameColor
            __addColorbar(res, ax, legendNameColor)
    
    if flagLegendSize and not flagCst:
        if nameSize is not None:
            if legendNameSize is None:
                legendNameSize = nameSize
            labels = lambda marker_size : (M - m)*(marker_size - sizmin)/(sizmax - sizmin) + m
            ax.legend(*res.legend_elements("sizes", num=6, func=labels), 
                      title=legendNameSize)
         
    return res

def literal(db, *args, **kwargs):
    '''
    Construct a layer for plotting a point data base, with literal variables
    
    ax: matplotlib.Axes (necessary when used as a method of the class)
    db: Db containing the variable to be plotted
    name: Name of the variable containing the label per sample
    nameCoorX: Name of the variable standing for X coordinate 
    nameCoorY: Name of the variable standing for Y coordinate 
    useSel : Boolean to indicate if the selection has to be considered
    flagLegend: Flag for representing the Color Bar
    legendName: title of the Legend (set to 'name' if not defined
    posX: rank of the first coordinate
    posY: rank of the second coordinate
    **kwargs : arguments passed to matplotllib.pyplot.scatter
    '''
    ax = __getNewAxes(None, 1)
    return __ax_literal(ax, db, *args, **kwargs)
    
def __ax_literal(ax, db, name=None, nameCoorX=None, nameCoorY=None, 
                 useSel=True, flagLegend=True, legendName=None, 
                 posX=0, posY=1, **kwargs):
    name = __defaultVariable(db, name)
    
    if len(ax.get_title()) <= 0:
        ax.decoration(title = db.getName(name)[0])
    
    # Read the coordinates
    tabx, taby = __readCoorPoint(db, nameCoorX, nameCoorY, useSel, posX, posY)
    
    labval = __getDefinedValues(db, name, 0, 1, None, useSel, 
                                compress=False, asGrid=False, 
                                flagConvertNanToZero=False)
    valid = ~np.isnan(labval)

    # We pass 'labval' to scatter function in order to mask undefined values
    res = ax.scatter(x = tabx[valid], y = taby[valid], **kwargs)
    
    for i, txt in enumerate(labval):
        if not np.isnan(txt):
            ax.annotate(round(txt,2), (tabx[i], taby[i]))
  
    if legendName is None:
        legendName = name
        
    return res

def gradient(db, *args, **kwargs):
    '''
    Construct a layer for plotting the gradient information of a data base
    
    ax: matplotlib.Axes (necessary when used as a method of the class)
    db: Db containing the variable to be plotted
    nameCoorX: Name of the variable standing for X coordinate 
    nameCoorY: Name of the variable standing for Y coordinate 
    useSel : Boolean to indicate if the selection has to be considered
    '''
    ax = __getNewAxes(None, 1)
    return __ax_gradient(ax, db, *args, **kwargs)
    
def __ax_gradient(ax, db, nameCoorX=None, nameCoorY=None, useSel=True, 
                  posX=0, posY=1, **kwargs):
    if db.getLocNumber(gl.ELoc.G) <= 0:
        return None
    
    # Extracting coordinates
    tabx, taby = __readCoorPoint(db, nameCoorX, nameCoorY, useSel, posX, posY)
    
    # Reading the Gradient components
    if db.getNDim() > 1:
        tabgx = db.getGradient(0,useSel)
        tabgy = db.getGradient(1,useSel)
    else:
        tabgy = -db.getGradient(0,useSel)
        tabgx = -np.ones(len(tabgy))

    if len(tabx) <= 0 or len(taby) <= 0 or len(tabgx) <= 0 or len(tabgy) <= 0:
        return None
    
    res = ax.quiver(tabx, taby, -tabgx, -tabgy, angles='xy', **kwargs)
            
    return res

def tangent(db, *args, **kwargs):
    '''
    Construct a layer for plotting a tangent data base
    
    ax: matplotlib.Axes (necessary when used as a method of the class)
    db: Db containing the variable to be plotted
    nameCoorX: Name of the variable standing for X coordinate 
    nameCoorY: Name of the variable standing for Y coordinate 
    useSel : Boolean to indicate if the selection has to be considered
    '''
    ax = __getNewAxes(None, 1)
    return __ax_tangent(ax, db, *args, **kwargs)

def __ax_tangent(ax, db, nameCoorX=None, nameCoorY=None, useSel=True, 
                 posX=0, posY=1, **kwargs):
    if db.getLocNumber(gl.ELoc.TGTE) <= 0:
        return None

    # Extracting coordinates
    tabx, taby = __readCoorPoint(db, nameCoorX, nameCoorY, useSel, posX, posY)

    # Extract Tangent information
    tabtx = db.getTangent(0,useSel)
    tabty = db.getTangent(1,useSel)

    if len(tabx) <= 0 or len(taby) <= 0 or len(tabtx) <= 0 or len(tabty) <= 0:
        return None
    
    res = ax.quiver(tabx, taby, -tabtx, -tabty, **kwargs)
    res = ax.quiver(tabx, taby,  tabtx,  tabty, **kwargs)
            
    return res

def point(db, *args, **kwargs):
    '''
    Construct a figure for plotting a point data base
    
    ax: matplotlib.Axes
    db: Db containing the variable to be plotted
    nameColor: Name of the variable containing the color per sample
    nameSize: Name of the variable containing the size per sample
    nameLabel: Name of the variable containing the label per sample
    nameCoorX: Name of the variable standing for X coordinate 
    nameCoorY: Name of the variable standing for Y coordinate 
    useSel : Boolean to indicate if the selection has to be considered
    color: Constant color (used if 'nameColor' is not defined)
    size: Constant size (used if 'nameSize' is not defined)
    sizmin: Size corresponding to the smallest value (used if 'nameSize' is defined)
    sizmax: Size corresponding to the largest value (used if 'nameSize' is defined)
    flagAbsSize: Represent the Absolute value in Size representation
    flagCst: When True, the size is kept constant (equel to 's')
    cmap: Optional Color scale
    flagGradient: Draw Gradient (if Gradients are defined)
    colorGradient: Color attached to the Gradient representation
    scaleGradient: Scale of the Gradient representation
    flagTangent: Draw Tangent (if Tangents are defined)
    colorTangent: Color attached to the Tangent representation
    scaleTangent: Scale of the Tangent representation
    flagLegendColor: Flag for representing the Color Legend (only if name-color is defined)
    flagLegendSize: Flag for representing the Size legend (only if nameSize is defined)
    flagLegendLabel: Flag for representing the Label Legend (only if nameLabel is defined)
    legendNameColor: Title for the Color Legend (set to 'nameColor' if not defined)
    legendNameSize: Title for the Size legend (set to 'nameSize' if not defined)
    legendNameLabel: Title for the Label Legend (set to 'nameLabel' if not defined)
    posX: rank of the first coordinate
    posY: rank of the second coordinate

    **kwargs : arguments passed to matplotllib.pyplot.scatter
    '''
    ax = __getNewAxes(None, 1)
    return __ax_point(ax, db, *args, **kwargs)
    
def __ax_point(ax, db, 
               nameColor=None, nameSize=None, nameLabel=None,
               nameCoorX=None, nameCoorY=None, useSel=True, 
               color='r', size=20, sizmin=10, sizmax=200, cmap=None,
               flagAbsSize=False, flagCst=False,
               flagGradient=False, colorGradient='black', scaleGradient=20,
               flagTangent=False, colorTangent='black', scaleTangent=20,
               flagLegendColor=False, flagLegendSize=False, flagLegendLabel=False, 
               legendNameColor=None, legendNameSize=None, legendNameLabel=None,
               posX=0, posY=1, **kwargs):

    if __isNotCorrect(object=db, types=["Db", "DbGrid"]):
        return None

    if (nameColor is None) and (nameSize is None) and (nameLabel is None):
        nameSize = __defaultVariable(db, None)
        flagCst = True

    title = ""
    if (nameColor is not None) or (nameSize is not None):
        pt = __ax_symbol(ax, db, nameColor=nameColor, nameSize=nameSize, 
                         nameCoorX=nameCoorX, nameCoorY=nameCoorY, useSel=useSel, 
                         c=color, s=size, sizmin=sizmin, sizmax=sizmax, 
                         flagAbsSize=flagAbsSize, flagCst=flagCst,cmap=cmap, 
                         flagLegendColor=flagLegendColor, flagLegendSize=flagLegendSize,
                         legendNameColor=legendNameColor, legendNameSize=legendNameSize,
                         posX=posX, posY=posY, 
                         **kwargs)
        if nameColor is not None:
            title = title + nameColor +  " (Color) "
        if nameSize is not None:
            if flagCst:
                title = title + " Sample Locations "
            else:
                title = title + nameSize + " (Size) "
    
    if nameLabel is not None:
        tx = __ax_literal(ax, db, name=nameLabel, 
                          nameCoorX=nameCoorX, nameCoorY=nameCoorY, useSel=useSel, 
                          flagLegend=flagLegendLabel, legendName=legendNameLabel,
                          posX=posX, posY=posY, **kwargs)
        title = title + nameLabel + " (Label) "
        
    if flagGradient:
        gr = __ax_gradient(ax, db, nameCoorX=nameCoorX, nameCoorY=nameCoorY, 
                           useSel=useSel, color=colorGradient, scale=scaleGradient,
                           posX=posX, posY=posY)
        title = title + " (Gradient) "

    if flagTangent:
        tg = __ax_tangent(ax, db, nameCoorX=nameCoorX, nameCoorY=nameCoorY,
                          useSel=useSel, color=colorTangent, scale=scaleTangent,
                          posX=posX, posY=posY)
        title = title + " (Tangent) "
    
    ax.decoration(title = title)
    
    return ax

def modelOnGrid(model, db, *args, **kwargs):
    ax = __getNewAxes(None, 1)
    return __ax_modelOnGrid(ax, model, db=db, *args, **kwargs)
    
def __ax_modelOnGrid(ax, model, db, useSel=True, icov=0, color='black', flagOrtho=True, **kwargs):
    '''
    Display the Model characteristics on a Grid
    This makes sense when the model contains some non-stationarity
    '''
    # Extracting coordinates
    tabx = db.getCoordinates(0,useSel)
    taby = db.getCoordinates(1,useSel)
    if len(tabx) <= 0 or len(taby) <= 0:
        return None
    
    gl.db_model_nostat(db, model, icov)
    tabR1 = db.getColumn("Nostat.Range-1", useSel)
    tabR2 = db.getColumn("Nostat.Range-2", useSel)
    tabA  = db.getColumn("Nostat.Angle-1", useSel)
    if len(tabR1) <= 0 or len(tabR2) <= 0 or len(tabA) <= 0:
        return None
    
    if flagOrtho:
        tabA = 90 + tabA
    ax.quiver(tabx, taby, tabR2, tabR2, angles=tabA, color=color, **kwargs)
            
    return ax
    
def polygon(poly, *args, **kwargs):
    '''
    Construct a Figure for plotting a polygon
    ax: matplotlib.Axes
    **kwargs: arguments passed to matplotlib.fill
    '''
    ax = __getNewAxes(None, 1)
    return __ax_polygon(ax, poly, *args, **kwargs)

def __ax_polygon(ax, poly, facecolor='yellow', edgecolor = 'blue', 
                 colorPerSet = False, flagEdge=True, flagFace=False, 
                 **kwargs):
    if __isNotCorrect(object=poly, types=["Polygons"]):
        return None
    
    npol = poly.getPolyElemNumber()
    cols = getColorMap(npol)
    
    for ipol in range(npol):
        x = poly.getX(ipol)
        y = poly.getY(ipol)
        
        facecolor_local = 'none'
        if flagFace:
            facecolor_local = facecolor
            if colorPerSet:
                facecolor_local = cols(ipol)
                
        edgecolor_local = 'none'
        if flagEdge:
            edgecolor_local = edgecolor
            if colorPerSet:
                edgecolor_local = cols(ipol)

        ax.fill(x, y, facecolor=facecolor_local, edgecolor=edgecolor_local,
                **kwargs)
        
    return ax

def __readGrid(dbgrid, name, useSel=True, posX=0, posY=1, corner=None, shading = "flat"):
    
    x0 = dbgrid.getX0(posX)
    y0 = dbgrid.getX0(posY)
    nx = dbgrid.getNX(posX)
    ny = dbgrid.getNX(posY)
    dx = dbgrid.getDX(posX)
    dy = dbgrid.getDX(posY)
    angle = 0
    if posX==0 and posY==1:
        angle = dbgrid.getAngle(posX)
    
    data = __getDefinedValues(dbgrid, name, posX, posY, corner, useSel, 
                              compress=False, asGrid=True)
    data = np.reshape(data, (ny,nx))

    tr = transform.Affine2D().rotate_deg_around(x0,y0,angle)
    
    if shading == "nearest":
        X = np.linspace(x0, x0 + (nx-1)*dx, nx)
        Y = np.linspace(y0, y0 + (ny-1)*dy, ny)
    elif shading == "flat":
        X = np.linspace(x0, x0 + nx*dx, nx+1)
        Y = np.linspace(y0, y0 + ny*dy, ny+1)
    else:
        print("The argument shading should be either 'nearest' for cells centered on (x,y)"
              " or 'flat' for cells with low-left corner in (x,y)")
        
    return x0, y0, X, Y, data, tr

def cell(dbgrid, *args, **kwargs):
    '''
    Plotting the cell edges from a DbGrid 

    ax: matplotlib.Axes (necessary when used as a method of the class)
    dbgrid: DbGrid containing the variable to be plotted
    **kwargs : arguments passed to matplotlib.pyplot.pcolormesh
    '''
    ax = __getNewAxes(None, 1)
    return __ax_cell(ax, dbgrid, *args, **kwargs)

def __ax_cell(ax, dbgrid, posX=0, posY=1, step=1, **kwargs):

    indices = np.zeros(dbgrid.getNDim())
    shift = np.ones(dbgrid.getNDim()) * (-1)
    for i in range(0,dbgrid.getNX(posX)+1,step):
        indices[posX] = i
        indices[posY] = 0
        tab1 = dbgrid.getCoordinatesByIndice(indices, True, shift)
        indices[posY] = dbgrid.getNX(posY)
        tab2 = dbgrid.getCoordinatesByIndice(indices, True, shift)
        ax.plot([tab1[0],tab2[0]],[tab1[1],tab2[1]], **kwargs)
    for i in range(0,dbgrid.getNX(posY)+1,step):
        indices[posX] = 0
        indices[posY] = i
        tab1 = dbgrid.getCoordinatesByIndice(indices, True, shift)
        indices[posX] = dbgrid.getNX(posX)
        tab2 = dbgrid.getCoordinatesByIndice(indices, True, shift)
        ax.plot([tab1[0],tab2[0]],[tab1[1],tab2[1]], **kwargs)
    return

def __ax_box(ax, dbgrid, posX=0, posY=1, **kwargs):
    indices = np.zeros(dbgrid.getNDim())
    shift = np.ones(dbgrid.getNDim()) * (-1)
    step = dbgrid.getNX(posX)
    for i in range(0,dbgrid.getNX(posX)+1,step):
        indices[posX] = i
        indices[posY] = 0
        tab1 = dbgrid.getCoordinatesByIndice(indices, True, shift)
        indices[posY] = dbgrid.getNX(posY)
        tab2 = dbgrid.getCoordinatesByIndice(indices, True, shift)
        ax.plot([tab1[0],tab2[0]],[tab1[1],tab2[1]], **kwargs)
    
    step = dbgrid.getNX(posY)
    for i in range(0,dbgrid.getNX(posY)+1,step):
        indices[posX] = 0
        indices[posY] = i
        tab1 = dbgrid.getCoordinatesByIndice(indices, True, shift)
        indices[posX] = dbgrid.getNX(posX)
        tab2 = dbgrid.getCoordinatesByIndice(indices, True, shift)
        ax.plot([tab1[0],tab2[0]],[tab1[1],tab2[1]], **kwargs)
    return

def raster(dbgrid, *args, **kwargs):
    '''
    Plotting a variable from a DbGrid in Raster

    ax: matplotlib.Axes (necessary when used as a method of the class)
    dbgrid: DbGrid containing the variable to be plotted
    name: Name of the variable to be represented (by default, the first Z locator, or the last field)
    useSel : Boolean to indicate if the selection has to be considered
    flagLegend: Flag for representing the Color Bar
    legendName: Name given to the Legend (set to 'name' if not defined)
    **kwargs : arguments passed to matplotlib.pyplot.pcolormesh
    '''
    ax = __getNewAxes(None, 1)
    return __ax_raster(ax, dbgrid, *args, **kwargs)

def __ax_raster(ax, dbgrid, name=None, useSel = True, posX=0, posY=1, corner=None, 
                flagLegend=False, legendName=None, **kwargs):
    name = __defaultVariable(dbgrid, name)
            
    if len(ax.get_title()) <= 0:
        ax.decoration(title = dbgrid.getName(name)[0])
    
    x0, y0, X, Y, data, tr = __readGrid(dbgrid, name, useSel, posX=posX, posY=posY, corner=corner)
    
    res = ax.pcolormesh(X, Y, data, transform=tr + ax.transData, **kwargs)
    
    if flagLegend:
        if legendName is None:
            legendName = name
        __addColorbar(res, ax, legendName)
    
    return res
        
def isoline(dbgrid, *args, **kwargs):
    '''
    Plotting a variable (referred by its name) with isoline representation from a DbGrid

    ax: matplotlib.Axes (necessary when used as a method of the class)
    dbgrid: DbGrid containing the variable to be plotted
    name: Name of the variable to be represented (by default, the first Z locator, or the last field)
    useSel : Boolean to indicate if the selection has to be considered
    levels: Vector of isovalues to be represented
    flagLegend: Flag for representing the Color Bar (not represented if alpha=0)
    legendName: Name given to the Legend (set to 'name' if not defined)
    ax: Reference for the plot within the figure
    
    **kwargs : arguments passed to matplotlib.pyplot.contour
    '''
    ax = __getNewAxes(None, 1)
    return __ax_isoline(ax, dbgrid, *args, **kwargs)

def __ax_isoline(ax, dbgrid, name=None, useSel = True, 
                 posX=0, posY=1, corner=None, levels=None,
                 flagLegend=False, legendName=None, **kwargs):
    name = __defaultVariable(dbgrid, name)
        
    if len(ax.get_title()) <= 0:
        ax.decoration(title = dbgrid.getName(name)[0])
    
    x0, y0, X, Y, data, tr = __readGrid(dbgrid, name, useSel, posX=posX, posY=posY, 
                                        corner=corner, shading="nearest")
    trans_data = tr + ax.transData
    
    res = ax.contour(X, Y, data, levels, **kwargs)
    
    if flagLegend:
        h1,l1 = res.legend_elements()
        if legendName is None:
            legendName = name
        ax.legend([h1[0]], [legendName])
        
    return res

def grid(dbgrid, *args, **kwargs):
    '''
    Plotting a variable (referred by its name) informed in a DbGrid

    dbgrid: DbGrid containing the variable to be plotted
    nameRaster: Name of the variable to be represented as raster
    nameContour: Name of the variable tp be represented as contours
    useSel : Boolean to indicate if the selection has to be considered
    flagCell: When True, the edge of the grid cells are represented
    flagBox: when True, the bounding box of the Grid is represented
    flagLegendRaster: Flag for representing the Raster Legend
    flagLegendContour: Flag for representing the Contour Legend
    legendNameRaster: Title for the Raster Legend (set to 'nameRaster' if not defined)
    legendNameSize: Title for the Contour Legend (set to 'nameContour' if not defined)
    **kwargs : arguments passed to matplotlib.pyplot.pcolormesh
    '''
    ax = __getNewAxes(None, 1)
    return __ax_grid(ax, dbgrid, *args, **kwargs)

def __ax_grid(ax, dbgrid, nameRaster = None, nameContour = None, useSel = True, 
              posX=0, posY=1, corner=None, flagCell=False, flagBox=False,
              flagLegendRaster=False, flagLegendContour=False,
              legendNameRaster=None, legendNameContour=None,
              levels=None, **kwargs):
    if __isNotCorrect(object=dbgrid, types=["DbGrid"]):
        return None

    if (nameRaster is None) and (nameContour is None) and (not flagCell):
        nameRaster = __defaultVariable(dbgrid, None)

    title = ""
    if nameRaster is not None:
        rs = __ax_raster(ax, dbgrid = dbgrid, name = nameRaster, useSel = useSel,  
                         posX=posX, posY=posY, corner=corner, 
                         flagLegend=flagLegendRaster, legendName=legendNameRaster,
                         **kwargs)
        title = title + nameRaster + " (Raster) "
    
    if nameContour is not None:
        ct = __ax_isoline(ax, dbgrid = dbgrid, name = nameContour, useSel = useSel, 
                          posX=posX, posY=posY, corner=corner, levels=levels, 
                          flagLegend=flagLegendContour, legendName=legendNameContour,
                          **kwargs)
        title = title + nameContour + " (Isoline) "
    
    if flagCell:
        cl = __ax_cell(ax, dbgrid, posX=posX, posY=posY, **kwargs)
    
    if flagBox:
        cl = __ax_box(ax, dbgrid, posX=posX, posY=posY, **kwargs)
        

    ax.decoration(title = title)
    
    return ax

def grid1D(dbgrid, *args, **kwargs):
    '''
    Plotting a variable (referred by its name) informed in a DbGrid

    dbgrid: DbGrid containing the variable to be plotted
    name: Name of the variable to be represented (by default, the first Z locator, or the last field)
    useSel : Boolean to indicate if the selection has to be considered
    flagLegend: Flag for representing the Legend
    **kwargs : arguments passed to matplotlib.pyplot.curve
    '''
    ax = __getNewAxes(None, 1)
    return __ax_grid1D(ax, dbgrid, *args, **kwargs)

def __ax_grid1D(ax, dbgrid, name = None, useSel = True,
                color='black',flagLegend=False, label='curve',
                **kwargs):
    if dbgrid.getNDim() != 1:
        print("This function is dedicated to 1-D Grid")
        return None
    
    if not(dbgrid.isGrid()):
        print("This function is dedicated to Grid Db and cannot be used here")
        return None
    
    if name is None:
        if dbgrid.getLocNumber(gl.ELoc.Z) > 0:
            name = dbgrid.getNameByLocator(gl.ELoc.Z,0) # select locator z1, prints an error if no Z locator
        else : # if no Z locator, choose the last field
            name = dbgrid.getLastName()
    x0 = dbgrid.getX0(0)
    nx = dbgrid.getNX(0)
    dx = dbgrid.getDX(0)
    
    tabx = dbgrid.getColumnByLocator(gl.ELoc.X, 0, useSel)
    data = __getDefinedValues(dbgrid, name, 0, 1, None, useSel, 
                              compress=False, asGrid=True)

    __ax_curve(ax, data1=tabx, data2=data, color=color, flagLegend=flagLegend, 
               **kwargs)

    ax.decoration(title = dbgrid.getName(name)[0])
        
    return ax

def histogram(db, *args, **kwargs):
    '''
    Plotting the histogram of a variable contained in a Db
    ax: matplotlib.Axes
    kwargs : arguments passed to matplotlib.pyplot.hist
    '''
    ax = __getNewAxes(None, 0)
    return __ax_histogram(ax, db=db, *args, **kwargs)
    
def __ax_histogram(ax, db, name, useSel=True, **kwargs):
    
    if __isNotCorrect(object=db, types=["Db", "DbGrid"]):
        return None

    db.useSel = useSel
    val = db[name]
    if len(val) == 0:
        return None
    
    ax.hist(val, **kwargs)
    
    ax.decoration(title = db.getName(name)[0], xlabel="Values", ylabel="Count")
        
    return ax

def sortedcurve(tabx, taby, *args, **kwargs):
    '''
    Plotting a set of points after they have been sorted in increasing X
    '''
    ax = __getNewAxes(None, 0)
    return __ax_sortedcurve(ax, tabx=tabx, taby=taby, *args, **kwargs)

def __ax_sortedcurve(ax, tabx, taby, color='black', flagLegend=False,
                     *args, **kwargs):
    # Account for possible 'nan'  values
    mask = np.logical_and(np.isfinite(tabx), np.isfinite(taby))
    stabx = tabx[mask]
    staby = taby[mask]
    
    # Indices of the sorted elements of stabx
    indices = np.argsort(stabx)
    return __ax_curve(ax, data1=stabx[indices], data2=staby[indices], color=color, 
                      flagLegend=flagLegend, *args, **kwargs)
    
def curve(data1, data2=None, *args, **kwargs):
    '''
    Plotting the curve of an array (argument 'data1')
        if data1 is a tuple, it should contain x=data1[0] and y=data1[1]
        or
        'data1' and 'data2' are provided
        otherwise:
        icas=1 when 'data1' contains the abscissa and ordinates are regular
        icas=2 when 'data1' contains the ordinate and abscissa are regular
    **kwargs : arguments passed to matplotlib.pyplot.plot
    '''
    ax = __getNewAxes(None, 0)
    return __ax_curve(ax, data1=data1, data2=data2, *args, **kwargs)

def __ax_curve(ax, data1, data2=None, icas=1, color0='black',flagLegend=False, 
               **kwargs):
    color = kwargs.setdefault('color', color0)
    label = kwargs.setdefault('label', 'curve')
    
    if len(data1) == 0:
        return None

    filetype = type(data1).__name__
    if filetype == "tuple":
        tabx = data1[0]
        taby = data1[1]
    elif filetype == "ndarray" and len(data1) == 2:
        tabx = data1[0]
        taby = data1[1]
    else:
        nbpoint = len(data1)
        if data2 is not None:
            if len(data2) != nbpoint:
                print("Arrays 'data1' and 'data2' should have same dimensions")
                return None
            tabx = data1
            taby = data2
        else:
            regular = [i for i in range(nbpoint)]
            if icas == 1:
                tabx = data1
                taby = regular
            else:
                tabx = regular
                taby = data1
    
    ax.plot(tabx, taby, **kwargs)
    
    if flagLegend:
        ax.legend()
        
    return ax

def multisegments(center, data, *args, **kwargs):
    ax = __getNewAxes(None, 1)
    return __ax_multisegments(ax, center=center, data=data, **kwargs)

def __ax_multisegments(ax, center, data, color='black',flagLegend=False, 
                       label="segments", *args, **kwargs):
    '''
    Plotting a set of segments joining 'center' to any of vertices
    stored in 'data'.
    **kwargs : arguments passed to matplotlib.pyplot.plot
    '''
    color = kwargs.setdefault('color', color)
    label = kwargs.setdefault('label', label)
        
    if len(data) == 0:
        return None
    
    nseg = len(data[0])
    
    for iseg in range(nseg):
        ax.plot([center[0],data[0][iseg]], [center[1],data[1][iseg]], **kwargs)
    
    if flagLegend:
        ax.legend()
        
    return ax

def fault(faults, *args, **kwargs):
    '''
    Plotting a Fault system.
    **kwargs : arguments passed to matplotlib.pyplot.plot
    '''
    ax = __getNewAxes(None, 1)
    return __ax_fault(ax, faults=faults, *args, **kwargs)

def __ax_fault(ax, faults, color='black', flagLegend=False, label="segments",
               **kwargs):
    color = kwargs.setdefault('color', color)
    label = kwargs.setdefault('label', label)
        
    nfaults = faults.getNFaults()
    for ifault in range(nfaults):

        fault = faults.getFault(ifault);
        npoints = fault.getNPoints() - 1
        xtab = fault.getX()
        ytab = fault.getY()
        ax.plot(xtab, ytab, **kwargs)
    
    if flagLegend:
        ax.legend()
        
    return ax

def XY(xtab, ytab, *args, **kwargs):
    ax = __getNewAxes(None, 0)
    return __ax_XY(ax, xtab=xtab, ytab=ytab, *args, **kwargs)

def __ax_XY(ax, xtab, ytab, flagAsPoint=False, flagLegend=False, 
            color='blue', marker='o', markersize=5, linestyle='-', label='data', 
            **kwargs):

    kwargs.setdefault('label', label)
    kwargs.setdefault('color', color)
    if flagAsPoint:
        kwargs.setdefault('markersize', markersize)
        kwargs.setdefault('marker', marker)
        kwargs.setdefault('linestyle',linestyle)
    else:
        kwargs.setdefault('linestyle',linestyle)

    if not len(ytab) == len(xtab):
        print("Arrays 'xtab' and 'ytab' should have same dimensions")
        return None;
    
    ax.plot(xtab, ytab, **kwargs)
            
    if flagLegend:
        ax.legend()
        
    return ax

def sample(sampleobj, *args, **kwargs):
    ax = __getNewAxes(None, 1)
    return __ax_sample(ax, sampleobj=sampleobj, *args, **kwargs)

def __ax_sample(ax, sampleobj, color='black', marker='o', markersize=10,
                linestyle=' ', flagLegend=False, label='data', 
                **kwargs):
    
    ax.plot(sampleobj[0], sampleobj[1], marker=marker, markersize=markersize, color=color,
            linestyle=linestyle, label=label, **kwargs)
            
    if flagLegend:
        ax.legend()
        
    return ax
    
def rule(ruleobj, *args, **kwargs):
    ax = __getNewAxes(None, 0)
    return __ax_rule(ax, ruleobj=ruleobj, *args, **kwargs)

def __ax_rule(ax, ruleobj, proportions=[],cmap=None, maxG=3.):
    if __isNotCorrect(object=ruleobj, types=["Rule"]):
        return None
    
    nfac = ruleobj.getFaciesNumber()
    ruleobj.setProportions(proportions)
    
    cols = getColorMap(nfac, cmap)

    for ifac in range(nfac):
        bds = ruleobj.getThresh(ifac+1)
        rect = ptc.Rectangle((bds[0],bds[2]),bds[1]-bds[0], bds[3]-bds[2], 
                              color=cols(ifac))
        ax.add_patch(rect)

    ax.geometry(xlim=[-maxG,+maxG], ylim=[-maxG,+maxG]) 
       
    return ax

def table(tableobj, ranks=None, *args, **kwargs):
    '''
    Plotting the contents of a Table (argument 'table')
    ax: matplotlib.Axes
    icols: designates the ranks of the variable (0: ordinate; 1: abscissae [or regular]) 
    fmt: designates [marker][line][color] information
    **kwargs
    '''
    ax = __getNewAxes(None, 0)
    return __ax_table(ax, tableobj=tableobj, icols=ranks, *args, **kwargs)
    
def __ax_table(ax, tableobj, icols, fmt='ok', flagLegend=False, **kwargs):
    if __isNotCorrect(object=tableobj, types=["Table"]):
        return None
    
    if len(icols) == 0:
        datay = tableobj.getColumn(0)
        datax = [i for i in range(tableobj.getNRows())]
    elif len(icols) == 1:
        datay = tableobj.getColumn(int(icols[0]))
        datax = [i for i in range(tableobj.getNRows())]
    else:
        datay = tableobj.getColumn(int(icols[0]))
        datax = tableobj.getColumn(int(icols[1]))
    
    data = np.stack((np.array(datax), np.array(datay)))
    data = data[:, ~np.isnan(data).any(axis=0)]

    ax.plot(data[0,:], data[1,:], **kwargs)
    
    if flagLegend:
        ax.legend()
        
    return ax

def mesh(meshobj, *args, **kwargs):
    """
    Plotting the contents of a Mesh
    **kwargs : arguments passed to matplotlib.pyplot.fill
    """
    ax = __getNewAxes(None, 1) 
    return __ax_mesh(ax, meshobj=meshobj, *args, **kwargs)

def __ax_mesh(ax, meshobj, 
              flagEdge=True, flagFace=False, flagApex=False, 
              facecolor="yellow", edgecolor="blue", linewidth=1,
              **kwargs):
    if __isNotCorrect(object=meshobj, types=["Mesh","MeshETurbo","MeshEStandardExt"]):
        return None
    
    if flagFace:
        kwargs.setdefault('facecolor', facecolor)
    else:
        kwargs.setdefault('facecolor', "white")
       
    if flagEdge:
        kwargs.setdefault('edgecolor', edgecolor) 
        kwargs.setdefault('linewidth', linewidth)

    nmesh = meshobj.getNMeshes()
    
    for imesh in range(nmesh):
        tabx = meshobj.getCoordinatesPerMesh(imesh, 0, True)
        taby = meshobj.getCoordinatesPerMesh(imesh, 1, True)
        ax.fill(tabx, taby, **kwargs)
        
        if flagApex:
            ax.scatter(tabx, taby, color='black')

    return ax

def correlation(db, namex, namey, *args, **kwargs):
    '''
    Plotting the scatter plot between two variables contained in a Db
    
    kwargs: additional arguments used in hist2d or scatter
    '''
    ax = __getNewAxes(None, 0)
    return __ax_correlation(ax, db=db, namex=namex, namey=namey, *args, **kwargs)
    
def __ax_correlation(ax, db, namex, namey, db2=None, 
                     asPoint = False,  flagSameAxes=False,
                     diagLine=False, diagColor="black", diagLineStyle='-',
                     bissLine=False, bissColor="red", bissLineStyle='-',
                     regrLine=False, regrColor="blue", regrLineStyle='-',
                     **kwargs):
    if __isNotCorrect(object=db, types=["Db", "DbGrid"]):
        return None
        
    if db2 is None:
        db2 = db
   
    if db.getSampleNumber() != db2.getSampleNumber():
        print("Db and Db2 should have the same number of samples")
        return None

    res = gl.correlationPairs(db, db, namex, namey)
    tabx = db.getValuesByNames(res[0], [namex])
    if len(tabx) == 0:
        return None
    taby = db2.getValuesByNames(res[0], [namey])
    if len(taby) == 0:
        return None
    
    xmin = np.min(tabx)
    xmax = np.max(tabx)
    ymin = np.min(taby)
    ymax = np.max(taby)
    
    range = None
    if flagSameAxes:
        xmin = ymin = min(xmin, ymin)
        xmax = ymax = max(xmax, ymax)
        ax.geometry(xlim=[xmin, xmax], ylim=[ymin, ymax])

    if asPoint:
        ax.scatter(tabx, taby, **kwargs)
    else:
        range = [[xmin,xmax],[ymin,ymax]]
        ax.hist2d(tabx, taby, range=range, **kwargs)

    if diagLine:
        u=[xmin, xmax]
        v=[ymin, ymax]
        ax.plot(u,v,color=diagColor,linestyle=diagLineStyle)
        
    if bissLine:
        bmin = min(xmin, ymin)
        bmax = max(xmax, ymax)
        u=[bmin, bmax]
        ax.plot(u,u,color=bissColor,linestyle=bissLineStyle)

    if regrLine:
        regr = gl.regression(db2, namey, [namex], flagCst=True)
        if regr.getNvar() == 0:
            return None
        a = regr.getCoeff(0)
        b = regr.getCoeff(1)
        u=[xmin, xmax]
        v=[a+b*xmin, a+b*xmax]
        ax.plot(u,v,color=regrColor,linestyle=regrLineStyle)
        
    ax.decoration(xlabel = db.getName(namex)[0], ylabel = db.getName(namey)[0])

    return ax

def hscatter(db, namex, namey, varioparam, ipas, *args, **kwargs):
    '''
    Plotting the scatter plot between two variables contained in a Db
    
    kwargs: additional arguments used in hist2d or scatter
    '''
    ax = __getNewAxes(None, 0)
    return __ax_hscatter(ax, db=db, namex=namex, namey=namey, varioparam=varioparam, ipas=ipas, 
                         *args, **kwargs)
    
def __ax_hscatter(ax, db, namex, namey, varioparam, ipas=0, idir=0, 
                     asPoint = False,  flagSameAxes=False,
                     diagLine=False, diagColor="black", diagLineStyle='-',
                     bissLine=False, bissColor="red", bissLineStyle='-',
                     **kwargs):
    if __isNotCorrect(object=db, types=["Db", "DbGrid"]):
        return None
        
    res = gl.hscatterPairs(db, namex, namey, varioparam, ipas, idir)
    tabx = db.getValuesByNames(res[0], [namex])
    if len(tabx) == 0:
        return None
    taby = db.getValuesByNames(res[0], [namey])
    if len(taby) == 0:
        return None
    
    xmin = np.min(tabx)
    xmax = np.max(tabx)
    ymin = np.min(taby)
    ymax = np.max(taby)
    
    range = None
    if flagSameAxes:
        xmin = ymin = min(xmin, ymin)
        xmax = ymax = max(xmax, ymax)
        ax.geometry(xlim=[xmin, xmax], ylim=[ymin, ymax])

    if asPoint:
        ax.scatter(tabx, taby, **kwargs)
    else:
        range = [[xmin,xmax],[ymin,ymax]]
        ax.hist2d(tabx, taby, range=range, **kwargs)

    if diagLine:
        u=[xmin, xmax]
        v=[ymin, ymax]
        ax.plot(u,v,color=diagColor,linestyle=diagLineStyle)
        
    if bissLine:
        bmin = min(xmin, ymin)
        bmax = max(xmax, ymax)
        u=[bmin, bmax]
        ax.plot(u,u,color=bissColor,linestyle=bissLineStyle)

    ax.decoration(xlabel = db.getName(namex)[0], ylabel = db.getName(namey)[0])

    return ax

def anam(anamobj, *args, **kwargs):
    ax = __getNewAxes(None, 0)
    return __ax_anam(ax, anamobj=anamobj, *args, **kwargs)
    
def __ax_anam(ax, anamobj, color='blue', linestyle='-', flagLegend=False):
    
    if __isNotCorrect(object=anamobj, types=["Anam","AnamHermite"]):
        return None

    res = anamobj.sample()
    
    ax = __ax_XY(ax, res.getY(), res.getZ(),
                 flagLegend=flagLegend, color=color, linestyle=linestyle,
                 label='Anamorphosis')
    ax.geometry(xlim = res.getAylim(), ylim=res.getAzlim())
    ax.decoration(xlabel="Gaussian values", ylabel="Raw values")
    
    return ax

def neigh(ax, neigh, grid, node=0, flagCell=False, flagZoom=False, **kwargs):
    
    # Identify target location
    target = grid.getSampleCoordinates(node)
    
    # Represent the target location
    __ax_sample(ax, target, **kwargs)
    
    # Represent the edge of the target (if block)
    if flagCell and grid.isGrid():
        __ax_curve(ax, grid.getCellEdges(node), **kwargs)
    
    # Represent the Neighborhood Ellipsoid
    if neigh.getType() == gl.ENeigh.MOVING:
        __ax_curve(ax, neigh.getEllipsoid(target), **kwargs)
    
        # Represent the Angular sectors
        if neigh.getFlagSector():
            __ax_multisegments(ax, target, neigh.getSectors(target), **kwargs)
        
        # Zoom to the Maximum radius circle (optional)
        if flagZoom:
            limits = neigh.getZoomLimits(target)
            ax.set_xlim(limits[0])
            ax.set_ylim(limits[1])
    
def neighWeights(ax, res, flagWeights=True, 
                 horizontalalignment='center',
                 verticalalignment='bottom',
                 **kwargs):
    # Number of neighboring samples
    nech = res.nech
    
    # Get the coordinates of the neighborhoods
    X = res.xyz[0]
    Y = res.xyz[1]
    ax.XY(X, Y, flagAsPoint=True, linestyle='')
    
    # Annotate the weights
    if flagWeights:
        for i in range(nech):
            ax.annotate(round(100.*res.wgt.getValue(i,0),2), (X[i], Y[i]), 
                        horizontalalignment=horizontalalignment,
                        verticalalignment=verticalalignment,
                        **kwargs)
    
def drawCircles(m,M,middle = False):
    x = np.linspace(-m,m,100)
    plt.plot(x, np.sqrt(m**2-x**2),c="g")
    plt.plot(x,-np.sqrt(m**2-x**2),c="g")
    
    x = np.linspace(-M,M,100)
    plt.plot(x, np.sqrt(M**2-x**2),c="g")
    plt.plot(x,-np.sqrt(M**2-x**2),c="g")
    if middle:
        mid = .5 * (m+M)
        x = np.linspace(-mid,mid,100)
        plt.plot(x, np.sqrt(mid**2-x**2),c="r")
        plt.plot(x,-np.sqrt(mid**2-x**2),c="r")
        
def drawDir(angle,col,R=10000):
    a = np.deg2rad(angle)
    u = R*np.array([0,np.cos(a)])
    v = R*np.array([0,np.sin(a)])
    plt.plot( u, v,c=col)
    plt.plot(-u,-v,c=col)
    
def drawCylrad(angle,cylrad,R=10000,col="purple"):
    a = np.deg2rad(angle)
    x = cylrad/np.sin(a-np.pi/2)
    
    u=R*np.array([0,np.cos(a)])
    v=R*np.array([0,np.sin(a)])
    
    plt.plot( u, (v+x),c=col)
    plt.plot(-u,-(v+x),c=col)
    plt.plot( u, (v-x),c=col)
    plt.plot(-u,-(v-x),c=col)
    
# Function to retrieve the limits of a lag
def lagDefine(i, lag, tol=0):
    center = mini = maxi = i*lag
    if tol>0:
        reltol = tol * lag / 100
        if i>0:
            mini = center - reltol
        maxi = maxi + reltol
    return mini, center, maxi  

def plot(object, name1=None, name2=None, ranks=None, **kwargs):
    
    filetype = type(object).__name__

    if filetype == "Db":
        if name2 is None:
            point(object, **kwargs)
        else:
            correlation(object, namex=name1, namey=name2, **kwargs)
            
    elif filetype == "DbGrid":
        grid(object, name1, **kwargs)
    
    elif filetype == "Vario":
        variogram(object, **kwargs)
    
    elif filetype == "Model":
        model(object, **kwargs)
    
    elif filetype == "Mesh":
        mesh(object, **kwargs)
    
    elif filetype == "Rule":
        rule(object, **kwargs)
    
    elif filetype == "Table":
        table(object, ranks, **kwargs)

    elif filetype == "Polygons":
        poly(object, **kwargs)
        
    else:
        print("Unknown type:",filetype)

def plotFromNF(filename, name1=None, name2=None, ranks=None, **kwargs):
    filetype = gl.ASerializable.getFileIdentity(filename)
    if filetype == "":
        exit()

    if filetype == "Db":
        db = gl.Db.createFromNF(filename,False)
        plot(db, name1, name2, **kwargs)
            
    elif filetype == "DbGrid":
        dbgrid = gl.DbGrid.createFromNF(filename,False)
        plot(dbgrid, name1, **kwargs)
    
    elif filetype == "Vario":
        vario_item = gl.Vario.createFromNF(filename,False)
        plot(vario_item, **kwargs)
    
    elif filetype == "Model":
        model_item = gl.Model.createFromNF(filename,False)
        plot(model_item, **kwargs)
    
    elif filetype == "Rule":
        rule_item = gl.Rule.createFromNF(filename,False)
        plot(rule_item, **kwargs)
    
    elif filetype == "Table":
        table_item = gl.Table.createFromNF(filename,False)
        plot(table_item,ranks, **kwargs)

    elif filetype == "Polygons":
        polygon_item = gl.Polygons.createFromNF(filename,False)
        plot(polygon_item, **kwargs)
        
    else:
        print("Unknown type:",filetype)

# Select data on interactive figures with matplotlib

from matplotlib.widgets import PolygonSelector
from matplotlib.path import Path

class PointSelection:
    """
    Select indices from a matplotlib collection using point selector.
    Left click on data points for selecting it, and right click to remove selection on the point.
    Press 'escape' to remove the current selection and start a new one.

    Selected indices are saved in the `ind` attribute, and the mask of the selection in the
    'mask' attribute. If mydb is provided, a new variable "interactive_selection" is created. 
    This tool changes color (to red by default) for the selected points.

    Note that this tool selects collection objects based on their *origins*
    (i.e., `offsets`).

    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        Axes to interact with.
    collection : `matplotlib.collections.Collection` subclass
        Collection you want to select from.
        At least one of 'ax' or 'collection' must be provided.
    mydb : gstlearn.Db or DbGrid. If provided, a new variable "interactive_selection" is 
        created (or modified if already existing)
    pickradius : precision of the picker in points (default is 7)
    color : new color for the selected points (default is red)
    """
    def __init__(self, ax=None, collection=None, mydb=None, pickradius=7, 
                 color='r', verbose=False):
        self.ax = ax
        if ax is None and collection is None:
            raise ValueError("ax and collection cannot be None at the same time,"
                             " at least one must be given.")
        elif collection is None:
            self.collection = ax.collections[0]
        else:
            self.collection = collection
            if ax is None:
                self.ax = self.collection.axes
            
        self.fig = self.collection.get_figure()
        self.collection.set_picker(True)
        self.collection.set_pickradius(pickradius)
        self.color = mcolors.to_rgba(color)
        self.verbose = verbose
            
        self.cid = self.fig.canvas.mpl_connect('pick_event', self.on_pick)
        self.cid_esc = self.fig.canvas.mpl_connect("key_press_event", self.onkeypress)
        
        self.data = self.collection.get_offsets()
        self.Ndata = len(self.data)
        
        self.collection.update_scalarmappable()
        colors = self.collection.get_facecolor()
        if len(colors)==1:
            colors = np.empty((self.Ndata,4))
            colors[:,:] = self.collection.get_facecolors()[0]
        self.initial_colors = np.copy(colors)
        
        self.list_clicks = []
        self.ind = []
        self.mask = np.zeros(self.Ndata)
        self.mydb = mydb
        if mydb is not None:
            mydb["interactive_selection"] = self.mask
        
        print("Select points on the plot: left click for selecting, right click to remove selection, "
              "'escape' for deleting current selection and starting a new one")
        
    def on_pick(self,event):
        self.event=event
        if event.mouseevent.button in (1,3):
            xmouse, ymouse = event.mouseevent.xdata, event.mouseevent.ydata
            ind = event.ind
            if self.verbose:
                print( 'Artist picked:', event.artist)
                print( '{} vertices picked'.format(len(ind)))
                print( 'Vertices picked:',ind)
                print( 'x, y of mouse: {:.2f},{:.2f}'.format(xmouse, ymouse))
                print( 'Data point:', self.data[ind[0]])
            for i in range(len(ind)):
                self.list_clicks.append(ind[i])
                if event.mouseevent.button == 1: #left click = select
                    self.mask[ind[i]] = 1
                else: #right click = unselect
                    self.mask[ind[i]] = 0
                self.ind = np.where(self.mask == 1)[0]
                
            if self.mydb is not None:
                self.mydb["interactive_selection"] = self.mask
            
            self.collection.update_scalarmappable()
            self.collection.set_array(None)
            colors = np.copy(self.initial_colors)
            colors[self.ind] = np.array(self.color)
            self.collection.set_facecolors(colors)
            self.fig.canvas.draw()

    def onkeypress(self,event):
        """Reinitialize when press ESC"""
        if event.key == "escape":
            self.list_clicks = []
            self.ind = []
            self.mask = np.zeros(self.Ndata)
            if self.mydb is not None:
                # self.mydb.deleteColumn("interactive_selection")
                self.mydb["interactive_selection"] = self.mask
            self.collection.set_facecolors(self.initial_colors)
            self.fig.canvas.draw()
             
    def disconnect(self):
        """Disconnect matplotlib events (ends interactive selection), 
        and deletes variable "interactive_selection" in 'mydb' if provided."""
        self.fig.canvas.mpl_disconnect(self.cid)
        self.fig.canvas.mpl_disconnect(self.cid_esc)
        self.mydb.deleteColumn("interactive_selection")


class PolygonSelection:
    """
    Select indices from a matplotlib collection using `PolygonSelector`.
    Draw polygon to select point inside a region. Press 'escape' to remove the polygon and start a new one.

    Selected indices are saved in the `ind` attribute, and the mask of the selection in the
    'mask' attribute. If mydb is provided, a new variable "interactive_selection" is created. 
    This tool fades out the points that are not part of the selection (i.e., reduces their alpha
    values). If your collection has alpha < 1, this tool will permanently alter the alpha values.

    Note that this tool selects collection objects based on their *origins*
    (i.e., `offsets`).

    Parameters
    ----------
    ax : `~matplotlib.axes.Axes`
        Axes to interact with.
    collection : `matplotlib.collections.Collection` subclass
        Collection you want to select from.
        At least one of 'ax' or 'collection' must be provided.
    mydb : gstlearn.Db or DbGrid. If provided, a new variable "interactive_selection" is 
        created (or modified if already existing)
    alpha_other : 0 <= float <= 1
        To highlight a selection, this tool sets all selected points to an
        alpha value of 1 and non-selected points to *alpha_other*.
    """

    def __init__(self, ax=None, collection=None, mydb=None, alpha_other=0.3):
        self.ax = ax
        if ax is None and collection is None:
            raise ValueError("ax and collection cannot be None at the same time,"
                             " at least one must be given.")
        elif collection is None:
            self.collection = ax.collections[0]
        else:
            self.collection = collection
            if ax is None:
                self.ax = self.collection.axes
        self.canvas = self.ax.figure.canvas
        self.alpha_other = alpha_other

        self.xys = self.collection.get_offsets()
        self.Ndata = len(self.xys)

        # Ensure that we have separate colors for each object
        self.collection.update_scalarmappable()
        self.fc = self.collection.get_facecolors()
        if len(self.fc) == 0:
            raise ValueError('Collection must have a facecolor')
        elif len(self.fc) == 1:
            self.fc = np.tile(self.fc, (self.Ndata, 1))
        self.initial_xlim = self.ax.get_xlim()
        
        self.poly = PolygonSelector(self.ax, self.onselect)
        self.ind = []
        self.mask = np.zeros(self.Ndata)
        self.mydb = mydb
        if self.mydb is not None:
            self.mydb["interactive_selection"] = self.mask
        
        print("Draw a polygon on the plot to select points inside the polygon."
              "Press 'escape' for deleting current polygon and starting a new one.")

    def onselect(self, verts):
        path = Path(verts)
        self.ind = np.nonzero(path.contains_points(self.xys))[0]
        self.mask = np.zeros(self.Ndata)
        self.mask[self.ind] = 1
        if self.mydb is not None:
            self.mydb["interactive_selection"] = self.mask
        
        self.collection.update_scalarmappable()
        self.collection.set_array(None)
        self.fc[:, -1] = self.alpha_other
        self.fc[self.ind, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()

    def disconnect(self):
        """Disconnect matplotlib events (ends interactive selection), 
        and deletes variable "interactive_selection" in 'mydb' if provided."""
        self.poly.disconnect_events()
        self.fc[:, -1] = 1
        self.collection.set_facecolors(self.fc)
        self.canvas.draw_idle()
        self.mydb.deleteColumn("interactive_selection")

## Add plot functions as methods of the class ##
## ------------------------------------------ ##
import gstlearn.plot         as gp

# Functions called using the generic *plot* function, based on the object recognition
setattr(gl.Db,               "plot",             gp.point)
setattr(gl.DbGrid,           "plot",             gp.grid)
setattr(gl.Polygons,         "plot",             gp.polygon)
setattr(gl.Rule,             "plot",             gp.rule)
setattr(gl.Faults,           "plot",             gp.fault)
setattr(gl.AnamHermite,      "plot",             gp.anam)
setattr(gl.Vario,            "plot",             gp.variogram)
setattr(gl.Model,            "plot",             gp.model)

setattr(gl.Table,            "plot",             gp.table)
setattr(gl.MeshETurbo,       "plot",             gp.mesh)
setattr(gl.Db,               "histogram",        gp.histogram)
setattr(gl.Db,               "correlation",      gp.correlation)
setattr(gl.Db,               "hscatter",         gp.hscatter)
setattr(gl.Db,               "grid1D"   ,        gp.grid1D)
setattr(gl.Vario,            "varmod",           gp.varmod)

# New style attribute setting functions
setattr(plt.Axes, "decoration",    gp.decoration)
setattr(plt.Axes, "geometry",      gp.geometry)

# Functions considered as members f the Axis class
# The name "grid" must not be used as confusing for matplotlib
setattr(plt.Axes, "gstgrid",       gp.__ax_grid)
setattr(plt.Axes, "gstpoint",      gp.__ax_point)

# Functions considered as members f the Axis class
setattr(plt.Axes, "polygon",       gp.__ax_polygon)
setattr(plt.Axes, "rule",          gp.__ax_rule)
setattr(plt.Axes, "fault",         gp.__ax_fault)
setattr(plt.Axes, "anam",          gp.__ax_anam)
setattr(plt.Axes, "grid1D"   ,     gp.__ax_grid1D)
setattr(plt.Axes, "curve",         gp.__ax_curve)
setattr(plt.Axes, "sortedcurve",   gp.__ax_sortedcurve)
setattr(plt.Axes, "multisegments", gp.__ax_multisegments)
setattr(plt.Axes, "histogram",     gp.__ax_histogram)
setattr(plt.Axes, "correlation",   gp.__ax_correlation)
setattr(plt.Axes, "hscatter",      gp.__ax_hscatter)
setattr(plt.Axes, "table",         gp.__ax_table)

setattr(plt.Axes, "model",         gp.__ax_model)
setattr(plt.Axes, "mesh",          gp.__ax_mesh)
setattr(plt.Axes, "variogram",     gp.__ax_variogram)

setattr(plt.Axes, "neigh",         gp.neigh)
setattr(plt.Axes, "neighWeights",  gp.neighWeights)

setattr(gl.Db,    "point",         gp.point)

setattr(gl.Db,    "symbol",        gp.symbol)
setattr(gl.DbGrid,"symbol",        gp.symbol)
setattr(plt.Axes, "symbol",        gp.__ax_symbol)

setattr(gl.Db,    "literal",       gp.literal)
setattr(gl.DbGrid,"literal",       gp.literal)
setattr(plt.Axes, "literal",       gp.__ax_literal)

setattr(gl.Db,    "gradient",      gp.gradient)
setattr(gl.DbGrid,"gradient",      gp.gradient)
setattr(plt.Axes, "gradient",      gp.__ax_gradient)

setattr(gl.Db,    "tangent",       gp.tangent)
setattr(gl.DbGrid,"tangent",       gp.tangent)
setattr(plt.Axes, "tangent",       gp.__ax_tangent)

setattr(gl.DbGrid,"raster",        gp.raster)
setattr(plt.Axes, "raster",        gp.__ax_raster)

setattr(gl.DbGrid,"isoline",       gp.isoline)
setattr(plt.Axes, "isoline",       gp.__ax_isoline)

setattr(gl.DbGrid,"cell",          gp.cell)
setattr(plt.Axes, "cell",          gp.__ax_cell)

setattr(plt.Axes, "XY",            gp.__ax_XY)

setattr(plt.Figure, "varmod",      gp.__ax_varmod)
setattr(plt.Figure, "variogram",   gp.__ax_variogram)
setattr(plt.Figure, "model",       gp.__ax_model)

setattr(plt.Figure, "decoration",  gp.decoration)
setattr(plt.Figure, "geometry",    gp.geometry)
