################################################################################
#                                                                              #
#                         gstlearn Python package                              #
# Copyright (c) (2023) MINES Paris / ARMINES                                   #
# Authors: gstlearn Team                                                       #
# License: BSD 3-clause                                                        #
#                                                                              #
################################################################################
try:
    import matplotlib
    import matplotlib.pyplot     as plt
    import matplotlib.patches    as ptc
    import matplotlib.transforms as transform
    import matplotlib.colors     as mcolors
except ModuleNotFoundError as ex:
    msg = ("Python dependency 'matplotlib' not found.\n"
          "To install it alongside gstlearn, please run `pip install gstlearn[plot]'")
    raise ModuleNotFoundError(msg) from ex

import numpy                 as np
import numpy.ma              as ma
import gstlearn              as gl
import gstlearn.plot         as gp
# import gstlearn.proj         as prj
import math

from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy                   import shape
from pandas.io               import orc
from matplotlib.pyplot       import axes

#Set of global values
defaultDims = [[5,5], [8,8]]
defaultXlim = [ None, None ]
defaultYlim = [ None, None ]
defaultAspect = [ 'auto', 1 ]

def getColorMap(n, cmap=None):
    '''
    Returns a resampled Matplotlib colormap for a given number of colors
    
    n: requested number of different colors
    cmap: name of a listed matplotlib colormap, or an instance of Colormap or None.
    '''
    if isinstance(cmap, matplotlib.colors.Colormap):
        return cmap.resampled(n)
    name = cmap
    if name is None:
        name = 'viridis'
    return plt.colormaps[name].resampled(n)
    
def _selectiItems(nvalues, sitem=-1):
    outs = range(0, nvalues)
    nout = nvalues
    if sitem >= 0:
        outs = range(sitem, sitem+1)
        nout = 1
    return outs, nout

def _getArgument(arg, default, *args, **kwargs):
    if arg in kwargs:
        return kwargs[arg]
    elif arg in args:
        return arg
    return default

def _isNotCorrect(object, types):
    if object is None:
        print("Argument 'object' must be provided")
        return True
    filetype = type(object).__name__
    if filetype in types:
        return False
    
    print("Argument 'object' (",filetype,") must be a valid type among",types)
    return True

def _getDefaultVariableName(db, name):
    if name is None:
        if db.getNLoc(gl.ELoc.Z) > 0:
            name = db.getNameByLocator(gl.ELoc.Z,0)
        else : # if no Z locator, choose the last field
            name = db.getLastName()
    else:
        if db.getUID(name) < 0:
            name = db.getLastName()
    return name

def geometry(dims=None, xlim=None, ylim=None, aspect=None):
    '''
    Define the geometry for the current graphic
    If Axes is defined, the arguments are passed to _ax_geometry
    Otherwize if Figure is defined, the arguments are passed to _ax_geometry for all axes.
    Otherwise, this statement is ignored.

    dims: Extension of ALL graphic Axes
    xlim: Range of values along the X-axis
    ylim: Range of values along the Y-axis
    aspect: Y/X ratio
    '''
    fig = plt.gcf()
    ax  = _getCurrentAx(fig)
    if ax is not None:
        _ax_geometry(ax, dims=dims, xlim=xlim, ylim=ylim, aspect=aspect)
        return
    
    if fig is not None:
        axes = fig.get_axes()
        for ax in axes:
            _ax_geometry(ax, dims=dims, xlim=xlim, ylim=ylim, aspect=aspect)
        return
    
def _ax_geometry(ax, dims=None, xlim=None, ylim=None, aspect=None):
    '''
    Set the default values for the geometric parameters for one Axes
    
    ax: matplotlib.Axes
    dims: Extension of ALL graphic Axes (even if attached to one Axes in particular)
    xlim: Range of values along the X-axis
    ylim: Range of values along the Y-axis
    aspect: Y/X ratio
    '''
    if dims is not None:
        ax.figure.set_size_inches(dims[0], dims[1])
    if xlim is not None:
        ax.set_xlim(left = xlim[0], right = xlim[1])
    if ylim is not None:
        ax.set_ylim(bottom = ylim[0], top = ylim[1])
    if aspect is not None:
        ax.set_aspect(aspect)
            
def decoration(title=None, xlabel=None, ylabel=None, **kwargs):
    '''
    Attach a decoration to the current graphic (either current Figure or current Axes)

    title: Title assigned to the Figure
    xlabel: Label on the X-axis (only for the current Axes)
    ylabel: Label on the Y-axis (only for the current Axes)
    '''
    fig = plt.gcf()
    if fig is None:
        return
    
    if _isMultiAxes(fig):
        _fig_decoration(fig, title=title, **kwargs)
        return
    else:
        ax  = _getCurrentAx(fig)
        if ax is not None:
            _ax_decoration(ax, title=title, xlabel=xlabel, ylabel=ylabel, **kwargs)
        return
    
def _fig_decoration(fig, title=None, **kwargs):
    '''
    Add the decoration to a Figure.
    fig: matplotlib.Figure
    title: title contents
    '''
    if title is not None:
        fig.suptitle(title, **kwargs)
    
def _ax_decoration(ax, title=None, xlabel=None, ylabel=None, **kwargs):
    '''
    Add the decoration to a Axes.
    ax: matplotlib.Axes
    title: title contents (for the main for a collection of Axes)
    xlabel: label along the horizontal axis
    ylabel: label along the vertical axis
    '''
    if title is not None:
        ax.set_title(title, **kwargs)
    if xlabel is not None:
        ax.set_xlabel(xlabel)
    if ylabel is not None:
        ax.set_ylabel(ylabel)

def _isMultiAxes(fig):
    axs = fig.get_axes()

    if len(axs) <= 1:
        return False
    
    # Find the number of Axes in the initial subplots definition of the Figure
    grid_spec = axs[0].get_subplotspec().get_gridspec()
    nrows, ncols = grid_spec.get_geometry()
    return nrows * ncols > 1

def _getAxesFromFigureSubplots(fig):
    axs = fig.get_axes()

    # Find the number of Axes in the initial subplots definition of the Figure
    grid_spec = axs[0].get_subplotspec().get_gridspec()
    nrows, ncols = grid_spec.get_geometry()

    # If nrows * ncols > 1, turn 'axs' back to the numpy array
    if nrows * ncols > 1:
        axs = np.array(axs[:nrows*ncols]).reshape(nrows, ncols)
    return axs

def _getCurrentAx(fig):
    '''
    Returns the current Axes.
    Note that if an Axes has been dedicated to Legend, the previous one is returned
    '''
    axs = fig.get_axes()
    ax  = plt.gca()
    index = axs.index(ax)

    # Find the initial number of subplots in the figure
    grid_spec = axs[0].get_subplotspec().get_gridspec()
    nrows, ncols = grid_spec.get_geometry()
    ninit = nrows * ncols

    if index >= ninit:
        return axs[ninit-1]
    else:
        return ax

def _getNewAxes(nx=1, ny=1):
    ''' 
    Creates a new figure (possibly containing multiple subplots)
    nx, ny: Number of subplots along X and Y

    Remarks
        If 'ax' does not exist, a new figure is created. 
        Otherwise, the input argument is simply returned.
    '''
    if len(plt.get_fignums()) <= 0:
        fig, ax = plt.subplots(nx, ny, squeeze=False)
    else:
        fig = plt.gcf()
        ax = _getAxesFromFigureSubplots(fig)
    
    if not _isMultiAxes(fig):
        ax = _getFirstElement(ax)
    return ax

def _getFirstElement(tab):
    '''
    Returns the first element of an Axes or an array of Axes.
    tab:  Argument to be checked
    '''
    if hasattr(tab, '__getitem__') and len(tab) > 0:  # Check that `tab` is subscriptable and not empty
        first = tab[0]
        if hasattr(first, '__getitem__') and len(first) > 0:  # Check that `tab[0]` is also subscriptable
            return first[0]
        return first
    return None  # Return None if `tab` is empty or not subscriptable.

def close():
    ''' 
    General procedure for closing a graphic with matplotlib.pyplot.
    It allows suppressing the dependency to matplotlib in your scripts
    '''
    plt.show()

def init(nx=1, ny=1, figsize=None, flagEqual=False):
    ''' 
    General procedure for initializing a new graphic with matplotlib.pyplot.
    It allows suppressing the dependency to matplotlib in your scripts
    nx, ny:     Number of subplots along X and Y
    figsize: When defined, this dictates the dimension of all Axes 
    flagEqual: When True, all subsequent Axes have an 'aspect' set to 1
    '''
    fig, axs = plt.subplots(nrows=nx, ncols=ny)

    aspect = None
    if flagEqual:
        aspect = 1

    if nx * ny > 1:
        for ax in axs.flat:
            _ax_geometry(ax, dims=figsize, aspect=aspect)
    else:
        _ax_geometry(axs, dims=figsize, aspect=aspect)

    return fig, axs

def _legendContinuous(ax, im, legendName = None):
    '''
    Attach the Legend of a continuous variable

    ax: matplotlib.Axes
    im: matplolib element giving the color map
    legendName: Name given to the Legend
    '''
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, ax=ax, cax=cax)
    if legendName is not None:
        cbar.ax.set_title(legendName)
    return cbar

def _legendDiscrete(ax, sizmin, sizmax, sizvmin, sizvmax, color, 
                    legendName = None, loc='upper right', num = 5, ndec=3):

    indices = np.arange(num)
    sizes = sizmin  + indices * (sizmax - sizmin)
    values = sizvmin + indices * (sizvmax - sizvmin)
    size_handles = [plt.scatter([], [], s=sizes[i], color=color, alpha=0.6, edgecolors='k', 
                                label=f"{round(values[i],ndec)}") for i in indices]
    ax.legend(handles=size_handles, title=legendName, loc=loc)

def _getCoordinates(db, nameCoorX=None, nameCoorY=None, useSel=True, posX=0, posY=1):
    '''
    Returns a set of two vectors containing the coordinates

    db: current Data Base
    nameCoorX: Name of a variable that may be used to given the coordinates along X
    nameCoorY: Name of a variable that may be used to given the coordinates along X
    useSel: True if the selection must be taken into account
    posX: rank of the first coordinate
    posY: rank of the second coordinate

    Remarks:
    - if 'nameCoorX' is defined, this variable provides the coordinates.
      Otherwise, we use the coordinates (Locator:X) for the 'posX' item.
      - resp for 'nameCoorY'.
    '''
    if nameCoorX is not None:
        tabx = db.getColumn(nameCoorX, useSel)
    else:
        if db.getNDim() > 0:
            tabx = db.getOneCoordinate(posX,useSel)
            
    if nameCoorY is not None:
        taby = db.getColumn(nameCoorY, useSel)
    else:
        if db.getNDim() > 1:
            taby = db.getOneCoordinate(posY,useSel)

    if len(tabx) <= 0 or len(taby) <= 0:
        return None
    
    return tabx, taby

def _getGradients(db, useSel=True):
    '''
    Returns a set of two vectors containing the two gradient components

    db: current Data Base
    useSel: True if the selection must be taken into account
    '''
    if db.getNLoc(gl.ELoc.G) <= 0:
        return None, None
    if db.getNDim() > 1:
        tabgx = db.getGradient(0,useSel)
        tabgy = db.getGradient(1,useSel)
    else:
        tabgy = -db.getGradient(0,useSel)
        tabgx = -np.ones(len(tabgy))

    return tabgx, tabgy

def _getTangents(db, useSel = True):
    '''
    Returns a set of two vectors containing the two tangent components

    db: current Data Base
    useSel: True if the selection must be taken into account
    '''
    if db.getNLoc(gl.ELoc.TGTE) <= 0:
        return None, None
        
    # Extract Tangent information
    tabtx = db.getTangent(0,useSel)
    tabty = db.getTangent(1,useSel)

    return tabtx, tabty

def _getVariable(db, name, posX=0, posY=1, corner=None, useSel=True, asGrid=True):

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

    return tab

def _getGridVariable(dbgrid, name, useSel=True, posX=0, posY=1, corner=None, shading = "flat"):
    
    x0 = dbgrid.getX0(posX)
    y0 = dbgrid.getX0(posY)
    nx = dbgrid.getNX(posX)
    ny = dbgrid.getNX(posY)
    dx = dbgrid.getDX(posX)
    dy = dbgrid.getDX(posY)
    angle = 0
    if posX==0 and posY==1:
        angle = dbgrid.getAngle(posX)
    
    data = _getVariable(dbgrid, name, posX, posY, corner, useSel, True)
    data = np.reshape(data, (ny,nx))

    tr = transform.Affine2D().rotate_deg_around(x0,y0,angle)
    Xrot = dbgrid.getOneSliceForCoordinate(posX, posX, posY, corner, useSel)
    Xrot = np.reshape(Xrot, (ny,nx))
    Yrot = dbgrid.getOneSliceForCoordinate(posY, posX, posY, corner, useSel)
    Yrot = np.reshape(Yrot, (ny,nx))

    if shading == "nearest":
        X = np.linspace(x0, x0 + (nx-1)*dx, nx)
        Y = np.linspace(y0, y0 + (ny-1)*dy, ny)
    elif shading == "flat":
        X = np.linspace(x0, x0 + nx*dx, nx+1)
        Y = np.linspace(y0, y0 + ny*dy, ny+1)
    else:
        print("The argument shading should be either 'nearest' for cells centered on (x,y)"
              " or 'flat' for cells with low-left corner in (x,y)")
    
    return x0, y0, X, Y, Xrot, Yrot, data, tr

def __ax_varioElem(ax, vario, ivar=0, jvar=0, idir=0, hmax=None, showPairs = False,
                   varColor='black', varLinestyle='dotted', 
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

def varmod(vario=None, model=None, ivar=-1, jvar=-1, *args, **kwargs):
    """
    Construct a figure for plotting experimental variogram(s) and model. Both of them are optional
    
    vario : experimental variogram to be represented
    model : optional, variogram model
    ivar, jvar : Indices of the variables for the variogram to be represented. If -1 (default), all 
                 variables are selected and all the simple and crossed variograms are represented.
    idir : Index of the direction of the variogram to be represented. If -1 (default) all available
           directions are selected and multi-directional variograms are represented.
    **kwargs : arguments passed to _ax_varmod.
    """
    nvar = 0
    if vario is not None:
        nvar = vario.getNVar()
    if model is not None:
        nvar = model.getNVar()
    ivarUtil, ivarN = _selectiItems(nvar, ivar)
    jvarUtil, jvarN = _selectiItems(nvar, jvar)
    axs = _getNewAxes(nx=ivarN, ny=jvarN)

    return _ax_varmod(axs, vario=vario, model=model, ivar=ivar, jvar=jvar, 
                       *args, **kwargs)

def _fig_varmod(fig, vario=None, model=None, **kwargs):
    axs = _getAxesFromFigureSubplots(fig)
    return _ax_varmod(axs, vario=vario, model=model, **kwargs)

def _ax_varmod(axs, vario=None, model=None, ivar=-1, jvar=-1, idir=-1,
               nh = 100, hmax = None, codir=None,
               showPairs=False, asCov=False, 
               flagDrawVariance = True,
               varioLinestyle = 'dashed', modelLinestyle = 'solid',
               varColor='black', varLinestyle="dotted",
               envColor='black', envLinestyle="dotted",
               cmap=None, flagLegend=False,
               **kwargs):
    """
    Construct a figure for plotting experimental variogram(s) and model. Both of them are optional
    
    axs : matplotlib.Axes
    vario : experimental variogram to be represented
    model : optional, variogram model
    ivar, jvar : Indices of the variables for the variogram to be represented. If -1 (default), all 
                 variables are selected and all the simple and crossed variograms are represented.
    idir : Index of the direction of the variogram to be represented. If -1 (default) all available
           directions are selected and multi-directional variograms are represented.
    nh : number of points between 0 and hmax where the model variogram is calculated (default is 100).
    hmax : Maximum distance to be represented.
    showPairs : True to show the number of pairs per lag
    flagDrawVariance : Flag to add the variance (default is True)  
    varioLinestyle: Linestyle for representing the experimental variogram
    modelLinestyle: Linestyle for representing the Model
    varColor, varLinestyle: parameters for representing variance-covariance line
    envColor, envLinestyle: parameters for representing coregionalization envelop
    cmap : Optional Color scale
    flagLegend : Flag to display the axes legend.
    **kwargs : arguments passed to matplotlib.pyplot.plot for all variograms plotted (not models!)
    """
    flagDef = False
    if vario is not None:
        if _isNotCorrect(object=vario, types=["Vario"]):
            return None
        flagDef = True
    if model is not None:
        if _isNotCorrect(object=model, types=["Model"]):
            return None
        flagDef = True
    if not flagDef:
        print("You must define either 'vario' or 'model' or both")
        return None

    color_in_kwargs = 'color' in kwargs

    ndir = 1
    if vario is not None:
        ndir = vario.getNDir()
    nvar = 0
    if vario is not None:
        nvar = vario.getNVar()
    if model is not None:
        nvar = model.getNVar() 
    cols = getColorMap(ndir,cmap)
    hmax = _getHmax(hmax, vario, model)
    
    ndirUtil, ivarD = _selectiItems(ndir, idir)
    ivarUtil, ivarN = _selectiItems(nvar, ivar)
    jvarUtil, jvarN = _selectiItems(nvar, jvar)
            
    if ndir > 1:
        flagLabelDir = True
    else:
        flagLabelDir = False
        
    # Loop on the variables 
    iaxes = 0 
    for iv in ivarUtil:
        for jv in jvarUtil:
            
            if type(axs) == np.ndarray:
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
                codirLoc = _getCodir(codir, idirUtil, vario, model)
                
                # Plotting the Variogram (optional)
                if vario is not None:
                    __ax_varioElem(ax, vario, iv, jv, idirUtil, 
                                   showPairs=showPairs, hmax=hmax,
                                   flagDrawVariance = flagDrawVariance,
                                   varColor=varColor, varLinestyle=varLinestyle,  
                                   flagLabelDir=flagLabelDir, flagLegend=flagLegend, 
                                   **kwargs)

                # Plotting the Model (optional)
                if model is not None:
                    if modelLinestyle is not None:
                        kwargs.update({'linestyle': modelLinestyle})
                    _ax_modelElem(ax, model, ivar=iv, jvar=jv, codir=codirLoc, 
                                  hmax=hmax, nh=nh, asCov=asCov,
                                  envColor=envColor, envLinestyle=envLinestyle, 
                                  flagLabelDir=flagLabelDir, flagLegend=flagLegend, 
                                  **kwargs)

            ax.autoscale(True)
            
            if vario is not None:
                if vario.drawOnlyPositiveX(iv, jv):
                    ax.set_xlim(left=0)
                if vario.drawOnlyPositiveY(iv, jv):
                    ax.set_ylim(bottom=0)

            # Add the point (0,0) to fix the scale of the graphic
            # This should be valid for any representation of Variogram and/or Model
            ax.plot(0., 0.)
    
    return axs

def variogram(vario, ivar=0, jvar=0, *args, **kwargs):
    """
    Plot experimental variogram(s) (can be multidirectional and multivariable or selected ones).
    
    vario : experimental variogram to be represented (gstlearn.Vario).
    ivar, jvar : Indices of the variables for the variogram to be represented. If -1 (default), all 
                 variables are selected and all the simple and crossed variograms are represented.
    **kwargs : arguments passed to _ax_variogram.
    """
    nvar = vario.getNVar()
    ivarUtil, ivarN = _selectiItems(nvar, ivar)
    jvarUtil, jvarN = _selectiItems(nvar, jvar)
    axs = _getNewAxes(nx=ivarN, ny=jvarN)
    
    return _ax_variogram(axs, vario, ivar=ivar, jvar=jvar, *args, **kwargs)
    
def _getCodir(codir, idir, vario, model):
    if codir is not None:
        return codir
    ndim = 0
    if vario is not None:
        ndim = vario.getNDim()
        return vario.getCodirs(idir)
    if model is not None:
        ndim = model.getNDim()
    if ndim <= 0:
        print("You must define either 'vario' or 'model' or both")
        return None
    else:
        codir = np.zeros(ndim)
        codir[0] = 1
        return codir
        
def _getHmax(hmax, vario, model):
    if hmax is not None:
        return hmax
    
    hmax = 0
    if model is not None:
        for icova in range(model.getNCov()):
            range_max = np.max(model.getCovAniso(icova).getRanges())
            if 3*range_max > hmax:
                hmax = 3*range_max
    if vario is not None:
        hmax = vario.getHmax()
    if hmax == 0: # if the model has no range defined
        hmax = 1
    return hmax

def _fig_variogram(fig, vario, **kwargs):
    axs = _getAxesFromFigureSubplots(fig)
    return _ax_variogram(axs, vario=vario, **kwargs)

def _ax_variogram(axs, vario, ivar=0, jvar=0, idir=0, *args, **kwargs):
    """
    Plot experimental variogram(s) (can be multidirectional and multivariable or selected ones).
    
    vario : experimental variogram to be represented (gstlearn.Vario).
    ivar, jvar : Indices of the variables for the variogram to be represented. If -1 (default), all 
                 variables are selected and all the simple and crossed variograms are represented.
    idir : Index of the direction of the variogram to be represented. If -1 (default) all available
           directions are selected and multidirectional variograms are represented.
    **kwargs : arguments passed to _ax_varmod
    """
    varioLineStyle = _getArgument("varioLineStyle", "solid", *args, **kwargs)
    return _ax_varmod(axs, vario, ivar=ivar, jvar=jvar, idir=idir, 
                      varioLinestyle = varioLineStyle, 
                      *args, **kwargs)

def _ax_modelElem(ax, modelobj, ivar=0, jvar=0, codir=None, vario=None, idir=0,
                  nh = 100, hmax = None, asCov=False,
                  envColor='black', envLinestyle='dashed',
                  label=None, flagLabelDir=False, flagEnvelop = True, flagLegend=False, 
                  **kwargs):
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
    if _isNotCorrect(object=modelobj, types=["Model"]):
        return None
             
    if label is None:
        if flagLabelDir:
            angles = gl.GeometryHelper.rotationGetAngles(codir,True)
            label = "model dir={}".format(np.round(angles,3))
        else:
            label = "model"

    istart = 0
    for i in range(modelobj.getNCov()):
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
        ggm = modelobj.envelop(hh, ivar, jvar, -1, codir, mode)
        ax.plot(hh[istart:], ggm[istart:], c = envColor, linestyle = envLinestyle)
    
    # Representation bounds
    if ivar == jvar:
        ax.set_xlim(left=0)
    ax.set_ylim(bottom=0)

    # Draw the Legend (optional)
    if flagLegend:
        ax.legend()
        
    return res

def model(modelobj, ivar=0, jvar=0, *args, **kwargs):
    """
    Construct a figure for plotting a model.
    
    modelobj : the variogram model
    ivar, jvar : Indices of the variables for the variogram to be represented. If -1 (default), all 
                 variables are selected and all the simple and crossed variograms are represented.
    **kwargs : arguments passed to _ax_model.
    """
    nvar = modelobj.getNVar()
    ivarUtil, ivarN = _selectiItems(nvar, ivar)
    jvarUtil, jvarN = _selectiItems(nvar, jvar)
    axs = _getNewAxes(nx=ivarN, ny=jvarN)
    return _ax_model(axs, modelobj, ivar=ivar, jvar=jvar, *args, **kwargs)
    
def _fig_model(fig, modelobj, **kwargs):
    axs = _getAxesFromFigureSubplots(fig)
    return _ax_model(axs, modelobj=modelobj, **kwargs)

def _ax_model(axs, modelobj, ivar=0, jvar=0, *args, **kwargs):
    """
    Construct a figure for plotting a Model.
    
    axs : current matplotlib.Axes
    modelobj : variogram model
    ivar, jvar : Indices of the variables for the variogram to be represented. If -1 (default), all 
                 variables are selected and all the simple and crossed variograms are represented.
    *args, **kwargs : arguments passed to _ax_varmod
    """
    return _ax_varmod(axs, model=modelobj, ivar=ivar, jvar=jvar, *args, **kwargs)

def symbol(db, nameColor=None, nameSize=None, *args, **kwargs):
    '''
    Construct a Layer for plotting a point data base, with optional color and size variables
    
    db: Db containing the variable to be plotted
    nameColor: Name of the variable containing the color per sample
    nameSize: Name of the variable containing the size per sample
    *args, **kwargs : arguments passed to _ax_symbol
    '''
    ax = _getNewAxes()
    return _ax_symbol(ax, db, nameColor=nameColor, nameSize=nameSize, 
                       *args, **kwargs)
    
def _ax_symbol(ax, db, nameColor=None, nameSize=None, 
               nameCoorX=None, nameCoorY=None, useSel=True, 
               c='r', s=20, sizmin=10, sizmax=200, flagAbsSize=False, flagCst=False,
               sizvmin=None, sizvmax = None, 
               flagLegendColor=False, flagLegendSize=False, 
               legendNameColor=None, legendNameSize=None, posX=0, posY=1, **kwargs):
    '''
    Construct a Layer for plotting a point data base, with optional color and size variables
    
    db: Db containing the variable to be plotted
    nameColor: Name of the variable containing the color per sample
    nameSize: Name of the variable containing the size per sample
    useSel : Boolean to indicate if the selection has to be considered
    c: Constant color (used if 'nameColor' is not defined)
    s: Constant size (used if 'nameSize' is not defined)
    sizmin: Size corresponding to the smallest value (used if 'nameSize' is defined)
    sizmax: Size corresponding to the largest value (used if 'nameSize' is defined)
    sizvmin: Minimum value below which samples are not represented with proportional symbols
    sizvmax: Maximum value above which samples are not represented with proportional symbols
    flagAbsSize: Represent the Absolute value in Size representation
    flagCst: when True, the symbol size is constant and equel to 's'
    flagLegendColor: Flag for representing the Color Legend
    flagLegendSize: Flag for representing the Size Legend
    legendNameColor: Title for the Color Legend (set to 'nameColor' if not defined)
    legendNameSize: Title for the Size Legend (set to 'size_name' if not defined)
    posX: rank of the first coordinate
    posY: rank of the second coordinate
    **kwargs : arguments passed to matplotllib.pyplot.scatter
    '''
    tabx, taby = _getCoordinates(db, nameCoorX, nameCoorY, useSel, posX, posY)
    if len(tabx) <= 0:
        return
    
    nb = len(tabx)
    valid = np.full(nb, True)
    name = ''
    
    # Color of symbol
    colval = c
    if nameColor is not None:
        name = name + ' ' + nameColor
        colval = _getVariable(db, nameColor, posX, posY, None, useSel, False)
        valid = np.logical_and(valid, ~np.isnan(colval))
    else:
        colval = np.full(nb, c)

    # Size of symbol
    sizval = s
    if nameSize is not None:
        name = name + ' ' + nameSize
        sizval = _getVariable(db, nameSize, posX, posY, None, useSel, False)
        valid = np.logical_and(valid, ~np.isnan(sizval))

        if flagCst:
            sizval = np.full(nb, s)
        else:
            if flagAbsSize:
                sizval = np.absolute(sizval)
            if sizvmin is None:
                sizvmin = np.nanmin(sizval)
            if sizvmax is None:
                sizvmax = np.nanmax(sizval)
            if sizvmax > sizvmin:
                sizval = np.where((sizval >= sizvmin) & (sizval <= sizvmax), 
                        np.interp(sizval, [sizvmin, sizvmax], [sizmin, sizmax]), s)
    else:
        sizval = np.full(nb, s)

    res = ax.scatter(x = tabx[valid], y = taby[valid], 
                     s = sizval[valid], c = colval[valid], **kwargs)

    if len(ax.get_title()) <= 0:
        ax.decoration(title = name)

    if flagLegendColor and nameColor is not None:
        _legendContinuous(ax, res, legendNameColor)
    
    if flagLegendSize and nameSize is not None:
        _legendDiscrete(ax, sizmin, sizmax, sizvmin, sizvmax, c, legendNameSize)
         
    return res

def literal(db, name=None, *args, **kwargs):
    '''
    Construct a layer for plotting a point data base, with literal variables
    
    db: Db containing the variable to be plotted
    name: Name of the variable containing the label per sample
    *args, **kwargs : arguments passed to _ax_literal
    '''
    ax = _getNewAxes()
    return _ax_literal(ax, db, name=name, *args, **kwargs)
    
def _ax_literal(ax, db, name=None, 
                nameCoorX=None, nameCoorY=None, useSel=True, 
                flagLegend=True, legendName=None, 
                posX=0, posY=1, fontsize=10, **kwargs):
    '''
    Construct a layer for plotting a point data base, with literal variables
    
    ax: matplotlib.Axes (necessary when used as a method of the class)
    db: Db containing the variable to be plotted
    name: Name of the variable containing the label per sample
    useSel : Boolean to indicate if the selection has to be considered
    flagLegend: Flag for representing the Color Bar
    legendName: title of the Legend (set to 'name' if not defined
    posX: rank of the first coordinate
    posY: rank of the second coordinate
    **kwargs : arguments passed to matplotllib.pyplot.scatter
    '''
    tabx, taby = _getCoordinates(db, nameCoorX, nameCoorY, useSel, posX, posY)
    if len(tabx) <= 0: 
        return
    
    labval = _getVariable(db, name, posX, posY, None, useSel, False)
    valid = ~np.isnan(labval)

    res = ax.scatter(x = tabx[valid], y = taby[valid], **kwargs)
    
    for i, txt in enumerate(labval):
        if not np.isnan(txt):
            ax.annotate(round(txt,2), (tabx[i], taby[i]), fontsize=fontsize)
  
    name = _getDefaultVariableName(db, name)
    
    if len(ax.get_title()) <= 0:
        ax.decoration(title = db.getName(name)[0])

    return res

def gradient(db, *args, **kwargs):
    '''
    Construct a layer for plotting the gradient information of a data base
    
    db: Db containing the variable to be plotted
    *args, **kwargs : arguments passed to _ax_gradient
    '''
    ax = _getNewAxes()
    return _ax_gradient(ax, db, *args, **kwargs)
    
def _ax_gradient(ax, db, 
                 nameCoorX=None, nameCoorY=None, useSel=True, 
                 posX=0, posY=1, **kwargs):
    '''
    Construct a layer for plotting the gradient information of a data base
    
    ax: matplotlib.Axes (necessary when used as a method of the class)
    db: Db containing the variable to be plotted
    useSel : Boolean to indicate if the selection has to be considered
    posX: rank of the first coordinate
    posY: rank of the second coordinate
    **kwargs : arguments passed to quiver.
    '''    
    tabx, taby = _getCoordinates(db, nameCoorX, nameCoorY, useSel, posX, posY)
    if len(tabx) <= 0: 
        return
    
    tabgx, tabgy = _getGradients(db, useSel)
    if len(tabgx) <= 0:
        return

    res = ax.quiver(tabx, taby, -tabgx, -tabgy, angles='xy', **kwargs)
            
    return res

def tangent(db, *args, **kwargs):
    '''
    Construct a layer for plotting a tangent data base
    
    db: Db containing the variable to be plotted
    *args, **kwargs : arguments passed to _ax_tangent
    '''
    ax = _getNewAxes()
    return _ax_tangent(ax, db, *args, **kwargs)

def _ax_tangent(ax, db, nameCoorX=None, nameCoorY=None, useSel=True, 
                posX=0, posY=1, **kwargs):
    '''
    Construct a layer for plotting a tangent data base
    
    ax: matplotlib.Axes (necessary when used as a method of the class)
    db: Db containing the variable to be plotted
    useSel : Boolean to indicate if the selection has to be considered
    posX: rank of the first coordinate
    posY: rank of the second coordinate
    **kwargs : arguments passed to quiver.
    '''
    tabx, taby = _getCoordinates(db, nameCoorX, nameCoorY, useSel, posX, posY)
    if len(tabx) <= 0:
        return

    tabtx, tabty = _getTangents(db, useSel)
    if len(tabtx) <= 0:
        return

    res = ax.quiver(tabx, taby, -tabtx, -tabty, **kwargs)
    res = ax.quiver(tabx, taby,  tabtx,  tabty, **kwargs)
            
    return res

def covaOnGrid(cova, db, *args, **kwargs):
    '''
    Display the Model characteristics on a Grid
    This makes sense when the model contains some non-stationarity
    '''
        
    ax = _getNewAxes()
    return _ax_covaOnGrid(ax, cova, db=db, *args, **kwargs)
    
def _ax_covaOnGrid(ax, cova, db, useSel=True, color='black', flagOrtho=True, **kwargs):

    tab = db.getAllCoordinates(useSel)
    if len(tab) <= 0 :
        return None
    
    tabR1 = cova.informCoords(tab,gl.EConsElem.RANGE,0)
    tabR2 = cova.informCoords(tab,gl.EConsElem.RANGE,1)
    tabA  = cova.informCoords(tab,gl.EConsElem.ANGLE,0)

    if len(tabR1) <= 0 or len(tabR2) <= 0 or len(tabA) <= 0:
        return None
    
    if flagOrtho:
        tabA = 90 + tabA
    ax.quiver(tab[0,:], tab[1,:], tabR2, tabR2, angles=tabA, color=color, **kwargs)
            
    return ax
    
def polygon(poly, *args, **kwargs):
    '''
    Construct a Figure for plotting a polygon
    poly: Polygon class
    *args, **kwargs: arguments passed to _ax_polygon
    '''
    ax = _getNewAxes()
    return _ax_polygon(ax, poly, *args, **kwargs)

def _ax_polygon(ax, poly, facecolor='yellow', edgecolor = 'blue', 
                 colorPerSet = False, flagEdge=True, flagFace=False, 
                 **kwargs):
    '''
    Construct a Figure for plotting a polygon
    ax: matplotlib.Axes
    poly: Polygon class
    facecolor: Color assigned to each polygon area
    edgecolor: Color assigned to the polygon edges
    colorPerSet: when True, each polygon is represented with a different color
    flagEdge: when True, the polygon edges are represented
    flagFace: when True, the polygon edges are represented
    **kwargs: arguments passed to matplotlib.fill
    '''
    if _isNotCorrect(object=poly, types=["Polygons"]):
        return None
    
    npol = poly.getNPolyElem()
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

def cell(dbgrid, *args, **kwargs):
    '''
    Plotting the cell edges from a DbGrid 

    dbgrid: DbGrid containing the variable to be plotted
    *args, **kwargs : arguments passed to _ax_cell
    '''
    ax = _getNewAxes()
    return _ax_cell(ax, dbgrid, *args, **kwargs)

def _ax_cell(ax, dbgrid, posX=0, posY=1, corner=None, step=1, **kwargs):
    '''
    Plotting the cell edges from a DbGrid 

    ax: matplotlib.Axes (necessary when used as a method of the class)
    dbgrid: DbGrid containing the variable to be plotted
    posX: rank of the first coordinate
    posY: rank of the second coordinate
    step: step for representing the cell edge every 'step' values
    **kwargs : arguments passed to subsequent functions
    '''
    shift = np.ones(dbgrid.getNDim()) * (-1)
    if corner is None:
        corner = np.zeros(dbgrid.getNDim())
    indices = corner
        
    for i in range(0,dbgrid.getNX(posX)+1,step):
        indices[posX] = i
        indices[posY] = 0
        tab1 = dbgrid.getCoordinatesByIndice(indices, True, shift)
        indices[posY] = dbgrid.getNX(posY)
        tab2 = dbgrid.getCoordinatesByIndice(indices, True, shift)
        ax.plot([tab1[posX],tab2[posX]],[tab1[posY],tab2[posY]], **kwargs)
    for i in range(0,dbgrid.getNX(posY)+1,step):
        indices[posX] = 0
        indices[posY] = i
        tab1 = dbgrid.getCoordinatesByIndice(indices, True, shift)
        indices[posX] = dbgrid.getNX(posX)
        tab2 = dbgrid.getCoordinatesByIndice(indices, True, shift)
        ax.plot([tab1[posX],tab2[posX]],[tab1[posY],tab2[posY]], **kwargs)
    return

def _ax_box(ax, dbgrid, posX=0, posY=1, step=1, **kwargs):
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

def raster(dbgrid, name=None, *args, **kwargs):
    '''
    Plotting a variable from a DbGrid in Raster

    dbgrid: DbGrid containing the variable to be plotted
    name: Name of the variable to be represented (by default, the first Z locator, or the last field)
    *args, **kwargs : arguments passed to _ax_raster
    '''
    ax = _getNewAxes()
    return _ax_raster(ax, dbgrid, name=name, *args, **kwargs)

def _ax_raster(ax, dbgrid, name=None, useSel = True, posX=0, posY=1, corner=None, 
               flagLegend=False, legendName=None, **kwargs):
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
    name = _getDefaultVariableName(dbgrid, name)

    x0, y0, X, Y, Xrot, Yrot, data, tr = _getGridVariable(dbgrid, name, useSel, posX=posX, posY=posY, corner=corner)
    
    res = ax.pcolormesh(X, Y, data, transform=tr + ax.transData, **kwargs)
    
    if flagLegend:
        _legendContinuous(ax, res, legendName)
            
    if len(ax.get_title()) <= 0:
        ax.decoration(title = dbgrid.getName(name)[0])

    return res
        
def isoline(dbgrid, name=None, *args, **kwargs):
    '''
    Plotting a variable (referred by its name) with isoline representation from a DbGrid

    dbgrid: DbGrid containing the variable to be plotted
    name: Name of the variable to be represented (by default, the first Z locator, or the last field)
    *args, **kwargs : arguments passed to _ax_isoline
    '''
    ax = _getNewAxes()
    return _ax_isoline(ax, dbgrid, name=name, *args, **kwargs)

def _ax_isoline(ax, dbgrid, name=None, useSel = True, 
                 posX=0, posY=1, corner=None, levels=None, nlevel=5,
                 flagLegend=False, legendName=None, **kwargs):
    '''
    Plotting a variable (referred by its name) with isoline representation from a DbGrid

    ax: matplotlib.Axes (necessary when used as a method of the class)
    dbgrid: DbGrid containing the variable to be plotted
    name: Name of the variable to be represented (by default, the first Z locator, or the last field)
    useSel : Boolean to indicate if the selection has to be considered
    levels: Vector of isovalues to be represented
    nlevel: Number of levels for automatic generation of 'labels' (if not provided)
    flagLegend: Flag for representing the Color Bar (not represented if alpha=0)
    legendName: Name given to the Legend (set to 'name' if not defined)
    ax: Reference for the plot within the figure
    
    **kwargs : arguments passed to matplotlib.pyplot.contour
    '''    
    x0, y0, X, Y, Xrot, Yrot, data, tr = _getGridVariable(dbgrid, name, useSel, posX=posX, posY=posY, corner=corner, shading="nearest")
    ax.set_xlim(np.nanmin(Xrot), np.nanmax(Xrot))
    ax.set_ylim(np.nanmin(Yrot), np.nanmax(Yrot))

    trans_data = tr + ax.transData
    
    if levels is None:
        levels = np.linspace(np.nanmin(data), np.nanmax(data), nlevel)

    res = ax.contour(Xrot, Yrot, data, levels, **kwargs)
    
    name = _getDefaultVariableName(dbgrid, name)
        
    if len(ax.get_title()) <= 0:
        ax.decoration(title = dbgrid.getName(name)[0])

    if flagLegend:
        h1,l1 = res.legend_elements()
        ax.legend([h1[0]], [legendName])
        
    return res

def line(dbline, *args, **kwargs):
    '''
    Plotting a variable (referred by its name) informed in a DbLine

    dbline: DbLine containing the variable to be plotted
    *args, **kwargs : arguments passed to _ax_line.
    '''
    ax = _getNewAxes()
    return _ax_line(ax, dbline, *args, **kwargs)

def _ax_line(ax, dbline, color = 'blue', colorPoint='black', colorHeader='red', 
              flagHeader=True, flagSample=False, flagAnnotateHeader=False, offset=[-1.0,0.5],
              **kwargs):
    '''
    Plotting a variable (referred by its name) informed in a DbLine

    dbline: DbLine containing the variable to be plotted
    name: Name of the variable to be represented
    useSel : Boolean to indicate if the selection has to be considered
    **kwargs : arguments passed to ...
    '''
    if _isNotCorrect(object=dbline, types=["DbLine"]):
        return None
    
    if dbline.getNDim() != 2:
        return None
    
    nbline = dbline.getNLine()
    
    for iline in range(nbline):
        x = dbline.getCoordinatesPerLine(iline, 0)
        y = dbline.getCoordinatesPerLine(iline, 1)
        
        ax.plot(x, y, color=color, **kwargs)

        if flagHeader:
            ax.plot(x[0], y[0], marker='D', color=colorHeader)
            if flagAnnotateHeader:
                ax.text(x[0]+offset[0], y[0]+offset[1], "L#"+str(iline+1))

        if flagSample:
            ax.plot(x, y, marker='.', color=colorPoint, linestyle='None')
        
    return ax

def graphO(dbgraphO, *args, **kwargs):
    '''
    Plotting a variable (referred by its name) informed in a DbGraphO

    dbgraphO: DbGraphO containing the variable to be plotted
    *args, **kwargs : arguments passed to _ax_graphO.
    '''
    ax = _getNewAxes()
    return _ax_graphO(ax, dbgraphO, *args, **kwargs)

def _ax_graphO(ax, dbgraphO, color = 'blue', colorPoint='black', flagSample=False, flagAnnotate=False,
                flagByRank=False, ndec=2, **kwargs):
    '''
    Plotting a variable (referred by its name) informed in a DbGraphO

    dbgraphO: DbGraphO containing the variable to be plotted
    useSel : Boolean to indicate if the selection has to be considered
    **kwargs : arguments passed to subsequent functions.
    '''
    if _isNotCorrect(object=dbgraphO, types=["DbGraphO"]):
        return None
    if dbgraphO.getNDim() != 2:
        return None
    
    if flagSample:
        x = dbgraphO.getOneCoordinate(0)
        y = dbgraphO.getOneCoordinate(1)
        ax.plot(x, y, marker='.', color=colorPoint, linestyle='None')
    
    narcs = dbgraphO.getNArc()
    for iarc in range(narcs):
        x = dbgraphO.getArc(iarc, 0)
        y = dbgraphO.getArc(iarc, 1)
        value = dbgraphO.getArcValue(iarc)
        if value > 0:
            xmid = (x[0]+x[1])/2.
            ymid = (y[0]+y[1])/2.
            rotation = np.arctan2(y[1]-y[0],x[1]-x[0])*180/np.pi
        
            ax.plot(x, y, color=color, **kwargs)

            if flagAnnotate:
                if flagByRank:
                    ax.text(xmid, ymid, str(iarc+1), ha="center", va="bottom", rotation=rotation)
                else:
                    ax.text(xmid, ymid, str(round(value,ndec)), ha="center", va="bottom", 
                            rotation=rotation)

    return ax

def grid1D(dbgrid, name=None, *args, **kwargs):
    '''
    Plotting a variable (referred by its name) informed in a DbGrid

    dbgrid: DbGrid containing the variable to be plotted
    name: Name of the variable to be represented (by default, the first Z locator, or the last field)
    *args, **kwargs : arguments passed to _ax_grid1D
    '''
    ax = _getNewAxes()
    return _ax_grid1D(ax, dbgrid, name=name, *args, **kwargs)

def _ax_grid1D(ax, dbgrid, name = None, useSel = True,
                color='black',flagLegend=False, label='curve',
                **kwargs):
    '''
    Plotting a variable (referred by its name) informed in a DbGrid

    dbgrid: DbGrid containing the variable to be plotted
    name: Name of the variable to be represented (by default, the first Z locator, or the last field)
    useSel : Boolean to indicate if the selection has to be considered
    flagLegend: Flag for representing the Legend
    **kwargs : arguments passed to matplotlib.pyplot.curve
    '''
    if dbgrid.getNDim() != 1:
        print("This function is dedicated to 1-D Grid")
        return None
    
    if not(dbgrid.isGrid()):
        print("This function is dedicated to Grid Db and cannot be used here")
        return None
    
    if name is None:
        if dbgrid.getNLoc(gl.ELoc.Z) > 0:
            name = dbgrid.getNameByLocator(gl.ELoc.Z,0) # select locator z1, prints an error if no Z locator
        else : # if no Z locator, choose the last field
            name = dbgrid.getLastName()
    x0 = dbgrid.getX0(0)
    nx = dbgrid.getNX(0)
    dx = dbgrid.getDX(0)
    
    tabx = dbgrid.getColumnByLocator(gl.ELoc.X, 0, useSel)
    data = _getVariable(dbgrid, name, 0, 1, None, useSel, True)

    _ax_curve(ax, data1=tabx, data2=data, color=color, flagLegend=flagLegend, 
               **kwargs)

    ax.decoration(title = dbgrid.getName(name)[0])
        
    return ax

def histogram(db, name=None, *args, **kwargs):
    '''
    Plotting the histogram of a variable contained in a Db

    db: Db containing the variable to be plotted
    name: Name of the variable to be used for histogram
    *args, **kwargs : arguments passed to _ax_histogram
    '''
    ax = _getNewAxes()
    return _ax_histogram(ax, db=db, name=name, *args, **kwargs)
    
def _ax_histogram(ax, db, name, useSel=True, **kwargs):
    '''
    Plotting the histogram of a variable contained in a Db
    ax: matplotlib.Axes
    db: Db containing the variable to be plotted
    name: Name of the variable to be used for histogram
    *args, **kwargs : arguments passed to matplotlib.pyplot.hist
    '''    
    if _isNotCorrect(object=db, types=["Db", "DbGrid"]):
        return None

    db.useSel = useSel
    val = db[name]
    if len(val) == 0:
        return None
    
    ax.hist(val, **kwargs)
    
    ax.decoration(title = db.getName(name)[0], xlabel="Values", ylabel="Count")
    return ax

def sortedCurve(tabx, taby, *args, **kwargs):
    '''
    Plotting a set of points after they have been sorted in increasing X
    '''
    ax = _getNewAxes()
    return _ax_sortedCurve(ax, tabx=tabx, taby=taby, *args, **kwargs)

def _ax_sortedCurve(ax, tabx, taby, color='black', flagLegend=False,
                     *args, **kwargs):
    # Account for possible 'nan'  values
    mask = np.logical_and(np.isfinite(tabx), np.isfinite(taby))
    stabx = tabx[mask]
    staby = taby[mask]
    
    # Indices of the sorted elements of stabx
    indices = np.argsort(stabx)
    return _ax_curve(ax, data1=stabx[indices], data2=staby[indices], color=color, 
                      flagLegend=flagLegend, *args, **kwargs)
    
def curve(data1, data2=None, *args, **kwargs):
    '''
    Plotting the curve of an array (argument 'data1')
    **kwargs : arguments passed to _ax_curve
    '''
    ax = _getNewAxes()
    return _ax_curve(ax, data1=data1, data2=data2, *args, **kwargs)

def _ax_curve(ax, data1, data2=None, icas=1, color0='black',flagLegend=False, 
               **kwargs):
    '''
    Plotting the curve of an array (argument 'data1')
        if data1 is a tuple, it should contain x=data1[0] and y=data1[1]
        or
        'data1' and 'data2' are provided
        otherwise:
        icas=1 when 'data1' contains the abscissa and ordinates are regular
        icas=2 when 'data1' contains the ordinate and abscissa are regular
    **kwargs : arguments passed to matplotlib.plot
    '''
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

def multiSegments(center, data, *args, **kwargs):
    ax = _getNewAxes()
    return _ax_multiSegments(ax, center=center, data=data, **kwargs)

def _ax_multiSegments(ax, center, data, color='black',flagLegend=False, 
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
    **kwargs : arguments passed to _ax_faults
    '''
    ax = _getNewAxes()
    return _ax_faults(ax, faults=faults, *args, **kwargs)

def _ax_faults(ax, faults, color='black', flagLegend=False, label="segments",
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
    ax = _getNewAxes()
    return _ax_XY(ax, xtab=xtab, ytab=ytab, *args, **kwargs)

def _ax_XY(ax, xtab, ytab, flagAsPoint=False, flagLegend=False, 
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
    ax = _getNewAxes()
    return _ax_sample(ax, sampleobj=sampleobj, *args, **kwargs)

def _ax_sample(ax, sampleobj, color='black', marker='o', markersize=10,
                linestyle=' ', flagLegend=False, label='data', 
                **kwargs):
    
    ax.plot(sampleobj[0], sampleobj[1], marker=marker, markersize=markersize, color=color,
            linestyle=linestyle, label=label, **kwargs)
            
    if flagLegend:
        ax.legend()
        
    return ax
    
def rule(ruleobj, *args, **kwargs):
    ax = _getNewAxes()
    return _ax_rule(ax, ruleobj=ruleobj, *args, **kwargs)

def _ax_rule(ax, ruleobj, proportions=[],cmap=None, maxG=3.):
    if _isNotCorrect(object=ruleobj, types=["Rule"]):
        return None
    
    nfac = ruleobj.getNFacies()
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
    tableobj: object containing the Table
    ranks: designates the ranks of the variable (0: ordinate; 1: abscissae [or regular]) 
    *args, **kwargs: arguments passed to _ax_table
    '''
    ax = _getNewAxes()
    return _ax_table(ax, tableobj=tableobj, icols=ranks, *args, **kwargs)
    
def _ax_table(ax, tableobj, icols, fmt='ok', flagLegend=False, **kwargs):
    '''
    Plotting the contents of a Table (argument 'table')
    ax: matplotlib.Axes
    icols: designates the ranks of the variable (0: ordinate; 1: abscissae [or regular]) 
    fmt: designates [marker][line][color] information
    **kwargs
    '''
    if _isNotCorrect(object=tableobj, types=["Table"]):
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
    **kwargs : arguments passed to _ax_mesh
    """
    ax = _getNewAxes() 
    return _ax_mesh(ax, meshobj=meshobj, *args, **kwargs)

def _ax_mesh(ax, meshobj, 
              flagEdge=True, flagFace=False, flagApex=False, 
              facecolor="yellow", edgecolor="blue", linewidth=1,
              **kwargs):
    if _isNotCorrect(object=meshobj, types=["Mesh", "MeshETurbo", "MeshEStandardExt", "DbMeshTurbo", "DbMeshStandard"]):
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

def baseMap(db, crsFrom="EPSG:4326", crsTo="EPSG:3857", 
            box=None, flagProj = False, color='blue', size=10, 
            *args, **kwargs):
    '''
    Plotting a variable from a Db using BaseMap
    '''
    ax = _getNewAxes()

    return _ax_baseMap(ax, db, crsFrom, crsTo, box,   flagProj, 
                        color, size, *args, **kwargs)

def _ax_baseMap(ax, db, crsFrom="EPSG:4326", crsTo="EPSG:3857", 
                 box=None, flagProj = False, color='blue', size=10,
                 *args, **kwargs):
    '''
    Plotting a variable from a Db using BaseMap

    db: Db defining the data to be plotted
    crsFrom Input projection characteristics
    crsTo Output projection charactieristics
    box Optional VVD for bounds (Dimension: [ndim][2])
    flagProj True if projection should be involved (using crsX)
    color Color used for displaying the data samples
    size Size used for displaying the data samples
    **kwargs : arguments passed to matplotlib.pyplot.pcolormesh

    Note: to add a basemap, simply add a sentence such as:
        import contextily as ctx
          ctx.add_basemap(ax1, source=ctx.providers.OpenStreetMap.Mapnik)
        on the Axis returned by this function

    Note: the dependency to projection (and therefore to geopanda) is conditional
    '''
    # Draw the data points
    if box is not None:
        pts = db.getAllCoordinatesMat(box).toTL()
    else:
        pts = db.getAllCoordinatesMat().toTL()

    if len(pts) > 0:
        if flagProj:
            import gstlearn.proj as prj
            from shapely.geometry import Point
            points = [Point(i) for i in pts]
            data = prj.projGP(points, crsFrom, crsTo)
            data.plot(ax=ax, color=color, markersize=size)
        else:
            plt.scatter(pts[:,0], pts[:,1], c=color, s=size)
    #         if literal:
    #             plt.annotate()
    # for i, txt in enumerate(labval):
    #     if not np.isnan(txt):
    #         ax.annotate(round(txt,2), (tabx[i], taby[i]))

    # Display bounding points (optional)
    if box is not None:
        extPoints = np.array([[box[0,0], box[1,0]],
                              [box[0,1], box[1,1]]])
        if flagProj:
            import gstlearn.proj as prj
            from shapely.geometry import Point
            geometry = [Point(xy) for xy in extPoints]
            gdf = prj.projGP(geometry, crsFrom, crsTo)
            gdf.plot(ax=ax, color='black', markersize=0.1)
        else:
            plt.scatter(extPoints[:,0], extPoints[:,1], c="white", s=0.1)

def correlation(db, namex, namey, *args, **kwargs):
    '''
    Plotting the scatter plot between two variables contained in a Db
    
    **kwargs: additional arguments passed to _ax_scatter
    '''
    ax = _getNewAxes()
    return _ax_correlation(ax, db=db, namex=namex, namey=namey, *args, **kwargs)

def _ax_correlation(ax, db, namex, namey, db2=None, 
                     asPoint = False,  flagSameAxes=False,
                     diagLine=False, diagColor="black", diagLineStyle='-',
                     bissLine=False, bissColor="red", bissLineStyle='-',
                     regrLine=False, regrColor="blue", regrLineStyle='-',
                     horizLine=False, horizColor="blue", horizLineStyle='-', hValue=0.,
                     vertLine=False, vertColor="blue", vertLineStyle='-', vValue=0.,
                     **kwargs):
    if _isNotCorrect(object=db, types=["Db", "DbGrid"]):
        return None
        
    if db2 is None:
        db2 = db
   
    if db.getNSample() != db2.getNSample():
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

    if horizLine:
        u=[xmin, xmax]
        v=[hValue, hValue]
        ax.plot(u,v,color=horizColor,linestyle=horizLineStyle)
        
    if vertLine:
        u=[vValue, vValue]
        v=[ymin, ymax]
        ax.plot(u,v,color=vertColor,linestyle=vertLineStyle)
        
    ax.decoration(xlabel = db.getName(namex)[0], ylabel = db.getName(namey)[0])

    return ax

def hscatter(db, namex, namey, varioparam, ilag, *args, **kwargs):
    '''
    Plotting the scatter plot between two variables contained in a Db
    
    **kwargs: additional arguments passed to _ax_hscatter
    '''
    ax = _getNewAxes()
    return _ax_hscatter(ax, db=db, namex=namex, namey=namey, varioparam=varioparam, ilag=ilag, 
                         *args, **kwargs)
    
def _ax_hscatter(ax, db, namex, namey, varioparam, ilag=0, idir=0, 
                     asPoint = False,  flagSameAxes=False,
                     diagLine=False, diagColor="black", diagLineStyle='-',
                     bissLine=False, bissColor="red", bissLineStyle='-',
                     **kwargs):
    if _isNotCorrect(object=db, types=["Db", "DbGrid"]):
        return None
        
    res = gl.hscatterPairs(db, namex, namey, varioparam, ilag, idir)
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
    ax = _getNewAxes()
    return _ax_anam(ax, anamobj=anamobj, *args, **kwargs)
    
def _ax_anam(ax, anamobj, color='blue', linestyle='-', flagLegend=False):
    
    if _isNotCorrect(object=anamobj, types=["Anam","AnamHermite"]):
        return None

    res = anamobj.sample()
    
    ax = _ax_XY(ax, res.getY(), res.getZ(),
                 flagLegend=flagLegend, color=color, linestyle=linestyle,
                 label='Anamorphosis')
    ax.geometry(xlim = res.getAylim(), ylim=res.getAzlim())
    ax.decoration(xlabel="Gaussian values", ylabel="Raw values")
    
    return ax

def neigh(neighobj, grid=None, node=0, *args, **kwargs):
    ax = _getNewAxes()
    return _ax_neigh(ax, neighobj=neighobj, grid=grid, node=node, *args, **kwargs)

def _ax_neigh(ax, neighobj, grid=None, node=0, flagCell=False, flagZoom=False, **kwargs):

    if _isNotCorrect(object=neighobj, types=["NeighMoving", "NeighCell"]):
        return None
    if _isNotCorrect(object=grid, types=["Db", "DbGrid"]):
        return None

    if grid is None:
        return None

    # Identify target location
    target = grid.getSampleCoordinates(node)
    
    # Represent the target location
    _ax_sample(ax, target, **kwargs)
    
    # Represent the edge of the target (if block)
    if flagCell and grid.isGrid():
        _ax_curve(ax, grid.getCellEdges(node), **kwargs)
    
    # Represent the Neighborhood Ellipsoid
    if neighobj.getType() == gl.ENeigh.MOVING:
        _ax_curve(ax, neighobj.getEllipsoid(target), **kwargs)
    
        # Represent the Angular sectors
        if neighobj.getFlagSector():
            _ax_multiSegments(ax, target, neighobj.getSectors(target), **kwargs)
        
        # Zoom to the Maximum radius circle (optional)
        if flagZoom:
            limits = neighobj.getZoomLimits(target)
            ax.set_xlim(limits[0])
            ax.set_ylim(limits[1])
    
def _ax_neighWeights(ax, res, flagWeights=True, 
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
    '''
    Generic Plot function which can be used whatever its first argument 'object'.
    '''
    res = None
    filetype = type(object).__name__

    if filetype == "Db":
        if name2 is None:
            res = symbol(object, **kwargs)
        else:
            res = correlation(object, namex=name1, namey=name2, **kwargs)
            
    elif filetype == "DbGrid":
        res = raster(object, name1, **kwargs)
    
    elif filetype == "DbLine":
        res = line(object, **kwargs)

    elif filetype == "DbGraphO":
        res = graphO(object, **kwargs)

    elif filetype == "DbMeshTurbo" or filetype == "DbMeshStandard" or filetype == "Mesh" or filetype == "MeshETurbo":
        res = mesh(object, **kwargs)
 
    elif filetype == "Vario":
        res = variogram(object, **kwargs)
    
    elif filetype == "Model":
        res = model(object, **kwargs)
    
    elif filetype == "Rule":
        res = rule(object, **kwargs)
    
    elif filetype == "AnamHermite":
        res = anam(object, **kwargs)
    
    elif filetype == "Table":
        res = table(object, ranks, **kwargs)

    elif filetype == "Faults":
        res = fault(object, **kwargs)

    elif filetype == "Polygons":
        res = polygon(object, **kwargs)

    elif filetype == "NeighMoving":
        res = neigh(object, **kwargs)
        
    else:
        print("Unknown type (in function 'plot'):",filetype)
    return res

def plotFromNF(filename, name1=None, name2=None, ranks=None, **kwargs):
    '''
    Generic function to plot the contents of any NF function
    '''
    res = None
    filetype = gl.ASerializable.getFileIdentity(filename)
    if filetype == "":
        exit()

    if filetype == "Db":
        db = gl.Db.createFromNF(filename,False)
        res = symbol(db, **kwargs)
            
    elif filetype == "DbGrid":
        dbgrid = gl.DbGrid.createFromNF(filename,False)
        res = raster(dbgrid, **kwargs)
    
    elif filetype == "Vario":
        vario_item = gl.Vario.createFromNF(filename,False)
        res = variogram(vario_item, **kwargs)
    
    elif filetype == "Model":
        model_item = gl.Model.createFromNF(filename,False)
        res = model(model_item, **kwargs)
    
    elif filetype == "Rule":
        rule_item = gl.Rule.createFromNF(filename,False)
        res = rule(rule_item, **kwargs)
    
    elif filetype == "Table":
        table_item = gl.Table.createFromNF(filename,False)
        res = table(table_item,ranks, **kwargs)

    elif filetype == "Polygons":
        polygon_item = gl.Polygons.createFromNF(filename,False)
        res = polygon(polygon_item, **kwargs)
        
    else:
        print("Unknown type:",filetype)

    return res

## ------------------------------------------ ##
## Add plot functions as methods of the class ##
## ------------------------------------------ ##
import gstlearn.plot as gp

# New style attribute setting functions
setattr(plt.Axes, "decoration",    gp._ax_decoration)
setattr(plt.Axes, "geometry",      gp._ax_geometry)

# Functions considered as members of the Axis class
setattr(plt.Axes, "polygon",       gp._ax_polygon)
setattr(plt.Axes, "rule",          gp._ax_rule)
setattr(plt.Axes, "fault",         gp._ax_faults)
setattr(plt.Axes, "anam",          gp._ax_anam)
setattr(plt.Axes, "grid1D"   ,     gp._ax_grid1D)
setattr(plt.Axes, "curve",         gp._ax_curve)
setattr(plt.Axes, "sortedcurve",   gp._ax_sortedCurve)
setattr(plt.Axes, "multiSegments", gp._ax_multiSegments)
setattr(plt.Axes, "histogram",     gp._ax_histogram)
setattr(plt.Axes, "correlation",   gp._ax_correlation)
setattr(plt.Axes, "hscatter",      gp._ax_hscatter)
setattr(plt.Axes, "table",         gp._ax_table)

setattr(plt.Axes, "model",         gp._ax_model)
setattr(plt.Axes, "mesh",          gp._ax_mesh)
setattr(plt.Axes, "variogram",     gp._ax_variogram)
setattr(plt.Axes, "varmod",        gp._ax_varmod)

setattr(plt.Axes, "neigh",         gp._ax_neigh)
setattr(plt.Axes, "neighWeights",  gp._ax_neighWeights)

setattr(plt.Axes, "symbol",        gp._ax_symbol)
setattr(plt.Axes, "literal",       gp._ax_literal)
setattr(plt.Axes, "gradient",      gp._ax_gradient)
setattr(plt.Axes, "tangent",       gp._ax_tangent)
setattr(plt.Axes, "baseMap",       gp._ax_baseMap)
setattr(plt.Axes, "raster",        gp._ax_raster)
setattr(plt.Axes, "isoline",       gp._ax_isoline)
setattr(plt.Axes, "cell",          gp._ax_cell)
setattr(plt.Axes, "box",           gp._ax_box)
setattr(plt.Axes, "XY",            gp._ax_XY)

setattr(plt.Figure, "varmod",      gp._fig_varmod)
setattr(plt.Figure, "variogram",   gp._fig_variogram)
setattr(plt.Figure, "model",       gp._fig_model)
setattr(plt.Figure, "decoration",  gp._fig_decoration)

setattr(gl.DbGrid, "plot",         gp.raster)
