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

#Set of global values
default_dims = [[5,5], [8,8]]
default_xlim = [ None, None ]
default_ylim = [ None, None ]
default_aspect = [ 'auto', 1 ]

def setDefault(mode, dims=None, xlim=None, ylim=None, aspect=None):
    global default_dims
    global default_xlim
    global default_ylim
    global default_aspect
        
    if dims is not None:
        default_dims[mode] = dims
    if xlim is not None:
        default_xlim[mode] = xlim
    if ylim is not None:
        default_ylim[mode] = ylim
    if aspect is not None:
        default_aspect[mode] = aspect

def printDefault():
    for mode in range(2):
        if mode == 0:
            print("Non geographical defaults:")
        else:
            print("Geographical defaults:")
            
        if default_dims[mode] is not None:
            print("- Figure dimensions =", default_dims[mode])
        else:
            print("- Figure dimensions (not defined)")
        if default_xlim[mode] is not None:
            print("- Limits along X =",default_xlim[mode])
        else:
            print("- Limits along X (not defined)")
        if default_ylim[mode] is not None:
            print("- Limits along Y =",default_ylim[mode])
        else:
            print("- Limits along Y (not defined)")
        if default_aspect[mode] is not None:
            print("- Aspect =",default_aspect[mode])
        else:
            print("- Aspect (not defined)")
        
def get_cmap(n, name='gist_rainbow'):
    '''
    Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
    RGB color; the keyword argument name must be a standard mpl colormap name.
    '''
    return plt.cm.get_cmap(name, n)
    
def selectItems(nvalues, sitem=-1):
    outs = range(0, nvalues)
    nout = nvalues
    if sitem >= 0:
        outs = range(sitem, sitem+1)
        nout = 1
    return outs, nout

def geometry(ax=None, dims=None, xlim=None, ylim=None, aspect=None):
    '''
    Set the default values for the geometrical parameters for one or a set of Axes
    
    ax: matplotlib.Axes (necessary when used as a method of the class)
    dims: Extension of graphic Axes
    xlim: Range of values along the X-axis
    ylim: Range of values along the Y-axis
    aspect: Y/X ratio
    
    Remark: When 'ax' designates a set of Axes, parameters are applied to all of them.
    '''
    
    if ax is None:
        print("This method is only available in 'overlay' mode: 'ax' must be provided")
        return None
    
    if dims is not None:
        if is_array(dims, 2):
            if type(ax) == np.ndarray:
                ax[0,0].figure.set_size_inches(dims[0]*ax.shape[0], 
                                               dims[1]*ax.shape[1])
            else:
                ax.figure.set_size_inches(dims[0], dims[1])
        else:
            print("'dims' should be [a,b]. Ignored")
        
    if xlim is not None:
        if is_array(xlim, 2):
            if type(ax) == np.ndarray:
                ax[ix,iy].set_xlim(left = xlim[0], right = xlim[1])
            else:
                ax.set_xlim(left = xlim[0], right = xlim[1])
        else:
            print("'xlim' should be [a,b] or [None,b] or [a,None]. Ignored")
    
    if ylim is not None:
        if is_array(ylim, 2):
            if type(ax) == np.ndarray:
                for ix in range(ax.shape[0]):
                    for iy in range(ax.shape[1]):
                        ax[ix,iy].set_ylim(bottom = ylim[0], top = ylim[1])
            else:
                ax.set_ylim(bottom = ylim[0], top = ylim[1])
        else:
           print("'ylim' should be [a,b] or [None,b] or [a,None]. Ignored")
        
    if aspect is not None:
        if type(ax) == np.ndarray:
            for ix in range(ax.shape[0]):
                for iy in range(ax.shape[1]):
                    ax[ix,iy].set_aspect(aspect)
        else:
            ax.set_aspect(aspect)

def decoration(ax=None, xlabel=None, ylabel=None, title=None, **kwargs):
    '''
    Add the decoration to a figure.
    
    Parameters
    ----------
    The procedure depends whether 'ax' is already defined or not.
    ax: matplotlib.Axes (necessary when used as a method of the class)
    xlabel: label along the horizontal axis
    ylabel: label along the vertical axis
    title: title contents (for the main for a collection of Axes)
    '''
    if ax is None:
        print("This method is only available in 'overlay' mode: 'ax' must be provided")
        return None

    if title is not None:
        if type(ax) == np.ndarray:
            ax[0,0].figure.suptitle(title, **kwargs)
        else:
            ax.set_title(title, **kwargs)

    if xlabel is not None:
        if type(ax) == np.ndarray:
            print("decoration() cannot be used when 'ax' is an array. Ignored")
        else:
            ax.set_xlabel(xlabel)
    if ylabel is not None:
        if type(ax) == np.ndarray:
            print("decoration() cannot be used when 'ax' is an array. Ignored")
        else:
            ax.set_ylabel(ylabel)

def getNewAxes(ax=None, mode=0, nx=1, ny=1, sharex=False, sharey=False):
    ''' Creates a new figure (possibly containing multiple subplots)
    
        Parameters
        ----------
        ax: Axes description (or array of them). See remarks.
        mode: Type of Global value (used for default assignment). 0 for Non_Geographical; 1 for Geographical
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
        if len(plt.get_fignums()) == 0:
            
            # Axes is None and no Figure already exists. Create it
            fig, ax = plt.subplots(nx, ny, squeeze=False, sharex=sharex, sharey=sharey)
            
            if is_array(ax, 1):
                ax = ax[0,0]
                
            # Apply the Global Geometry parameters (when defined)
            geometry(ax,
                     dims = default_dims[mode], 
                     xlim = default_xlim[mode], 
                     ylim = default_ylim[mode], 
                     aspect = default_aspect[mode])
            
        else:
            # Axes is None but a figure already exists, return the (last) Axes of Figure
            ax = plt.gca()
    
    if is_array(ax, 1):
        ax = ax[0,0]
        
    return ax

def shape_Nsubplots(N):
    '''
    Calculates numbers of lines and columns for N subplots. 
    If N is a perfect square, then nlines = ncols = sqrt(N). 
    Else, the number of columns increase first.
    '''
    nlines = np.floor(np.sqrt(N))
    ncols = np.ceil(N/nlines)
    return int(nlines), int(ncols)

def is_array(tab, ndim=None):
    '''
    Check if the input argument is an array (of dimension 'ndim' when defined) or a scalar
    tab:  Argument to be checked
    ndim: Requested dimension. The dimension check is not performed if None 
    '''
    if not hasattr(tab, "__len__"):
        return False
    
    if (ndim is not None) and (len(tab) is not ndim):
        return False
    
    return True

def addColorbar(im, ax):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, ax=ax, cax=cax)
    return cbar

def getDefinedValues(db, name, posX=0, posY=1, corner=None, usesel=True, 
                     compress=False, asGrid=True, 
                     flagConvertNanToZero=False):

    if db.isGrid() and asGrid:
        if corner is None:
            corner = np.zeros(db.getNDim())
        
        if db.getNDim() == 1:
            tab = db.getColumn(name, usesel, False)
        else:
            tab = db.getOneSlice(name, posX, posY, corner, usesel)
    else:
        tab = db.getColumn(name, usesel, compress)
    tab = np.array(tab).transpose()

    if flagConvertNanToZero:
        tab[np.isnan(tab)] = 0
    else:
        tab = ma.array(tab,mask=np.isnan(tab))
    
    if compress:
        tab = tab[np.logical_not(np.isnan(tab))]
        
    return tab

def getBiDefinedValues(db1, name1, name2, db2, usesel=True):
    tabx = db1.getColumn(name1, usesel)
    tabx = np.array(tabx).transpose()
    
    taby = db2.getColumn(name2, usesel)
    taby = np.array(taby).transpose()
    
    sel  = np.logical_not(np.logical_or(np.isnan(tabx), np.isnan(taby)))
    tabx = tabx[sel]
    taby = taby[sel]
    return tabx, taby

def getFileIdentity(filename):
    type = gl.Aserializable.getFileIdentity(filename)
    print(type)
    return type

def varioElem(ax=None, vario=None, ivar=0, jvar=0, idir=0, hmax=None, show_pairs = False,
              var_color='black', var_linestyle='dashed', 
              flagDrawVariance = True, flagLabelDir=False, flagLegend=False, 
              label=None, **kwargs):
    """
    Plot a single experimental variogram (one direction and fixed pair of variable(s)).
    
    Parameters
    ----------
    ax: matplotlib.Axes (necessary when used as a method of the class)
    vario : experimental variogram to be represented (gstlearn.Vario).
    ivar, jvar : Indices of the variables for the variogram to be represented (the default is 0).
    idir : Index of the direction of the variogram to be represented (the default is 0).
    hmax : Maximum distance to be represented.
    show_pairs : Flag for annotating the number of pairs for each lag on the plot (the default is False).
    var_color : color of the horizontal line representing the sill (the default is 'black').
    var_linestyle : linestyle of the horizontal line representing the sill (the default is 'dashed').
    label : Label to be drawn (constructed if not provided)
    flagDrawVariance : Flag to add the variance (default is True)    
    flagLabelDir : Encode the direction in the label (when constructed)
    flagLegend : Flag to display the axes legend.
    **kwargs : arguments passed to matplotlib.pyplot.plot

    Returns
    -------
    ax : axes where the variogram is represented
    """
    if vario is None:
        print("'vario' is compulsory")
        return None

    ax = getNewAxes(ax,0)

    if label is None:
        if flagLabelDir:
            label = "vario dir={}".format(np.round(vario.getCodirs(idir),3))
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
        ax.hlines(vario.getVar(ivar,jvar), 0, hmax, var_color, var_linestyle)
    
    # Representing the number of pairs (optional)
    if show_pairs:
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

def varmod(vario, model=None, ivar=-1, jvar=-1, idir=-1,
           nh = 100, hmax = None, show_pairs=False, asCov=False, 
           var_color='black', var_linestyle="dotted",
           env_color='black', env_linestyle="dotted",
           cmap=None, flagLegend=False, axs=None, 
           **kwargs):
    """
    Construct a figure for plotting experimental variogram(s) and model.
    
    Parameters
    ----------
    vario : experimental variogram to be represented
    model : optional, variogram model
    ivar, jvar : Indices of the variables for the variogram to be represented. If -1 (default), all 
                 variables are selected and all the simple and crossed variograms are represented.
    idir : Index of the direction of the variogram to be represented. If -1 (default) all available
           directions are selected and multidirectional variograms are represented.
    var_color, var_linestyle: parameters for representing variance-covariance line
    env_color, env_linestyle: parameters for representing coregionalization envelop
    nh : number of points between 0 and hmax where the model variogram is calculated (default is 100).
    hmax : Maximum distance to be represented.
    cmap : Optional Color scale
    flagLegend : Flag to display the axes legend.
    axs : Reference for the plot(s) within the figure. If None (default),
          it creates a new figure (with multiple axes for multivariate variograms).

    **kwargs : arguments passed to matplotlib.pyplot.plot for all variograms plotted (not models!)
    
    Returns
    -------
    ax : axes where the variograms are represented
    """
    color_in_kwargs = 'color' in kwargs
    
    if hmax is None:
        hmax = vario.getHmax(ivar, jvar, idir)
        
    ndir = vario.getDirectionNumber()
    nvar = vario.getVariableNumber()
    cols = get_cmap(ndir,cmap)
    
    ndirUtil, ivarD = selectItems(ndir, idir)
    ivarUtil, ivarN = selectItems(nvar, ivar)
    jvarUtil, jvarN = selectItems(nvar, jvar)
        
    axs = getNewAxes(axs, 0, nx=ivarN, ny=jvarN)
        
    if ndir > 1:
        flagLabelDir = True
    else:
        flagLabelDir = False
        
    # Loop on the variables  
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
                
                varioElem(ax, vario, iv, jv, idirUtil, 
                          show_pairs=show_pairs, hmax=hmax,
                          var_color=var_color, var_linestyle=var_linestyle,  
                          flagLabelDir=flagLabelDir, flagLegend=flagLegend, **kwargs)

                # Plotting the Model (optional)
                if model is not None:
                    codir = vario.getCodirs(idirUtil)
                    modelElem(ax, model, ivar=iv, jvar=jv, codir=codir, 
                              hmax=hmax, nh=nh, asCov=asCov,
                              env_color=env_color, env_linestyle=env_linestyle, 
                              flagLabelDir=flagLabelDir, flagLegend=flagLegend, **kwargs)

            ax.autoscale(True)
            
            if vario.drawOnlyPositiveX(iv, jv):
                ax.set_xlim(left=0)
            if vario.drawOnlyPositiveY(iv, jv):
                ax.set_ylim(bottom=0)
    
    if is_array(axs, 1):
        return axs[0,0]
    else:
        return axs

def vario(vario, ivar=0, jvar=0, idir=0,
          var_color='black', var_linestyle='dashed', hmax=None,  
          cmap = None, flagLegend=False, 
          axs = None, **kwargs):
    """
    Plot experimental variogram(s) (can be multidirectional and multivariable or selected ones).
    
    Parameters
    ----------
    vario : experimental variogram to be represented (gstlearn.Vario).
    ivar, jvar : Indices of the variables for the variogram to be represented. If -1 (default), all 
                 variables are selected and all the simple and crossed variograms are represented.
    idir : Index of the direction of the variogram to be represented. If -1 (default) all available
           directions are selected and multidirectional variograms are represented.
    var_color, var_linestyle : parameters for representing variance-covariance line
    hmax : Maximum distance to be represented.
    cmap : Optional Color scale
    flagLegend : Flag to display the axes legend.
    axs : Reference for the plot(s) within the figure. If None (default),
          it creates a new figure (with multiple axes for multivariate variograms).
    **kwargs : arguments passed to matplotlib.pyplot.plot for all variograms plotted

    Returns
    -------
    ax : axes where the variograms are represented
    """
    axs = varmod(vario, ivar=ivar, jvar=jvar, idir=idir, 
                 var_color=var_color, var_linestyle=var_linestyle, 
                 hmax=hmax, cmap=cmap, 
                 flagLegend=flagLegend, axs=axs, 
                 **kwargs)
    
    return axs

def modelElem(ax = None, model = None, ivar=0, jvar=0, codir=None, vario=None, idir=0,
              nh = 100, hmax = None, asCov=False,
              env_color='black', env_linestyle='dashed',
              label=None, flagLabelDir=False, flagEnvelop = True, flagLegend=False, 
              **kwargs):
    """
    Construct a Layer for plotting a model
    
    Parameters
    ----------
    ax: matplotlib.Axes (necessary when used as a method of the class)
    model : variogram model to be represented (gstlearn.Model).
    ivar, jvar : Indices of the variables for the variogram to be represented (the default is 0).
    codir : Vector of the direction of the variogram to be represented. The default is the unit 
            vector in the first space dimension.
    vario, idir: Vario information used to set the direction (when codir is not provided)
    env_color, env_linestyle : color and linestyle for correlation envelop 
    nh : number of points between 0 and hmax where the model variogram is calculated (default is 100).
    flagEnv : flag for representing the correlation envelop (the default is True)
    hmax : Maximum distance to be represented. By default: 3 times the maximum range of the
           basic structures, or 1 if no range is defined.
    asCov : Present the Model as a Covariance (rather than as a Variogram)
    label: Label displayed in the Legend (constructed if not provided)
    flagLabelDir : Add the direction vector to the label (if constructed)
    flagEnvelop: Represent the coregionalization envelop (in multivariate case only)
    flagLegend : Flag to display the axes legend.
    """
    if model is None:
        print("'model' is compulsory")
        return None

    if codir is None:
        if vario is None:
            codir = [0] * model.getDimensionNumber()
            codir[0] = 1
        else:
            codir = vario.getCodirs(idir)
            
    # if hmax not specified = 3*maximum range of the model's basic structures
    if hmax is None:
        hmax = 0
        for icova in range(model.getCovaNumber()):
            range_max = np.max(model.getCova(icova).getRanges())
            if 3*range_max > hmax:
                hmax = 3*range_max
    if hmax == 0: # if the model has no range defined
        hmax = 1
            
    ax = getNewAxes(ax, 0)
    
    if label is None:
        if flagLabelDir:
            label = "model dir={}".format(np.round(codir,3))
        else:
            label = "model"

    istart = 0
    for i in range(model.getCovaNumber()):
        if model.getCovName(i) == 'Nugget Effect':
            istart = 1 # do not plot the first lag (h=0) for nugget effect (discontinuity)
     
    # Represent the Model 
    hh = np.linspace(0, hmax, nh+1)
    gg = model.sample(hh, ivar, jvar, codir, 0, asCov=asCov)
    res = ax.plot(hh[istart:], gg[istart:], label=label, **kwargs)
    
    # Represent the coregionalization envelop (optional)
    if ivar != jvar and flagEnvelop:
        ggp = model.sample(hh, ivar, jvar, codir, 1, asCov=asCov)
        ax.plot(hh[istart:], ggp[istart:], c = env_color, linestyle = env_linestyle)
        ggm = model.sample(hh, ivar, jvar, codir,-1, asCov=asCov)
        ax.plot(hh[istart:], ggm[istart:], c = env_color, linestyle = env_linestyle)
    
    # Draw the Legend (optional)
    if flagLegend:
        ax.legend()
        
    return res

def model(model = None, ivar=0, jvar=0, codir=None, vario=None, idir=0,
          nh = 100, hmax = None, asCov=False,
          env_color='black', env_linestyle='dashed',
          label=None, flagLabelDir=False, flagEnvelop = True, flagLegend=False, 
          ax = None, **kwargs):
    
    ax = getNewAxes(ax, 0)
    
    modelElem(ax, model = model, ivar=ivar, jvar=jvar, codir=codir, 
              vario=vario, idir=idir,
              nh = nh, hmax = hmax, asCov=asCov,
              env_color=env_color, env_linestyle=env_linestyle,
              label=label, flagLabelDir=flagLabelDir, flagEnvelop = flagEnvelop, 
              flagLegend=flagLegend, 
              **kwargs)
    
    return ax

def readCoorPoint(db, coorX_name=None, coorY_name=None, 
                  usesel=True, posX=0, posY=1):
    
    # Extracting coordinates
    if coorX_name is not None:
        tabx = db.getColumn(coorX_name, usesel)
    else:
        if db.getNDim() > 0:
            tabx = db.getCoordinates(posX,usesel)
            
    if coorY_name is not None:
        taby = db.getColumn(coorY_name, usesel)
    else:
        if db.getNDim() > 1:
            taby = db.getCoordinates(posY,usesel)
    
    if len(tabx) <= 0 or len(taby) <= 0:
        return None
    
    return tabx, taby
    
def pointSymbol(ax=None, db=None, name_color=None, name_size=None, 
                coorX_name=None, coorY_name=None, usesel=True, 
                c='r', s=20, sizmin=10, sizmax=200, flagAbsSize=False, flagCst=False,
                flagLegend=True, legendName=None, posX=0, posY=1, **kwargs):
    '''
    Construct a Layer for plotting a point data base, with optional color and size variables
    
    ax: matplotlib.Axes (necessary when used as a method of the class)
    db: Db containing the variable to be plotted
    name_color: Name of the variable containing the color per sample
    name_size: Name of the variable containing the size per sample
    coorX_name: Name of the variable standing for X coordinate 
    coorY_name: Name of the variable standing for Y coordinate 
    usesel : Boolean to indicate if the selection has to be considered
    c: Constant color (used if 'name_color' is not defined)
    s: Constant size (used if 'name_size' is not defined)
    sizmin: Size corresponding to the smallest value (used if 'name_size' is defined)
    sizmax: Size corresponding to the largest value (used if 'name_size' is defined)
    flagAbsSize: Represent the Absolute value in Size representation
    flagCst: When True, the size is kept constant (equal to 's')
    flagLegend: Flag for representing the Legend
    legendName: Title for the Legend
    posX: rank of the first coordinate
    posY: rank of the second coordinate
    **kwargs : arguments passed to matplotllib.pyplot.scatter
    '''
    if db is None:
        print("'db' is compulsory")
        return None
    
    ax = getNewAxes(ax, 1)

    name = ''
    
    # Read the coordinates
    tabx, taby = readCoorPoint(db, coorX_name, coorY_name, usesel, posX, posY)
    
    # Color of symbol
    if name_color is not None:
        colval = getDefinedValues(db, name_color, 0, 1, None, usesel, 
                                  compress=True, asGrid=False, 
                                  flagConvertNanToZero=True)
        name = name + ' ' + name_color
    else:
        colval = c

    # Size of symbol
    if name_size is not None:
        sizval = getDefinedValues(db, name_size, 0, 1, None, usesel, 
                                  compress=True, asGrid=False, 
                                  flagConvertNanToZero=True)
        if not flagCst:
            if flagAbsSize:
                sizval = np.absolute(sizval)
            m = np.nanmin(np.absolute(sizval))
            M = np.nanmax(np.absolute(sizval))
            if M > m:
                sizval = (sizmax - sizmin) * (np.absolute(sizval) - m) / (M-m) + sizmin
                
            name = name + ' ' + name_size
        else:
            sizval = s
    else:
        sizval = s

    if len(ax.get_title()) <= 0:
        ax.decoration(title = name)
    
    res = ax.scatter(x = tabx, y = taby, s = sizval, c = colval, **kwargs)

    if flagLegend:
        if name_color is not None:
            addColorbar(res, ax)
    
        if name_size is not None:
            labels = lambda marker_size : (M - m)*(marker_size - sizmin)/(sizmax - sizmin) + m
            ax.legend(*res.legend_elements("sizes", num=6, func=labels), 
                      title=legendName)
         
    return res

def pointLabel(ax=None, db=None, name=None, coorX_name=None, coorY_name=None, 
               usesel=True, flagLegend=True, legendName=None, 
               posX=0, posY=1, **kwargs):
    '''
    Construct a layer for plotting a point data base, with optional color and size variables
    
    ax: matplotlib.Axes (necessary when used as a method of the class)
    db: Db containing the variable to be plotted
    name: Name of the variable containing the label per sample
    coorX_name: Name of the variable standing for X coordinate 
    coorY_name: Name of the variable standing for Y coordinate 
    usesel : Boolean to indicate if the selection has to be considered
    flagLegend: Flag for representing the Color Bar
    legendName: title of the Legend
    posX: rank of the first coordinate
    posY: rank of the second coordinate
    **kwargs : arguments passed to matplotllib.pyplot.scatter
    '''
    if (db is None) or (name is None):
        print("'db' and 'name' are compulsory")
        return None
    
    ax = getNewAxes(ax, 1)

    if len(ax.get_title()) <= 0:
        ax.decoration(title = db.getName(name)[0])
    
    # Read the coordinates
    tabx, taby = readCoorPoint(db, coorX_name, coorY_name, usesel, posX, posY)
    
    labval = getDefinedValues(db, name, 0, 1, None, usesel, 
                              compress=True, asGrid=False, 
                              flagConvertNanToZero=True)

    res = ax.scatter(x = tabx, y = taby, **kwargs)
    
    for i, txt in enumerate(labval):
        ax.annotate(round(txt,2), (tabx[i], taby[i]))
  
    return res

def pointGradient(ax=None, db=None, coorX_name=None, coorY_name=None, usesel=True, 
                  posX=0, posY=1, **kwargs):
    '''
    Construct a layer for plotting a gradient data base
    
    ax: matplotlib.Axes (necessary when used as a method of the class)
    db: Db containing the variable to be plotted
    coorX_name: Name of the variable standing for X coordinate 
    coorY_name: Name of the variable standing for Y coordinate 
    usesel : Boolean to indicate if the selection has to be considered
    '''
    if db is None:
        print("'db' is compulsory")
        return None
    
    ax = getNewAxes(ax, 1)

    if db.getLocNumber(gl.ELoc.G) <= 0:
        return None
    
    # Extracting coordinates
    tabx, taby = readCoorPoint(db, coorX_name, coorY_name, usesel, posX, posY)
    
    # Reading the Gradient components
    if db.getNDim() > 1:
        tabgx = db.getGradient(0,usesel)
        tabgy = db.getGradient(1,usesel)
    else:
        tabgy = -db.getGradient(0,usesel)
        tabgx = -np.ones(len(tabgy))

    if len(tabx) <= 0 or len(taby) <= 0 or len(tabgx) <= 0 or len(tabgy) <= 0:
        return None
    
    res = ax.quiver(tabx, taby, -tabgx, -tabgy, angles='xy', **kwargs)
            
    return res

def pointTangent(ax=None, db=None, coorX_name=None, coorY_name=None, usesel=True, 
                 posX=0, posY=1, **kwargs):
    '''
    Construct a layer for plotting a tangent data base
    
    ax: matplotlib.Axes (necessary when used as a method of the class)
    db: Db containing the variable to be plotted
    coorX_name: Name of the variable standing for X coordinate 
    coorY_name: Name of the variable standing for Y coordinate 
    usesel : Boolean to indicate if the selection has to be considered
    '''
    if db is None:
        print("'db' is compulsory")
        return None
    
    ax = getNewAxes(ax, 1)

    if db.getLocNumber(gl.ELoc.TGTE) <= 0:
        return None

    # Extracting coordinates
    tabx, taby = readCoorPoint(db, coorX_name, coorY_name, usesel, posX, posY)

    # Extract Tangent information
    tabtx = db.getTangent(0,usesel)
    tabty = db.getTangent(1,usesel)

    if len(tabx) <= 0 or len(taby) <= 0 or len(tabtx) <= 0 or len(tabty) <= 0:
        return None
    
    res = ax.quiver(tabx, taby, -tabtx, -tabty, **kwargs)
    res = ax.quiver(tabx, taby,  tabtx,  tabty, **kwargs)
            
    return res

def point(db, 
          name_color=None, name_size=None, name_label=None,
          coorX_name=None, coorY_name=None, usesel=True, 
          color='r', size=20, sizmin=10, sizmax=200, cmap=None,
          flagAbsSize=False, flagCst=False,
          flagGradient=False, colorGradient='black', scaleGradient=20,
          flagTangent=False, colorTangent='black', scaleTangent=20,
          flagLegendSymbol=False, legendSymbolName=None,
          flagLegendLabel=False, legendLabelName=None,
          posX=0, posY=1, ax=None, **kwargs):
    '''
    Construct a figure for plotting a point data base
    
    db: Db containing the variable to be plotted
    name_color: Name of the variable containing the color per sample
    name_size: Name of the variable containing the size per sample
    name_label: Name of the variable containing the label per sample
    coorX_name: Name of the variable standing for X coordinate 
    coorY_name: Name of the variable standing for Y coordinate 
    usesel : Boolean to indicate if the selection has to be considered
    color: Constant color (used if 'name_color' is not defined)
    size: Constant size (used if 'name_size' is not defined)
    sizmin: Size corresponding to the smallest value (used if 'name_size' is defined)
    sizmax: Size corresponding to the largest value (used if 'name_size' is defined)
    flagAbsSize: Represent the Absolute value in Size representation
    flagCst: When True, the size is kept constant (equel to 's')
    cmap: Optional Color scale
    flagGradient: Draw Gradient (if Gradients are defined)
    colorGradient: Color attached to the Gradient representation
    scaleGradient: Scale of the Gradient representation
    flagTangent: Draw Tangent (if Tangents are defined)
    colorTangent: Color attached to the Gradient representation
    scaleTangent: Scale of the Gradient representation
    flagLegendSymbol: Flag for representing the Color Bar (only if name_color is defined)
    legendSymbolName: Title for the Symbol Legend
    flagLegendLabel: Flag for representing the Legend for marker size (only if name_size is defined)
    legendLabelName: Title for the Label Legend
    posX: rank of the first coordinate
    posY: rank of the second coordinate

    **kwargs : arguments passed to matplotllib.pyplot.scatter
    '''
    ax = getNewAxes(ax, 1)
    
    # If no variable is defined, use the default variable for Symbol(size) representation
    # The default variable is the first Z-locator one, or the last variable in the file
    if (name_color is None) and (name_size is None) and (name_label is None):
        if db.getLocNumber(gl.ELoc.Z) > 0:
            name_size = db.getNameByLocator(gl.ELoc.Z,0)
        else : # if no Z locator, choose the last field
            name_size = db.getLastName()
            flagCst = True

    if (name_color is not None) or (name_size is not None):
        pt = pointSymbol(ax, db, name_color=name_color, name_size=name_size, 
                         coorX_name=coorX_name, coorY_name=coorY_name, usesel=usesel, 
                         c=color, s=size, sizmin=sizmin, sizmax=sizmax, 
                         flagAbsSize=flagAbsSize, flagCst=flagCst,
                         cmap=cmap, 
                         flagLegend=flagLegendSymbol, legendName=legendSymbolName,
                         posX=posX, posY=posY, 
                         **kwargs)
    
    if name_label is not None:
        tx = pointLabel(ax, db, name=name_label, 
                        coorX_name=coorX_name, coorY_name=coorY_name, 
                        usesel=usesel, 
                        flagLegend=flagLegendLabel, legendName=legendLabelName,
                        posX=posX, posY=posY, **kwargs)
        
    if flagGradient:
        gr = pointGradient(ax, db, coorX_name=coorX_name, coorY_name=coorY_name, 
                           usesel=usesel, color=colorGradient, scale=scaleGradient,
                           posX=posX, posY=posY)

    if flagTangent:
        tg = pointTangent(ax, db, coorX_name=coorX_name, coorY_name=coorY_name,
                          usesel=usesel, color=colorTangent, scale=scaleTangent,
                          posX=posX, posY=posY)
        
    return ax

def modelOnGrid(model, db, usesel=True, icov=0, color='black', scale=1,
                ax=None):
    '''
    Display the Model characteristics on a Grid
    This makes sense when the model contains some non-stationarity
    '''
    ax = getNewAxes(ax, 1)

    # Extracting coordinates
    tabx = db.getCoordinates(0,usesel)
    taby = db.getCoordinates(1,usesel)
    if len(tabx) <= 0 or len(taby) <= 0:
        return None
    
    gl.db_model_nostat(db, model, icov)
    tabR1 = db.getColumn("Nostat.Range-1", usesel)
    tabR2 = db.getColumn("Nostat.Range-2", usesel)
    tabA  = db.getColumn("Nostat.Angle-1", usesel)
    if len(tabR1) <= 0 or len(tabR2) <= 0 or len(tabA) <= 0:
        return None
    
    ax.quiver(tabx, taby, tabR2, tabR2, angles=tabA, color=color)
            
    return ax
    
def polygon(poly, facecolor='yellow', edgecolor = 'blue', 
            colorPerSet = False, flagEdge=True, flagFace=False, 
            ax=None, **kwargs):
    '''
    Construct a Figure for plotting a polygon
    **kwargs: arguments passed to matplotlib.fill
    '''
    ax = getNewAxes(ax, 1)
    
    npol = poly.getPolySetNumber()
    cols = get_cmap(npol)
    
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

def readGrid(dbgrid, name, usesel=True, 
             posx=0, posy=1, corner=None, shading = "nearest"):
    
    x0 = dbgrid.getX0(posx)
    y0 = dbgrid.getX0(posy)
    nx = dbgrid.getNX(posx)
    ny = dbgrid.getNX(posy)
    dx = dbgrid.getDX(posx)
    dy = dbgrid.getDX(posy)
    angles = dbgrid.getAngles()
    
    data = getDefinedValues(dbgrid, name, posx, posy, corner, usesel, 
                            compress=False, asGrid=True)
    data = np.reshape(data, (ny,nx))

    tr = transform.Affine2D().rotate_deg_around(x0,y0,angles[0])
    
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

def gridRaster(ax=None, dbgrid=None, name=None, usesel = True, posx=0, posy=1, corner=None, 
               flagLegend=True, **kwargs):
    '''
    Plotting a variable from a DbGrid in Raster

    ax: matplotlib.Axes (necessary when used as a method of the class)
    dbgrid: DbGrid containing the variable to be plotted
    name: Name of the variable to be represented (by default, the first Z locator, or the last field)
    usesel : Boolean to indicate if the selection has to be considered
    flagLegend: Flag for representing the Color Bar
    **kwargs : arguments passed to matplotlib.pyplot.pcolormesh
    '''
    if (dbgrid is None) or (name is None):
        print("'dbgrid' and 'name' are compulsory")
        return None
    
    ax = getNewAxes(ax, 1)
        
    if len(ax.get_title()) <= 0:
        ax.decoration(title = dbgrid.getName(name)[0])
    
    x0, y0, X, Y, data, tr = readGrid(dbgrid, name, usesel, 
                                      posx=posx, posy=posy, corner=corner)
    trans_data = tr + ax.transData
    
    res = ax.pcolormesh(X, Y, data, **kwargs)
    res.set_transform(trans_data)
    
    x1, x2, y1, y2 = x0, X[-1], y0, Y[-1]
    ax.plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], marker='', linestyle='', 
            transform=trans_data)
   
    if flagLegend:
        addColorbar(res, ax)
    
    return res
        
def gridContour(ax=None, dbgrid=None, name=None, usesel = True, 
                posx=0, posy=1, corner=None, levels=None,
                flagLegend=True, **kwargs):
    '''
    Plotting a variable (referred by its name) informed in a DbGrid

    ax: matplotlib.Axes (necessary when used as a method of the class)
    dbgrid: DbGrid containing the variable to be plotted
    name: Name of the variable to be represented (by default, the first Z locator, or the last field)
    usesel : Boolean to indicate if the selection has to be considered
    levels: Vector of isovalues to be represented
    flagLegend: Flag for representing the Color Bar (not represented if alpha=0)
    ax: Reference for the plot within the figure
    
    **kwargs : arguments passed to matplotlib.pyplot.contour
    '''
    if (dbgrid is None) or (name is None):
        print("'dbgrid' and 'name' are compulsory")
        return None
    
    ax = getNewAxes(ax, 1)

    if len(ax.get_title()) <= 0:
        ax.decoration(title = dbgrid.getName(name)[0])
    
    x0, y0, X, Y, data, tr = readGrid(dbgrid, name, usesel, 
                                      posx=posx, posy=posy, corner=corner)
    trans_data = tr + ax.transData
    
    res = ax.contour(X, Y, data, levels, **kwargs)
    
    if flagLegend:
        h1,l1 = res.legend_elements()
        ax.legend([h1[0]], ["Contour"])
        
    return res

def grid(dbgrid, name_raster = None, name_contour = None, usesel = True, 
         posx=0, posy=1, corner=None, 
         flagLegendRaster=False, flagLegendContour=False,
         levels=None, ax=None, **kwargs):
    '''
    Plotting a variable (referred by its name) informed in a DbGrid

    dbgrid: DbGrid containing the variable to be plotted
    name_raster: Name of the variable to be represented as raster
    name_contour: Name of the variable tp be represented as contours
    usesel : Boolean to indicate if the selection has to be considered
    flagLegendColor: Flag for representing the Color Bar (not represented if alpha=0)
    **kwargs : arguments passed to matplotlib.pyplot.pcolormesh
    '''
    if not(dbgrid.isGrid()):
        print("This function is dedicated to Grid Db and cannot be used here")
        return None;
    
    ax = getNewAxes(ax, 1)
    
    # If no variable is defined, use the default variable for Raster representation
    # The default variable is the first Z-locator one, or the last variable in the file
    if (name_raster is None) and (name_contour is None):
        if dbgrid.getLocNumber(gl.ELoc.Z) > 0:
            name_raster = dbgrid.getNameByLocator(gl.ELoc.Z,0)
        else : # if no Z locator, choose the last field
            name_raster = dbgrid.getLastName()

    if name_raster is not None:
        rs = gridRaster(ax, dbgrid = dbgrid, name = name_raster, usesel = usesel,  
                        posx=posx, posy=posy, corner=corner, 
                        flagLegend=flagLegendRaster,
                        **kwargs)
    
    if name_contour is not None:
        ct = gridContour(ax, dbgrid = dbgrid, name = name_contour, usesel = usesel, 
                         posx=posx, posy=posy, corner=corner, levels=levels, 
                         flagLegend=flagLegendContour, 
                         **kwargs)
    
    return ax

def grid1D(dbgrid, name = None, usesel = True, flagLegendColor=True,
           color='black',flagLegend=False, label='curve',
           ax=None, **kwargs):
    '''
    Plotting a variable (referred by its name) informed in a DbGrid

    dbgrid: DbGrid containing the variable to be plotted
    name: Name of the variable to be represented (by default, the first Z locator, or the last field)
    usesel : Boolean to indicate if the selection has to be considered
    flagLegendColor: Flag for representing the Color Bar
    ax: Reference for the plot within the figure
    
    **kwargs : arguments passed to matplotlib.pyplot.curve
    '''
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
    
    ax = getNewAxes(ax, 1)
        
    x0 = dbgrid.getX0(0)
    nx = dbgrid.getNX(0)
    dx = dbgrid.getDX(0)
    
    tabx = dbgrid.getColumnByLocator(gl.ELoc.X, 0, usesel)
    data = getDefinedValues(dbgrid, name, 0, 1, None, usesel, 
                            compress=False, asGrid=True)

    curve(data1=tabx, data2=data, color=color, flagLegend=flagLegend, 
          ax=ax, **kwargs)

    ax.decoration(title = dbgrid.getName(name)[0])
        
    return ax

def hist_tab(val, ax = None, **kwargs):
    '''
    Plotting the histogram of an array (argument 'val')
    
    kwargs : arguments passed to matplotlib.pyplot.hist
    '''
    ax = getNewAxes(ax, 0)
        
    ax.hist(val, **kwargs)
    
    return ax
    
def hist(db, name, ax=None, usesel=True, **kwargs):
    '''
    Plotting the histogram of a variable contained in a Db
    
    kwargs : arguments passed to matplotlib.pyplot.hist
    '''
    db.useSel = usesel
    val = db[name]
    if len(val) == 0:
        return None
    
    ax = hist_tab(val, ax=ax, **kwargs)
    
    ax.decoration(title = db.getName(name)[0])
        
    return ax

def sortedcurve(tabx, taby, color='black', flagLegend=False,
                ax=None, **kwargs):
    '''
    Plotting a set of points after they have been sorted in increasing X
    '''
    # Account for possible 'nan'  values
    mask = np.logical_and(np.isfinite(tabx), np.isfinite(taby))
    stabx = tabx[mask]
    staby = taby[mask]
    
    # Indices of the sorted elements of stabx
    indices = np.argsort(stabx)
    ax = curve(stabx[indices], staby[indices], color=color, 
               flagLegend=flagLegend, ax=ax, **kwargs)
    
    return ax
    
def curve(data1, data2=None, icas=1, color0='black',flagLegend=False, 
          ax=None, **kwargs):
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
    color = kwargs.setdefault('color', color0)
    label = kwargs.setdefault('label', 'curve')
    
    if len(data1) == 0:
        return None

    ax = getNewAxes(ax, 0)
    
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

def multisegments(center, data, color='black',flagLegend=False, label="segments",
                  ax=None, **kwargs):
    '''
    Plotting a set of segments joining 'center' to any of vertices
    stored in 'data'.
    **kwargs : arguments passed to matplotlib.pyplot.plot
    '''
    color = kwargs.setdefault('color', color)
    label = kwargs.setdefault('label', label)
        
    if len(data) == 0:
        return None
    
    ax = getNewAxes(ax, 0)
    
    nseg = len(data[0])
    
    for iseg in range(nseg):
        ax.plot([center[0],data[0][iseg]], [center[1],data[1][iseg]], **kwargs)
    
    if flagLegend:
        ax.legend()
        
    return ax

def fault(faults, color='black', flagLegend=False, label="segments",
          ax=None, **kwargs):
    '''
    Plotting a Fault system.
    **kwargs : arguments passed to matplotlib.pyplot.plot
    '''
    color = kwargs.setdefault('color', color)
    label = kwargs.setdefault('label', label)
        
    ax = getNewAxes(ax, 1)
    
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

def XY(xtab, ytab, flagAsPoint=False, flagLegend=False, 
       color='blue', marker='o', markersize=10, linestyle='-',
       ax=None, label='data', **kwargs):

    kwargs.setdefault('label', label)
    kwargs.setdefault('color', color)
    if flagAsPoint:
        kwargs.setdefault('markersize', markersize)
        kwargs.setdefault('marker', marker)
    else:
        kwargs.setdefault('linestyle',linestyle)

    if not len(ytab) == len(xtab):
        print("Arrays 'xtab' and 'ytab' should have same dimensions")
        return None;
    
    ax = getNewAxes(ax, 0)
        
    ax.plot(xtab, ytab, **kwargs)
            
    if flagLegend:
        ax.legend()
        
    return ax

def sample(sample, color='black', marker='o', markersize=10,
           flagLegend=False,
           ax=None, label='data', 
           **kwargs):
    
    ax = getNewAxes(ax, 0)
    
    ax.plot(sample[0], sample[1], marker=marker, markersize=markersize, color=color,
            label=label, **kwargs)
            
    if flagLegend:
        ax.legend()
        
    return ax
    
def rule(rule, proportions=[],cmap=None, maxG=3., ax=None):

    ax = getNewAxes(ax, 0)
    
    ax.geometry(xlim=[-maxG,+maxG], ylim=[-maxG,+maxG])    
    nfac = rule.getFaciesNumber()
    rule.setProportions(proportions)
    
    cols = get_cmap(nfac, cmap)

    for ifac in range(nfac):
        bds = rule.getThresh(ifac+1)
        rect = ptc.Rectangle((bds[0],bds[2]),bds[1]-bds[0], bds[3]-bds[2], 
                              color=cols(ifac))
        ax.add_patch(rect)

    return ax


def table(table, icols, fmt='ok', flagLegend=False, ax=None, **kwargs):
    '''
    Plotting the contents of a Table (argument 'table')
        icols designates the ranks of the variable (0: ordinate; 1: abscissae [or regular]) 
        fmt designates [marker][line][color] information
    **kwargs
    '''
    ax = getNewAxes(ax, 0)
    
    if len(icols) == 1:
        datay = table.getColumn(int(icols[0]))
        datax = [i for i in range(table.getRowNumber())]
    else:
        datay = table.getColumn(int(icols[0]))
        datax = table.getColumn(int(icols[1]))
    
    data = np.stack((np.array(datax), np.array(datay)))
    data = data[:, ~np.isnan(data).any(axis=0)]

    ax.plot(data[0,:], data[1,:], **kwargs)
    
    if flagLegend:
        ax.legend()
        
    return ax

def mesh(mesh, 
         flagEdge=True, flagFace=False, flagApex=False, 
         facecolor="yellow", edgecolor="blue", linewidth=1,
         ax=None, **kwargs):
    """
    Plotting the contents of a Mesh
    **kwargs : arguments passed to matplotlib.pyplot.fill
    """
    if flagFace:
        kwargs.setdefault('facecolor', facecolor)
    else:
        kwargs.setdefault('facecolor', "white")
       
    if flagEdge:
        kwargs.setdefault('edgecolor', edgecolor) 
        kwargs.setdefault('linewidth', linewidth)

    ax = getNewAxes(ax, 1)   

    nmesh = mesh.getNMeshes()
    
    for imesh in range(nmesh):
        tabx = mesh.getCoordinatesPerMesh(imesh, 0, True)
        taby = mesh.getCoordinatesPerMesh(imesh, 1, True)
        ax.fill(tabx, taby, **kwargs)
        
        if flagApex:
            ax.scatter(tabx, taby, color='black')

    return ax

def correlation(db, namex, namey, db2=None, usesel=True, 
                asPoint = False,  flagSameAxes=False,
                diagLine=False, diagColor="black", diagLineStyle='-',
                bissLine=False, bissColor="red", bissLineStyle='-',
                regrLine=False, regrColor="blue", regrLineStyle='-',
                ax=None, **kwargs):
    '''
    Plotting the scatter plot between two variables contained in a Db
    
    kwargs: additional arguments used in hist2d or scatter
    '''
    ax = getNewAxes(ax, 0)
        
    if db2 is None:
        db2 = db
   
    if db.getSampleNumber() != db2.getSampleNumber():
        print("Db and Db2 should have the same number of samples")
        return None

    tabx, taby = getBiDefinedValues(db, namex, namey, db2, usesel)
    if len(tabx) == 0:
        return None
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
        regr = gl.regression(db2, namey, [namex], flagCste=True)
        if regr.nvar == 0:
            return None
        a = regr.coeffs[0]
        b = regr.coeffs[1]
        u=[xmin, xmax]
        v=[a+b*xmin, a+b*xmax]
        ax.plot(u,v,color=regrColor,linestyle=regrLineStyle)
        
    ax.decoration(xlabel = db.getName(namex)[0], ylabel = db.getName(namey)[0])

    return ax

def anam(anam, color='blue', linestyle='-', flagLegend=False, ax=None):
    
    res = anam.sample()
    ax = XY(res.getY(), res.getZ(),
            flagLegend=flagLegend, color=color, linestyle=linestyle,
            label='Anamorphosis', ax=ax)
    ax.geometry(xlim = res.getAylim(), ylim=res.getAzlim())
    
    return ax

def plot(object, name1=None, name2=None, ranks=None, **kwargs):
    filetype = type(object).__name__

    if filetype == "Db":
        if name1 is None:
            name1 = object.getLastName()
        flagDb = True
        if name2 is not None:
            flagDb = False
        if flagDb:
            point(object, name1, **kwargs)
        else:
            correlation(object, name1, name2, **kwargs)
            
    elif filetype == "DbGrid":
        if name1 is None:
            name1 = object.getLastName()
        grid(object, name1, **kwargs)
    
    elif filetype == "Vario":
        vario(object, **kwargs)
    
    elif filetype == "Model":
        model(object, **kwargs)
    
    elif filetype == "Rule":
        rule(object, **kwargs)
    
    elif filetype == "Table":
        table(object,ranks, **kwargs)

    elif filetype == "Polygon":
        polygon(object,colorPerSet=True,flagFace=True, **kwargs)
        
    else:
        print("Unknown type")

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

    elif filetype == "Polygon":
        polygon_item = gl.Polygons.createFromNF(filename,False)
        plot(polygon_item, **kwargs)
        
    else:
        print("Unknown type")


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
    def __init__(self, ax=None, collection=None, mydb=None, pickradius=7, color='r', verbose=False):
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

setattr(gl.Db,"plot", gp.point)
setattr(gl.Db,"plot_correlation", gp.correlation)
setattr(gl.Db,"plot_hist", gp.hist)

setattr(gl.DbGrid,"plot", gp.grid)
setattr(gl.DbGrid,"plot_point", gp.point)

# plot_correlation and plot_hist are already inherited from the parent class Db

setattr(gl.Vario,            "plot", gp.vario)
setattr(gl.Model,            "plot", gp.model)
setattr(gl.Rule,             "plot", gp.rule)
setattr(gl.Table,            "plot", gp.table)
setattr(gl.Faults,           "plot", gp.fault)
setattr(gl.Polygons,         "plot", gp.polygon)
setattr(gl.AnamHermite,      "plot", gp.anam)
setattr(gl.MeshEStandardExt, "plot", gp.mesh)
setattr(gl.MeshETurbo,       "plot", gp.mesh)

setattr(plt.Axes, "decoration",    gp.decoration)
setattr(plt.Axes, "geometry",      gp.geometry)

setattr(plt.Axes, "pointSymbol",   gp.pointSymbol)
setattr(plt.Axes, "pointLabel",    gp.pointLabel)
setattr(plt.Axes, "pointGradient", gp.pointGradient)
setattr(plt.Axes, "pointTangent",  gp.pointTangent)

setattr(plt.Axes, "gridRaster",    gp.gridRaster)
setattr(plt.Axes, "gridContour",   gp.gridContour)

setattr(plt.Axes, "varioElem",     gp.varioElem)
setattr(plt.Axes, "modelElem",     gp.modelElem)