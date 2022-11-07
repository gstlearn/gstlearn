import matplotlib.pyplot     as plt
import matplotlib.patches    as ptc
import matplotlib.transforms as transform
import matplotlib.colors     as mcolors
import numpy                 as np
import numpy.ma              as ma
import gstlearn              as gl
from mpl_toolkits.axes_grid1 import make_axes_locatable
from numpy import shape

def get_cmap(n, name='gist_rainbow'):
    '''Returns a function that maps each index in 0, 1, ..., n-1 to a distinct 
        RGB color; the keyword argument name must be a standard mpl colormap name.'''
    return plt.cm.get_cmap(name, n)
    
def selectItems(nvalues, sitem=-1):
    outs = range(0, nvalues)
    nout = nvalues
    if sitem >= 0:
        outs = range(sitem, sitem+1)
        nout = 1
    return outs, nout

def newFigure(figsize = (8,8), xlim = None, ylim = None, nx=1, ny=1, ylimnodiag = None,
              sharex=False, sharey=False):
    ''' Creates a new figure (possibly containing multiple subplots)
    
        Parameters
        ----------
        figsize:    Vector of dimensions along X and Y (per subplot), default is (8,8).
        xlim, ylim: Limits along X and Y (this applies to all subplots of the figure)
        ylimnodiag: Same as ylim for non-diagonal subplot
        nx, ny:     Number of subplots along X and Y
        sharex, sharey: if the subplots should all share respectively X and Y axis (default are False)
        
        Returns
        -------
        Tuple composed of a figure and ax description'''
        
    if figsize is not None:
        figsize = [figsize[0]*nx, figsize[1]*ny]
        
    fig, ax = plt.subplots(nx, ny, figsize=figsize, squeeze=False, sharex=sharex, sharey=sharey)
    
    if xlim is not None:
        for ix in range(nx):
            for iy in range(ny):
                ax[ix,iy].set_xlim(xlim)

    if ylim is not None:
        for ix in range(nx):
            for iy in range(ny):
                if ix != iy:
                    if ylimnodiag is not None:
                        ax[ix,iy].set_ylim(ylimnodiag)
                else:
                    ax[ix,iy].set_ylim(ylim)
        
    if nx * ny == 1:
        ax = ax[0,0]
    return fig, ax

def shape_Nsubplots(N):
    """
    Calculates numbers of lines and columns for N subplots. 
    If N is a perfect square, then nlines = ncols = sqrt(N). 
    Else, the number of columns increase first.
    """
    nlines = np.floor(np.sqrt(N))
    ncols = np.ceil(N/nlines)
    return int(nlines), int(ncols)

def drawDecor(ax=None, xlabel=None, ylabel=None, aspect=None, title=None, flagLegend=False):
    if ax is None:
        if xlabel is not None:
            plt.xlabel(xlabel)
        if ylabel is not None:
            plt.ylabel(ylabel)
        if title is not None:
            plt.title(title)
    else:
        if xlabel is not None:
            ax.set_xlabel(xlabel)
        if ylabel is not None:
            ax.set_ylabel(ylabel)
        if title is not None:
            ax.set_title(title)
        if aspect is not None:
            ax.set_aspect(aspect)
        if flagLegend:
            ax.legend()

    return

def addColorbar(im, ax):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, ax=ax, cax=cax)
    return cbar

def getDefinedValues(db, name, posx=0, posy=1, corner=None, usesel=True, 
                     compress=False, asGrid=True, flagConvertNanToZero=False):

    if db.isGrid() and asGrid:
        if db.getNDim() == 2:
            posx = 0
            posy = 1
            
        if corner is None:
            corner = np.zeros(db.getNDim())
        
        if db.getNDim() == 1:
            tabx = db.getColumn(name, usesel)
        else:
            tabx = db.getOneSlice(name, posx, posy, corner, usesel)
    else:
        tabx = db.getColumn(name, usesel)
    tabx = np.array(tabx).transpose()

    if flagConvertNanToZero:
        tabx[np.isnan(tabx)] = 0
    else:
        tabx = ma.array(tabx,mask=np.isnan(tabx))
    
    if compress:
        tabx = tabx[np.logical_not(np.isnan(tabx))]
        
    return tabx

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

def update_xylim(ax, xlim=None, ylim=None):
    """Update x and y limits by keeping the maximum extent between initial limits and input."""
    if xlim is not None:
        xlim0 = ax.get_xlim()
        if xlim[0] is not None and xlim[0] < xlim0[0]:
            ax.set_xlim(left=xlim[0])
        if xlim[1] is not None and xlim[1] > xlim0[1]:
            ax.set_xlim(right=xlim[1])
    if ylim is not None:
        ylim0 = ax.get_ylim()
        if ylim[0] is not None and ylim[0] < ylim0[0]:
            ax.set_ylim(bottom=ylim[0])
        if ylim[1] is not None and ylim[1] > ylim0[1]:
            ax.set_ylim(top=ylim[1])

def varioElem(vario, ivar=0, jvar=0, idir=0, color0='black', 
              linestyle='solid', linestyle0='dashed', hmax=None, gmax=None, show_pairs = False,
              flagLabelDir=False, flagLegend=False, flagLabelSill=False,
              title=None, xlabel=None, ylabel=None, label=None,
              ax=None, figsize=None, end_plot = False, flagDrawVariance = True,
              **plot_args):
    """Plot a single and unidirectional experimental variogram (one direction and fixed variable(s)).
    
    Parameters
    ----------
    vario : experimental variogram to be represented (gstlearn.Vario).
    ivar, jvar : Indices of the variables for the variogram to be represented (the default is 0).
    idir : Index of the direction of the variogram to be represented (the default is 0).
    color0 : color of the horizontal line representing the sill (the default is 'black').
    linestyle0 : linestyle of the horizontal line representing the sill (the default is 'dashed').
    hmax : Maximum distance to be represented.
    gmax : Maximum variogram value to be represented.
    show_pairs : Flag for annotating the number of pairs for each lag on the plot (the default is False).
    label : Label to be drawn if flagLegend is switched ON
    flagLabelDir : Flag to add the direction vector of the variogram in the label of the line. 
    flagLegend : Flag to display the axes legend (The default is False).
    flagLabelSill : Flag to define the label for the sill (The default is True).
    title : Optional title for the axes.
    ax : Reference for the plot within the figure. If None (default), it creates a new figure.
    figsize : (if ax is None) Sizes (width, height) of figure (in inches).
    end_plot : Flag for closing the graphics (the default is False).
    flagDrawVariance : Flag to add the variance (default is True)
    **plot_args : arguments passed to matplotlib.pyplot.plot

    Returns
    -------
    ax : axes where the variogram is represented

    """
    color = plot_args.setdefault('color', color0)
    linestyle = plot_args.setdefault('linestyle', linestyle)    
    
    if ax is None:
        fig, ax = newFigure(figsize, None, None)

    if label is None:
        label = "vario"
    if flagLabelDir:
        label = "vario dir={}".format(np.round(vario.getCodir(idir),3))
    
    # Plotting the experimental variogram
    gg = vario.getGgVec(idir,ivar,jvar)
    hh = vario.getHhVec(idir,ivar,jvar)
    if len(hh) == 0:
        return ax
    
    ax.plot(hh, gg, label = label, **plot_args)
    hmax = np.nanmax(hh)
    
    # Plotting relevant control lines
    if ivar != jvar:
        ax.hlines(0, 0, hmax, colors="black", linewidth=0.5)
    labsill = None
    if flagLabelSill:
        labsill = "sill"

    if flagDrawVariance:
        ax.hlines(vario.getVar(ivar,jvar), 0, hmax, color0, linestyle0, label = labsill)
    
    if show_pairs:
        pairs = vario.getSwVec(idir,ivar,jvar)
        for i in range(len(hh)):
            ax.annotate(str(int(pairs[i])), (hh[i],gg[i]), xytext=(0,5), xycoords = 'data',
                        textcoords = 'offset points', ha='center')
    
    drawDecor(ax, xlabel=xlabel, ylabel=ylabel, title=title, flagLegend=flagLegend)
    
    if vario.drawOnlyPositiveX(ivar, jvar):
        ax.set_xlim(left=0)
    if vario.drawOnlyPositiveY(ivar, jvar):
        ax.set_ylim(bottom=0)
        
    if end_plot:
        plt.show()
    
    return ax

def varioDir(vario, ivar=0, jvar=0,
             color0='black', linestyle0='dashed', hmax=None, gmax=None, 
             show_pairs=False, cmap=None, flagLegend=False, title=None, 
             xlabel=None, ylabel=None, label=None, ax=None, figsize=None, 
             end_plot=False, **plot_args):
    """Plot a single directional experimental variogram (all available directions, for fixed variable(s)).
    
    Calls the function varioElem for each direction, and labels are automatically set with direction vectors.
    
    Parameters
    ----------
    vario : experimental variogram to be represented (gstlearn.Vario).
    ivar, jvar : Indices of the variables for the variogram to be represented (the default is 0).
    color0 : color of the horizontal line representing the sill (the default is 'black').
    linestyle0 : linestyle of the horizontal line representing the sill (the default is 'dashed').
    hmax : Maximum distance to be represented.
    gmax : Maximum variogram value to be represented.
    cmap : Optional Color scale
    flagLegend : Flag to display the axes legend (The default is False).
    title : Optional title for the axes.
    ax : Reference for the plot within the figure. If None (default), it creates a new figure.
    figsize : (if ax is None) Sizes (width, height) of figure (in inches).
    end_plot : Flag for closing the graphics (the default is False).

    **plot_args : arguments passed to matplotlib.pyplot.plot for all directions plotted
    
    Returns
    -------
    ax : axes where the variogram is represented

    """
    
    if hmax is None:
        hmax = vario.getHmax(ivar, jvar, -1)
    if gmax is None:
        gmax = vario.getGmax(ivar, jvar, -1, True, True) * 1.1

    ndir = vario.getDirectionNumber()
    ndirUtil, ivarD = selectItems(ndir)
    cols = get_cmap(vario.getDirectionNumber(),cmap)
    
    if ax is None:
        fig, ax = newFigure(figsize, None, None)
    
    for idirUtil in ndirUtil:
        flagLabelSill = idirUtil == 0
        varioElem(vario, ivar=ivar, jvar=jvar, idir=idirUtil, color=cols(idirUtil), 
                  color0=color0, linestyle0=linestyle0, show_pairs=show_pairs,
                  ax=ax, hmax=hmax, gmax=gmax, figsize=figsize, flagLabelDir=True, 
                  flagLabelSill=flagLabelSill, label=label,
                  **plot_args)
        
    drawDecor(ax, xlabel=xlabel, ylabel=ylabel, title=title, flagLegend=flagLegend)
    
    ax.autoscale(True)
    
    if vario.drawOnlyPositiveX(ivar, jvar):
        ax.set_xlim(left=0)
    if vario.drawOnlyPositiveY(ivar, jvar):
        ax.set_ylim(bottom=0)
    
    if end_plot:
        plt.show()
        
    return ax

def varmod(vario, mymodel=None, ivar=-1, jvar=-1, idir=-1,
           linestyle='solid', linestylem="dashed", color0='black', linestyle0="dotted",
           nh = 100, hmax = None, gmax = None, show_pairs=False, asCov=False,
           cmap=None, flagLegend=False, title=None, axs=None, figsize=None, end_plot=False, 
           **plot_args):
    """Plot experimental variogram(s) and model (can be multidirectional and multivariable or selected ones).
    
    Same as vario plus the possible model.
    
    Parameters
    ----------
    vario : experimental variogram to be represented (gstlearn.Vario).
    model : optional, variogram model (gstlearn.Model).
    ivar, jvar : Indices of the variables for the variogram to be represented. If -1 (default), all 
                 variables are selected and all the simple and crossed variograms are represented.
    idir : Index of the direction of the variogram to be represented. If -1 (default) all available
           directions are selected and multidirectional variograms are represented.
    linestylem : linestyle of the model lines (the default is 'dashed').
    color0 : color of the horizontal lines representing the experimental sills (the default is 'black').
    linestyle0 : linestyle of the horizontal lines representing the sills (the default is 'dotted').
    nh : number of points between 0 and hmax where the model variogram is calculated (default is 100).
    hmax : Maximum distance to be represented.
    gmax : Maximum variogram value to be represented.
    cmap : Optional Color scale
    flagLegend : Flag to display the axes legend (The default is False).
    title : Optional title for the figure (suptitle).
    axs : Reference for the plot(s) within the figure. If None (default),
          it creates a new figure (with multiple axes for multivariate variograms).
    figsize : (if ax is None) Sizes (width, height) of figure (in inches).
    end_plot : Flag for closing the graphics (the default is False).

    **plot_args : arguments passed to matplotlib.pyplot.plot for all variograms plotted (not models!)
    
    Returns
    -------
    ax : axes where the variograms are represented

    """
    
    if hmax is None:
        hmax = vario.getHmax(ivar, jvar, idir)
    if gmax is None:
        gmax = vario.getGmax(ivar, jvar, idir)

    ylimnodiag = [-gmax, gmax]
        
    ndir = vario.getDirectionNumber()
    nvar = vario.getVariableNumber()
    cols = get_cmap(ndir,cmap)
    
    ndirUtil, ivarD = selectItems(ndir, idir)
    ivarUtil, ivarN = selectItems(nvar, ivar)
    jvarUtil, jvarN = selectItems(nvar, jvar)
        
    # Create a new figure
    if axs is None:
        fig, axs = newFigure(figsize, None, None, ivarN, jvarN, ylimnodiag)   
        
    # if several directions, label with the direction vectors
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
                varioElem(vario, iv, jv, idirUtil, 
                          color=cols(idirUtil), linestyle=linestyle,
                          color0=color0, linestyle0=linestyle0, show_pairs=show_pairs,
                          ax=ax, hmax=hmax, gmax=None, 
                          flagLabelDir=flagLabelDir, flagLegend=flagLegend, **plot_args)

                # Plotting the Model (optional)
                if mymodel is not None:
                    codir = vario.getCodir(idirUtil)
                    model(mymodel, ivar=iv, jvar=jv, codir=codir, 
                          color=cols(idirUtil), linestyle=linestylem, 
                          color0=color0, linestyle0=linestyle0, ax=ax,
                          hmax=hmax, gmax=None, nh=nh, asCov=asCov,
                          flagLabelDir=flagLabelDir, flagLegend=flagLegend)

            ax.autoscale(True)
            
            if vario.drawOnlyPositiveX(iv, jv):
                ax.set_xlim(left=0)
            if vario.drawOnlyPositiveY(iv, jv):
                ax.set_ylim(bottom=0)
    
    if title is not None:
        plt.suptitle(title)
        
    if end_plot:
        plt.show()
        
    return axs

def vario(vario, ivar=-1, jvar=-1, idir=-1,
          linestyle='solid', color0='black', linestyle0='dashed', hmax=None, gmax=None, 
          cmap = None, flagLegend=False, 
          title = None, axs = None, figsize = None, end_plot=False, **plot_args):
    """Plot experimental variogram(s) (can be multidirectional and multivariable or selected ones).
    
    Parameters
    ----------
    vario : experimental variogram to be represented (gstlearn.Vario).
    ivar, jvar : Indices of the variables for the variogram to be represented. If -1 (default), all 
                 variables are selected and all the simple and crossed variograms are represented.
    idir : Index of the direction of the variogram to be represented. If -1 (default) all available
           directions are selected and multidirectional variograms are represented.
    color0 : color of the horizontal lines representing the sills (the default is 'black').
    linestyle0 : linestyle of the horizontal lines representing the sills (the default is 'dashed').
    hmax : Maximum distance to be represented.
    gmax : Maximum variogram value to be represented.
    cmap : Optional Color scale
    flagLegend : Flag to display the axes legend (The default is False).
    title : Optional title for the figure (suptitle).
    axs : Reference for the plot(s) within the figure. If None (default),
          it creates a new figure (with multiple axes for multivariate variograms).
    figsize : (if ax is None) Sizes (width, height) of figure (in inches).
    end_plot : Flag for closing the graphics (the default is False).

    **plot_args : arguments passed to matplotlib.pyplot.plot for all variograms plotted

    Returns
    -------
    ax : axes where the variograms are represented

    """
    axs = varmod(vario, mymodel=None, ivar=ivar, jvar=jvar, idir=idir, 
                 linestyle=linestyle, color0=color0, linestyle0=linestyle0, 
                 hmax=hmax, gmax=gmax, cmap=cmap, 
                 flagLegend=flagLegend, title=title, axs=axs, figsize=figsize, 
                 end_plot=end_plot, **plot_args)
    
    return axs

def model(model, ivar=0, jvar=0, codir=None, color0='black', linestyle0='dashed',
          nh = 100, flagEnv = True, hmax = None, gmax = None, 
          flagLabelDir=False, flagLegend=False, asCov=False,
          title=None, xlabel=None, ylabel=None, ax=None, 
          figsize = None, end_plot =False, **plot_args):
    """Plot a single and unidirectional variogram model (one direction and fixed variable(s)).
    
    Parameters
    ----------
    model : variogram model to be represented (gstlearn.Model).
    ivar, jvar : Indices of the variables for the variogram to be represented (the default is 0).
    codir : Vector of the direction of the variogram to be represented. The default is the unit 
            vector in the first space dimension.
    color0, linestyle0 : (if flagEnv = True) color and linestyle of the correlation envelop 
                         (the defaults are 'black' and 'dashed')
    nh : number of points between 0 and hmax where the model variogram is calculated (default is 100).
    flagEnv : flag for representing the correlation envelop (the default is True)
    hmax : Maximum distance to be represented. The default is 3 times the maximum range of the
           basic structures, or 1 if no range is defined.
    gmax : Maximum variogram value to be represented.
    flagLabelDir : Flag to add the direction vector codir in the label of the line. 
                   The default label is "model" (default is False).
    flagLegend : Flag to display the axes legend (The default is False).
    asCov : Present the Model as a Covariance (rather than as a Variogram)
    title : Optional title for the axes.
    ax : Reference for the plot within the figure. If None (default), it creates a new figure.
    figsize : (if ax is None) Sizes (width, height) of figure (in inches).
    end_plot : Flag for closing the graphics (the default is False).

    Returns
    -------
    ax : axes where the variogram is represented

    """
    color = plot_args.setdefault('color', color0)    
    
    if codir is None:
        ndim = model.getDimensionNumber()
        codir = [0] * ndim
        codir[0] = 1

    # if hmax not specified = 3*maximum range of the model's basic structures
    if hmax is None:
        hmax = 0
        for icova in range(model.getCovaNumber()):
            range_max = np.max(model.getCova(icova).getRanges())
            if 3*range_max > hmax:
                hmax = 3*range_max
    if hmax == 0: # if the model has no range defined
        hmax = 1
            
    hh = np.linspace(0, hmax, nh+1)
    gg = model.sample(hmax, nh, ivar, jvar, codir, asCov=asCov, addZero=True)
    
    if ax is None:
        fig, ax = newFigure(figsize, None, None)
        
    label = "model"
    if flagLabelDir:
        label = "model dir={}".format(codir)
    
    istart = 0
    for i in range(model.getCovaNumber()):
        if model.getCovName(i) == 'Nugget Effect':
            istart = 1 # do not plot the first lag (h=0) for nugget effect (discontinuity)
            
    ax.plot(hh[istart:], gg[istart:], label=label, **plot_args)
    
    if ivar != jvar and flagEnv:
        ggp = model.sample(hmax, nh, ivar, jvar, codir, 1, asCov=asCov, addZero=True)
        ax.plot(hh[istart:], ggp[istart:], color = color0, linestyle = linestyle0, label="plus")
        ggm = model.sample(hmax, nh, ivar, jvar, codir,-1, asCov=asCov, addZero=True)
        ax.plot(hh[istart:], ggm[istart:], color = color0, linestyle = linestyle0, label="minus")
    
    drawDecor(ax, xlabel=xlabel, ylabel=ylabel, title=title, flagLegend=flagLegend)
    
    if end_plot:
        plt.show()

    return ax

def point(db, 
          color_name=None, size_name=None, elev1D_name=None, label_name=None, usesel=True, 
          color='r', size=20, sizmin=10, sizmax=200, edgecolors=None,
          xlim=None, ylim=None, directColor=False, flagAbsSize=False,
          cmap=None, flagColorBar=True, flagSizeLegend=True, aspect=None,
          posX=0, posY=1, xlabel=None, ylabel=None,
          title=None, ax=None, figsize=None, end_plot=False, **scatter_args):
    '''Function for plotting a point data base, with optional color and size variables
    
    db: Db containing the variable to be plotted
    color_name: Name of the variable containing the color per sample
    size_name: Name of the variable containing the size per sample
    elev1D_name: Name of the variable standing for Y coordinate in 1-D case
    label_name: Name of the variable containing the label per sample
    usesel : Boolean to indicate if the selection has to be considered
    color: Constant color (used if 'color_name' is not defined)
    size: Constant size (used if 'size_name' is not defined)
    sizmin: Size corresponding to the smallest value (used if 'size_name' is defined)
    sizmax: Size corresponding to the largest value (used if 'size_name' is defined)
    xlim: Bounds defined along the first axis
    ylim: Bounds defined along the second axis
    directColor: True if the value of the field directly indicates the Rank in the Color Scale
    flagAbsSize: Represent the Absolute value in Size representation
    cmap: Optional Color scale
    flagColorBar: Flag for representing the Color Bar (not represented if color_name=None)
    flagSizeLegend: Flag for representing the Legend for marker size (not represented if size_name=None)
    aspect: aspect ratio of the axes scaling, i.e. y/x-scale. 
    posX: rank of the first coordinate
    posY: rank of the second coordinate
    title: Title given to the plot
    ax: Reference for the plot within the figure
    figsize: (if ax is None) Sizes (width, height) of figure (in inches)
    end_plot: Flag for closing the graphics

    **scatter_args : arguments passed to matplotllib.pyplot.scatter
    '''
    
    edgecolors = scatter_args.setdefault('edgecolors', edgecolors)
    
    if ax is None:
        fig, ax = newFigure(figsize, xlim, ylim)

    # Extracting coordinates
    tabx = db.getCoordinates(posX,usesel)
    if db.getNDim() > 1:
        taby = db.getCoordinates(posY,usesel)
    else:
        taby = db.getColumn(elev1D_name, usesel)
    if len(tabx) <= 0 or len(taby) <= 0:
        return
    
    # Color of symbol
    if color_name is not None:
        colval = getDefinedValues(db, color_name, 0, 1, None, usesel, 
                                  compress=False, asGrid=False, flagConvertNanToZero=True)
        if (directColor and cmap is not None):
            colval = colval.astype(int)
            colval = cmap(colval)
    else:
        colval = color

    # Size of symbol
    if size_name is not None:
        sizval = getDefinedValues(db, size_name, 0, 1, None, usesel, 
                                  compress=False, asGrid=False, flagConvertNanToZero=True)
        
        if flagAbsSize:
            sizval = np.absolute(sizval)
        m = np.nanmin(np.absolute(sizval))
        M = np.nanmax(np.absolute(sizval))
        sizval = (sizmax - sizmin) * (np.absolute(sizval) - m) / (M-m) + sizmin
    else:
        sizval = size

    im = ax.scatter(x = tabx, y = taby, s = sizval, c = colval, cmap=cmap, **scatter_args)

    if label_name is not None:
        labval = getDefinedValues(db, label_name, 0, 1, None, usesel, 
                                  compress=False, asGrid=False, flagConvertNanToZero=True)
        for i in range(len(labval)):
            ax.text(tabx[i], taby[i], round(labval[i],2))
            
    drawDecor(ax, xlabel=xlabel, ylabel=ylabel, title=title, aspect=aspect)
        
    if flagColorBar and (color_name is not None):
        addColorbar(im, ax)
    
    if flagSizeLegend and (size_name is not None):
        labels = lambda marker_size : (M - m)*(marker_size - sizmin)/(sizmax - sizmin) + m
        # labels is the inverse transformation from marker sizes to variable values
        ax.legend(*im.legend_elements("sizes", num=5, color=cmap, func=labels))
         
    if end_plot:
        plt.show()

    return ax

def gradient(db, elev1D_name=None, usesel=True, color='black', scale=20, 
             xlim=None, ylim=None, aspect=None,
             title=None, ax=None, figsize=None, end_plot=False):
    '''Function for plotting a gradient data base
    
    db: Db containing the variable to be plotted
    usesel : Boolean to indicate if the selection has to be considered
    color: Constant color 
    scale: Constant scale
    xlim: Bounds defined along the first axis
    ylim: Bounds defined along the second axis
    aspect: aspect ratio of the axes scaling, i.e. y/x-scale. 
    title: Title given to the plot
    ax: Reference for the plot within the figure
    figsize: (if ax is None) Sizes (width, height) of figure (in inches)
    end_plot: Flag for closing the graphics
    '''
    
    if ax is None:
        fig, ax = newFigure(figsize, xlim, ylim)

    # Extracting coordinates
    tabx  = db.getCoordinates(0,usesel)
    if db.getNDim() > 1:
        taby  = db.getCoordinates(1,usesel)
    else:
        taby = db.getColumn(elev1D_name, usesel)

    if db.getNDim() > 1:
        tabgx = db.getGradients(0,usesel)
        tabgy = db.getGradients(1,usesel)
    else:
        tabgy = -db.getGradients(0,usesel)
        tabgx = -np.ones(len(tabgy))

    if len(tabx) <= 0 or len(taby) <= 0 or len(tabgx) <= 0 or len(tabgy) <= 0:
        return
    
    ax.quiver(tabx, taby, -tabgx, -tabgy, angles='xy', color=color, scale=scale)
            
    drawDecor(ax, title=title, aspect=aspect)
        
    if end_plot:
        plt.show()

    return ax


def tangent(db, usesel=True, color='black', scale=20, 
            xlim=None, ylim=None, aspect=None,
            title=None, ax=None, figsize=None, end_plot=False):
    '''Function for plotting a tangent data base
    
    db: Db containing the variable to be plotted
    usesel : Boolean to indicate if the selection has to be considered
    color: Constant color 
    scale: Constant scale
    xlim: Bounds defined along the first axis
    ylim: Bounds defined along the second axis
    aspect: aspect ratio of the axes scaling, i.e. y/x-scale. 
    title: Title given to the plot
    ax: Reference for the plot within the figure
    figsize: (if ax is None) Sizes (width, height) of figure (in inches)
    end_plot: Flag for closing the graphics
    '''
    
    if ax is None:
        fig, ax = newFigure(figsize, xlim, ylim)

    # Extracting coordinates
    tabx  = db.getCoordinates(0,usesel)
    taby  = db.getCoordinates(1,usesel)
    tabtx = db.getTangents(0,usesel)
    tabty = db.getTangents(1,usesel)

    if len(tabx) <= 0 or len(taby) <= 0 or len(tabtx) <= 0 or len(tabty) <= 0:
        return
    
    ax.quiver(tabx, taby, -tabtx, -tabty, color=color, scale=scale)
    ax.quiver(tabx, taby,  tabtx,  tabty, color=color, scale=scale)
            
    drawDecor(ax, title=title, aspect=aspect)
        
    if end_plot:
        plt.show()

    return ax

def polygon(poly, faceColor='yellow', edgeColor = 'blue', aspect=None,
            colorPerSet = False, flagEdge=True, flagFace=False, linewidth=2,
            title= None, ax=None, figsize=None, end_plot=False, **fill_args):
    '''Function to display a polygon
    
    **fill_args: arguments passed to ax.fill
    '''
    
    if ax is None:
        fig, ax = newFigure(figsize)
    
    npol = poly.getPolySetNumber()
    cols = get_cmap(npol)
    
    for ipol in range(npol):
        x = poly.getX(ipol)
        y = poly.getY(ipol)
        
        faceColor_local = 'none'
        if flagFace:
            faceColor_local = faceColor
            if colorPerSet:
                faceColor_local = cols(ipol)
                
        edgeColor_local = 'none'
        if flagEdge:
            edgeColor_local = edgeColor
            if colorPerSet:
                edgeColor_local = cols(ipol)

        ax.fill(x, y, facecolor=faceColor_local, edgecolor=edgeColor_local,
                linewidth=linewidth, **fill_args)
        
    drawDecor(ax, title=title, aspect=aspect)
        
    if end_plot:
        plt.show()
        
    return ax
        
def grid(dbgrid, name = None, usesel = True, flagColorBar=True, aspect=None,
         xlim=None, ylim=None, posx=0, posy=1, corner=None, 
         levels=None, colorL='black', linestyleL = 'solid', 
         flagRaster=True,
         title = None, ax=None, figsize = None, end_plot=False, 
         **plot_args):
    '''
    Function for plotting a variable (referred by its name) informed in a DbGrid

    dbgrid: Db, organized as a Grid, containing the variable to be plotted
    name: Name of the variable to be represented (by default, the first Z locator, or the last field)
    usesel : Boolean to indicate if the selection has to be considered
    flagColorBar: Flag for representing the Color Bar (not represented if alpha=0)
    aspect: aspect ratio of the axes scaling, i.e. y/x-scale.
    xlim: Bounds defined along the first axis
    ylim: Bounds defined along the second axis
    title: Title given to the plot
    ax: Reference for the plot within the figure
    figsize: (if ax is None) Sizes (width, height) of figure (in inches)
    end_plot: Flag for closing the graphics
    
    **plot_args : arguments passed to matplotlib.pyplot.pcolormesh
    '''
    clip_on = plot_args.setdefault('clip_on', True)
    shading = plot_args.setdefault('shading', 'nearest')
    
    if not(dbgrid.isGrid()):
        print("This function is dedicated to Grid Db and cannot be used here")
        return None;
    
    if name is None:
        if dbgrid.getVariableNumber() > 0:
            name = dbgrid.getNameByLocator(gl.ELoc.Z,0) # select locator z1, prints an error if no Z locator
        else : # if no Z locator, choose the last field
            name = dbgrid.getLastName()
    
    if ax is None:
        fig, ax = newFigure(figsize, xlim, ylim)
        
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
    trans_data = tr + ax.transData
    
    if shading == "nearest":
        X = np.linspace(x0, x0 + (nx-1)*dx, nx)
        Y = np.linspace(y0, y0 + (ny-1)*dy, ny)
    elif shading == "flat":
        X = np.linspace(x0, x0 + nx*dx, nx+1)
        Y = np.linspace(y0, y0 + ny*dy, ny+1)
    else:
        print("The argument shading should be either 'nearest' for cells centered on (x,y)"
              " or 'flat' for cells with low-left corner in (x,y)")

    if flagRaster:
        im = ax.pcolormesh(X, Y, data, **plot_args)
        im.set_transform(trans_data)
        
        if flagColorBar:
            addColorbar(im, ax)
    
    if levels is not None:
        ax.contour(X, Y, data, levels, colors=colorL, linestyles=linestyleL)
        
    x1, x2, y1, y2 = x0, X[-1], y0, Y[-1]
    ax.plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], marker='', linestyle='', 
            transform=trans_data)
    
    update_xylim(ax, xlim=xlim, ylim=ylim) 
    
    if title is None:
        title = dbgrid.getNames(name)[0]
        
    drawDecor(ax, title=title, aspect=aspect)
    
    if end_plot:
        plt.show()
    
    return ax

def grid1D(dbgrid, name = None, usesel = True, flagColorBar=True, aspect=None,
         xlim=None, ylim=None,
         color='black',flagLegend=False, label='curve',
         title = None, ax=None, figsize = None, end_plot=False, 
         **plot_args):
    '''
    Function for plotting a variable (referred by its name) informed in a DbGrid

    dbgrid: Db, organized as a Grid, containing the variable to be plotted
    name: Name of the variable to be represented (by default, the first Z locator, or the last field)
    usesel : Boolean to indicate if the selection has to be considered
    flagColorBar: Flag for representing the Color Bar (not represented if alpha=0)
    aspect: aspect ratio of the axes scaling, i.e. y/x-scale.
    xlim: Bounds defined along the first axis
    ylim: Bounds defined along the second axis
    title: Title given to the plot
    ax: Reference for the plot within the figure
    figsize: (if ax is None) Sizes (width, height) of figure (in inches)
    end_plot: Flag for closing the graphics
    
    **plot_args : arguments passed to matplotlib.pyplot.pcolormesh
    '''
    if dbgrid.getNDim() != 1:
        print("This function is dedicated to 1-D Grid")
        return None
    
    if not(dbgrid.isGrid()):
        print("This function is dedicated to Grid Db and cannot be used here")
        return None
    
    if name is None:
        if dbgrid.getVariableNumber() > 0:
            name = dbgrid.getNameByLocator(gl.ELoc.Z,0) # select locator z1, prints an error if no Z locator
        else : # if no Z locator, choose the last field
            name = dbgrid.getLastName()
    
    if ax is None:
        fig, ax = newFigure(figsize, xlim, ylim)
        
    x0 = dbgrid.getX0(0)
    nx = dbgrid.getNX(0)
    dx = dbgrid.getDX(0)
    
    tabx = dbgrid.getColumnByLocator(gl.ELoc.X, 0, usesel)
    data = getDefinedValues(dbgrid, name, 0, 1, None, usesel, 
                            compress=False, asGrid=True)

    curve(data1=tabx, data2=data, color=color, flagLegend=flagLegend, 
          label=label,
          title=title, ax=ax, figsize = figsize, end_plot=end_plot,
          **plot_args)

    return ax

def hist_tab(val, xlabel=None, ylabel=None, nbins=30, color='yellow', edgecolor='red',
             title = None, ax = None, figsize=None, end_plot=False, **hist_args):
    '''Function for plotting the histogram of an array (argument 'val')
    
    hist_args : arguments passed to matplotlib.pyplot.hist
    '''
    color     = hist_args.setdefault("color", color)
    edgecolor = hist_args.setdefault("edgecolor", edgecolor)
    bins      = hist_args.setdefault("bins", nbins)
    
    if ax is None:
        fig, ax = newFigure(figsize)
        
    ax.hist(val, **hist_args)
    
    drawDecor(ax, xlabel=xlabel, ylabel=ylabel, title=title)
        
    if end_plot:
        plt.show()
        
    return ax
    
def hist(db, name, xlabel=None, ylabel=None, title = None, ax=None,
         figsize=None, end_plot=False, usesel=True, **hist_args):
    '''Function for plotting the histogram of a variable contained in a Db
    
    hist_args : arguments passed to matplotlib.pyplot.hist'''
    
    db.useSel = usesel
    val = db[name]
    if len(val) == 0:
        return None
    
    if title is None:
        title = db.getNames(name)[0]
    ax = hist_tab(val, title=title, xlabel=xlabel, ylabel=ylabel, ax=ax, figsize=figsize, end_plot=end_plot, **hist_args)
    
    return ax

def sortedcurve(tabx, taby, color='black', flagLegend=False,
                label='curve', xlabel=None, ylabel=None,
                title=None, ax=None, figsize=None, end_plot=False, 
                **plot_args):
    '''
    Function for plotting a set of points after they have been sorted in increasing X
    '''
        # Account for possible 'nan'  values
    mask = np.logical_and(np.isfinite(tabx), np.isfinite(taby))
    stabx = tabx[mask]
    staby = taby[mask]
    
    # Indices of the sorted elements of stabx
    indices = np.argsort(stabx)
    ax = curve(stabx[indices], staby[indices], color=color, 
          flagLegend=flagLegend, label=label, xlabel=xlabel, ylabel=ylabel,
          title=title, ax=ax, figsize=figsize, end_plot=end_plot)
    
    return ax
    
def curve(data1, data2=None, icas=1, color='black',flagLegend=False, 
          label='curve', xlabel=None, ylabel=None, 
          title=None, ax=None, figsize = None, end_plot=False, **plot_args):
    '''
    Function for plotting the curve of an array (argument 'data1')
        if data1 is a tuple, it should contain x=data1[0] and y=data1[1]
        or
        'data1' and 'data2' are provided
        otherwise:
        icas=1 when 'data1' contains the abscissa and ordinates are regular
        icas=2 when 'data1' contains the ordinate and abscissa are regular
    **plot_args : arguments passed to matplotlib.pyplot.plot
    '''
    color = plot_args.setdefault('color', color)
    label = plot_args.setdefault('label', label)
    
    if len(data1) == 0:
        return None

    if ax is None:
        fig, ax = newFigure(figsize)
    
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
    
    ax.plot(tabx, taby, **plot_args)
    
    drawDecor(ax, xlabel=xlabel, ylabel=ylabel, title=title, 
              flagLegend=flagLegend)
    
    if end_plot:
        plt.show()
        
    return ax

def multisegments(center, data, color='black',flagLegend=False, label="segments",
                  title=None, ax=None, figsize = None, end_plot=False, **plot_args):
    '''
    Function for plotting a set of segments joining 'center' to any of vertices
    stored in 'data'.
    **plot_args : arguments passed to matplotlib.pyplot.plot
    '''
    color = plot_args.setdefault('color', color)
    label = plot_args.setdefault('label', label)
        
    if len(data) == 0:
        return None
    
    if ax is None:
        fig, ax = newFigure(figsize)
    
    nseg = len(data[0])
    
    for iseg in range(nseg):
        ax.plot([center[0],data[0][iseg]], [center[1],data[1][iseg]], **plot_args)
    
    drawDecor(ax, title=title, flagLegend=flagLegend)
    
    if end_plot:
        plt.show()
        
    return ax

def fault(faults, color='black',flagLegend=False, label="segments",
          title=None, ax=None, figsize = None, end_plot=False, **plot_args):
    '''
    Function for plotting a Fault system.
    **plot_args : arguments passed to matplotlib.pyplot.plot
    '''
    color = plot_args.setdefault('color', color)
    label = plot_args.setdefault('label', label)
        
    if ax is None:
        fig, ax = newFigure(figsize)
    
    nfaults = faults.getNFaults()
    for ifault in range(nfaults):

        fault = faults.getFault(ifault);
        npoints = fault.getNPoints() - 1
        xtab = fault.getX()
        ytab = fault.getY()
        ax.plot(xtab, ytab, **plot_args)
    
    drawDecor(ax, title=title, flagLegend=flagLegend)
    
    if end_plot:
        plt.show()
        
    return ax

def XY(xtab, ytab, flagAsPoint=False, xlim=None, ylim=None, flagLegend=False, 
       color='blue', marker='o', markersize=10, linestyle='-',
       title=None, ax=None, figsize = None, label='data', end_plot=False, **plot_args):
    
    plot_args.setdefault('label', label)
    plot_args.setdefault('color', color)
    if flagAsPoint:
        plot_args.setdefault('markersize', markersize)
        plot_args.setdefault('marker', marker)
    else:
        plot_args.setdefault('linestyle',linestyle)

    if not len(ytab) == len(xtab):
        print("Arrays 'xtab' and 'ytab' should have same dimensions")
        return None;
    
    if ax is None:
        fig, ax = newFigure(figsize, xlim, ylim)
        
    ax.plot(xtab, ytab, **plot_args)
            
    drawDecor(ax, title=title, flagLegend=flagLegend)
    
    if end_plot:
        plt.show()
    
    return ax

def sample(sample, xlim=None, ylim=None, aspect=None,
           color='black', marker='o', markersize=10,
           flagLegend=False,
           title=None, ax=None, figsize = None, label='data', end_plot=False, **plot_args):
    
    if ax is None:
        fig, ax = newFigure(figsize, xlim, ylim)
    
    ax.plot(sample[0], sample[1], marker=marker, markersize=markersize, color=color,
            label=label, **plot_args)
            
    drawDecor(ax, title=title, aspect=aspect, flagLegend=flagLegend)
    
    if end_plot:
        plt.show()

    return ax
    
def rule(rule, proportions=[],cmap=None, 
         title=None, xlim=[-5,+5], ylim=[-5,+5], ax=None, figsize=None, end_plot=False):
    
    if ax is None:
        fig, ax = newFigure(figsize, xlim=xlim, ylim=ylim)
        
    nfac = rule.getFaciesNumber()
    rule.setProportions(proportions)
    
    cols = get_cmap(nfac, cmap)

    for ifac in range(nfac):
        bds = rule.getThresh(ifac+1)
        rect = ptc.Rectangle((bds[0],bds[2]),bds[1]-bds[0], bds[3]-bds[2], 
                              color=cols(ifac))
        ax.add_patch(rect)

    drawDecor(ax, title=title)
        
    if end_plot:
        plt.show()

    return ax


def table(table, icols, fmt='ok', xlim=None, ylim=None, flagLegend=False,
          color0='b', linestyle0='-', marker0='', label='table',
          title=None, ax=None, figsize=None, end_plot=False, **plot_args):
    '''
    Function for plotting the contents of a Table (argument 'tablr')
        icols designates the ranks of the variable (0: ordinate; 1: abscissae [or regular]) 
        fmt designates [marker][line][color] information
    **plot_args
    '''
    plot_args.setdefault('label', label)
    
    if ax is None:
        fig, ax = newFigure(figsize, xlim, ylim)
    
    if len(icols) == 1:
        datay = table.getColumn(int(icols[0]))
        datax = [i for i in range(table.getRowNumber())]
    else:
        datay = table.getColumn(int(icols[0]))
        datax = table.getColumn(int(icols[1]))
    
    data = np.stack((np.array(datax), np.array(datay)))
    data = data[:, ~np.isnan(data).any(axis=0)]

    ax.plot(data[0,:], data[1,:], color=color0, linestyle=linestyle0, marker=marker0, 
            **plot_args)
    
    drawDecor(ax, title=title, flagLegend=flagLegend)
    
    if end_plot:
        plt.show()
        
    return ax

def mesh(mesh, 
         flagEdge=True, flagFace=False, flagApex=False, aspect=None,
         xlim=None, ylim=None, facecolor="yellow", edgecolor="blue", linewidth=1,
         title=None, ax=None, figsize = None, end_plot =False, **plot_args):
    """
    **plot_args : arguments passed to matplotlib.pyplot.fill
    """
    if flagFace:
        plot_args.setdefault('facecolor', facecolor)
    else:
        plot_args.setdefault('facecolor', "white")
       
    if flagEdge:
        plot_args.setdefault('edgecolor', edgecolor) 
        plot_args.setdefault('linewidth', linewidth)

    if ax is None:
        fig, ax = newFigure(figsize, xlim, ylim)   

    nmesh = mesh.getNMeshes()
            
    for imesh in range(nmesh):
        tabx = mesh.getCoordinatesPerMesh(imesh, 0, True)
        taby = mesh.getCoordinatesPerMesh(imesh, 1, True)
        ax.fill(tabx, taby, **plot_args)
        if flagApex:
            ax.scatter(tabx, taby, color='black')

    drawDecor(ax, title=title, aspect=aspect)
    
    if end_plot:
        plt.show()

    return ax

def correlation(db, namex, namey, db2=None, bins=50, xlim=None, ylim=None, usesel=True, 
                asPoint = False, flagAxisLabel = True,
                diagLine=False, diagColor="black", diagLineStyle='-',
                bissLine=False, bissColor="red", bissLineStyle='-',
                regrLine=False, regrColor="blue", regrLineStyle='-',
                xlabel=None, ylabel=None, aspect=None, 
                title = None, ax=None, figsize=None, end_plot=False):
    '''Function for plotting the scatter plot between two variables contained in a Db'''
 
    if ax is None:
        fig, ax = newFigure(figsize, xlim, ylim)
        
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
    
    if asPoint:
        ax.scatter(tabx, taby)
    else:
        ax.hist2d(tabx, taby, bins, cmin=1)

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
        
    if flagAxisLabel:
        if xlabel is None:
            xlabel = db.getNames(namex)[0]
        if ylabel is None:
            ylabel = db.getNames(namey)[0]

    drawDecor(ax, xlabel=xlabel, ylabel=ylabel, title=title, aspect=aspect)
    
    if end_plot:
        plt.show()
        
    return ax

def anam(anam, xlim=None, ylim=None, 
         color='blue', linestyle='-', flagLegend=False,
         xlabel=None, ylabel=None, title = None, ax=None, 
         figsize=None, end_plot=False):
    
    res = anam.sample()
    ax = XY(res.getY(), res.getZ(),
            xlim=res.getAylim(), ylim=res.getAzlim(),
            flagLegend=flagLegend, color=color, linestyle=linestyle,
            label='Anamorphosis', title=title,
            ax=ax, figsize=figsize)
    
    if end_plot:
        plt.show()

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
            point(object, name1, end_plot=True, **kwargs)
        else:
            correlation(object, name1, name2, end_plot=True, **kwargs)
            
    elif filetype == "DbGrid":
        if name1 is None:
            name1 = object.getLastName()
        grid(object, name1, end_plot=True, **kwargs)
    
    elif filetype == "Vario":
        vario(object,end_plot=True, **kwargs)
    
    elif filetype == "Model":
        model(object,end_plot=True, **kwargs)
    
    elif filetype == "Rule":
        rule(object,end_plot=True, **kwargs)
    
    elif filetype == "Table":
        table(object,ranks,end_plot=True, **kwargs)

    elif filetype == "Polygon":
        polygon(object,colorPerSet=True,flagFace=True,end_plot=True, **kwargs)
        
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


### Plot several variables at once

def grids(dbgrid, names = None, usesel = True, flagColorBar=True, aspect=None,
         xlim=None, ylim=None, norm=None,
         title = None, axs=None, figsize = None, end_plot=False, **plot_args):
    '''
    Function for plotting several variables (referred by their names) informed in a grid Db in subplots (Nsubplots=Nvariables)

    dbgrid: Db, organized as a Grid, containing the variable to be plotted
    names: Name of the variables to be represented (by default all the Z locators, or the last field)
    usesel : Boolean to indicate if the selection has to be considered
    flagColorBar: Flag for representing the Color Bar (not represented if alpha=0)
    aspect: aspect ratio of the axes scaling, i.e. y/x-scale.
    norm: Optional norm. It can be either an instance of matplotlib.colors.Normalize when using a unique norm for 
          all variables (e.g. matplotlib.colors.LogNorm()), or a class of matplotlib.colors.Normalize when using 
          a type of normalization scaled independently for each variable (e.g. matplotlib.colors.LogNorm). 
    xlim: Bounds defined along the first axis
    ylim: Bounds defined along the second axis
    title: Title given to the figure (each subplot is titled with the name of the variable represented)
    axs: References for the subplots within the figure: list or array (1D or 2D) containing at least Nvar Axes.
    figsize: (if ax is None) Sizes (width, height) of figure (in inches)
    end_plot: Flag for closing the graphics
    
    **plot_args : arguments passed to matplotlib.pyplot.pcolormesh for every grid plots
    '''
    
    if names is None:
        names = dbgrid.getNamesByLocator(gl.ELoc.Z) # all Z locators
        if names == () : # if no Z locator, choose the last field
            names = dbgrid.getLastName()
    else:
        names = dbgrid.getNames(names)
    names = np.atleast_1d(names)
    
    Nplots = len(names)
    if Nplots == 1:
        ax = grid(dbgrid, names[0], usesel=usesel, flagColorBar=flagColorBar, 
                  aspect=aspect, xlim=xlim, ylim=ylim, 
                  title=title, ax=axs, figsize=figsize, end_plot=end_plot, **plot_args)
        return ax
    elif Nplots == 0:
        print("There is no variable in the dbgrid corresponding to the name given.")
        return None
    
    if axs is None:
        nlines, ncols = shape_Nsubplots(Nplots)
        fig, axs = newFigure(figsize, xlim, ylim, nx=nlines, ny=ncols, sharex=True, sharey=True)
    
    if title is not None:
        fig.suptitle(title)
    
    for i,ax in enumerate(axs.flat):
        if i >= Nplots:
            ax.set_visible(False)
            continue;
        
        norm_i = None
        if norm is not None:
            if isinstance(norm, type): #independent norms for each subplot
                norm_i = norm()
            else: # a unique norm for all subplots
                norm_i = norm
            
        grid(dbgrid, names[i], ax=ax, usesel=usesel, flagColorBar=flagColorBar, 
             aspect=aspect, xlim=xlim, ylim=ylim, norm=norm_i, **plot_args)
    if end_plot:
        plt.show()
        
    return axs

def color_plots(db, names = None, usesel = True, flagColorBar=True, aspect='auto',
         xlim=None, ylim=None, size=20, cmap=None,
         title = None, axs=None, figsize = None, end_plot=False, **plot_args):
    '''
    Function for plotting several variables (referred by their names) informed in a Db in subplots (Nsubplots=Nvariables).
    Variables are represented with color plots, i.e. each data point is represented with a color corresponding to the data value.

    db: Db containing the variable to be plotted
    names: Name of the variables to be represented (by default all the Z locators, or the last field)
    usesel : Boolean to indicate if the selection has to be considered
    flagColorBar: Flag for representing the Color Bar (not represented if alpha=0)
    aspect: Aspect ratio of the axes scaling, i.e. y/x-scale.
    xlim: Bounds defined along the first axis
    ylim: Bounds defined along the second axis
    size: Size of the data points (default 20)
    cmap: Optional color scale
    title: Title given to the figure (each subplot is titled with the name of the variable represented)
    axs: References for the subplots within the figure: list or array (1D or 2D) containing at least Nvar Axes.
    figsize: (if ax is None) Sizes (width, height) of figure (in inches)
    end_plot: Flag for closing the graphics
    
    **plot_args : arguments passed to matplotlib.pyplot.pcolormesh for every grid plots
    '''
    
    if names is None:
        names = db.getNamesByLocator(gl.ELoc.Z) # all Z locators
        if names == () : # if no Z locator, choose the last field
            names = db.getLastName()
    else:
        names = db.getNames(names)
    names = np.atleast_1d(names)
    
    Nplots = len(names)
    if Nplots == 1:
        axs = point(db, color_name=names[0], usesel=usesel, flagColorBar=flagColorBar, xlim=xlim, ylim=ylim, 
                   size=size, cmap=cmap, aspect=aspect,
                 title=title, ax=axs, figsize=figsize, end_plot=end_plot, **plot_args)

    elif Nplots == 0:
        print("There is no variable in the dbgrid corresponding to the name given.")
        axs = None
    
    else:
        if axs is None:
            nlines, ncols = shape_Nsubplots(Nplots)
            fig, axs = newFigure(figsize, xlim, ylim, nx=nlines, ny=ncols, sharex=True, sharey=True)
        
        if title is not None:
            fig.suptitle(title)
        
        for i,ax in enumerate(axs.flat):
            if i >= Nplots:
                ax.set_visible(False)
                continue;
            
            point(db, color_name=names[i], ax=ax, usesel=usesel, flagColorBar=flagColorBar,
                  aspect=aspect, cmap=cmap, size=size, xlim=xlim, ylim=ylim, 
                  title=names[i], **plot_args)
            
        if end_plot:
            plt.show()
        
    return axs

def size_plots(db, names = None, usesel = True, flagColorBar=True, aspect=None,
               xlim=None, ylim=None, color='r', sizmin=20, sizmax=200,
               title = None, axs=None, figsize = None, end_plot=False, **plot_args):
    '''
    Function for plotting several variables (referred by their names) informed in a Db in subplots (Nsubplots=Nvariables)
    Variables are represented with size plots, i.e. each data point is represented with a size corresponding to the data value.

    db: Db containing the variable to be plotted
    names: Name of the variables to be represented (by default all the Z locators, or the last field)
    usesel : Boolean to indicate if the selection has to be considered
    flagColorBar: Flag for representing the Color Bar (not represented if alpha=0)
    aspect: Aspect ratio of the axes scaling, i.e. y/x-scale.
    xlim: Bounds defined along the first axis
    ylim: Bounds defined along the second axis
    color: Color of data points
    sizmin: Size corresponding to the smallest value
    sizmax: Size corresponding to the largest value
    title: Title given to the figure (each subplot is titled with the name of the variable represented)
    axs: References for the subplots within the figure: list or array (1D or 2D) containing at least Nvar Axes.
    figsize: (if ax is None) Sizes (width, height) of figure (in inches)
    end_plot: Flag for closing the graphics
    
    **plot_args : arguments passed to matplotlib.pyplot.pcolormesh for every grid plots
    '''
    
    if names is None:
        names = db.getNamesByLocator(gl.ELoc.Z) # all Z locators
        if names == () : # if no Z locator, choose the last field
            names = db.getLastName()
    else:
        names = db.getNames(names)
    names = np.atleast_1d(names)
    
    Nplots = len(names)
    if Nplots == 1:
        axs = point(db, size_name=names[0], usesel=usesel, flagColorBar=flagColorBar, xlim=xlim, ylim=ylim, 
                   sizmin=sizmin, sizmax=sizmax, color=color, aspect=aspect,
                 title=title, ax=axs, figsize=figsize, end_plot=end_plot, **plot_args)

    elif Nplots == 0:
        print("There is no variable in the dbgrid corresponding to the name given.")
        axs = None
    
    else:
        if axs is None:
            nlines, ncols = shape_Nsubplots(Nplots)
            fig, axs = newFigure(figsize, xlim, ylim, nx=nlines, ny=ncols, sharex=True, sharey=True)
        
        if title is not None:
            fig.suptitle(title)
        
        for i,ax in enumerate(axs.flat):
            if i >= Nplots:
                ax.set_visible(False)
                continue;
            
            point(db, size_name=names[i], ax=ax, usesel=usesel, flagColorBar=flagColorBar, 
                  aspect=aspect, sizmin=sizmin, sizmax=sizmax, color=color, 
                  xlim=xlim, ylim=ylim, title=names[i], **plot_args)
            
        if end_plot:
            plt.show()
        
    return axs



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

import gstlearn.plot as gp

setattr(gl.Db,"plot", gp.point)
setattr(gl.Db,"plot_correlation", gp.correlation)
setattr(gl.Db,"plot_hist", gp.hist)
setattr(gl.Db,"color_plots", gp.color_plots)
setattr(gl.Db,"size_plots", gp.size_plots)

setattr(gl.DbGrid,"plot", gp.grid)
setattr(gl.DbGrid,"plot_grids", gp.grids)
setattr(gl.DbGrid,"plot_point", gp.point)
# plot_correlation and plot_hist are already inherited from the parent class Db

setattr(gl.Vario,"plot", gp.vario)
setattr(gl.Vario,"plot_varioElem", gp.varioElem)
setattr(gl.Vario,"plot_varioDir", gp.varioDir)
setattr(gl.Vario,"plot_varmod", gp.varmod)

setattr(gl.Model,"plot", gp.model)

setattr(gl.Rule,"plot", gp.rule)

setattr(gl.Table,"plot", gp.table)

setattr(gl.Polygons,"plot", gp.polygon)