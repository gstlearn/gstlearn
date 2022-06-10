import matplotlib.pyplot     as plt
import matplotlib.patches    as ptc
import matplotlib.transforms as transform
import matplotlib.colors     as mcolors
import numpy                 as np
import gstlearn              as gl
from mpl_toolkits.axes_grid1 import make_axes_locatable

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

def drawDecor(ax=None, xlab=None, ylab=None, title=None):
    if ax is None:
        if xlab is not None:
            plt.xlabel(xlab)
        if ylab is not None:
            plt.ylabel(ylab)
        if title is not None:
            plt.title(title)
    else:
        if xlab is not None:
            ax.set_xlabel(xlab)
        if ylab is not None:
            ax.set_ylabel(ylab)
        if title is not None:
            ax.set_title(title)
    return

def addColorbar(im, ax):
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = plt.colorbar(im, ax=ax, cax=cax)
    return cbar

def getDefinedValues(db, name, posx=0, posy=1, corner=None, usesel=True, compress=False):

    if db.isGrid():
        if db.getNDim() == 2:
            posx = 0
            posy = 1
        if corner is None:
            corner = gl.ut_vector_int(db.getNDim(),0)
        tabx = db.getOneSlice(name, posx, posy, corner, usesel)
    else:
        tabx = db.getColumn(name, usesel)
    tabx = np.array(tabx).transpose()
    tabx[tabx == gl.getTEST()] = np.nan
    
    if compress:
        tabx = tabx[np.logical_not(np.isnan(tabx))]
        
    return tabx

def getBiDefinedValues(db, name1, name2, usesel=True):
    tabx = db.getColumn(name1, usesel)
    tabx = np.array(tabx).transpose()
    tabx[tabx == gl.getTEST()] = np.nan
    
    taby = db.getColumn(name2, usesel)
    taby = np.array(taby).transpose()
    taby[taby == gl.getTEST()] = np.nan
    
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
              linestyle0='dashed', hmax=None, gmax=None, show_pairs = False,
              flagLabelDir=False, flagLegend=False, flagLabelSill=False,
              title=None, xlabel=None, ylabel=None, label=None,
              ax=None, figsize=None, end_plot = False, **plot_args):
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
    
    **plot_args : arguments passed to matplotlib.pyplot.plot

    Returns
    -------
    ax : axes where the variogram is represented

    """
    color = plot_args.setdefault('color', 'black')
    
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
        
    ax.hlines(vario.getVar(ivar,jvar), 0, hmax, color0, linestyle0, label = labsill)
    
    if show_pairs:
        pairs = vario.getSwVec(idir,ivar,jvar)
        for i in range(len(hh)):
            ax.annotate(str(int(pairs[i])), (hh[i],gg[i]), xytext=(0,5), xycoords = 'data',
                        textcoords = 'offset points', ha='center')
    
    drawDecor(ax, xlabel, ylabel, title)
    
    ax.set_xlim(left=0)
    if ivar == jvar:
        ax.set_ylim(bottom=0)
        
    if flagLegend:
        ax.legend()
        
    if end_plot:
        plt.show()
    
    return ax

def varioDir(vario, ivar=0, jvar=0,
             color0='black', linestyle0='dashed', hmax=None, gmax=None, 
             cmap=None, flagLegend=False, title=None, xlabel=None, ylabel=None, label=None,
             ax=None, figsize=None, 
             end_plot=False, **plot_args):
    """Plot a single directional experimental variogram (all avalaible directions, for fixed variable(s)).
    
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
    
    # xlim = [0,hmax]
    # gmin = 0
    # if ivar != jvar:
    #     gmin = -gmax
    # ylim = [gmin,gmax]
        
    # if ax is None:
    #     fig, ax = newFigure(figsize, xlim, ylim)
    # update_xylim(ax, xlim, ylim)
    if ax is None:
        fig, ax = newFigure(figsize, None, None)
    
    for idirUtil in ndirUtil:
        flagLabelSill = idirUtil == 0
        varioElem(vario, ivar=ivar, jvar=jvar, idir=idirUtil, color=cols(idirUtil), color0=color0, linestyle0=linestyle0, 
                  ax=ax, hmax=hmax, gmax=gmax, figsize=figsize, flagLabelDir=True, flagLabelSill=flagLabelSill, label=label,
                  **plot_args)
        
    drawDecor(ax, xlabel, ylabel, title)
    
    ax.autoscale(True)
    ax.set_xlim(left=0)
    if ivar == jvar:
        ax.set_ylim(bottom=0)
    
    if flagLegend:
        ax.legend()
        
    if end_plot:
        plt.show()
        
    return ax

def varmod(vario, mymodel=None, ivar=-1, jvar=-1, idir=-1,
           linestylem="dashed", color0='black', linestyle0="dotted",
           nh = 100, hmax = None, gmax = None, 
           cmap=None, flagLegend=False, title=None, axs=None, figsize=None, end_plot=False, **plot_args):
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

    # xlim = [0, hmax]
    # ylim = [0, gmax]
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
                          color=cols(idirUtil),
                          color0=color0, linestyle0=linestyle0, 
                          ax=ax, hmax=hmax, gmax=None, 
                          flagLabelDir=flagLabelDir, flagLegend=flagLegend, **plot_args)

                # Plotting the Model (optional)
                if mymodel is not None:
                    codir = vario.getCodir(idirUtil)
                    model(mymodel, ivar=iv, jvar=jv, codir=codir, 
                          color=cols(idirUtil), linestyle=linestylem, color0=color0, linestyle0=linestyle0, ax=ax,
                          hmax=hmax, gmax=None, nh=nh,
                          flagLabelDir=flagLabelDir, flagLegend=flagLegend)

    ax.autoscale(True)
    ax.set_xlim(left=0)
    if ivar == jvar:
        ax.set_ylim(bottom=0)
    
    if title is not None:
        plt.suptitle(title)
        
    if end_plot:
        plt.show()
        
    return axs

def vario(vario, ivar=-1, jvar=-1, idir=-1,
          color0='black', linestyle0='dashed', hmax=None, gmax=None, 
          cmap = None, flagLegend=False, title = None, axs = None, figsize = None, end_plot=False, **plot_args):
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
                 color0=color0, linestyle0=linestyle0, hmax=hmax, gmax=gmax, cmap=cmap, 
                 flagLegend=flagLegend, title=title, axs=axs, figsize=figsize, end_plot=end_plot, **plot_args)
    
    return axs

def model(model, ivar=0, jvar=0, codir=None, color0='black', linestyle0='dashed',
          nh = 100, flagEnv = True, hmax = None, gmax = None, 
          flagLabelDir=False, flagLegend=False, title=None, xlabel=None, ylabel=None, ax=None, 
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
    title : Optional title for the axes.
    ax : Reference for the plot within the figure. If None (default), it creates a new figure.
    figsize : (if ax is None) Sizes (width, height) of figure (in inches).
    end_plot : Flag for closing the graphics (the default is False).

    Returns
    -------
    ax : axes where the variogram is represented

    """
    color = plot_args.setdefault('color', 'black')    
    
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
    gg = model.sample(hmax, nh, ivar, jvar, codir, addZero=True)
    
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
        ggp = model.sample(hmax, nh, ivar, jvar, codir, 1, addZero=True)
        ax.plot(hh[istart:], ggp[istart:], color = color0, linestyle = linestyle0, label="plus")
        ggm = model.sample(hmax, nh, ivar, jvar, codir,-1, addZero=True)
        ax.plot(hh[istart:], ggm[istart:], color = color0, linestyle = linestyle0, label="minus")
    
    drawDecor(ax, xlabel, ylabel, title)
    
    ax.autoscale(True)
    ax.set_xlim(left=0)
    if ivar == jvar:
        ax.set_ylim(bottom=0)
    
    if flagLegend:
        ax.legend()
        
    if end_plot:
        plt.show()

    return ax

def point(db, 
          color_name=None, size_name=None, usesel=True, 
          color='r', size=20, sizmin=10, sizmax=200, 
          xlim=None, ylim=None, directColor=False,
          cmap=None, flagColorBar=True, flagSizeLegend=True, aspect='auto',
          title=None, ax=None, figsize = None, end_plot =False, **scatter_args):
    '''Function for plotting a point data base, with optional color and size variables
    
    db: Db containing the variable to be plotted
    color_name: Name of the variable containing the color per sample
    size_name: Name of the variable containing the size per sample
    usesel : Boolean to indicate if the selection has to be considered
    color: Constant color (used if 'color_name' is not defined)
    size: Constant size (used if 'size_name' is not defined)
    sizmin: Size corresponding to the smallest value (used if 'size_name' is defined)
    sizmax: Size corresponding to the largest value (used if 'size_name' is defined)
    xlim: Bounds defined along the first axis
    ylim: Bounds defined along the second axis
    directColor: True if the value of the field directly indicates the Rank in the Color Scale
    cmap: Optional Color scale
    flagColorBar: Flag for representing the Color Bar (not represented if color_name=None)
    flagSizeLegend: Flag for representing the Legend for marker size (not represented if size_name=None)
    title: Title given to the plot
    ax: Reference for the plot within the figure
    figsize: (if ax is None) Sizes (width, height) of figure (in inches)
    end_plot: Flag for closing the graphics

    **scatter_args : arguments passed to matplotllib.pyplot.scatter
    '''
    
    if ax is None:
        fig, ax = newFigure(figsize, xlim, ylim)
    else:
        ax.set_xlim(xlim)
        ax.set_ylim(ylim)

    # Extracting coordinates
    tabx = db.getCoordinates(0,usesel)
    taby = db.getCoordinates(1,usesel)
    if len(tabx) <= 0 or len(taby) <= 0:
        return
    
    # Color of symbol
    if color_name is not None:
        colval = getDefinedValues(db, color_name, 0, 1, None, usesel, False)
        if (directColor and cmap is not None):
            colval = colval.astype(int)
            colval = cmap(colval)
    else:
        colval = color

    # Size of symbol
    if size_name is not None:
        sizval = getDefinedValues(db, size_name, 0, 1, None, usesel, False)
        m = np.nanmin(np.absolute(sizval))
        M = np.nanmax(np.absolute(sizval))
        sizval = (sizmax - sizmin) * (np.absolute(sizval) - m) / (M-m) + sizmin
        mask = ~np.isnan(sizval)
        sizval = sizval[mask]
        tabx = np.array(tabx)[mask]
        taby = np.array(taby)[mask]
        if color_name is not None:
            colval = colval[mask]
            color = cmap
    else:
        sizval = size

    im = ax.scatter(x = tabx, y = taby, s = sizval, c = colval, cmap=cmap, **scatter_args)
    ax.set_aspect(aspect)
    
    if flagColorBar and (color_name is not None):
        addColorbar(im, ax)
    
    if flagSizeLegend and (size_name is not None):
        labels = lambda marker_size : (M - m)*(marker_size - sizmin)/(sizmax - sizmin) + m
        # labels is the inverse transformation from marker sizes to variable values
        ax.legend(*im.legend_elements("sizes", num=5, color=cmap, func=labels))
         
    if title is not None:
        ax.set_title(title)
    
    if end_plot:
        plt.show()

    return ax

def polygon(poly, faceColor='yellow', edgeColor = 'blue', 
            colorPerSet = False, flagEdge=True, flagFace=False, linewidth=2,
            title= None, ax=None, figsize=None, end_plot=False, **fill_args):
    '''Function to display a polygon
    
    **fill_args: arguments passed to ax.fill'''
    
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
        
    if title is not None:
        ax.set_title(title)
        
    if end_plot:
        plt.show()
        
    return ax
        
def grid(dbgrid, name = None, usesel = True, flagColorBar=True, aspect='equal',
         xlim=None, ylim=None, posx=0, posy=1, corner=None,
         title = None, ax=None, figsize = None, end_plot=False, **plot_args):
    '''
    Function for plotting a variable (referred by its name) informed in a DbGrid

    dbgrid: Db, organized as a Grid, containing the variable to be plotted
    name: Name of the variable to be represented (by default, the first Z locator, or the last field)
    usesel : Boolean to indicate if the selection has to be considered
    flagColorBar: Flag for representing the Color Bar (not represented if alpha=0)
    xlim: Bounds defined along the first axis
    ylim: Bounds defined along the second axis
    title: Title given to the plot
    ax: Reference for the plot within the figure
    figsize: (if ax is None) Sizes (width, height) of figure (in inches)
    end_plot: Flag for closing the graphics
    
    **plot_args : arguments passed to matplotlib.pyplot.pcolormesh
    '''
    # define some default values (add them to hist_args if not already set)
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
    # else:
    #     ax.set_xlim(xlim)
    #     ax.set_ylim(ylim)
        
    x0 = dbgrid.getX0(posx)
    y0 = dbgrid.getX0(posy)
    nx = dbgrid.getNX(posx)
    ny = dbgrid.getNX(posy)
    dx = dbgrid.getDX(posx)
    dy = dbgrid.getDX(posy)
    angles = dbgrid.getAngles()
    
    data = getDefinedValues(dbgrid, name, posx, posy, corner, usesel, False)
    data   = np.reshape(data, (ny,nx))

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

    im = ax.pcolormesh(X, Y, data, **plot_args)
    
    ax.set_aspect(aspect)
    im.set_transform(trans_data)
    
    x1, x2, y1, y2 = x0, X[-1], y0, Y[-1]
    ax.plot([x1, x2, x2, x1, x1], [y1, y1, y2, y2, y1], "red", transform=trans_data)
    
    update_xylim(ax, xlim=xlim, ylim=ylim) 
    
    if flagColorBar:# and alpha > 0:
        addColorbar(im, ax)
    
    if title is not None:
        ax.set_title(title)
    else:
        ax.set_title(dbgrid.getNames(name)[0])
    
    if end_plot:
        plt.show()
    
    return ax

def grids(dbgrid, names = None, usesel = True, flagColorBar=True, aspect='equal',
         xlim=None, ylim=None, norm=None,
         title = None, axs=None, figsize = None, end_plot=False, **plot_args):
    '''
    Function for plotting several variables (referred by their names) informed in a grid Db in subplots (Nsubplots=Nvariables)

    dbgrid: Db, organized as a Grid, containing the variable to be plotted
    names: Name of the variables to be represented (by default all the Z locators, or the last field)
    usesel : Boolean to indicate if the selection has to be considered
    flagColorBar: Flag for representing the Color Bar (not represented if alpha=0)
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
        ax = grid(dbgrid, name=names[0], usesel=usesel, flagColorBar=flagColorBar, xlim=xlim, ylim=ylim, 
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
            
        grid(dbgrid, name=names[i], ax=ax, usesel=usesel, flagColorBar=flagColorBar, aspect=aspect,
                 xlim=xlim, ylim=ylim, norm=norm_i, **plot_args)
    if end_plot:
        plt.show()
        
    return axs

def hist_tab(val, xlab=None, ylab=None, 
             title = None, ax = None, figsize=None, end_plot=False, **hist_args):
    '''Function for plotting the histogram of an array (argument 'val')
    
    hist_args : arguments passed to matplotlib.pyplot.hist'''
    
    # define some default values (add them to hist_args if not already set)
    color     = hist_args.setdefault("color", "yellow")
    edgecolor = hist_args.setdefault("edgecolor", "red")
    bins      = hist_args.setdefault("bins", 30)
    
    if ax is None:
        fig, ax = newFigure(figsize)
        
    ax.hist(val, **hist_args)
    
    drawDecor(ax, xlab=xlab, ylab=ylab, title=title)
        
    if end_plot:
        plt.show()
        
    return ax
    
def hist(db, name, xlab=None, ylab=None, title = None, ax=None,
         figsize=None, end_plot=False, usesel=True, **hist_args):
    '''Function for plotting the histogram of a variable contained in a Db
    
    hist_args : arguments passed to matplotlib.pyplot.hist'''
    
    db.useSel = usesel
    #val = db.getColumn(name)
    val = db[name]
    if len(val) == 0:
        return None
    
    if title is None:
        title = db.getNames(name)[0]
    ax = hist_tab(val, title=title, xlab=xlab, ylab=ylab, ax=ax, figsize=figsize, end_plot=end_plot, **hist_args)
    
    return ax

def curve(data, icas=1, color='black', 
          title=None, ax=None, figsize = None, end_plot=False, **plot_args):
    '''
    Function for plotting the curve of an array (argument 'data')
        icas=1 when 'data' represent the abscissa and 2 when 'data' represents the ordinate
    **plot_args : arguments passed to matplotlib.pyplot.plot
    '''
    color = plot_args.setdefault('color', 'black')
    label = plot_args.setdefault('label', 'model')
    
    if ax is None:
        fig, ax = newFigure(figsize)
        
    nbpoint = len(data)
    regular = [i for i in range(nbpoint)]
    
    if icas == 1:
        ax.plot(data, regular, **plot_args)
    else:
        ax.plot(regular, data, **plot_args)
    
    drawDecor(ax, title=title)
    ax.legend()
        
    if end_plot:
        plt.show()
        
    return ax

def XY(xtab, ytab, flagAsPoint=False, xlim=None, ylim=None, 
       title=None, ax=None, figsize = None, end_plot=False, **plot_args):
    """Plot Y against X.
    **plot_args : arguments passed to matplotlib.pyplot.plot"""
    
    if not len(ytab) == len(xtab):
        print("Arrays 'xtab' and 'ytab' should have same dimensions")
        return None;
    
    if ax is None:
        fig, ax = newFigure(figsize, xlim, ylim)
    
    if flagAsPoint:
        plot_args.setdefault('label', 'point')
        plot_args.setdefault('markersize', 10)
        plot_args.setdefault('color', 'blue')
        plot_args.setdefault('marker', 'o')
        ax.plot(xtab, ytab, **plot_args)
    else:
        plot_args.setdefault('label', 'model')
        #plot_args.setdefault('fmt', 'g-') # TODO : to be restored ? (doesn't work with python 3.9.2 / mpl 3.5.1)
        ax.plot(xtab, ytab, **plot_args)
    
    drawDecor(ax, title=title)
    ax.legend()
    
    if end_plot:
        plt.show()
    
    return ax

def rule(rule, proportions=[], 
         title=None, ax=None, figsize=None, end_plot=False):
    
    if ax is None:
        fig, ax = newFigure(figsize, [-10.,10.], [-10.,10.])
        
    nfac = rule.getFaciesNumber()
    rule.setProportions(proportions)
    cols = get_cmap(nfac)

    for ifac in range(nfac):
        bds = rule.getThresh(ifac+1)
        rect = ptc.Rectangle((bds[0],bds[2]),bds[1]-bds[0], bds[3]-bds[2], 
                              color=cols(ifac))
        ax.add_patch(rect)

    drawDecor(ax, title=title)
        
    if end_plot:
        plt.show()

    return ax


def table(table, icols, fmt='ok', xlim=None, ylim=None,
          title=None, ax=None, figsize=None, end_plot=False, **plot_args):
    '''
    Function for plotting the contents of a Table (argument 'tablr')
        icols designates the ranks of the variable (0: ordinate; 1: abscissae [or regular]) 
        fmt designates [marker][line][color] information
    **plot_args
    '''
    
    if ax is None:
        fig, ax = newFigure(figsize, xlim, ylim)
    
    if len(icols) == 1:
        datay = table.getCol(int(icols[0]))
        datax = [i for i in range(table.getRowNumber())]
    else:
        datay = table.getCol(int(icols[0]))
        datax = table.getCol(int(icols[1]))
    
    plot_args.setdefault('label', 'table')
    ax.plot(datax, datay, fmt, **plot_args)
    
    drawDecor(ax, title=title)
    ax.legend()
    
    if end_plot:
        plt.show()
        
    return ax

def mesh(mesh, 
         #color='r', size=20, sizmin=10, sizmax=200,#arguments not used ?
         flagEdge=True, flagFace=False,
         xlim=None, ylim=None, 
         title=None, ax=None, figsize = None, end_plot =False, **plot_args):
    """
    **plot_args : arguments passed to matplotlib.pyplot.fill
    """
    
    if ax is None:
        if xlim is not None:
            plt.xlim(mesh.getExtrema(0))
        if ylim is not None:
            plt.ylim(mesh.getExtrema(1))
        fig, ax = newFigure(figsize, xlim, ylim)   

    nmesh = mesh.getNMeshes()
    for imesh in range(nmesh):
        tabx = mesh.getCoordinatesPerMesh(imesh, 0, True)
        taby = mesh.getCoordinatesPerMesh(imesh, 1, True)
            

    if flagFace:
        plot_args.setdefault('facecolor', 'yellow')
                
    if flagEdge:
        plot_args.setdefault('edgecolor', 'blue')

    plot_args.setdefault('linewidth', 1)
    ax.fill(tabx, taby, **plot_args)

    drawDecor(ax, title=title)
    
    if end_plot:
        plt.show()

    return ax

def correlation(db, namex, namey, bins=50, xlim=None, ylim=None, usesel=True, asPoint = False,
                diagLine = False, 
                xlab=None, ylab=None, title = None, ax=None, figsize=None, end_plot=False):
    '''Function for plotting the scatter plot between two variables contained in a Db'''
 
    if ax is None:
        fig, ax = newFigure(figsize, xlim, ylim)
   
    tabx, taby = getBiDefinedValues(db, namex, namey, usesel)
    if len(tabx) == 0:
        return None
    if len(taby) == 0:
        return None
    
    if asPoint:
        ax.scatter(tabx, taby)
    else:
        ax.hist2d(tabx, taby, bins, cmin=1)


    if diagLine:
        u=[np.min(tabx),np.min(taby)]
        v=[np.max(tabx),np.max(taby)]
        ax.plot(u,v,c="r")

    drawDecor(ax, xlab=xlab, ylab=ylab, title=title)
    
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
        


