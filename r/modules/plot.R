#Set of global values
default_working_mode <<- FALSE # False for old and True for new

default_size <<- c(c(8,8), c(8,8))
default_xlim <<- c( NA, NA )
default_ylim <<- c( NA, NA )
default_sameLim <<- c( FALSE, FALSE )
default_aspect <- c('auto', 1 )

ensure_dependencies <- function() {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("The ggplot2 package must be installed to use R plots functionality")
  }
  if (!requireNamespace("ggpubr", quietly = TRUE)) {
    stop("The ggpubr package must be installed to use R plots functionality")
  }
}

end.func <- function(p, end.plot=TRUE)
{
  padd = p
  if (end.plot)
  {
    print(padd)
    padd = NULL
    invisible()
  } else {
    padd
  }
}

get.colors <- function()
{
  c("blue", "red", "green", "brown", "orange", "purple", "yellow")
}

getFigure <- function(padd = NULL)
{
  ensure_dependencies()
  if (length(padd) > 0)
    p <- padd
  else
    p <- ggplot2::ggplot()
  p
}

decor <- function(p, xlab = "", ylab = "", asp = NULL, title = "")
{
  ensure_dependencies()
  if (xlab != "")
    p <- p + ggplot2::labs(x = xlab)
  if (ylab != "")
    p <- p + ggplot2::labs(y = ylab)
  if (title != "")
    p <- p + ggplot2::ggtitle(title) + ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  if (! is.null(asp))
    suppressWarnings(suppressMessages(p <- p + ggplot2::coord_fixed(ratio = asp)))
  p
}

# Function for representing a Model
plot.model <- function(model, vario=NULL, hmax=1, codir=NULL, 
					   ivar=0, jvar=0, idir=0, asCov=FALSE, 
                       xlab = "", ylab = "",title="", nh=100, padd=NULL, end.plot=TRUE)
{
  ensure_dependencies()
  if (! is.null(vario))
  {
    codir = vario$getCodirs(idir)
    hmax = vario$getHmax(ivar, jvar, idir)
  }
  else
  {
  if (is.null(codir))
    {
      ndim = model$getDimensionNumber()
      codir = rep(0,ndim)
      codir[1] = 1
    }
  }

  p <- getFigure(padd)
  
  hh = seq(from=0, to=hmax, length.out=nh)
  gg = model$sample(hmax, nh, ivar, jvar, codir, asCov=asCov)
  df = data.frame(cbind(hh,gg))
  
  p <- p + geom_line(data = df, aes(x=hh,y=gg), na.rm=TRUE)
  
  p <- decor(p, xlab = xlab, ylab = ylab, title = title)
  
  end.func(p, end.plot)
}
setMethod("plot", signature(x="_p_Model"), function(x,y="missing",...) plot.model(x,...))

# Function for representing the Experimental Variogram together with the Model (optional)

plot.varmod <- function(vario, model=NULL, ivar=-1, jvar=-1, idir=-1,
                        asCov=FALSE, nh=100, draw_psize=FALSE, draw_plabels=FALSE, 
                        color_psize="black", ratio_psize=10,
                        color_plabel="black", size_plabel=2, nudge_y=0.1,
                        title="", show.legend=FALSE, end.plot=TRUE, ...)
{
  ensure_dependencies()
  ndir = vario$getDirectionNumber()
  nvar = vario$getVariableNumber()
  cols = get.colors()
  idirUtil = seq(0,ndir-1)
  if (idir >= 0) idirUtil = idir
  ivarUtil = seq(0, nvar-1)
  if (ivar >= 0) ivarUtil = ivar
  ivarN = length(ivarUtil)
  jvarUtil = seq(0, nvar-1)
  if (jvar >= 0) jvarUtil = jvar
  jvarN = length(jvarUtil)
  
  # Loop on the variables
  
  plot_lst <- vector("list", length = ivarN * jvarN)
  index = 0
  for (iv in ivarUtil)
    for (jv in jvarUtil)
      {
        index = index + 1
        g = ggplot2::ggplot()
        
        if (iv < jv) 
        {
          plot_lst[[index]] <- g
          next
        }
        
        xlim = c(0,0)
        ylim = c(0,0)
        for (id in idirUtil)
        {
          sill = vario$getVar(iv,jv)
          nlag = vario$getLagNumber(id)
          sw = vario$getSwVec(id,iv,jv)
          gg = vario$getGgVec(id,iv,jv, asCov=asCov)
          hh = vario$getHhVec(id,iv,jv)
          hmax = max(hh)
          gmax = max(abs(gg))
          if (abs(sill) > gmax) gmax = abs(sill)
          gmax = gmax * 1.1

          # Bounds
          if (hmax > xlim[2]) xlim[2] = hmax 
          gmin = 0
          if (iv != jv) gmin = -gmax
          if (gmin < ylim[1]) ylim[1] = gmin
          if (gmax > ylim[2]) ylim[2] = gmax
                
          # Plotting the experimental variogram
          df = data.frame(cbind(hh,gg))
          g <- g + ggplot2::geom_line(data = df, ggplot2::aes(x=hh,y=gg), color=cols[id+1], na.rm=TRUE)
          
          if (draw_psize)
            g <- g + ggplot2::geom_point(data = df, ggplot2::aes(x=hh, y=gg), 
                                         size=sw/ratio_psize, color=color_psize, 
                                         na.rm=TRUE, show.legend=show.legend)
            
          if (draw_plabels)
             g <- g + ggplot2::geom_text(data = df, ggplot2::aes(x=hh, y=gg, label=as.character(sw)),
                                         color=color_plabel, size=size_plabel, nudge_y=nudge_y, 
                                         show.legend=show.legend, check_overlap=TRUE)
   
          # Plotting the Model (optional)
          if (! is.null(model))
          {
            hh = seq(0, hmax, hmax/nh)
            nhh = length(hh)
            codir = vario$getCodirs(id)
            gg = model$sample(hmax, nhh, iv, jv, codir, asCov=asCov)
            dfg = data.frame(cbind(hh,gg))
            g <- g + ggplot2::geom_line(data = dfg, ggplot2::aes(x=hh,y=gg), color=cols[id+1], size=1, na.rm=TRUE) 

            if (iv != jv)
            {
              ggp = model$sample(hmax, nhh, iv, jv, codir, 1)
              dfg = data.frame(cbind(hh,ggp))
              g <- g + ggplot2::geom_line(data = dfg, ggplot2::aes(x=hh,y=ggp), color=cols[id+1],
                                          linetype = 'twodash', na.rm=TRUE) 
              ggm = model$sample(hmax, nhh, iv, jv, codir,-1)
              dfm = data.frame(cbind(hh,ggp))
              g <- g + ggplot2::geom_line(data = dfm, ggplot2::aes(x=hh,y=ggm), color=cols[id+1],
                                          linetype = 'twodash', na.rm=TRUE) 
            }
          } 
        } # End of loop on Directions
        
        g <- g + ggplot2::scale_x_continuous("Distance", limits = xlim, expand = c(0,0))
        if (iv == jv)
          g <- g + ggplot2::scale_y_continuous("Variogram", limits = ylim, expand = c(0,0))
        else
          g <- g + ggplot2::scale_y_continuous("Cross-Variogram", limits = ylim, expand = c(0,0))
                
        # Plotting relevant control lines
        if (iv != jv)
          g <- g + ggplot2::geom_hline(yintercept = 0.)
        g <- g + ggplot2::geom_hline(yintercept = sill, linetype = 'longdash')
        plot_lst[[index]] <- g
      }
  p = ggplot2::ggarrange(plotlist=plot_lst, nrow=ivarN, ncol = jvarN)
  
  p <- decor(p, title = title)
  
  end.func(p, end.plot)
}
setMethod("plot", signature(x="_p_Vario"), function(x,y,...) plot.varmod(x,...))

# Function for plotting a point data base, with optional color and size variables
plot.point <- function(db, color_name=NULL, size_name=NULL, label_name=NULL,
                       color0='red', size0=0.2, color_label="black", nudge_y=0.1,
                       sizmin=10, sizmax=100, flagAbsSize = FALSE, 
                       show.legend.color=FALSE, legend.name.color="P-Color",
                       show.legend.size =FALSE, legend.name.size ="P-Size",
                       show.legend.label=FALSE, legend.name.label="P-Label",
                       asp=1, xlab="", ylab="", title="", 
                       padd = NULL, end.plot=TRUE, ...) 
{  
  ensure_dependencies()
  # Creating the necessary data frame
  np   = db$getSampleNumber(TRUE)
  xtab = db$getCoordinates(0,TRUE)
  ytab = db$getCoordinates(1,TRUE)
    
  # Color of symbol
  if (! is.null(color_name))
  {
    colval  = db$getColumn(color_name)
  }
  else
  {
    colval = rep(color0,np)
  }

  # Size of symbol
  if (! is.null(size_name))
  {
  	reduction = 100
    sizval  = db$getColumn(size_name)
    if (flagAbsSize) sizval = abs(sizval)
    m = min(sizval,na.rm=TRUE)
    M = max(sizval,na.rm=TRUE)
    sizval = (sizmax * (sizval - m) / (M-m) + sizmin) / reduction
  }
  else
  {
    sizval = rep(size0,np)
  }

  # Label of symbols
  if (! is.null(label_name))
  {
    label_round = 2
    labval  = round(db$getColumn(label_name,TRUE),label_round)
  }
  else
  {
    labval = rep(0,np)
  }
  df = data.frame(xtab,ytab)
  
  p <- getFigure(padd)
  
  p <- p + geom_point(data=df, aes(x=xtab,y=ytab), color=colval, size=sizval,
           na.rm=TRUE)
  
  if (! is.null(label_name)) 
  {
    p <- p + geom_text(data = df, aes(x=x, y=y), label=as.character(labval),
               nudge_y = nudge_y, color=color_label, check_overlap=TRUE)
               
    if (show.legend.label) {
      p <- p + ggplot2::guides(label = ggplot2::guide_legend(title = legend.name.label))
    } else {
      p <- p + ggplot2::guides(label = "none")
    }
  }
  
  if (show.legend.color) {
    p <- p + ggplot2::guides(color = ggplot2::guide_legend(title = legend.name.color))
  } else {
    p <- p + ggplot2::guides(color = "none")
  }
    
  if (show.legend.size) {
    p <- p + ggplot2::guides(size = ggplot2::guide_legend(title = legend.name.size))
  } else {
    p <- p + ggplot2::guides(size = "none")
  }
      
  p <- decor(p, xlab = xlab, ylab = ylab, asp = asp, title = title)
  
  end.func(p, end.plot)
}

# Function for plotting a variable (referred by its name) informed in a grid Db
#
# option Indicates the color map (from "A", "B", "C", "D", "E", "F", "G", "H")
plot.grid <- function(dbgrid, name=NULL, na.color = "white", asp=1,
      option="B", zlim = NULL, useSel = TRUE,
      show.legend=TRUE, legend.name="",
      xlab="", ylab="", title="", 
      padd=NULL, end.plot=TRUE)
{
  ensure_dependencies()
  if (! dbgrid$isGrid())
  {
    cat("This function is restricted to Grid Db and cannot be used here")
    stop()
  }

  if (is.null(name))
  {
    if (dbgrid$getLocatorNumber(ELoc_Z()) > 0) 
      name = dbgrid$getNameByLocator(ELoc_Z())
    else
      name = dbgrid$getLastName()
  }

  # Building the necessary data frame
  x = dbgrid$getColumnByLocator(ELoc_X(),0, FALSE, FALSE)
  y = dbgrid$getColumnByLocator(ELoc_X(),1, FALSE, FALSE)
  data = dbgrid$getColumn(name, useSel, FALSE)
  if (length(data) <= 0)
  {
    cat("Variable",name,"does not exist\n")
    stop()
  }
  
  p <- getFigure(padd)
  
  # Define the contents
  if (dbgrid$getAngles()[1] == 0)
  {
    df = data.frame(x,y,data)
  	p <- p + geom_tile(data = df, aes(x = x, y = y, fill = data))
  }
  else
  {
    ids = seq(1, dbgrid$getNTotal())
    coords = dbgrid$getAllCellsEdges()
    positions = data.frame(id = rep(ids, each=4),x=coords[[1]],y=coords[[2]])
    values = data.frame(id = ids, value = data)
    df <- merge(values, positions, by = c("id"))
    p <- p + ggplot2::geom_polygon(data = df, ggplot2::aes(x = x, y = y, fill = value, group = id))
  }
  
  # Define the color scale
  p = p + ggplot2::scale_fill_viridis_c(option = option, na.value = na.color, limits=zlim)
  
  if (show.legend) {
    p <- p + ggplot2::guides(fill = ggplot2::guide_colorbar(title=legend.name, reverse=FALSE))
  } else {
    p <- p + ggplot2::guides(fill = "none")
  }
       
  p <- decor(p, xlab = xlab, ylab = ylab, asp=asp, title = title)
  
  end.func(p, end.plot)
}

plot.db <- function(db, padd=NULL, end.plot=TRUE, ...)
{
  ensure_dependencies()
  if (db$isGrid())
    p = plot.grid(db, padd=padd, end.plot=end.plot, ...)
  else
    p = plot.point(db, padd=padd, end.plot=end.plot, ...)
  p
}

setMethod("plot", signature(x="_p_Db"), function(x,padd=NULL,...) plot.db(x,padd,...))

# Function to display a polygon (not tested)

plot.polygon <- function(poly, xlab="", ylab="", title="", color="black", 
		fill=NA, asp=1, padd = NULL, end.plot=TRUE)
{
  ensure_dependencies()
  npol = poly$getPolySetNumber()
  
  p <- getFigure(padd)
  
  for (ipol in 1:npol)
  {
    xtab = poly$getX(ipol-1)
    ytab = poly$getY(ipol-1)
    rp = data.frame(xtab, ytab)
    p <- p + geom_polygon(data = rp, aes(x=xtab,y=ytab), color=color, fill=fill)
  }  
  
  p <- decor(p, xlab = xlab, ylab = ylab, asp=asp, title = title)
  
  end.func(p, end.plot)
}
setMethod("plot", signature(x="_p_Polygons"), function(x,y=missing,...) plot.polygon(x,...))
        
# Function for plotting the histogram of a variable
plot.hist <- function(db, name, nbins=30, col='grey', fill='yellow',
                      xlab="", ylab="", title="", 
                      padd = NULL, end.plot=TRUE)
{
  ensure_dependencies()
  val  = dbg$etColumn(name)
  rp = data.frame(val)
    
  p <- getFigure(padd)
     
  p <- p + ggplot2::geom_histogram(data=rp, ggplot2::aes(x=val), bins=nbins, color=col, fill=fill,
                                   na.rm=TRUE) 
  
  p <- decor(p, xlab = xlab, ylab = ylab, title = title)
  
  end.func(p, end.plot)
}

# Function for plotting histogram for a table of values
plot.hist_tab <- function(val, nbins=30, xlab="", ylab="", title="", 
                          padd=FALSE, end.plot=TRUE)
{
  ensure_dependencies()
  rp = data.frame(val)
  
  p <- getFigure(padd)
     
  p <- p + ggplot2::geom_histogram(data = rp, ggplot2::aes(x=val), bins=nbins, color='grey', fill='yellow') 

  p <- decor(p, xlab = xlab, ylab = ylab, title = title)

  end.func(p, end.plot)
}

# Function for plotting a curve of regularly sampled values
plot.curve <- function(data, color="black", xlab="", ylab="", title="", 
                       padd=NULL, end.plot=TRUE)
{
  ensure_dependencies()
  nbpoint = length(data)
  absc = seq(1,nbpoint)
  rp = data.frame(absc,data)
  
  p <- getFigure(padd)
    
  p <- p + ggplot2::geom_line(data = rp, ggplot2::aes(x=absc,y=data), color=color, na.rm=TRUE)
  
  p <- decor(p, xlab = xlab, ylab = ylab, title = title)
  
  end.func(p, end.plot)
}

# Function for representing a line between points provided as arguments
plot.XY <-function(xtab, ytab, join=TRUE,
                   color="black", linetype="solid", shape=20,
                   flagDiag = FALSE, 
                   diag_color = "red", diag_line = "solid",
                   xlim="", ylim="", xlab="", ylab="", title="", 
                   padd=NULL, end.plot=TRUE)
{
  ensure_dependencies()
  if (length(ytab) != length(xtab))
  {
    cat("Arrays 'xtab' and 'ytab' should have same dimensions")
    stop()
  }
  rp = data.frame(xtab, ytab)

  p <- getFigure(padd)
     
  if (is.numeric(xlim) && length(xlim) == 2)
    p <- p + ggplot2::scale_x_continuous(limits = xlim, expand = c(0,0))
  if (is.numeric(ylim) && length(ylim) == 2)
    p <- p + ggplot2::scale_y_continuous(limits = ylim, expand = c(0,0))
  
  if (flagDiag)
  {
    u = min(xtab, ytab, na.rm=TRUE)
    v = max(xtab, ytab, na.rm=TRUE)
    p <- p + ggplot2::geom_segment(ggplot2::aes(x=u,y=u,xend=v,yend=v),
                                   linetype = diag_line, color = diag_color, na.rm=TRUE)
  }
  
  if (join)
    p <- p + ggplot2::geom_line(data = rp, ggplot2::aes(x=xtab,y=ytab), 
                                linetype = linetype, color=color, na.rm=TRUE)
  else 
    p <- p + ggplot2::geom_point(data = rp, ggplot2::aes(x=xtab,y=ytab), 
                                 shape=shape, color=color, na.rm=TRUE)
  
  p <- decor(p, xlab = xlab, ylab = ylab, title = title)
  
  end.func(p, end.plot)
}

# Function for representing an anamorphosis
plot.anam <- function(anam, ndisc=100, aymin=-10, aymax=10, 
                      color="black", linetype="solid",
                      xlim="", ylim="", xlab="Y", ylab="Z", title="", 
                      padd=NULL, end.plot=TRUE)
{
  ensure_dependencies()
  res = anam$sample(ndisc, aymin, aymax)
  valY = res$getY()
  valZ = res$getZ()
  
  p = plot.XY(valY, valZ, join=TRUE, flagDiag = FALSE,
              color=color, linetype=linetype, 
              xlim=res$getAylim(), ylim=res$getAzlim(), xlab=xlab, ylab=ylab, 
              title=title, 
              padd=padd, end.plot=FALSE)
  
  end.func(p, end.plot)
}

# Function for representing a scatter plot
plot.correlation <- function(db1, name1, name2, db2=NULL, useSel=FALSE,
							 flagDiag = FALSE,
                             color="black", linetype = "solid",
                             diag_color = "red", diag_line = "solid",
                             xlim="", ylim="", xlab="", ylab="", title="", 
                             padd=NULL, end.plot=TRUE)
{
  ensure_dependencies()
  if (is.null(db2)) db2 = db1
  val1 = db1$getColumn(name1, useSel)
  val2 = db2$getColumn(name2, useSel)
  p = plot.XY(val1, val2, join=FALSE, flagDiag=flagDiag, 
              color = color, linetype = linetype, 
              diag_color = diag_color, diag_line = diag_line,
              xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, title=title, 
              padd=padd, end.plot=FALSE)
  
  end.func(p, end.plot) 
}

# Representing a Lithotype rule
plot.rule <- function(rule, proportions=NULL, xlab="", ylab="", title="",
                      padd=NULL, end.plot=TRUE)
{
  ensure_dependencies()
  nrect = rule$getFaciesNumber()
  if (! is.null(proportions)) 
    rule$setProportions(proportions)
  else
    rule$setProportions()
  cols = get.colors()

  df = data.frame(xmin=rep(0,nrect),xmax=rep(0,nrect),
            ymin=rep(0,nrect),ymax=rep(0,nrect),
            colors=cols[1:nrect])
  for (ifac in 1:nrect)
  {
    rect = rule$getThresh(ifac)
    df$xmin[ifac] = rect[1]
    df$xmax[ifac] = rect[2]
    df$ymin[ifac] = rect[3]
    df$ymax[ifac] = rect[4]
  }
  
  p <- getFigure(padd)
     
  p <- p + ggplot2::geom_rect(data = df, ggplot2::aes(xmin = xmin, xmax = xmax, 
                              ymin = ymin, ymax = ymax, fill = colors))
  
  p <- decor(p, xlab = xlab, ylab = ylab, title = title)
  
  end.func(p, end.plot)
}
 
 
# Function to display a polygon (not tested)
plot.mesh <- function(mesh, 
                      flagEdge=TRUE, flagFace=FALSE, flagApex=FALSE, asp=1,
                      xlim="", ylim="", facecolor="yellow", edgecolor="blue", linewidth=1,
                      show.legend = FALSE, xlab="", ylab="", title="", 
                      padd = NULL, end.plot=TRUE)
{
  p <- getFigure(padd)
  
  if (! flagFace) facecolor = "white"
  if (! flagEdge) edgecolor = facecolor
  
  nmesh = mesh$getNMeshes()
  for (imesh in 1:nmesh)
  {
    xtab = mesh$getCoordinatesPerMesh(imesh-1, 0, TRUE)
    ytab = mesh$getCoordinatesPerMesh(imesh-1, 1, TRUE)
    rp = data.frame(xtab, ytab)
    p <- p + ggplot2::geom_polygon(data = rp, ggplot2::aes(x=xtab,y=ytab), 
                                   linewidth=linewidth, 
                                   fill=facecolor, 
                                   color=edgecolor, show.legend=show.legend)
  if (flagApex)
    p <- p + ggplot2::geom_point(data = rp, ggplot2::aes(x=xtab, y=ytab))
  }  
  
  p <- decor(p, xlab = xlab, ylab = ylab, asp=asp, title = title)
  
  end.func(p, end.plot)
}
setMethod("plot", signature(x="_p_AMesh"), function(x,y=missing,...) plot.mesh(x,...))
 