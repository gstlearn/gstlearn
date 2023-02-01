#Set of global values

plot.initialize <- function() 
{
	plot.default_size <<- list(c(8,8), c(8,8))
	plot.default_xlim <<- c( NA, NA )
	plot.default_ylim <<- c( NA, NA )
	plot.default_sameLim <<- c( FALSE, FALSE )
	plot.default_asp <<- c(NA, 1 )
	invisible()
}

isNotDef <- function(arg)
{
	if (length(arg) == 1)
	{
		if (is.na(arg)) return (TRUE)
	}
	else
	{
		for (i in 1:length(arg))
		{
			if (is.na(arg[i])) return (TRUE)
		}
	}
	return (FALSE)
}

plot.setDefault <- function(icas=1, size=NA, xlim=NA, ylim=NA, sameLim=NA, asp=NA)
{
    if (!isNotDef(size))
        plot.default_size[icas] <<- size
    if (!isNotDef(xlim))
        plot.default_xlim[icas] <<- xlim
    if (!isNotDef(ylim))
        plot.default_ylim[icas] <<- ylim
    if (!isNotDef(sameLim))
        plot.default_sameLim[icas] <<- sameLim
    if (!isNotDef(asp))
        plot.default_asp[icas] <<- asp
}

plot.printDefault <- function()
{
    for (icas in 1:2)
    {
        if (icas == 1)
            cat("Non geographical defaults:\n")
        else
            cat("Geographical defaults:\n")
            
        if (!isNotDef(plot.default_size[[icas]]))
             cat("- Figure size =", plot.default_size[[icas]],"\n")
        if (!isNotDef(plot.default_xlim[icas]))
            cat("- Limits along X =",plot.default_xlim[icas],"\n")
        if (!isNotDef(plot.default_ylim[icas]))
            cat("- Limits along Y =",plot.default_ylim[icas],"\n")
        if (!isNotDef(plot.default_sameLim[icas]))
        {
        	if (plot.default_sameLim[icas])
	            cat("- Limits are the same on both axes\n")
	    }
        if (!isNotDef(plot.default_asp[icas]))
            cat("- Aspect =",plot.default_asp[icas],"\n")
 	}    
}

get.colors <- function()
{
  c("blue", "red", "green", "brown", "orange", "purple", "yellow")
}

getFigure <- function(padd = NULL)
{
  if (length(padd) > 0)
    p <- padd
  else
    p <- ggplot()
  p
}

is_array <- function(arg, ndim=NA)
{
	if (length(arg) <= 1) return (FALSE)
	
	if (!isNotDef(ndim) && length(arg) != ndim) return (FALSE)
	
	TRUE
}

plot.geometry <- function(ax, icas=1, size=NA, xlim=NA, ylim=NA, asp=NA, sameLim=NA)
{
    if (isNotDef(size[1]))
    	size = plot.default_size[[icas]]
    if (! isNotDef(size[1]))
    {
        if (is_array(size, 2))
        {
            if (is_array(ax, 2))
            {
                for (ix in 1:dim(ax)[1])
                    for (iy in 1:dim(ax)[2])
                        ax[ix,iy] <- ax[ix,iy] + ggsave(ax, width=size[1], height=size[2])
            }
            else
            {
            	options(repr.ax.width  = size[1], repr.ax.height = size[2])
            }
        }
        else
            cat("'size' should be [a,b]. Ignored\n")
    }
    if (isNotDef(sameLim))
        sameLim = plot.default_sameLim[icas]
        
    if (isNotDef(xlim[1]))
        xlim = plot.default_xlim[icas]
    if (!isNotDef(xlim[1]))
    {
        if (is_array(xlim, 2))
        {
            if (is_array(ax, 2))
            {
                for (ix in 1:dim(ax)[1])
                    for (iy in 1:dim(ax)[2])
                        ax[ix,iy] <- ax[ix,iy] + xlim(xlim)
 			}
            else
            {
                ax <- ax + xlim(xlim)
            }
		}
        else
            cat("'xlim' should be [a,b] or [None,b] or [a,None]. Ignored\n")
    }
    
    if (isNotDef(ylim[1]))
        ylim = plot.default_ylim[icas]
    if (!isNotDef(ylim[1]))
    {
        if (is_array(ylim, 2))
        {
            if (is_array(ax, 2))
            {
               for (ix in 1:dim(ax)[1])
                    for (iy in 1:dim(ax)[2])
                        ax[ix,iy] <- ax[ix,iy] + ylim(ylim)
 			}
            else
                ax <- ax + ylim(ylim)
        }
        else
           cat("'ylim' should be [a,b] or [None,b] or [a,None]. Ignored\n")
    }
        
    if (isNotDef(asp))
        asp = plot.default_asp[icas]
    if (!isNotDef(asp))
    {
        if (is_array(ax, 2))
        {
        	for (ix in 1:dim(ax)[1])
            	for (iy in 1:dim(ax)[2])
                    ax[ix,iy] <- ax[ix,iy] + coor_fixed(asp)
        }
        else
            ax = ax + coord_fixed(asp)
	}
	ax
}

plot.decoration <- function(p, xlab = NA, ylab = NA, title = NA)
{
  if (!isNotDef(xlab))
    p <- p + labs(x = xlab)
  if (!isNotDef(ylab))
    p <- p + labs(y = ylab)
  if (!isNotDef(title))
    p <- p + ggtitle(title) #+ theme(plot.title = element_text(hjust = 0.5))

  p
}

# Function for representing a Model
plot.model <- function(model, vario=NULL, hmax=1, codir=NULL, 
					   ivar=0, jvar=0, idir=0, asCov=FALSE, nh=100, padd=NULL)
{
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
  
  p  
}
setMethod("plot", signature(x="_p_Model"), function(x,y="missing",...) plot.model(x,...))

# Function for representing the Experimental Variogram together with the Model (optional)

plot.varmod <- function(vario, model=NULL, ivar=-1, jvar=-1, idir=-1,
                        asCov=FALSE, nh=100, draw_psize=FALSE, draw_plabels=FALSE, 
                        color_psize="black", ratio_psize=10,
                        color_plabel="black", size_plabel=2, nudge_y=0.1,
                        show.legend=FALSE, ...)
{
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
        g = ggplot()
        
        if (iv < jv) 
        {
          plot_lst[[index]] <- g
          next
        }
        
        xvlim = c(0,0)
        yvlim = c(0,0)
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
          if (hmax > xvlim[2]) xvlim[2] = hmax 
          gmin = 0
          if (iv != jv) gmin = -gmax
          if (gmin < yvlim[1]) yvlim[1] = gmin
          if (gmax > yvlim[2]) yvlim[2] = gmax
                
          # Plotting the experimental variogram
          df = data.frame(cbind(hh,gg))
          g <- g + geom_line(data = df, aes(x=hh,y=gg), color=cols[id+1], na.rm=TRUE)
          
          if (draw_psize)
            g <- g + geom_point(data = df, aes(x=hh, y=gg), 
                                         size=sw/ratio_psize, color=color_psize, 
                                         na.rm=TRUE, show.legend=show.legend)
            
          if (draw_plabels)
             g <- g + geom_text(data = df, aes(x=hh, y=gg, label=as.character(sw)),
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
            g <- g + geom_line(data = dfg, aes(x=hh,y=gg), color=cols[id+1], size=1, na.rm=TRUE) 

            if (iv != jv)
            {
              ggp = model$sample(hmax, nhh, iv, jv, codir, 1)
              dfg = data.frame(cbind(hh,ggp))
              g <- g + geom_line(data = dfg, aes(x=hh,y=ggp), color=cols[id+1],
                                          linetype = 'twodash', na.rm=TRUE) 
              ggm = model$sample(hmax, nhh, iv, jv, codir,-1)
              dfm = data.frame(cbind(hh,ggp))
              g <- g + geom_line(data = dfm, aes(x=hh,y=ggm), color=cols[id+1],
                                          linetype = 'twodash', na.rm=TRUE) 
            }
          } 
        } # End of loop on Directions
        
        g <- g + scale_x_continuous("Distance", limits = xvlim, expand = c(0,0))
        if (iv == jv)
          g <- g + scale_y_continuous("Variogram", limits = yvlim, expand = c(0,0))
        else
          g <- g + scale_y_continuous("Cross-Variogram", limits = yvlim, expand = c(0,0))
                
        # Plotting relevant control lines
        if (iv != jv)
          g <- g + geom_hline(yintercept = 0.)
        g <- g + geom_hline(yintercept = sill, linetype = 'longdash')
        plot_lst[[index]] <- g
      }
  p = ggarrange(plotlist=plot_lst, nrow=ivarN, ncol = jvarN)
  
  p
}
setMethod("plot", signature(x="_p_Vario"), function(x,y,...) plot.varmod(x,...))

read.pointCoor <- function(db)
{
  xtab = db$getCoordinates(0,TRUE)
  ytab = db$getCoordinates(1,TRUE)
  df = data.frame(xtab,ytab)
  df
}

# Function for plotting a point data base, with optional color and size variables
plot.pointSymbol <- function(p, db, color_name=NULL, size_name=NULL,
                      		 color0='red', size0=0.2, 
                      		 sizmin=10, sizmax=100, flagAbsSize = FALSE, 
                      		 show.legend.color=FALSE, legend.name.color = "P-Color", 
                      		 show.legend.size=FALSE,  legend.name.size="P-Size", 
                      		 ...) 
{  
  # Creating the necessary data frame
  df = read.pointCoor(db)
  np = dim(df)[2]
    
  # Color of symbol
  if (! is.null(color_name)) {
    colval  = db$getColumn(color_name)
  } else {
    colval = rep(color0,np)
  }


  # Size of symbol
  if (! is.null(size_name)) {
  	reduction = 100
    sizval  = db$getColumn(size_name)
    if (flagAbsSize) sizval = abs(sizval)
    m = min(sizval,na.rm=TRUE)
    M = max(sizval,na.rm=TRUE)
    sizval = (sizmax * (sizval - m) / (M-m) + sizmin) / reduction
  } else {
    sizval = rep(size0,np)
  }

  p <- p + geom_point(data = df, aes(x=xtab,y=ytab), color=colval, size=sizval,
                      na.rm=TRUE)
  
  if (show.legend.color && ! is.null(color_name)) {
  	p <- p + guides(color = guide_legend(title = legend.name.color))
  } else {
 	p <- p + guides(color = "none")
  }
    
  if (show.legend.size && ! is.null(size_name)) {
	p <- p + guides(size = guide_legend(title = legend.name.size))
  } else {
   	p <- p + guides(size = "none")
  }
  
  p
}

# Function for plotting a point data base, with optional color and size variables
plot.pointLabel <- function(p, db, label_name=NULL,
                  	  	    color="black", nudge_y=0.1, label_round=2,
                   	 		show.legend=FALSE, legend.name.label = "P-Label", ...) 
{  
  # Creating the necessary data frame
  df = read.pointCoor(db)
  np = dim(df)[2]
    
  # Label of symbols
  labval  = round(db$getColumn(label_name,TRUE),label_round)
  
  p <- p + geom_text(data = df, aes(x=xtab, y=ytab), label=as.character(labval),
            		   nudge_y = nudge_y, color=color, check_overlap=TRUE)
               
  if (show.legend) {
      p <- p + guides(label = guide_legend(title = legend.name.label))
  } else {
      p <- p + guides(label = "none")
  }
  p
}

# Function for plotting a point data base, with optional color and size variables
plot.point <- function(db, color_name=NULL, size_name=NULL, label_name=NULL,
                       color0='red', 
                       size0=0.2, sizmin=10, sizmax=100, flagAbsSize = FALSE,  
                       color_label="black", nudge_y=0.1, label_round=2,
                       show.legend.color=FALSE, legend.name.color="P-Color",
                       show.legend.size=FALSE,  legend.name.size="P-Size",
                       show.legend.label=FALSE, legend.name.label="P-Label",
                       padd = NULL, ...) 
{ 
  p <- getFigure(padd)
  
  if (! is.null(color_name) || ! is.null(size_name))
	  p <- p + plot.pointSymbol(p, db, color_name=color_name, size_name=size_name,
                      	  		color0=color0, size0=size0, 
                      			sizmin=sizmin, sizmax=sizmax, flagAbsSize = flagAbsSize, 
                      			show.legend.symbol=show.legend.symbol, 
                      			legend.name.color = legend.name.color,
                      			legend.name.size = legend.name.size, ...)
  
  if (! is.null(label_name)) 
  	p <- p + plot.pointLabel(p, db, label_name=label_name,
                    	     color=color_label, nudge_y=nudge_y, label_round=label_round,
                      	   	 show.legend=show.legend.label, ...) 
  p
}

#
# Function for plotting a variable (referred by its name) informed in a grid Db
#
# option Indicates the color map (from "A", "B", "C", "D", "E", "F", "G", "H")
plot.grid <- function(dbgrid, name=NULL, na.color = "white", 
      option="B", zlim = NULL, useSel = TRUE,
      show.legend=TRUE, legend.name="", padd=NULL)
{
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
    p <- p + geom_polygon(data = df, aes(x = x, y = y, fill = value, group = id))
  }
  
  # Define the color scale
  p = p + scale_fill_viridis_c(option = option, na.value = na.color, limits=zlim)
  
  if (show.legend) {
    p <- p + guides(fill = guide_colorbar(title=legend.name, reverse=FALSE))
  } else {
    p <- p + guides(fill = "none")
  }
       
  p
}

plot.db <- function(db, padd=NULL, ...)
{
  if (db$isGrid())
    p = plot.grid(db, padd=padd, ...)
  else
    p = plot.point(db, padd=padd, ...)
  p
}

setMethod("plot", signature(x="_p_Db"), function(x,padd=NULL,...) plot.db(x,padd,...))

# Function to display a polygon (not tested)

plot.polygon <- function(poly, color="black", fill=NA, padd = NULL)
{
  npol = poly$getPolySetNumber()
  
  p <- getFigure(padd)
  
  for (ipol in 1:npol)
  {
    xtab = poly$getX(ipol-1)
    ytab = poly$getY(ipol-1)
    rp = data.frame(xtab, ytab)
    p <- p + geom_polygon(data = rp, aes(x=xtab,y=ytab), color=color, fill=fill)
  }  
  
  p
}
setMethod("plot", signature(x="_p_Polygons"), function(x,y=missing,...) plot.polygon(x,...))
        
# Function for plotting the histogram of a variable
plot.hist <- function(db, name, nbins=30, col='grey', fill='yellow', padd = NULL)
{
  val  = db$getColumn(name)
  rp = data.frame(val)
    
  p <- getFigure(padd)
     
  p <- p + geom_histogram(data=rp, aes(x=val), bins=nbins, color=col, fill=fill,
                                   na.rm=TRUE) 
  
  p
}

# Function for plotting histogram for a table of values
plot.hist_tab <- function(val, nbins=30, padd=FALSE)
{
  rp = data.frame(val)
  
  p <- getFigure(padd)
     
  p <- p + geom_histogram(data = rp, aes(x=val), bins=nbins, color='grey', fill='yellow') 

  p
}

# Function for plotting a curve of regularly sampled values
plot.curve <- function(data, color="black", padd=NULL)
{
  nbpoint = length(data)
  absc = seq(1,nbpoint)
  rp = data.frame(absc,data)
  
  p <- getFigure(padd)
    
  p <- p + geom_line(data = rp, aes(x=absc,y=data), color=color, na.rm=TRUE)
  
  p
}

# Function for representing a line between points provided as arguments
plot.XY <-function(xtab, ytab, join=TRUE,
                   color="black", linetype="solid", shape=20,
                   flagDiag = FALSE, 
                   diag_color = "red", diag_line = "solid", padd=NULL)
{
  if (length(ytab) != length(xtab))
  {
    cat("Arrays 'xtab' and 'ytab' should have same dimensions")
    stop()
  }
  rp = data.frame(xtab, ytab)

  p <- getFigure(padd)
     
  if (flagDiag)
  {
    u = min(xtab, ytab, na.rm=TRUE)
    v = max(xtab, ytab, na.rm=TRUE)
    p <- p + geom_segment(aes(x=u,y=u,xend=v,yend=v),
                          linetype = diag_line, color = diag_color, na.rm=TRUE)
  }
  
  if (join)
    p <- p + geom_line(data = rp, aes(x=xtab,y=ytab), 
   		               linetype = linetype, color=color, na.rm=TRUE)
  else 
    p <- p + geom_point(data = rp, aes(x=xtab,y=ytab), 
                        shape=shape, color=color, na.rm=TRUE)
  
  p
}

# Function for representing an anamorphosis
plot.anam <- function(anam, ndisc=100, aymin=-10, aymax=10, 
                      color="black", linetype="solid", padd=NULL)
{
  res = anam$sample(ndisc, aymin, aymax)
  
  p = plot.XY(res$getY(), res$getZ(), join=TRUE, flagDiag = FALSE,
              color=color, linetype=linetype, padd=padd)
              
  p <- p + plot.geometry(p, xlim=res$getAylim(), ylim=res$getAzlim())
  
  p <- p + plot.decoration(p, xlab = "X", ylab = "Z")
  
  p
}

# Function for representing a scatter plot
plot.correlation <- function(db1, name1, name2, db2=NULL, useSel=FALSE,
							 flagDiag = FALSE,
                             color="black", linetype = "solid",
                             diag_color = "red", diag_line = "solid", padd=NULL)
{
  if (is.null(db2)) db2 = db1
  val1 = db1$getColumn(name1, useSel)
  val2 = db2$getColumn(name2, useSel)
  p = plot.XY(val1, val2, join=FALSE, flagDiag=flagDiag, 
              color = color, linetype = linetype, 
              diag_color = diag_color, diag_line = diag_line, padd=padd)
  p 
}

# Representing a Lithotype rule
plot.rule <- function(rule, proportions=NULL, padd=NULL)
{
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
     
  p <- p + geom_rect(data = df, aes(xmin = xmin, xmax = xmax, 
                              ymin = ymin, ymax = ymax, fill = colors))
  
  p
}
 
 
# Function to display a polygon (not tested)
plot.mesh <- function(mesh, 
                      flagEdge=TRUE, flagFace=FALSE, flagApex=FALSE, 
                      facecolor="yellow", edgecolor="blue", linewidth=1,
                      show.legend = FALSE, padd = NULL)
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
    p <- p + geom_polygon(data = rp, aes(x=xtab,y=ytab), 
                                   linewidth=linewidth, 
                                   fill=facecolor, 
                                   color=edgecolor, show.legend=show.legend)
  if (flagApex)
    p <- p + geom_point(data = rp, aes(x=xtab, y=ytab))
  }  
  
  p
}
setMethod("plot", signature(x="_p_AMesh"), function(x,y=missing,...) plot.mesh(x,...))

#setMethod('decoration', 'ggplot', plot.decoration)
