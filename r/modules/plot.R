library(ggplot2)
library(ggpubr)

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

decor <- function(p, xlab = "", ylab = "", asp = NULL, title = "")
{
  if (xlab != "")
	  p <- p + labs(x = xlab)
  if (ylab != "")
	  p <- p + labs(y = ylab)
  if (title != "")
	  p <- p + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  if (! is.null(asp))
  	  suppressWarnings(suppressMessages(p <- p + coord_fixed(ratio = asp)))
  p
}

# Function for representing a Model
plot.model <- function(model, hmax, codir=NULL, ivar=0, jvar=0, 
                       title="", nh=100, padd=NULL)
{
  if (is.null(codir))
  {
    ndim = model$getDimensionNumber()
    codir = rep(0,ndim)
    codir[1] = 1
  }

  p <- getFigure(padd)
  
  hh = seq(from=0, to=hmax, length.out=nh)
  gg = model$sample(hmax, nh, ivar, jvar, codir)
  df = data.frame(cbind(hh,gg))
  
  plot(hh, gg, type="l")
  p <- p + geom_line(data = df, aes(x=hh,y=gg), na.rm=TRUE)
  
  p <- decor(p, xlab = xlab, ylab = ylab, asp=as, title = title)
  
  p
}
setMethod("plot", signature(x="_p_Model"), function(x,y="missing",...) plot.model(x,...))


# Function for representing the Experimental Variogram together with the Model (optional)

plot.varmod <- function(vario, model=NULL, ivar=-1, jvar=-1, idir=-1,
                        nh=100, draw_psize=FALSE, draw_plabels=FALSE, 
                        color_psize="black", ratio_psize=3,
                        color_plabel="black", size_plabel=2, nudge_y=0.1,
                        title="", ...)
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
        
        xlim = c(0,0)
        ylim = c(0,0)
        for (id in idirUtil)
        {
          sill = vario$getVar(iv,jv)
          nlag = vario$getLagNumber(id)
          sw = vario$getSwVec(id,iv,jv)
          gg = vario$getGgVec(id,iv,jv)
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
          g <- g + geom_line(data = df, aes(x=hh,y=gg), color=cols[id+1], na.rm=TRUE)
          
          if (draw_psize)
         	 g <- g + geom_point(data = df, aes(x=hh, y=gg), 
         	 	size=sw/ratio_psize, color=color_psize, 
         	 	na.rm=TRUE, show.legend=FALSE)
         	 
          if (draw_plabels)
          	 g <- g + geom_text(data = df, aes(x=hh, y=gg, label=as.character(sw)),
          	 	color=color_plabel, size=size_plabel, nudge_y=nudge_y, 
          	 	show.legend=FALSE, check_overlap=TRUE)
 	
          # Plotting the Model (optional)
          if (! is.null(model))
          {
            hh = seq(0, hmax, hmax/nh)
            nhh = length(hh)
            codir = vario$getCodir(id)
            gg = model$sample(hmax, nhh, iv, jv, codir)
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
        
        g <- g + scale_x_continuous("Distance", limits = xlim, expand = c(0,0))
        if (iv == jv)
          g <- g + scale_y_continuous("Variogram", limits = ylim, expand = c(0,0))
        else
          g <- g + scale_y_continuous("Cross-Variogram", limits = ylim, expand = c(0,0))
                
        # Plotting relevant control lines
        if (iv != jv)
          g <- g + geom_hline(yintercept = 0.)
        g <- g + geom_hline(yintercept = sill, linetype = 'longdash')
        plot_lst[[index]] <- g
      }
	p = ggarrange(plotlist=plot_lst, nrow=ivarN, ncol = jvarN)
	
	p <- decor(p, title = title)
	
	p
}
setMethod("plot", signature(x="_p_Vario"), function(x,y,...) plot.varmod(x,...))

# Function for plotting a point data base, with optional color and size variables
plot.point <- function(db, color_name=NULL, size_name=NULL, label_name=NULL,
              color0='red', size0=0.2, color_label="black", nudge_y=0.1,
              sizmin=10, sizmax=100, flagAbsSize = FALSE, 
              show.legend.color=FALSE, name.legend.color="P-Color",
              show.legend.size =FALSE, name.legend.size ="P-Size",
              show.legend.label=FALSE, name.legend.label="P-Label",
              asp=1, xlab="", ylab="", title="", padd = NULL, ...) 
{  
  # Creating the necessary data frame
  np   = db$getSampleNumber(TRUE)
  tabx = db$getCoordinates(0,TRUE)
  taby = db$getCoordinates(1,TRUE)
    
  # Color of symbol
  if (! is.null(color_name))
  {
    colval  = Db_getColumn(db,color_name,TRUE)
  }
  else
  {
    colval = rep(color0,np)
  }

  # Size of symbol
  reduction = 100
  if (! is.null(size_name))
  {
    sizval  = Db_getColumn(db,size_name,TRUE)
    if (flagAbsSize) sizval = abs(sizval)
    m = min(sizval,na.rm=TRUE)
    M = max(sizval,na.rm=TRUE)
    sizval = (sizmax * (sizval - m) / (M-m) + sizmin) / reduction
  }
  else
  {
    sizval = rep(size0,np)
  }

  # Label of sylbols
  if (! is.null(label_name))
  {
 	label_round = 2
    labval  = round(Db_getColumn(db,label_name,TRUE),label_round)
  }
  else
  {
    labval = rep(0,np)
  }
  
  df = data.frame(tabx,taby,colval,sizval,labval)
  
  p <- getFigure(padd)
     
  p <- p + geom_point(data=df, aes(x=tabx,y=taby,color=colval,size=sizval),
  		na.rm=TRUE)
  
  if (! is.null(label_name)) 
  {
 	p <- p + geom_text(data = df, aes(x=x, y=y, label=as.character(labval)),
          	 	nudge_y=nudge_y, color=color_label, check_overlap=TRUE)
	if (show.legend.label) {
	  p <- p + guides(label = guide_legend(title = name.legend.label))
	} else {
	  p <- p + guides(label = "none")
	}
  }
  
  if (show.legend.color) {
  	p <- p + guides(color = guide_legend(title = name.legend.color))
  } else {
    p <- p + guides(color = "none")
  }
  	
  if (show.legend.size) {
  	p <- p + guides(size = guide_legend(title = name.legend.size))
  } else {
    p <- p + guides(size = "none")
  }
  		
  p <- decor(p, xlab = xlab, ylab = ylab, asp = asp, title = title)
  
  p
}

# Function for plotting a variable (referred by its name) informed in a grid Db
plot.grid <- function(dbgrid, name=NULL, color_NA = "white", asp=1,
			show.legend=TRUE, name_legend="G-Fill",
			xlab="", ylab="", title="", padd=NULL)
{
  if (! dbgrid$isGrid())
  {
    cat("This function is restricted to Grid Db and cannot be used here")
    return
  }

  # Building the necessary data frame
  x = dbgrid$getColumnByLocator(ELoc_X(),0)
  y = dbgrid$getColumnByLocator(ELoc_X(),1)
  
  if (is.null(name))
  {
  	if (dbgrid$getLocatorNumber(ELoc_Z()) > 0) 
		name = dbgrid$getNameByLocator(ELoc_Z())
	else
		name = dbgrid$getLastName()
  }
  data = Db_getColumn(dbgrid, name)
  df = data.frame(x,y,data)
  
  p <- getFigure(padd)
   
  p <- p + geom_raster(data = df, aes(x = x, y = y, fill = data), 
  			show.legend=show.legend) + 
       scale_fill_viridis_c(option = "inferno", na.value = color_NA)
       
  if (show.legend) {
  	p <- p + guides(fill = guide_legend(title=name_legend))
  } else {
  	p <- p + guides(fill = "none")
  }
       
  p <- decor(p, xlab = xlab, ylab = ylab, asp=asp, title = title)
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

setMethod('$', '_p_Limits', function(x, name)
setMethod("plot", signature(x="_p_Db"), function(x,padd=NULL,...) plot.db(x,padd,...))

# Function to display a polygon (not tested)
plot.polygon <- function(poly, xlab="", ylab="", title="", padd = NULL)
{    
  npol = poly$getPolySetNumber()
  cols = get.colors()
  
  ids = seq(1,npol)
  values = data.frame(
    id = ids,
    value = cols[ids]
  )
   
  p <- getFigure(padd)
  
  for (ipol in 1:npol)
  {
    x = poly$getX(ipol)
    y = poly$getY(ipol)
    p <- p + plt.fill(x, y, color)
  }  
  
  p <- decor(p, xlab = xlab, ylab = ylab, asp=asp, title = title)
  p
}
setMethod("plot", signature(x="_p_Polygons"), function(x,y=missing,...) plot.polygon(x,...))
        
# Function for plotting the histogram of a variable
plot.hist <- function(db, name, nbins=30, col='grey', fill='yellow',
            xlab="", ylab="", title="", padd = NULL)
{    
  val  = Db_getColumn(db,name)
  rp = data.frame(val)
    
  p <- getFigure(padd)
     
  p <- p + geom_histogram(data=rp, aes(x=val), bins=nbins, color=col, fill=fill,
  						  na.rm=TRUE) 
  
  p <- decor(p, xlab = xlab, ylab = ylab, title = title)
  
  p
}

# Function for plotting histogram for a table of values
plot.hist_tab <- function(val, nbins=30, xlab="", ylab="", title="", padd=FALSE)
{
  rp = data.frame(val)
  
  p <- getFigure(padd)
     
  p <- p + geom_histogram(data = rp, aes(x=val), bins=nbins, color='grey', fill='yellow') 

  p <- decor(p, xlab = xlab, ylab = ylab, title = title)

  p
}

# Function for plotting a curve of regularly sampled values
plot.curve <- function(data, color="black", xlab="", ylab="", title="", padd=NULL)
{
  nbpoint = length(data)
  absc = seq(1,nbpoint)
  rp = data.frame(absc,data)
  
  p <- getFigure(padd)
    
  p <- p + geom_line(data = rp, aes(x=absc,y=data), color=color, na.rm=TRUE)
  
  p <- decor(p, xlab = xlab, ylab = ylab, title = title)
  
  p
}

# Function for representing a line between points provided as arguments
plot.XY <-function(xtab, ytab, join=TRUE,
	               color="black", linetype="solid", shape=20,
	               flagDiag = FALSE, 
	               diag_color = "red", diag_line = "solid",
	               xlim="", ylim="", xlab="", ylab="", title="", padd=NULL)
{    
  if (length(ytab) != length(xtab))
  {
    cat("Arrays 'xtab' and 'ytab' should have same dimensions")
    return
  }
  rp = data.frame(xtab, ytab)

  p <- getFigure(padd)
     
  if (is.numeric(xlim) && length(xlim) == 2)
    p <- p + scale_x_continuous(limits = xlim, expand = c(0,0))
  if (is.numeric(ylim) && length(ylim) == 2)
    p <- p + scale_y_continuous(limits = ylim, expand = c(0,0))
  
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
  
  p <- decor(p, xlab = xlab, ylab = ylab, title = title)
  
  p
}

# Function for representing an anamorphosis
plot.anam <- function(anam, ndisc=100, aymin=-10, aymax=10, 
    color="black", linetype="solid",
    xlim="", ylim="", xlab="Y", ylab="Z", title="", padd=NULL)
{
  res = anam$sample(ndisc, aymin, aymax)
  valY = res$getY()
  valZ = res$getZ()
  
  p = plot.XY(valY, valZ, join=TRUE, flagDiag = FALSE,
  		      color=color, linetype=linetype, 
	          xlim=res$getAylim(), ylim=res$getAzlim(), xlab=xlab, ylab=ylab, title=title, 
              padd=padd)
  p
}

# Function for representing a scatter plot
plot.correlation <- function(db1, name1, name2, db2=NULL, flagDiag = FALSE,
		color="black", linetype = "solid",
 		diag_color = "red", diag_line = "solid",
 		xlim="", ylim="", xlab="", ylab="", title="", 
  		padd=NULL)
{
  if (is.null(db2)) db2 = db1
  val1 = Db_getColumn(db1,name1)
  val2 = Db_getColumn(db2,name2)
  p = plot.XY(val1, val2, join=FALSE, flagDiag=flagDiag, 
  		      color=color, linetype = linetype, 
  		      diag_color = diag_color, diag_line = diag_line,
     	      xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, title=title, 
    	      padd=padd)
  p 
}

# Representing a Lithotype rule
plot.rule <- function(rule, proportions=NULL, xlab="", ylab="", title="", padd=NULL)
{
  nrect = rule$getFaciesNumber()
  if (! is.null(proportions)) 
    rule$setProportions(proportions)
  else
    rule$setProportions()
  cols = get.colors()

  df = data.frame(xmin=rep(0,nrect),xmax=rep(0,nrect),ymin=rep(0,nrect),ymax=rep(0,nrect),
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
  
  p <- decor(p, xlab = xlab, ylab = ylab, title = title)
  
  p
}
 