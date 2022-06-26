library(ggplot2)
library(gridExtra)

get.colors <- function()
{
  c("blue", "red", "green", "brown", "orange", "purple", "yellow")
}

decor <- function(p, xlab = xlab, ylab = ylab, title = title)
{
  if (xlab != "")
	  p <- p + labs(x = xlab)
  if (ylab != "")
	  p <- p + labs(y = ylab)
  if (title != "")
	  p <- p + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  p
}

# Function for representing a Model

plot.model <- function(model, hmax, codir=NA, ivar=0, jvar=0, title="", nh=100, padd=NULL)
{
  if (is.na(codir))
  {
    ndim = model$getDimensionNumber()
    codir = rep(0,ndim)
    codir[1] = 1
  }

  hh = seq(0, hmax, hmax/nh)
  gg = model$sampleModel(hmax, nh, ivar, jvar, codir)
  df = data.frame(cbind(hh,gg))
  
  if (length(padd) > 0)
    p = padd
  else
    p <- ggplot()
  
  plot(hh, gg, type="l")
  p <- p + geom_line(data = df, aes(x=hh,y=gg)) + ggtitle(title)
  p
}

# Function for representing the Experimental Variogram together with the Model (optional)

plot.varmod <- function(vario, model=NULL, ivar=-1, jvar=-1, idir=-1,
                        nh=100, title="", ...)
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
          sw = vario$getSwVec(iv,jv,id)
          gg = vario$getGgVec(iv,jv,id)
          hh = vario$getHhVec(iv,jv,id)
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
          g <- g + geom_line(data = df, aes(x=hh,y=gg), color=cols[id+1]) 
 
          # Plotting the Model (optional)
          if (! is.null(model))
          {
            hh = seq(0, hmax, hmax/nh)
            nhh = length(hh)
            codir = vario$getCodir(id)
            gg = model$sampleModel(hmax, nhh, iv, jv, codir)
            dfg = data.frame(cbind(hh,gg))
            g <- g + geom_line(data = dfg, aes(x=hh,y=gg), color=cols[id+1], size=1) 

            if (iv != jv)
            {
              ggp = model$sampleModel(hmax, nhh, iv, jv, codir, 1)
              dfg = data.frame(cbind(hh,ggp))
              g <- g + geom_line(data = dfg, aes(x=hh,y=ggp), color=cols[id+1],
                        linetype = 'twodash') 
              ggm = model$sampleModel(hmax, nhh, iv, jv, codir,-1)
              dfm = data.frame(cbind(hh,ggp))
              g <- g + geom_line(data = dfm, aes(x=hh,y=ggm), color=cols[id+1],
                        linetype = 'twodash') 
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
  all.plots <- marrangeGrob(plot_lst, nrow = ivarN, ncol = jvarN, top = title)
  all.plots
}

# Function for plotting a point data base, with optional color and size variables

plot.point <- function(db, color_name=NA, size_name=NA,
              col0='red', cex0=0.2, sizmin=10, sizmax=200, asp=1, pch0=19, 
              xlab="", ylab="", title="", padd = NULL, ...) 
{    
  # Extracting coordinates
  tabx = db$getCoordinates(0,TRUE)
  taby = db$getCoordinates(1,TRUE)
  np   = length(tabx)
    
  # Color of symbol
  if (! is.na(color_name))
  {
    colval  = Db_getColumn(db,color_name,TRUE)
    colval[colval == getTEST()] = NA
  }
  else
  {
    colval = rep(col0,np)
  }

  # Size of symbol
  reduction = 100
  if (! is.na(size_name))
  {
    sizval  = Db_getColumn(db,size_name,TRUE)
    sizval[sizval == getTEST()] = NA
    m = min(abs(sizval))
    M = max(abs(sizval))
    sizval = (sizmax * (abs(sizval) - m) / (M-m) + sizmin) / reduction
  }
  else
  {
    sizval = rep(cex0,np)
  }

  df = data.frame(tabx,taby,colval,sizval)
  
  if (length(padd) > 0)
    p = padd
  else
    p <- ggplot()
  p <- p + geom_point(data=df,aes(x=tabx,y=taby,size=sizval,color=colval))

  p <- p + theme(legend.position = "none")
  
  p <- decor(p, xlab = xlab, ylab = ylab, title = title)
  p
}

# Function to display a polygon (not tested)

plot.polygon <- function(poly, title="")
{    
  npol = poly$getPolySetNumber()
  cols = get.colors()
  
  ids = seq(1,npol)
  values = data.frame(
    id = ids,
    value = cols[ids]
  )
    
  for (ipol in 1:npol)
  {
    x = poly$getX(ipol)
    y = poly$getY(ipol)
    plt.fill(x, y, color)
  }  
  ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
}
        
# Function for plotting a variable (referred by its name) informed in a grid Db

plot.grid <- function(dbgrid, name, xlab="", ylab="", title = "", padd=NULL)
{
  if (! dbgrid$isGrid())
  {
    cat("This function is restricted to Grid Db and cannot be used here")
    return
  }

  x = dbgrid$getColumnByLocator(ELoc_X(),0)
  y = dbgrid$getColumnByLocator(ELoc_X(),1)
  data = Db_getColumn(dbgrid, name)
  data[data == getTEST()] = NA
  df = data.frame(x,y,data)

  if (length(padd) > 0)
   p <- padd
  else
   p <- ggplot()
  p <- p + geom_raster(data = df, aes(x = x, y = y, fill = data)) + 
       scale_fill_viridis_c(option = "inferno", na.value = 'white')
       
  p <- decor(p, xlab = xlab, ylab = ylab, title = title)
  p
}

# Function for plotting the histogram of a variable

plot.hist <- function(db, name, nbins=30, col='grey', fill='yellow',
            xlab="", ylab="", title="", padd = NULL)
{    
  val  = Db_getColumn(db,name)
  rp = data.frame(val)
    
  if (length(padd) > 0)
    p <- padd
  else
    p <- ggplot() 
  p <- p + geom_histogram(data=rp, aes(x=val), bins=nbins, color=col, fill=fill) 
  p <- decor(p, xlab = xlab, ylab = ylab, title = title)
  p
}

# Function for plotting histogram for a table of values
plot.hist_tab <- function(val, nbins=30, xlab="", ylab="", title="", padd=FALSE)
{
  rp = data.frame(val)
  
  if (length(padd) > 0)
    p <- padd
  else
    p <- ggplot() 
  p <- p + geom_histogram(data = rp, aes(x=val), bins=nbins, color='grey', fill='yellow') 
  p <- decor(p, xlab = xlab, ylab = ylab, title = title)
  p
}

# Function for plotting a curve of regularly sampled values
plot.curve <- function(data, color="black",
    xlab="", ylab="", title="", padd=NULL)
{
  nbpoint = length(data)
  absc = seq(1,nbpoint)
  rp = data.frame(absc,data)
  
  if (length(padd) > 0)
    p <- padd
  else
    p <- ggplot() 
  p <- p + geom_line(data = rp, aes(x=absc,y=data), color=color)
  p <- decor(p, xlab = xlab, ylab = ylab, title = title)
  p
}

# Function for representing a line between points provided as arguments
plot.XY <-function(xtab, ytab, join=TRUE,
	color="black", linetype="solid", shape=20,
	flagDiag = FALSE, diag_color = "red", diag_line = "solid",
	xlim="", ylim="", xlab="", ylab="", title="", padd=NULL)
{    
  if (length(ytab) != length(xtab))
  {
    cat("Arrays 'xtab' and 'ytab' should have same dimensions")
    return
  }
  rp = data.frame(xtab, ytab)

  if (length(padd) > 0)
    p <- padd
  else
    p <- ggplot() 
    
  if (is.numeric(xlim) && length(xlim) == 2)
    p <- p + scale_x_continuous(limits = xlim, expand = c(0,0))
  if (is.numeric(ylim) && length(ylim) == 2)
    p <- p + scale_y_continuous(limits = ylim, expand = c(0,0))
  
  if (flagDiag)
  {
  	u = min(xtab, ytab)
  	v = max(xtab, ytab)
  	p <- p + geom_segment(aes(x=u,y=u,xend=v,yend=v),
  		linetype = diag_line, color = diag_color)
  }
  
  if (join)
 	 p <- p + geom_line(data = rp, aes(x=xtab,y=ytab), 
 	 	linetype = linetype, color=color)
  else 
    p <- p + geom_point(data = rp, aes(x=xtab,y=ytab), 
    	shape=shape, color=color)
  
  p <- decor(p, xlab = xlab, ylab = ylab, title = title)
  p
}

plot.anam <- function(anam, ndisc=100, aymin=-10, aymax=10, 
    color="black", linetype="solid",
    xlim="", ylim="", xlab="Y", ylab="Z", title="", padd=NULL)
{
  res = anam$sample(ndisc, aymin, aymax)
  valY = res$getY()
  valZ = res$getZ()
  p = plot.XY(valY, valZ, join=TRUE, flagDiag = FALSE,
  		color=color, linetype=linetype, 
	    xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, title=title, 
        padd=padd)
  p
}

plot.correlation <- function(db1, name1, name2, db2, flagDiag = FALSE,
  color="black", linetype = "solid",
  diag_color = "red", diag_line = "dashed",
  xlim="", ylim="", xlab="", ylab="", title="", 
  padd=NULL)
{
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
plot.rule <- function(rule, proportions=NA, xlab="", ylab="", title="", padd=NULL)
{
  nrect = rule$getFaciesNumber()
  if (! is.na(proportions)) 
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
  
  if (length(padd) > 0)
    p <- padd
  else
    p <- ggplot() 
  p <- p + geom_rect(data = df, aes(xmin = xmin, xmax = xmax, 
                                    ymin = ymin, ymax = ymax, fill = colors))
  
  p <- decor(p, xlab = xlab, ylab = ylab, title = title)
  p
}
 