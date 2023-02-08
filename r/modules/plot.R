#
# This is a set of functions working with ggplot2() which enable performing
# plots easily 
#

# Define the global values
plot.initialize <- function() 
{
	plot.default_dims <<- list(c(8,8), c(8,8))
	plot.default_xlim <<- list(c(NA,NA), c(NA,NA))
	plot.default_ylim <<- list(c(NA,NA), c(NA,NA))
	plot.default_asp  <<- c(0, 1 )
	invisible()
}

isNotDef <- function(arg)
{
	if (is.null(arg)) return (TRUE)
	warn.old = options("warn")
	options(warn = -1)
	
	if (length(arg) == 1)
	{
		if (is.na(arg)) 
		{
			options(warn.old)
			return (TRUE)
		}
	}
	else
	{
		for (i in 1:length(arg))
		{
			if (is.na(arg[i])) 
			{
				options(warn.old)
				return (TRUE)
			}
		}
	}
	options(warn.old)
	return (FALSE)
}

plot.setDefault <- function(mode=1, dims=NA, xlim=NA, ylim=NA, asp=NA)
{
    if (!isNotDef(dims))
        plot.default_dims[[mode]] <<- dims
    if (!isNotDef(xlim))
        plot.default_xlim[[mode]] <<- xlim
    if (!isNotDef(ylim))
        plot.default_ylim[[mode]] <<- ylim    
    if (!isNotDef(asp))
	    plot.default_asp[[mode]] <<- asp
}

plot.printDefault <- function()
{
    for (mode in 1:2)
    {
        if (mode == 1)
            cat("Non geographical defaults (mode=1):\n")
        else
            cat("Geographical defaults (mode=2):\n")
            
        if (!isNotDef(plot.default_dims[[mode]]))
             cat("- Figure dimensions =", plot.default_dims[[mode]],"\n")
        else
        	 cat("- Figure dimensions (not defined)\n")
        	 
        if (!isNotDef(plot.default_xlim[[mode]]))
            cat("- Limits along X =",plot.default_xlim[[mode]],"\n")
        else
        	cat("- Limits along X (not defined)\n")
        	
        if (!isNotDef(plot.default_ylim[[mode]]))
            cat("- Limits along Y =",plot.default_ylim[[mode]],"\n")
        else
        	cat("- Limits along Y (not defined)\n")
        	
        if (plot.default_asp[mode] != 0)
            cat("- Aspect =",plot.default_asp[mode],"\n")
        else
        	cat("- Aspect (automatic)\n")
 	}    
}

get.colors <- function()
{
  c("blue", "red", "green", "brown", "orange", "purple", "yellow")
}

getNewFigure <- function(padd = NULL, mode = 1)
{
  if (length(padd) > 0)
  {
    p <- padd
  }
  else
  {
    p <- ggplot()
    
    p <- plot.geometry(p,
    				   dims=plot.default_dims[[mode]], 
    			  	   xlim=plot.default_xlim[[mode]], 
    				   ylim=plot.default_ylim[[mode]], 
    				   asp=plot.default_asp[mode])
  }
  p
}

is_array <- function(arg, ndim=NA)
{
	if (length(arg) <= 1) return (FALSE)
	
	if (!isNotDef(ndim) && length(arg) != ndim) return (FALSE)
	
	TRUE
}

plot.decoration <- function(p, xlab = NA, ylab = NA, title = NA)
{
  if (! is.ggplot(p))
  {
  	print("This function is (currently) limited to ggplot(s)")
  	return(p)
  }
  if (!isNotDef(xlab))
    p <- p + labs(x = xlab)
  if (!isNotDef(ylab))
    p <- p + labs(y = ylab)
  if (!isNotDef(title))
    p <- p + ggtitle(title) + theme(plot.title = element_text(hjust = 0.5))
  p
}

plot.empty <- function(p)
{
	p = p +
  	theme(
   	 	panel.background = element_rect(fill='transparent'), #transparent panel bg
   	 	plot.background = element_rect(fill='transparent', color=NA), #transparent plot bg
    	panel.grid.major = element_blank(), #remove major gridlines
    	panel.grid.minor = element_blank(), #remove minor gridlines
   		legend.background = element_rect(fill='transparent'), #transparent legend bg
    	legend.box.background = element_rect(fill='transparent') #transparent legend panel
    )
    p
}

# Set the Geometry for the plot 'p'
# asp: Specify a value of "0" for an automatic aspect ratio
#
plot.geometry <- function(p, dims=NA, xlim=NA, ylim=NA, asp=NA, expand=waiver())
{
    if (! isNotDef(dims[1]))
    {
        if (is_array(dims, 2))
        {
            if (is_array(p, 2))
            {
                for (ix in 1:dim(p)[1])
                    for (iy in 1:dim(p)[2])
                        p[ix,iy] <- p[ix,iy] + ggsave(p, width=dims[1], height=dims[2])
            }
            else
            {
            	options(repr.p.width  = dims[1], repr.p.height = dims[2])
            }
        }
        else
        {
            cat("'dims' should be [a,b]. Ignored\n")
        }
    }
    
    if (!isNotDef(xlim[1]))
    {
        if (is_array(xlim, 2))
        {
            if (is_array(p, 2))
            {
                for (ix in 1:dim(p)[1])
                    for (iy in 1:dim(p)[2])
                        p[ix,iy] <- p[ix,iy] + scale_x_continuous(limits=xlim, expand=expand)
 			}
            else
            {
                p <- p + scale_x_continuous(limits=xlim, expand=expand)
            }
		}
        else
            cat("'xlim' should be [a,b] or [NA,b] or [a,NA]. Ignored\n")
    }
    
    if (!isNotDef(ylim[1]))
    {
        if (is_array(ylim, 2))
        {
            if (is_array(p, 2))
            {
               for (ix in 1:dim(p)[1])
                    for (iy in 1:dim(p)[2])
                        p[ix,iy] <- p[ix,iy] + scale_y_continuous(limits=ylim, expand=expand)
 			}
            else
                p <- p + scale_y_continuous(limits=ylim, expand=expand)
        }
        else
           cat("'ylim' should be [a,b] or [NA,b] or [a,NA]. Ignored\n")
    }
    
    if (!isNotDef(asp))
    {     
    	if (is_array(p, 2))
    	{
	   		for (ix in 1:dim(p)[1])
    	   		for (iy in 1:dim(p)[2])
    	   		{
    	   			if (asp != 0)
	       	    		p[ix,iy] <- p[ix,iy] + coor_fixed(asp)
	       	    	else
	       	    		p[ix,iy] <- p[ix,iy] + theme(aspect.ratio = 1)
	       	    }
    	}
    	else
    	{
    		if (asp != 0)
	    		p = p + coord_fixed(asp)
	    	else
	    		p = p + theme(aspect.ratio = 1)
	    }
    }

	p
}

# Function for representing a Model
plot.model <- function(model, ivar=0, jvar=0, codir=NA, vario=NA, idir=0,
                       nh = 100, hmax = NA, asCov=FALSE,
          			   env_color='black', env_linetype='dashed', env_size=0.5,
          			   flag.envelop = TRUE, padd=NULL, ...)
{
    p = getNewFigure(padd, 1)
    
    # Calculate the Direction if not defined (possibly using the Variogram)
    if (isNotDef(codir))
    {
        if (isNotDef(vario))
        {
            codir = rep(0, model$getDimensionNumber())
            codir[0] = 1
        }
        else
        {
            codir = vario$getCodirs(idir)
        }
    }
    
    p = modelElem(p, model, ivar=ivar, jvar=jvar, codir=codir,
                  nh = nh, hmax = hmax, asCov=asCov,  flag.envelop=flag.envelop,
                  env_color = env_color, env_linetype= env_linetype, env_size=env_size,
                  ...)
   p  
}

# Function for representing the Experimental Variogram together with the Model (optional)

plot.vario <- function(vario, ivar=0, jvar=0, idir=0, hmax=NA,
    			       var_color='black', var_linetype='dashed', var_size=0.5,
    			       draw_psize = FALSE, draw_plabel = FALSE, 
          			   show.legend=FALSE, padd=NULL, ...)
{
    p = plot.varmod(vario, ivar=ivar, jvar=jvar, idir=idir, hmax=hmax,
                    var_color=var_color, var_linetype=var_linetype, var_size=var_size,
                    draw_psize=draw_psize, draw_plabel=draw_plabel, 
                    show.legend=show.legend, padd=padd, ...)
    
    p
}

selectItems <- function(nvalues, sitem=-1)
{
    if (sitem >= 0)
        outs = sitem
    else
    	outs = seq(0, nvalues-1)
    outs
}

varioLayer <- function(vario, ivar=0, jvar=0, idir=0, mode=0, ...)
{
     # Plotting the experimental variogram
    
    gg = vario$getGgVec(idir,ivar,jvar) 
    hh = vario$getHhVec(idir,ivar,jvar)
    sw = vario$getSwVec(idir,ivar,jvar)
    df = data.frame(gg = gg, hh = hh, sw = sw)
    
    # Dispatch on the type of representation
    
    if (mode == 0)
    {
    	# Representing the Experimental variogram
    	layer = geom_line(data = df, mapping=aes(x=hh, y=gg),  ...)
    }
    else if (mode == 1) 
	{   
    	# Representing the number of pairs (by size)
	    layer <- geom_point(data = df, mapping=aes(x=hh, y=gg, size=sw), ...)
	}
	else
	{
	    # Representing the number of pairs (by label)
        layer <- geom_text(data = df, mapping=aes(x=hh, y=gg, label=as.character(sw)), ...)
    } 			
    layer
}

varioElem <- function(p, vario, ivar=0, jvar=0, idir=0, 
          			  var_color='black', var_linetype="dashed", var_size=0.5, 
					  draw_variance = TRUE, draw_psize = FALSE, draw_plabel = FALSE, 
					  label=NULL, ...)
{
	p = p + varioLayer(vario, ivar=ivar, jvar=jvar, idir=idir, mode=0, ...)
	
	if (draw_psize)
			p = p + varioLayer(vario, ivar=ivar, jvar=jvar, idir=idir, mode=1, ...)
			                   
	if (draw_plabel)
			p = p + varioLayer(vario, ivar=ivar, jvar=jvar, idir=idir, mode=2, ...)
			                   
	# Constructing the label for Legend
    if (length(label) <= 0)
	   	label = paste("vario dir=", paste(vario$getCodirs(idir), collapse=' '))
    
    # Adding the vertical axis at X=0
    p = p + geom_vline(xintercept = 0., color='black', size=0.5)
                     
    # Adding the horizontal axis at Y=0   	     
    p = p + geom_hline(yintercept = 0., color="black", size=0.5)

   	# Drawing the variance-covariance reference line (optional)
   	if (draw_variance)
        p = p + geom_hline(yintercept=vario$getVar(ivar,jvar), 
        				   color=var_color, linetype=var_linetype, size=var_size)
	
	# Tuning the bounds of graphics
    if (vario$drawOnlyPositiveX(ivar, jvar))
		p = plot.geometry(p, xlim = c(0, NA))
    if (vario$drawOnlyPositiveY(ivar, jvar))
		p = plot.geometry(p, ylim = c(0, NA))
	
	p
}

modelLayer <- function(model, ivar=0, jvar=0, codir=NA,
            	       nh = 100, hmax = NA, asCov=FALSE, nostd=0, ...)
{
    istart = 0
    for (icova in 1:model$getCovaNumber())
    {
        if (model$getCovName(icova-1) == 'Nugget Effect')
            istart = 1 # do not plot the first lag (h=0) for nugget effect (discontinuity)
    }
     
    # Calculating distances 
    hh = seq(0., hmax, length.out=nh+1)
    
    # Representing the Model
	gg = model$sample(hmax=hmax, nh=nh, ivar=ivar, jvar=jvar, codir=codir, 
					  nostd=nostd, asCov=asCov, addZero=TRUE)
    df = data.frame(gg = gg[istart:nh], hh = hh[istart:nh])
    layer = geom_line(data = df, mapping=aes(x=hh, y=gg), na.rm=TRUE, ...)
    
    layer
}

modelElem <- function(p, model, ivar=0, jvar=0, codir=NA, 
            	      nh = 100, hmax = NA, asCov=FALSE, flag.envelop = TRUE, 
            	      env_color='black', env_linetype="dashed", env_size=0.5,
            	      ...)
{
    # if hmax not specified = 3*maximum range of the model's basic structures
    if (isNotDef(hmax))
    {
        hmax = 0
        for (icova in 1:model$getCovaNumber())
        {
            range_max = max(model$getCova(icova-1)$getRanges())
            if (3*range_max > hmax)
                hmax = 3*range_max
        }
    }
    if (hmax == 0) hmax = 1.
            
    istart = 0
    for (icova in 1:model$getCovaNumber())
    {
        if (model$getCovName(icova-1) == 'Nugget Effect')
            istart = 1 # do not plot the first lag (h=0) for nugget effect (discontinuity)
    }
     
    # Represent the Model
    p = p + modelLayer(model, ivar=ivar, jvar=jvar, codir=codir, 
    				   nh=nh, hmax=hmax, asCov=asCov, nostd=0, ...)
    
    # Represent the coregionalization envelop
    if (ivar != jvar && flag.envelop)
    {
    	p = p + modelLayer(model, ivar=ivar, jvar=jvar, codir=codir, 
    				       nh=nh, hmax=hmax, asCov=asCov, nostd=-1, 
    				       color = env_color, linetype = env_linetype, size=env_size)
    	
    	p = p + modelLayer(model, ivar=ivar, jvar=jvar, codir=codir, 
    				       nh=nh, hmax=hmax, asCov=asCov, nostd=+1, 
    				       color = env_color, linetype = env_linetype, size=env_size)
    }
    
    p
}

plot.varmod <- function(vario, model=NA, ivar=-1, jvar=-1, idir=-1,
          				nh = 100, hmax = NA, draw_psize=FALSE, draw_plabel=FALSE, 
          				asCov=FALSE, draw_variance = TRUE, flag.envelop=TRUE, 
          				var_color='black', var_linetype="dashed", var_size=0.5, 
          				env_color='black', env_linetype="dashed", env_size=0.5,
          				show.legend=FALSE, label=NULL, draw.vario=TRUE, padd=NULL, ...)
{
  dots = list(...)
  has_color = "color" %in% names(dots)
  
  ndir = vario$getDirectionNumber()
  nvar = vario$getVariableNumber()
  cols = get.colors()
  
  idirUtil = selectItems(ndir, idir)
  ivarUtil = selectItems(nvar, ivar)
  jvarUtil = selectItems(nvar, jvar)
  ivarN = length(ivarUtil)
  jvarN = length(jvarUtil)

  empty_figure = length(padd) <= 0
  multi_figure = (ivarN * jvarN > 1)

  if (isNotDef(hmax))
	hmax = vario$getHmax(ivar, jvar, idir)

  # Loop on the variables
  
  index = 0
  plot_lst <- vector("list", length = ivarN * jvarN)

  for (ivar in ivarUtil)
    for (jvar in jvarUtil)
      {

        # Define the current plot
        index = index + 1
        if (multi_figure)
        {
        	if (empty_figure)
        		g = getNewFigure(NULL, 1)
        	else
        		g = padd[[index]]
        }
        else
        {
        	g = getNewFigure(padd, 1)
        }
        
        if (ivar < jvar) 
        {
          plot_lst[[index]] <- plot.empty(g)
          next
        }
        
        for (idir in idirUtil)
        {
        	color = NA
       		if (! has_color) color=cols[idir+1]
        	g = varioElem(g, vario, ivar, jvar, idir, 
        			      var_color=var_color, var_linetype=var_linetype, var_size=var_size,
   				          draw_variance=draw_variance, draw_psize=draw_psize, draw_plabel=draw_plabel,
        			      label=label, color=color, ...)

            # Plotting the Model (optional)
            if (! isNotDef(model))
            {
            	codir = vario$getCodirs(idir)
            	g = modelElem(g, model, ivar, jvar, codir=codir, 
            			      nh = nh, hmax = hmax, asCov=asCov, flag.envelop=flag.envelop,
            			      env_color = env_color, env_linetype = env_linetype, 
            			      env_size=env_size, color = color, ...)
            }
        }
        
		# Adding some decoration
		g = plot.decoration(g, xlab = "Distance", ylab = "Variogram")
		
		# Add this plot to the list
		plot_lst[[index]] <- g
    } 

	if (multi_figure)
		p = ggarrange(plotlist=plot_lst, nrow=ivarN, ncol = jvarN)
  	else
  		p = plot_lst[[1]]
  		
    p
}

readPointCoor <- function(db)
{
  x = db$getCoordinates(0,TRUE)
  y = db$getCoordinates(1,TRUE)
  df = data.frame(x,y)
  df
}

readGridCoor <- function(dbgrid, name, usesel= FALSE)
{
  x = dbgrid$getColumnByLocator(ELoc_X(),0, FALSE, FALSE)
  y = dbgrid$getColumnByLocator(ELoc_X(),1, FALSE, FALSE)
  data = dbgrid$getColumn(name, usesel, FALSE)
  if (length(data) != length(x))
  {
    cat("Variable",name,"does not exist or does not have correction dimension\n")
    stop()
  }
  df = data.frame(x,y,data)
  df
}

# Function for plotting a point data base, with optional color and size variables
pointSymbol <- function(db, name_color=NULL, name_size=NULL,
                      	sizmin=10, sizmax=100, flagAbsSize = FALSE, flagCst=FALSE,
                      	...) 
{  
  # Creating the necessary data frame
  df = readPointCoor(db)
   
  # Color of symbol
  colval = NULL
  if (! is.null(name_color)) {
    colval  = db$getColumn(name_color, TRUE)
  }
  df["colval"] = colval 

  # Size of symbol
  sizval = NULL
  if (! is.null(name_size)) {
  	if (! flagCst)
   	{
	  	reduction = 100
  	  	sizval  = db$getColumn(name_size, TRUE)
  		if (flagAbsSize) sizval = abs(sizval)
    	m = min(sizval,na.rm=TRUE)
    	M = max(sizval,na.rm=TRUE)
    	sizval = (sizmax * (sizval - m) / (M-m) + sizmin) / reduction
    }
  }
  df["sizval"] = sizval

  layer <- geom_point(data = df, mapping = aes(x=x, y=y, color=colval, size=sizval), 
  		   na.rm=TRUE, ...)
  layer
}

# Function for plotting a point data base, with optional color and size variables
pointLabel <- function(db, name, digit=2, ...) 
{  
  # Creating the necessary data frame
  df = readPointCoor(db)
    
  # Label of symbols
  labval  = round(db$getColumn(name,TRUE),digit)
  df["labval"] = as.character(labval)
  
  layer <- geom_text(data = df, mapping=aes(x=x, y=y, label=labval), ...)
               
  layer
}

# Function for plotting a point data base, with optional color and size variables
plot.point <- function(db, name_color=NULL, name_size=NULL, name_label=NULL,
                       sizmin=10, sizmax=100, flagAbsSize = FALSE, flagCst=FALSE,
                       color_label="black", nudge_y=0.1, digit_label=2,
                       show.legend.symbol=FALSE, legend.name.color="P-Color",
                       legend.name.size="P-Size",
                       show.legend.label=FALSE, legend.name.label="P-Label",
                       show.title = TRUE, padd = NULL, ...) 
{ 
  p <- getNewFigure(padd, 2)
  title = ""
  
# If no variable is defined, use the default variable for Symbol(size) representation
# The default variable is the first Z-locator one, or the last variable in the file
	if (is.null(name_color) && is.null(name_size) && is.null(name_label))
	{
    	if (db$getVariableNumber() > 0)
            name_size = db$getNameByLocator(ELoc_Z(),0)
        else 
        {
        	# if no Z locator, choose the last field
            name_size = db$getLastName()
            flagCst = TRUE
		}
	}
  
  if (! is.null(name_color) || ! is.null(name_size))
  {
	  p <- p + pointSymbol(db, name_color=name_color, name_size=name_size,
                      	   sizmin=sizmin, sizmax=sizmax, flagAbsSize = flagAbsSize, flagCst=flagCst,
                      	   show.legend = show.legend.symbol,
                      	   ...)
 	  
 	  # Set the default title
	  if (! is.null(name_color))
	 	  title = paste(title, name_color)
	  if (! is.null(name_size))
	  	  title = paste(title, name_size)
	  
	  # Set the Legend
	  p <- p + labs(color = legend.name.color)
	  p <- p + labs(size = legend.name.size)
  }
  
  if (! is.null(name_label))
  {
  	p <- p + pointLabel(db, name=name_label, digit=digit_label, ...)
                    
    # Set the title  	   	 
   	title = paste(title, name_label)
    
    # Set the legend
   	p <- p + labs(label = legend.name.label)
  }
  
  # Decoration
  if (show.title)
	  p <- plot.decoration(p, title = title)

  p
}

gridRaster <- function(dbgrid, name, usesel = TRUE, ...)
{
  # Reading the Grid information
  df = readGridCoor(dbgrid, name, usesel)
  
  # Define the contents
  if (dbgrid$getAngles()[1] == 0)
  {
  	layer <- geom_tile(data = df, mapping=aes(x = x, y = y, fill = data), ...)
  }
  else
  {
    ids = seq(1, dbgrid$getNTotal())
    coords = dbgrid$getAllCellsEdges()
    positions = data.frame(id = rep(ids, each=4), x=coords[[1]], y=coords[[2]])
    values = data.frame(id = ids, value = df$data)
    df <- merge(values, positions, by = c("id"))
    layer <- geom_polygon(data = df, mapping=aes(x = x, y = y, fill = value, group = id), ...)
  }
  layer
}

gridContour <- function(dbgrid, name, usesel = TRUE, ...)
{
  # Reading the Grid information
  df = readGridCoor(dbgrid, name, usesel)
  
  layer <- geom_contour(data = df, mapping=aes(x = x, y = y, z = data), ...)

  layer
}

#
# Function for plotting a variable informed in a grid Db
#
plot.grid <- function(dbgrid, name_raster=NULL, name_contour=NULL,
					  usesel = TRUE, 
    				  option="B", na.color = "white", zlim = NULL, 
    				  bins = 10, line.color="black",
      				  show.legend.raster=FALSE, legend.name.raster="G-Raster", 
      				  show.legend.contour=FALSE, legend.name.contour="G-contour", 
      				  show.title = TRUE, 
      				  padd=NULL, ...)
{
  if (! dbgrid$isGrid())
  {
    cat("This function is restricted to Grid Db and cannot be used here\n")
    stop()
  }

  p <- getNewFigure(padd, 2)
  title = ""
  
  # If no variable is defined, use the default variable for Raster representation
  # The default variable is the first Z-locator one, or the last variable in the file
  if (is.null(name_raster) && is.null(name_contour))
  {
	if (dbgrid$getVariableNumber() > 0)
       name_raster = dbgrid$getNameByLocator(ELoc_Z(),0)
    else
    # if no Z locator, choose the last field
       name_raster = dbgrid$getLastName()
  }
  
  # Raster representation
   
  if (! is.null(name_raster))
  {
  	p = p + gridRaster(dbgrid, name=name_raster, usesel=usesel, ...)
  	p <- p + scale_color_viridis_c(option = option, na.value = na.color, limits=zlim)
  
  	# Set the title
  	title = paste(title,name_raster)
  	
	# Set the Legend
	if (show.legend.raster)
	  	p <- p + guides(fill = guide_colorbar(title=legend.name.raster, reverse=FALSE))
  }
  
  # Contour representation
  
  if (! is.null(name_contour))
  {
  	p = p + gridContour(dbgrid, name=name_contour, usesel=usesel, ...)
  	
  	# Set the title					
  	title = paste(title, name_contour)
  		
  	# Set the Legend
  	if (show.legend.contour)	 
	  	p <- p + labs(contour = legend.name.contour)
  }  
  
  # Decoration
  if (show.title)
	  p <- plot.decoration(p, title = title)
  
  p
}

# Function to display a polygon (not tested)

plot.polygon <- function(poly, padd = NULL, ...)
{
  npol = poly$getPolySetNumber()
  
  p <- getNewFigure(padd, 2)
  
  for (ipol in 1:npol)
  {
    df = data.frame(x = poly$getX(ipol-1), y = poly$getY(ipol-1))
    p <- p + geom_polygon(data = df, mapping=aes(x=x,y=y),  ...)
  }  
  p
}
        
# Function for plotting the histogram of a variable
plot.hist <- function(db, name, usesel=TRUE, nbins=30, col='grey', fill='yellow', padd = NULL)
{
  p <- getNewFigure(padd, 1)
  
  val  = db$getColumn(name, usesel)
  df = data.frame(val)
    
  p <- p + geom_histogram(data=df, 
  				mapping=aes(x=val), bins=nbins, color=col, fill=fill, na.rm=TRUE) 
  				
  p <- plot.decoration(p, title=name)
  
  p
}

# Function for plotting histogram for a table of values
plot.hist_tab <- function(val, nbins=30, padd=NULL)
{
  p <- getNewFigure(padd, 1)
  
  df = data.frame(val)
     
  p <- p + geom_histogram(data = df, 
  						  mapping=aes(x=val), bins=nbins, color='grey', fill='yellow') 

  p
}

# Function for plotting a curve of regularly sampled values
plot.curve <- function(data, color="black", padd=NULL)
{
  p <- getNewFigure(padd, 1)
  
  absc = seq(1,length(data))
  df = data.frame(absc,data)
  
  p <- p + geom_line(data = df, 
  					 mapping=aes(x=absc,y=data), color=color, na.rm=TRUE)
  
  p
}

# Function for representing a line between points provided as arguments
plot.XY <-function(x, y, join=TRUE,
                   color="black", linetype="solid", shape=20,
                   flagDiag = FALSE, 
                   diag_color = "red", diag_line = "solid", padd=NULL)
{
  if (length(y) != length(x))
  {
    cat("Arrays 'x' and 'y' should have same dimensions")
    stop()
  }
  
  p <- getNewFigure(padd, 1)
    
  df = data.frame(x, y)
     
  if (flagDiag)
  {
    u = min(x, y, na.rm=TRUE)
    v = max(x, y, na.rm=TRUE)
    p <- p + geom_segment(aes(x=u,y=u,xend=v,yend=v),
                          linetype = diag_line, color = diag_color, na.rm=TRUE)
  }
  
  if (join)
    p <- p + geom_line(data = df, mapping=aes(x=x,y=y), 
   		               linetype = linetype, color=color, na.rm=TRUE)
  else 
    p <- p + geom_point(data = df, mapping=aes(x=x,y=y), 
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
              
  p <- plot.geometry(p, xlim=res$getAylim(), ylim=res$getAzlim())
  
  p <- plot.decoration(p, xlab = "X", ylab = "Z")
  
  p
}

# Function for representing a scatter plot
plot.correlation <- function(db1, name1, name2, db2=NULL, usesel=FALSE,
							 flagDiag = FALSE,
                             color="black", linetype = "solid",
                             diag_color = "red", diag_line = "solid", padd=NULL)
{
  if (is.null(db2)) db2 = db1
  val1 = db1$getColumn(name1, usesel)
  val2 = db2$getColumn(name2, usesel)
  p = plot.XY(val1, val2, join=FALSE, flagDiag=flagDiag, 
              color = color, linetype = linetype, 
              diag_color = diag_color, diag_line = diag_line, padd=padd)
  
  p = plot.decoration(p, xlab=name1, ylab=name2)
  
  p 
}

# Representing a Lithotype rule
plot.rule <- function(rule, proportions=NULL, padd=NULL)
{
  p <- getNewFigure(padd, 1)
  
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
    rect = rule$getThresh(ifac-1)
    df$xmin[ifac] = rect[1]
    df$xmax[ifac] = rect[2]
    df$ymin[ifac] = rect[3]
    df$ymax[ifac] = rect[4]
  }
     
  p <- p + geom_rect(data = df, mapping=aes(xmin = xmin, xmax = xmax, 
                     ymin = ymin, ymax = ymax, fill = colors))
  
  p
}
 
 
# Function to display a polygon (not tested)
plot.mesh <- function(mesh, 
                      flagEdge=TRUE, flagFace=FALSE, flagApex=FALSE, 
                      facecolor="yellow", edgecolor="blue", linewidth=1,
                      show.legend = FALSE, padd = NULL)
{
  p <- getNewFigure(padd, 2)
  
  if (! flagFace) facecolor = "white"
  if (! flagEdge) edgecolor = facecolor
  
  nmesh = mesh$getNMeshes()
  for (imesh in 1:nmesh)
  {
    x = mesh$getCoordinatesPerMesh(imesh-1, 0, TRUE)
    y = mesh$getCoordinatesPerMesh(imesh-1, 1, TRUE)
    df = data.frame(x, y)
    p <- p + geom_polygon(data = df, mapping=aes(x=x,y=y), 
                          linewidth=linewidth, fill=facecolor, 
                          color=edgecolor, show.legend=show.legend)
  if (flagApex)
    p <- p + geom_point(data = df, mapping=aes(x=x, y=y))
  }  
  p
}

setMethod("plot", signature(x="_p_AMesh"), function(x,y=missing,...)   plot.mesh(x,...))
setMethod("plot", signature(x="_p_DbGrid"), function(x,y=missing,...)  plot.grid(x,...))

setMethod("plot", signature(x="_p_Db"), function(x,y=missing,...) plot.point(x,...))
setMethod("plot", signature(x="_p_Polygons"), function(x,y=missing,...) plot.polygon(x,...))

setMethod("plot", signature(x="_p_Vario"), function(x,y=missing,...) plot.vario(x,...))
setMethod("plot", signature(x="_p_Model"), function(x,y="missing",...) plot.model(x,...))
