args = commandArgs(trailingOnly=TRUE)
libdir = args[1]
dyn.load(paste0(libdir,"rgstlearn.so"))
cacheMetaData(1)
rgeostats <- new.env()
eval(parse(paste0(libdir,"rgstlearn.R")), envir=rgeostats)
attach(rgeostats)

# Add default show S4 methods to all objects (assuming that they inherits from AStringable)
cl = showMethods('$', printTo=FALSE)
cl = cl[grep("^x=\"_p_", cl)]
cl = gsub("^x=\"","",cl)
cl = gsub("\"","",cl)
invisible(lapply(cl,function(x) { 
  cmd = paste0("setMethod(f = \"show\", signature = \"",x,"\"",
               ", definition = function(object){ ",
               "AStringable_display(object) })")
  eval(parse(text=cmd))
}))
