################################################
#              THE CRAZY INDEXER
################################################
Rebuilding indexer is often necessary
But if you still have errors in the Problems view:

1- Close eclipse

2- Remove markers from the .metadata directory :
  
  .plugins/org.eclipse.core.resources/.projects/gstlearn/.markers

3- Open eclipse



################################################
#                 DEBUG TIPS
################################################

# Attach to process
===================

How to debug gstlearn whatever the target language?

1. Build debug version of gstlearn and make python_install DEBUG=1 or make r_install DEBUG=1
2. Launch your favorite software and import or load gstlearn package
3. Put a breakpoint in the C++ code where you want to stop
3. Run an "attach to process" debug session under eclipse
4. Connect to the following processes (check that the prompt is frozen at debug startup):

Python:
- Console : python3
- Spyder : python3 -m spyder_kernels.console -f XXX.json
- Jupyter : python3 -m ipykernel_launcher -f XXX.json

R:
- Console : exec/R
- rstudio : rsession

5. Continue the debug session (F5), prompt should become available
6. Execute the problematic command in the target language 
 