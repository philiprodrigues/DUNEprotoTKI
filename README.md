# DUNEprotoTKI: tools for TKI analyses on protoDUNE data

How-to's:

1) * The symbolic link in 
     anaData/input/
     needs to be unlinked and then relinked (ln -s) to the actual path of the data. (Do not push root files to repository!)
   * A print out of the data structure can be found in doc.
2) * Setup ROOT and then do
     source setup.sh
   * Go to style/ and compile by executing do.sh only if style is changed
3) * Go to anaData/, drawTKI/, drawTracking/ and executing do.sh
   - for truth signal distributions: anaData with kTruth=1 + drawTKI
   - for efficiency, momentum resolution, and dEdx: anaData with kTruth=0 -> drawTracking
   - for event selection: anaData with kTruth=0



