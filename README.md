# DUNEprotoTKI: tools for TKI analyses on protoDUNE data

## How-to's:

1. The symbolic links in 
   ```
   anaData/input/
   ```
   need to be unlinked (unlink) and then relinked (ln -s) to the actual paths of the data. (Do not push root files to repository!)

   - Printouts of the data structure can be found in doc/.

2. Setup ROOT and then do
   ```
   source setup.sh
   ```

3. Go to style/ and compile 
   ```
   ./do.sh
   ```
   - style/ is the lower level tool for all other programmes, needed for histogram manipulation.

4. Go to anaData/, and then drawTKI/, drawTracking/ if needed:
   - for truth signal distributions: anaData/ with kTruth=1 + drawTKI/
   - for efficiency, momentum resolution, and dEdx: anaData/ with kTruth=1 + drawTracking/
   - for event selection: anaData/ with kTruth=0 alone
   run
   ```
   ./do.sh
   ```
   and all output will be in the ./output/ subdirectories with all printout in see*.log.





