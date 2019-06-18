# LEMming2

This repo is mainly for archival purposes. It contains the code used and described by Ward (2019) Lithosphere. Given the inputs in that paper it should generate the same results.

However it is GPL software so use it as you wish. No warranty, expressed or implied. Validate against your use case before considering the results.

The code can't be considered to have a software architecture - it is more like the formless, unconsolidated nest of a software dove.
Mostly it is one long script, fed by global variables and broken up across files, but it generally follows DRY. 

Input is basically a namelist that sets all the parameters, and a couple examples are included as "LEMming_Input*.m". They're pretty well-commented but some comments may be cryptic.

This version is set (hacked, really) for detachment-limited erosion only.

You might need to re-MEX the following components: "get_neighbors.c" and "LocalQ.c". MEX-binaries are included for Mac 64-bit systems.

If you have a question feel free to email the developer: dylan.ward [at] uc.edu.