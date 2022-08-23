# RealFOAM

Extensions for newer OpenFOAM versions to include the common real gas equations of state. Implementation of the common cubic equations of state Soave-Redlich-Kwong, Peng-Robinson and Redlich-Kwong-Peng-Robinson using a generalised formulation of cubic equations of state. On branch `esi` for OpenFOAM-v2112 and OpenFOAM-v2206 (see tags) and on branch `org` for OpenFOAM 10.  

**Check out the desired branch `org` or `esi`.**

## Requirements
- Download OpenFOAM version - Important: Please download the source code in order to be able to later modify it (e.g. OpenFOAM v2112 can be downloaded under https://dl.openfoam.com/source/v2112/)
- Compile it on your laptop or workstation (**Attention: This takes about 4 hours**)

After successful compilation:
- Add the files from the corresponding repository, link them and recompile. 
- Run the test case tutorials.

## References
For further details on the thermodynamic model and a validation of the implementations, we refer to the corresponding paper _Trummler et al., Thermodynamic modeling for numerical simulations based on the generalized cubic equation of state._

## Citation & License information
If you use the OpenFOAM extensions provided here for your work, cite the corresponding paper: _Trummler, T., Glatzle, M., Doehring, A., Urban, N., Klein, M.,Thermodynamic modeling for numerical simulations based on the generalized cubic equation of state._

See LICENSE for information on the terms & conditions for usage of this software, and a DISCLAIMER OF ALL WARRANTIES.