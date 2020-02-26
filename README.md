![cnem](https://github.com/isolat/cnem/blob/master/logo.png)
# Description
This package contains an implementation of natural neighbour interpolation (exact Sibson and Laplace) in 3d (2d coming soon), liceneced under GNU Lesser General Public License.

This implementation is platform independent, and for python and Matlab users, the main lib is written in c++ and parallelized (shared memory).

Authors : amran.illoul(at)ensam.eu , philippe.lorong(at)ensam.eu

# Dependency
* **TETGEN** : To build the initial 3d constrained delaunay tetradedrisation we use TETGEN devellopped by Hang SI (open source library, AGPL v.3, http://www.tetgen.org/ ). TETGEN 1.5.1 is supplied in the package.

* **TBB** : for parallelization, open source library, gpl v2 licence, http://www.threadingbuildingblocks.org/ 

[![View cnem on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://fr.mathworks.com/matlabcentral/fileexchange/74351-cnem)
