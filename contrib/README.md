OpenMS contributing libraries
=============

This directory contains the OpenMS contrib package.

If you downloaded a stable release of OpenMS,
detailed installation instructions and other information
can be found in the documentation (`OpenMS/doc/index.html`).

If this is a development version obtained via Git or as
nightly snapshot, the documentation is not contained
in this package. Please refer to the OpenMS website
https://www.openms.de for installation instructions.

Please note that an internet connection is required to
build the contrib since the actual source packages are
only downloaded when built.

You can find prebuilt versions and the used sources
(e.g. for manual download) on our [archive](https://abibuilder.informatik.uni-tuebingen.de/archive/openms/contrib).

Build
============

To see all available build types, execute 

`$ cmake -DBUILD_TYPE=LIST .`

Using the `-DBUILD_TYPE` switch on the command line, you can customize which
external libraries should be built (depending on your system, you might already
have some of them installed and may not want to re-install them). 
CMake will then download the requested libraries and compile them for you.

Notes
============
Note that in case you do a partial build, it is wise to install BZIP2 and ZLIB first, since other libraries (e.g. BOOST) depend on it.
