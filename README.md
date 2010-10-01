sot-dynamic
===========

This software provides robot dynamic computation for dynamic-graph
by using jrl-dynamics.

Setup
-----

To compile this package, it is recommended to create a separate build
directory:

    mkdir _build
    cd _build
    cmake [OPTIONS] ..
    make install

Please note that CMake produces a `CMakeCache.txt` file which should
be deleted to reconfigure a package from scratch.


### Dependencies

The matrix abstract layer depends on several packages which
have to be available on your machine.

 - Libraries:
   - [jrl-dynamics][jrl-dynamics] (>= 1.16.1)
   - [dynamic-graph][dynamic-graph] (>= 1.0.0)
   - sot-core (>= 1.0.0)
   - [jrl-mal][jrl-mal] (>= 1.8.0)
 - Closed source libraries:
   - hrp2Dynamics (>= 1.3.0)
   - hrp2-10-optimized (>= 1.3.0) [optional]
 - System tools:
   - CMake (>=2.6)
   - pkg-config
   - usual compilation tools (GCC/G++, make, etc.)


[jrl-dynamics]: http://github.com/jrl-umi3128/jrl-dynamics
[dynamic-graph]: http://github.com/jrl-umi3128/dynamic-graph
[jrl-mal]: http://github.com/jrl-umi3128/jrl-mal
