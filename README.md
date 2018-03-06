### Windows

First clone and build libOTe, libPoly. libOTe, libPoly and libPSU should share the same parent directory. Then clone this library and open the solution in Visaul Studio.

### Linux


libOTe, libPoly and libPSU should share the same parent directory.

```

git clone git@github.com:nitrieu/psu_impl.git
cd psu_impl
[libOTe clone build steps](https://github.com/nitrieu/libOTe)
[libPoly clone build steps] git@github.com:nitrieu/libPoly.git
cmake .
make