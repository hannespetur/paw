I[![Build Status](https://travis-ci.org/hannespetur/paw.svg?branch=master)](https://travis-ci.org/hannespetur/paw)

## paw
Paw libraries are header-only libraries made to be both simple and convenient to use.


### paw::Parser
paw::Parser is a library for parsing command-line arguments. Its most distinctive feature that it automatically parses argument values to their respective types and creates a neat help page based on the available option of the program.

**Dependencies**: C++11 support (GCC >= 4.8.1, Clang >= 3.3)

### paw::Station
paw::Station is a threadpool library built on top of the `std::thread` class. Each "station" manages the work to available threads.

**Dependencies**: C++11 support (GCC >= 4.8.1, Clang >= 3.3), multi-threading library (compile with the `-pthread` flag on Linux)

### paw::Align
paw::Align is a pairwise alignment library. The alignments are SIMD optimized and the library includes backtracing. The library is compiled with various CPU extension and the optimal one is selected at runtime.

**Dependencies**: C++11 support (GCC >= 4.8.1, Clang >= 3.3)

### More libraries to come...
Hopefully.

### Usage
Paw libraries have separated their implementations from their declarations to allow compiling implementations only once, regardless how often its header files are included in your program. Users of paw should put `#define IMPLEMENT_PAW` in **EXCACTLY** one of their `.cc/.cpp` file to achieve this. Without the definitions, all implementations will be ignored which greatly reduces compilation times. A good practice is to have one file with nothing but these two lines:
```cpp
#define IMPLEMENT_PAW
#include "paw.hpp" // Or only the header file(s) you need
```
Paw implementations will only need to be recompiled if this file is changed or the paw library is updated, which should likely happen very rarely.

Another possibility is build a library with all the paw implementations.
```sh
mkdir build
cd build
cmake ..
make shared # Builds lib/libpaw.so in the build directory
# Or "make static" to build a static library
```
Then, you no longer need a single file with `#define IMPLEMENT_PAW` but instead you can simply link your program to `libpaw.so` (or `libpaw.a`). Use whichever method that suits you and enjoy using paw!

### Examples
See the `examples` directory.

### Author
Hannes P Eggertsson

### License
GNU GPLv3
