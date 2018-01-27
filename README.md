[![Build Status](https://travis-ci.org/hannespetur/paw.svg?branch=master)](https://travis-ci.org/hannespetur/paw)

## Paw
The paw library is a header-only C++11 library. Its main design goal is to be both simple and convenient to use.

Paw separates its implementations from its declarations to allow compiling implementations only once, regardless how often its header files are included in your program. Users of paw should put `#define IMPLEMENT_PAW` in **EXCACTLY** one of their `.cc/.cpp` file to achieve this. Without the definitions, all implementations will be ignored and thus greatly reduce compilation times. A good practice is to have one file with nothing but these two lines:

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

### paw::parser
paw::parser is an utility for parsing command-line arguments. Its most distinctive feature that it automatically parses argument values to their respective types (see example below) and can automatically creates a neat help page. Examples on how to use paw::parser are in the `examples` directory.

### More libraries to come...
Hopefully.

### Author
Hannes P Eggertsson

### License
GNU GPLv3
