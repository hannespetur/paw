[![Build Status](https://travis-ci.org/hannespetur/paw.svg?branch=master)](https://travis-ci.org/hannespetur/paw)

## Paw
The paw library is a header-only C++11 library. Its main design goal is to be both simple and convenient to use.

Paw separates its implementations from its declarations to allow compiling implementations only once, regardless how often the header file is including in your program. Users of paw should put `#define IMPLEMENT_PAW` in **EXCACTLY** one `.cc/.cpp` file that includes paw to achieve this. Without the definition all implementations will be ignored, greatly reducing compilation times. A good practice is to have one file with nothing but these two lines:

```cpp
#define IMPLEMENT_PAW
#include "paw.hpp" // Or only the header file(s) you need
```

Paw implementations will only need to be recompiled if this file is changed or the library is update, which should likely happen very rarely.

Another possibility is build a shared (or a static) library with all the paw implementations.

```sh
mkdir build
cd build
cmake ..
make shared # Builds lib/libpaw.so in the build directory
#"make static" to build a static library
```

Then, you no longer need a single file with `#define IMPLEMENT_PAW` but instead you can simply link your program to `libpaw.so`. Use whichever method that suits you and enjoy using paw!

### paw::parser
paw::parser is an utility for parsing command-line arguments. Its most distinctive feature that it automatically parses argument values to their respective types (see example below) and can automatically creates a neat help page. Examples on how to use paw::parser are in the `examples` directory.

### More libraries to come...
Hopefully.

### Author
Hannes P Eggertsson

### License
GNU GPLv3
