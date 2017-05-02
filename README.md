## Paw
The paw library is a header-only C++11 library.
Its main design goal is to be very simple and convenient to use.

Currently, the library only contains an argument parser but in the future more functionality will be added.

### Parser
`parser.hpp` is a small header-only file for parsing command-line arguments.
The most distinctive feature of the parser is that it automatically parses argument values to their respective types using the `std::istream <<` operator.

Note that paw parser is still under development, so its API is subject to change.

#### TODO
 - Generate the "help" page
 - Write more unit tests
 - Write an example

### Author
Hannes P Eggertsson

### License
GNU GPLv3
