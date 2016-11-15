#TinyLEGO Implementation

This is the C++14 implementation of the modified TinyLEGO [1] protocol appearing in the paper "Constant Round Maliciously Secure 2PC with Function-independent Preprocessing using LEGO" [2]. The code has been authored by Roberto Trifiletti and builds on the [OTExtension](https://github.com/encryptogroup/OTExtension) implementation of the ALSZ15 protocol [3] and uses the PRG code of [APRICOT](https://github.com/bristolcrypto/apricot). Also, the code of the OTExtension library in this project has been modified slightly to support extension of globally correlated random OTs (Delta-ROTs).

The codebase also contains a (more or less) stand-alone implementation of the recent additively homomorphic commitment scheme of [4].


##Description
The implementation was written with the purpose of exploring the practical efficiency of the TinyLEGO protocol for general secure two-party computation (2PC). It should therefore be treated as a prototype and we make no claim about actual security guarantees for any real world use cases. The code includes main and test functions that exhibit how to use the underlying commitment and 2PC code can be used. Included is also circuit representations of common cryptographic functions (available from https://www.cs.bris.ac.uk/Research/CryptographySecurity/MPC/) that are typically used for benchmarking MPC protocols.
At this time the implementation has a few limitations compared to the full potential of the general protocol:
* Due to simplicity the preprocessing, circuit building and evaluation phases can only be invoked once for the lifetime of the program, but we stress this is not due to a restriction from the TinyLEGO protocol. This has the effect that the caller needs to know a priori how many AND gates the final functionality requires as he cannot later produce more with the current codebase.
* The implementation uses no disk I/O whatsoever and the complexity of the desired secure function (# AND gates) is therefore bounded by the amount of RAM on the current machine.
* In light of the above restrictions, pipelining the evaluation of the final garbled circuit is not implemented, but if anyone wants to extend the code to handle this (or in any other way) you are very welcome to.

##Installation
The code has been tested to work on MAC OSX 10.11 (El Capitan), Amazon Linux 2016.03, Ubuntu 14.04, and Ubuntu 16.04.

###Requirements
* C/C++ compiler with C++14 support. The code has been successfully built with GCC 5.3.1, GCC 6.1 and CLANG 3.8.
* Required packages:
    * make
    * cmake
    * libssl-dev ([OpenSSL](https://www.openssl.org) for OTX code)
    * libgmp-dev ([GMP](https://gmplib.org) for arbitrary-precision arithmetic)
    * libzmq3-dev ([ZeroMQ](http://zeromq.org) for networking)

The above can be installed via a package manager, i.e. apt on Ubuntu, yum on Amazon Linux and Homebrew on OSX are typical.

After installing the above prerequisites the code can be built with the provided script _cmake-release.sh_. To verify everything is working the script _test/test-all.sh_ can afterwards be run. If all tests succeed you are good to go.

##Running the main files
Two main files are produced during compilation, build/release/Tinyconst and build/release/Tinyeval. An example run of the two clients on different machines could be
* [Machine A] ./build/release/Tinyconst -n 100 -c aes -e 8,4,2 -ip [A's IP] -p [port_num]
* [Machine B] ./build/release/Tinyeval -n 100 -c aes -e 8,4,2 -ip [A's IP] -p [port_num]

The above code precomputes enough AND gates for 100 secure computations of AES-128 (including key-expansion). The -e parameters specifies how many parallel executions this should be split into for the independent preprocessing, dependent preprocessing and online phase, respectively. In this example 8 parallel threads will produce the preprocessing required for 100 AES computations, 4 parallel threads will then build 25 AES circuits each using soldering and finally 2 parallel threads will evaluate 50 AES each. In order to measure latency for sequential evaluations the last argument of e should be 1. For information about which execution arguments were used for the results provided in [2] we refer to the timing_scripts folder as the performance benefit of the chosen execution parameters is highly platform and network dependent.

##References
* [1] T. K. Frederiksen, T. P. Jakobsen, J. B. Nielsen, R. Trifiletti, “TinyLEGO: An Interactive Garbling Scheme for Maliciously Secure Two-Party Computation,” IACR Cryptology ePrint Archive, vol. 2015, p. 309, 2015. [Online]. Available: http://eprint.iacr.org/2015/309.

* [2] J. B. Nielsen, T. Schneider, R. Trifiletti, "Constant Round Maliciously Secure 2PC with Function-independent Preprocessing using LEGO", to appear at The Network and Distributed System Security Symposium (NDSS) 2017.

* [3] G. Asharov, Y. Lindell, T. Schneider, M. Zohner, "More Efficient Oblivious Transfer Extensions with Security for Malicious Adversaries", in EUROCRYPT 2015, Part I, ser. LNCS, E. Oswald and M. Fischlin, Eds., vol. 9056. Springer, Apr. 2015, pp. 673–701.

* [4] T. K. Frederiksen, T. P. Jakobsen, J. B. Nielsen, R. Trifiletti, “On the Complexity of Additively Homomorphic UC commitments,” in TCC 2016-A, Part I, ser. LNCS, E. Kushilevitz and T. Malkin, Eds., vol. 9562. Springer, Jan. 2016, pp. 542–565.