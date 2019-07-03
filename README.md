# Private Set Union
This is the implementation of our paper: **Scalable Private Set Union from Symmetric-Key Techniques **[[ePrint](https://eprint.iacr.org/2019/xxx.pdf)]. 

Evaluating on a single server (`2 36-cores Intel Xeon CPU E5-2699 v3 @ 2.30GHz and 256GB of RAM`) with a single thread per party, each party has `2^20` items, our protocol requires  `239` seconds and `144` MB.

## Installations

### Required libraries
 C++ compiler with C++14 support. There are several library dependencies including [`Boost`](https://sourceforge.net/projects/boost/), [`Miracl`](https://github.com/miracl/MIRACL), [`NTL`](http://www.shoup.net/ntl/) with GMP, and [`libOTe`](https://github.com/osu-crypto/libOTe). For `libOTe`, it requires CPU supporting `PCLMUL`, `AES-NI`, and `SSE4.1`. Optional: `nasm` for improved SHA1 performance.   Our code has been tested on both Windows (Microsoft Visual Studio) and Linux. To install the required libraries: 
  * For building boost, miracl and libOTe, please follow the more instructions at [`libOTe`](https://github.com/osu-crypto/libOTe)
  * For NTL with GMP, `cd ./thirdparty`, and run `gmp.get` and `ntl.get`.   


### Building the Project
After cloning project from git,
##### Windows:
1. build cryptoTools,libOTe, and libPSU projects in order.
3. run frontend project
 
##### Linux:
1. make (requirements: `CMake`, `Make`, `g++` or similar)
2. for test:
	./bin/frontend.exe -t


## Running the code

	./bin/frontend.exe -r 1 
	& ./bin/frontend.exe -r 0 

		
## Help
For any questions on building or running the library, please contact [`Ni Trieu`](http://people.oregonstate.edu/~trieun/) at trieun at oregonstate dot edu

