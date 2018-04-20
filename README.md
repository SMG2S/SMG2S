# SMG2S
Sparse Matrix Generator with Given Spectrum


===============================================================================


It is underdevelopment by @[Xinzhe Wu](https://brunowu.github.io)...

The SpMV functionality implemented inside is a hypergraph based version which aims at avoiding the communicationi based on the paper [[Chen2014](https://link.springer.com/chapter/10.1007/978-3-319-17353-5_1)].

## Installation

### Binary Executable file
In the main directory:

```bash
cmake .  -DCMAKE_INSTALL_PREFIX=${INSTALL_DIRECTORY}
```

Install

```bash
make
```

Execution

```bash
mpirun -np ${PROCS} ./smg2s.exe -SIZE ${MAT_SIZE} -L ${LOW_BANDWIDTH} -C ${CONTINUOUS_ONES}
```

### Include files

Install the binary file and includes files into ${INSTALL_DIRECTORY}
```bash
make install
```

Include the SMG2S header file
```cpp
#include <smg2s/smg2s.h>
```

Include and Compile
```bash
mpicxx example.cpp -I${INSTALL_DIRECTORY}/include
```

## Example
### Creation
Generate the Nilpotent Matrix for operation:
```cpp
Nilpotency<int> nilp;
nilp.NilpType1(length,probSize);
```
Create the parallel Sparse Matrix Object Mt:
```cpp
parMatrixSparse<std::complex<float>,int> *Mt;
```
Generate a new matrix:
```cpp
Mt = smg2s<std::complex<float>,int>(probSize, nilp, lbandwidth);

```

![Matrix Generation Pattern](figure/matgen.png)

## Verification

![Comparison of generated spectrum with given spectrum](figure/vector.png)

