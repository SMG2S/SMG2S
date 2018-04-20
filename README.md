# SMG2S
Sparse Matrix Generator with Given Spectrum

===============================================================================


Author [Xinzhe Wu](https://brunowu.github.io) @ [Maison de la Simulation](http://www.maisondelasimulation.fr), France.


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

Include header file

```cpp
#include <smg2s/smg2s.h>
```

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

## Interface
The cmake will check if PETSc is installed in the platfrom, if yes, header file to interface will also be copied to ${INSTALL_DIRECTORY}/include when installing SMG2S.

### Interface to PETSc

Include header file

```cpp
#include <interface/PETSc/petsc_interface.h>
```

Create parMatrixSparse type matrix

```cpp
parMatrixSparse<std::complex<double>,int> *Mt;
```

Restore this matrix into CSR format

```cpp
Mt->Loc_ConvertToCSR();
```

Create PETSc MAT type
```cpp
MatCreate(PETSC_COMM_WORLD,&A);
```

Convert to PETSc MAT format
```cpp
A = ConvertToPETSCMat(Mt);
```

More information: [PETSc GMRES example](https://github.com/brunowu/SMG2S/tree/master/example/gmres) and [PETSc Arnoldi example](https://github.com/brunowu/SMG2S/tree/master/example/arnoldi)

### Interface to Trilinos and other libraries

Coming soon.


![Matrix Generation Pattern](figure/matgen.png)

## Verification

![Comparison of generated spectrum with given spectrum](figure/vector.png)

