# SMG2S
Sparse Matrix Generator with Given Spectrum

-------------------------------------------------------------------------------

Author [Xinzhe Wu](https://brunowu.github.io) @ [Maison de la Simulation](http://www.maisondelasimulation.fr), France.

![Matrix Generation Pattern](figure/matgen.png)

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
mpirun -np ${PROCS} ./smg2s.exe -SIZE ${MAT_SIZE} -L ${LOW_BANDWIDTH} -C ${CONTINUOUS_ONES} -SPTR ${GIVEN_SPECTRUM_FILE}
```

If ${GIVEN_SPECTRUM_FILE} is not given, SMG2S will use the internal eigenvalue generation method to generate a default spectrum.

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
Mt = smg2s<std::complex<float>,int>(probSize, nilp, lbandwidth, spectrum);

```

### Given Spectra file in pseudo-Matrix Market Vector Format

#### Complex values
For the complex values, the given spectrum is stored in three columns, the first column is the coordinates, the second column is the real part of complex values, and the third column is the imaginary part of complex values.

    %%MatrixMarket matrix coordinate complex general
    3 3 3
    1 10 6.5154
    2 10.6288 3.4790
    3 10.7621 5.0540

#### Real Values
For the real values, the given spectrum is stored in two columns, the first column is the coordinates, the second column is related values.

    %%MatrixMarket matrix coordinate real general
    3 3
    1 10
    2 10.6288    
    3 10.7621

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


## Verification

![Comparison of generated spectrum with given spectrum](figure/vector.png)

