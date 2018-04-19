# SMG2S
Sparse Matrix Generator with Given Spectrum


===============================================================================


It is underdevelopment by @[Xinzhe Wu](https://brunowu.github.io)...

The SpMV functionality implemented inside is a hypergraph based version which aims at avoiding the communicationi based on the paper [[Chen2014](https://link.springer.com/chapter/10.1007/978-3-319-17353-5_1)].

## Installation
### Source Codes
### Shared Library

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

