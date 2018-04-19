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

```cpp
Nilpotency<int> nilp;

nilp.NilpType1(length,probSize);

parMatrixSparse<std::complex<float>,int> *Mt;

Mt = smg2s<std::complex<float>,int>(probSize, nilp, lbandwidth);

```

![Matrix Generation Pattern](figure/matgen.eps)

## Verification

