#include <C/c-smg2s.h>
#include <mpi.h>
#include <math.h>

int main(int argc, char* argv[]) {

	MPI_Init(&argc, &argv);

    int world_size;
    int world_rank;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int probSize = 7;
    int span, lower_b, upper_b;

    int diagl = -5;
    int diagu = -3;
    int offset = 1;
    int nbOne = 2;
    double sparsity = 0.5;
    /* construct a nilpotent object for generation */
    nilp_t *nilp = newNilp_2(nbOne, offset, probSize);

    span = (int)ceil((double)probSize/(double)world_size);


    if(world_rank == world_size - 1){
        lower_b = world_rank * span;
        upper_b = probSize - 1 + 1;
    }else{
        lower_b = world_rank * span;
        upper_b = (world_rank + 1) * span - 1 + 1;
    }

    /* construct a initMat object for SMG2S*/
    initMatrix_t *initMat = newInitMatrix_3(diagl, diagu, sparsity);  
    /* construct a parVecMap object which determines the distribution scheme of vectors and matrices*/
    parVecMap_t *p = newParVecMap(MPI_COMM_WORLD, lower_b, upper_b);

    /* example 1, generation of a non-Hermtian matrix */
    // 1. generate the spectrum on the fly
    zs_parVec_t *spec1 = new_zs_ParVec_2(p);
    for(int i = lower_b; i < upper_b; i++){
        dcomplex_t v = {i+1, i+2};
        zs_parVecSetVal(spec1, i, v);
    }
    // 2. generation 
    zs_parMatSparse_t *mat1 = zs_nonherm_2(probSize, nilp, initMat, spec1);    
    zs_parMatSparse_destory(mat1);

    /* example 2, generation of a non-Symmetric matrix with real eigenvalues */
    // 1. generate the spectrum on the fly
	ds_parVec_t *spec2 = new_ds_ParVec_2(p);
    for(int i = lower_b; i < upper_b; i++){
        ds_parVecSetVal(spec2, i, i + 1);
    }
    // 2. generation 
    ds_parMatSparse_t *mat2 = ds_nonsymm_2(probSize, nilp, initMat, spec2);    
    ds_parMatSparse_destory(mat2);

    /* example 3, generation of a non-Symmetric matrix with conjugated eigenvalues */
    // 1. generate the spectrum on the fly
    zs_parVec_t *spec3 = new_zs_ParVec_2(p);
    for(int i = lower_b; i < upper_b; i++){
        if(i % 2 == 0){
            dcomplex_t v = {i/2 + 1, i/2 + 2};
            zs_parVecSetVal(spec3, i, v);
        }else{
            dcomplex_t v = {i/2 + 1, -i/2 - 2};
            zs_parVecSetVal(spec3, i, v);
        }
        if(i == probSize - 1){
            dcomplex_t v = {i + 1, 0};
            zs_parVecSetVal(spec3, i, v);
        }
    }
    // 2. generation 
    ds_parMatSparse_t *mat3 = ds_nonsymmconj_2(probSize, nilp, initMat, spec3);    
    ds_parMatSparse_destory(mat3);

    initMatrix_destory(initMat);
    zs_parVec_destory(spec1);
    ds_parVec_destory(spec2);
    zs_parVec_destory(spec3);
    nilp_destory (nilp);

	MPI_Finalize();
}