#include <C/c-smg2s.h>
#include <stdio.h>
#include <stdlib.h>
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
    double sparsity = 0.5;

    nilp_t *nilp = newNilp_1(2, 7);
    nilp_show(nilp);

    span = (int)ceil((double)probSize/(double)world_size);


    if(world_rank == world_size - 1){
        lower_b = world_rank * span;
        upper_b = probSize - 1 + 1;
    }else{
        lower_b = world_rank * span;
        upper_b = (world_rank + 1) * span - 1 + 1;
    }

    initMatrix_t *initMat = newInitMatrix_3(diagl, diagu, sparsity);  

    parVecMap_t *p = newParVecMap(MPI_COMM_WORLD, lower_b, upper_b);

    //non-symm case 1
	ds_parVec_t *spec = new_ds_ParVec_2(p);

    for(int i = lower_b; i < upper_b; i++){
        ds_parVecSetVal(spec, i, i + 1);
    }

    ds_parVecView(spec);

    //initMatrix_show (initMat);
    ds_parMatSparse_t *mat = ds_nonsymm_2(probSize, nilp, initMat, spec);    

    ds_parMatSparse_View(mat);  
    ds_parMatSparse_destory(mat);

    //non-symm case 2 with conjugate eigenvalues
    zs_parVec_t *spec2 = new_zs_ParVec_2(p);
    for(int i = lower_b; i < upper_b; i++){
        if(i % 2 == 0){
            dcomplex_t v = {i/2 + 1, i/2 + 2};
            zs_parVecSetVal(spec2, i, v);
        }else{
            dcomplex_t v = {i/2 + 1, -i/2 - 2};
            zs_parVecSetVal(spec2, i, v);
        }
        if(i == probSize - 1){
            dcomplex_t v = {i + 1, 0};
            zs_parVecSetVal(spec2, i, v);
        }
    }

    zs_parVecView(spec2);
    ds_parMatSparse_t *mat2 = ds_nonsymmconj_2(probSize, nilp, initMat, spec2);    
    ds_parMatSparse_View(mat2);  
    ds_parMatSparse_destory(mat2);

    //non-herm case 
    zs_parVec_t *spec3 = new_zs_ParVec_2(p);
    for(int i = lower_b; i < upper_b; i++){
        dcomplex_t v = {i+1, i+2};
        zs_parVecSetVal(spec3, i, v);
    }

    zs_parVecView(spec3);

    zs_parMatSparse_t *mat3 = zs_nonherm_2(probSize, nilp, initMat, spec3);    
    zs_parMatSparse_View(mat3);  
    zs_parMatSparse_destory(mat3);
    initMatrix_destory(initMat);

    ds_parVec_destory(spec);
    zs_parVec_destory(spec2);
    zs_parVec_destory(spec3);

    nilp_destory (nilp);

	MPI_Finalize();
}