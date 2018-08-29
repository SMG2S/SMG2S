#!/usr/bin/env python
"""
SMG2S Hello World using Python
"""

from mpi4py import MPI
import smg2s
import sys

size = MPI.COMM_WORLD.Get_size()
rank = MPI.COMM_WORLD.Get_rank()
name = MPI.Get_processor_name()

sys.stdout.write(
    "Hello, World! I am process %d of %d on %s.\n"
    % (rank, size, name))

if rank == 0:
	print ('INFO ]> Starting ...')
	print("INFO ]> The MPI Comm World Size is %d" %size)

lbandwidth = 3

nilp = smg2s.NilpotencyInt()

nilp.NilpType1(2,10)

if rank == 0:
	print("Nilptency matrix continuous one nb = %d" %nilp.nbOne)

Mt = smg2s.parMatrixSparseRealDoubleInt()

Mt=smg2s.smg2sRealDoubleInt(10,nilp,lbandwidth," ", MPI.COMM_WORLD)

Mt.LOC_MatView()


