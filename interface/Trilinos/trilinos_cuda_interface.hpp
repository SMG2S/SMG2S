/*
   This file is part of SMG2S.
   Author(s): Xinzhe WU <xinzhe.wu@ed.univ-lille1.fr or xinzhe.wu1990@gmail.com>
        Date: 2018-04-20
   Copyright (C) 2018-     Xinzhe WU

   SMG2S is free software: you can redistribute it and/or modify
   it under the terms of the GNU Lesser General Public License as published
   by the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.
   SMG2S is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU Lesser General Public License for more details.
   You should have received a copy of the GNU Lesser General Public License
   along with SMG2S.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef __TRILINOS_INTERFACE_H__
#define __TRILINOS_INTERFACE_H__

#include "../../parMatrix/parMatrixSparse.h"
#include <math.h>
#include <complex.h>
#include "../../utils/MPI_DataType.h"
#ifdef __APPLE__
#include <sys/malloc.h>
#else
#include <malloc.h>
#endif

#include <Teuchos_CommandLineProcessor.hpp>

#include <Teuchos_GlobalMPISession.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_CrsMatrix.hpp>

using Tpetra::CrsMatrix;
using Tpetra::global_size_t;
using Teuchos::Array;
using Teuchos::ArrayView;
using Teuchos::ArrayRCP;
using Teuchos::arcp;
using Teuchos::tuple;
using Tpetra::Map;
using std::vector;

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::tuple;
using std::cout;
using std::endl;

typedef Kokkos::Compat::KokkosDeviceWrapperNode<Kokkos::Cuda> node_type;

#if defined(__USE_COMPLEX__) && defined(__USE_DOUBLE__) && defined (__USE_64BIT__)
typedef std::complex<double>               ST;
typedef Teuchos::ScalarTraits<ST>          SCT;
typedef int 							   local_ordinal_type;
typedef long 							   global_ordinal_type;

#elif defined(__USE_COMPLEX__) && defined(__USE_DOUBLE__)
typedef std::complex<double>               ST;
typedef Teuchos::ScalarTraits<ST>          SCT;
typedef int 							   local_ordinal_type;
typedef int 							   global_ordinal_type;

#elif defined (__USE_COMPLEX__) && defined(__USE_64BIT__)
typedef std::complex<float>               ST;
typedef Teuchos::ScalarTraits<ST>          SCT;
typedef int 							   local_ordinal_type;
typedef long 							   global_ordinal_type;

#elif defined (__USE_COMPLEX__)
typedef std::complex<float>               ST;
typedef Teuchos::ScalarTraits<ST>          SCT;
typedef int 							   local_ordinal_type;
typedef int 							   global_ordinal_type;

#elif defined (__USE_DOUBLE__) && defined(__USE_64BIT__)
typedef double               ST;
typedef Teuchos::ScalarTraits<ST>          SCT;
typedef int 							   local_ordinal_type;
typedef long 							   global_ordinal_type;

#elif defined (__USE_DOUBLE__)
typedef double               ST;
typedef Teuchos::ScalarTraits<ST>          SCT;
typedef int 							   local_ordinal_type;
typedef int 							   global_ordinal_type;

#elif defined (__USE_64BIT__)
typedef float               ST;
typedef Teuchos::ScalarTraits<ST>          SCT;
typedef int 							   local_ordinal_type;
typedef long 							   global_ordinal_type;

#else
typedef float               ST;
typedef Teuchos::ScalarTraits<ST>          SCT;
typedef int 							   local_ordinal_type;
typedef int 							   global_ordinal_type;

#endif

RCP<CrsMatrix<ST, local_ordinal_type, global_ordinal_type,node_type> > ConvertToTrilinosCudaMat(parMatrixSparse<ST, global_ordinal_type > *M){

	RCP<const Teuchos::Comm<int> > comm =	
    	Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();

	parVectorMap<global_ordinal_type> *pmap = M->GetYMap();

	global_ordinal_type gsize = pmap->GetGlobalSize();

	RCP<const Map<> > map = rcp (new Map<> (gsize, 0, comm));

  RCP<CrsMatrix<ST, local_ordinal_type, global_ordinal_type,node_type> > K = rcp (new CrsMatrix<ST, local_ordinal_type, global_ordinal_type,node_type> (map, 0,  Tpetra::DynamicProfile));

  const local_ordinal_type numMyElements = map->getNodeNumElements();
  
  ArrayView<const global_ordinal_type> myGlobalElements = map->getNodeElementList ();

  global_ordinal_type i;

  global_ordinal_type col;

 	ST val;

 	typename std::map<int ,ST>::iterator it;

 	std::map<global_ordinal_type,ST> *dynloc;
 	dynloc = M->GetDynMatLoc();
 	for(i = 0; i < numMyElements; i++){
 		for(it = dynloc[i].begin(); it != dynloc[i].end(); ++it){
 			col = it->first;
 			val = it->second;

	 		K->insertGlobalValues(myGlobalElements[i],tuple<global_ordinal_type>(col),tuple<ST>(val));
 		}
 	}

 	//K->fillComplete ();

 	return K;
}

#endif