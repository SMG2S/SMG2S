/*

MIT License

Copyright (c) 2019 Xinzhe WU

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
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

template<typename T, typename S>
RCP<CrsMatrix<T> > ConvertToTrilinosMat(parMatrixSparse<T, S > *M){


  typedef T               ST;
  typedef Teuchos::ScalarTraits<ST>          SCT;
  typedef int                  local_ordinal_type;
  typedef S                 global_ordinal_type;

	RCP<const Teuchos::Comm<int> > comm =	
  Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();

	parVectorMap<global_ordinal_type> *pmap = M->GetYMap();

	global_ordinal_type gsize = pmap->GetGlobalSize();

	RCP<const Map<> > map = rcp (new Map<> (gsize, 0, comm));

  RCP<CrsMatrix<ST> > K = rcp (new CrsMatrix<ST> (map, 0,  Tpetra::DynamicProfile));

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

// 	K->fillComplete ();

 	return K;
}

#endif
