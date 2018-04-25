%module smg2s
%{
#include <iostream>
#include <ctime>
#include <cstdlib>
#include <map>
#include "../../utils/utils.h"
%}

%include <stl.i>
%include "../../utils/utils.h"

%template(NilpotencyInt) Nilpotency<int>;



