#include "c_wrapper.h"
#include <stdio.h>

int main(int argc, char* argv[]) {

	struct NilpotencyInt *n;
	n = newNilpotencyInt();
	NilpType1(n, 2, 10);
	showNilpotencyInt(n);
	ReleaseNilpotencyInt(&n);

	return 0;
}
