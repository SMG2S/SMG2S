/*
MIT License
Copyright (c) 2019 Xinzhe WU @ Maison de la Simulation, France
Copyright (c) 2019-2022, Xinzhe Wu @ Simulation and Data Laboratory Quantum 
									 Materials,  Forschungszentrum Juelich GmbH.
									 
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

#ifndef __LOGO_H__
#define __LOGO_H__

#include <iostream>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

void border_print(void)
{
	printf(
	"==================================================================="
	"=======================================\n");
}

void border_print2(void)
{
	printf(
	"\n-------------------------------------------------------------------"
	"---------------------------------------\n\n");
}

void center_print(const char *s, int width)
{
	int length = strlen(s);
	int i;
	for (i=0; i<=(width-length)/2; i++) {
		fputs(" ", stdout);
	}
	fputs(s, stdout);
	fputs("\n", stdout);
}


void logo(float version)
{
	border_print();

	printf("\n");
	printf("\n");
	printf("\n");
	center_print("SMG2S: Scalable Matrix Generator with Given Spetrum", 100);
	char v[100];
    sprintf(v, "Version: %.1f", version);
	center_print(v, 100);
	printf("\n");
	printf("\n");
	printf("\n");
	border_print();
    center_print("Developed by Xinzhe WU at Maison de la Simulation, France", 100);
	border_print();

	printf("\n\n\n");

}

static void show_usage(std::string name)
{
    std::cerr << "Usage: mpirun -np ${PROCS} " << name << " -SIZE ${MATRIX SIZE} -L ${LOWER BANDWIDTH} -C ${NB ONES FOR NILPOTENT MATRIX}\n"
              << "Options:\n"
              << "\t-h,--help\t\tShow this HELP message\n\n"
              << std::endl;
}


#endif


