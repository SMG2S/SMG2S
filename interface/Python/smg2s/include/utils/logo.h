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
	"=============\n");
}

void border_print2(void)
{
	printf(
	"-------------------------------------------------------------------"
	"-------------\n");
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
/*
	printf(
	"	 ________   ____  ____    _________   ________   ________                     \n"
	"	||         ||   ||   ||  ||       ||         || ||                 \n"
	"	||______   ||   ||   ||  ||            ______|| ||______     \n"    
	"	|--------| ||   ||   ||  ||      ___ |------- | |--------|\n"
	"	        || ||   --   ||  ||       || ||                 ||     				\n"
	"	||______|| ||        ||  ||_______|| ||______|| ||______||\n"
	"			\n"
	);
*/
	printf("\n");
	printf("\n");
	printf("\n");
	center_print("SMG2S: Scalable Matrix Generator with Given Spetrum", 79);
	char v[100];
    sprintf(v, "Version: %.1f", version);
	center_print(v, 79);
	printf("\n");
	printf("\n");
	printf("\n");
	border_print();
    center_print("Developed by Xinzhe WU at Maison de la Simulation, France", 79);
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


