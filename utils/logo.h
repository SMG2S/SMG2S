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
	printf(
	"	 ________   ____  ____    _________   ________   ________                     \n"
	"	||         ||   ||   ||  ||       ||         || ||                 \n"
	"	||______   ||   ||   ||  ||            ______|| ||______     \n"    
	"	|--------| ||   ||   ||  ||      ___ |------- | |--------|\n"
	"	        || ||   --   ||  ||       || ||                 ||     				\n"
	"	||______|| ||        ||  ||_______|| ||______|| ||______||\n"
	"			\n"
	);

    border_print();
    center_print("Developed by Xinzhe WU at Maison de la Simulation, France", 79);
    char v[100];
    sprintf(v, "Version: %.1f", version);
	center_print(v, 79);
	border_print();

	printf("\n\n\n");

}

#endif


