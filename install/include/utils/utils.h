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

#ifndef __UTILS_H__
#define __UTILS_H__

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <map>

template<class T>
T random_unint(T min, T max)
{
/*
	std::default_random_engine e(time(NULL));
	std::uniform_real_distribution<T> u(min, max);
*/	
	return (min + static_cast<T>(max * rand() / static_cast<T>(RAND_MAX + 1)));

}

template<class T, class S>
T random(S min, S max)
{
	return (min + (rand()%(max - min + 1)));

}

template<class S>
S factorial(S start, S end)
{
	S i;
	S value;

	if(start > end){
		value = 1;
	}
	else{
		value = start;
		for(i = start + 1; i <= end; i++){
			value *= i;
		}
	}

	return value;
}

template<typename S>
struct Nilpotency
{
	S	diagPosition; //off-diagonal offset
	S	nbOne; //continuous 1 number of nilpotent matrix
	S   matrix_size; //matrix size
	S   nilpotency;
	bool setup;

	Nilpotency()
	{
		diagPosition = 0;
		nbOne = 0;
		matrix_size = 0;
		nilpotency = 0;
		setup = false;
	};

	Nilpotency(S offset, S num, S size, S nil)
	{
		diagPosition = offset;
		nbOne = num;
		matrix_size = size;
		nilpotency = nil;
		setup = true;
	};

	void NilpType1(S num, S size)
	{
		diagPosition = 2;
		nbOne = num;
		matrix_size = size;
		nilpotency = num+1;
		setup = true;
	}

	void NilpType2(S num, S size)
	{
		if(num%2 == 0){
			diagPosition = 3;
			nbOne = num;
			matrix_size = size;
			nilpotency = num+1;
			setup = true;
		}
		else{
			setup = false;
			printf("Please choose the right nb of continuous 1, for NilpType2, it should be divisible by 2 \n");	
		} 
	}

	~Nilpotency()
	{

	};

};

#endif
