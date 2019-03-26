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

#ifndef __UTILS_H__
#define __UTILS_H__

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <map>
#include "../config/config.h"

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

	return static_cast<T>(min + (rand()%(max - min + 1)));

}


__int64_t factorial(__int64_t start, __int64_t end)
{
        __int64_t i;
        __int64_t value;

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

	void NilpType3(S diagP, S num, S size)
	{

		if (2*(diagP - 1)*num - 1 > size){
			setup = false;
			printf("Please choose the right parameter diaP, num and Size to satisfy the relation: 2* 2*(diagP - 1)*num - 1 <= size\n");
		}
		else if(num % (diagP - 1) != 0){
			setup = false;
			printf("Please choose the right parameter num = %d, for NilpType3 with diagP = %d, it should be divisible by diagP - 1 = %d \n", num, diagP, diagP - 1);
		}
		else{
			diagPosition = diagP;
			nbOne = num;
			matrix_size = size;
			nilpotency = num+1;
			setup = true;
		}
	}

	~Nilpotency()
	{

	};

};

#endif
