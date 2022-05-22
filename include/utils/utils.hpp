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

#ifndef __UTILS_H__
#define __UTILS_H__

#include <iostream>
#include <ctime>
#include <cstdlib>
#include <map>
#include<set>

template <class Q>
struct Base_Class {
  typedef Q type;
};

template <class Q>
struct Base_Class<std::complex<Q> > {
  typedef Q type;
};

template <typename Q>
using Base = typename Base_Class<Q>::type;


template<class T>
T random_unint(T min, T max)
{
/*
	std::default_random_engine e(time(NULL));
	std::uniform_real_distribution<T> u(min, max);
*/	
	return (min + static_cast<T>(max * rand() / static_cast<T>(RAND_MAX + 1)));

}

template<class T>
int pw(T x)
{
	int stop = 0;
	int count = 0;

	while (!stop){
		x = x / 10;
		if(x >= 1){
			count++;
		}
		else {
			stop = 1;
		}
	}
	return count;
}

template<class T>
int demical(T x)
{
	int stop = 0;
	int demical;
	int count = 0;

	while (!stop){
		demical = (int)x;
		x = x / 10;
		if(x >= 1){
			count++;
		}
		else {
			stop = 1;
		}
	}
	return demical;
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


template<class T, class S>
int distinct(T* arr, S len)
{
    std::set<T> set;
    for (auto i = 0; i < len; i++) {
        set.insert(arr[i]);
    }
    return set.size(); 
}


template<class T, class S>
int distinct(T* arr, S len, T val)
{
    
    std::set<T> set;
    for (auto i = 0; i < len; i++) {
    	if(arr[i] != val){
            set.insert(arr[i]);
	}
    }

    return set.size(); 
}



#endif
