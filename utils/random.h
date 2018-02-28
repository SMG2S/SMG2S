//#include <random>
#include <iostream>
#include <ctime>
#include <stdlib.h>

template<class T>
T random_unint(T min, T max)
{
/*
	std::default_random_engine e(time(NULL));
	std::uniform_real_distribution<T> u(min, max);
*/	
	return (min + static_cast<T>(max * rand() / static_cast<T>(RAND_MAX + 1)));

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
