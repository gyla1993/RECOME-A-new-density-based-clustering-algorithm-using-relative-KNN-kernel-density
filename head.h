#ifndef _HEAD_
#define _HEAD_


#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <algorithm>
#include <utility>
#include <queue>
#include <cmath>
#include <cstring>
#include <map>
#include <sstream>
#include <queue>
#include <vector>
#include <set>
#include <ctime>
#pragma comment(linker, "/STACK:62400000,62400000")


template<typename T>
inline bool checkmax(T &a, const T &b);

template<typename T>
inline bool checkmin(T &a, const T &b);

template<typename T>
inline T ABS(T a);

const int MAXN(21000);  //max #samples
const int MAXK(150);	//max K for KNN
const int MAXD(18);		//max #dimensions
const double LIM(1e30);
const double EPS(1e-8);
const double PI(acos(-1.0));
const double EC(2.7182818284);


//declearation of edge structure
struct E {
	int u, v;
	double w;
	E *next;
	bool operator < (const E &t);
};
//declearation of graph structure
struct G {
	E **h, *e, *re;
	void init(int n, int en);
	void add(int u, int v, double w = 0);
	void destroy();
};

//FJDD algorithm, which can extract all jump discontinuity values of parameter alpha, stored in alphaArr
void getAlpha(int n, G &sg, std::vector<double> &alphaArr, double *rhoStar);

#endif 
