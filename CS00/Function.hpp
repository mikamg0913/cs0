#pragma once
#include <functional>
#include <valarray>
#include <vector>

class Func {
public:
	template <class T> Func(T f, double min, double max) : f(f), min(min), max(max) { }

	double min;
	double max;
	double operator()(const std::valarray<double>& x) const { return f(x); }

private:
	std::function<double(const std::valarray<double>&)> f;
};

class FuncSet {
public:
	FuncSet();
	const Func& operator[](size_t index) { return ff[index]; }
private:
	std::vector<Func> ff;
};

double fit_calc(double);

// なし, 多峰性,0
double rastrigin(const std::valarray<double>& x);//0

											  //なし, 多峰性,-418.98*D
double schwefel(const std::valarray<double>& x);//1

											 //あり, 単峰性,0
double rosenbrock(const std::valarray<double>& x);//2

											   //あり, 多峰性,0
double griewank(const std::valarray<double>& x);//3

											 //あり, 単峰性,0
double ridge(const std::valarray<double>& xx);//4

										   //多峰性関数。大域的最適解の周辺に多数の局所解をもつ。,0
double ackley(const std::valarray<double>& xx);//5

											//基本的な関数,0
double sphere(const std::valarray<double>& xx);//6

											//2変数,-1
double easom(const std::valarray<double>& xx);//7

										   //なし,多峰,0
double xin_she_yang(const std::valarray<double>& xx);//8

												  //多峰,次元により最適解の位置が違う,(2,-1.8013)(5,-4.687658),(10,-9.66015)
double michalewicz(const std::valarray<double>& xx);//9

												 //2変数,あり,多峰,-1
double eggholder(const std::valarray<double>& xx);//10


