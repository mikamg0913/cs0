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

// �Ȃ�, ������,0
double rastrigin(const std::valarray<double>& x);//0

											  //�Ȃ�, ������,-418.98*D
double schwefel(const std::valarray<double>& x);//1

											 //����, �P����,0
double rosenbrock(const std::valarray<double>& x);//2

											   //����, ������,0
double griewank(const std::valarray<double>& x);//3

											 //����, �P����,0
double ridge(const std::valarray<double>& xx);//4

										   //�������֐��B���I�œK���̎��ӂɑ����̋Ǐ��������B,0
double ackley(const std::valarray<double>& xx);//5

											//��{�I�Ȋ֐�,0
double sphere(const std::valarray<double>& xx);//6

											//2�ϐ�,-1
double easom(const std::valarray<double>& xx);//7

										   //�Ȃ�,����,0
double xin_she_yang(const std::valarray<double>& xx);//8

												  //����,�����ɂ��œK���̈ʒu���Ⴄ,(2,-1.8013)(5,-4.687658),(10,-9.66015)
double michalewicz(const std::valarray<double>& xx);//9

												 //2�ϐ�,����,����,-1
double eggholder(const std::valarray<double>& xx);//10


