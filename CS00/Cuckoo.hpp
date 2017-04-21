#pragma once
#include "Function.hpp"
#include <random>
#include <functional>
#include <fstream>
#include <valarray>
#include <numeric>

class Levy {
private:
	double beta, sigma;
	std::normal_distribution<> n_d_u, n_d_v;

public:
	int roop_num;
	Levy();
	Levy(double beta);
	template <class Engine> double Make(Engine& engine) {
		//return pow(roop_num, -beta);
		//return (n_d_u(eng)*sigma)/ pow(abs(n_d_v(eng)), 1 / beta);
		return (n_d_u(engine)*sigma) / pow(abs(n_d_v(engine)), 1 / beta);
		//return n_d_u(eng);
		return 0.1;
	}
	template <class Engine> double Make_s(Engine& engine, double x, double best) {
		double step = 0.01*Make(engine)*(best - x);
		//return step*Make();
		//return step;
		return step*n_d_v(engine);
		//return 0.1;
	}
};

class Bird {
public:
	Bird() { }
	Bird(const Bird& src) : x(src.x), fval(src.fval) { }
	Bird(Bird&& src) : x(std::move(src.x)), fval(src.fval) { }

	std::valarray<double> x;

	Bird& operator =(const Bird& src) {
		x = src.x;
		fval = src.fval;
		return *this;
	}
	Bird& operator =(Bird&& src) {
		x = std::move(src.x);
		fval = src.fval;
		return *this;
	}
	void Initialize(const Func& func, std::mt19937& rng, int dimension) {
		x.resize(dimension);
		for (auto& it : x) {
			it = std::uniform_real_distribution<double>(func.min, func.max)(rng);
		}
		fval = func(x);
	}
	double f() const { return fval; }
	void update(const Bird& source, const std::valarray<double> levyFlight, const Func& func) {
		x = wrap_values(source.x + levyFlight, func);
		fval = func(x);
	}
	void worst_update(const Bird& source, double _F, const std::valarray<bool>& wnest, const Bird& r1, const Bird& r2, const Func& func) {
		x = wrap_values(source.x + _F * mask_value(wnest, r1.x - r2.x), func);
		fval = func(x);
	}

private:
	double fval;

	template <class T> static std::valarray<T> mask_value(const std::valarray<bool>& mask, const std::valarray<T>& target) {
		std::valarray<T> result(mask.size());
		for (size_t i = 0; i < result.size(); ++i) {
			result[i] = mask[i] ? target[i] : static_cast<T>(0);
		}
		return result;
	}
	//折り返し
	static double wrap_value(double xi, const Func& func) {
		auto min = func.min;
		auto max = func.max;
		if (xi <= max && xi >= min) {
			return xi;
		}
		auto am = (xi - min) - (int)((xi - min) / (max - min))*(max - min);
		return (xi > max ? max : min) - am;
	}
	static std::valarray<double> wrap_values(const std::valarray<double>& x, const Func& func) {
		std::valarray<double> result(x.size());
		for (size_t i = 0; i < result.size(); ++i) {
			result[i] = wrap_value(x[i], func);
		}
		return result;
	}
};

class Cuckoo {
private:
	static const double alpha;
	static const double pa;
	const Func& func;
	int roop_num;
	std::mt19937& rng;

	//値
	std::vector<Bird> xx;
	Bird best;

	std::vector<std::valarray<bool>> wnest;
	std::vector<int> sh_num, sh_num2;

	Levy levy{ 1.5 };
	//Levy levy;

public:
	double best_f() const { return best.f(); }

	Cuckoo(int N, int SN, std::mt19937& rng, const Func& funcc);
	Cuckoo(const Cuckoo&) = delete;
	Cuckoo& operator =(const Cuckoo&) = delete;

	//ランダムウォーク含めた更新
	void x_update();
	void x_update_s();
	void x_update_s(int);

	//ソートするだけ
	void sort_update();
	//悪巣を排除、新しいの作る
	void worst_update();

	//全体の更新
	void update();
};