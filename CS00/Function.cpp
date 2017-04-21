#include "Function.hpp"
#define _USE_MATH_DEFINES
#include <math.h>

FuncSet::FuncSet() {
	ff.emplace_back(rastrigin,           -5.12,         5.12 );
	ff.emplace_back(schwefel,          -512,          512    );
	ff.emplace_back(rosenbrock,          -2.048,        2.048);
	ff.emplace_back(griewank,          -512,          512    );
	ff.emplace_back(ridge,              -64,           64    );
	ff.emplace_back(ackley,             -32.768,       32.768);
	ff.emplace_back(sphere,       -10000000,     10000000    );
	ff.emplace_back(easom,             -100,          100    );
	ff.emplace_back(xin_she_yang, -2 * M_PI,     2 * M_PI    );
	ff.emplace_back(michalewicz,          0,         M_PI    );
	ff.emplace_back(eggholder,         -512,          512    );
}

double fit_calc(double f) {
	if (f < 0)
		return 1 + fabs(f);
	else
		return 1 / (1 + fabs(f));
}

double rastrigin(const std::valarray<double>& xx) {
	int A = 10, n = 0;
	double ans = 0;
	for (auto x : xx) {
		ans += x * x - A * cos(2 * M_PI*x);
		n++;
	}
	ans += A * n;
	return ans;
}

double schwefel(const std::valarray<double>& xx) {
	double ans = 0;
	int i = 0;
	for (auto x : xx) {
		ans -= x*sin(sqrt(fabs(x)));
		i++;
	}
	return ans + i*418.9829;
}

double rosenbrock(const std::valarray<double>& xx) {
	double ans = 0;
	for (unsigned i = 0; i < xx.size() - 1; i++) {
		auto s = xx[i + 1] - xx[i] * xx[i];
		ans += 100 * s * s + (1 - xx[i]) * (1 - xx[i]);
	}
	return ans;
}

double griewank(const std::valarray<double>& xx) {
	double ans = 1, dm = 1;
	int i = 1;
	for (auto x : xx) {
		ans += x * x / 4000.0;
		dm *= cos(x / sqrt(i));
	}
	return ans - dm;
}

double ridge(const std::valarray<double>& xx) {
	double ans = 0;
	for (unsigned i = 1; i <= xx.size(); i++) {
		for (int j = 1; j <= i; j++) {
			ans += pow(xx[j - 1], 2);
		}
	}
	return ans;
}

double ackley(const std::valarray<double>& xx) {
	double ans = 20, x1 = 0, x2 = 0;
	int n = 0;
	for (auto x : xx) {
		x1 += x*x;
		x2 += cos(2 * M_PI*x);
		n++;
	}
	return ans - 20 * exp(-0.2*pow(x1 / n*1.0, 0.5)) + exp(1) - exp(x2 / n);
}

double sphere(const std::valarray<double>& xx) {
	double ans = 0;
	for (auto x : xx) {
		ans += x*x;
	}
	return ans;
}

double easom(const std::valarray<double>& xx) {
	double ans = (-1.0)*cos(xx[0])*cos(xx[1])*exp(-pow(xx[0] - M_PI, 2.0) - pow(xx[1] - M_PI, 2.0));
	return ans;
}

double xin_she_yang(const std::valarray<double>& xx) {
	double ans = 0, x0 = 0;
	for (auto x : xx) {
		ans += abs(x);
		x0 += sin(x*x);
	}
	ans *= exp(-x0);
	return ans;
}

double michalewicz(const std::valarray<double>& xx) {
	double x0 = 0;
	int i = 1;
	for (auto x : xx) {
		x0 += sin(x)*pow(sin(i*x*x / M_PI), 20);
		i++;
	}
	return -x0;
}

double eggholder(const std::valarray<double>& xx) {

	return -(xx[1] + 47)*sin(pow(abs(xx[1] + xx[0] / 2.0 + 47), 0.5)) - xx
		[0] * sin(pow(abs(xx[0] - (xx[1] + 47)), 0.5));
}



