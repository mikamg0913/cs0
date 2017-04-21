#include <iostream>
#include <valarray>
#include <algorithm>
#include "Cuckoo.hpp"

void update(int N, int SN, std::mt19937& rng, const Func& _func, int _Loop, int tryy);

int main(void) {

	int N;//次元数
	int try_num = 10;//試行回数
	int loop_num = 20;//ループ回数
	int SN;//
	int func_num = 0;//使う関数

	std::cout << "関数:";
	std::cin >> func_num;
	std::cout << "次元数:";
	std::cin >> N;
	std::cout << "個体数:";
	std::cin >> SN;
	std::cout << "イテレーション:";
	std::cin >> loop_num;
	std::cout << "試行回数:";
	std::cin >> try_num;

	FuncSet fset;
	std::random_device rnd;
	std::mt19937 rng(rnd());

	update(N, SN, rng, fset[func_num], loop_num, try_num);
	/*Cuckoo cs;
	cs.move_test();*/

	std::cout << "end";
	return 0;
}

void update(int N, int SN, std::mt19937& rng, const Func& _func, int _Loop, int tryy) {
	std::valarray<double> avgres(_Loop);
	std::valarray<double> bestres(std::numeric_limits<double>::infinity(), _Loop);
	std::valarray<double> worstres(-std::numeric_limits<double>::infinity(), _Loop);
	for (int i = 0; i < tryy; i++) {
		Cuckoo cs(N, SN, rng, _func);
		for (int i = 0; i < _Loop; i++) {
			cs.update();
			avgres[i] += cs.best_f();
			bestres[i] = std::min(bestres[i], cs.best_f());
			worstres[i] = std::max(worstres[i], cs.best_f());
		}
		std::cout << cs.best_f() << std::endl;
	}
	std::ofstream out("out.txt");
	avgres /= tryy;
	for (int i = 0; i < _Loop; i++) {
		out << i + 1 << " " << avgres[i] << " " << bestres[i] << " " << worstres[i] << std::endl;
	}
}
