#include "Cuckoo.hpp"
#include <algorithm>

const double Cuckoo::alpha = 1;
const double Cuckoo::pa = 0.25;

Levy::Levy() : Levy(0.5) { }

Levy::Levy(double _beta) {
	roop_num = 1;
	beta = _beta;
	sigma = pow(tgamma(beta + 1) * sin(std::_Pi * beta / 2) /
		tgamma((beta + 1) / 2) * beta * pow(2, (beta - 1) / 2), 1 / beta);

	n_d_u = std::normal_distribution<>(0, sigma);
	n_d_v = std::normal_distribution<>(0, 1.0);
}

Cuckoo::Cuckoo(int N, int SN, std::mt19937& rng, const Func& funcc) : rng(rng), func(funcc) {
	//levy = Levy(betaa);
	roop_num = 0;

	xx.resize(SN);
	for (auto&& it : xx) {
		it.Initialize(funcc, rng, N);
	}
	

	//wnestÇÃèâä˙âª
	wnest.resize(SN);
	for (auto&& it : wnest) {
		it.resize(N);
	}

	sh_num.resize(SN);
	std::iota(sh_num.begin(), sh_num.end(), 0);
	sh_num2 = sh_num;

	sort_update();
}

void Cuckoo::x_update() {
	int x_rand = 1;
	//int x_rand = data->rd_make_int(sn);
	int j_rand = 1;
	//int j_rand = data->rd_make_int(sn);
	Bird dammy;
	std::valarray<double> lf(xx[x_rand].x.size());
	for (auto& it : lf) {
		//it = 0.1;
		it = Cuckoo::alpha * levy.Make(rng);
	}
	dammy.update(xx[x_rand], lf, func);
	if (xx[j_rand].f() > dammy.f()) {
		xx[j_rand] = dammy;
	}
}

void Cuckoo::x_update_s() {
	x_update_s(std::uniform_int_distribution<int>(0, xx.size())(rng));
}

void Cuckoo::x_update_s(int s) {
	//int x_rand = 1;
	//int x_rand = data->rd_make_int(sn);
	Bird dammy;
	std::valarray<double> lf(xx[s].x.size());
	for (unsigned i = 0; i < xx[s].x.size(); ++i) {
		//lf[i] = 0.1;
		lf[i] = levy.Make_s(rng, xx[s].x[i], best.x[i]);
		//lf[i] = levy.Make();
	}
	dammy.update(xx[s], lf, func);
	if (xx[s].f() > dammy.f()) {
		xx[s] = dammy;
	}
}

void Cuckoo::sort_update() {
	best = *std::min_element(xx.cbegin(), xx.cend(), [](const Bird& a, const Bird& b) { return a.f() < b.f(); });
}

void Cuckoo::worst_update() {
	/*for (int i = 0; i < sn*pa; i++)
	xx.pop_back();
	for (int i = 0; i < sn*pa; i++) {
	xx.push_back(x_make());
	}*/
	//wnestê∂ê¨
	for (unsigned i = 0; i < wnest.size(); i++) {
		for (unsigned j = 0; j < wnest[i].size(); j++) {
			wnest[i][j] = std::uniform_real_distribution<double>(0, 1)(rng) >= pa;
		}
	}

	//F
	double _F = std::uniform_real_distribution<double>(0, 1)(rng);

	auto x_dm = xx;

	std::shuffle(sh_num.begin(), sh_num.end(), rng);
	std::shuffle(sh_num2.begin(), sh_num2.end(), rng);

	//dm=nest+F*wnestÅõ(nest_per-nest_per)
	for (unsigned j = 0; j < xx.size(); j++) {
		Bird dammy;
		dammy.worst_update(xx[j], _F, wnest[j], x_dm[sh_num[j]], x_dm[sh_num2[j]], func);
		//std::cout << "å≥:" << xx[j].f << " dm:" << dammy.f << std::endl;
		if (xx[j].f() > dammy.f()) {
			//xx[j] = dammy;
			xx[j] = std::move(dammy);
		}
	}
}

void Cuckoo::update() {
	for (unsigned i = 0; i < xx.size(); i++) {
		//x_update();
		x_update_s(i);
	}
	//x_update();
	sort_update();
	worst_update();
	roop_num++;
	levy.roop_num = roop_num;
}
