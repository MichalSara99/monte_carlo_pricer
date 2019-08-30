#pragma once
#if !defined(_PAYOFF_T_H_)
#define _PAYOFF_T_H_

#include"payoff.h"
#include<random>
using namespace payoff;
using namespace std::placeholders;

void payoff1() {
	// Binding the already existing strategy:
	double call_strike{ 20.0 };
	PlainCallStrategy call_s{ call_strike };
	auto call_payoff = std::bind(&PlainCallStrategy::payoff, call_s, _1);
	Payoff<double> pay_call{ call_payoff };
	std::cout << "Enter price of underlying: ";
	double S;
	std::cin >> S;
	std::cout << "Plain call payoff ( K = " << call_strike << "): " << pay_call.payoff(S);
	std::cout << "\n";
	std::cout << "Template type: " << typeid(decltype(pay_call)::underlyingType).name() << "\n";
}

void payoff2() {
	// Binding the already existing strategy:
	std::vector<double> prices(static_cast<std::size_t>(1000));
	std::random_device rd;
	std::normal_distribution<double> normals;
	double init{ 20.0 };
	std::generate(prices.begin(), prices.end(), [&]() {return init * std::exp(normals(rd)); });
	double asian_strike{ 21.0 };
	AsianAvgCallStrategy asian_s{ asian_strike };
	auto call_payoff = std::bind(&AsianAvgCallStrategy::payoff, asian_s, _1);
	auto asian_pay = Payoff<PathValuesType<double>>{ call_payoff };
	std::cout << "Average asian call payoff ( K = " << asian_strike << "): " << asian_pay.payoff(prices);
	std::cout << "\n";
	std::cout << "Template type: " << typeid(decltype(asian_pay)::underlyingType).name() << "\n";
}

void payoff3() {
	// Create the strategy on the fly:
	std::vector<double> prices(static_cast<std::size_t>(10));
	std::random_device rd;
	std::normal_distribution<double> normals;
	double init{ 100.0 };
	std::generate(prices.begin(), prices.end(), [&]() {return init * std::exp(normals(rd)); });
	double asian_strike{ 1.23 };

	auto asianGeoCallStrategy = [&](PathValuesType<double> const &underlying) {
		std::size_t N{ underlying.size() };
		auto logsum = std::accumulate(underlying.begin(), underlying.end(), 1.0, [](double x, double y) {
			return (std::log(x) + std::log(y));
		});
		auto geo = std::exp(logsum / static_cast<double>(N));
		return std::max(0.0, geo - asian_strike);
	};

	Payoff<PathValuesType<double>> asian_pay{ asianGeoCallStrategy };
	std::cout << "Geometric asian call payoff ( K = " << asian_strike << "): " << asian_pay.payoff(prices);
	std::cout << "\n";
	std::cout << "Template type: " << typeid(decltype(asian_pay)::underlyingType).name() << "\n";


}

#endif //_PAYOFF_T_H_