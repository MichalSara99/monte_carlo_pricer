#pragma once
#if !defined(_EXAMPLES_H_)
#define _EXAMPLES_H_

#include<iostream>
#include<string>
#include"payoff.h"
#include"payoff_strategy.h"
#include"fdm.h"
#include"sde_builder.h"

using namespace finite_difference_method;
using namespace sde_builder;

// Pricing european options 
// using paths from geometric brownian motion  
void europeanOptionsGBMEuler() {
	
	// First generate paths using GBM 
	double rate{ 0.001 };
	double sigma{ 0.005 };
	double s{ 100.0 };
	double maturityInYears{ 1.0 };
	std::size_t numberSteps{ 720 }; // two times a day
	std::size_t simuls{ 70'000 };

	// Construct the model:
	GeometricBrownianMotion<> gbm{ rate,sigma,s };
	std::cout << "Model: " << gbm.name() << "\n";
	std::cout << "Factors: "<< GeometricBrownianMotion<>::FactorCount << "\n";
	// Construct the engine: 
	Fdm<GeometricBrownianMotion<>::FactorCount, double> fdm_gbm{ gbm.model(),maturityInYears,numberSteps };
	auto start = std::chrono::system_clock::now();
	auto paths_euler = fdm_gbm(simuls, FDMScheme::EulerScheme);
	auto end = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
	std::cout << "Euler scheme for GBM<1> took: " << end << " seconds.\n";

	// last values:
	for (std::size_t t = 0; t < 30; ++t) {
		std::cout << t << " paths: \n";
		std::cout << paths_euler[t][719] << ", ";
		std::cout << "\n";
	}
	std::cout << "==============================\n";


	// Construct call payoff of the option:
	double call_strike{ 100.0 };
	PlainCallStrategy<> call_strategy{ call_strike };
	auto call_payoff = std::bind(&PlainCallStrategy<>::payoff, &call_strategy, std::placeholders::_1);
	// Construct put payoff of the option:
	double put_strike{ 100.0 };
	PlainPutStrategy<> put_strategy{ put_strike };
	auto put_payoff = std::bind(&PlainPutStrategy<>::payoff, &put_strategy, std::placeholders::_1);

	double lastPrice{0.0};
	double call_sum{0.0};
	double put_sum{0.0};
	for (auto const &path : paths_euler) {
		lastPrice = path.back();
		call_sum += call_payoff(lastPrice);
		put_sum += put_payoff(lastPrice);
	}
	auto call_avg = (call_sum / static_cast<double>(paths_euler.size()));
	auto put_avg = (put_sum / static_cast<double>(paths_euler.size()));

	std::cout << "Call price: " << ((std::exp(-1.0*rate*maturityInYears))*call_avg) << "\n";
	std::cout << "Put price: " << ((std::exp(-1.0*rate*maturityInYears))*put_avg) << "\n";
	std::cout << "=========================================================\n";
}

// Pricing european options 
// using paths from geometric brownian motion  
void europeanOptionsGBMMilstein() {

	// First generate paths using GBM 
	double rate{ 0.001 };
	double sigma{ 0.005 };
	double s{ 100.0 };
	double maturityInYears{ 1.0 };
	std::size_t numberSteps{ 720 }; // two times a day
	std::size_t simuls{ 50000 };

	// Construct the model:
	GeometricBrownianMotion<> gbm{ rate,sigma,s };
	std::cout << "Model: " << gbm.name() << "\n";
	std::cout << "Factors: " << GeometricBrownianMotion<>::FactorCount << "\n";
	// Construct the engine: 
	Fdm<GeometricBrownianMotion<>::FactorCount, double> fdm_gbm{ gbm.model(),maturityInYears,numberSteps };
	auto start = std::chrono::system_clock::now();
	auto paths_euler = fdm_gbm(simuls, FDMScheme::MilsteinScheme);
	auto end = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
	std::cout << "Milstein scheme for GBM<1> took: " << end << " seconds.\n";

	// Construct call payoff of the option:
	double call_strike{ 100.0 };
	PlainCallStrategy<> call_strategy{ call_strike };
	auto call_payoff = std::bind(&PlainCallStrategy<>::payoff, &call_strategy, std::placeholders::_1);
	// Construct put payoff of the option:
	double put_strike{ 100.0 };
	PlainPutStrategy<> put_strategy{ put_strike };
	auto put_payoff = std::bind(&PlainPutStrategy<>::payoff, &put_strategy, std::placeholders::_1);

	double lastPrice{ 0.0 };
	double call_sum{ 0.0 };
	double put_sum{ 0.0 };
	for (auto const &path : paths_euler) {
		lastPrice = path.back();
		call_sum += call_payoff(lastPrice);
		put_sum += put_payoff(lastPrice);
	}
	auto call_avg = (call_sum / static_cast<double>(paths_euler.size()));
	auto put_avg = (put_sum / static_cast<double>(paths_euler.size()));

	std::cout << "Call price: " << ((std::exp(-1.0*rate*maturityInYears))*call_avg) << "\n";
	std::cout << "Put price: " << ((std::exp(-1.0*rate*maturityInYears))*put_avg) << "\n";
	std::cout << "=========================================================\n";
}

// Pricing european options 
// using paths from geometric brownian motion  
void europeanOptionsHestonEuler() {

	// First generate paths using Heston Model:
	double rate_d{ 0.005 };
	double rate_f{ 0.004 };
	double mu{ rate_d - rate_f };
	double sigma{ 0.005 };
	double theta{ 0.0155 };
	double kappa{ 0.14 };
	double etha{ 0.012 };
	double s_init{ 100.0 };
	double var_init{ 0.025 };
	double correlation{ 0.0 };
	double maturityInYears{ 1.0 };
	std::size_t numberSteps{ 720 }; // two times a day
	std::size_t simuls{ 50000 };

	// Construct the model:
	HestonModel<> heston{ mu,sigma,kappa,theta,etha,s_init,var_init };
	std::cout << "Model: " << heston.name() << "\n";
	std::cout << "Factors: " << HestonModel<>::FactorCount << "\n";
	// Construct the engine: 
	Fdm<HestonModel<>::FactorCount, double> fdm_heston{ heston.model(),maturityInYears,correlation,numberSteps };
	auto start = std::chrono::system_clock::now();
	auto paths_euler = fdm_heston(simuls, FDMScheme::EulerScheme);
	auto end = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
	std::cout << "Euler scheme for Heston<2> took: " << end << " seconds.\n";

	// Construct call payoff of the option:
	double call_strike{ 100.0 };
	PlainCallStrategy<> call_strategy{ call_strike };
	auto call_payoff = std::bind(&PlainCallStrategy<>::payoff, &call_strategy, std::placeholders::_1);
	// Construct put payoff of the option:
	double put_strike{ 100.0 };
	PlainPutStrategy<> put_strategy{ put_strike };
	auto put_payoff = std::bind(&PlainPutStrategy<>::payoff, &put_strategy, std::placeholders::_1);

	double lastPrice{ 0.0 };
	double call_sum{ 0.0 };
	double put_sum{ 0.0 };
	for (auto const &path : paths_euler) {
		lastPrice = path.back();
		call_sum += call_payoff(lastPrice);
		put_sum += put_payoff(lastPrice);
	}
	auto call_avg = (call_sum / static_cast<double>(paths_euler.size()));
	auto put_avg = (put_sum / static_cast<double>(paths_euler.size()));

	std::cout << "Call price: " << ((std::exp(-1.0*mu*maturityInYears))*call_avg) << "\n";
	std::cout << "Put price: " << ((std::exp(-1.0*mu*maturityInYears))*put_avg) << "\n";
	std::cout << "=========================================================\n";
}


// Pricing european options 
// using paths from geometric brownian motion  
void europeanOptionsHestonMilstein() {

	// First generate paths using Heston Model:
	double rate_d{ 0.005 };
	double rate_f{ 0.004 };
	double mu{ rate_d - rate_f };
	double sigma{ 0.005 };
	double theta{ 0.0155 };
	double kappa{ 0.14 };
	double etha{ 0.012 };
	double s_init{ 100.0 };
	double var_init{ 0.025 };
	double correlation{ 0.0 };
	double maturityInYears{ 1.0 };
	std::size_t numberSteps{ 720 }; // two times a day
	std::size_t simuls{ 50000 };

	// Construct the model:
	HestonModel<> heston{ mu,sigma,kappa,theta,etha,s_init,var_init };
	std::cout << "Model: " << heston.name() << "\n";
	std::cout << "Factors: " << HestonModel<>::FactorCount << "\n";
	// Construct the engine: 
	Fdm<HestonModel<>::FactorCount, double> fdm_heston{ heston.model(),maturityInYears,correlation,numberSteps };
	auto start = std::chrono::system_clock::now();
	auto paths_euler = fdm_heston(simuls, FDMScheme::MilsteinScheme);
	auto end = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
	std::cout << "Milstein scheme for Heston<2> took: " << end << " seconds.\n";

	// Construct call payoff of the option:
	double call_strike{ 100.0 };
	PlainCallStrategy<> call_strategy{ call_strike };
	auto call_payoff = std::bind(&PlainCallStrategy<>::payoff, &call_strategy, std::placeholders::_1);
	// Construct put payoff of the option:
	double put_strike{ 100.0 };
	PlainPutStrategy<> put_strategy{ put_strike };
	auto put_payoff = std::bind(&PlainPutStrategy<>::payoff, &put_strategy, std::placeholders::_1);

	double lastPrice{ 0.0 };
	double call_sum{ 0.0 };
	double put_sum{ 0.0 };
	for (auto const &path : paths_euler) {
		lastPrice = path.back();
		call_sum += call_payoff(lastPrice);
		put_sum += put_payoff(lastPrice);
	}
	auto call_avg = (call_sum / static_cast<double>(paths_euler.size()));
	auto put_avg = (put_sum / static_cast<double>(paths_euler.size()));

	std::cout << "Call price: " << ((std::exp(-1.0*mu*maturityInYears))*call_avg) << "\n";
	std::cout << "Put price: " << ((std::exp(-1.0*mu*maturityInYears))*put_avg) << "\n";
	std::cout << "=========================================================\n";
}


// Pricing asian options 
// using paths from geometric brownian motion  
void asianOptionsGBMEuler() {

	// First generate paths using GBM 
	double rate{ 0.001 };
	double sigma{ 0.005 };
	double s{ 100.0 };
	double maturityInYears{ 1.0 };
	std::size_t numberSteps{ 720 }; // two times a day
	std::size_t simuls{ 50000 };

	// Construct the model:
	GeometricBrownianMotion<> gbm{ rate,sigma,s };
	std::cout << "Model: " << gbm.name() << "\n";
	std::cout << "Factors: " << GeometricBrownianMotion<>::FactorCount << "\n";
	// Construct the engine: 
	Fdm<GeometricBrownianMotion<>::FactorCount, double> fdm_gbm{ gbm.model(),maturityInYears,numberSteps };
	auto start = std::chrono::system_clock::now();
	auto paths_euler = fdm_gbm(simuls, FDMScheme::EulerScheme);
	auto end = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
	std::cout << "Euler scheme for GBM<1> took: " << end << " seconds.\n";

	// Construct call payoff of the option:
	double call_strike{ 100.0 };
	AsianAvgCallStrategy<> call_strategy{ call_strike };
	auto call_payoff = std::bind(&AsianAvgCallStrategy<>::payoff, &call_strategy, std::placeholders::_1);
	// Construct put payoff of the option:
	double put_strike{ 100.0 };
	AsianAvgPutStrategy<> put_strategy{ put_strike };
	auto put_payoff = std::bind(&AsianAvgPutStrategy<>::payoff, &put_strategy, std::placeholders::_1);

	double call_sum{ 0.0 };
	double put_sum{ 0.0 };
	for (auto const &path : paths_euler) {
		call_sum += call_payoff(path);
		put_sum += put_payoff(path);
	}
	auto call_avg = (call_sum / static_cast<double>(paths_euler.size()));
	auto put_avg = (put_sum / static_cast<double>(paths_euler.size()));

	std::cout << "Asian call price: " << ((std::exp(-1.0*rate*maturityInYears))*call_avg) << "\n";
	std::cout << "Asian put price: " << ((std::exp(-1.0*rate*maturityInYears))*put_avg) << "\n";
	std::cout << "=========================================================\n";
}


// Pricing european options 
// using paths from geometric brownian motion  
void asianOptionsGBMMilstein() {

	// First generate paths using GBM 
	double rate{ 0.001 };
	double sigma{ 0.005 };
	double s{ 100.0 };
	double maturityInYears{ 1.0 };
	std::size_t numberSteps{ 720 }; // two times a day
	std::size_t simuls{ 50000 };

	// Construct the model:
	GeometricBrownianMotion<> gbm{ rate,sigma,s };
	std::cout << "Model: " << gbm.name() << "\n";
	std::cout << "Factors: " << GeometricBrownianMotion<>::FactorCount << "\n";
	// Construct the engine: 
	Fdm<GeometricBrownianMotion<>::FactorCount, double> fdm_gbm{ gbm.model(),maturityInYears,numberSteps };
	auto start = std::chrono::system_clock::now();
	auto paths_euler = fdm_gbm(simuls, FDMScheme::MilsteinScheme);
	auto end = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
	std::cout << "Milstein scheme for GBM<1> took: " << end << " seconds.\n";

	// Construct call payoff of the option:
	double call_strike{ 100.0 };
	AsianAvgCallStrategy<> call_strategy{ call_strike };
	auto call_payoff = std::bind(&AsianAvgCallStrategy<>::payoff, &call_strategy, std::placeholders::_1);
	// Construct put payoff of the option:
	double put_strike{ 100.0 };
	AsianAvgPutStrategy<> put_strategy{ put_strike };
	auto put_payoff = std::bind(&AsianAvgPutStrategy<>::payoff, &put_strategy, std::placeholders::_1);

	double call_sum{ 0.0 };
	double put_sum{ 0.0 };
	for (auto const &path : paths_euler) {
		call_sum += call_payoff(path);
		put_sum += put_payoff(path);
	}
	auto call_avg = (call_sum / static_cast<double>(paths_euler.size()));
	auto put_avg = (put_sum / static_cast<double>(paths_euler.size()));

	std::cout << "Call price: " << ((std::exp(-1.0*rate*maturityInYears))*call_avg) << "\n";
	std::cout << "Put price: " << ((std::exp(-1.0*rate*maturityInYears))*put_avg) << "\n";
	std::cout << "=========================================================\n";
}

// Pricing asian options 
// using paths from geometric brownian motion  
void asianOptionsHestonEuler() {

	// First generate paths using Heston Model:
	double rate_d{ 0.005 };
	double rate_f{ 0.004 };
	double mu{ rate_d - rate_f };
	double sigma{ 0.005 };
	double theta{ 0.0155 };
	double kappa{ 0.14 };
	double etha{ 0.012 };
	double s_init{ 100.0 };
	double var_init{ 0.025 };
	double correlation{ 0.0 };
	double maturityInYears{ 1.0 };
	std::size_t numberSteps{ 720 }; // two times a day
	std::size_t simuls{ 50000 };

	// Construct the model:
	HestonModel<> heston{ mu,sigma,kappa,theta,etha,s_init,var_init };
	std::cout << "Model: " << heston.name() << "\n";
	std::cout << "Factors: " << HestonModel<>::FactorCount << "\n";
	// Construct the engine: 
	Fdm<HestonModel<>::FactorCount, double> fdm_heston{ heston.model(),maturityInYears,correlation,numberSteps };
	auto start = std::chrono::system_clock::now();
	auto paths_euler = fdm_heston(simuls, FDMScheme::EulerScheme);
	auto end = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
	std::cout << "Euler scheme for Heston<2> took: " << end << " seconds.\n";

	// Construct call payoff of the option:
	double call_strike{ 100.0 };
	AsianAvgCallStrategy<> call_strategy{ call_strike };
	auto call_payoff = std::bind(&AsianAvgCallStrategy<>::payoff, &call_strategy, std::placeholders::_1);
	// Construct put payoff of the option:
	double put_strike{ 100.0 };
	AsianAvgPutStrategy<> put_strategy{ put_strike };
	auto put_payoff = std::bind(&AsianAvgPutStrategy<>::payoff, &put_strategy, std::placeholders::_1);

	double call_sum{ 0.0 };
	double put_sum{ 0.0 };
	for (auto const &path : paths_euler) {
		call_sum += call_payoff(path);
		put_sum += put_payoff(path);
	}
	auto call_avg = (call_sum / static_cast<double>(paths_euler.size()));
	auto put_avg = (put_sum / static_cast<double>(paths_euler.size()));

	std::cout << "Call price: " << ((std::exp(-1.0*mu*maturityInYears))*call_avg) << "\n";
	std::cout << "Put price: " << ((std::exp(-1.0*mu*maturityInYears))*put_avg) << "\n";
	std::cout << "=========================================================\n";
}

// Pricing asian options 
// using paths from geometric brownian motion  
void asianOptionsHestonMilstein() {

	// First generate paths using Heston Model:
	double rate_d{ 0.005 };
	double rate_f{ 0.004 };
	double mu{ rate_d - rate_f };
	double sigma{ 0.005 };
	double theta{ 0.0155 };
	double kappa{ 0.14 };
	double etha{ 0.012 };
	double s_init{ 100.0 };
	double var_init{ 0.025 };
	double correlation{ 0.0 };
	double maturityInYears{ 1.0 };
	std::size_t numberSteps{ 720 }; // two times a day
	std::size_t simuls{ 50000 };

	// Construct the model:
	HestonModel<> heston{ mu,sigma,kappa,theta,etha,s_init,var_init };
	std::cout << "Model: " << heston.name() << "\n";
	std::cout << "Factors: " << HestonModel<>::FactorCount << "\n";
	// Construct the engine: 
	Fdm<HestonModel<>::FactorCount, double> fdm_heston{ heston.model(),maturityInYears,correlation,numberSteps };
	auto start = std::chrono::system_clock::now();
	auto paths_euler = fdm_heston(simuls, FDMScheme::MilsteinScheme);
	auto end = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
	std::cout << "Milstein scheme for Heston<2> took: " << end << " seconds.\n";

	// Construct call payoff of the option:
	double call_strike{ 100.0 };
	AsianAvgCallStrategy<> call_strategy{ call_strike };
	auto call_payoff = std::bind(&AsianAvgCallStrategy<>::payoff, &call_strategy, std::placeholders::_1);
	// Construct put payoff of the option:
	double put_strike{ 100.0 };
	AsianAvgPutStrategy<> put_strategy{ put_strike };
	auto put_payoff = std::bind(&AsianAvgPutStrategy<>::payoff, &put_strategy, std::placeholders::_1);

	double call_sum{ 0.0 };
	double put_sum{ 0.0 };
	for (auto const &path : paths_euler) {
		call_sum += call_payoff(path);
		put_sum += put_payoff(path);
	}
	auto call_avg = (call_sum / static_cast<double>(paths_euler.size()));
	auto put_avg = (put_sum / static_cast<double>(paths_euler.size()));

	std::cout << "Call price: " << ((std::exp(-1.0*mu*maturityInYears))*call_avg) << "\n";
	std::cout << "Put price: " << ((std::exp(-1.0*mu*maturityInYears))*put_avg) << "\n";
	std::cout << "=========================================================\n";
}


#endif ///_EXAMPLES_H_