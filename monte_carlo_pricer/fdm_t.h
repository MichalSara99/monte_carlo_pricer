#pragma once
#if !defined(_FDM_H_)
#define _FDM_T_H_

#include"sde_builder.h"
#include"fdm.h"
#include<iostream>
#include<chrono>


using namespace finite_difference_method;
using namespace sde_builder;

void fdm_gbm() {

	double r{ 0.05 };
	double sigma{ 0.01 };
	double s{ 100.0 };

	GeometricBrownianMotion<> gbm{ r,sigma,s };
	std::cout << "Number of factors: " << GeometricBrownianMotion<>::FactorCount << "\n";
	auto sde = gbm.model();
	auto diffusion = gbm.diffusion();
	auto drift = gbm.drift();

	constexpr std::size_t factors = GeometricBrownianMotion<>::FactorCount;
	Fdm<factors,double> gbm_fdm{ sde,1.0,720 };
	auto times = gbm_fdm.timeResolution();

	std::cout << "timing: \n";
	auto start = std::chrono::system_clock::now();
	auto paths_euler = gbm_fdm(50000);
	auto end = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
	std::cout << "Euler took: " << end << " seconds\n";
	start = std::chrono::system_clock::now();
	auto paths_milstein = gbm_fdm(50000,FDMScheme::MilsteinScheme);
	end = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
	std::cout << "Milstein took: " << end << " seconds\n";
	std::cout << "\n";

	// To see few generated values:
	for (std::size_t t = 0; t < 30; ++t) {
		std::cout << t << " paths: \n";
		for (std::size_t v = 0; v < 10; ++v) {
			std::cout << paths_euler[t][v] << ", ";
		}
		std::cout << "\n";
	}

}

void fdm_heston() {

	float r_d{ 0.05f };
	float r_f{ 0.01f };
	float mu = r_d - r_f;
	float sigma{ 0.01f };
	float theta{ 0.015f };
	float kappa{ 0.12f };
	float etha{ 0.012f };
	float stock_init{ 100.0f };
	float var_init{ 0.025f };


	HestonModel<float> hest{ mu,sigma,kappa,theta,etha,stock_init,var_init };
	std::cout << "Number of factors: " << HestonModel<>::FactorCount << "\n";

	constexpr std::size_t factors = HestonModel<>::FactorCount;
	Fdm<HestonModel<>::FactorCount, float> heston_fdm{ hest.model(),1.0,0.8f };
	auto times = heston_fdm.timeResolution();

	std::cout << "timing: \n";
	auto start = std::chrono::system_clock::now();
	auto paths_euler = heston_fdm(50000);
	auto end = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
	std::cout << "Euler took: " << end << " seconds\n";
	start = std::chrono::system_clock::now();
	std::cout << "\n";

	// To see few generated values:
	for (std::size_t t = 0; t < 30; ++t) {
		std::cout << t << " paths: \n";
		for (std::size_t v = 0; v < 10; ++v) {
			std::cout << paths_euler[t][v] << ", ";
		}
		std::cout << "\n";
	}
}





#endif ///_FDM_T_H_

