#pragma once
#if !defined(_FDM_H_)
#define _FDM_T_H_

#include"sde_builder.h"
#include"fdm.h"
#include<iostream>
#include<ctime>

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
	auto start = clock();
	auto paths_euler = gbm_fdm(30000);
	auto end = clock() - start;
	std::cout << "Euler took: " << (end / static_cast<float>(CLOCKS_PER_SEC)) << "\n";
	start = clock();
	auto paths_milstein = gbm_fdm(30000,FDMScheme::MilsteinScheme);
	end = clock() - start;
	std::cout << "Milstein took: " << (end / static_cast<float>(CLOCKS_PER_SEC)) << "\n";
	std::cout << "\n";

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

	Fdm<HestonModel<>::FactorCount, float> heston_fdm{ hest.model(),1.0 };


}





#endif ///_FDM_T_H_

