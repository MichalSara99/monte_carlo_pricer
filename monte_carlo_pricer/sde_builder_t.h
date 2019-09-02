#pragma once
#if !defined(_SDE_T_H_)
#define _SDE_T_H_

#include"sde_builder.h"
#include<iostream>

using namespace sde_builder;


void gbm() {
	double r{ 0.05 };
	double sigma{ 0.01 };
	double s{ 100.0 };

	GeometricBrownianMotion<> gbm{ r,sigma,s };
	std::cout << "Number of factors: " << GeometricBrownianMotion<>::FactorCount << "\n";
	auto sde = gbm.model();
	auto diffusion = gbm.diffusion();
	auto drift = gbm.drift();

	std::cout << "drift(101.1,0.1): " << drift(0.1,101.1) << "\n";
	std::cout << "diffusion(101.1,0.1): " << diffusion(0.1,101.0) << "\n";
	std::cout << "====================================\n";
	std::cout << "drift(101.1,0.1): " << sde->drift(0.1,101.1) << "\n";
	std::cout << "diffusion(101.1,0.1): " << sde->diffusion(0.1,101.1) << "\n";
}


void cev() {
	float r{ 0.05f };
	float sigma{ 0.01f };
	float s{ 100.0f };
	float beta{ 0.25f };

	ConstantElasticityVariance<float> cev{ r,sigma,beta,s };
	std::cout << "Number of factors: " << ConstantElasticityVariance<>::FactorCount << "\n";
	auto sde = cev.model();
	auto diffusion = cev.diffusion();
	auto drift = cev.drift();

	std::cout << "drift(101.1,0.1): " << drift(0.1,101.1) << "\n";
	std::cout << "diffusion(101.1,0.1): " << diffusion(0.1,101.0) << "\n";
	std::cout << "====================================\n";
	std::cout << "drift(101.1,0.1): " << sde->drift(0.1,101.1) << "\n";
	std::cout << "diffusion(101.1,0.1): " << sde->diffusion(0.1,101.1) << "\n";
}


void heston() {
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
	auto stock_sde = std::get<0>(hest.model());
	auto var_sde = std::get<1>(hest.model());
	auto diffusion1 = hest.diffusion1();
	auto drift1 = hest.drift1();
	auto diffusion2 = hest.diffusion2();
	auto drift2 = hest.drift2();

	//std::cout << "drift1(101.1,0.1): " << drift1(101.1, 0.1) << "\n";
	//std::cout << "diffusion1(101.1,0.1): " << diffusion1(101.0, 0.1) << "\n";
	//std::cout << "drift1(101.1,0.1): " << drift2(101.1, 0.1) << "\n";
	//std::cout << "diffusion1(101.1,0.1): " << diffusion2(101.0, 0.1) << "\n";
	std::cout << "====================================\n";
	std::cout << "Number of factors: " << HestonModel<>::FactorCount << "\n";
	std::cout << "drift1(101.1,0.012,0.1): " << stock_sde->drift(0.1,101.1, 0.012) << "\n";
	std::cout << "diffusion1(101.1,0.012,0.1): " << stock_sde->diffusion(0.1,101.1,0.012) << "\n";
	std::cout << "drift2(101.1,0.012,0.1): " << var_sde->drift(0.1,101.1, 0.012) << "\n";
	std::cout << "diffusion2(101.1,0.012,0.1): " << var_sde->diffusion(0.1,101.1, 0.012) << "\n";
}






#endif ///_SDE_T_H_