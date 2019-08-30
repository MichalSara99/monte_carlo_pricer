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
	auto sde = gbm.model();
	auto diffusion = gbm.diffusion();
	auto drift = gbm.drift();

	std::cout << "drift(101.1,0.1): " << drift(101.1, 0.1) << "\n";
	std::cout << "diffusion(101.1,0.1): " << diffusion(101.0, 0.1) << "\n";
	std::cout << "====================================\n";
	std::cout << "drift(101.1,0.1): " << sde->drift(101.1, 0.1) << "\n";
	std::cout << "diffusion(101.1,0.1): " << sde->diffusion(101.1, 0.1) << "\n";
}


void cev() {
	float r{ 0.05f };
	float sigma{ 0.01f };
	float s{ 100.0f };
	float beta{ 0.25f };

	ConstantElasticityVariance<float> cev{ r,sigma,beta,s };
	auto sde = cev.model();
	auto diffusion = cev.diffusion();
	auto drift = cev.drift();

	std::cout << "drift(101.1,0.1): " << drift(101.1, 0.1) << "\n";
	std::cout << "diffusion(101.1,0.1): " << diffusion(101.0, 0.1) << "\n";
	std::cout << "====================================\n";
	std::cout << "drift(101.1,0.1): " << sde->drift(101.1, 0.1) << "\n";
	std::cout << "diffusion(101.1,0.1): " << sde->diffusion(101.1, 0.1) << "\n";
}







#endif ///_SDE_T_H_