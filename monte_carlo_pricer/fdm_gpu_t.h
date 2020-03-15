#pragma once
#if !defined(_FDM_GPU_T_H_)
#define _FDM_GPU_T_H_

#include"sde_builder_gpu.h"
#include"fdm_gpu.h"
#include<iostream>
#include<chrono>
#include<iomanip>


using namespace finite_difference_method_gpu;
using namespace sde_builder_gpu;

void fdm_gbm_gpu() {

	double r{ 0.05 };
	double sigma{ 0.01 };
	double s{ 100.0 };

	GeometricBrownianMotionGPU<> gbm{ r,sigma,s };
	std::cout << "Number of factors: " << GeometricBrownianMotionGPU<>::FactorCount << "\n";
	auto sde = gbm.sde();
	auto diffusion = gbm.diffusion();
	auto drift = gbm.drift();

	constexpr std::size_t factors = GeometricBrownianMotionGPU<>::FactorCount;
	FdmGPU<factors, double> gbm_fdm{1.0,2*360 };
	auto times = gbm_fdm.timeResolution();

	std::cout << "timing: \n";
	//braces here are to release the memory allocation before generating another collection of paths:
	{
		auto start = std::chrono::system_clock::now();
		auto paths_euler = gbm_fdm(sde, 70'000);
		auto end = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
		std::cout << "Euler took: " << end << " seconds\n";

	}
		auto start = std::chrono::system_clock::now();
		auto paths_milstein = gbm_fdm(sde, 70'000, FDMScheme::MilsteinScheme);
		auto end = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
		std::cout << "Milstein took: " << end << " seconds\n";
		std::cout << "\n";
	

	// To see few generated values:
	for (std::size_t t = 0; t < 30; ++t) {
		std::cout << t << " paths: \n";
		for (std::size_t v = 0; v < 10; ++v) {
			std::cout << paths_milstein[t][v] << ", ";
		}
		std::cout << "\n";
	}

}


void fdm_heston_gpu() {

	float r_d{ 0.05f };
	float r_f{ 0.01f };
	float mu = r_d - r_f;
	float sigma{ 0.01f };
	float theta{ 0.015f };
	float kappa{ 0.12f };
	float etha{ 0.012f };
	float stock_init{ 100.0f };
	float var_init{ 0.025f };


	HestonModelGPU<> hest{ mu,sigma,kappa,theta,etha,stock_init,var_init };
	std::cout << "Number of factors: " << HestonModel<>::FactorCount << "\n";

	constexpr std::size_t factors = HestonModel<>::FactorCount;
	FdmGPU<HestonModelGPU<>::FactorCount,double> heston_fdm{ 1.0f,0.8f,2 * 360 };
	auto times = heston_fdm.timeResolution();

	std::cout << "timing: \n";
	{
		auto start = std::chrono::system_clock::now();
		auto paths_euler = heston_fdm(hest.sdes(), 70'000);
		auto end = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
		std::cout << "Euler took: " << end << " seconds\n";
	}
	
		auto start = std::chrono::system_clock::now();
		auto paths_milstein = heston_fdm(hest.sdes(), 70'000, FDMScheme::MilsteinScheme);
		auto end = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
		std::cout << "Milstein took: " << end << " seconds\n";
		std::cout << "\n";
	

	

	// To see few generated values:
	for (std::size_t t = 0; t < 30; ++t) {
		std::cout << t << " paths: \n";
		for (std::size_t v = 0; v < 10; ++v) {
			std::cout <<std::setprecision(7)<< paths_milstein[t][v] << ", ";
		}
		std::cout << "\n";
	}
}





#endif ///_FDM_GPU_T_H_