#pragma once
#if !defined(_FDM_GPU_T_H_)
#define _FDM_GPU_T_H_

#include"sde_builder_gpu.h"
#include"fdm_gpu.h"
#include<iostream>
#include<chrono>


using namespace finite_difference_method_gpu;
using namespace sde_builder_gpu;

void fdm_gbm_gpu() {

	double r{ 0.05 };
	double sigma{ 0.01 };
	double s{ 100.0 };

	GeometricBrownianMotionGPU<> gbm{ r,sigma,s };
	std::cout << "Number of factors: " << GeometricBrownianMotionGPU<>::FactorCount << "\n";
	auto sde = gbm.model();
	auto diffusion = gbm.diffusion();
	auto drift = gbm.drift();

	constexpr std::size_t factors = GeometricBrownianMotionGPU<>::FactorCount;
	FdmGPU<factors, double> gbm_fdm{1.0,3*360 };
	auto times = gbm_fdm.timeResolution();

	std::cout << "timing: \n";
	auto start = std::chrono::system_clock::now();
	auto paths_euler = gbm_fdm(sde,100'000);
	auto end = std::chrono::duration<double>(std::chrono::system_clock::now() - start).count();
	std::cout << "Euler took: " << end << " seconds\n";
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






#endif ///_FDM_GPU_T_H_