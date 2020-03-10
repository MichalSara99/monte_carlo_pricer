#pragma once
#if !defined(_SDE_BUILDER_GPU_H_)
#define _SDE_BUILDER_GPU_H_

#include"mc_types.h"
#include<tuple>

namespace sde_builder_gpu {

	using mc_types::SdeModelType;

	template<std::size_t Factor,typename T>
	class SdeBuilderGPU{};

	// Abstract specialized temnplate class for on-factoir models
	template<typename T>
	class SdeBuilderGPU<1, T> {
	public:
		virtual ~SdeBuilderGPU(){}
		SdeModelType modelType()const { return SdeModelType::oneFactor; }
		enum { FactorCount = 1 };
		inline virtual std::string name() const = 0;

	};

	// Abstract specialized temnplate class for on-factoir models
	template<typename T>
	class SdeBuilderGPU<2, T> {
	public:
		virtual ~SdeBuilderGPU() {}
		SdeModelType modelType()const { return SdeModelType::twoFactor; }
		enum { FactorCount = 2; };
	};


	template<typename T = double,
		typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
		class GeometricBrownianMotionGPU : public SdeBuilderGPU<1, T> {
		private:
			T mu_;
			T sigma_;
			T init_;

		public:
			GeometricBrownianMotionGPU(T mu, T sigma, T initialCondition)
				:mu_{ mu }, sigma_{ sigma }, init_{ initialCondition } {}
			GeometricBrownianMotionGPU()
				:GeometricBrownianMotion{ 0.0,1.0,1.0 } {}

			GeometricBrownianMotionGPU(GeometricBrownianMotionGPU<T> const &copy)
				:mu_{ copy.mu_ }, sigma_{ copy.sigma_ }, init_{ copy.init_ } {}

			GeometricBrownianMotionGPU& operator=(GeometricBrownianMotionGPU<T> const &copy) {
				if (this != &copy) {
					mu_ = copy.mu_;
					sigma_ = copy.sigma_;
					init_ = copy.init_;
				}
				return *this;
			}

			inline T const &mu()const { return mu_; }
			inline T const &sigma()const { return sigma_; }
			inline T const &init()const { return init_; }

			inline std::string name() const override { return std::string{ "Geometric Brownian Motion" }; }

			auto drift()const {
				T mu = mu_;
				return [=](T time, T underlyingPrice) restrict(amp)->T {
					return mu * underlyingPrice;
				};
			}

			auto diffusion()const {
				T sigma = sigma_;
				return [=](T time, T underlyingPrice) restrict(amp)->T {
					return sigma * underlyingPrice;
				};
			}

			auto model()const {
				T init = init_;
				return std::make_tuple(drift(), diffusion(), init);
			}
	};



}




#endif ///_SDE_BUILDER_GPU_H_