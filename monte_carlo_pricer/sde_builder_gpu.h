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
		enum { FactorCount = 2 };
		inline virtual std::string name() const = 0;
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

			auto sde()const {
				T init = init_;
				return std::make_tuple(drift(), diffusion(), init);
			}
	};

	template<typename T = double,
		typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
		class ArithmeticBrownianMotionGPU :public SdeBuilder<1, T> {
		private:
			T mu_;
			T sigma_;
			T init_;

		public:
			ArithmeticBrownianMotionGPU(T mu, T sigma, T initialCondition)
				:mu_{ mu }, sigma_{ sigma }, init_{ initialCondition } {}
			ArithmeticBrownianMotionGPU()
				:ArithmeticBrownianMotionGPU{ 0.0,1.0,1.0 } {}

			ArithmeticBrownianMotionGPU(ArithmeticBrownianMotionGPU<T> const &copy)
				:mu_{ copy.mu_ }, sigma_{ copy.sigma_ }, init_{ copy.init_ } {}

			ArithmeticBrownianMotionGPU& operator=(ArithmeticBrownianMotionGPU<T> const &copy) {
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

			inline std::string name() const override { return std::string{ "Arithmetic Brownian Motion" }; }

			auto drift()const  {
				auto mu = this->mu_;
				return [=](T time, T underlyingPrice)restrict(amp)->T {
					return mu;
				};
			}

			auto diffusion()const  {
				auto sigma = this->sigma_;
				return [=](T time, T underlyingPrice)restrict(amp)->T {
					return sigma;
				};
			}

			auto model() const {
				T init = this->init_;
				return std::make_tuple(drift(), diffusion(), init);
			}
	};

	template<typename T = double,
		typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
		class ConstantElasticityVarianceGPU :public SdeBuilder<1, T> {
		private:
			T beta_;
			T mu_;
			T sigma_;
			T init_;

		public:
			ConstantElasticityVarianceGPU(T mu, T sigma, T beta, T initialCondition)
				:mu_{ mu }, sigma_{ sigma }, beta_{ beta }, init_{
				initialCondition
			} {}
			ConstantElasticityVarianceGPU()
				:ConstantElasticityVarianceGPU{ 0.0,1.0,0.5,1.0 } {}

			ConstantElasticityVarianceGPU(ConstantElasticityVarianceGPU<T> const &copy)
				:mu_{ copy.mu_ }, sigma_{ copy.sigma_ }, beta_{ copy.beta_ },
				init_{ copy.init_ } {}

			ConstantElasticityVarianceGPU& operator=(ConstantElasticityVarianceGPU<T> const &copy) {
				if (this != &copy) {
					mu_ = copy.mu_;
					sigma_ = copy.sigma_;
					beta_ = copy.beta_;
					init_ = copy.init_;
				}
				return *this;
			}

			inline T const &mu()const { return mu_; }
			inline T const &sigma()const { return sigma_; }
			inline T const &init()const { return init_; }
			inline T const &beta()const { return beta_; }

			inline std::string name() const override { return std::string{ "Constant Elasticity Variance" }; }

			auto drift()const  {
				auto mu = this->mu_;
				return [=](T time, T underlyingPrice)restrict(amp)->T {
					return mu * underlyingPrice;
				};
			}

			auto diffusion()const {
				auto beta = this->beta_;
				auto sigma = this->sigma_;
				return [=](T time, T underlyingPrice)restrict(amp)->T {
					return sigma * concurrency::fast_math::pow(underlyingPrice, beta);
				};
			}

			auto model() const {
				T init = this->init_;
				return std::make_tuple(drift(), diffusion(), init);
			}
	};


	template<typename T = double,
		typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
		class HestonModelGPU :public SdeBuilderGPU<2, T> {
		private:
			// for first underlying factor:
			T init1_;
			T mu_;
			T sigma_;
			// for second variance factor:
			T init2_;
			T kappa_;
			T theta_;
			T etha_;
			// correlation between factors:
			T rho_;

		public:
			HestonModelGPU(T mu, T sigma, T kappa, T theta, T etha,
				T init1, T init2, T rho = 0.0)
				:mu_{ mu }, sigma_{ sigma }, kappa_{ kappa }, theta_{ theta }, etha_{ etha },
				init1_{ init1 }, init2_{ init2 }, rho_{ rho } {}

			HestonModelGPU()
				:HestonModel{ 0.5,0.05,0.05,0.05,0.05,1.0,0.01 } {}

			HestonModelGPU(HestonModelGPU<T> const &copy)
				:mu_{ copy.mu_ }, sigma_{ copy.sigma_ }, kappa_{ copy.kappa_ },
				theta_{ copy.theta_ }, etha_{ copy.etha_ },
				init1_{ copy.init1_ }, init2_{ copy.init2_ },
				rho_{ copy.rho_ } {}

			HestonModelGPU& operator=(HestonModelGPU<T> const &copy) {
				if (this != &copy) {
					mu_ = copy.mu_;
					sigma_ = copy.sigma_;
					kappa_ = copy.kappa_;
					theta_ = copy.theta_;
					etha_ = copy.etha_;
					init1_ = copy.init1_;
					init2_ = copy.init2_;
					rho_ = copy.rho_;
				}
				return *this;
			}

			inline T const &mu()const { return mu_; }
			inline T const &sigma()const { return sigma_; }
			inline T const &kappa()const { return kappa_; }
			inline T const &theta()const { return theta_; }
			inline T const &etha()const { return etha_; }
			inline T const &init1()const { return init1_; }
			inline T const &init2()const { return init2_; }
			inline T const &rho()const { return rho_; }

			inline std::string name() const override { return std::string{ "Heston Model" }; }

			auto drift1()const {
				auto mu = this->mu_;
				return [=](T time, T underlyingPrice, T varianceProcess)restrict(amp)->T {
					return mu * underlyingPrice;
				};
			}

			auto diffusion1()const  {
				auto sigma = this->sigma_;
				return [=](T time, T underlyingPrice, T varianceProcess)restrict(amp)->T {
					return sigma * underlyingPrice  *concurrency::fast_math::sqrt(varianceProcess);
				};
			}

			auto drift2() const {
				auto kappa = this->kappa_;
				auto theta = this->theta_;
				return [=](T time, T underlyingPrice, T varianceProcess)restrict(amp)->T {
					return kappa * (theta - varianceProcess);
				};
			}

			auto diffusion2()const  {
				auto etha = this->etha_;
				return [=](T time, T underlyingPrice, T varianceProcess)restrict(amp)->T {
					return etha * concurrency::fast_math::sqrt(varianceProcess);
				};
			}

			auto sde1()const {
				T init1 = this->init1_;
				return std::make_tuple(drift1(), diffusion1(), init1);
			}

			auto sde2()const {
				T init2 = this->init2_;
				return std::make_tuple(drift2(), diffusion2(), init2);
			}

			auto sdes()const {
				return std::make_tuple(sde1(), sde2());
			}
	};

}




#endif ///_SDE_BUILDER_GPU_H_