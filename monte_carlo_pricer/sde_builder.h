#pragma once
#if !defined(_SDE_BUILDER_H_)
#define _SDE_BUILDER_H_

#include"sde.h"
#include<memory>

namespace sde_builder {

	using mc_types::SdeComponent;
	using mc_types::ISde;
	using sde::Sde;
	using mc_types::SdeModelType;


	template<std::size_t Factor,typename T,typename ...Ts>
	class SdeBuilder {
	};

	// Abstract specialized template class for one-factor models:
	template<typename T, typename ...Ts>
	class SdeBuilder<1,T,Ts...> {
	public:
		virtual ~SdeBuilder(){}
		virtual inline  std::string name()const = 0;
		virtual SdeComponent<T, Ts...> drift()const = 0;
		virtual SdeComponent<T, Ts...> diffusion()const = 0;
		virtual std::shared_ptr<Sde<T, Ts...>> model()const = 0;
		SdeModelType modelType()const { return SdeModelType::oneFactor; }
		enum {FactorCount = 1};
	};

	// Abstract specialized template class for two-factor models:
	template<typename T, typename ...Ts>
	class SdeBuilder<2, T, Ts...> {
	public:
		virtual ~SdeBuilder() {}
		virtual  inline std::string name()const = 0;
		virtual  SdeComponent<T, Ts...> drift1()const = 0;
		virtual  SdeComponent<T, Ts...> diffusion1()const = 0;
		virtual  SdeComponent<T, Ts...> drift2() const = 0;
		virtual  SdeComponent<T, Ts...> diffusion2()const = 0;
		virtual std::tuple<std::shared_ptr<Sde<T, Ts...>>,
							std::shared_ptr<Sde<T, Ts...>>> model() const = 0;
		SdeModelType modelType()const { return SdeModelType::twoFactor; }
		enum { FactorCount = 2 };
	};



	template<typename T = double,
				typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	class GeometricBrownianMotion: public SdeBuilder<1, T, T, T> {
	private:
		T mu_;
		T sigma_;
		T init_;

	public:
		GeometricBrownianMotion(T mu,T sigma,T initialCondition)
			:mu_{ mu }, sigma_{ sigma }, init_{initialCondition} {}
		GeometricBrownianMotion()
			:GeometricBrownianMotion{ 0.0,1.0,1.0 } {}

		GeometricBrownianMotion(GeometricBrownianMotion<T> const &copy)
			:mu_{copy.mu_},sigma_{copy.sigma_},init_{copy.init_}{}

		GeometricBrownianMotion& operator=(GeometricBrownianMotion<T> const &copy) {
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
		
		SdeComponent<T,T,T> drift()const override{
			return [this](T time,T underlyingPrice) {
				return mu_ * underlyingPrice;
			};
		}

		SdeComponent<T,T,T> diffusion()const override {
			return [this](T time,T underlyingPrice) {
				return sigma_ * underlyingPrice;
			};
		}
	
		std::shared_ptr<Sde<T,T,T>> model()const override {
			auto drift = this->drift();
			auto diff = this->diffusion();
			ISde<T,T,T> modelPair = std::make_tuple(drift, diff);
			return std::shared_ptr<Sde<T,T,T>>{ new Sde<T,T,T>{ modelPair,init_} };
		}

	};


	template<typename T=double,
			typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	class ArithmeticBrownianMotion:public SdeBuilder<1, T, T, T> {
	private:
		T mu_;
		T sigma_;
		T init_;

	public:
		ArithmeticBrownianMotion(T mu, T sigma, T initialCondition)
			:mu_{ mu }, sigma_{ sigma }, init_{ initialCondition } {}
		ArithmeticBrownianMotion()
			:ArithmeticBrownianMotion{ 0.0,1.0,1.0 } {}

		ArithmeticBrownianMotion(ArithmeticBrownianMotion<T> const &copy)
			:mu_{ copy.mu_ }, sigma_{ copy.sigma_ }, init_{ copy.init_ } {}

		ArithmeticBrownianMotion& operator=(ArithmeticBrownianMotion<T> const &copy) {
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

		SdeComponent<T,T,T> drift()const override {
			return [this](T time,T underlyingPrice) {
				return mu_;
			};
		}

		SdeComponent<T,T,T> diffusion()const override {
			return [this](T time,T underlyingPrice) {
				return sigma_;
			};
		}

		std::shared_ptr<Sde<T,T,T>> model() const override{
			auto drift = this->drift();
			auto diff = this->diffusion();
			ISde<T,T,T> modelPair = std::make_tuple(drift, diff);
			return std::shared_ptr<Sde<T,T,T>>{ new Sde<T,T,T>{ modelPair,init_ } };
		}
	};

	template<typename T=double,
		typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	class ConstantElasticityVariance:public SdeBuilder<1, T,T,T> {
	private:
		T beta_;
		T mu_;
		T sigma_;
		T init_;

	public:
		ConstantElasticityVariance(T mu, T sigma,T beta, T initialCondition)
			:mu_{ mu }, sigma_{ sigma }, beta_{beta}, init_ {
			initialCondition
		} {}
		ConstantElasticityVariance()
			:ConstantElasticityVariance{ 0.0,1.0,0.5,1.0 } {}

		ConstantElasticityVariance(ConstantElasticityVariance<T> const &copy)
			:mu_{ copy.mu_ }, sigma_{ copy.sigma_ }, beta_{copy.beta_}, 
			init_ {copy.init_} {}

		ConstantElasticityVariance& operator=(ConstantElasticityVariance<T> const &copy) {
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

		SdeComponent<T,T,T> drift()const override {
			return [this](T time,T underlyingPrice) {
				return mu_* underlyingPrice;
			};
		}

		SdeComponent<T,T,T> diffusion()const override {
			return [this](T time,T underlyingPrice) {
				return sigma_*std::pow(underlyingPrice,beta_);
			};
		}

		std::shared_ptr<Sde<T,T,T>> model()const override {
			auto drift = this->drift();
			auto diff = this->diffusion();
			ISde<T,T,T> modelPair = std::make_tuple(drift, diff);
			return std::shared_ptr<Sde<T,T,T>>{ new Sde<T,T,T>{ modelPair,init_ } };
		}
	};

	template<typename T = double,
			typename =typename std::enable_if<std::is_arithmetic<T>::value>::type>
	class HestonModel :public SdeBuilder<2, T,T,T,T> {
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
		HestonModel(T const &mu, T const &sigma, T const &kappa,T const &theta,T const &etha,
			T const &init1,T const &init2,T const &rho=0.0)
			:mu_{ mu }, sigma_{ sigma }, kappa_{ kappa }, theta_{theta},etha_{etha},
			init1_{ init1 }, init2_{ init2 }, rho_{rho} {}

		HestonModel()
			:HestonModel{ 0.5,0.05,0.05,0.05,0.05,1.0,0.01} {}

		HestonModel(HestonModel<T> const &copy)
			:mu_{ copy.mu_ }, sigma_{ copy.sigma_ }, kappa_{ copy.kappa_ },
			theta_{ copy.theta_ }, etha_{ copy.etha_ },
			init1_{ copy.init1_ }, init2_{ copy.init2_ },
			rho_{copy.rho_} {}

		HestonModel& operator=(HestonModel<T> const &copy) {
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

		SdeComponent<T,T,T,T> drift1()const override {
			return [this](T time, T underlyingPrice, T varianceProcess) {
				return mu_ * underlyingPrice;
			};
		}

		SdeComponent<T,T,T,T> diffusion1()const override {
			return [this](T time, T underlyingPrice,T varianceProcess) {
				return sigma_ * underlyingPrice  *std::sqrt(varianceProcess);
			};
		}

		SdeComponent<T,T,T,T> drift2() const override {
			return [this](T time,T underlyingPrice, T varianceProcess) {
				return kappa_*(theta_ - varianceProcess);
			};
		}

		SdeComponent<T,T,T,T> diffusion2()const override {
			return [this](T time, T underlyingPrice, T varianceProcess) {
				return etha_ * std::sqrt(varianceProcess);
			};
		}

		std::tuple<std::shared_ptr<Sde<T,T,T,T>>,std::shared_ptr<Sde<T,T,T,T>>> model()const override {
			auto drift1 = this->drift1();
			auto diff1 = this->diffusion1();
			auto drift2 = this->drift2();
			auto diff2 = this->diffusion2();
			ISde<T,T,T,T> modelPair1 = std::make_tuple(drift1, diff1);
			ISde<T,T,T,T> modelPair2 = std::make_tuple(drift2, diff2);
			return std::make_tuple(std::shared_ptr<Sde<T, T,T,T>>{ new Sde<T, T,T,T>{ modelPair1,init1_ } },
				std::shared_ptr<Sde<T, T,T,T>>{new Sde<T, T,T,T>{ modelPair2,init2_ }});
		}
	};



}


#endif ///_SDE_BUILDER_H_