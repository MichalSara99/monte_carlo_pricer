#pragma once
#if !defined(_SDE_BUILDER_H_)
#define _SDE_BUILDER_H_

#include"sde.h"
#include<memory>

namespace sde_builder {

	using mc_types::SdeComponent2Args;
	using mc_types::SdeComponent3Args;
	using mc_types::ISde1;
	using mc_types::ISde2;
	using sde::Sde;
	using mc_types::SdeModelType;


	template<std::size_t Factor,
		typename T = double,
		typename = typename std::enable_if<(std::is_arithmetic<T>::value) &&
											(Factor > 0)>::type>
	class SdeBuilder {
		enum { FactorCount = SdeBuilder<Factor,T>::FactorCount};
	};

	// Abstract specialized template class for one-factor models:
	template<typename T>
	class SdeBuilder<1,T> {
	public:
		virtual ~SdeBuilder(){}
		virtual SdeComponent2Args<T> drift() = 0;
		virtual SdeComponent2Args<T> diffusion() = 0;
		virtual std::shared_ptr<Sde<1,T>> model() = 0;
		SdeModelType modelType()const { return SdeModelType::oneFactor; }
		enum {FactorCount = 1};
	};

	// Abstract specialized template class for two-factor models:
	template<typename T>
	class SdeBuilder<2, T> {
	public:
		virtual ~SdeBuilder() {}
		virtual std::tuple<SdeComponent3Args<T>, SdeComponent3Args<T>> drifts() = 0;
		virtual std::tuple<SdeComponent3Args<T>, SdeComponent3Args<T>> diffusions() = 0;
		virtual std::tuple<std::shared_ptr<Sde<2,T>>, std::shared_ptr<Sde<2,T>>> models() = 0;
		SdeModelType modelType()const { return SdeModelType::twoFactor; }
		enum { FactorCount = 2 };
	};




	template<typename T = double,
				typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	class GeometricBrownianMotion: public SdeBuilder<1,T> {
	private:
		T mu_;
		T sigma_;
		T init_;

	public:
		GeometricBrownianMotion(T const &mu,T const &sigma,T const &initialCondition)
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
		
		SdeComponent2Args<T> drift() {
			return [this](T const &underlyingPrice,T const &time) {
				return mu_ * underlyingPrice;
			};
		}

		SdeComponent2Args<T> diffusion() {
			return [this](T const &underlyingPrice, T const &time) {
				return sigma_ * underlyingPrice;
			};
		}
	
		std::shared_ptr<Sde<1,T>> model() {
			auto drift = this->drift();
			auto diff = this->diffusion();
			ISde1<T> modelPair = std::make_tuple(drift, diff);
			return std::shared_ptr<Sde<1,T>>{ new Sde<1,T>{ modelPair,init_} };
		}

	};


	template<typename T=double,
			typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	class ArithmeticBrownianMotion:public SdeBuilder<1, T> {
	private:
		T mu_;
		T sigma_;
		T init_;

	public:
		ArithmeticBrownianMotion(T const &mu, T const &sigma, T const &initialCondition)
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

		SdeComponent2Args<T> drift() {
			return [this](T const &underlyingPrice, T const &time) {
				return mu_;
			};
		}

		SdeComponent2Args<T> diffusion() {
			return [this](T const &underlyingPrice, T const &time) {
				return sigma_;
			};
		}

		std::shared_ptr<Sde<1,T>> model() {
			auto drift = this->drift();
			auto diff = this->diffusion();
			ISde1<T> modelPair = std::make_tuple(drift, diff);
			return std::shared_ptr<Sde<1,T>>{ new Sde<1,T>{ modelPair,init_ } };
		}
	};

	template<typename T=double,
		typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	class ConstantElasticityVariance:public SdeBuilder<1, T> {
	private:
		T beta_;
		T mu_;
		T sigma_;
		T init_;

	public:
		ConstantElasticityVariance(T const &mu, T const &sigma,T const &beta, T const &initialCondition)
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

		SdeComponent2Args<T> drift() {
			return [this](T const &underlyingPrice, T const &time) {
				return mu_* underlyingPrice;
			};
		}

		SdeComponent2Args<T> diffusion() {
			return [this](T const &underlyingPrice, T const &time) {
				return sigma_*std::pow(underlyingPrice,beta_);
			};
		}

		std::shared_ptr<Sde<1,T>> model() {
			auto drift = this->drift();
			auto diff = this->diffusion();
			ISde1<T> modelPair = std::make_tuple(drift, diff);
			return std::shared_ptr<Sde<1,T>>{ new Sde<1,T>{ modelPair,init_ } };
		}
	};

	template<typename T = double,
			typename =typename std::enable_if<std::is_arithmetic<T>::value>::type>
	class HestonModel :public SdeBuilder<2, T> {
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

	public:
		HestonModel(T const &mu, T const &sigma, T const &kappa,T const &theta,T const &etha,
			T const &init1,T const &init2)
			:mu_{ mu }, sigma_{ sigma }, kappa_{ kappa }, theta_{theta},etha_{etha},
			init1_ {init1}, init2_{init2} {}

		HestonModel()
			:HestonModel{ 0.5,0.05,0.05,0.05,0.05,1.0,0.01} {}

		HestonModel(HestonModel<T> const &copy)
			:mu_{ copy.mu_ }, sigma_{ copy.sigma_ }, kappa_{ copy.kappa_ },
			theta_{ copy.theta_ }, etha_{ copy.etha_ },
			init1_{ copy.init1_ }, init2_{ copy.init2_ } {}

		HestonModel& operator=(HestonModel<T> const &copy) {
			if (this != &copy) {
				mu_ = copy.mu_;
				sigma_ = copy.sigma_;
				kappa_ = copy.kappa_;
				theta_ = copy.theta_;
				etha_ = copy.etha_;
				init1_ = copy.init1_;
				init2_ = copy.init2_;
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

		SdeComponent3Args<T> drift1() {
			return [this](T const &underlyingPrice, T const &varianceProcess, T const &time) {
				return mu_ * underlyingPrice;
			};
		}

		SdeComponent3Args<T> diffusion1() {
			return [this](T const &underlyingPrice,T const &varianceProcess, T const &time) {
				return sigma_ * underlyingPrice  *std::sqrt(varianceProcess);
			};
		}

		SdeComponent3Args<T> drift2() {
			return [this](T const &underlyingPrice, T const &varianceProcess, T const &time) {
				return kappa_*(theta_ - varianceProcess);
			};
		}

		SdeComponent3Args<T> diffusion2() {
			return [this](T const &underlyingPrice, T const &varianceProcess, T const &time) {
				return etha_ * std::sqrt(varianceProcess);
			};
		}

		std::tuple<std::shared_ptr<Sde<2,T>>,std::shared_ptr<Sde<2,T>>> model() {
			auto drift1 = this->drift1();
			auto diff1 = this->diffusion1();
			auto drift2 = this->drift2();
			auto diff2 = this->diffusion2();
			ISde2<T> modelPair1 = std::make_tuple(drift1, diff1);
			ISde2<T> modelPair2 = std::make_tuple(drift2, diff2);
			return std::make_tuple(std::shared_ptr<Sde<2, T>>{ new Sde<2, T>{ modelPair1,init1_ } },
				std::shared_ptr<Sde<2, T>>{new Sde<2, T>{ modelPair2,init2_ }});
		}
	};



}


#endif ///_SDE_BUILDER_H_