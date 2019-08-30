#pragma once
#if !defined(_SDE_BUILDER_H_)
#define _SDE_BUILDER_H_

#include"sde.h"
#include<memory>

namespace sde_builder {

	using mc_types::SdeComponent;
	using mc_types::ISde;
	using sde::Sde;


	template<typename T = double,
		typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	class SdeBuilder {
	public:
		virtual SdeComponent<T> drift() = 0;
		virtual SdeComponent<T> diffusion() = 0;
		virtual std::shared_ptr<Sde<T>> model() = 0;
	};



	template<typename T = double,
				typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	class GeometricBrownianMotion: public SdeBuilder<T> {
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
		
		SdeComponent<T> drift() {
			return [this](T const &underlyingPrice,T const &time) {
				return mu_ * underlyingPrice;
			};
		}

		SdeComponent<T> diffusion() {
			return [this](T const &underlyingPrice, T const &time) {
				return sigma_ * underlyingPrice;
			};
		}
	
		std::shared_ptr<Sde<T>> model() {
			auto drift = this->drift();
			auto diff = this->diffusion();
			ISde<T> modelPair = std::make_tuple(drift, diff);
			return std::shared_ptr<Sde<T>>{ new Sde<T>{ modelPair,init_} };
		}

	};


	template<typename T=double,
			typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	class ArithmeticBrownianMotion:public SdeBuilder<T> {
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

		SdeComponent<T> drift() {
			return [this](T const &underlyingPrice, T const &time) {
				return mu_;
			};
		}

		SdeComponent<T> diffusion() {
			return [this](T const &underlyingPrice, T const &time) {
				return sigma_;
			};
		}

		std::shared_ptr<Sde<T>> model() {
			auto drift = this->drift();
			auto diff = this->diffusion();
			ISde<T> modelPair = std::make_tuple(drift, diff);
			return std::shared_ptr<Sde<T>>{ new Sde<T>{ modelPair,init_ } };
		}
	};

	template<typename T=double,
		typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	class ConstantElasticityVariance:public SdeBuilder<T> {
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

		SdeComponent<T> drift() {
			return [this](T const &underlyingPrice, T const &time) {
				return mu_* underlyingPrice;
			};
		}

		SdeComponent<T> diffusion() {
			return [this](T const &underlyingPrice, T const &time) {
				return sigma_*std::pow(underlyingPrice,beta_);
			};
		}

		std::shared_ptr<Sde<T>> model() {
			auto drift = this->drift();
			auto diff = this->diffusion();
			ISde<T> modelPair = std::make_tuple(drift, diff);
			return std::shared_ptr<Sde<T>>{ new Sde<T>{ modelPair,init_ } };
		}
	};




}


#endif ///_SDE_BUILDER_H_