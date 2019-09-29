#pragma once
#if !defined(_FDM_SCHEME_H_)
#define _FDM_SCHEME_H_

#include"mc_types.h"
#include"mc_utilities.h"
#include<random>

namespace finite_difference_method {


	using mc_types::TimePointsType;
	using mc_utilities::PartialCentralDifference;

	template<typename T>
	using kernel = std::function<void(PathValuesType<T> &)>;

	template<std::size_t FactorCount,typename T,typename ...Ts>
	class SchemeBuilder {
	public:
		virtual ~SchemeBuilder(){}
	};

	// Scheme builder for one-factor models
	template<typename T,typename ...Ts>
	class SchemeBuilder<1,T,Ts...> {
	protected:
		std::size_t numberSteps_;
		T delta_;
		TimePointsType<T> timePoints_;
		std::shared_ptr<Sde<T, Ts...>> model_;

	public:
		SchemeBuilder(std::shared_ptr<Sde<T, Ts...>> const &model,
					T const &delta,std::size_t numberSteps)
					:model_{model},delta_ {delta}, 
					numberSteps_{ numberSteps } {}

		SchemeBuilder(std::shared_ptr<Sde<T, Ts...>> const &model, 
					TimePointsType<T> const &timePoints, std::size_t numberSteps)
					:model_{model}, timePoints_{ timePoints },
					numberSteps_{ numberSteps } {}

		virtual void operator()(PathValuesType<T> &values) = 0;
	};

	// Scheme builder for two-factor models:
	template<typename T, typename ...Ts>
	class SchemeBuilder<2, T, Ts...> {
	protected:
		std::size_t numberSteps_;
		T delta_;
		T correlation_;
		TimePointsType<T> timePoints_;
		std::tuple<std::shared_ptr<Sde<T, Ts...>>, std::shared_ptr<Sde<T, Ts...>>> model_;

	public:
		SchemeBuilder(std::tuple<std::shared_ptr<Sde<T, Ts...>>, std::shared_ptr<Sde<T, Ts...>>> const &model,
			T const &correlation,T const &delta, std::size_t numberSteps)
			:model_{ model }, delta_{ delta }, correlation_{ correlation },
			numberSteps_ {numberSteps} {}

		SchemeBuilder(std::tuple<std::shared_ptr<Sde<T, Ts...>>, std::shared_ptr<Sde<T, Ts...>>> const &model,
			T const &correlation,TimePointsType<T> const &timePoints, std::size_t numberSteps)
			:model_{ model }, timePoints_{ timePoints }, correlation_{correlation},
			numberSteps_ {numberSteps} {}
	};


	template<std::size_t FactorCount,
			typename T,
			typename = typename std::enable_if<(std::is_arithmetic<T>::value) && 
												(FactorCount > 0)>::type>
	class EulerScheme{};


	template<typename T>
	class EulerScheme<1, T> :public SchemeBuilder<1, T, T, T> {
	private:
		std::normal_distribution<T> normal_;
		std::mt19937 mt_;

	public:
		EulerScheme(std::shared_ptr<Sde<T,T,T>> const &model,
			T const &delta, std::size_t numberSteps)
			:SchemeBuilder<1,T,T,T>{model,delta,numberSteps}{
		}

		EulerScheme(std::shared_ptr<Sde<T, T, T>> const &model,
			TimePointsType<T> const &timePoints, std::size_t numberSteps)
			:SchemeBuilder<1, T, T, T>{ model,timePoints,numberSteps } {
		}

		void operator()(PathValuesType<T> &values) {
			values[0] = this->model_->initCondition();
			for (std::size_t i = 1; i < values.size(); ++i) {
				values[i] = values[i - 1] +
					this->model_->drift((i-1)*(this->delta_), values[i - 1])*(this->delta_) +
					this->model_->diffusion((i - 1)*(this->delta_), values[i - 1]) * 
					std::sqrt((this->delta_)) * normal_(mt_);
			}
		}
	};


	template<std::size_t FactorCount,
			typename T,
			typename  = typename std::enable_if<(std::is_arithmetic<T>::value) && 
												(FactorCount > 0)>::type>
	class MilsteinScheme{};


	template<typename T>
	class MilsteinScheme<1, T> :public SchemeBuilder<1, T, T, T> {
	private:
		std::normal_distribution<T> normal_;
		std::mt19937 mt_;

	public:
		MilsteinScheme(std::shared_ptr<Sde<T, T, T>> const &model,
			T const &delta, std::size_t numberSteps)
			:SchemeBuilder<1, T, T, T>{model,delta,numberSteps}{}

		MilsteinScheme(std::shared_ptr<Sde<T, T, T>> const &model,
			TimePointsType<T> const &timePoints, std::size_t numberSteps)
			:SchemeBuilder<1, T, T, T>{model,timePoints,numberSteps}{}

		void operator()(PathValuesType<T> &values) {
			values[0] = this->model_->initCondition();
			T z{};
			for (std::size_t i = 1; i < values.size(); ++i) {
				z = normal_(mt_);
				auto diffusion_proc = std::bind(&Sde<T,T,T>::diffusion,*(this->model_),
												(i - 1) * (this->delta_),std::placeholders::_1);
				auto diffusion_prime = PartialCentralDifference<1, T>()(std::forward<decltype(diffusion_proc)>(diffusion_proc));
				values[i] = values[i - 1] +
					this->model_->drift((i - 1)*(this->delta_), values[i - 1])*(this->delta_) +
					this->model_->diffusion((i - 1)*(this->delta_), values[i - 1])*
					std::sqrt((this->delta_)) * z +
					0.5*this->model_->diffusion((i - 1) * (this->delta_), values[i - 1]) * diffusion_prime(values[i - 1])*
					((std::sqrt(this->delta_)*z)*(std::sqrt(this->delta_)*z) - (this->delta_));
			}
		}


	};


}



#endif ///_FDM_SCHEME_H_