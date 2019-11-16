#pragma once
#if !defined(_FDM_SCHEME_H_)
#define _FDM_SCHEME_H_

#include"mc_types.h"
#include"mc_utilities.h"
#include<random>
#include<cassert>

namespace finite_difference_method {


	using mc_types::TimePointsType;
	using mc_utilities::PartialCentralDifference;
	using mc_utilities::withRespectTo;


	template<typename T,typename ...Ts>
	using asyncKernel = std::function<PathValuesType<T>(Ts...)>;

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
		std::shared_ptr<Sde<T, Ts...>> model_;

	public:
		SchemeBuilder(std::shared_ptr<Sde<T, Ts...>> const &model,
					T const &delta,std::size_t numberSteps)
					:model_{model},delta_ {delta}, 
					numberSteps_{ numberSteps } {}
		SchemeBuilder(std::shared_ptr<Sde<T,Ts...>> const &model):
					model_{model}{}


		virtual PathValuesType<T> simulate(std::random_device::result_type seed) = 0;
		virtual PathValuesType<T> simulateWithTimePoints(std::random_device::result_type seed,
			TimePointsType<T> &timePoints) = 0;
	};

	// Scheme builder for two-factor models:
	template<typename T, typename ...Ts>
	class SchemeBuilder<2, T, Ts...> {
	protected:
		std::size_t numberSteps_;
		T delta_;
		T correlation_;
		std::tuple<std::shared_ptr<Sde<T, Ts...>>, std::shared_ptr<Sde<T, Ts...>>> model_;

	public:
		SchemeBuilder(std::tuple<std::shared_ptr<Sde<T, Ts...>>, std::shared_ptr<Sde<T, Ts...>>> const &model,
			T correlation,T const &delta, std::size_t numberSteps)
			:model_{ model }, delta_{ delta }, correlation_{ correlation },
			numberSteps_ {numberSteps} {}

		SchemeBuilder(std::tuple<std::shared_ptr<Sde<T, Ts...>>, std::shared_ptr<Sde<T, Ts...>>> const &model,
			T correlation)
			:model_{ model }, correlation_{ correlation } {}

		virtual PathValuesType<T> simulate(std::random_device::result_type seed) = 0;
		virtual PathValuesType<T> simulateWithTimePoints(std::random_device::result_type seed,
			TimePointsType<T> &timePoints) = 0;

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

		EulerScheme(std::shared_ptr<Sde<T,T,T>> const &model):
			SchemeBuilder<1, T, T, T>{model}{}


		PathValuesType<T> simulate(std::random_device::result_type seed) override {
			this->mt_.seed(seed );
			PathValuesType<T> path(this->numberSteps_);
			path[0] = this->model_->initCondition();
			auto spot = this->model_->initCondition();
			T spotNew{};
			for (std::size_t i = 1; i < path.size(); ++i) {
				spotNew = spot +
					this->model_->drift((i - 1)*(this->delta_), spot)*(this->delta_) +
					this->model_->diffusion((i - 1)*(this->delta_), spot) *
					std::sqrt((this->delta_)) * normal_(mt_);
				path[i] = spotNew;
				spot = spotNew;
			}
			return path;
		}

		PathValuesType<T> simulateWithTimePoints(std::random_device::result_type seed,TimePointsType<T> &timePoints)override {
			this->mt_.seed(seed);
			assert(!timePoints.empty());
			PathValuesType<T> path(timePoints.size());
			path[0] = this->model_->initCondition();
			auto spot = this->model_->initCondition();
			T spotNew{};
			for (std::size_t i = 1; i < path.size(); ++i) {
				spotNew = spot +
					this->model_->drift(timePoints[i - 1], spot)*(timePoints[i] - timePoints[i - 1]) +
					this->model_->diffusion(timePoints[i - 1], spot) *
					std::sqrt((timePoints[i] - timePoints[i - 1])) * normal_(mt_);
				path[i] = spotNew;
				spot = spotNew;
			}
			return path;
		}

	};

	template<typename T>
	class EulerScheme<2, T> :public SchemeBuilder<2, T, T, T, T> {
	private:
		std::normal_distribution<T> normal1_;
		std::normal_distribution<T> normal2_;
		std::mt19937 mt_;

	public:
		EulerScheme(std::tuple<std::shared_ptr<Sde<T, T, T, T>>, std::shared_ptr<Sde<T, T,T,T>>> const &model,
			T correlation, T const &delta, std::size_t numberSteps)
			:SchemeBuilder<2,T,T,T,T>{model,correlation,delta,numberSteps}{}

		EulerScheme(std::tuple<std::shared_ptr<Sde<T,T,T,T>>,std::shared_ptr<Sde<T,T,T,T>>> const &model,
			T correlation):
			SchemeBuilder<2,T,T,T,T>{model,correlation}{}

		PathValuesType<T> simulate(std::random_device::result_type seed)override {
			this->mt_.seed(seed);
			T z1{};
			T z2{};
			PathValuesType<T> path(this->numberSteps_);
			auto firstModel = std::get<0>(this->model_);
			auto secondModel = std::get<1>(this->model_);
			path[0] = firstModel->initCondition();
			auto firstSpot = firstModel->initCondition();
			T firstSpotNew{};
			auto secondSpot = secondModel->initCondition();
			T secondSpotNew{};

			for (std::size_t i = 1; i < path.size(); ++i) {
				z1 = normal1_(mt_);
				z2 = normal2_(mt_);
				firstSpotNew = firstSpot +
					firstModel->drift((i - 1)*(this->delta_), firstSpot,secondSpot)*(this->delta_) +
					firstModel->diffusion((i - 1)*(this->delta_), firstSpot, secondSpot) *
					std::sqrt((this->delta_)) * z1;
				secondSpotNew = secondSpot +
					secondModel->drift((i - 1)*(this->delta_), firstSpot, secondSpot)*(this->delta_) +
					secondModel->diffusion((i - 1)*(this->delta_), firstSpot, secondSpot) *
					std::sqrt((this->delta_)) * 
					(this->correlation_ * z1 + std::sqrt(1.0 - (this->correlation_ * this->correlation_)) * z2);
				path[i] = firstSpotNew;
				firstSpot = firstSpotNew;
				secondSpot = secondSpotNew;
			}
			return path;
		}

		PathValuesType<T> simulateWithTimePoints(std::random_device::result_type seed, TimePointsType<T> &timePoints)override {
			this->mt_.seed(seed);
			T z1{};
			T z2{};
			assert(!timePoints.empty());
			PathValuesType<T> path(timePoints.size());
			auto firstModel = std::get<0>(this->model_);
			auto secondModel = std::get<1>(this->model_);
			path[0] = firstModel->initCondition();
			auto firstSpot = firstModel->initCondition();
			T firstSpotNew{};
			auto secondSpot = secondModel->initCondition();
			T secondSpotNew{};

			for (std::size_t i = 1; i < path.size(); ++i) {
				z1 = normal1_(mt_);
				z2 = normal2_(mt_);
				firstSpotNew = firstSpot +
					firstModel->drift(timePoints[i - 1], firstSpot, secondSpot)*(timePoints[i] - timePoints[i - 1]) +
					firstModel->diffusion(timePoints[i - 1], firstSpot, secondSpot) *
					std::sqrt(timePoints[i] - timePoints[i - 1]) * z1;
				secondSpotNew = secondSpot +
					secondModel->drift(timePoints[i - 1], firstSpot, secondSpot)*(timePoints[i] - timePoints[i - 1]) +
					secondModel->diffusion(timePoints[i - 1], firstSpot, secondSpot) *
					std::sqrt(timePoints[i] - timePoints[i - 1]) *
					(this->correlation_ * z1 + std::sqrt(1.0 - (this->correlation_*this->correlation_)) * z2);
				path[i] = firstSpotNew;
				firstSpot = firstSpotNew;
				secondSpot = secondSpotNew;
			}
			return path;
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
		T step_ = 10e-6;

	public:
		MilsteinScheme(std::shared_ptr<Sde<T, T, T>> const &model,
			T const &delta, std::size_t numberSteps)
			:SchemeBuilder<1, T, T, T>{model,delta,numberSteps}{}

		MilsteinScheme(std::shared_ptr<Sde<T, T, T>> const &model) :
			SchemeBuilder<1, T, T, T>{ model } {};

		inline T diffusionPrime(T time,T price) {
			auto fun = pcd_(std::bind(&Sde<T, T, T>::diffusion, *(this->model_), std::placeholders::_1, std::placeholders::_2),
				withRespectTo::secondArg);
			return fun(time, price);
		}

		PathValuesType<T> simulate(std::random_device::result_type seed) override {
			this->mt_.seed(seed);
			PathValuesType<T> path(this->numberSteps_);
			path[0] = this->model_->initCondition();
			auto spot = this->model_->initCondition();
			T spotNew{};
			T z{};
			for (std::size_t i = 1; i < path.size(); ++i) {
				z = normal_(mt_);
				spotNew = spot +
					this->model_->drift((i - 1)*(this->delta_), spot)*(this->delta_) +
					this->model_->diffusion((i - 1)*(this->delta_), spot) *
					std::sqrt((this->delta_)) * z +
					0.5*this->model_->diffusion((i - 1) * (this->delta_), spot) *
					((this->model_->diffusion((i - 1) * (this->delta_), spot + 0.5*(this->step_)) - 
						this->model_->diffusion((i - 1) * (this->delta_), spot - 0.5*(this->step_)))/(this->step_))*
					((std::sqrt(this->delta_)*z)*(std::sqrt(this->delta_)*z) - (this->delta_));
				path[i] = spotNew;
				spot = spotNew;
			}
			return path;
		}

		PathValuesType<T> simulateWithTimePoints(std::random_device::result_type seed,TimePointsType<T> &timePoints)override {
			this->mt_.seed(seed);
			assert(!timePoints.empty());
			PathValuesType<T> path(timePoints.size());
			path[0] = this->model_->initCondition();
			auto spot = this->model_->initCondition();
			T spotNew{};
			T z{};
			for (std::size_t i = 0; i < path.size(); ++i) {
				z = normal_(mt_);
				spotNew = spot +
					this->model_->drift(timePoints[i - 1], spot)*(timePoints[i] - timePoints[i - 1]) +
					this->model_->diffusion(timePoints[i - 1], spot) *
					std::sqrt(timePoints[i] - timePoints[i - 1]) * z +
					0.5*this->model_->diffusion(timePoints[i - 1], spot)*
					((this->model_->diffusion(timePoints[i - 1], spot + 0.5*(this->step_)) -
						this->model_->diffusion(timePoints[i - 1], spot - 0.5*(this->step_))) / (this->step_))*
						((std::sqrt(timePoints[i] - timePoints[i - 1])*z)*(std::sqrt(timePoints[i] - timePoints[i - 1])*z) - 
					(timePoints[i] - timePoints[i - 1]));
				path[i] = spotNew;
				spot = spotNew;
			}
			return path;
		}
	};


	template<typename T>
	class MilsteinScheme<2, T> :public SchemeBuilder<2, T, T, T, T> {
	private:
		std::normal_distribution<T> normal1_;
		std::normal_distribution<T> normal2_;
		std::mt19937 mt_;
		T step_ = 10e-6;

	public:
		MilsteinScheme(std::tuple<std::shared_ptr<Sde<T, T,T,T>>, std::shared_ptr<Sde<T, T,T,T>>> const &model,
			T correlation, T const &delta, std::size_t numberSteps):
			SchemeBuilder<2,T,T,T,T>{model,correlation,delta,numberSteps}{
		}

		MilsteinScheme(std::tuple<std::shared_ptr<Sde<T,T,T,T>>,std::shared_ptr<Sde<T,T,T,T>>> const &model,
			T correlation):
			SchemeBuilder<2,T,T,T,T>{model,correlation}{}

		PathValuesType<T> simulate(std::random_device::result_type seed) override {
			mt_.seed(seed);
			T z1{};
			T z2{};

			PathValuesType<T> path(this->numberSteps_);
			auto firstModel = std::get<0>(this->model_);
			auto secondModel = std::get<1>(this->model_);
			path[0] = firstModel->initCondition();
			auto firstSpot = firstModel->initCondition();
			T firstSpotNew{};
			auto secondSpot = secondModel->initCondition();
			T secondSpotNew{};

			for (std::size_t i = 1; i < path.size(); ++i) {
				z1 = normal1_(mt_);
				z2 = normal2_(mt_);

				firstSpotNew = firstSpot +
					firstModel->drift((i - 1)*(this->delta_), firstSpot, secondSpot)*(this->delta_) +
					firstModel->diffusion((i - 1)*(this->delta_), firstSpot, secondSpot) *
					std::sqrt((this->delta_)) * z1 +
					0.5*firstModel->diffusion((i - 1)*(this->delta_), firstSpot, secondSpot) *
					((firstModel->diffusion((i - 1)*(this->delta_), firstSpot + 0.5*(this->step_), secondSpot) -
						firstModel->diffusion((i - 1)*(this->delta_), firstSpot - 0.5*(this->step_), secondSpot)) / (this->step_))*
						(this->delta_)*((z1)*(z1)-1.0) +
					0.5*(this->correlation_)*secondModel->diffusion((i - 1)*(this->delta_), firstSpot, secondSpot) *
					((firstModel->diffusion((i - 1)*(this->delta_), firstSpot, secondSpot + 0.5*(this->step_)) -
						firstModel->diffusion((i - 1)*(this->delta_), firstSpot, secondSpot - 0.5*(this->step_))) / (this->step_)) *
						(this->delta_)*((z1)*(z1)-1.0) +
					std::sqrt(1.0 - (this->correlation_)*(this->correlation_))*
					secondModel->diffusion((i - 1)*(this->delta_), firstSpot, secondSpot)*
					((firstModel->diffusion((i - 1)*(this->delta_), firstSpot, secondSpot + 0.5*(this->step_)) -
						firstModel->diffusion((i - 1)*(this->delta_), firstSpot, secondSpot - 0.5*(this->step_))) / (this->step_)) *
					(this->delta_)*z1*z2;

				secondSpotNew = secondSpot +
					secondModel->drift((i - 1)*(this->delta_), firstSpot, secondSpot)*(this->delta_) +
					secondModel->diffusion((i - 1)*(this->delta_), firstSpot, secondSpot) *
					std::sqrt((this->delta_)) *
					(this->correlation_ * z1 + std::sqrt(1.0 - (this->correlation_ * this->correlation_)) * z2) +
					0.5*(this->correlation_)*firstModel->diffusion((i - 1)*(this->delta_), firstSpot, secondSpot)*
					((secondModel->diffusion((i - 1)*(this->delta_), firstSpot + 0.5*(this->step_), secondSpot) -
						secondModel->diffusion((i - 1)*(this->delta_), firstSpot - 0.5*(this->step_), secondSpot)) / (this->step_))*
						(this->delta_) * ((z1)*(z1)-1.0) +
					0.5 * secondModel->diffusion((i - 1)*(this->delta_), firstSpot, secondSpot)*
					((secondModel->diffusion((i - 1)*(this->delta_), firstSpot, secondSpot + 0.5*(this->step_)) -
						secondModel->diffusion((i - 1)*(this->delta_), firstSpot, secondSpot - 0.5*(this->step_))) / (this->step_))*
						(this->delta_) * (((this->correlation_)*z1 + std::sqrt(1.0 - (this->correlation_)*(this->correlation_))*z2)*
					((this->correlation_)*z1 + std::sqrt(1.0 - (this->correlation_)*(this->correlation_))*z2) - 1.0) +
					std::sqrt(1.0 - (this->correlation_)*(this->correlation_))*
					firstModel->diffusion((i - 1)*(this->delta_), firstSpot, secondSpot)*
					((secondModel->diffusion((i - 1)*(this->delta_), firstSpot + 0.5*(this->step_), secondSpot) -
						secondModel->diffusion((i - 1)*(this->delta_), firstSpot - 0.5*(this->step_), secondSpot)) / (this->step_)) *
					(this->delta_)*z1*z2;

				path[i] = firstSpotNew;
				firstSpot = firstSpotNew;
				secondSpot = secondSpotNew;
			}
			return path;
		}

		PathValuesType<T> simulateWithTimePoints(std::random_device::result_type seed,TimePointsType<T> &timePoints)override{
			this->mt_.seed(seed);
			T z1{};
			T z2{};
			assert(!timePoints.empty());
			PathValuesType<T> path(timePoints.size());
			auto firstModel = std::get<0>(this->model_);
			auto secondModel = std::get<1>(this->model_);
			path[0] = firstModel->initCondition();
			auto firstSpot = firstModel->initCondition();
			T firstSpotNew{};
			auto secondSpot = secondModel->initCondition();
			T secondSpotNew{};

			for (std::size_t i = 1; i < path.size(); ++i) {
				z1 = normal1_(mt_);
				z2 = normal2_(mt_);

				firstSpotNew = firstSpot +
					firstModel->drift(timePoints[i - 1], firstSpot, secondSpot)*(timePoints[i] - timePoints[i - 1]) +
					firstModel->diffusion(timePoints[i - 1], firstSpot, secondSpot) *
					std::sqrt(timePoints[i] - timePoints[i - 1]) * z1 +
					0.5*firstModel->diffusion(timePoints[i - 1], firstSpot, secondSpot) *
					((firstModel->diffusion(timePoints[i - 1], firstSpot + 0.5*(this->step_), secondSpot) -
						firstModel->diffusion(timePoints[i - 1], firstSpot - 0.5*(this->step_), secondSpot)) / (this->step_))*
						(timePoints[i] - timePoints[i - 1])*((z1)*(z1)-1.0) +
					0.5*(this->correlation_)*secondModel->diffusion(timePoints[i - 1], firstSpot, secondSpot) *
					((firstModel->diffusion(timePoints[i - 1], firstSpot, secondSpot + 0.5*(this->step_)) -
						firstModel->diffusion(timePoints[i - 1], firstSpot, secondSpot - 0.5*(this->step_))) / (this->step_)) *
						(timePoints[i] -timePoints[i - 1])*((z1)*(z1)-1.0) +
					std::sqrt(1.0 - (this->correlation_)*(this->correlation_))*
					secondModel->diffusion(timePoints[i - 1], firstSpot, secondSpot)*
					((firstModel->diffusion(timePoints[i - 1], firstSpot, secondSpot + 0.5*(this->step_)) -
						firstModel->diffusion(timePoints[i - 1], firstSpot, secondSpot - 0.5*(this->step_))) / (this->step_)) *
						(timePoints[i] -timePoints[i - 1])*z1*z2;

				secondSpotNew = secondSpot +
					secondModel->drift(timePoints[i - 1], firstSpot, secondSpot)*(timePoints[i] - timePoints[i - 1]) +
					secondModel->diffusion(timePoints[i - 1], firstSpot, secondSpot) *
					std::sqrt((timePoints[i] - timePoints[i - 1])) *
					(this->correlation_ * z1 + std::sqrt(1.0 - (this->correlation_ * this->correlation_)) * z2) +
					0.5*(this->correlation_)*firstModel->diffusion(timePoints[i - 1], firstSpot, secondSpot)*
					((secondModel->diffusion(timePoints[i - 1], firstSpot + 0.5*(this->step_), secondSpot) -
						secondModel->diffusion(timePoints[i - 1], firstSpot - 0.5*(this->step_), secondSpot)) / (this->step_))*
						(timePoints[i] -timePoints[i - 1]) * ((z1)*(z1)-1.0) +
					0.5 * secondModel->diffusion(timePoints[i - 1], firstSpot, secondSpot)*
					((secondModel->diffusion(timePoints[i - 1], firstSpot, secondSpot + 0.5*(this->step_)) -
						secondModel->diffusion(timePoints[i - 1], firstSpot, secondSpot - 0.5*(this->step_))) / (this->step_))*
						(timePoints[i] -timePoints[i - 1]) * (((this->correlation_)*z1 + std::sqrt(1.0 - (this->correlation_)*(this->correlation_))*z2)*
					((this->correlation_)*z1 + std::sqrt(1.0 - (this->correlation_)*(this->correlation_))*z2) - 1.0) +
					std::sqrt(1.0 - (this->correlation_)*(this->correlation_))*
					firstModel->diffusion(timePoints[i - 1], firstSpot, secondSpot)*
					((secondModel->diffusion(timePoints[i - 1], firstSpot + 0.5*(this->step_), secondSpot) -
						secondModel->diffusion(timePoints[i - 1], firstSpot - 0.5*(this->step_), secondSpot)) / (this->step_)) *
						(timePoints[i] -timePoints[i - 1])*z1*z2;

				path[i] = firstSpotNew;
				firstSpot = firstSpotNew;
				secondSpot = secondSpotNew;
			}
			return path;
		}



	};

}



#endif ///_FDM_SCHEME_H_