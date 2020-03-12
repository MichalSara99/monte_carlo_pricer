#pragma once
#if !defined(_FDM_GPU_H_)
#define _FDM_GPU_H_

#include<random>
#include"mc_types.h"
#include"fdm_scheme_gpu.h"
#include<algorithm>


namespace finite_difference_method_gpu {


	using mc_types::TimePointsType;
	using mc_types::PathValuesType;
	using mc_types::FDMScheme;

	template<std::size_t FactorCount,typename T>
	class FdmBuilderGPU {
	public:
		virtual ~FdmBuilderGPU(){}
	};

	// Finite difference method builder for one-factor models:
	template<typename T>
	class FdmBuilderGPU<1, T> {
	protected:
		bool timePointsOn_;
		T terminationTime_;
		TimePointsType<T> timePoints_;
		std::size_t numberSteps_;

	public:
		explicit FdmBuilderGPU(T const &terminationTime, std::size_t numberSteps = 360)
			:terminationTime_{ terminationTime },
			numberSteps_{ numberSteps },
			timePointsOn_{false}{}

		explicit FdmBuilderGPU(TimePointsType<T> const &timePoints)
			:timePoints_{ timePoints }, timePointsOn_{true} {}

		inline TimePointsType<T> timeResolution()const {
			if (timePointsOn_ == false) {
				TimePointsType<T> points(numberSteps_ + 1);
				auto delta = (terminationTime_ / static_cast<T>(numberSteps_));
				std::size_t i = 0;
				std::generate(points.begin(), points.end(), [&i,&delta]() {return delta * i++; });
				return points;
			}
			return timePoints_;
		}
	};

	// Finite difference method builder for one-factor models:
	template<typename T>
	class FdmBuilderGPU<2, T> {
	protected:
		bool timePointsOn_;
		T terminationTime_;
		T correlation_;
		TimePointsType<T> timePoints_;
		std::size_t numberSteps_;

	public:
		explicit FdmBuilderGPU(T const &terminationTime, T correlation = 0.0, std::size_t numberSteps = 360)
			:terminationTime_{ terminationTime },
			numberSteps_{ numberSteps },
			correlation_{correlation},
			timePointsOn_{ false } {}

		explicit FdmBuilderGPU(TimePointsType<T> const &timePoints, T correlation = 0.0)
			:timePoints_{ timePoints },
			correlation_{ correlation },
			timePointsOn_ {true} {}

		inline TimePointsType<T> timeResolution()const {
			if (timePointsOn_ == false) {
				TimePointsType<T> points(numberSteps_ + 1);
				auto delta = (terminationTime_ / static_cast<T>(numberSteps_));
				std::size_t i = 0;
				std::generate(points.begin(), points.end(), [&i, &delta]() {return delta * i++; });
				return points;
			}
			return timePoints_;
		}
	};


	template<std::size_t FactorCount,
			typename T,
			typename =typename std::enable_if<std::is_floating_point<T>::value && (FactorCount>0)>::type>
	class FdmGPU{};

	template<typename T>
	class FdmGPU<1, T> :public FdmBuilderGPU<1, T> {
	private:
		std::random_device rd_;

	public:
		FdmGPU(T const &terminationTime,std::size_t numberSteps = 360)
			:FdmBuilderGPU<1,T>{terminationTime,numberSteps}{}

		FdmGPU(TimePointsType<T> const &timePoints)
			:FdmBuilderGPU<1,T>{timePoints}{}

		template<typename SDE>
		PathValuesType<PathValuesType<T>> operator()(SDE const &sde, std::size_t iterations,
			FDMScheme scheme = FDMScheme::EulerScheme) {

			T delta = (this->terminationTime_ / static_cast<T>(this->numberSteps_));
			std::size_t seed = 56489;
			switch (scheme) {
			case FDMScheme::EulerScheme:
			{
				if (this->timePointsOn_ == false) {
					EulerSchemeGPU<1, T> euler(delta, this->numberSteps_);
					return euler.simulate(sde, iterations, seed);
				}
				else {
					EulerSchemeGPU<1, T> euler;
					const TimePointsType<T> timePoints = this->timePoints_;
					return euler.simulateWithTimePoints(sde, iterations,timePoints, seed);
				}
			}
			break;
			case FDMScheme::MilsteinScheme:
			{
				if (this->timePointsOn_ == false) {
					MilsteinSchemeGPU<1, T> milstein(delta, this->numberSteps_);
					return milstein.simulate(sde, iterations, seed);
				}
				else {
					MilsteinSchemeGPU<1, T> milstein;
					const TimePointsType<T> timePoints = this->timePoints_;
					return milstein.simulateWithTimePoints(sde, iterations, timePoints, seed);
				}
			}
			break;
			default:
			{
				if (this->timePointsOn_ == false) {
					EulerSchemeGPU<1, T> euler(delta, this->numberSteps_);
					return euler.simulate(sde, iterations, seed);
				}
				else {
					EulerSchemeGPU<1, T> euler;
					const TimePointsType<T> timePoints = this->timePoints_;
					return euler.simulateWithTimePoints(sde, iterations, timePoints, seed);
				}
			}
			}
		}
	};

	template<typename T>
	class FdmGPU<2, T> :public FdmBuilderGPU<2, T> {
	private:
		std::random_device rd_;

	public:
		FdmGPU(T const &terminationTime, T correlation = 0.0, std::size_t numberSteps = 360)
			:FdmBuilderGPU<2, T>{ terminationTime,correlation,numberSteps } {}

		FdmGPU(TimePointsType<T> const &timePoints, T correlation = 0.0)
			:FdmBuilderGPU<2, T>{ timePoints,correlation} {}

		template<typename SDEs>
		PathValuesType<PathValuesType<T>> operator()(SDEs const &sdes, std::size_t iterations,
			FDMScheme scheme = FDMScheme::EulerScheme) {

			T delta = (this->terminationTime_ / static_cast<T>(this->numberSteps_));
			std::size_t seed = 56489;
			switch (scheme) {
			case FDMScheme::EulerScheme:
			{
				if (this->timePointsOn_ == false) {
					EulerSchemeGPU<2, T> euler(delta, this->numberSteps_, this->correlation_);
					return euler.simulate(sdes, iterations, seed);
				}
				else {
					EulerSchemeGPU<2, T> euler(this->correlation_);
					const TimePointsType<T> timePoints = this->timePoints_;
					return euler.simulateWithTimePoints(sdes, iterations, timePoints, seed);
				}
			}
			break;
			case FDMScheme::MilsteinScheme:
			{
				if (this->timePointsOn_ == false) {
					MilsteinSchemeGPU<2, T> milstein(delta, this->numberSteps_, this->correlation_);
					return milstein.simulate(sdes, iterations, seed);
				}
				else {
					MilsteinSchemeGPU<2, T> milstein(this->correlation_);
					const TimePointsType<T> timePoints = this->timePoints_;
					return milstein.simulateWithTimePoints(sdes, iterations, timePoints, seed);
				}
			}
			break;
			default:
			{
				if (this->timePointsOn_ == false) {
					EulerSchemeGPU<2, T> euler(delta, this->numberSteps_, this->correlation_);
					return euler.simulate(sdes, iterations, seed);
				}
				else {
					EulerSchemeGPU<2, T> euler(this->correlation_);
					const TimePointsType<T> timePoints = this->timePoints_;
					return euler.simulateWithTimePoints(sdes, iterations, timePoints, seed);
				}
			}
			}
		}

	};

}





#endif ///_FDM_GPU_H_