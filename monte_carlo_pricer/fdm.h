#pragma once
#if !defined(_FDM_H_)
#define _FDM_H_

#include"mc_types.h"
#include"fdm_scheme.h"
#include"sde.h"
#include<thread>
#include<future>

namespace finite_difference_method {

	using mc_types::ISde;
	using mc_types::PathValuesType;
	using mc_types::TimePointsType;
	using mc_types::FDMScheme;
	using sde::Sde;


	template<std::size_t FactorCount,typename T,typename ...Ts>
	class FdmBuilder {

	};


	// Finite Difference Method builder for one-factor models
	template<typename T,typename ...Ts>
	class FdmBuilder<1,T,Ts...> {
	protected:
		T terminationTime_;
		std::size_t numberSteps_;
		std::shared_ptr<Sde<T,Ts...>> model_;

	public:
		FdmBuilder(std::shared_ptr<Sde<T,Ts...>> const &model,T const &terminationTime,
			std::size_t numberSteps = 360)
			:model_{ model },terminationTime_ {terminationTime},numberSteps_{numberSteps}{}

		FdmBuilder(ISde<T,Ts...> const &isde,T const &init, T const &terminationTime,
			std::size_t numberSteps = 360)
			:model_{ new Sde<T,Ts...>{isde,init} }, terminationTime_{ terminationTime },
			numberSteps_{ numberSteps } {}

		TimePointsType<T> timeResolution()const {
			TimePointsType<T> points(numberSteps_ + 1);
			auto delta = (terminationTime_ / static_cast<T>(numberSteps_));
			std::size_t i = 0;
			std::generate(points.begin(), points.end(), [&i, &delta]() {return delta * i++; });
			return points;
		}

		virtual PathValuesType<PathValuesType<T>> operator()(std::size_t iterations,FDMScheme scheme = FDMScheme::EulerScheme)=0;

	};


	// Finite Difference Method builder for two-factor models:
	// Factor1 always depends on factor2 (process driven by factor2 enters factor 1)
	template<typename T,typename ...Ts>
	class FdmBuilder<2, T, Ts...> {
	protected:
		T terminationTime_;
		T correlation_;
		std::size_t numberSteps_;
		std::shared_ptr<Sde<T,Ts...>> factor1_;
		std::shared_ptr<Sde<T,Ts...>> factor2_;

	public:
		FdmBuilder(std::tuple<std::shared_ptr<Sde<T,Ts...>>, std::shared_ptr<Sde<T, Ts...>>> const &model,
			T const &terminationTime,T const &correlation = 0.0,std::size_t numberSteps = 360)
			:factor1_{ std::get<0>(model) },factor2_{std::get<1>(model)},
			terminationTime_{ terminationTime }, correlation_{ correlation },
			numberSteps_ {numberSteps} {}

		FdmBuilder(std::shared_ptr<Sde<T, Ts...>> const &factor1,
			std::shared_ptr<Sde<T, Ts...>> const &factor2,
			T const &terminationTime, T const &correlation = 0.0, std::size_t numberSteps = 360)
			:factor1_{factor1 }, factor2_{ factor2 },
			terminationTime_{ terminationTime },correlation_{correlation},
			numberSteps_{ numberSteps } {}

		FdmBuilder(ISde<T,Ts...> const &isde1, T const &init1,
			ISde<T, Ts...> const &isde2, T const &init2,
			T const &terminationTime,T const & correlation=0.0,std::size_t numberSteps = 360)
			:factor1_{ new Sde<T,Ts...>{ isde1,init1 } },
			factor2_{ new Sde<T,Ts...>{ isde2,init2 } },
			terminationTime_{ terminationTime },
			correlation_{correlation},
			numberSteps_{ numberSteps } {}

		TimePointsType<T> timeResolution()const {
			TimePointsType<T> points(numberSteps_ + 1);
			auto delta = (terminationTime_ / static_cast<T>(numberSteps_));
			std::size_t i = 0;
			std::generate(points.begin(), points.end(), [&i, &delta]() {return delta * i++; });
			return points;
		}

		virtual PathValuesType<PathValuesType<T>> operator()(std::size_t iterations,FDMScheme scheme = FDMScheme::EulerScheme) = 0;
	};



	template<std::size_t FactorCount,
		    typename T,
			typename =typename std::enable_if<std::is_arithmetic<T>::value && 
											(FactorCount > 0)>::type>
	class Fdm {
	};

	template<typename T>
	class Fdm<1, T> :public FdmBuilder<1, T, T,T> {
	private:
		std::random_device rd_;
	public:
		Fdm(std::shared_ptr<Sde<T, T, T>> const &model, T const &terminationTime,
			std::size_t numberSteps = 360)
			:FdmBuilder<1,T,T,T>{model,terminationTime,numberSteps}{}
		Fdm(ISde<T, T,T> const &isde, T const &init, T const &terminationTime,
			std::size_t numberSteps = 360)
			:FdmBuilder<1,T,T,T>{isde,init,terminationTime,numberSteps}{}

		PathValuesType<PathValuesType<T>> operator()(std::size_t iterations,
													FDMScheme scheme = FDMScheme::EulerScheme)override{

			asyncKernel<T> asyncGenerator = nullptr;
			T delta = (this->terminationTime_ / static_cast<T>(this->numberSteps_));

			switch (scheme) {
			case FDMScheme::EulerScheme:
			{
				EulerScheme<1, T> euler(this->model_, delta, this->numberSteps_);
				asyncGenerator = euler;
			}
			break;
			case FDMScheme::MilsteinScheme:
			{
				MilsteinScheme<1, T> milstein(this->model_, delta, this->numberSteps_);
				asyncGenerator = milstein;
			}
			break;
			}

			PathValuesType<std::future<PathValuesType<T>>> futures;
			futures.reserve(iterations);

			for (std::size_t i = 0; i < iterations; ++i) {
				futures.emplace_back(std::async(std::launch::async, asyncGenerator, rd_()));
			}

			PathValuesType<PathValuesType<T>> paths;
			paths.reserve(iterations);

			for (auto &path : futures) {
				paths.emplace_back(std::move(path.get()));
			}

			return paths;
		}

	};

	template<typename T>
	class Fdm<2, T> :public FdmBuilder<2, T, T, T,T> {
	public:
		Fdm(std::tuple<std::shared_ptr<Sde<T, T,T,T>>, std::shared_ptr<Sde<T, T,T,T>>> const &model,
			T const &terminationTime,T const &correlation = 0.0, std::size_t numberSteps = 360)
			:FdmBuilder<2, T, T,T,T>{ model,terminationTime,correlation,numberSteps } {}

		Fdm(std::shared_ptr<Sde<T, T, T, T>> const &factor1,
			std::shared_ptr<Sde<T, T, T, T>> const &factor2,
			T const &terminationTime, T const &correlation = 0.0,std::size_t numberSteps = 360)
			:FdmBuilder<2, T, T, T,T>{ factor1,factor2,terminationTime,correlation,numberSteps } {}

		Fdm(ISde<T, T,T,T> const &isde1, T const &init1,
			ISde<T, T,T,T> const &isde2, T const &init2,
			T const &terminationTime,T const &correlation = 0.0,
			std::size_t numberSteps = 360)
			:FdmBuilder<2, T, T, T,T>{ isde1,init1,isde2,init2,terminationTime,
			correlation,numberSteps } {}


		PathValuesType<PathValuesType<T>> operator()(std::size_t iterations,FDMScheme scheme = FDMScheme::EulerScheme) override {
			throw std::exception("Not yet implemented");
		}

	};


}





#endif ///_FDM_H_