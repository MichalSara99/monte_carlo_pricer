#pragma once
#if !defined(_FDM_H_)
#define _FDM_H_

#include"mc_types.h"
#include"fdm_scheme.h"
#include"sde.h"

namespace finite_difference_method {

	using mc_types::ISde;
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
			std::size_t numberSteps = 10000)
			:model_{ model },terminationTime_ {terminationTime},numberSteps_{numberSteps}{}

		FdmBuilder(ISde<T,Ts...> const &isde,T const &init, T const &terminationTime,
			std::size_t numberSteps = 10000)
			:model_{ new Sde<T,Ts...>{isde,init} }, terminationTime_{ terminationTime },
			numberSteps_{ numberSteps } {}
	};


	//Finite Difference Method builder for two-factor models
	template<typename T,typename ...Ts>
	class FdmBuilder<2, T, Ts...> {
	protected:
		T terminationTime_;
		std::size_t numberSteps_;
		std::shared_ptr<Sde<T,Ts...>> factor1_;
		std::shared_ptr<Sde<T,Ts...>> factor2_;

	public:
		FdmBuilder(std::tuple<std::shared_ptr<Sde<T,Ts...>>, std::shared_ptr<Sde<T, Ts...>>> const &model,
			T const &terminationTime,std::size_t numberSteps = 10000)
			:factor1_{ std::get<0>(model) },factor2_{std::get<1>(model)},
			terminationTime_{ terminationTime }, numberSteps_{ numberSteps } {}

		FdmBuilder(std::shared_ptr<Sde<T, Ts...>> const &factor1,
			std::shared_ptr<Sde<T, Ts...>> const &factor2,
			T const &terminationTime, std::size_t numberSteps = 10000)
			:factor1_{factor1 }, factor2_{ factor2 },
			terminationTime_{ terminationTime }, numberSteps_{ numberSteps } {}

		FdmBuilder(ISde<T,Ts...> const &isde1, T const &init1,
			ISde<T, Ts...> const &isde2, T const &init2,
			T const &terminationTime,std::size_t numberSteps = 10000)
			:factor1_{ new Sde<T,Ts...>{ isde1,init1 } },
			factor2_{ new Sde<T,Ts...>{ isde2,init2 } },
			terminationTime_{ terminationTime },
			numberSteps_{ numberSteps } {}
	};



	template<std::size_t FactorCount,
		    typename T,
			typename =typename std::enable_if<std::is_arithmetic<T>::value && 
											(FactorCount > 0)>::type>
	class Fdm {
	};

	template<typename T>
	class Fdm<1, T> :public FdmBuilder<1, T, T,T> {
	public:
		Fdm(std::shared_ptr<Sde<T, T, T>> const &model, T const &terminationTime,
			std::size_t numberSteps = 10000)
			:FdmBuilder<1,T,T,T>{model,terminationTime,numberSteps}{}
		Fdm(ISde<T, T,T> const &isde, T const &init, T const &terminationTime,
			std::size_t numberSteps = 10000)
			:FdmBuilder<1,T,T,T>{isde,init,terminationTime,numberSteps}{}



	};

	template<typename T>
	class Fdm<2, T> :public FdmBuilder<2, T, T, T,T> {
	public:
		Fdm(std::tuple<std::shared_ptr<Sde<T, T,T,T>>, std::shared_ptr<Sde<T, T,T,T>>> const &model,
			T const &terminationTime, std::size_t numberSteps = 10000)
			:FdmBuilder<2, T, T,T,T>{ model,terminationTime,numberSteps } {}

		Fdm(std::shared_ptr<Sde<T, T, T, T>> const &factor1,
			std::shared_ptr<Sde<T, T, T, T>> const &factor2,
			T const &terminationTime, std::size_t numberSteps = 10000)
			:FdmBuilder<2, T, T, T,T>{ factor1,factor2,terminationTime,numberSteps } {}

		Fdm(ISde<T, T,T,T> const &isde1, T const &init1,
			ISde<T, T,T,T> const &isde2, T const &init2,
			T const &terminationTime, std::size_t numberSteps = 10000)
			:FdmBuilder<2, T, T, T,T>{ isde1,init1,isde2,init2,terminationTime,numberSteps } {}


	};


}





#endif ///_FDM_H_