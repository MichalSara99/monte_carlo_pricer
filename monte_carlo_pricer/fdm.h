#pragma once
#if !defined(_FDM_H_)
#define _FDM_H_

#include"mc_types.h"
#include"fdm_scheme.h"
#include"sde.h"

namespace finite_difference_method {

	using mc_types::ISde;
	using sde::Sde;

	template<std::size_t Factor,
			typename T = double,
			typename =typename std::enable_if<(std::is_arithmetic<T>::value) && 
											(Factor > 0)>::type>
	class Fdm {

	};

	// Finite Difference Method for one-factor models
	template<typename T>
	class Fdm<1,T> {
	private:
		T terminationTime_;
		std::size_t numberSteps_;
		std::shared_ptr<Sde<1,T>> model_;

	public:
		Fdm(std::shared_ptr<Sde<1,T>> const &model,T const &terminationTime,
			std::size_t numberSteps = 10000)
			:model_{ model },terminationTime_ {terminationTime},numberSteps_{numberSteps}{}

		Fdm(ISde1<T> const &isde,T const &init, T const &terminationTime,
			std::size_t numberSteps = 10000)
			:model_{ new Sde<T>{isde,init} }, terminationTime_{ terminationTime },
			numberSteps_{ numberSteps } {}



	};

	//Finite Difference Method for two - factor models
	template<typename T>
	class Fdm<2, T> {
	private:
		T terminationTime_;
		std::size_t numberSteps_;
		std::shared_ptr<Sde<2,T>> model1_;
		std::shared_ptr<Sde<2,T>> model2_;

	public:
		Fdm(std::shared_ptr<Sde<2,T>> const &model1,
			std::shared_ptr<Sde<2, T>> const &model2,
			T const &terminationTime,
			std::size_t numberSteps = 10000)
			:model1_{ model1 }, model2_{model2},
			terminationTime_{ terminationTime }, numberSteps_{ numberSteps } {}

		Fdm(ISde2<T> const &isde1, T const &init1,
			ISde2<T> const &isde2, T const &init2,
			T const &terminationTime,
			std::size_t numberSteps = 10000)
			:model1_{ new Sde<2,T>{ isde1,init1 } },
			model2_{ new Sde<2,T>{ isde2,init2 } },
			terminationTime_{ terminationTime },
			numberSteps_{ numberSteps } {}




	};



}





#endif ///_FDM_H_