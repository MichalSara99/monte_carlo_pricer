#pragma once
#if !defined(_SDE_H_)
#define _SDE_H_

#include"mc_types.h"

namespace sde {

	using mc_types::ISde1;
	using mc_types::ISde2;
	using mc_types::SdeComponent2Args;
	using mc_types::SdeComponent3Args;

	template<std::size_t Factor,
			typename T = double,
			typename =typename std::enable_if<(std::is_arithmetic<T>::value) && 
											(Factor > 0)>::type>
	class Sde {

	};


	template<typename T>
	class Sde<1,T> {
	private:
		T initCond_;
		SdeComponent2Args<T> drift_;
		SdeComponent2Args<T> diffusion_;

	public:
		Sde() = default;
		Sde(ISde1<T> const &sdeComponents,T const &initialCondition = 0.0)
			:drift_{std::get<0>(sdeComponents)},diffusion_{std::get<1>(sdeComponents)},
			initCond_{ initialCondition}{}

		Sde(Sde<1,T> const &copy)
			:drift_{copy.drift_},diffusion_{copy.diffusion_},
			initCond_{ copy.initCond_ }{}

		T drift(T const &first,T const &second)const {
			return drift_(first, second);
		}

		T diffusion(T const &first, T const &second)const {
			return diffusion_(first, second);
		}
	};

	template<typename T>
	class Sde<2, T> {
	private:
		T initCond_;
		SdeComponent3Args<T> drift_;
		SdeComponent3Args<T> diffusion_;

	public:
		Sde() = default;
		Sde(ISde2<T> const &sdeComponents,  T const &initialCondition = 0.0)
			:drift_{ std::get<0>(sdeComponents) }, diffusion_{ std::get<1>(sdeComponents) },
			initCond_{ initialCondition } {}

		Sde(Sde<2, T> const &copy)
			:drift_{ copy.drift_ }, diffusion_{ copy.diffusion_ },
			initCond_{ copy.initCond_ } {}

		T drift(T const &first, T const &second, T const &third)const {
			return drift_(first, second, third);
		}

		T diffusion(T const &first, T const &second,T const &third)const {
			return diffusion_(first, second, third);
		}

	};


}



#endif ///_SDE_H_