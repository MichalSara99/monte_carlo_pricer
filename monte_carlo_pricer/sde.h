#pragma once
#if !defined(_SDE_H_)
#define _SDE_H_

#include"mc_types.h"

namespace sde {

	using mc_types::ISde;
	using mc_types::SdeComponent;

	template<typename T = double,
			typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	class Sde {
	private:
		T initCond_;
		SdeComponent<T> drift_;
		SdeComponent<T> diffusion_;

	public:
		Sde() = default;
		Sde(ISde<T> const &sdeComponents,T const &initialCondition = 0.0)
			:drift_{std::get<0>(sdeComponents)},diffusion_{std::get<1>(sdeComponents)},
			initCond_{ initialCondition}{}

		Sde(Sde<T> const &copy)
			:drift_{copy.drift_},diffusion_{copy.diffusion_},
			initCond_{ copy.initCond_ }{}

		T drift(T const &first,T const &second)const {
			return drift_(first, second);
		}

		T diffusion(T const &first, T const &second)const {
			return diffusion_(first, second);
		}


	};


}



#endif ///_SDE_H_