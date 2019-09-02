#pragma once
#if !defined(_SDE_H_)
#define _SDE_H_

#include"mc_types.h"

namespace sde {

	using mc_types::ISde;
	using mc_types::SdeComponent;

	template<typename T,typename ...Ts>
	class Sde {
	private:
		T initCond_;
		SdeComponent<T,Ts...> drift_;
		SdeComponent<T,Ts...> diffusion_;

	public:
		Sde(ISde<T,Ts...> const &sdeComponents,  T const &initialCondition = 0.0)
			:drift_{ std::get<0>(sdeComponents) }, diffusion_{ std::get<1>(sdeComponents) },
			initCond_{ initialCondition } {}

		Sde(Sde<T,Ts...> const &copy)
			:drift_{ copy.drift_ }, diffusion_{ copy.diffusion_ },
			initCond_{ copy.initCond_ } {}

		T drift(Ts...args)const {
			return drift_(args...);
		}

		T diffusion(Ts...args)const {
			return diffusion_(args...);
		}

	};


}



#endif ///_SDE_H_