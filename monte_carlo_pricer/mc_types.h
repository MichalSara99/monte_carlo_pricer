#pragma once
#if !defined(_MC_TYPES_H_)
#define _MC_TYPES_H_

#include<vector>
#include<type_traits>
#include<functional>
#include<tuple>

namespace mc_types {

	template<typename ReturnType,typename ArgType>
	using PayoffFunType = std::function<ReturnType(ArgType)>;

	template<typename T>
	using PathValuesType = std::vector<T>;

	template<typename T,typename ...Args>
	using SdeComponent = std::function<T(Args...)>;

	// For models
	// 0: drift lambda 1: diffusion lambda
	template<typename T,typename ...Ts>
	using ISde = std::tuple<SdeComponent<T,Ts...>, SdeComponent<T,Ts...>>;

	enum class SdeModelType { oneFactor, twoFactor };

}



#endif ///_MC_TYPES_H_