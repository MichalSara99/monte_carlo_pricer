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

	// For 2 parameter function
	template<typename T>
	using SdeComponent2Args = std::function<T(T const &, T const &)>;
	
	//For 3 parameter function
	template<typename T>
	using SdeComponent3Args = std::function<T(T const &, T const &, T const &)>;

	// For one-factor models
	// 0: drift lambda 1: diffusion lambda
	template<typename T>
	using ISde1 = std::tuple<SdeComponent2Args<T>, SdeComponent2Args<T>>;

	// For two-factor models
	// 0: drift lambda 1: diffusion lambda
	template<typename T>
	using ISde2 = std::tuple<SdeComponent3Args<T>, SdeComponent3Args<T>>;

	enum class SdeModelType { oneFactor, twoFactor };

}



#endif ///_MC_TYPES_H_