#pragma once
#if !defined(_MC_UTILITIES_H_)
#define _MC_UTILITIES_H_

#include<cassert>
#include<functional>
#include<amp.h>
#include<amp_math.h>

namespace mc_utilities {

	#define PI 3.14159265358979323846

	enum class withRespectTo {
		firstArg,
		secondArg,
	};


	template<std::size_t Order,
		typename T,
		typename ...Args>
		struct PartialCentralDifferenceBuilder {
	};


	template<typename T, typename ...Args>
	struct PartialCentralDifferenceBuilder<1, T, Args...> {
	protected:
		T step_ = 10e-6;
	public:
		inline void setStep(double step) { step_ = step; }
		inline double step()const { return step_; }

		virtual std::function<T(Args...)> operator()(std::function<T(Args...)> &&fun, withRespectTo withRespectEnum) = 0;
	};

	template<std::size_t Order,
		std::size_t ArgsCount,
		typename = typename std::enable_if<(Order > 0) && (ArgsCount > 0)>::type>
		struct PartialCentralDifference {

	};

	template<typename T>
	struct PartialCentralDifference<1, 1, T> : public PartialCentralDifferenceBuilder<1, T, T> {
	public:
		std::function<T(T)> operator()(std::function<T(T)> &&fun, withRespectTo withRespectEnum)override {
			return [=](T arg) {
				return ((fun(arg + 0.5 * (this->step_)) - fun(arg - 0.5 * (this->step_))) / (this->step_));
			};
		}
	};

	template<typename T>
	struct PartialCentralDifference<1, 2, T> : public PartialCentralDifferenceBuilder<1, T, T, T> {
	public:
		std::function<T(T,T)> operator()(std::function<T(T,T)> &&fun, withRespectTo withRespectEnum)override {
			if (withRespectEnum == withRespectTo::firstArg) {
				return [=](T arg1, T arg2) {
					return ((fun(arg1 + 0.5 * (this->step_),arg2) - fun(arg1 - 0.5 * (this->step_),arg2)) / (this->step_));
				};
			}
			else if (withRespectEnum == withRespectTo::secondArg) {
				return [=](T arg1, T arg2) {
					return ((fun(arg1,arg2 + 0.5 * (this->step_)) - fun(arg1,arg2 - 0.5 * (this->step_))) / (this->step_));
				};
			}
		}
	};


}



#endif ///_MC_UTILITIES_H_