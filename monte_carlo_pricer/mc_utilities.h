#pragma once
#if !defined(_MC_UTILITIES_H_)
#define _MC_UTILITIES_H_

#include<functional>

namespace mc_utilities{


	template<std::size_t Order,
			typename T,
			typename = typename std::enable_if<(std::is_arithmetic<T>::value) &&
												(Order > 0) >::type>
	struct PartialCentralDifference {
	};


	template<typename T>
	struct PartialCentralDifference<1,T> {
	private:
		T step_ = 10e-6;
	public:
		inline void setStep(double step) { step_ = step; }
		inline double step()const { return step_; }

		auto operator()(std::function<T(T)> &&fun) {
			return [=](T arg) {
				return ((fun(arg + 0.5 * (this->step_)) - fun(arg - 0.5 * (this->step_))) / (this->step_));
			};
		}
	};


}



#endif ///_MC_UTILITIES_H_