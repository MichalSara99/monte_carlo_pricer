#pragma once
#if !defined(_PAYOFF_H_)
#define _PAYOFF_H_

#include"payoff_strategy.h"

namespace payoff {

	using mc_types::PayoffFunType;

	template<typename UnderlyingType,
			typename = typename std::enable_if<std::is_scalar<UnderlyingType>::value||
												std::is_compound<UnderlyingType>::value>::type>
	class Payoff {
	private:
		PayoffFunType<double,UnderlyingType> payoff_;
	public:
		typedef UnderlyingType underlyingType;

		Payoff(PayoffFunType<double,UnderlyingType> &&fun)
			:payoff_{std::forward<PayoffFunType<double,UnderlyingType>>(fun)}{}

		virtual ~Payoff(){}

		Payoff(Payoff<UnderlyingType> const &copy)
			:payoff_{copy.payoff_}{}

		Payoff& operator=(Payoff<UnderlyingType> const &copy) {
			if (this != &copy ) {
				payoff_ = copy.payoff_;
			}
			return *this;
		}

		double payoff(UnderlyingType const &underlying)const {
			return payoff_(underlying);
		}




	};



}




#endif ///_PAYOFF_H_

