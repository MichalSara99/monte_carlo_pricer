#pragma once
#if !defined(_PAYOFF_STARTEGY_H_)
#define _PAYOFF_STARTEGY_H_

#include"mc_types.h"
#include<algorithm>
#include<numeric>

namespace payoff {

	using mc_types::PathValuesType;


	template<typename UnderlyingType,
			typename = typename std::enable_if<std::is_scalar<UnderlyingType>::value||
											  std::is_compound<UnderlyingType>::value>::type>
	class PayoffStrategy {
	public:
		typedef UnderlyingType underlyingType;
		virtual double payoff(UnderlyingType const &underlying)const = 0;
	};

	template<typename T=double>
	class PlainCallStrategy: public PayoffStrategy<T> {
	private:
		T strike_;
	public:
		PlainCallStrategy(T strike)
			:strike_{ strike }{}

		double payoff(T const &underlying)const override {
			return static_cast<double>(std::max(0.0, underlying - strike_));
		}

	};

	template<typename T=double>
	class PlainPutStrategy :public PayoffStrategy<T> {
	private:
		T strike_;
	public:
		PlainPutStrategy(T strike)
			:strike_{ strike } {}

		double payoff(T const &underlying)const override {
			return static_cast<double>(std::max(0.0, strike_ - underlying));
		}
	};

	template<typename T = double>
	class AsianAvgCallStrategy :public PayoffStrategy<PathValuesType<T>> {
	private:
		T strike_;
	public:
		AsianAvgCallStrategy(T strike)
			:strike_{strike}{}

		double payoff(PathValuesType<T> const &underlying)const {
			std::size_t N = underlying.size();
			auto sum = std::accumulate(underlying.begin(), underlying.end(), 0.0);
			auto avg = (sum / static_cast<double>(N));
			return static_cast<double>(std::max(0.0, avg - strike_));
		}
	};

	template<typename T = double>
	class AsianAvgPutStrategy :public PayoffStrategy<PathValuesType<T>> {
	private:
		T strike_;
	public:
		AsianAvgPutStrategy(T strike)
			:strike_{ strike } {}

		double payoff(PathValuesType<T> const &underlying)const {
			std::size_t N = underlying.size();
			auto sum = std::accumulate(underlying.begin(), underlying.end(), 0.0);
			auto avg = (sum / static_cast<double>(N));
			return static_cast<double>(std::max(0.0, strike_ - avg));
		}
	};


}



#endif ///_PAYOFF_STARTEGY_H_