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


	class PlainCallStrategy: public PayoffStrategy<double> {
	private:
		double strike_;
	public:
		PlainCallStrategy(double strike)
			:strike_{ strike }{}

		double payoff(double const &underlying)const override {
			return std::max(0.0, underlying - strike_);
		}

	};

	class PlainPutStrategy :public PayoffStrategy<double> {
	private:
		double strike_;
	public:
		PlainPutStrategy(double strike)
			:strike_{ strike } {}

		double payoff(double const &underlying)const override {
			return std::max(0.0, strike_ - underlying);
		}
	};

	class AsianAvgCallStrategy :public PayoffStrategy<PathValuesType<double>> {
	private:
		double strike_;
	public:
		AsianAvgCallStrategy(double strike)
			:strike_{strike}{}

		double payoff(PathValuesType<double> const &underlying)const {
			std::size_t N = underlying.size();
			auto sum = std::accumulate(underlying.begin(), underlying.end(), 0.0);
			auto avg = (sum / static_cast<double>(N));
			return std::max(0.0, avg - strike_);
		}
	};

	class AsianAvgPutStrategy :public PayoffStrategy<PathValuesType<double>> {
	private:
		double strike_;
	public:
		AsianAvgPutStrategy(double strike)
			:strike_{ strike } {}

		double payoff(PathValuesType<double> const &underlying)const {
			std::size_t N = underlying.size();
			auto sum = std::accumulate(underlying.begin(), underlying.end(), 0.0);
			auto avg = (sum / static_cast<double>(N));
			return std::max(0.0, strike_ - avg);
		}
	};


}



#endif ///_PAYOFF_STARTEGY_H_