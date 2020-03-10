#pragma once
#if !defined(_FDM_SCHEME_GPU_H_)
#define _FDM_SCHEME_GPU_H_

#include"mc_types.h"
#include"mc_utilities.h"
#include<amp.h>
#include<random>
#include<amp_math.h>
#include"amp_sobol_rng.h"

namespace finite_difference_method_gpu {

	using mc_types::TimePointsType;
	using mc_types::PathValuesType;
	using mc_utilities::PartialCentralDifference;
	using mc_utilities::withRespectTo;
	using mc_utilities::NormalVariate;



	template<std::size_t FactorCount, typename T>
	class SchemeBuilderGPU {
	public:
		virtual ~SchemeBuilderGPU(){}
	};

	// Scheme builder for one-factor models:
	template<typename T>
	class SchemeBuilderGPU<1, T> {
	protected:
		std::size_t numberSteps_;
		T delta_;
	public:
		SchemeBuilderGPU(){}
		SchemeBuilderGPU(T delta,std::size_t numberSteps)
			:delta_{delta},numberSteps_{ numberSteps }{}
	};


	// Scheme builder for two-factor models:
	template<typename T>
	class SchemeBuilderGPU<2, T> {
	protected:
		std::size_t numberSteps_;
		T delta_;
		T correlation_;

	public:
		SchemeBuilderGPU(T delta,std::size_t numberSteps,T correlation)
			:delta_{delta},
			correlation_{correlation},
			numberSteps_{numberSteps}{}

		SchemeBuilderGPU(T correlation)
			:correlation_{correlation}{}

	};


	template<std::size_t FactorCount,
			typename T,
			typename =typename std::enable_if<(std::is_floating_point<T>::value) && 
					(FactorCount>0)>::type>
	class EulerSchemeGPU{};

	template<typename T>
	class EulerSchemeGPU<1, T> :public SchemeBuilderGPU<1, T> {
	public:
		EulerSchemeGPU()
			:SchemeBuilderGPU<1,T>{}{}
		EulerSchemeGPU(T delta,std::size_t numberSteps)
			:SchemeBuilderGPU<1,T>{delta,numberSteps}{}

	private:
		template<typename SDE>
		PathValuesType<PathValuesType<T>> _generator(SDE const &sde,std::size_t iterations, std::size_t numberSteps,
			T delta, std::random_device::result_type seed)const {
			
			// retrieve the sde components:
			auto drift = std::get<0>(sde);
			auto diffusion = std::get<1>(sde);
			auto init = std::get<2>(sde);

			auto randExtent = concurrency::extent<1>(iterations);
			sobol_rng_collection<sobol_rng<1>, 1> sbRand(randExtent, seed);
			auto simExtent = concurrency::extent<2>(iterations, numberSteps);
			concurrency::array_view<T, 2> sims(simExtent);

			NormalVariate<T> normalRand(0.0, 1.0);

			concurrency::parallel_for_each(randExtent, [=](concurrency::index<1> simIdx)restrict(amp) {
				auto rnd = sbRand[simIdx];
				rnd.skip(sbRand.direction_numbers(), simIdx[0]);
				sims(simIdx[0], 0) = init;
				for (std::size_t t = 1; t < numberSteps; ++t) {
					sims(simIdx[0], t) = sims(simIdx[0], t - 1) +
						drift(t - 1, sims(simIdx[0], t - 1))*delta +
						diffusion(t - 1, sims(simIdx[0], t - 1))*concurrency::fast_math::sqrt(delta)*
						normalRand(rnd.get_single(1), rnd.get_single(1));
				}
			});
			sims.synchronize();

			PathValuesType<PathValuesType<T>> paths(iterations);
			for (std::size_t i = 0; i < iterations; ++i) {
				PathValuesType<T> path(numberSteps);
				for (std::size_t s = 0; s < numberSteps; ++s) {
					path[s] = std::move(sims(i, s));
				}
				paths[i] = std::move(path);
			}
			return paths;
		}

		template<typename SDE>
		PathValuesType<PathValuesType<T>> _generatorWithTimePoints(SDE const &sde, std::size_t iterations,
			TimePointsType<T> const &timePoints, std::random_device::result_type seed)const {

			// retrieve the sde components:
			auto drift = std::get<0>(sde);
			auto diffusion = std::get<1>(sde);
			auto init = std::get<2>(sde);


			std::size_t const numberSteps = timePoints.size();
			auto randExtent = concurrency::extent<1>(iterations);
			sobol_rng_collection<sobol_rng<1>, 1> sbRand(randExtent, seed);
			auto simExtent = concurrency::extent<2>(iterations, numberSteps);
			concurrency::array_view<T, 2> sims(simExtent);
			concurrency::array_view<const T, 1> tmps(numberSteps, timePoints);

			NormalVariate<T> normalRand(0.0, 1.0);

			concurrency::parallel_for_each(randExtent, [=](concurrency::index<1> simIdx)restrict(amp) {
				auto rnd = sbRand[simIdx];
				rnd.skip(sbRand.direction_numbers(), simIdx[0]);
				sims(simIdx[0], 0) = init;
				for (std::size_t t = 1; t < numberSteps; ++t) {
					sims(simIdx[0], t) = sims(simIdx[0], t - 1) +
						drift(tmps[t - 1], sims(simIdx[0], t - 1))*(tmps[t] - tmps[t - 1]) +
						diffusion(tmps[t - 1], sims(simIdx[0], t - 1))*concurrency::fast_math::sqrt(tmps[t] - tmps[t - 1])*
						normalRand(rnd.get_single(1), rnd.get_single(1));
				}
			});
			sims.synchronize();

			PathValuesType<PathValuesType<T>> paths(iterations);
			for (std::size_t i = 0; i < iterations; ++i) {
				PathValuesType<T> path(numberSteps);
				for (std::size_t s = 0; s < numberSteps; ++s) {
					path[s] = std::move(sims(i, s));
				}
				paths[i] = std::move(path);
			}
			return paths;
		}


	public:
		template<typename SDE>
		PathValuesType<PathValuesType<T>> simulate(SDE const &sde, std::size_t iterations, 
			std::random_device::result_type seed) const {
			auto numberSteps = this->numberSteps_;
			auto delta = this->delta_;
			return _generator(sde, iterations, numberSteps, delta, seed);
		}

		template<typename SDE>
		PathValuesType<PathValuesType<T>> simulateWithTimePoints(SDE const &sde, std::size_t iterations,
			TimePointsType<T> const &timePoints, std::random_device::result_type seed) const {
			return _generatorWithTimePoints(sde, iterations, timePoints, seed);
		}





	};
	





}



#endif ///_FDM_SCHEME_GPU_H_