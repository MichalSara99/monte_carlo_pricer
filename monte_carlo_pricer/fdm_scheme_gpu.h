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
				T z{};
				rnd.skip(sbRand.direction_numbers(), simIdx[0]);
				sims(simIdx[0], 0) = init;
				for (std::size_t t = 1; t < numberSteps; ++t) {
					z = normalRand(rnd.get_single(1), rnd.get_single(1));
					sims(simIdx[0], t) = sims(simIdx[0], t - 1) +
						drift((t - 1)*delta, sims(simIdx[0], t - 1))*delta +
						diffusion((t - 1)*delta, sims(simIdx[0], t - 1))*concurrency::fast_math::sqrt(delta)*z;
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
				T z{};
				rnd.skip(sbRand.direction_numbers(), simIdx[0]);
				sims(simIdx[0], 0) = init;
				for (std::size_t t = 1; t < numberSteps; ++t) {
					z = normalRand(rnd.get_single(1), rnd.get_single(1));
					sims(simIdx[0], t) = sims(simIdx[0], t - 1) +
						drift(tmps[t - 1], sims(simIdx[0], t - 1))*(tmps[t] - tmps[t - 1]) +
						diffusion(tmps[t - 1], sims(simIdx[0], t - 1))*concurrency::fast_math::sqrt(tmps[t] - tmps[t - 1])*z;
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
	


	template<typename T>
	class EulerSchemeGPU<2, T> :public SchemeBuilderGPU<2, T> {
	public:
		explicit EulerSchemeGPU(T delta,std::size_t numberSteps,T correlation)
			:SchemeBuilderGPU<2, T>{delta, numberSteps, correlation}
		{}

		explicit EulerSchemeGPU<2,T>(T correlation)
			: SchemeBuilderGPU<2, T>{correlation } {}
	private:
		template<typename SDEs>
		PathValuesType<PathValuesType<T>> _generator(SDEs const &sdes, std::size_t iterations, std::size_t numberSteps,
			T delta,T correlation, std::random_device::result_type seed)const {

			// retrieve the sde components:
			auto drift1 = std::get<0>(std::get<0>(sdes));
			auto diffusion1 = std::get<1>(std::get<0>(sdes));
			auto init1 = std::get<2>(std::get<0>(sdes));

			auto drift2 = std::get<0>(std::get<1>(sdes));
			auto diffusion2 = std::get<1>(std::get<1>(sdes));
			auto init2 = std::get<2>(std::get<1>(sdes));

			auto randExtent = concurrency::extent<1>(iterations);
			sobol_rng_collection<sobol_rng<1>, 1> sbRand(randExtent, seed);
			auto simExtent = concurrency::extent<2>(iterations, numberSteps);
			concurrency::array_view<T, 2> sims(simExtent);

			NormalVariate<T> normalRand(0.0, 1.0);

			concurrency::parallel_for_each(randExtent, [=](concurrency::index<1> simIdx)restrict(amp) {
				auto rnd = sbRand[simIdx];
				T firstSpot{ init1 };
				T secondSpot{ init2 };
				T firstSpotNew{};
				T secondSpotNew{};
				T z1{};
				T z2{};
				rnd.skip(sbRand.direction_numbers(), simIdx[0]);

				sims(simIdx[0], 0) = init1;
				for (std::size_t t = 1; t < numberSteps; ++t) {
					z1 = normalRand(rnd.get_single(1), rnd.get_single(1));
					z2 = normalRand(rnd.get_single(1), rnd.get_single(1));

					firstSpotNew = firstSpot +
						drift1((t - 1)*delta, firstSpot, secondSpot) * delta +
						diffusion1((t - 1)*delta, firstSpot, secondSpot) *
						concurrency::fast_math::sqrt(delta) * z1;

					secondSpotNew = secondSpot +
						drift2((t - 1)*delta, firstSpot, secondSpot)*delta +
						diffusion2((t - 1)*delta, firstSpot, secondSpot)*
						concurrency::fast_math::sqrt(delta)*
						(correlation*z1 + concurrency::fast_math::sqrt(1.0 - (correlation*correlation))*z2);
					sims(simIdx[0], t) = firstSpotNew;
					firstSpot = firstSpotNew;
					secondSpot = secondSpotNew;
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

		template<typename SDEs>
		PathValuesType<PathValuesType<T>> _generatorWithTimePoints(SDEs const &sdes, std::size_t iterations,
			TimePointsType<T> const &timePoints,T correlation, std::random_device::result_type seed)const {

			// retrieve the sde components:
			auto drift1 = std::get<0>(std::get<0>(sdes));
			auto diffusion1 = std::get<1>(std::get<0>(sdes));
			auto init1 = std::get<2>(std::get<0>(sdes));

			auto drift2 = std::get<0>(std::get<1>(sdes));
			auto diffusion2 = std::get<1>(std::get<1>(sdes));
			auto init2 = std::get<2>(std::get<1>(sdes));


			std::size_t const numberSteps = timePoints.size();
			auto randExtent = concurrency::extent<1>(iterations);
			sobol_rng_collection<sobol_rng<1>, 1> sbRand(randExtent, seed);
			auto simExtent = concurrency::extent<2>(iterations, numberSteps);
			concurrency::array_view<T, 2> sims(simExtent);
			concurrency::array_view<const T, 1> tmps(numberSteps, timePoints);

			NormalVariate<T> normalRand(0.0, 1.0);

			concurrency::parallel_for_each(randExtent, [=](concurrency::index<1> simIdx)restrict(amp) {
				auto rnd = sbRand[simIdx];
				T firstSpot{ init1 };
				T secondSpot{ init2 };
				T firstSpotNew{};
				T secondSpotNew{};
				T z1{};
				T z2{};
				rnd.skip(sbRand.direction_numbers(), simIdx[0]);

				sims(simIdx[0], 0) = init1;
				for (std::size_t t = 1; t < numberSteps; ++t) {
					z1 = normalRand(rnd.get_single(1), rnd.get_single(1));
					z2 = normalRand(rnd.get_single(1), rnd.get_single(1));

					firstSpotNew = firstSpot +
						drift1(tmps[t-1], firstSpot, secondSpot) * (tmps[t]- tmps[t - 1]) +
						diffusion1(tmps[t - 1], firstSpot, secondSpot) *
						concurrency::fast_math::sqrt(tmps[t]- tmps[t - 1]) * z1;

					secondSpotNew = secondSpot +
						drift2(tmps[t - 1], firstSpot, secondSpot)*(tmps[t]- tmps[t - 1]) +
						diffusion2(tmps[t - 1], firstSpot, secondSpot)*
						concurrency::fast_math::sqrt(tmps[t]- tmps[t - 1])*
						(correlation*z1 + concurrency::fast_math::sqrt(1.0 - (correlation*correlation))*z2);
					sims(simIdx[0], t) = firstSpotNew;
					firstSpot = firstSpotNew;
					secondSpot = secondSpotNew;
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
		template<typename SDEs>
		PathValuesType<PathValuesType<T>> simulate(SDEs const &sdes, std::size_t iterations,
			std::random_device::result_type seed) const {
			auto numberSteps = this->numberSteps_;
			auto delta = this->delta_;
			auto correlation = this->correlation_;
			return _generator(sdes, iterations, numberSteps, delta, correlation, seed);
		}

		template<typename SDEs>
		PathValuesType<PathValuesType<T>> simulateWithTimePoints(SDEs const &sdes, std::size_t iterations,
			TimePointsType<T> const &timePoints, std::random_device::result_type seed) const {
			auto correlation = this->correlation_;
			return _generatorWithTimePoints(sdes, iterations, timePoints,correlation, seed);
		}
	};

	template<std::size_t FactorCount,
			typename T,
			typename =typename std::enable_if<(std::is_floating_point<T>::value) &&
							(FactorCount > 0)>::type>
	class MilsteinSchemeGPU{};


	template<typename T>
	class MilsteinSchemeGPU<1, T> :public SchemeBuilderGPU<1, T> {
	private:
		T step_ = 10e-6;
	public:
		MilsteinSchemeGPU()
			:SchemeBuilderGPU<1, T>{} {}
		MilsteinSchemeGPU(T delta, std::size_t numberSteps)
			:SchemeBuilderGPU<1, T>{ delta,numberSteps } {}

	private:
		template<typename SDE>
		PathValuesType<PathValuesType<T>> _generator(SDE const &sde, std::size_t iterations, std::size_t numberSteps,
			T delta,T step, std::random_device::result_type seed)const {

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
				T z{};
				rnd.skip(sbRand.direction_numbers(), simIdx[0]);
				sims(simIdx[0], 0) = init;
				for (std::size_t t = 1; t < numberSteps; ++t) {
					z = normalRand(rnd.get_single(1), rnd.get_single(1));

					sims(simIdx[0], t) = sims(simIdx[0], t - 1) +
						drift((t - 1)*delta, sims(simIdx[0], t - 1))*delta +
						diffusion((t - 1)*delta, sims(simIdx[0], t - 1))*
						concurrency::fast_math::sqrt(delta)*z +
						0.5*diffusion((t - 1)*delta, sims(simIdx[0], t - 1))*
						((diffusion((t - 1)*delta, sims(simIdx[0], t - 1) + 0.5*step) -
							diffusion((t - 1)*delta, sims(simIdx[0], t - 1) - 0.5*step)) / step)*
							((concurrency::fast_math::sqrt(delta)*z)*
						(concurrency::fast_math::sqrt(delta)*z) - delta);
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
			TimePointsType<T> const &timePoints,T step, std::random_device::result_type seed)const {

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
				T z{};
				rnd.skip(sbRand.direction_numbers(), simIdx[0]);
				sims(simIdx[0], 0) = init;
				for (std::size_t t = 1; t < numberSteps; ++t) {
					z = normalRand(rnd.get_single(1), rnd.get_single(1));
					sims(simIdx[0], t) = sims(simIdx[0], t - 1) +
						drift(tmps[t - 1], sims(simIdx[0], t - 1))*(tmps[t]- tmps[t - 1]) +
						diffusion(tmps[t - 1], sims(simIdx[0], t - 1))*
						concurrency::fast_math::sqrt(tmps[t]- tmps[t - 1])*z +
						0.5*diffusion(tmps[t - 1], sims(simIdx[0], t - 1))*
						((diffusion(tmps[t - 1], sims(simIdx[0], t - 1) + 0.5*step) -
							diffusion(tmps[t - 1], sims(simIdx[0], t - 1) - 0.5*step)) / (step))*
							((concurrency::fast_math::sqrt(tmps[t]- tmps[t - 1])*z)*
						(concurrency::fast_math::sqrt(tmps[t]- tmps[t - 1])*z) - (tmps[t]- tmps[t - 1]));
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
			auto step = this->step_;
			return _generator(sde, iterations, numberSteps, delta, step, seed);
		}

		template<typename SDE>
		PathValuesType<PathValuesType<T>> simulateWithTimePoints(SDE const &sde, std::size_t iterations,
			TimePointsType<T> const &timePoints, std::random_device::result_type seed) const {
			auto step = this->step_;
			return _generatorWithTimePoints(sde, iterations, timePoints, step, seed);
		}
	};



	template<typename T>
	class MilsteinSchemeGPU<2, T> :public SchemeBuilderGPU<2, T> {
	private:
		T step_ = 10e-6;
	public:
		explicit MilsteinSchemeGPU(T delta, std::size_t numberSteps, T correlation)
			:SchemeBuilderGPU<2, T>{delta, numberSteps,correlation }
		{}

		explicit MilsteinSchemeGPU<2, T>(T correlation)
			:SchemeBuilderGPU<2, T>{correlation } {}
	private:
		template<typename SDEs>
		PathValuesType<PathValuesType<T>> _generator(SDEs const &sdes, std::size_t iterations, std::size_t numberSteps,
			T delta,T correlation,T step, std::random_device::result_type seed)const {

			// retrieve the sde components:
			auto drift1 = std::get<0>(std::get<0>(sdes));
			auto diffusion1 = std::get<1>(std::get<0>(sdes));
			auto init1 = std::get<2>(std::get<0>(sdes));

			auto drift2 = std::get<0>(std::get<1>(sdes));
			auto diffusion2 = std::get<1>(std::get<1>(sdes));
			auto init2 = std::get<2>(std::get<1>(sdes));

			auto randExtent = concurrency::extent<1>(iterations);
			sobol_rng_collection<sobol_rng<1>, 1> sbRand(randExtent, seed);
			auto simExtent = concurrency::extent<2>(iterations, numberSteps);
			concurrency::array_view<T, 2> sims(simExtent);

			NormalVariate<T> normalRand(0.0, 1.0);

			concurrency::parallel_for_each(randExtent, [=](concurrency::index<1> simIdx)restrict(amp) {
				auto rnd = sbRand[simIdx];
				T firstSpot{ init1 };
				T secondSpot{ init2 };
				T firstSpotNew{};
				T secondSpotNew{};
				T z1{};
				T z2{};
				rnd.skip(sbRand.direction_numbers(), simIdx[0]);
				sims(simIdx[0], 0) = init1;
				for (std::size_t t = 1; t < numberSteps; ++t) {
					z1 = normalRand(rnd.get_single(1), rnd.get_single(1));
					z2 = normalRand(rnd.get_single(1), rnd.get_single(1));

					firstSpotNew = firstSpot +
						drift1((t - 1)*delta, firstSpot, secondSpot)*delta +
						diffusion1((t - 1)*delta, firstSpot, secondSpot)*
						concurrency::fast_math::sqrt(delta)*z1 +
						0.5*diffusion1((t - 1)*delta, firstSpot, secondSpot) *
						((diffusion1((t - 1)*delta, firstSpot + 0.5*step, secondSpot) -
							diffusion1((t - 1)*delta, firstSpot - 0.5*step, secondSpot)) / (step))*
						delta*((z1)*(z1)-1.0) +
						0.5*correlation*diffusion2((t - 1)*delta, firstSpot, secondSpot) *
						((diffusion1((t - 1)*delta, firstSpot, secondSpot + 0.5*step) -
							diffusion1((t - 1)*delta, firstSpot, secondSpot - 0.5*step)) / (step)) *
						delta*((z1)*(z1)-1.0) +
						concurrency::fast_math::sqrt(1.0 - (correlation)*(correlation))*
						diffusion2((t - 1)*delta, firstSpot, secondSpot)*
						((diffusion1((t - 1)*delta, firstSpot, secondSpot + 0.5*step) -
							diffusion1((t - 1)*delta, firstSpot, secondSpot - 0.5*step)) / (step)) *
							delta*z1*z2;

					secondSpotNew = secondSpot +
						drift2((t - 1)*delta, firstSpot, secondSpot)*delta +
						diffusion2((t - 1)*delta, firstSpot, secondSpot) *
						concurrency::fast_math::sqrt(delta) *
						(correlation * z1 + concurrency::fast_math::sqrt(1.0 - (correlation*correlation)) * z2) +
						0.5*correlation*diffusion1((t - 1)*delta, firstSpot, secondSpot)*
						((diffusion2((t - 1)*delta, firstSpot + 0.5*step, secondSpot) -
							diffusion2((t - 1)*delta, firstSpot - 0.5*step, secondSpot)) / (step))*
							delta * ((z1)*(z1)-1.0) +
						0.5 * diffusion2((t - 1)*delta, firstSpot, secondSpot)*
						((diffusion2((t - 1)*delta, firstSpot, secondSpot + 0.5*step) -
							diffusion2((t - 1)*delta, firstSpot, secondSpot - 0.5*step)) / (step))*
							delta * ((correlation*z1 + concurrency::fast_math::sqrt(1.0 - (correlation*correlation))*z2)*
						(correlation*z1 + concurrency::fast_math::sqrt(1.0 - (correlation*correlation))*z2) - 1.0) +
						concurrency::fast_math::sqrt(1.0 - (correlation*correlation))*
						diffusion1((t - 1)*delta, firstSpot, secondSpot)*
						((diffusion2((t - 1)*delta, firstSpot + 0.5*step, secondSpot) -
							diffusion2((t - 1)*delta, firstSpot - 0.5*step, secondSpot)) / (step)) *
							delta*z1*z2;

					sims(simIdx[0], t) = firstSpotNew;
					firstSpot = firstSpotNew;
					secondSpot = secondSpotNew;
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

		template<typename SDEs>
		PathValuesType<PathValuesType<T>> _generatorWithTimePoints(SDEs const &sdes, std::size_t iterations,
			TimePointsType<T> const &timePoints,T correlation, T step, std::random_device::result_type seed)const {

			// retrieve the sde components:
			auto drift1 = std::get<0>(std::get<0>(sdes));
			auto diffusion1 = std::get<1>(std::get<0>(sdes));
			auto init1 = std::get<2>(std::get<0>(sdes));

			auto drift2 = std::get<0>(std::get<1>(sdes));
			auto diffusion2 = std::get<1>(std::get<1>(sdes));
			auto init2 = std::get<2>(std::get<1>(sdes));


			std::size_t const numberSteps = timePoints.size();
			auto randExtent = concurrency::extent<1>(iterations);
			sobol_rng_collection<sobol_rng<1>, 1> sbRand(randExtent, seed);
			auto simExtent = concurrency::extent<2>(iterations, numberSteps);
			concurrency::array_view<T, 2> sims(simExtent);
			concurrency::array_view<const T, 1> tmps(numberSteps, timePoints);

			NormalVariate<T> normalRand(0.0, 1.0);

			concurrency::parallel_for_each(randExtent, [=](concurrency::index<1> simIdx)restrict(amp) {
				auto rnd = sbRand[simIdx];
				T firstSpot{ init1 };
				T secondSpot{ init2 };
				T firstSpotNew{};
				T secondSpotNew{};
				T z1{};
				T z2{};
				rnd.skip(sbRand.direction_numbers(), simIdx[0]);
				sims(simIdx[0], 0) = init1;
				for (std::size_t t = 1; t < numberSteps; ++t) {
					firstSpotNew = firstSpot +
						drift1(tmps[t - 1], firstSpot, secondSpot)*(tmps[t]- tmps[t - 1]) +
						diffusion1(tmps[t - 1], firstSpot, secondSpot)*
						concurrency::fast_math::sqrt(tmps[t]- tmps[t - 1])*z1 +
						0.5*diffusion1(tmps[t - 1], firstSpot, secondSpot) *
						((diffusion1(tmps[t - 1], firstSpot + 0.5*step, secondSpot) -
							diffusion1(tmps[t - 1], firstSpot - 0.5*step, secondSpot)) / (step))*
						(tmps[t]- tmps[t - 1])*((z1)*(z1)-1.0) +
						0.5*correlation*diffusion2(tmps[t - 1], firstSpot, secondSpot) *
						((diffusion1(tmps[t - 1], firstSpot, secondSpot + 0.5*step) -
							diffusion1(tmps[t - 1], firstSpot, secondSpot - 0.5*step)) / (step)) *
						(tmps[t]- tmps[t - 1])*((z1)*(z1)-1.0) +
						concurrency::fast_math::sqrt(1.0 - (correlation)*(correlation))*
						diffusion2(tmps[t - 1], firstSpot, secondSpot)*
						((diffusion1(tmps[t - 1], firstSpot, secondSpot + 0.5*step) -
							diffusion1(tmps[t - 1], firstSpot, secondSpot - 0.5*step)) / (step)) *
						(tmps[t]- tmps[t - 1])*z1*z2;

					secondSpotNew = secondSpot +
						drift2(tmps[t - 1], firstSpot, secondSpot)*(tmps[t]- tmps[t - 1]) +
						diffusion2(tmps[t - 1], firstSpot, secondSpot) *
						concurrency::fast_math::sqrt(tmps[t]- tmps[t - 1]) *
						(correlation * z1 + concurrency::fast_math::sqrt(1.0 - (correlation*correlation)) * z2) +
						0.5*correlation*diffusion1(tmps[t - 1], firstSpot, secondSpot)*
						((diffusion2(tmps[t - 1], firstSpot + 0.5*step, secondSpot) -
							diffusion2(tmps[t - 1], firstSpot - 0.5*step, secondSpot)) / (step))*
						(tmps[t]- tmps[t - 1]) * ((z1)*(z1)-1.0) +
						0.5 * diffusion2(tmps[t - 1], firstSpot, secondSpot)*
						((diffusion2(tmps[t - 1], firstSpot, secondSpot + 0.5*step) -
							diffusion2(tmps[t - 1], firstSpot, secondSpot - 0.5*step)) / (step))*
						(tmps[t]- tmps[t - 1]) * ((correlation*z1 + concurrency::fast_math::sqrt(1.0 - (correlation*correlation))*z2)*
						(correlation*z1 + concurrency::fast_math::sqrt(1.0 - (correlation*correlation))*z2) - 1.0) +
						concurrency::fast_math::sqrt(1.0 - (correlation*correlation))*
						diffusion1(tmps[t - 1], firstSpot, secondSpot)*
						((diffusion2(tmps[t - 1], firstSpot + 0.5*step, secondSpot) -
							diffusion2(tmps[t - 1], firstSpot - 0.5*step, secondSpot)) / (step)) *
						(tmps[t]- tmps[t - 1])*z1*z2;

					sims(simIdx[0], t) = firstSpotNew;
					firstSpot = firstSpotNew;
					secondSpot = secondSpotNew;
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
		template<typename SDEs>
		PathValuesType<PathValuesType<T>> simulate(SDEs const &sdes, std::size_t iterations,
			std::random_device::result_type seed) const {
			auto numberSteps = this->numberSteps_;
			auto delta = this->delta_;
			auto correlation = this->correlation_;
			auto step = this->step_;
			return _generator(sdes, iterations, numberSteps, delta, correlation, step, seed);
		}

		template<typename SDEs>
		PathValuesType<PathValuesType<T>> simulateWithTimePoints(SDEs const &sdes, std::size_t iterations,
			TimePointsType<T> const &timePoints, std::random_device::result_type seed) const {
			auto correlation = this->correlation_;
			auto step = this->step_;
			return _generatorWithTimePoints(sdes, iterations, timePoints, correlation,step, seed);
		}
	};


}



#endif ///_FDM_SCHEME_GPU_H_