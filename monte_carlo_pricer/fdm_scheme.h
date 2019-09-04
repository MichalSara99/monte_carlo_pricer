#pragma once
#if !defined(_FDM_SCHEME_H_)
#define _FDM_SCHEME_H_

#include"mc_types.h"

namespace finite_difference_method {


	template<typename T = double,
			typename =typename std::enable_if<std::is_arithmetic<T>::value>::type>
	class Scheme {
	public:
		virtual ~Scheme(){}



	};

	template<typename T = double,
		typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	class EulerScheme:public Scheme<T> {

	};

	template<typename T = double,
		typename = typename std::enable_if<std::is_arithmetic<T>::value>::type>
	class MilsteinScheme:public Scheme<T> {


	};
}



#endif ///_FDM_SCHEME_H_