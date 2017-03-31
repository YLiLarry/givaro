/* -*- mode: C++; tab-width: 4; indent-tabs-mode: t; c-basic-offset: 4 -*- */
// vim:sts=4:sw=4:ts=4:noet:sr:cino=>s,f0,{0,g0,(0,\:0,t0,+0,=s
// ==========================================================================
// Copyright(c)'1994-2017 by The Givaro group
// This file is part of Givaro.
// Givaro is governed by the CeCILL-B license under French law
// and abiding by the rules of distribution of free software.
// see the COPYRIGHT file for more details.
// Authors: B. Grenet, R. Lebreton from existing files
// ==========================================================================


#ifndef __GIVARO_modular_implem_INL
#define __GIVARO_modular_implem_INL

#include <cmath>

#include "modular-defines.h"

#define SELF_T Modular_implem<Storage_t, Compute_t, Residu_t>

namespace Givaro {

	// ------------------------
	// ----- Classic arithmetic

	template<typename Storage_t, typename Compute_t, typename Residu_t>
	inline typename SELF_T::Element&
	SELF_T::mul
	(Element& r, const Element& a, const Element& b) const
	{
		return  __GIVARO_MODULAR_INTEGER_MUL(r,_p,a,b);
	}

	template<typename Storage_t, typename Compute_t, typename Residu_t>
	inline typename SELF_T::Element&
	SELF_T::sub
	(Element& r, const Element& a, const Element& b) const
	{
		return __GIVARO_MODULAR_INTEGER_SUB(r,_p,a,b);
	}

	template<typename Storage_t, typename Compute_t, typename Residu_t>
	inline typename SELF_T::Element&
	SELF_T::add
	(Element& r, const Element& a, const Element& b) const
	{
		__GIVARO_MODULAR_INTEGER_ADD(r,_p,a,b);
		return r;
	}

	template<typename Storage_t, typename Compute_t, typename Residu_t>
	inline typename SELF_T::Element&
	SELF_T::neg
	(Element& r, const Element& a) const
	{
		return __GIVARO_MODULAR_INTEGER_NEG(r,_p,a);
	}

	template<typename Storage_t, typename Compute_t, typename Residu_t>
	inline typename SELF_T::Element&
	SELF_T::inv
	(Element& r, const Element& a) const
	{
		invext(r, a, typename SELF_T::Element(_p));
		return (r < 0)? r += typename SELF_T::Element(_p) : r;
	}

	template<typename Storage_t, typename Compute_t, typename Residu_t>
	inline typename SELF_T::Element&
	SELF_T::div
	(Element& r, const Element& a, const Element& b) const
	{
		return mulin(inv(r, b), a);
	}

	template<typename Storage_t, typename Compute_t, typename Residu_t>
	inline typename SELF_T::Element&
	SELF_T::mulin
	(Element& r, const Element& a) const
	{
		return __GIVARO_MODULAR_INTEGER_MULIN(r, _p, a);
	}

	template<typename Storage_t, typename Compute_t, typename Residu_t>
	inline typename SELF_T::Element&
	SELF_T::divin
	(Element& r, const Element& a) const
	{
		Element ia;
		return mulin(r, inv(ia, a));
	}

	template<typename Storage_t, typename Compute_t, typename Residu_t>
	inline typename SELF_T::Element&
	SELF_T::addin
	(Element& r, const Element& a) const
	{
		typename SELF_T::Element tmp = r;
		__GIVARO_MODULAR_INTEGER_ADDIN(tmp,_p, a);
		return r = Element(tmp);
	}

	template<typename Storage_t, typename Compute_t, typename Residu_t>
	inline typename SELF_T::Element&
	SELF_T::subin
	(Element& r, const Element& a) const
	{
		typename SELF_T::Element tmp = r;
		__GIVARO_MODULAR_INTEGER_SUBIN(tmp,_p, a);
		return r = (typename SELF_T::Element)tmp;
	}

	template<typename Storage_t, typename Compute_t, typename Residu_t>
	inline typename SELF_T::Element&
	SELF_T::negin
	(Element& r) const
	{
		return __GIVARO_MODULAR_INTEGER_NEGIN(r,_p);
	}

	template<typename Storage_t, typename Compute_t, typename Residu_t>
	inline typename SELF_T::Element&
	SELF_T::invin
	(Element& r) const
	{
		return inv(r, r);
	}

//	template<typename Storage_t, typename Compute_t, typename Residu_t>
//	inline typename SELF_T::Element&
//	SELF_T::axpy
//	(Element& r, const Element& a, const Element& b, const Element& c) const
//	{
//		return __GIVARO_MODULAR_INTEGER_MULADD(r, _p, a, b, c);
//	}

//	template<typename Storage_t, typename Compute_t, typename Residu_t>
//	inline typename SELF_T::Element&
//	SELF_T::axpyin
//	(Element& r, const Element& a, const Element& b) const
//	{
//		return __GIVARO_MODULAR_INTEGER_MULADDIN(r, _p, a, b);
//	}

//	template<typename Storage_t, typename Compute_t, typename Residu_t>
//	inline typename SELF_T::Element&
//	SELF_T::maxpy
//	(Element& r, const Element& a, const Element& b, const Element& c) const
//	{
//		int32_t tmp;
//		__GIVARO_MODULAR_INTEGER_MUL(tmp, _p, a, b);
//		__GIVARO_MODULAR_INTEGER_SUB(r, _p, c, tmp);
//		return r;
//	}

//	template<typename Storage_t, typename Compute_t, typename Residu_t>
//	inline typename SELF_T::Element&
//	SELF_T::axmy
//	(Element& r, const Element& a, const Element& b, const Element& c) const
//	{
//		int32_t tmp;
//		__GIVARO_MODULAR_INTEGER_MULSUB(tmp, _p, a, b, c);
//		return r = (Modular<int32_t, COMP>::Element)tmp;
//	}

//	template<typename Storage_t, typename Compute_t, typename Residu_t>
//	inline typename SELF_T::Element&
//	SELF_T::maxpyin
//	(Element& r, const Element& a, const Element& b) const
//	{
//		__GIVARO_MODULAR_INTEGER_SUBMULIN(r, _p, a, b );
//		return r;
//	}

//	template<typename Storage_t, typename Compute_t, typename Residu_t>
//	inline typename SELF_T::Element&
//	SELF_T::axmyin
//	(Element& r, const Element& a, const Element& b) const
//	{
//		return __GIVARO_MODULAR_INTEGER_MULSUB(r, _p, a, b, r );
//	}

//	template<typename Storage_t, typename Compute_t, typename Residu_t>
//	inline bool
//	SELF_T::isUnit(const Element& a) const
//	{
//		Element u,d;
//		invext(u,d,a,int32_t(_p));
//		return isOne(d) || isMOne(d);
//	}

}

#endif // __GIVARO_modular_implem_INL
