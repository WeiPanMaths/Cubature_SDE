#pragma once
#ifndef BASKETSPREADQUANTLIB_H__
#define BASKETSPREADQUANTLIB_H__

#include <ql/pricingengines/basket/mceuropeanbasketengine.hpp>
#ifdef BOOST_MSVC
#  include <ql/auto_link.hpp>
#endif

using namespace QuantLib;

namespace MyBasket {
    enum BasketType { SpreadBasket };
	
    struct BasketOptionTwoData {
        BasketType basketType;
        Option::Type type;
        Real strike;
        Real s1;
        Real s2;
        Rate q1;
        Rate q2;
        Rate r;
        Time t; // years
        Volatility v1;
        Volatility v2;
        Real rho;
        Real result;
        Real tol;
    };

    class MyBasketSolver {
        BasketOptionTwoData _data;
        ext::shared_ptr<PricingEngine> analyticEngine;
        ext::shared_ptr<PricingEngine> mcEngine;
        ext::shared_ptr<PricingEngine> fdEngine;
        ext::shared_ptr<BasketOption> basketOption;

    public:
        MyBasketSolver(BasketOptionTwoData data);
        Real CalculationMC() const;
        Real CalculationFD() const;
        Real CalculationKirk() const;

    };
} // namespace MyBasket

#endif // !_BASKETSPREADQUANTLIB_H__

