#include <iostream>

#include "BasketSpreadQuantLib.h"

#include <ql/quotes/simplequote.hpp>
#include <ql/time/daycounters/actual360.hpp>
#include <ql/instruments/basketoption.hpp>

//#include <ql/pricingengines/basket/stulzengine.hpp>

#include <ql/pricingengines/basket/kirkengine.hpp>
#include <ql/pricingengines/basket/mceuropeanbasketengine.hpp>

//#include <ql/pricingengines/basket/mcamericanbasketengine.hpp>

#include <ql/processes/blackscholesprocess.hpp>
#include <ql/processes/stochasticprocessarray.hpp>
#include <ql/termstructures/yield/flatforward.hpp>
#include <ql/termstructures/volatility/equityfx/blackconstantvol.hpp>

//#include <ql/termstructures/volatility/equityfx/hestonblackvolsurface.hpp>

#include <ql/pricingengines/basket/fd2dblackscholesvanillaengine.hpp>
//#include <ql/utilities/dataformatters.hpp>
//#include <ql/functional.hpp>
#include <ql/time/calendars/nullcalendar.hpp>

#ifdef BOOST_MSVC
#  include <ql/auto_link.hpp>
#endif

using namespace QuantLib;

//#define LENGTH(a) (sizeof(a)/sizeof(a[0]))

namespace MyBasket{

    //enum BasketType { MinBasket, MaxBasket, SpreadBasket };
    //enum BasketType { SpreadBasket };
    
    Real relativeError(Real x1, Real x2, Real reference) {
        if (reference != 0.0)
            return std::fabs(x1 - x2) / reference;
        else
            // fall back to absolute error
            return std::fabs(x1 - x2);
    }

    ext::shared_ptr<BasketPayoff> basketTypeToPayoff(
        BasketType basketType,
        const ext::shared_ptr<Payoff>& p) {
        switch (basketType) {
        //case MinBasket:
        //    return ext::shared_ptr<BasketPayoff>(new MinBasketPayoff(p));
        //case MaxBasket:
        //    return ext::shared_ptr<BasketPayoff>(new MaxBasketPayoff(p));
        case SpreadBasket:
            return ext::shared_ptr<BasketPayoff>(new SpreadBasketPayoff(p));
        }
        QL_FAIL("unknown basket option type");
    }

    //struct BasketOptionTwoData {
    //    BasketType basketType;
    //    Option::Type type;
    //    Real strike;
    //    Real s1;
    //    Real s2;
    //    Rate q1;
    //    Rate q2;
    //    Rate r;
    //    Time t; // years
    //    Volatility v1;
    //    Volatility v2;
    //    Real rho;
    //    Real result;
    //    Real tol;
    //};

    ext::shared_ptr<YieldTermStructure>
        flatRate(const Date& today,
            const ext::shared_ptr<Quote>& forward,
            const DayCounter& dc) {
        return ext::shared_ptr<YieldTermStructure>(
            new FlatForward(today, Handle<Quote>(forward), dc));
    }

    ext::shared_ptr<YieldTermStructure>
        flatRate(const Date& today, Rate forward, const DayCounter& dc) {
        return flatRate(
            today, ext::shared_ptr<Quote>(new SimpleQuote(forward)), dc);
    }

    ext::shared_ptr<YieldTermStructure>
        flatRate(const ext::shared_ptr<Quote>& forward,
            const DayCounter& dc) {
        return ext::shared_ptr<YieldTermStructure>(
            new FlatForward(0, NullCalendar(), Handle<Quote>(forward), dc));
    }

    ext::shared_ptr<YieldTermStructure>
        flatRate(Rate forward, const DayCounter& dc) {
        return flatRate(ext::shared_ptr<Quote>(new SimpleQuote(forward)),
            dc);
    }


    ext::shared_ptr<BlackVolTermStructure>
        flatVol(const Date& today,
            const ext::shared_ptr<Quote>& vol,
            const DayCounter& dc) {
        return ext::shared_ptr<BlackVolTermStructure>(new
            BlackConstantVol(today, NullCalendar(), Handle<Quote>(vol), dc));
    }

    ext::shared_ptr<BlackVolTermStructure>
        flatVol(const Date& today, Volatility vol,
            const DayCounter& dc) {
        return flatVol(today,
            ext::shared_ptr<Quote>(new SimpleQuote(vol)),
            dc);
    }

    ext::shared_ptr<BlackVolTermStructure>
        flatVol(const ext::shared_ptr<Quote>& vol,
            const DayCounter& dc) {
        return ext::shared_ptr<BlackVolTermStructure>(new
            BlackConstantVol(0, NullCalendar(), Handle<Quote>(vol), dc));
    }

    ext::shared_ptr<BlackVolTermStructure>
        flatVol(Volatility vol,
            const DayCounter& dc) {
        return flatVol(ext::shared_ptr<Quote>(new SimpleQuote(vol)), dc);
    }

    //class MyBasketSolver {
    //    BasketOptionTwoData _data;
    //    ext::shared_ptr<PricingEngine> analyticEngine;
    //    ext::shared_ptr<PricingEngine> mcEngine;
    //    ext::shared_ptr<PricingEngine> fdEngine;
    //    ext::shared_ptr<BasketOption> basketOption;

    //public:
    MyBasketSolver::MyBasketSolver(BasketOptionTwoData data) : _data(data) 
    {
        DayCounter dc = Actual360();
        Date today = Date::todaysDate();

        ext::shared_ptr<SimpleQuote> spot1(new SimpleQuote(0.0));
        ext::shared_ptr<SimpleQuote> spot2(new SimpleQuote(0.0));

        ext::shared_ptr<SimpleQuote> qRate1(new SimpleQuote(0.0));
        ext::shared_ptr<SimpleQuote> qRate2(new SimpleQuote(0.0));

        ext::shared_ptr<SimpleQuote> rRate(new SimpleQuote(0.0));
        ext::shared_ptr<YieldTermStructure> rTS = MyBasket::flatRate(today, rRate, dc);

        ext::shared_ptr<SimpleQuote> vol1(new SimpleQuote(0.0));
        ext::shared_ptr<BlackVolTermStructure> volTS1 = MyBasket::flatVol(today, vol1, dc);
        ext::shared_ptr<SimpleQuote> vol2(new SimpleQuote(0.0));
        ext::shared_ptr<BlackVolTermStructure> volTS2 = MyBasket::flatVol(today, vol2, dc);

        ext::shared_ptr<PlainVanillaPayoff> payoff(new
            PlainVanillaPayoff(_data.type, _data.strike));

        Date exDate = today + Integer(_data.t * 360 + 0.5);
        ext::shared_ptr<Exercise> exercise(new EuropeanExercise(exDate));

        spot1->setValue(_data.s1);
        spot2->setValue(_data.s2);
        qRate1->setValue(_data.q1);
        qRate2->setValue(_data.q2);
        rRate->setValue(_data.r);
        vol1->setValue(_data.v1);
        vol2->setValue(_data.v2);

        ext::shared_ptr<GeneralizedBlackScholesProcess> p1, p2;
        switch (_data.basketType) {
        case MyBasket::SpreadBasket:
            p1 = ext::shared_ptr<GeneralizedBlackScholesProcess>(
                new BlackProcess(Handle<Quote>(spot1),
                    Handle<YieldTermStructure>(rTS),
                    Handle<BlackVolTermStructure>(volTS1)));
            p2 = ext::shared_ptr<GeneralizedBlackScholesProcess>(
                new BlackProcess(Handle<Quote>(spot2),
                    Handle<YieldTermStructure>(rTS),
                    Handle<BlackVolTermStructure>(volTS2)));

            analyticEngine = ext::shared_ptr<PricingEngine>(
                new KirkEngine(ext::dynamic_pointer_cast<BlackProcess>(p1),
                    ext::dynamic_pointer_cast<BlackProcess>(p2),
                    _data.rho));
            break;
        default:
            QL_FAIL("unknown basket type");
        }

        std::vector<ext::shared_ptr<StochasticProcess1D> > procs;
        procs.push_back(p1);
        procs.push_back(p2);

        Matrix correlationMatrix(2, 2, _data.rho);
        for (Integer j = 0; j < 2; j++) {
            correlationMatrix[j][j] = 1.0;
        }

        ext::shared_ptr<StochasticProcessArray> process(
            new StochasticProcessArray(procs, correlationMatrix));

        mcEngine =
            MakeMCEuropeanBasketEngine<PseudoRandom, Statistics>(process)
            .withStepsPerYear(1)
            .withSamples(100000)  // 10000
            .withSeed(42);

        fdEngine = ext::shared_ptr<PricingEngine>(
            new Fd2dBlackScholesVanillaEngine(p1, p2, _data.rho,
                50, 50, 15));  // 50 50 15

        basketOption = ext::shared_ptr<BasketOption>(new BasketOption(basketTypeToPayoff(_data.basketType,
            payoff),
            exercise));

    }

    Real MyBasketSolver::CalculationMC() const
    {
        const Real mcRelativeErrorTolerance = 0.01;
        // mc engine
        basketOption->setPricingEngine(mcEngine);
        Real calculated = basketOption->NPV();
        Real expected = _data.result;
        Real relError = MyBasket::relativeError(calculated, expected, expected);
        //if (relError > mcRelativeErrorTolerance) {
        //    std::cout << "didn't make (mc) tolerance " << _data.tol << std::endl;
        //}
        //std::cout << relError << "\n";
        return calculated;
    }

    Real MyBasketSolver::CalculationFD() const
    {
        const Real fdRelativeErrorTolerance = 0.01;
        // fd engine
        basketOption->setPricingEngine(fdEngine);
        Real calculated = basketOption->NPV();
        Real expected = _data.result;
        Real relError = MyBasket::relativeError(calculated, expected, expected);
        //if (relError > fdRelativeErrorTolerance) {
        //    std::cout << "didn't make (fd) tolerance " << _data.tol << std::endl;
        //}
        //std::cout << relError << ", MC: ";
        return calculated;
    }

    Real MyBasketSolver::CalculationKirk() const
    {            
        // analytic engine
        basketOption->setPricingEngine(analyticEngine);
        Real calculated = basketOption->NPV();
        Real expected = _data.result;
        Real error = std::fabs(calculated - expected);
        //if (error > _data.tol) {
        //    std::cout << "didn't make (analystic) tolerance " << _data.tol << std::endl;
        //}
        //std::cout << "analytic: " << error << ", FD: ";
        return calculated;
    }

    //};
} // namespace MyBasket 

# if 0
void Test()
{
    /*
        Data from:
        Excel spreadsheet www.maths.ox.ac.uk/~firth/computing/excel.shtml
        and
        "Option pricing formulas", E.G. Haug, McGraw-Hill 1998 pag 56-58
        European two asset max basket options
    */
    MyBasket::BasketOptionTwoData values[] = {
        //      basketType,   optionType, strike,    s1,    s2,   q1,   q2,    r,    t,   v1,   v2,  rho, result, tol
        // data from http://www.maths.ox.ac.uk/~firth/computing/excel.shtml

        //      basketType,   optionType, strike,    s1,    s2,   q1,   q2,    r,    t,   v1,   v2,  rho,  result, tol
        // data from "Option pricing formulas" VB code + spreadsheet

        /* "Option pricing formulas", E.G. Haug, McGraw-Hill 1998 pag 59-60
            Kirk approx. for a european spread option on two futures*/
#if 1
        {MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.1, 0.20, 0.20, -0.5, 4.7530, 1.0e-3},
        {MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.1, 0.20, 0.20,  0.0, 3.7970, 1.0e-3},
        {MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.1, 0.20, 0.20,  0.5, 2.5537, 1.0e-3},
        {MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.1, 0.25, 0.20, -0.5, 5.4275, 1.0e-3},
        {MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.1, 0.25, 0.20,  0.0, 4.3712, 1.0e-3},
        {MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.1, 0.25, 0.20,  0.5, 3.0086, 1.0e-3},
        {MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.1, 0.20, 0.25, -0.5, 5.4061, 1.0e-3},
        {MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.1, 0.20, 0.25,  0.0, 4.3451, 1.0e-3},
        {MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.1, 0.20, 0.25,  0.5, 2.9723, 1.0e-3},
        {MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.5, 0.20, 0.20, -0.5,10.7517, 1.0e-3},
        {MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.5, 0.20, 0.20,  0.0, 8.7020, 1.0e-3},
        {MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.5, 0.20, 0.20,  0.5, 6.0257, 1.0e-3},
        {MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.5, 0.25, 0.20, -0.5,12.1941, 1.0e-3},
        {MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.5, 0.25, 0.20,  0.0, 9.9340, 1.0e-3},
        {MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.5, 0.25, 0.20,  0.5, 7.0067, 1.0e-3},
        {MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.5, 0.20, 0.25, -0.5,12.1483, 1.0e-3},
        {MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.5, 0.20, 0.25,  0.0, 9.8780, 1.0e-3},
#endif
        {MyBasket::SpreadBasket, Option::Call, 3.0,  122.0, 120.0, 0.0, 0.0, 0.10,  0.5, 0.20, 0.25,  0.5, 6.9284, 1.0e-3}
    };

    for (Size i = 0; i < LENGTH(values); i++) 
    {
        MyBasket::MyBasketSolver solver(values[i]);
        solver.CalculationKirk();
        solver.CalculationFD();
        solver.CalculationMC();
    }
#if 0
    DayCounter dc = Actual360();
    Date today = Date::todaysDate();

    ext::shared_ptr<SimpleQuote> spot1(new SimpleQuote(0.0));
    ext::shared_ptr<SimpleQuote> spot2(new SimpleQuote(0.0));

    ext::shared_ptr<SimpleQuote> qRate1(new SimpleQuote(0.0));
    //ext::shared_ptr<YieldTermStructure> qTS1 = MyBasket::flatRate(today, qRate1, dc);
    ext::shared_ptr<SimpleQuote> qRate2(new SimpleQuote(0.0));
    //ext::shared_ptr<YieldTermStructure> qTS2 = MyBasket::flatRate(today, qRate2, dc);

    ext::shared_ptr<SimpleQuote> rRate(new SimpleQuote(0.0));
    ext::shared_ptr<YieldTermStructure> rTS = MyBasket::flatRate(today, rRate, dc);

    ext::shared_ptr<SimpleQuote> vol1(new SimpleQuote(0.0));
    ext::shared_ptr<BlackVolTermStructure> volTS1 = MyBasket::flatVol(today, vol1, dc);
    ext::shared_ptr<SimpleQuote> vol2(new SimpleQuote(0.0));
    ext::shared_ptr<BlackVolTermStructure> volTS2 = MyBasket::flatVol(today, vol2, dc);

    const Real mcRelativeErrorTolerance = 0.01;
    const Real fdRelativeErrorTolerance = 0.01;

    for (Size i = 0; i < LENGTH(values); i++) {

        ext::shared_ptr<PlainVanillaPayoff> payoff(new
            PlainVanillaPayoff(values[i].type, values[i].strike));

        Date exDate = today + Integer(values[i].t * 360 + 0.5);
        ext::shared_ptr<Exercise> exercise(new EuropeanExercise(exDate));

        spot1->setValue(values[i].s1);
        spot2->setValue(values[i].s2);
        qRate1->setValue(values[i].q1);
        qRate2->setValue(values[i].q2);
        rRate->setValue(values[i].r);
        vol1->setValue(values[i].v1);
        vol2->setValue(values[i].v2);

        ext::shared_ptr<PricingEngine> analyticEngine;
        ext::shared_ptr<GeneralizedBlackScholesProcess> p1, p2;
        switch (values[i].basketType) {
        case MyBasket::SpreadBasket:
            p1 = ext::shared_ptr<GeneralizedBlackScholesProcess>(
                new BlackProcess(Handle<Quote>(spot1),
                    Handle<YieldTermStructure>(rTS),
                    Handle<BlackVolTermStructure>(volTS1)));
            p2 = ext::shared_ptr<GeneralizedBlackScholesProcess>(
                new BlackProcess(Handle<Quote>(spot2),
                    Handle<YieldTermStructure>(rTS),
                    Handle<BlackVolTermStructure>(volTS2)));

            analyticEngine = ext::shared_ptr<PricingEngine>(
                new KirkEngine(ext::dynamic_pointer_cast<BlackProcess>(p1),
                    ext::dynamic_pointer_cast<BlackProcess>(p2),
                    values[i].rho));
            break;
        default:
            QL_FAIL("unknown basket type");
        }

        std::vector<ext::shared_ptr<StochasticProcess1D> > procs;
        procs.push_back(p1);
        procs.push_back(p2);

        Matrix correlationMatrix(2, 2, values[i].rho);
        for (Integer j = 0; j < 2; j++) {
            correlationMatrix[j][j] = 1.0;
        }

        ext::shared_ptr<StochasticProcessArray> process(
            new StochasticProcessArray(procs, correlationMatrix));

        ext::shared_ptr<PricingEngine> mcEngine =
            MakeMCEuropeanBasketEngine<PseudoRandom, Statistics>(process)
            .withStepsPerYear(1)
            .withSamples(10000)  // 10000
            .withSeed(42);

        ext::shared_ptr<PricingEngine> fdEngine(
            new Fd2dBlackScholesVanillaEngine(p1, p2, values[i].rho,
                50, 50, 15));  // 50 50 15

        BasketOption basketOption(basketTypeToPayoff(values[i].basketType,
            payoff),
            exercise);

        // analytic engine
        basketOption.setPricingEngine(analyticEngine);
        Real calculated = basketOption.NPV();
        Real expected = values[i].result;
        Real error = std::fabs(calculated - expected);

        std::cout << "analytic: " << error << ", FD: ";

        // fd engine
        basketOption.setPricingEngine(fdEngine);
        calculated = basketOption.NPV();
        Real relError = MyBasket::relativeError(calculated, expected, expected);
        std::cout << relError << ", MC: ";

        // mc engine
        basketOption.setPricingEngine(mcEngine);
        calculated = basketOption.NPV();
        relError = MyBasket::relativeError(calculated, expected, expected);
        std::cout << relError << "\n";
    }
#endif
}

#endif
//
//int main()
//{
//    Test();
//
//	return 0;
//}

