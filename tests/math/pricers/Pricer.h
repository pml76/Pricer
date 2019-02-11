/*
 *
 * (c) 2019, by Peter Lennartz  // peter.lennartz@gmail.com
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */


#pragma once

#include <arb.h>
#include <arb_hypgeom.h>

#include <iostream>
#include <sstream>
#include <string>
#include <exception>
#include <vector>

#include <cstdlib>
#include <execinfo.h>
#include <cmath>

#define THROW_EXCEPTION(e, t) ({std::stringstream ss;void *buffer[1000];int size = backtrace(buffer, 1000); char **bt = backtrace_symbols(buffer, size);for( int i = 0; i < size; i++ ) { ss << bt[i] << std::endl;} free(bt); ss << std::endl << "********************************" << std::endl << " E X C E P T I O N " << std::endl << "********************************" << std::endl << t << std::endl << "********************************" << std::endl; throw e(ss.str().data()); });


class not_normal_number_exception : public std::exception {
public:
    not_normal_number_exception(const char *s) {
        str = s;
    }

    virtual const char *what() const noexcept override {
        return str.data();
    }

private:

    std::string str;
};

class unprecise_calculation_exception : public std::exception {
public:
    unprecise_calculation_exception(const char *s) : str(s) {}

    virtual const char *what() const noexcept override {
        return str.data();
    }

private:

    std::string str;
};


/**
 * operator <<
 *
 * For debugging we it is good to be able to print arb-numbers.
 */
inline
std::ostream &operator<<(std::ostream &o, const arb_t &t) {
    o << "[" << arf_get_d(arb_midref(t), ARF_RND_NEAR) << "+-" << mag_get_d(arb_radref(t)) << "]";
    return o;
}


/**
 * Pricer routines.
 *
 * The class pricer contains all the routines that are used to price options.
 * All the computations are done using arb-routines. These routines offer us
 * the possibility to do computations at a given level of precision. While
 * this is not needed for every set of input variables, there is occasionally
 * a case where we have to increase the woring-precision to reach a given
 * level of precision in the target numbers.
 */
class Pricer {

private:

    mutable arb_t sigma,    /**< the vola */
            t,        /**< time to maturity in years */
            tau,    /**< time of the averaging-period in years*/
            r,        /**< interrest-rate */
            s,        /**< price of underlying */
            x;        /**< strike */

    mutable arb_t tmp1, tmp2,
            tmp3, tmp4,
            tmp5;

    /**
     * storage of intermediata results
     */
    mutable arb_t d1, d2,
            sigmaA,
            M,
            emrt,
            sigmaA2T2,
            sigmaAsqrtT;

    mutable arb_t price; /**< final price of the option */

    mutable arb_t msqrt2;

    /**
     * for precision-mgmt.
     */
    mutable slong prec, default_prec, target_prec;

    // compute_d_values
    //
    // factored out code. used by each of the two price computation for either
    // put or call.
    void compute_d_values() const {

        const unsigned int additional_prec = 17;
        unsigned int iter = 0;
        while (iter++ < 100) {

            if (arf_is_pos_inf(arb_midref(x)) || arf_is_neg_inf(arb_midref(x))) {
                arb_neg_inf(d1);
                arb_neg_inf(d2);
                break;
            } else {
                arb_div(tmp2, s, x, prec);

                arb_log(tmp2, tmp2, prec); // tmp2 = ln(s/x);

                arb_add(tmp1, tmp2, sigmaA2T2, prec); // tmp1 = ln(s/x) + sigmaA^2*T/2

                arb_div(d1, tmp1, sigmaAsqrtT, prec);
                arb_sub(d2, d1, sigmaAsqrtT, prec);

                arb_div(d1, d1, msqrt2, prec); // divide both by -sqrt(2) since we only
                arb_div(d2, d2, msqrt2, prec); // have the erfc function in stock ...

                if (arb_rel_accuracy_bits(d1) < target_prec + additional_prec ||
                        arb_rel_accuracy_bits(d2) < target_prec + additional_prec) {
                    prec *= 2;
                    init();
                    continue;
                }

                break;
            }

            if (iter == 100) THROW_EXCEPTION(unprecise_calculation_exception,
                                             "Cannot do computation to required precision");
        }
    }

public:
    Pricer(slong target_prec = 50, slong default_prec = 100)
            : prec(default_prec), default_prec(default_prec), target_prec(target_prec) {

        arb_init(sigma); // the sigmatility
        arb_init(t);     // the time until maturity
        arb_init(tau);   // the time until the begin of the averaging period
        arb_init(r);     // the interest rate
        arb_init(s);     // the future price of the stock
        arb_init(x);     // the strike of the option

        // only allocate the temporal variables, but do not initialize them.
        arb_init(tmp1);
        arb_init(tmp2);
        arb_init(tmp3);
        arb_init(tmp4);
        arb_init(tmp5);

        // thiese will hold temporal data, therefore we do not
        // initialize them here
        arb_init(d1);
        arb_init(d2);
        arb_init(sigmaA);
        arb_init(M);
        arb_init(emrt);
        arb_init(sigmaA2T2);
        arb_init(sigmaAsqrtT);

        arb_init(price);

        arb_init(msqrt2);

    }

    void set_underlying(const double &s_) {
        arb_set_d(s, s_);

        if (!arf_is_normal(arb_midref(s))) {
            THROW_EXCEPTION(not_normal_number_exception,
                            "The s variable has to be a normal number, not equal to zero.");
        }

    }


    void set_market_data(const double &sigma_, const double &t_,
                         const double &tau_, const double &r_, const double &s_) {


        reset();

        arb_set_d(sigma, sigma_);
        arb_set_d(t, t_);
        arb_set_d(tau, tau_);
        arb_set_d(r, r_);
        arb_set_d(s, s_);

        if (!arf_is_normal(arb_midref(sigma))) {
            THROW_EXCEPTION(not_normal_number_exception,
                            "The sigma variable has to be a normal number, not equal to zero.");
        }

        if (!arf_is_normal(arb_midref(t))) {
            THROW_EXCEPTION(not_normal_number_exception,
                            "The t variable has to be a normal number, not equal to zero.");
        }

        if (!arf_is_normal(arb_midref(tau))) {
            THROW_EXCEPTION(not_normal_number_exception,
                            "The tau variable has to be a normal number, not equal to zero.");
        }

        if (arf_is_inf(arb_midref(r)) || arf_is_nan(arb_midref(r))) {
            THROW_EXCEPTION(not_normal_number_exception,
                            "The r variable has to be a normal number, not equal to zero.");
        }

        if (!arf_is_normal(arb_midref(s))) {
            THROW_EXCEPTION(not_normal_number_exception,
                            "The s variable has to be a normal number, not equal to zero.");
        }

        init();
    }

    ~Pricer() {
        arb_clear(x);
        arb_clear(s);
        arb_clear(r);
        arb_clear(tau);
        arb_clear(t);
        arb_clear(sigma);

        arb_clear(tmp1);
        arb_clear(tmp2);
        arb_clear(tmp3);
        arb_clear(tmp4);
        arb_clear(tmp5);

        arb_clear(d1);
        arb_clear(d2);
        arb_clear(sigmaA);
        arb_clear(M);
        arb_clear(emrt);
        arb_clear(sigmaA2T2);
        arb_clear(sigmaAsqrtT);

        arb_clear(price);
    }

    // init
    //
    // to be called once market data and period is fixed.
    //
    // does not do any strike dependend computations
    void init() const {
        unsigned int iter = 0;
        prec /= 2;
        const slong additional_prec = 17;

        while (iter++ < 100) {

            prec *= 2;

            arb_set_ui(msqrt2, 2);
            arb_sqrt(msqrt2, msqrt2, prec);
            arb_neg(msqrt2, msqrt2); // set the const -sqrt(2)

            if (arb_rel_accuracy_bits(msqrt2) < target_prec + additional_prec)
                continue;

            // compute M
            arb_mul(tmp1, sigma, sigma, prec); // tmp1 = sigma^2

            arb_mul(tmp2, tmp1, t, prec);
            arb_exp(tmp2, tmp2, prec);
            arb_mul_ui(tmp2, tmp2, 2, prec); // tmp2 = 2exp(sigma^2*T)

            arb_sub(tmp4, t, tau, prec); // tmp4 = T - tau

            arb_mul(tmp3, tmp1, tmp4, prec);
            arb_exp(tmp3, tmp3, prec);
            arb_mul_ui(tmp3, tmp3, 2, prec); // tmp3 = 2exp(sigma^2*(t-tau))

            arb_mul(tmp5, tau, tmp1, prec);
            arb_add_ui(tmp5, tmp5, 1, prec); // tmp5 = 1 + sigma^2(tau)

            arb_mul(tmp3, tmp3, tmp5,
                    prec); // tmp3 = 2exp(sigma^2*(t-tau))(1+sigma^2(tau))

            arb_sub(tmp2, tmp2, tmp3, prec); // tmp2 = 2exp(sigma^2*t) -
            // 2exp(sigma^2*(t-tau))(1+sigma^2(tau))

            arb_mul(tmp1, tmp1, tmp1, prec);
            arb_mul(tmp4, tau, tau, prec);
            arb_mul(tmp1, tmp1, tmp4, prec); // tmp1 = sigma^4(tau)^2

            arb_div(M, tmp2, tmp1, prec);

            if (arb_rel_accuracy_bits(M) < target_prec + additional_prec)
                continue;

            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            arb_log(tmp1, M, prec);
            arb_div(tmp1, tmp1, t, prec);
            arb_sqrt(sigmaA, tmp1, prec); // sigmaA = sqrt( ln(M)/T )

            if (arb_rel_accuracy_bits(sigmaA) < target_prec + additional_prec)
                continue;

            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            arb_mul(tmp1, sigmaA, sigmaA, prec);
            arb_mul(tmp1, tmp1, t, prec);
            arb_div_ui(sigmaA2T2, tmp1, 2, prec); // sigma2T2 = sigmaA^2*T/2

            if (arb_rel_accuracy_bits(sigmaA2T2) < target_prec + additional_prec)
                continue;

            arb_sqrt(tmp2, t, prec);
            arb_mul(sigmaAsqrtT, tmp2, sigmaA,
                    prec); // sigmaAsqrtT = sigmaA * sqrt(T)

            if (arb_rel_accuracy_bits(sigmaAsqrtT) < target_prec + additional_prec)
                continue;

            //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            arb_mul(tmp1, r, t, prec);
            arb_neg(tmp1, tmp1);
            arb_exp(emrt, tmp1, prec);
            arb_div_ui(emrt, emrt, 2, prec); // emrt = exp(-rT)/2

            if (arb_rel_accuracy_bits(emrt) < target_prec + additional_prec)
                continue;

            break;
        } // iter-loop

        if (iter > 100) THROW_EXCEPTION(unprecise_calculation_exception,
                                        "cannot do computations to required precision");
    }

    void reset() { prec = default_prec; }

    double compute_put_price(double strike) const {


        arb_set_d(x, strike);

        if (arf_is_nan(arb_midref(x))) {
            std::stringstream ss;

            ss << "Variable x equals nan in compute_put_price" << std::endl
                    << " s         = " << s << std::endl
                    << " x         = " << x << std::endl
                    << " sigma     = " << sigma << std::endl
                    << " r         = " << r << std::endl
                    << " t         = " << t << std::endl
                    << " tau       = " << tau << std::endl
                    << std::endl;
            THROW_EXCEPTION(not_normal_number_exception, ss.str().data());
        }


        double result = 1. / 0.;

        if (prec != default_prec) {
            prec = default_prec;
            init();
        }

        unsigned int k;
        for (k = 0; k < 100; k++) {
            compute_d_values();

            arb_neg(tmp1, d2);
            arb_neg(tmp2, d1);

            arb_hypgeom_erfc(tmp1, tmp1, prec);
            arb_hypgeom_erfc(tmp2, tmp2, prec);

            arb_mul(tmp1, tmp1, x, prec);
            arb_mul(tmp2, tmp2, s, prec);

            arb_sub(tmp1, tmp1, tmp2, prec);

            arb_mul(price, tmp1, emrt,
                    prec); // price = exp(-rT)(x*erfc(-d2)-s*erfc(-d1))/2

            // check if we have enough precision in our computation
            if (arb_rel_accuracy_bits(price) >= target_prec) {
                // in case we do extract a double, save the result, and exit the loop
                result = arf_get_d(arb_midref(price), ARF_RND_NEAR);
                break;
            } else {
                // in case we don't, increase the working precision by a
                // factor of two and redo the complete computation
                prec *= 2;
                init();
            }
        }

        return result;
    }

    double get_sigmaA() const {

        double result;
        while (true) {
            if (arb_rel_accuracy_bits(sigmaA) >= target_prec) {
                // in case we do extract a double, save the result, and exit the loop
                result = arf_get_d(arb_midref(sigmaA), ARF_RND_NEAR);
                break;
            } else {
                // in case we don't, increase the working precision by a
                // factor of two and redo the complete computation
                prec *= 2;
                init();
            }
        }

        return result;
    }

    unsigned int get_number_of_settlement_levels() { return 61; }

    std::vector<double> compute_settlement_levels() {
        double v = arf_get_d(arb_midref(sigmaAsqrtT), ARF_RND_NEAR);
        double s0 = arf_get_d(arb_midref(s), ARF_RND_NEAR);

        double s = -3.0 * v;
        v /= 10.;

        std::vector<double> ret;
        for (unsigned int i = 0; i < 61; i++) {
            ret.push_back(s0 * std::exp(s));
            s += v;
        }

        return ret;
    }


    double compute_cdf_value_for_settlement_at(const double settlement) const {

        double result;

        while (true) {
            arb_set_d(tmp1, settlement);
            arb_div(tmp1, tmp1, s, prec);
            arb_log(tmp1, tmp1, prec);

            arb_set(tmp2, sigmaA);
            arb_set_d(tmp3, 2);
            arb_mul(tmp2, tmp2, tmp2, prec);
            arb_mul(tmp2, tmp2, t, prec);
            arb_div(tmp2, tmp2, tmp3, prec);

            arb_add(tmp1, tmp1, tmp2, prec);
            arb_div(x, tmp1, sigmaAsqrtT, prec);

            arb_sqrt(tmp4, tmp3, prec);
            arb_div(x, x, tmp4, prec);
            arb_neg(x, x);
            arb_hypgeom_erfc(x, x, prec);
            arb_div(tmp1, x, tmp3, prec);

            // check if we have enough precision in our computation
            if (arb_rel_accuracy_bits(tmp1) >= target_prec) {
                // in case we do extract a double, save the result, and exit the loop
                result = arf_get_d(arb_midref(tmp1), ARF_RND_NEAR);
                break;
            } else {
                // in case we don't, increase the working precision by a
                // factor of two and redo the complete computation
                prec *= 2;
                init();
            }
        }

        return result;

    }


    double compute_pdf_value_for_settlement_at(const double settlement) const {

        double result;

        while (true) {
            arb_set_d(tmp1, settlement);
            arb_div(tmp1, tmp1, s, prec);
            arb_log(tmp1, tmp1, prec);

            arb_set(tmp2, sigmaA);
            arb_set_d(tmp3, 2);
            arb_mul(tmp2, tmp2, tmp2, prec);
            arb_mul(tmp2, tmp2, t, prec);
            arb_div(tmp2, tmp2, tmp3, prec);

            arb_add(tmp1, tmp1, tmp2, prec);
            arb_div(x, tmp1, sigmaAsqrtT, prec);

            arb_mul(tmp1, x, x, prec);
            arb_div(tmp1, tmp1, tmp3, prec);
            arb_neg(tmp1, tmp1);
            arb_exp(tmp1, tmp1, prec);

            arb_const_pi(tmp2, prec);
            arb_mul(tmp2, tmp2, tmp3, prec);
            arb_sqrt(tmp2, tmp2, prec);
            arb_div(tmp1, tmp1, tmp2, prec);

            // check if we have enough precision in our computation
            mp_limb_signed_t n = arb_rel_accuracy_bits(tmp1);
            if (n >= target_prec) {
                // in case we do extract a double, save the result, and exit the loop
                result = arf_get_d(arb_midref(tmp1), ARF_RND_NEAR);
                break;
            } else {
                // in case we don't, increase the working precision by a
                // factor of two and redo the complete computation
                prec *= 2;
                init();
            }
        }

        return result;
    }

    double compute_call_price(double strike) const {
        arb_set_d(x, strike);

        if (arf_is_nan(arb_midref(x))) {
            std::stringstream ss;
            ss << "variable x equals nan in compute_put_price()" << std::endl
                    << " s         = " << s << std::endl
                    << " x         = " << x << std::endl
                    << " sigma     = " << sigma << std::endl
                    << " r         = " << r << std::endl
                    << " t         = " << t << std::endl
                    << " tau       = " << tau << std::endl
                    << std::endl;
            THROW_EXCEPTION(not_normal_number_exception, ss.str().data());
        }

        double result = 1. / 0.;

        if (prec != default_prec) {
            prec = default_prec;
            init();
        }

        unsigned int k;
        for (k = 0; k < 10; k++) {
            compute_d_values();

            arb_hypgeom_erfc(tmp1, d1, prec);
            arb_hypgeom_erfc(tmp2, d2, prec);

            arb_mul(tmp1, tmp1, s, prec);
            arb_mul(tmp2, tmp2, x, prec);

            arb_sub(tmp1, tmp1, tmp2, prec);

            arb_mul(price, tmp1, emrt,
                    prec); // price exp(-rT)(s*erfc(d1)-x*erfc(d2))/2

            // check if we have enough precision in our computation
            if (arb_rel_accuracy_bits(price) >= target_prec) {
                // in case we do extract a double, save the result, and exit the loop
                result = arf_get_d(arb_midref(price), ARF_RND_NEAR);
                break;
            } else {
                // in case we don't, increase the working precision by a
                // factor of two and redo the complete computation
                prec *= 2;
                init();
            }
        }

        return result;
    }

    double get_emrt() const {
        return arf_get_d(arb_midref(emrt), ARF_RND_NEAR);
    }

};
