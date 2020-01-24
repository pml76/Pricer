/*
 *  (c) by Peter Lennartz
 * 
 *   This program is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   This program is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with this program.  If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef PRICER_MEMORY_H
#define PRICER_MEMORY_H

#include <include/pricer-base.h>

#include "allocate.h"

namespace Pricer {


#define DEFINE_VARIABLE(type,x) private:type m__ ## x; public:inline type& INLINE get_ ## x() {return m__ ## x;}


    /**
     * class pricer_context. Memory management for simple option-pricer.
     *
     * The memory managed by this class handles the in- and output of the simple option-pricer. See variables below
     * to find out which variables are input- and which are output variables.
     *
     * When computing prices the workflow is
     *    1. to call init_tw_pricer(),
     *    2. to fill the necessary arrays contained in this class,
     *    3. to call prep_tw_pricer(),
     *    4. do call tw_pricer().
     *
     * After modifying data the user must ensure that she calls prep_tw_pricer() before calling tw_pricer() again.
     */
    class pricer_context {
        friend class ddx_pricer_context;
        friend class compute_prices_of_instruments_context;
    public:
        pricer_context(uint64_t n_max = 1024 * 64) :
                m__n_act(0), m__n_max(0) {
                alloc_mem(n_max);
        }

        virtual ~pricer_context() {
            dealloc_mem();
        };
        

        inline void copy_entry(uint64_t i, pricer_context *p, uint64_t j ) {
            
            get_x()[j]           = p->get_x()[i];
            get_s()[j]           = p->get_s()[i];
            get_sigma()[j]       = p->get_sigma()[i];
            get_r()[j]           = p->get_r()[i];
            get_t()[j]           = p->get_t()[i];
            get_tau()[j]         = p->get_tau()[i];
            get_put_call()[j]    = p->get_put_call()[i];
            get_long_short()[j]  = p->get_long_short()[i];
            get_tradedate()[j]   = p->get_tradedate()[i];

            get_sigmaA()[j]      = p->get_sigmaA()[i];
            get_sigmaA2T2()[j]   = p->get_sigmaA2T2()[i];
            get_sigmaAsqrtT()[j] = p->get_sigmaAsqrtT()[i];
            get_emrt()[j]        = p->get_emrt()[i];
            get_tradedate()[j]   = p->get_tradedate()[i];
            get_prices()[j]      = p->get_prices()[i];
            get_d1()[j]          = p->get_d1()[i];
            get_d2()[j]          = p->get_d2()[i];


        }

        inline
        uint64_t add_entry(double tradedate_p, double x_p, double s_p, double sigma_p, double r_p, double t_p, double tau_p, double put_call_p, double long_short_p) {

            if(UNLIKELY(get_n_act() >= get_n_max())) {
                realloc_mem(2*get_n_max());
            }

            get_x()[get_n_act()] = x_p;
            get_s()[get_n_act()] = s_p;
            get_sigma()[get_n_act()] = sigma_p;
            get_r()[get_n_act()] = r_p;
            get_t()[get_n_act()] = t_p;
            get_tau()[get_n_act()] = tau_p;
            get_put_call()[get_n_act()] = put_call_p;
            get_long_short()[get_n_act()] = long_short_p;
            get_tradedate()[get_n_act()] = tradedate_p;

            return get_n_act()++;

        }

    private:

        virtual void alloc_mem(uint64_t n);
        virtual void realloc_mem(uint64_t n_p);
        virtual void dealloc_mem();

    public:

        DEFINE_VARIABLE(Real_Ptr, tradedate);
        DEFINE_VARIABLE(uint64_t, n_act);

        DEFINE_VARIABLE(uint64_t, n_max) /// the maximum number of entries in the arrays managed by this class.
                                         ///
                                         /// NOTE: the getter routines declared here return mostly pointers to the
                                         /// beginnings of arrays. When writing to the arrays it is in the
                                         /// responsibility of the routine that performs the writing to ensure that
                                         /// it does not write beyond the n_max-th entry of the arrays.
                                         ///

        DEFINE_VARIABLE(Real_Ptr, prices)   /// [output]. Contains the prices of the options after the call
                                            /// to tw_pricer().

        DEFINE_VARIABLE(Real_Ptr, d1)       /// [output]. Needed for computation of derivatives
        DEFINE_VARIABLE(Real_Ptr, d2)       /// [output]. Needed for computation of derivatives
        DEFINE_VARIABLE(Real_Ptr, long_short)   /// [input]. When equal to 1 than a long leg is priced,
                                                /// When equal to -1 than a short leg is priced
        DEFINE_VARIABLE(Real_Ptr, put_call)     /// [input]. When equal to 1 a call is priced. When
                                                /// equal to -1 a put is pricedj
        DEFINE_VARIABLE(Real_Ptr, x)        /// [input]. The strike of the option.
        DEFINE_VARIABLE(Real_Ptr, s)        /// [input]. The price of the underlying (in the term-month)
        DEFINE_VARIABLE(Real_Ptr, sigma)    /// [input]. The volatility of the underlying in the term-month
        DEFINE_VARIABLE(Real_Ptr, t)        /// [input]. The time to maturity
        DEFINE_VARIABLE(Real_Ptr, tau)      /// [input]. The length of the averaging period of the option
        DEFINE_VARIABLE(Real_Ptr, r)        /// [input]. The interest rate at time t
        DEFINE_VARIABLE(Real_Ptr, sigmaA)           /// [internal use].
        DEFINE_VARIABLE(Real_Ptr, sigmaA2T2)        /// [internal use].
        DEFINE_VARIABLE(Real_Ptr, sigmaAsqrtT)      /// [internal use].
        DEFINE_VARIABLE(Real_Ptr, emrt)             /// [internal use].
        // DEFINE_VARIABLE(Real_Ptr, d2dx2_prep)       /// [internal use].

    };



    /**
     * class ddx_pricer_context. Memory management for for computing option prices and their first
     * partial derivatives with respect to the strikes.
     *
     * The memory managed by this class contains the input and the output of the calculation of the option price
     * as well as the calculation of the first partial derivative with respect to the strike. See the class
     * pricer_context and below to find out which fields are input, which fields are output, and which are for
     * internal use.
     *
     * NOTE: A caller has to adhere to the following workflow:
     *      1. call init_tw_pricer()
     *      2. create a dd_pricer_context opject and fill the necessary arrays.
     *      3. call prep_tw_pricer()
     *      4. call tw_pricer()
     *      5. call ddx_tw_pricer()
     *
     * After modifying the input data, steps 3. to 5. have to be repeated to ensure that everything
     * in upto date.
     */
    class ddx_pricer_context : public pricer_context {
        friend class d2dx2_pricer_context;
    public:
        ddx_pricer_context(uint64_t n_max  = 1024 * 64) : pricer_context( n_max ){
            alloc_mem( n_max);
        }

        virtual ~ddx_pricer_context() {
            dealloc_mem();
        };

        inline
        void copy_entry(uint64_t i, ddx_pricer_context *p, uint64_t j ) {
            pricer_context::copy_entry(i, p, j);
            get_ddx_price()[j] = p->get_ddx_price()[i];
        }


    private:

        virtual void alloc_mem(uint64_t n) override ;
        virtual void realloc_mem(uint64_t n_p) override ;
        virtual void dealloc_mem() override ;

    public:

        DEFINE_VARIABLE(Real_Ptr, ddx_price)  /// [output]. The value of the first partial derivative
                                              /// of the price function with respect to the strike.

    };


    /**
     * class d2dx2_pricer_context. Memory management for computing option prices and their first and second
     * partial derivatives with respect to the strikes.
     *
     * The memory managed by this class contains the input and the output of the calculation of the option price
     * as well as the calculation of the first and second partial derivative of the price function with respect to
     * the strike. See the class ddx_pricer_context and below to find out which fields are input,
     * which fields are output, and which are for internal use.
     *
     * NOTE: Computation of the first derivative is optional. The computation of the price, however, is mandatory
     * in order to compute the values fo the second derivative.
     *
     * NOTE: A caller has to adhere to the following workflow:
     *      1. call init_tw_pricer()
     *      2. create a dd_pricer_context opject and fill the necessary arrays.
     *      3. call prep_tw_pricer()
     *      4. call tw_pricer()
     *      5. optionally call ddx_tw_pricer()
     *      6. call d2dx2_tw_pricer()
     *
     * After modifying the input data, steps 3. to 6. have to be repeated to ensure that everything
     * in upto date.
     */
    class d2dx2_pricer_context : public ddx_pricer_context {

    public:
        d2dx2_pricer_context(uint64_t n_max = 1024 * 64) : ddx_pricer_context( n_max ){
            alloc_mem( n_max);
        }

        virtual ~d2dx2_pricer_context() {
            dealloc_mem();
        };

        inline
        void copy_entry(uint64_t i, d2dx2_pricer_context *p, uint64_t j ) {
            ddx_pricer_context::copy_entry(i, p, j);
            get_d2dx2_price()[j] = p->get_d2dx2_price()[i];
        }

    private:

        virtual void alloc_mem(uint64_t n) override ;
        virtual void realloc_mem(uint64_t n_p) override ;
        virtual void dealloc_mem() override ;


    public:

        DEFINE_VARIABLE(Real_Ptr, d2dx2_price) /// [output]. Holds the values of the second derivatives of
                                               /// the price-function after the call to d2dx2_tw_pricer()

    };



    /**
     * Memory for m_max instruments with a total of at most n…max legs.
     */
    class compute_prices_of_instruments_context : public pricer_context {
        friend class compute_instrument_strikes_from_premiums_context;
    public:
        compute_prices_of_instruments_context(uint64_t n_max  = 1024 * 64, uint64_t m_max = 1024 * 64) : pricer_context( n_max ){
                alloc_mem( n_max, m_max );
                get_m_act() = 0;
        }

        virtual ~compute_prices_of_instruments_context() {
            dealloc_mem();
        };


        void shrink_n_max_to_actual() {
            if(get_n_act() % 64 == 0 ) {
                get_n_max() = get_n_act();
            } else {
                get_n_max() = ((get_n_act() / 64)+1)*64;
            }
        }

        void shrink_m_max_to_actual() {
            if(get_m_act() % 64 == 0 ) {
                get_m_max() = get_m_act();
            } else {
                get_m_max() = ((get_m_act() / 64)+1)*64;
            }
        }

        uint64_t add_structure(double x_p ) {
            if(UNLIKELY(get_m_act() >= get_m_max() - 1)) {
                realloc_mem(get_n_max(), 2*get_m_max());
            }

            m__x_[m__m_act] = x_p;
            return m__m_act++;

        }

        inline
        void copy_structure(uint64_t i, compute_prices_of_instruments_context *p, uint64_t j ) {
            get_x_()[j] = p->get_x_()[i];
            get_instrument_prices()[j] = p->get_instrument_prices()[i];
        }


        inline
        void copy_leg(uint64_t i, compute_prices_of_instruments_context *p, uint64_t j ) {

            get_to_structure()[j] = p->get_to_structure()[i];
            get_offsets()[j]      = p->get_offsets()[i];
            pricer_context::copy_entry(i, p, j);

        }


        inline
        uint64_t  add_leg(double tradedate_p, double x_p, double s_p, double sigma_p, double r_p, double t_p, double tau_p, double put_call_p, double long_short_p, double offset_p, uint64_t structure_p) {


            get_to_structure()[get_n_act()] = structure_p;
            get_offsets()[get_n_act()] = offset_p;
            return pricer_context::add_entry(tradedate_p, x_p, s_p, sigma_p, r_p, t_p, tau_p, put_call_p, long_short_p);
        }

    private:

        virtual void alloc_mem(uint64_t n, uint64_t m);
        virtual void realloc_mem(uint64_t n_p, uint64_t m_p);
        virtual void dealloc_mem() override;

    public:
        DEFINE_VARIABLE(Int32_Ptr, to_structure)     /// [input]
        DEFINE_VARIABLE(Real_Ptr, offsets)           /// [input]
        DEFINE_VARIABLE(Real_Ptr, x_)                /// [input]

        DEFINE_VARIABLE(uint64_t, m_act)
        DEFINE_VARIABLE(uint64_t, m_max)             ///
        DEFINE_VARIABLE(Real_Ptr, instrument_prices) /// [output]
    };



    class compute_instrument_strikes_from_premiums_context : public compute_prices_of_instruments_context {

    public:
        compute_instrument_strikes_from_premiums_context(uint64_t n_max = 1024 * 64, uint64_t m_max = 1024 * 64)
            : compute_prices_of_instruments_context(n_max, m_max) {
            alloc_mem( n_max, m_max );
        }

        virtual ~compute_instrument_strikes_from_premiums_context() {
            dealloc_mem();
        };


        uint64_t add_structure(double x_p, double premium_p ) {
            compute_prices_of_instruments_context::add_structure(x_p);

            get_premiums()[get_m_act()] = premium_p;
            return get_m_act()++;

        }


        inline
        void copy_structure(uint64_t i, compute_instrument_strikes_from_premiums_context *p, uint64_t j ) {

            compute_prices_of_instruments_context::copy_structure(i, p, j);

            get_premium_modifier()[j]    = p->get_premium_modifier()[i];
            get_premiums()[j]            = p->get_premiums()[i];
            get_instrument_pricesl()[j]  = p->get_instrument_pricesl()[i];
            get_instrument_pricesh()[j]  = p->get_instrument_pricesh()[i];
            get_xl_()[j]                 = p->get_xl_()[i];
            get_xh_()[j]                 = p->get_xh_()[i];

        }



    private:

        virtual void alloc_mem(uint64_t n, uint64_t m) override ;
        virtual void realloc_mem(__uint64_t n_p, uint64_t m_p) override ;
        virtual void dealloc_mem() override ;


    public:


        DEFINE_VARIABLE(Real_Ptr, premium_modifier)    /// [internal use]
        DEFINE_VARIABLE(Real_Ptr, premiums)            /// [input]
        DEFINE_VARIABLE(Real_Ptr, instrument_pricesl)  /// [internal use]
        DEFINE_VARIABLE(Real_Ptr, instrument_pricesh)  /// [internal use]
        DEFINE_VARIABLE(Real_Ptr, xl_)                 /// [internal use]
        DEFINE_VARIABLE(Real_Ptr, xh_)                 /// [internal use]
    };


}


#endif //PRICER_MEMORY_H
