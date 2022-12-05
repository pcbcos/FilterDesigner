//
// Created by wengz on 2022/11/11.
//

#ifndef ELLIPTICFILTER_FILTERDESIGN_H
#define ELLIPTICFILTER_FILTERDESIGN_H

#include <vector>
#include <array>
#include <tuple>
#include "mpreal_ex.h"
#include "EllipticFunction.h"

namespace AF {
    enum filter_band_type {
        lowpass,
        highpass,
        bandpass,
        bandstop
    };
    using namespace mpfr;
    namespace detail {
        auto elliptic_lp_prototype(uint32_t N, const mpreal &Ap,
                                   const mpreal &As) -> std::tuple<std::vector<mpcomplex>, std::vector<mpcomplex>, mpreal, std::vector<std::array<mpreal, 3>>, std::vector<std::array<mpreal, 3>>>;


        auto elliptic_lp_prototype(mpreal wp, mpreal ws, mpreal Ap,
                                   mpreal As) -> std::tuple<std::vector<mpcomplex>, std::vector<mpcomplex>, mpreal, std::vector<std::array<mpreal, 3>>, std::vector<std::array<mpreal, 3>>>;
    }

    auto ellipitic_filter_order(mpreal wp, mpreal ws, mpreal Ap, mpreal As,
                                filter_band_type type = lowpass) -> std::tuple<uint32_t, mpreal>;

    auto
    ellipitic_filter_order2(mpfr::mpreal wpu, mpfr::mpreal wpl, mpfr::mpreal wsu, mpfr::mpreal wsl, mpfr::mpreal Ap,
                            mpfr::mpreal As,
                            AF::filter_band_type type = filter_band_type::bandpass) -> std::tuple<uint32_t, mpreal, mpreal>;

    auto ellipitic_filter(mpreal wp, mpreal ws, mpreal Ap, mpreal As,
                          filter_band_type type = lowpass) -> std::tuple<std::vector<mpcomplex>, std::vector<mpcomplex>, mpreal, std::vector<std::array<mpreal, 3>>, std::vector<std::array<mpreal, 3>>>;

    auto ellipitic_filter(mpfr::mpreal wpu, mpfr::mpreal wpl, mpfr::mpreal wsu, mpfr::mpreal wsl, mpfr::mpreal Ap,
                          mpfr::mpreal As,
                          AF::filter_band_type type = filter_band_type::bandpass) -> std::tuple<std::vector<mpcomplex>, std::vector<mpcomplex>, mpreal, std::vector<std::array<mpreal, 3>>, std::vector<std::array<mpreal, 3>>>;
}

namespace DF {


}


#endif //ELLIPTICFILTER_FILTERDESIGN_H
