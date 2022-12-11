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

enum filter_band_type {
    lowpass,
    highpass,
    bandpass,
    bandstop
};

using design_res = std::tuple<std::vector<mpcomplex>, std::vector<mpcomplex>, mpfr::mpreal, std::vector<std::array<mpfr::mpreal, 3>>, std::vector<std::array<mpfr::mpreal, 3>>>;
using zeros = std::vector<mpcomplex>;
using poles = std::vector<mpcomplex>;
using zero_and_pole = std::tuple<zeros, poles>;

namespace AF {
    using namespace mpfr;
    namespace detail {
        auto elliptic_lp_prototype(uint32_t N, const mpreal &Ap,
                                   const mpreal &As) -> design_res;

        auto elliptic_lp_prototype(const mpreal &wp, const mpreal &ws, const mpreal &Ap,
                                   const mpreal &As) -> design_res;

        auto ellipitic_filter_order(const mpreal &wp, const mpreal &ws, const mpreal &Ap, const mpreal &As,
                                    filter_band_type type = lowpass) -> std::tuple<uint32_t, mpreal>;

        auto
        ellipitic_filter_order(const mpreal &wpu, const mpreal &wpl, const mpreal &wsu, const mpreal &wsl, const mpreal &Ap,
                               const mpreal &As,
                               filter_band_type type = filter_band_type::bandpass) -> std::tuple<uint32_t, mpreal, mpreal>;
    }


    auto ellipitic_filter(const mpreal &Wp, const mpreal &Ws, const mpreal &Ap, const mpreal &As,
                          filter_band_type type = lowpass) -> design_res;

    auto ellipitic_filter(const mpreal &Wpu, const mpreal &Wpl, const mpreal &Wsu, const mpreal &Wsl, const mpreal &Ap,
                          const mpreal &As,
                          filter_band_type type = filter_band_type::bandpass) -> design_res;
}

namespace DF {
    using namespace mpfr;
    namespace detail {

    }

    auto ellipitic_filter(const mpreal &wp, const mpreal &ws, const mpreal &Ap, const mpreal &As, filter_band_type type = lowpass) -> design_res;

    auto ellipitic_filter(const mpreal &wpu, const mpreal &wpl, const mpreal &wsu, const mpreal &wsl, const mpreal &Ap,
                          const mpreal &As,
                          filter_band_type type = filter_band_type::bandpass) -> design_res;

}


#endif //ELLIPTICFILTER_FILTERDESIGN_H
