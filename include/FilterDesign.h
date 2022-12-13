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

enum filter_type {
    butterworth,
    chebyshev1,
    chebyshev2,
    elliptic,
    bessel
};

using design_res = std::tuple<std::vector<mpcomplex>, std::vector<mpcomplex>, mpfr::mpreal, std::vector<std::array<mpfr::mpreal, 3>>, std::vector<std::array<mpfr::mpreal, 3>>>;
using zeros = std::vector<mpcomplex>;
using poles = std::vector<mpcomplex>;
//using zero_and_pole = std::tuple<zeros, poles>;

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
        ellipitic_filter_order(const mpreal &wpu, const mpreal &wpl, const mpreal &wsu, const mpreal &wsl,
                               const mpreal &Ap,
                               const mpreal &As,
                               filter_band_type type = filter_band_type::bandpass) -> std::tuple<uint32_t, mpreal, mpreal>;

        //从零极点设计出低通滤波器
        auto zp_trans(zeros &zs, poles &ps, const mpreal &Gp) -> design_res;

        //低通转带通
        auto lp2bp(zeros &zs, poles &ps, const std::vector<std::array<mpreal, 3>> &B,
                   const std::vector<std::array<mpreal, 3>> &A,
                   const mpreal &w0, const mpreal &Gp = 1_mpr) -> design_res;

        //低通转带阻
        auto lp2bs(zeros &zs, poles &ps, const std::vector<std::array<mpreal, 3>> &B,
                   const std::vector<std::array<mpreal, 3>> &A,
                   const mpreal &w0, const mpreal &Gp = 1_mpr) -> design_res;
    }


    auto ellipitic_filter(const mpreal &Wp, const mpreal &Ws, const mpreal &Ap, const mpreal &As,
                          filter_band_type type = lowpass) -> design_res;

    auto ellipitic_filter(const mpreal &Wpu, const mpreal &Wpl, const mpreal &Wsu, const mpreal &Wsl, const mpreal &Ap,
                          const mpreal &As,
                          filter_band_type type = filter_band_type::bandpass) -> design_res;

    auto butterworth_filter(const mpreal &Wp, const mpreal &Ws, const mpreal &Ap, const mpreal &As,
                            filter_band_type type = lowpass) -> design_res;

    auto
    butterworth_filter(const mpreal &Wpu, const mpreal &Wpl, const mpreal &Wsu, const mpreal &Wsl, const mpreal &Ap,
                       const mpreal &As,
                       filter_band_type type = filter_band_type::bandpass) -> design_res;

    auto chebyshev1_filter(const mpreal &Wp, const mpreal &Ws, const mpreal &Ap, const mpreal &As,
                           filter_band_type type = lowpass) -> design_res;

    auto chebyshev1_filter(const mpreal &Wpu, const mpreal &Wpl, const mpreal &Wsu, const mpreal &Wsl, const mpreal &Ap,
                           const mpreal &As,
                           filter_band_type type = filter_band_type::bandpass) -> design_res;

    auto chebyshev2_filter(const mpreal &Wp, const mpreal &Ws, const mpreal &Ap, const mpreal &As,
                           filter_band_type type = lowpass) -> design_res;

    auto chebyshev2_filter(const mpreal &Wpu, const mpreal &Wpl, const mpreal &Wsu, const mpreal &Wsl, const mpreal &Ap,
                           const mpreal &As,
                           filter_band_type type = filter_band_type::bandpass) -> design_res;

}

namespace DF {
    using namespace mpfr;
    namespace detail {
        auto af2df(const zeros &zs, const poles &ps, const mpreal &c0, int q, const mpreal &Gp) -> design_res;
    }

    auto ellipitic_filter(const mpreal &wp, const mpreal &ws, const mpreal &Ap, const mpreal &As,
                          filter_band_type type = lowpass) -> design_res;

    auto ellipitic_filter(const mpreal &wpu, const mpreal &wpl, const mpreal &wsu, const mpreal &wsl, const mpreal &Ap,
                          const mpreal &As,
                          filter_band_type type = filter_band_type::bandpass) -> design_res;

    auto butterworth_filter(const mpreal &wp, const mpreal &ws, const mpreal &Ap, const mpreal &As,
                            filter_band_type type = lowpass) -> design_res;

    auto
    butterworth_filter(const mpreal &wpu, const mpreal &wpl, const mpreal &wsu, const mpreal &wsl, const mpreal &Ap,
                       const mpreal &As,
                       filter_band_type type = filter_band_type::bandpass) -> design_res;

    auto chebyshev1_filter(const mpreal &wp, const mpreal &ws, const mpreal &Ap, const mpreal &As,
                           filter_band_type type = lowpass) -> design_res;

    auto chebyshev1_filter(const mpreal &wpu, const mpreal &wpl, const mpreal &wsu, const mpreal &wsl, const mpreal &Ap,
                           const mpreal &As,
                           filter_band_type type = filter_band_type::bandpass) -> design_res;

    auto chebyshev2_filter(const mpreal &wp, const mpreal &ws, const mpreal &Ap, const mpreal &As,
                           filter_band_type type = lowpass) -> design_res;

    auto chebyshev2_filter(const mpreal &wpu, const mpreal &wpl, const mpreal &wsu, const mpreal &wsl, const mpreal &Ap,
                           const mpreal &As,
                           filter_band_type type = filter_band_type::bandpass) -> design_res;

}


#endif //ELLIPTICFILTER_FILTERDESIGN_H
