//
// Created by wengz on 2022/12/5.
//
#include <utility>
#include <vector>
#include "EllipticFunction.h"
#include "FilterDesign.h"

using mpfr::mpreal;

auto DF::ellipitic_filter(const mpreal &wp, const mpreal &ws, const mpreal &Ap, const mpreal &As,
                          filter_band_type type) -> design_res {
    if (type == lowpass) {
        //转换技术指标DF->AF
        mpreal Wp = tan(wp / 2_mpr);
        mpreal Ws = tan(ws / 2_mpr);
        mpreal Gp = pow(10_mpr, -Ap / 20_mpr);
        //设计模拟低通滤波器
        auto [z0, p0, H0, B0, A0] = AF::ellipitic_filter(Wp, Ws, Ap, As, lowpass);
        //反转换AF->DF
        return DF::detail::af2df(z0, p0, 1_mpr, 1, Gp);
    } else if (type == highpass) {
        //转换技术指标
        mpreal Wp = cot(wp / 2_mpr);
        mpreal Ws = cot(ws / 2_mpr);
        mpreal Gp = pow(10_mpr, -Ap / 20_mpr);
        //设计模拟低通滤波器
        auto [z0, p0, H0, B0, A0] = AF::ellipitic_filter(Wp, Ws, Ap, As, lowpass);
        //反转换
        return DF::detail::af2df(z0, p0, -1_mpr, 1, Gp);
    }
    return {};
}

auto DF::ellipitic_filter(const mpreal &wpu, const mpreal &wpl, const mpreal &wsu, const mpreal &wsl, const mpreal &Ap,
                          const mpreal &As,
                          filter_band_type type) -> design_res {
    mpreal Gp = pow(10_mpr, -Ap / 20_mpr);
    mpreal c0, Wsl, Wsu, Ws, Wp, Wpl, Wpu;
    int q;
    if (type == bandpass) {
        //转换技术指标//TODO:还有一种技术指标的转换方式
        c0 = sin(wpl + wpu) / (sin(wpl) + sin(wpu));
        Wsl = (c0 - cos(wsl)) / sin(wsl);
        Wsu = (c0 - cos(wsu)) / sin(wsu);
        Ws = std::min(abs(Wsl), abs(Wsu));
        Wp = tan((wpu - wpl) / 2_mpr);
        q = 1;
//        mpreal c0=sin(wsl+wsu)/(sin(wsl)+sin(wsu));
//        mpreal Wpl=(c0-cos(wpl))/sin(wpl);
//        mpreal Wpu=(c0-cos(wpu))/sin(wpu);
//        mpreal Ws=tan((wsu-wsl)/2_mpr);
//        mpreal Wp=std::max(abs(Wpl),abs(Wpu));
    } else if (type == bandstop) {
        //转换技术指标//TODO:还有一种技术指标的转换方式
//        c0 = sin(wpl + wpu) / (sin(wpl) + sin(wpu));
//        Wsl = sin(wsl) / (c0 - cos(wsl));
//        Wsu = sin(wsu) / (c0 - cos(wsu));
//        Ws = std::min(abs(Wsl), abs(Wsu));
//        Wp = cot((wpu - wpl) / 2_mpr);
        c0 = sin(wsl + wsu) / (sin(wsl) + sin(wsu));
        Wpl = sin(wpl) / (c0 - cos(wpl));
        Wpu = sin(wpu) / (c0 - cos(wpu));
        Wp = std::max(abs(Wpl), abs(Wpu));
        Ws = cot((wsu - wsl) / 2_mpr);
        q = -1;
    }
    //设计模拟低通滤波器
    auto [z0, p0, H0, B0, A0] = AF::ellipitic_filter(Wp, Ws, Ap, As);
    //反转换
    return DF::detail::af2df(z0, p0, c0, q, Gp);
}

auto
DF::butterworth_filter(const mpfr::mpreal &wp, const mpfr::mpreal &ws, const mpfr::mpreal &Ap, const mpfr::mpreal &As,
                       filter_band_type type) -> design_res {
    if (type == lowpass) {
        //转换技术指标DF->AF
        mpreal Wp = tan(wp / 2_mpr);
        mpreal Ws = tan(ws / 2_mpr);
        mpreal Gp = 1_mpr;
        //设计模拟低通滤波器
        auto [z0, p0, H0, B0, A0] = AF::butterworth_filter(Wp, Ws, Ap, As, lowpass);
        //反转换AF->DF
        return DF::detail::af2df(z0, p0, 1_mpr, 1, Gp);
    } else if (type == highpass) {
        //转换技术指标
        mpreal Wp = cot(wp / 2_mpr);
        mpreal Ws = cot(ws / 2_mpr);
        mpreal Gp = 1_mpr;
        //设计模拟低通滤波器
        auto [z0, p0, H0, B0, A0] = AF::butterworth_filter(Wp, Ws, Ap, As, lowpass);
        //反转换
        return DF::detail::af2df(z0, p0, -1_mpr, 1, Gp);
    }
    return {};
}

auto DF::butterworth_filter(const mpfr::mpreal &wpu, const mpfr::mpreal &wpl, const mpfr::mpreal &wsu,
                            const mpfr::mpreal &wsl, const mpfr::mpreal &Ap, const mpfr::mpreal &As,
                            filter_band_type type) -> design_res {
    mpreal Gp = 1_mpr;
    mpreal c0, Wsl, Wsu, Ws, Wp, Wpl, Wpu;
    int q;
    if (type == bandpass) {
        //转换技术指标//TODO:还有一种技术指标的转换方式
        c0 = sin(wpl + wpu) / (sin(wpl) + sin(wpu));
        Wsl = (c0 - cos(wsl)) / sin(wsl);
        Wsu = (c0 - cos(wsu)) / sin(wsu);
        Ws = std::min(abs(Wsl), abs(Wsu));
        Wp = tan((wpu - wpl) / 2_mpr);
        q = 1;
//        mpreal c0=sin(wsl+wsu)/(sin(wsl)+sin(wsu));
//        mpreal Wpl=(c0-cos(wpl))/sin(wpl);
//        mpreal Wpu=(c0-cos(wpu))/sin(wpu);
//        mpreal Ws=tan((wsu-wsl)/2_mpr);
//        mpreal Wp=std::max(abs(Wpl),abs(Wpu));
    } else if (type == bandstop) {
        //转换技术指标//TODO:还有一种技术指标的转换方式
//        c0 = sin(wpl + wpu) / (sin(wpl) + sin(wpu));
//        Wsl = sin(wsl) / (c0 - cos(wsl));
//        Wsu = sin(wsu) / (c0 - cos(wsu));
//        Ws = std::min(abs(Wsl), abs(Wsu));
//        Wp = cot((wpu - wpl) / 2_mpr);
        c0 = sin(wsl + wsu) / (sin(wsl) + sin(wsu));
        Wpl = sin(wpl) / (c0 - cos(wpl));
        Wpu = sin(wpu) / (c0 - cos(wpu));
        Wp = std::max(abs(Wpl), abs(Wpu));
        Ws = cot((wsu - wsl) / 2_mpr);
        q = -1;
    }
    //设计模拟低通滤波器
    auto [z0, p0, H0, B0, A0] = AF::butterworth_filter(Wp, Ws, Ap, As);
    //反转换
    return DF::detail::af2df(z0, p0, c0, q, Gp);
}

auto
DF::chebyshev1_filter(const mpfr::mpreal &wp, const mpfr::mpreal &ws, const mpfr::mpreal &Ap, const mpfr::mpreal &As,
                      filter_band_type type) -> design_res {
    return design_res();
}

auto DF::chebyshev1_filter(const mpfr::mpreal &wpu, const mpfr::mpreal &wpl, const mpfr::mpreal &wsu,
                           const mpfr::mpreal &wsl, const mpfr::mpreal &Ap, const mpfr::mpreal &As,
                           filter_band_type type) -> design_res {
    return design_res();
}

auto
DF::chebyshev2_filter(const mpfr::mpreal &wp, const mpfr::mpreal &ws, const mpfr::mpreal &Ap, const mpfr::mpreal &As,
                      filter_band_type type) -> design_res {
    return design_res();
}

auto DF::chebyshev2_filter(const mpfr::mpreal &wpu, const mpfr::mpreal &wpl, const mpfr::mpreal &wsu,
                           const mpfr::mpreal &wsl, const mpfr::mpreal &Ap, const mpfr::mpreal &As,
                           filter_band_type type) -> design_res {
    return design_res();
}

auto DF::detail::af2df(const zeros &zs, const poles &ps, const mpreal &c0, int q, const mpreal &Gp) -> design_res {
    std::vector<mpcomplex> z;
    std::vector<mpcomplex> p;
    std::vector<std::array<mpreal, 3>> B1;//分子
    std::vector<std::array<mpreal, 3>> A1;//分母
    uint32_t r = ps.size() % 2;
    uint32_t L = ps.size() / 2;

    std::vector<mpcomplex> z_hat;
    std::vector<mpcomplex> p_hat;

    mpreal H0 = Gp;
    for (auto &z_ai: zs) {
        if (mpfr::isinf(z_ai.real())) {
            z_hat.push_back(-1_mpr);
        } else {
            z_hat.push_back((1_mpr + z_ai) / (1_mpr - z_ai));
        }
    }
    for (auto &p_ai: ps) {
        p_hat.push_back((1_mpr + p_ai) / (1_mpr - p_ai));
    }
    if (r == 1) {
        auto p_h0 = p_hat[0];
        mpreal G0 = (1_mpr - p_h0).real() / 2_mpr;
        mpcomplex p01 = (c0 * (1_mpr + mpfr::mpreal(q) * p_h0) +
                         sqrt(c0 * c0 * (1_mpr + mpreal(q) * p_h0) * (1_mpr + mpreal(q) * p_h0) - 4_mpr * q * p_h0)) /
                        2_mpr;
        mpcomplex p02 = (c0 * (1_mpr + mpfr::mpreal(q) * p_h0) -
                         sqrt(c0 * c0 * (1_mpr + mpreal(q) * p_h0) * (1_mpr + mpreal(q) * p_h0) - 4_mpr * q * p_h0)) /
                        2_mpr;
        p.push_back(p01);
        p.push_back(p02);
        mpcomplex z01 = (q == 1) ? mpcomplex(1, 0) : exp(mpcomplex(0, acos(c0)));
        mpcomplex z02 = (q == 1) ? mpcomplex(-1, 0) : exp(mpcomplex(0, -acos(c0)));
        z.push_back(z01);
        z.push_back(z02);
        if (q == 1) {
            B1.push_back({G0, 0, -G0});
        } else {
            B1.push_back({G0, -2 * c0 * G0, G0});
        }
        A1.push_back({1, -(p01 + p02).real(), (p01 * p02).real()});
        H0 *= G0;
    } else {
        B1.push_back({Gp, 0, 0});
        A1.push_back({1, 0, 0});
    }
    for (uint32_t i = r; i < 2 * L + r; i += 2) {
        auto p_hi = p_hat[i];
        auto z_hi = z_hat[i - r];
        mpreal Gi_abs = abs((1_mpr - p_hi)) / abs((1_mpr - z_hi));
        mpreal Gi_square = Gi_abs * Gi_abs;
        mpcomplex zi1, zi1s, zi2, zi2s;
        auto getsign = [](const mpreal &c) {
            return mpfr::signbit(c) ? -1_mpr : 1_mpr;
        };
        //第一组
        mpcomplex pi1 = (c0 * (1_mpr + mpreal(q) * p_hi) - getsign(c0) *
                                                           sqrt(c0 * c0 * (1_mpr + mpreal(q) * p_hi) *
                                                                (1_mpr + mpreal(q) * p_hi) -
                                                                4_mpr * mpreal(q) * p_hi)) / 2_mpr;
        mpcomplex pi1s = conj(pi1);

        //第二组
        //mpcomplex pi2 = (c0 - pi1) / (1_mpr - c0 * pi1);
        mpcomplex pi2 = c0 * (1_mpr + mpreal(q) * p_hi) - pi1;
        mpcomplex pi2s = conj(pi2);

        if (mpfr::isinf(zs[i - r].real())) {
            zi1 = (q == 1) ? mpcomplex(1, 0) : exp(mpcomplex(0, acos(c0)));
            zi2 = (q == 1) ? mpcomplex(-1, 0) : exp(mpcomplex(0, -acos(c0)));
            zi1s = conj(zi1);
            zi2s = conj(zi2);
        } else {
            zi1 = (c0 * (1_mpr + mpreal(q) * z_hi) -getsign(c0)*
                   sqrt(c0 * c0 * (1_mpr + mpreal(q) * z_hi) * (1_mpr + mpreal(q) * z_hi) -
                        4_mpr * mpreal(q) * z_hi)) / 2_mpr;
            zi1s = conj(zi1);
            zi2 = c0 * (1_mpr + mpreal(q) * z_hi) - zi1;;
            zi2s = conj(zi2);
        }

        p.push_back(pi1);
        p.push_back(pi1s);
        z.push_back(zi1);
        z.push_back(zi1s);

        p.push_back(pi2);
        p.push_back(pi2s);
        z.push_back(zi2);
        z.push_back(zi2s);
        if (q == 1 && (c0 == 1_mpr) || (c0 == -1_mpr)) {//LP HP特判
            if (mpfr::isinf(zs[i - r].real())) {//Butterworth Chebyshev1 特判
                B1.push_back({Gi_square, 2*Gi_square*c0, Gi_square});
            } else {
                B1.push_back({Gi_square, -2_mpr * Gi_square * zi1.real(), Gi_square * abs(zi1) * abs(zi1)});
            }
            A1.push_back({1_mpr, -2_mpr * pi1.real(), abs(pi1) * abs(pi1)});
        } else {

            B1.push_back({Gi_abs, -2_mpr * Gi_abs * zi1.real(), Gi_abs * abs(zi1) * abs(zi1)});
            B1.push_back({Gi_abs, -2_mpr * Gi_abs * zi2.real(), Gi_abs * abs(zi2) * abs(zi2)});

            //B1.push_back({Gi_abs, -2_mpr * Gi_abs * zi1.real(), Gi_abs * abs(zi1) * abs(zi1)});
            //B1.push_back({Gi_abs, -2_mpr * Gi_abs * zi2.real(), Gi_abs * abs(zi2) * abs(zi2)});

            A1.push_back({1_mpr, -2_mpr * pi1.real(), abs(pi1) * abs(pi1)});
            A1.push_back({1_mpr, -2_mpr * pi2.real(), abs(pi2) * abs(pi2)});
        }

        H0 *= Gi_square;
    }
    return {z, p, H0, B1, A1};

}