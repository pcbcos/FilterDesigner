//
// Created by wengz on 2022/12/5.
//
#include <utility>
#include <vector>
#include "EllipticFunction.h"
#include "FilterDesign.h"

auto DF::ellipitic_filter(mpfr::mpreal wp, mpfr::mpreal ws, mpfr::mpreal Ap, mpfr::mpreal As,
                          filter_band_type type) -> std::tuple<std::vector<mpcomplex>, std::vector<mpcomplex>, mpreal, std::vector<std::array<mpreal, 3>>, std::vector<std::array<mpreal, 3>>> {
    if (type == lowpass) {
        //转换技术指标
        mpreal Wp = tan(wp / 2_mpr);
        mpreal Ws = tan(ws / 2_mpr);
        mpreal Gp = pow(10_mpr, -Ap / 20_mpr);
        //设计模拟低通滤波器
        auto [z0, p0, H0, B0, A0] = AF::ellipitic_filter(Wp, Ws, std::move(Ap), std::move(As), lowpass);
        //反转换
        std::vector<mpcomplex> z;
        std::vector<mpcomplex> p;
        std::vector<std::array<mpreal, 3>> B;//分子
        std::vector<std::array<mpreal, 3>> A;//分母
        uint32_t r = p0.size() % 2;
        uint32_t L = p0.size() / 2;
        if (r == 1) {
            auto p_a0 = p0[0];
            auto p_0 = (1_mpr + p_a0) / (1_mpr - p_a0);
            mpreal G0 = (1_mpr - p_0).real() / 2_mpr;
            H0 *= G0;
            p.push_back(p_0);
            B.push_back({H0 * G0, H0 * G0, 0});
            A.push_back({1, -p_0.real(), 0});
        } else {
            B.push_back({Gp, 0, 0});
            A.push_back({1, 0, 0});
        }
        for (uint32_t i = r; i < 2 * L + r; i += 2) {
            auto z_ai = z0[i - r];
            auto p_ai = p0[i];

            auto zi1 = (1_mpr + z_ai) / (1_mpr - z_ai);
            auto pi1 = (1_mpr + p_ai) / (1_mpr - p_ai);

            z.push_back(zi1);
            z.push_back(conj(zi1));
            p.push_back(pi1);
            p.push_back(conj(pi1));

            mpreal Gi_abs = (1_mpr - pi1).real() / (1_mpr - zi1).real();
            mpreal Gi_square = Gi_abs * Gi_abs;
            H0 *= Gi_square;
            B.push_back({Gi_square, -2_mpr * Gi_square * zi1.real(), Gi_square * abs(zi1) * abs(zi1)});
            A.push_back({1, -2_mpr * pi1.real(), abs(pi1) * abs(pi1)});
        }
        return {z, p, H0, B, A};

    } else if (type == highpass) {
        //转换技术指标
        mpreal Wp = cot(wp / 2_mpr);
        mpreal Ws = cot(ws / 2_mpr);
        mpreal Gp = pow(10_mpr, -Ap / 20_mpr);
        //设计模拟低通滤波器
        auto [z0, p0, H0, B0, A0] = AF::ellipitic_filter(Wp, Ws, std::move(Ap), std::move(As), lowpass);
        //反转换
        std::vector<mpcomplex> z;
        std::vector<mpcomplex> p;
        std::vector<std::array<mpreal, 3>> B;//分子
        std::vector<std::array<mpreal, 3>> A;//分母
        uint32_t r = p0.size() % 2;
        uint32_t L = p0.size() / 2;
        if (r == 1) {
            auto p_a0 = p0[0];
            auto p_i0 = -((1_mpr + p_a0) / (1_mpr - p_a0));
            mpreal G0 = (1_mpr + p_i0.real()) / 2_mpr;
            H0 *= G0;
            p.push_back(p_i0);
            B.push_back({G0, -G0, 0});
            A.push_back({1, -p_i0.real(), 0});
        } else {
            B.push_back({Gp, 0, 0});
            A.push_back({1, 0, 0});
        }
        for (uint32_t i = r; i < 2 * L + r; i += 2) {
            auto z_ai = z0[i - r];
            auto p_ai = p0[i];

            auto zi1 = -(1_mpr + z_ai) / (1_mpr - z_ai);
            auto pi1 = -(1_mpr + p_ai) / (1_mpr - p_ai);

            z.push_back(zi1);
            z.push_back(conj(zi1));
            p.push_back(pi1);
            p.push_back(conj(pi1));

            mpreal Gi_abs = (1_mpr + pi1).real() / (1_mpr + zi1).real();
            mpreal Gi_square = Gi_abs * Gi_abs;
            H0 *= Gi_square;
            B.push_back({Gi_square, -2_mpr * Gi_square * zi1.real(), Gi_square * abs(zi1) * abs(zi1)});
            A.push_back({1, -2_mpr * pi1.real(), abs(pi1) * abs(pi1)});
        }
        return {z, p, H0, B, A};
    }

    return {};
}

auto DF::ellipitic_filter(mpfr::mpreal wpu, mpfr::mpreal wpl, mpfr::mpreal wsu, mpfr::mpreal wsl, mpfr::mpreal Ap,
                          mpfr::mpreal As,
                          filter_band_type type) -> std::tuple<std::vector<mpcomplex>, std::vector<mpcomplex>, mpreal, std::vector<std::array<mpreal, 3>>, std::vector<std::array<mpreal, 3>>> {
    if (type == bandpass) {
        //转换技术指标//TODO:还有一种技术指标的转换方式
        mpreal c0 = sin(wpl + wpu) / (sin(wpl) + sin(wpu));
        mpreal Wsl = (c0 - cos(wsl)) / sin(wsl);
        mpreal Wsu = (c0 - cos(wsu)) / sin(wsu);
        mpreal Ws = std::min(abs(Wsl), abs(Wsu));
        mpreal Wp = tan((wpu - wpl) / 2_mpr);

//        mpreal c0=sin(wsl+wsu)/(sin(wsl)+sin(wsu));
//        mpreal Wpl=(c0-cos(wpl))/sin(wpl);
//        mpreal Wpu=(c0-cos(wpu))/sin(wpu);
//        mpreal Ws=tan((wsu-wsl)/2_mpr);
//        mpreal Wp=std::max(abs(Wpl),abs(Wpu));
        //设计模拟低通滤波器
        auto [z0, p0, H0, B0, A0] = AF::ellipitic_filter(Wp, Ws, std::move(Ap), std::move(As));
        //反转换

    } else if (type == bandstop) {
        //转换技术指标//TODO:还有一种技术指标的转换方式
        mpreal c0 = sin(wpl + wpu) / (sin(wpl) + sin(wpu));
        mpreal Wsl = sin(wsl) / (c0 - cos(wsl));
        mpreal Wsu = sin(wsu) / (c0 - cos(wsu));
        mpreal Ws = std::min(abs(Wsl), abs(Wsu));
        mpreal Wp = cot((wpu - wpl) / 2_mpr);
        //设计模拟低通滤波器
        auto [z0, p0, H0, B0, A0] = AF::ellipitic_filter(Wp, Ws, std::move(Ap), std::move(As));
        //反转换

    }
    return {};
}
