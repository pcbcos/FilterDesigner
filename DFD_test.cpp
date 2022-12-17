//
// Created by wengz on 2022/12/6.
//
#include <iostream>
#include <gmpxx.h>
#include <mpreal.h>
#include "TransferFunction.h"
#include "mpreal_ex.h"
#include "FilterDesign.h"

using mpfr::mpreal;
using mpcomplex = std::complex<mpreal>;

auto main() -> int {
    mpreal::set_default_prec(mpfr::digits2bits(50));
    mpreal k;

    std::cout.precision(20);
    std::cout << "Digital_Filter_Designer" << '\n';
    int ftype;
    std::cout << "Filter type(0:butterworth\t1:chebyshev1\t2:chebyshev2\t3:elliptic):";
    std::cin >> ftype;
    std::cout << '\n';
    std::cout << "band_type(0:lowpass\t1:highpass\t2:bandpass\t3:bandstop):";
    int type;
    std::cin >> type;
    mpreal wp, wpl, wpu;
    mpreal ws, wsl, wsu;

    if ((filter_band_type(type) == lowpass) || (filter_band_type(type) == highpass)) {
        std::cout << "wp=";
        std::cin >> wp;

        std::cout << "ws=";
        std::cin >> ws;

    } else {
        std::cout << "wpl=";
        std::cin >> wpl;

        std::cout << "wpu=";
        std::cin >> wpu;

        std::cout << "wsl=";
        std::cin >> wsl;

        std::cout << "wsu=";
        std::cin >> wsu;

    }

    std::cout << "Ap(dB)=";
    mpreal Ap;
    std::cin >> Ap;

    std::cout << "As(dB)=";
    mpreal As;
    std::cin >> As;

    zeros z;
    poles p;
    mpreal H0;
    std::vector<std::array<mpfr::mpreal, 3>> B;
    std::vector<std::array<mpfr::mpreal, 3>> A;
    if ((filter_band_type(type) == lowpass) || (filter_band_type(type) == highpass)) {
        if (ftype == 0) {
            auto [z1, p1, H01, B1, A1] = DF::butterworth_filter(wp, ws, Ap, As, (filter_band_type) type);
            z = z1;
            p = p1;
            H0 = H01;
            B = std::move(B1);
            A = std::move(A1);
        } else if (ftype == 1) {
            auto [z1, p1, H01, B1, A1] = DF::chebyshev1_filter(wp, ws, Ap, As, (filter_band_type) type);
            z = z1;
            p = p1;
            H0 = H01;
            B = std::move(B1);
            A = std::move(A1);
        } else if (ftype == 2) {
            auto [z1, p1, H01, B1, A1] = DF::chebyshev2_filter(wp, ws, Ap, As, (filter_band_type) type);
            z = z1;
            p = p1;
            H0 = H01;
            B = std::move(B1);
            A = std::move(A1);
        } else {
            auto [z1, p1, H01, B1, A1] = DF::ellipitic_filter(wp, ws, Ap, As, (filter_band_type) type);
            z = z1;
            p = p1;
            H0 = H01;
            B = std::move(B1);
            A = std::move(A1);
        }
    } else {
        if (ftype == 0) {
            auto [z1, p1, H01, B1, A1] = DF::butterworth_filter(wpu, wpl, wsu, wsl, Ap, As, (filter_band_type) type);
            z = z1;
            p = p1;
            H0 = H01;
            B = std::move(B1);
            A = std::move(A1);
        } else if (ftype == 1) {
            auto [z1, p1, H01, B1, A1] = DF::chebyshev1_filter(wpu, wpl, wsu, wsl, Ap, As, (filter_band_type) type);
            z = z1;
            p = p1;
            H0 = H01;
            B = std::move(B1);
            A = std::move(A1);
        } else if (ftype == 2) {
            auto [z1, p1, H01, B1, A1] = DF::chebyshev2_filter(wpu, wpl, wsu, wsl, Ap, As, (filter_band_type) type);
            z = z1;
            p = p1;
            H0 = H01;
            B = std::move(B1);
            A = std::move(A1);
        } else {
            auto [z1, p1, H01, B1, A1] = DF::ellipitic_filter(wpu, wpl, wsu, wsl, Ap, As, (filter_band_type) type);
            z = z1;
            p = p1;
            H0 = H01;
            B = std::move(B1);
            A = std::move(A1);
        }
    }

    std::cout << "\n零点z" << '\n';
    for (auto &zi: z) {
        std::cout << zi << '\n';
    }
    std::cout << "\n极点p" << '\n';
    for (auto &pi: p) {
        std::cout << pi << '\n';
    }
    std::cout << "\n直流增益H0" << '\n';
    std::cout << H0 << '\n';
    std::cout << "H(s)" << '\n';
    for (int i = 0; i < B.size(); i++) {
        std::cout << TransferFunction(B[i], A[i]).toLaTeX(TransferFunction::DF) << '\n';
    }

    std::vector<TransferFunction> v;
    for (int i = 0; i < B.size(); i++) {
        std::cout << TransferFunction(B[i], A[i]).toMMA() << '*';
        v.emplace_back(B[i], A[i]);
    }
    for(int i=0;i<=314;i++){
        mpreal w0=i/100_mpr;
        mpreal H=1_mpr;
        for(auto& tf:v){
            H*= tf.gain(w0, TransferFunction::DF);
        }
        if(!mpfr::isnan(H)){
            std::cout<<"w0="<<w0<<"  Gain="<<H<<'\n';
        }

    }

    return 0;
}