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
    std::cout << "elliptic_analog_design" << '\n';

    std::cout << "type:";
    int type;
    std::cin >> type;
    mpreal wp, wpl, wpu;
    mpreal ws, wsl, wsu;

    if ((filter_band_type(type) == lowpass) || (filter_band_type(type) == highpass)) {
        std::cout << "Wp=";
        std::cin >> wp;

        std::cout << "Ws=";
        std::cin >> ws;

    } else {
        std::cout << "Wpl=";
        std::cin >> wpl;

        std::cout << "Wpu=";
        std::cin >> wpu;

        std::cout << "Wsl=";
        std::cin >> wsl;

        std::cout << "Wsu=";
        std::cin >> wsu;

    }

    std::cout << "Ap(dB)=";
    mpreal Ap;
    std::cin >> Ap;

    std::cout << "As(dB)=";
    mpreal As;
    std::cin >> As;

    auto [z, p, H0, B, A] = AF::ellipitic_filter(wp, ws, Ap, As, (filter_band_type) type);
    //auto [z, p, H0, B, A] = AF::detail::elliptic_lp_prototype(N, Ap, As);
    if ((filter_band_type(type) == lowpass) || (filter_band_type(type) == highpass)) {
        //auto [z1, p1, H01, B1, A1] = AF::ellipitic_filter(wp, ws, Ap, As, (filter_band_type) type);
        auto [z1, p1, H01, B1, A1] = AF::chebyshev2_filter(wp, ws, Ap, As, (filter_band_type) type);
        z = z1;
        p = p1;
        H0 = H01;
        B = std::move(B1);
        A = std::move(A1);
    } else {
        auto [z1, p1, H01, B1, A1] = AF::chebyshev2_filter(wpu, wpl, wsu, wsl, Ap, As, (filter_band_type) type);
        z = z1;
        p = p1;
        H0 = H01;
        B = std::move(B1);
        A = std::move(A1);
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
        std::cout << TransferFunction(B[i], A[i]).toLaTex() << '\n';
    }

    for (int i = 0; i < B.size(); i++) {
        std::cout << B[i][0] << ',' << B[i][1] << ',' << B[i][2] << ',' << A[i][0] << ',' << A[i][1] << ',' << A[i][2]
                  << '\n';
    }
    for (int i = 0; i < B.size(); i++) {
        std::cout << TransferFunction(B[i], A[i]).toMMA() << '*';
    }

    return 0;
}