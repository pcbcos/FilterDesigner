//
// Created by wengz on 2022/11/11.
//

#include <vector>
#include "EllipticFunction.h"

namespace Ellipitic {
    auto __EllipticK(mpreal &__k) {
        //const mpreal kmin("1e-10");
        const mpreal kmin = std::numeric_limits<mpreal>::min();
        auto kmax = mpreal(sqrt(1 - kmin * kmin));
        if (__k > 0 && __k < kmax) {
            mpreal K = mpfr::const_pi() / 2;
            mpreal k_prime = sqrt(1 - __k * __k);
            while (__k > mpfr::machine_epsilon()) {
                __k = (__k / (1 + k_prime)) * (__k / (1 + k_prime));  //k_{n+1}=...
                k_prime = sqrt(1 - __k * __k);                      //k'_{n+1}=\sqrt{1-k^2_{n+1}}
                K *= (1 + __k);
            }
            return K;
        } else {
            mpreal k_prime = sqrt(1 - __k * __k);
            mpreal L = (-log(mpreal(k_prime) / 4));
            return L + (L - 1) * k_prime * k_prime / 2;
        }
    }

    auto EllipticK(const mpreal &k) -> std::pair<mpreal, mpreal> {
        mpreal __k = k;
        mpreal __k_prime = sqrt(1 - k * k);
        return std::make_pair(__EllipticK(__k), __EllipticK(__k_prime));
    }

    auto cd(const mpcomplex &u, mpreal k)->mpcomplex {
        /*calc w=cd(uK,k)  u is given
     * */
        std::vector<mpreal> kn;
        kn.push_back(k);
        const mpreal kmin = std::numeric_limits<mpreal>::min();
        auto kmax = mpreal(sqrt(1 - kmin * kmin));
        if (k > 0 && k < kmax) {
            mpreal K = mpfr::const_pi() / 2;
            mpreal k_prime = sqrt(1 - k * k);
            while (k > mpfr::machine_epsilon()) {
                k = (k / (1 + k_prime)) * (k / (1 + k_prime));  //k_{n+1}=...
                kn.push_back(k);
                k_prime = sqrt(1 - k * k);                      //k'_{n+1}=\sqrt{1-k^2_{n+1}}
                K *= (1 + k);
            }
        }
        mpcomplex w = cos(u * mpfr::const_pi() / mpreal(2));
        //std::cout<<w<<'\n';
        for (int m = kn.size() - 1; m >= 1; m--) {
            w = (1 + kn[m]) / (mpreal(1) / w + kn[m] * w);
        }
        return w;
    }

    auto sn(const mpcomplex &u, mpreal k)->mpcomplex {
        /*calc w=cd(uK,k)  u is given
     * */
        std::vector<mpreal> kn;
        kn.push_back(k);
        const mpreal kmin = std::numeric_limits<mpreal>::min();
        auto kmax = mpreal(sqrt(1 - kmin * kmin));
        if (k > 0 && k < kmax) {
            mpreal K = mpfr::const_pi() / 2;
            mpreal k_prime = sqrt(1 - k * k);
            while (k > mpfr::machine_epsilon()) {
                k = (k / (1 + k_prime)) * (k / (1 + k_prime));  //k_{n+1}=...
                kn.push_back(k);
                k_prime = sqrt(1 - k * k);                      //k'_{n+1}=\sqrt{1-k^2_{n+1}}
                K *= (1 + k);
            }
        }
        mpcomplex w = sin(u * mpfr::const_pi() / 2_mpr);
        for (int m = kn.size() - 1; m >= 1; m--) {
            w = (1 + kn[m]) / (mpreal(1) / w + kn[m] * w);
        }
        return w;
    }

    auto acd(mpcomplex w, mpreal k) ->mpcomplex{
        std::vector<mpreal> kn;
        kn.push_back(k);
        const mpreal kmin = std::numeric_limits<mpreal>::min();
        auto kmax = mpreal(sqrt(1 - kmin * kmin));
        if (k > 0 && k < kmax) {
            mpreal K = mpfr::const_pi() / 2;
            mpreal k_prime = sqrt(1 - k * k);
            while (k > mpfr::machine_epsilon()) {
                k = (k / (1 + k_prime)) * (k / (1 + k_prime));  //k_{n+1}=...
                kn.push_back(k);
                k_prime = sqrt(1 - k * k);                      //k'_{n+1}=\sqrt{1-k^2_{n+1}}
                K *= (1 + k);
            }
        }
        for (int n = 1; n <= kn.size() - 1; n++) {
            w = 2_mpr * w / (1 + kn[n]) / (1_mpr + sqrt(1_mpr - kn[n - 1] * kn[n - 1] * w * w));
        }
        mpcomplex u = 2 / mpfr::const_pi() * acos(w);
        while (u.real() > 2) {
            u -= 2;
        }
        while (u.real() < 0) {
            u += 2;
        }
        auto [K, Kp] = EllipticK(k);
        mpreal temp = Kp / K;
        while (u.imag() > temp) {
            u -= mpcomplex(0, 1) * temp;
        }
        while (u.imag() < -temp) {
            u += mpcomplex(0, 1) * temp;
        }
        return u;
    }

    auto asn(mpcomplex w, mpreal k)->mpcomplex {
        std::vector<mpreal> kn;
        kn.push_back(k);
        const mpreal kmin = std::numeric_limits<mpreal>::min();
        auto kmax = mpreal(sqrt(1 - kmin * kmin));
        if (k > 0 && k < kmax) {
            mpreal K = mpfr::const_pi() / 2;
            mpreal k_prime = sqrt(1 - k * k);
            while (k > mpfr::machine_epsilon()) {
                k = (k / (1 + k_prime)) * (k / (1 + k_prime));  //k_{n+1}=...
                kn.push_back(k);
                k_prime = sqrt(1 - k * k);                      //k'_{n+1}=\sqrt{1-k^2_{n+1}}
                K *= (1 + k);
            }
        }
        for (int n = 1; n <= kn.size() - 1; n++) {
            w = 2_mpr * w / (1 + kn[n]) / (1_mpr + sqrt(1_mpr - kn[n - 1] * kn[n - 1] * w * w));
        }
        mpcomplex u = 2 / mpfr::const_pi() * asin(w);
        while (u.real() > 1) {
            u -= 2;
        }
        while (u.real() < -1) {
            u += 2;
        }
        auto [K, Kp] = EllipticK(k);
        mpreal temp = Kp / K;
        while (u.imag() > temp) {
            u -= mpcomplex(0, 1) * temp;
        }
        while (u.imag() < -temp) {
            u += mpcomplex(0, 1) * temp;
        }
        return u;
    }

    mpreal deg1(uint32_t N, mpreal k) {
        uint32_t L = N / 2;
        uint32_t r = N % 2;
        mpreal k1 = pow(k, N);
        for (int i = 1; i <= L; i++) {
            mpreal ui = (2_mpr * i - 1_mpr) / N;
            k1 *= pow(sn(ui, k).real(), 4);
        }
        return k1;
    }

    mpreal deg(uint32_t N, mpreal k1) {
        uint32_t L = N / 2;
        uint32_t r = N % 2;
        mpreal kc = mpfr::sqrt(1 - k1 * k1);
        mpreal kp = pow(kc, N);
        for (int i = 1; i <= L; i++) {
            mpreal ui = (2_mpr * i - 1_mpr) / N;
            kp *= pow(sn(ui, kc).real(), 4);//supposed to be a real nuber
        }
        return sqrt(1 - kp * kp);
    }

}

