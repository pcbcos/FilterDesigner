//
// Created by wengz on 2022/11/11.
//

#ifndef ELLIPTICFILTER_ELLIPTICFUNCTION_H
#define ELLIPTICFILTER_ELLIPTICFUNCTION_H

#include "mpreal_ex.h"

namespace Ellipitic {
    using mpfr::mpreal;
    using mpcomplex = std::complex<mpreal>;

    auto __EllipticK(mpreal &__k);

    auto EllipticK(const mpreal &k)->std::pair<mpreal, mpreal>;

    auto cd(const mpcomplex &u, mpreal k)->mpcomplex;

    auto sn(const mpcomplex &u, mpreal k)->mpcomplex;

    auto acd(mpcomplex w, mpreal k)->mpcomplex;

    auto asn(mpcomplex w, mpreal k)->mpcomplex;

    mpreal deg1(uint32_t N, mpreal k);

    mpreal deg(uint32_t N, mpreal k1);

}

#endif //ELLIPTICFILTER_ELLIPTICFUNCTION_H
