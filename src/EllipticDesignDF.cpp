//
// Created by wengz on 2022/12/5.
//
#include <vector>
#include "EllipticFunction.h"
#include "FilterDesign.h"

auto DF::ellipitic_filter(mpfr::mpreal wp, mpfr::mpreal ws, mpfr::mpreal Ap, mpfr::mpreal As,
                          filter_band_type type) -> std::tuple<std::vector<mpcomplex>, std::vector<mpcomplex>, mpreal, std::vector<std::array<mpreal, 3>>, std::vector<std::array<mpreal, 3>>> {
    return {};
}

auto DF::ellipitic_filter(mpfr::mpreal wpu, mpfr::mpreal wpl, mpfr::mpreal wsu, mpfr::mpreal wsl, mpfr::mpreal Ap,
                          mpfr::mpreal As,
                          filter_band_type type) -> std::tuple<std::vector<mpcomplex>, std::vector<mpcomplex>, mpreal, std::vector<std::array<mpreal, 3>>, std::vector<std::array<mpreal, 3>>> {
    return {};
}
