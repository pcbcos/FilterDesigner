//
// Created by wengz on 2022/11/10.
//

#ifndef ELLIPTICFILTER_TRANSFERFUNCTION_H
#define ELLIPTICFILTER_TRANSFERFUNCTION_H

#include <utility>
#include <vector>
#include <mpreal.h>
#include <array>

using mpfr::mpreal;
using mpcomplex = std::complex<mpreal>;

class TransferFunction {
public:
    enum ftype{
        AF,
        DF
    };
    TransferFunction() = default;

    TransferFunction(std::array<mpreal, 3> B, std::array<mpreal, 3> A) : B(std::move(B)), A(std::move(A)) {

    };

    inline auto get_B() {
        return B;
    }

    inline auto get_A() {
        return A;
    }

    inline std::string toLaTex(ftype type=AF) {
        std::string s;
        auto get_flag = [](const mpreal &x) {
            return x > 0 ? std::string("+") : std::string("-");
        };

        if(type==AF){
            s += "\\frac{" + B[0].toString(5) + get_flag(B[1]) + abs(B[1]).toString(5) + "s" +
                 get_flag(B[2]) +
                 abs(B[2]).toString(5) + "s^2}{";
            s += A[0].toString(5) + get_flag(A[1]) + abs(A[1]).toString(5) + "s" + get_flag(A[2]) + abs(A[2]).toString(5) +
                 "s^2}";
        }else{
            s += "\\frac{" + B[0].toString(5) + get_flag(B[1]) + abs(B[1]).toString(5) + "z^{-1}" +
                 get_flag(B[2]) +
                 abs(B[2]).toString(5) + "z^{-2}}{";
            s += A[0].toString(5) + get_flag(A[1]) + abs(A[1]).toString(5) + "z^{-1}" + get_flag(A[2]) + abs(A[2]).toString(5) +
                 "z^{-2}}";
        }

        return s;
    }

private:
    std::array<mpreal, 3> B;
    std::array<mpreal, 3> A;
};

#endif //ELLIPTICFILTER_TRANSFERFUNCTION_H
