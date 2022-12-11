//
// Created by wengz on 2022/11/11.
//

#ifndef ELLIPTICFILTER_MPREAL_EX_H
#define ELLIPTICFILTER_MPREAL_EX_H

#include <mpreal.h>


static mpfr::mpreal operator ""_mpr(const char *s) {
    return {s};
}

using mpcomplex = std::complex<mpfr::mpreal>;

static std::ostream &operator<<(std::ostream &os, const mpcomplex &v) {
    os << v.real();
    if (v.imag() >= 0) {
        os << "+j" << v.imag();
    } else {
        os << "-j" << -v.imag();
    }
    return os;
}

#ifdef _MSC_VER
template <>
class std::_Ctraits<mpfr::mpreal> {
    using _Ty = mpfr::mpreal;
public:
    static _Ty _Flt_eps() { // get epsilon
        //return numeric_limits<_Ty>::epsilon();
        return mpfr::machine_epsilon();
    }

    static _Ty _Flt_max() {
        //return (numeric_limits<_Ty>::max)();
        return mpfr::maxval();
    }

    static _Ty _Flt_norm_min() {
        return mpfr::minval() > 0 ? mpfr::minval() : 0;
    }

    static _Ty _Abs(_Ty _Left) {
        return static_cast<_Ty>(_Signbit(_Left) ? -_Left : _Left);
    }

    static _Ty _Cosh(_Ty _Left, _Ty _Right) { // return cosh(_Left) * _Right
        return static_cast<_Ty>(_CSTD _Cosh(static_cast<double>(_Left), static_cast<double>(_Right)));
    }

    static _Ty _Copysign(_Ty _Magnitude, _Ty _Sign) {
        return static_cast<_Ty>(_Signbit(_Sign) ? -_Abs(_Magnitude) : _Abs(_Magnitude));
    }

    static short _Exp(_Ty* _Pleft, _Ty _Right, short _Exponent) { // compute exp(*_Pleft) * _Right * 2 ^ _Exponent
        double _Tmp = static_cast<double>(*_Pleft);
        short _Ans = _CSTD _Exp(&_Tmp, static_cast<double>(_Right), _Exponent);
        *_Pleft = static_cast<_Ty>(_Tmp);
        return _Ans;
    }

    static _Ty _Infv() { // return infinity
        return numeric_limits<_Ty>::infinity();
    }

    static bool _Isinf(_Ty _Left) { // test for infinity
        const auto _Tmp = static_cast<double>(_Left);
        const auto _Uint = _Bit_cast<uint64_t>(_Tmp);
        return (_Uint & 0x7fffffffffffffffU) == 0x7ff0000000000000U;
    }

    static bool _Isnan(_Ty _Left) {
        const auto _Tmp = static_cast<double>(_Left);
        const auto _Uint = _Bit_cast<uint64_t>(_Tmp);
        return (_Uint & 0x7fffffffffffffffU) > 0x7ff0000000000000U;
    }

    static _Ty _Nanv() { // return NaN
        return numeric_limits<_Ty>::quiet_NaN();
    }

    static bool _Signbit(_Ty _Left) {
        return (_STD signbit)(static_cast<double>(_Left));
    }

    static _Ty _Sinh(_Ty _Left, _Ty _Right) { // return sinh(_Left) * _Right
        return static_cast<_Ty>(_CSTD _Sinh(static_cast<double>(_Left), static_cast<double>(_Right)));
    }

    static _Ty asinh(_Ty _Left) {
        if (_Left == 0 || _Isnan(_Left) || _Isinf(_Left)) {
            return _Left;
        }

        _Ty _Ln2 = 0.69314718055994530941723212145817658_mpr;

        const _Ty _Old_left = _Left;
        _Ty _Ans;

        _Left = _Abs(_Left);

        if (_Left < 2 / _Flt_eps()) {
            _Ans = log1p(_Left + _Left * _Left / (1 + sqrt(_Left * _Left + 1)));
        }
        else {
            _Ans = log(_Left) + _Ln2;
        }

        return _Copysign(_Ans, _Old_left);
    }

    static _Ty atan2(_Ty _Yval, _Ty _Xval) { // return atan(_Yval / _Xval)
        return static_cast<_Ty>(_CSTD atan2(static_cast<double>(_Yval), static_cast<double>(_Xval)));
    }

    static _Ty cos(_Ty _Left) {
        return static_cast<_Ty>(_CSTD cos(static_cast<double>(_Left)));
    }

    static _Ty exp(_Ty _Left) {
        return static_cast<_Ty>(_CSTD exp(static_cast<double>(_Left)));
    }

    static _Ty ldexp(_Ty _Left, int _Exponent) { // return _Left * 2 ^ _Exponent
        return static_cast<_Ty>(_CSTD ldexp(static_cast<double>(_Left), _Exponent));
    }

    static _Ty log(_Ty _Left) {
        return static_cast<_Ty>(_CSTD log(static_cast<double>(_Left)));
    }

    static _Ty log1p(_Ty _Left) { // return log(1 + _Left)
        if (_Left < -1) {
            return _Nanv();
        }
        else if (_Left == 0) {
            return _Left;
        }
        else { // compute log(1 + _Left) with fixup for small _Left
            _Ty _Leftp1 = 1 + _Left;
            return log(_Leftp1) - ((_Leftp1 - 1) - _Left) / _Leftp1;
        }
    }

    static _Ty pow(_Ty _Left, _Ty _Right) {
        return static_cast<_Ty>(_CSTD pow(static_cast<double>(_Left), static_cast<double>(_Right)));
    }

    static _Ty sin(_Ty _Left) {
        return static_cast<_Ty>(_CSTD sin(static_cast<double>(_Left)));
    }

    static _Ty sqrt(_Ty _Left) {
        return static_cast<_Ty>(_CSTD sqrt(static_cast<double>(_Left)));
    }

    static _Ty tan(_Ty _Left) {
        return static_cast<_Ty>(_CSTD tan(static_cast<double>(_Left)));
    }

    static _Ty hypot(_Ty _Left, _Ty _Right) {
        return static_cast<_Ty>(_CSTD hypot(static_cast<double>(_Left), static_cast<double>(_Right)));
    }
};

template <>
mpfr::mpreal std::_Fabs<mpfr::mpreal>(const complex<mpfr::mpreal>& _Left, int* _Pexp) { // Used by sqrt(), return magnitude and scale factor.
                                                   // Returns a non-zero even integer in *_Pexp when _Left is finite
                                                   // and non-zero.
                                                   // Returns 0 in *_Pexp when _Left is zero, infinity, or NaN.
    using _Ty = mpfr::mpreal;
    *_Pexp = 0;
    _Ty _Av = _Ctraits<_Ty>::_Abs(_STD real(_Left));
    _Ty _Bv = _Ctraits<_Ty>::_Abs(_STD imag(_Left));

    if (_Ctraits<_Ty>::_Isinf(_Av) || _Ctraits<_Ty>::_Isinf(_Bv)) {
        return _Ctraits<_Ty>::_Infv(); // at least one component is INF
    }
    else if (_Ctraits<_Ty>::_Isnan(_Av)) {
        return _Av; // real component is NaN
    }
    else if (_Ctraits<_Ty>::_Isnan(_Bv)) {
        return _Bv; // imaginary component is NaN
    }
    else { // neither component is NaN or INF
        if (_Av < _Bv) { // ensure that |_Bv| <= |_Av|
            _STD swap(_Av, _Bv);
        }

        if (_Av == 0) {
            return _Av; // |0| == 0
        }

        if (1 <= _Av) {
            *_Pexp = 4;
            //_Av = _Av * static_cast<_Ty>(0.0625);
            //_Bv = _Bv * static_cast<_Ty>(0.0625);
            _Av = _Av * 0.0625_mpr;
            _Bv = _Bv * 0.0625_mpr;
        }
        else {
            _Ty _Flt_eps = _Ctraits<_Ty>::_Flt_eps();
            // TRANSITION, workaround for non-floating-point _Ty
            _Ty _Leg_tiny = _Flt_eps == 0 ? _Ty{ 0 } : 2 * _Ctraits<_Ty>::_Flt_norm_min() / _Flt_eps;

            if (_Av < _Leg_tiny) {
                int _Exponent = -2 * numeric_limits<_Ty>::digits(_Av);
                //int _Exponent = -2 * mpfr::mpreal::digits(_Av);
                *_Pexp = _Exponent;
                _Av = _Ctraits<_Ty>::ldexp(_Av, -_Exponent);
                _Bv = _Ctraits<_Ty>::ldexp(_Bv, -_Exponent);
            }
            else {
                *_Pexp = -2;
                _Av = _Av * 4;
                _Bv = _Bv * 4;
            }
        }

        const _Ty _Tmp = _Av - _Bv;
        if (_Tmp == _Av) {
            return _Av; // _Bv unimportant
        }
        else if (_Bv < _Tmp) { // use simple approximation
            const _Ty _Qv = _Av / _Bv;
            return _Av + _Bv / (_Qv + _Ctraits<_Ty>::sqrt(_Qv * _Qv + 1));
        }
        else { // use 1 1/2 precision to preserve bits
            _Ty _Root2 = 1.4142135623730950488016887242096981_mpr;
            _Ty _Oneplusroot2high = 10125945.0_mpr / 4194304.0_mpr; // exact if prec >= 24 bits
            _Ty _Oneplusroot2low = 1.4341252375973918872420969807856967e-7_mpr;

            const _Ty _Qv = _Tmp / _Bv;
            const _Ty _Rv = (_Qv + 2) * _Qv;
            const _Ty _Sv = _Rv / (_Root2 + _Ctraits<_Ty>::sqrt(_Rv + 2)) + _Oneplusroot2low + _Qv + _Oneplusroot2high;
            return _Av + _Bv / _Sv;
        }
    }
}

template <>
_NODISCARD std::complex<mpfr::mpreal> std::acos<mpfr::mpreal>(const std::complex<mpfr::mpreal>& _Left) {
    using _Ty = mpfr::mpreal;
    const _Ty _Arcbig = static_cast<_Ty>(0.25) * _Ctraits<_Ty>::sqrt(_Ctraits<_Ty>::_Flt_max());
    _Ty _Pi = 3.1415926535897932384626433832795029_mpr;

    const _Ty _Re = real(_Left);
    const _Ty _Im = imag(_Left);
    _Ty _Ux;
    _Ty _Vx;

    if (_Ctraits<_Ty>::_Isnan(_Re) || _Ctraits<_Ty>::_Isnan(_Im)) { // at least one NaN
        _Ux = _Ctraits<_Ty>::_Nanv();
        _Vx = _Ux;
    }
    else if (_Ctraits<_Ty>::_Isinf(_Re)) { // (+/-Inf, not NaN)
        if (_Ctraits<_Ty>::_Isinf(_Im)) {
            if (_Re < 0) {
                _Ux = static_cast<_Ty>(0.75) * _Pi; // (-Inf, +/-Inf)
            }
            else {
                _Ux = static_cast<_Ty>(0.25) * _Pi; // (+Inf, +/-Inf)
            }
        }
        else if (_Re < 0) {
            _Ux = _Pi; // (-Inf, finite)
        }
        else {
            _Ux = 0; // (+Inf, finite)
        }
        _Vx = -_Ctraits<_Ty>::_Copysign(_Ctraits<_Ty>::_Infv(), _Im);
    }
    else if (_Ctraits<_Ty>::_Isinf(_Im)) { // (finite, finite)
        _Ux = static_cast<_Ty>(0.50) * _Pi; // (finite, +/-Inf)
        _Vx = -_Im;
    }
    else { // (finite, finite)
        const complex<_Ty> _Wx = sqrt(complex<_Ty>(1 + _Re, -_Im));
        const complex<_Ty> _Zx = sqrt(complex<_Ty>(1 - _Re, -_Im));
        const _Ty _Wr = real(_Wx);
        const _Ty _Wi = imag(_Wx);
        const _Ty _Zr = real(_Zx);
        const _Ty _Zi = imag(_Zx);
        _Ty _Alfa;
        _Ty _Beta;

        _Ux = 2 * _Ctraits<_Ty>::atan2(_Zr, _Wr);

        if (_Arcbig < _Wr) { // real parts large
            _Alfa = _Wr;
            _Beta = _Zi + _Wi * (_Zr / _Alfa);
        }
        else if (_Arcbig < _Wi) { // imag parts large
            _Alfa = _Wi;
            _Beta = _Wr * (_Zi / _Alfa) + _Zr;
        }
        else if (_Wi < -_Arcbig) { // imag part of w large negative
            _Alfa = -_Wi;
            _Beta = _Wr * (_Zi / _Alfa) - _Zr;
        }
        else { // shouldn't overflow
            _Alfa = 0;
            _Beta = _Wr * _Zi + _Wi * _Zr; // Im(w * z)
        }

        _Vx = _Ctraits<_Ty>::asinh(_Beta);
        if (_Alfa != 0) {
            if (0 <= _Vx) {
                _Vx += _Ctraits<_Ty>::log(_Alfa);
            }
            else {
                _Vx -= _Ctraits<_Ty>::log(_Alfa); // asinh(a*b) = asinh(a)+log(b)
            }
        }
    }
    return complex<_Ty>(_Ux, _Vx);
};

template <>
_NODISCARD std::complex<mpfr::mpreal> std::asinh(const std::complex<mpfr::mpreal>& _Left) {
    using _Ty = mpfr::mpreal;
    const _Ty _Arcbig = static_cast<_Ty>(0.25) * _Ctraits<_Ty>::sqrt(_Ctraits<_Ty>::_Flt_max());
    _Ty _Pi = 3.1415926535897932384626433832795029_mpr;

    const _Ty _Re = real(_Left);
    _Ty _Im = imag(_Left);
    _Ty _Ux;
    _Ty _Vx;

    if (_Ctraits<_Ty>::_Isnan(_Re) || _Ctraits<_Ty>::_Isnan(_Im)) { // at least one NaN/Inf
        _Ux = _Ctraits<_Ty>::_Nanv();
        _Vx = _Ux;
    }
    else if (_Ctraits<_Ty>::_Isinf(_Re)) { // (+/-Inf, not NaN)
        _Ux = _Ctraits<_Ty>::_Infv();

        if (_Ctraits<_Ty>::_Isinf(_Im)) { // (+/-Inf, +/-Inf)
            _Ux = _Re;
            _Vx = _Ctraits<_Ty>::_Copysign(static_cast<_Ty>(0.25) * _Pi, _Im);
        }
        else { // (+/-Inf, finite)
            _Ux = _Re;
            _Vx = _Ctraits<_Ty>::_Copysign(_Ty{ 0 }, _Im);
        }
    }
    else if (_Ctraits<_Ty>::_Isinf(_Im)) { // (finite, +/-Inf)
        _Ux = _Ctraits<_Ty>::_Copysign(_Ctraits<_Ty>::_Infv(), _Re);
        _Vx = _Ctraits<_Ty>::_Copysign(static_cast<_Ty>(0.50) * _Pi, _Im);
    }
    else { // (finite, finite)
        const complex<_Ty> _Wx = sqrt(complex<_Ty>(1 - _Im, _Re));
        const complex<_Ty> _Zx = sqrt(complex<_Ty>(1 + _Im, -_Re));
        const _Ty _Wr = real(_Wx);
        const _Ty _Wi = imag(_Wx);
        const _Ty _Zr = real(_Zx);
        const _Ty _Zi = imag(_Zx);
        _Ty _Alfa;
        _Ty _Beta;

        if (_Arcbig < _Wr) { // real parts large
            _Alfa = _Wr;
            _Beta = _Wi * (_Zr / _Alfa) - _Zi;
        }
        else if (_Arcbig < _Wi) { // imag parts large
            _Alfa = _Wi;
            _Beta = _Zr - _Wr * (_Zi / _Alfa);
        }
        else if (_Wi < -_Arcbig) { // imag part of w large negative
            _Alfa = -_Wi;
            _Beta = -_Zr - _Wr * (_Zi / _Alfa);
        }
        else { // shouldn't overflow
            _Alfa = 0;
            _Beta = _Wi * _Zr - _Wr * _Zi; // Im(w * conj(z))
        }

        _Ux = _Ctraits<_Ty>::asinh(_Beta);
        if (_Alfa != 0) {
            if (0 <= _Ux) {
                _Ux += _Ctraits<_Ty>::log(_Alfa);
            }
            else {
                _Ux -= _Ctraits<_Ty>::log(_Alfa); // asinh(a*b) = asinh(a)+log(b)
            }
        }

        _Vx = _Ctraits<_Ty>::atan2(_Im, real(_Wx * _Zx));
    }

    return complex<_Ty>(_Ux, _Vx);
}

template <>
_NODISCARD std::complex<mpfr::mpreal> std::asin<mpfr::mpreal>(const std::complex<mpfr::mpreal>& _Left) {
    using _Ty = mpfr::mpreal;
    complex<_Ty> _Asinh = _STD asinh(complex<_Ty>(-_Left.imag(),_Left.real()));

    return complex<_Ty>(imag(_Asinh), -real(_Asinh));
}
#endif

#endif //ELLIPTICFILTER_MPREAL_EX_H
