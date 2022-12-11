//
// Created by wengz on 2022/11/11.
//

#include "FilterDesign.h"
using mpfr::mpreal;
auto AF::detail::elliptic_lp_prototype(uint32_t N, const mpfr::mpreal &Ap,
                                       const mpfr::mpreal &As) -> design_res {
    //N:阶数
    //Ap:通带最大衰减
    //As:阻带最小衰减
    mpreal Gp = pow(10_mpr, -Ap / 20_mpr);
    mpreal ep = sqrt(pow(10_mpr, Ap / 10_mpr) - 1);//通带波纹系数
    mpreal es = sqrt(pow(10_mpr, As / 10_mpr) - 1);//阻带波纹系数
    mpreal k1 = ep / es;
    mpreal k = Ellipitic::deg(N, k1);

    uint32_t L = N / 2;
    uint32_t r = N % 2;

    std::vector<std::array<mpreal, 3>> B;//分子
    std::vector<std::array<mpreal, 3>> A;//分母
    auto v0 = -mpcomplex(0, 1) * Ellipitic::asn(mpcomplex(0, 1) / ep, k1) / mpreal(N);

    std::vector<mpcomplex> z;
    std::vector<mpcomplex> p;
    if (r == 0) {
        B.push_back({Gp, 0, 0});
        A.push_back({1, 0, 0});
    } else {
        auto p0 = mpcomplex(0, 1) * Ellipitic::sn(mpcomplex(0, 1) * v0, k);
        p.push_back(p0);
        B.push_back({1, 0, 0});
        A.push_back({1, -(1_mpr / p0).real(), 0});
    }

    for (uint32_t i = 1; i <= L; i++) {
        auto ui = (2_mpr * i - 1_mpr) / N;
        auto zeta_i = Ellipitic::cd(ui, k);
        auto zi = mpcomplex(0, 1) / (k * zeta_i);
        auto pi = mpcomplex(0, 1) * Ellipitic::cd(ui - mpcomplex(0, 1) * v0, k);

        std::cout << k << '\n';
        std::cout << v0 << '\n';

        p.push_back(pi);
        p.push_back(conj(pi));

        z.push_back(zi);
        z.push_back(conj(zi));

        B.push_back({1, -2 * (1_mpr / zi).real(), 1_mpr / abs(zi) / abs(zi)});
        A.push_back({1, -2 * (1_mpr / pi).real(), 1_mpr / abs(pi) / abs(pi)});
    }

    auto H0 = r ? 1 : Gp;
    return std::make_tuple(z, p, H0, B, A);
}

auto AF::detail::elliptic_lp_prototype(const mpfr::mpreal &wp, const mpfr::mpreal &ws, const mpfr::mpreal &Ap,
                                       const mpfr::mpreal &As) -> design_res {
    //用四个指标确定模拟椭圆低通滤波器
    mpreal Gp = pow(10_mpr, -Ap / 20_mpr);
    mpreal ep = sqrt(pow(10_mpr, Ap / 10_mpr) - 1);//通带波纹系数
    mpreal es = sqrt(pow(10_mpr, As / 10_mpr) - 1);//阻带波纹系数
    mpreal k1 = ep / es;
    auto [N, wp0] = AF::detail::ellipitic_filter_order(wp, ws, Ap, As, lowpass);
    wp0 = wp;//保持通带截止频率不变,这是一种策略
    mpreal k = Ellipitic::deg(N, k1);
    uint32_t L = N / 2;
    uint32_t r = N % 2;
    auto v0 = -mpcomplex(0, 1) * Ellipitic::asn(mpcomplex(0, 1) / ep, k1) / mpreal(N);
    std::vector<mpcomplex> z;
    std::vector<mpcomplex> p;
    if (r == 1) {
        auto p0 = wp0 * mpcomplex(0, 1) * Ellipitic::sn(mpcomplex(0, 1) * v0, k);
        p.push_back(p0);
    }
    for (uint32_t i = 1; i <= L; i++) {
        auto ui = (2_mpr * i - 1_mpr) / N;
        auto zeta_i = Ellipitic::cd(ui, k);
        auto zi = wp0 * mpcomplex(0, 1) / (k * zeta_i);
        auto pi = wp0 * mpcomplex(0, 1) * Ellipitic::cd(ui - mpcomplex(0, 1) * v0, k);

        p.push_back(pi);
        p.push_back(conj(pi));

        z.push_back(zi);
        z.push_back(conj(zi));
    }

    return detail::zp_trans(z,p,Gp);
}


auto AF::detail::ellipitic_filter_order(const mpreal &wp, const mpreal &ws, const mpreal &Ap, const mpreal &As,
                                        filter_band_type type) -> std::tuple<uint32_t, mpfr::mpreal> {
    using mpfr::mpreal;
    mpreal wn0;
    uint32_t N;
    if (type == filter_band_type::lowpass) {
        mpreal k = wp / ws;
        mpreal eps_p = sqrt(pow(10, Ap / 10) - 1);
        mpreal eps_s = sqrt(pow(10, As / 10) - 1);
        mpreal k1 = eps_p / eps_s;
        auto [K1, K1p] = Ellipitic::EllipticK(k1);
        auto [K, Kp] = Ellipitic::EllipticK(k);
        mpreal N_exact = (K1p / K1) / (Kp / K);//确定滤波器阶数
        N = ceil(N_exact).toULong();
        k = Ellipitic::deg(N, k1);
        wn0 = ws * k;//新的通带边界频率,也可以计算新的阻带
    } else if (type == filter_band_type::highpass) {
        //进行频率变换,lambda_p=1
        mpreal lambda_p = 1_mpr;
        mpreal lambda_s = wp / ws;
        mpreal k = lambda_p / lambda_s;
        mpreal eps_p = sqrt(pow(10, Ap / 10) - 1);
        mpreal eps_s = sqrt(pow(10, As / 10) - 1);
        mpreal k1 = eps_p / eps_s;
        auto [K1, K1p] = Ellipitic::EllipticK(k1);
        auto [K, Kp] = Ellipitic::EllipticK(k);
        mpreal N_exact = (K1p / K1) / (Kp / K);//确定滤波器阶数
        N = ceil(N_exact).toULong();
        k = Ellipitic::deg(N, k1);
        mpreal lambda_s0 = lambda_p / k;//新的归一化阻带边界频率,也可以计算新的阻带
        wn0 = wp / lambda_s0;
    }
    return std::make_tuple(N, wn0);
}

auto
AF::detail::ellipitic_filter_order(const mpreal &wpu, const mpreal &wpl, const mpreal &wsu, const mpreal &wsl,
                                   const mpreal &Ap,
                                   const mpreal &As,
                                   filter_band_type type) -> std::tuple<uint32_t, mpfr::mpreal, mpfr::mpreal> {
    //TODO:添加policy,使这个调整策略可以选择
    if (type == filter_band_type::bandpass) {
        //按课本写的可能会有问题,到时候看看英文版
        mpreal BW = wpu - wpl;
        mpreal w0 = sqrt(wpl * wpu);
        mpreal wsl_p = wsl - w0 * w0 / wsl;
        mpreal wsu_p = wsu - w0 * w0 / wsu;
        mpreal ws_p = min(abs(wsl_p), abs(wsu_p));
        const mpreal& wp_p = BW;
        auto [N, wp0] = ellipitic_filter_order(wp_p, ws_p, Ap, As, lowpass);
        mpreal wpl0 = (sqrt(wp0 * wp0 + 4 * w0 * w0) - wp0) / 2_mpr;
        mpreal wpu0 = (sqrt(wp0 * wp0 + 4 * w0 * w0) + wp0) / 2_mpr;
        return {N, wpl0, wpu0};
    } else if (type == filter_band_type::bandstop) {
        mpfr::mpreal w0 = sqrt(wsl * wsu);//中心频率
        mpreal BW = wsu - wsl;
        mpreal wpl_p = 1_mpr / (wpl - w0 * w0 / wpl);
        mpreal wpu_p = 1_mpr / (wpu - w0 * w0 / wpu);
        mpreal ws_p = 1_mpr / BW;
        mpreal wp_p = max(abs(wpl_p), abs(wpu_p));
        auto [N, wp0] = ellipitic_filter_order(wp_p, ws_p, Ap, As, lowpass);
        mpreal wpl0 = (sqrt(1_mpr / (wp0 * wp0) + 4 * w0 * w0) - 1_mpr / wp0) / 2_mpr;
        mpreal wpu0 = (sqrt(1_mpr / (wp0 * wp0) + 4 * w0 * w0) + 1_mpr / wp0) / 2_mpr;
        return {N, wpl0, wpu0};
    } else {
        return {};
    }
}

auto AF::detail::zp_trans(zeros &zs, poles &ps, const mpreal& Gp) -> design_res {
    auto N=ps.size();
    auto L=N/2;
    auto r=N%2;
    std::vector<std::array<mpreal, 3>> B;//分子
    std::vector<std::array<mpreal, 3>> A;//分母
    if(r==0){
        B.push_back({Gp,0,0});
        A.push_back({1,0,0});
    }else{
        auto p0=ps[0];
        B.push_back({1, 0, 0});
        A.push_back({1, -(1_mpr / p0).real(), 0});
    }
    for(uint32_t i=r;i<=2*L+r;i+=2){
        auto zi=zs[i-r];
        auto pi=ps[i];
        B.push_back({1, -2 * (1_mpr / zi).real(), 1_mpr / abs(zi) / abs(zi)});
        A.push_back({1, -2 * (1_mpr / pi).real(), 1_mpr / abs(pi) / abs(pi)});
    }
    auto H0 = r ? 1_mpr : Gp;
    return std::make_tuple(zs, ps, H0, B, A);
}

auto AF::ellipitic_filter(const mpreal &Wp, const mpreal &Ws, const mpreal &Ap, const mpreal &As,
                          filter_band_type type) -> design_res {
    if (type == lowpass) {
        return AF::detail::elliptic_lp_prototype(Wp, Ws, Ap, As);
    } else if (type == highpass) {
        //频率变换
        mpreal wp_p = 1_mpr / Wp;
        mpreal ws_p = 1_mpr / Ws;
        //算等效低通
        auto [z, p, H0, B, A] = AF::detail::elliptic_lp_prototype(wp_p, ws_p, Ap, As);
        //频率逆变换
        for (auto &li: B) {
            std::reverse(li.begin(), li.end());
        }
        for (auto &li: A) {
            std::reverse(li.begin(), li.end());
        }
        for (auto &zi: z) {
            zi = 1_mpr / zi;
        }
        for (auto &pi: p) {
            pi = 1_mpr / pi;
        }

        return {z, p, H0, B, A};
    }
    return {};
}

auto AF::ellipitic_filter(const mpreal &Wpu, const mpreal &Wpl, const mpreal &Wsu, const mpreal &Wsl, const mpreal &Ap,
                          const mpreal &As,
                          filter_band_type type) -> design_res {
    if (type == bandpass) {

        //频率变换
        mpreal BW = Wpu - Wpl;
        mpreal w0 = sqrt(Wpl * Wpu);
        mpreal wsl_p = Wsl - w0 * w0 / Wsl;
        mpreal wsu_p = Wsu - w0 * w0 / Wsu;
        mpreal ws_p = min(abs(wsl_p), abs(wsu_p));
        const mpreal& wp_p = BW;

        //算等效低通
        auto [z0, p0, H0, B0, A0] = AF::detail::elliptic_lp_prototype(wp_p, ws_p, Ap, As);

        //频率逆变换
        //从B和A入手有一点点难度,考虑从零极点的变换入手
        std::vector<mpcomplex> z;
        std::vector<mpcomplex> p;
        std::vector<std::array<mpreal, 3>> B;//分子
        std::vector<std::array<mpreal, 3>> A;//分母
        uint32_t r = p0.size() % 2;
        uint32_t L = p0.size() / 2;
        mpreal Gp = pow(10_mpr, -Ap / 20_mpr);
        if (r == 1) {
            auto p_first = p0[0];
            A.push_back({1, 1_mpr / (A0[0][1] * w0 * w0), 1_mpr / (w0 * w0)});
            B.push_back({0, 1_mpr / (A0[0][1] * w0 * w0), 0});
            p.push_back((p_first + sqrt(p_first * p_first - 4_mpr * w0 * w0)) / 2_mpr);
            p.push_back((p_first - sqrt(p_first * p_first - 4_mpr * w0 * w0)) / 2_mpr);
        } else {
            B.push_back({Gp, 0, 0});
            A.push_back({1, 0, 0});
        }
        for (uint32_t i = r; i < 2 * L + r; i += 2) {
            auto z0i = z0[i - r];
            //auto z0i_conj = conj(z0i);
            auto p0i = p0[i];
            //auto p0i_conj = conj(p0i);

            //第一组
            auto zi1 = (z0i + sqrt(z0i * z0i - 4_mpr * w0 * w0)) / 2_mpr;
            z.push_back(zi1);
            z.push_back(conj(zi1));

            auto pi1 = (p0i + sqrt(p0i * p0i - 4_mpr * w0 * w0)) / 2_mpr;
            p.push_back(pi1);
            p.push_back(conj(pi1));
            auto B2 = B0[i / 2 + 1][2];
            auto A2 = A0[(i - r) / 2 + 1][2];

            B.push_back({B2, -2 * B2 * (1_mpr / zi1).real(), B2 / abs(zi1) / abs(zi1)});
            A.push_back({A2, -2 * A2 * (1_mpr / pi1).real(), A2 / abs(pi1) / abs(pi1)});

            //第二组
            auto zi2 = (z0i - sqrt(z0i * z0i - 4_mpr * w0 * w0)) / 2_mpr;
            z.push_back(zi2);
            z.push_back(conj(zi2));

            auto pi2 = (p0i - sqrt(p0i * p0i - 4_mpr * w0 * w0)) / 2_mpr;
            p.push_back(pi2);
            p.push_back(conj(pi2));

            //B2 = 1_mpr / abs(zi2) / abs(zi2);
            //A2 = 1_mpr / abs(pi2) / abs(pi2);

            B.push_back({1_mpr, -2 * (1_mpr / zi2).real(), 1_mpr / abs(zi2) / abs(zi2)});
            A.push_back({1_mpr, -2 * (1_mpr / pi2).real(), 1_mpr / abs(pi2) / abs(pi2)});
        }
        return {z, p, H0, B, A};


    } else if (type == bandstop) {

        //频率变换
        mpfr::mpreal w0 = sqrt(Wsl * Wsu);//中心频率
        mpreal BW = Wsu - Wsl;
        mpreal wpl_p = 1_mpr / (Wpl - w0 * w0 / Wpl);
        mpreal wpu_p = 1_mpr / (Wpu - w0 * w0 / Wpu);
        mpreal ws_p = 1_mpr / BW;
        mpreal wp_p = max(abs(wpl_p), abs(wpu_p));

        //算等效低通
        auto [z0, p0, H0, B0, A0] = AF::detail::elliptic_lp_prototype(wp_p, ws_p, Ap, As);

        //频率逆变换
        std::vector<mpcomplex> z;
        std::vector<mpcomplex> p;
        std::vector<std::array<mpreal, 3>> B;//分子
        std::vector<std::array<mpreal, 3>> A;//分母
        uint32_t r = p0.size() % 2;
        uint32_t L = p0.size() / 2;
        mpreal Gp = pow(10_mpr, -Ap / 20_mpr);
        if (r == 1) {
            auto p_first = p0[0];
            p_first = 1_mpr / p_first;
            A.push_back({1, 1_mpr / (A0[0][1] * w0 * w0), 1_mpr / (w0 * w0)});
            B.push_back({0, 1_mpr / (A0[0][1] * w0 * w0), 0});
            p.push_back((p_first + sqrt(p_first * p_first - 4_mpr * w0 * w0)) / 2_mpr);
            p.push_back((p_first - sqrt(p_first * p_first - 4_mpr * w0 * w0)) / 2_mpr);
        } else {
            B.push_back({Gp, 0, 0});
            A.push_back({1, 0, 0});
        }
        for (uint32_t i = r; i < 2 * L + r; i += 2) {
            auto z0i = z0[i - r];
            z0i = 1_mpr / z0i;
            //auto z0i_conj = conj(z0i);
            auto p0i = p0[i];
            p0i = 1_mpr / p0i;
            //auto p0i_conj = conj(p0i);

            //第一组
            auto zi1 = (z0i + sqrt(z0i * z0i - 4_mpr * w0 * w0)) / 2_mpr;
            z.push_back(zi1);
            z.push_back(conj(zi1));

            auto pi1 = (p0i + sqrt(p0i * p0i - 4_mpr * w0 * w0)) / 2_mpr;
            p.push_back(pi1);
            p.push_back(conj(pi1));

            B.push_back({1_mpr, -2 * (1_mpr / zi1).real(), 1_mpr / abs(zi1) / abs(zi1)});
            A.push_back({1_mpr, -2 * (1_mpr / pi1).real(), 1_mpr / abs(pi1) / abs(pi1)});

            //第二组
            auto zi2 = (z0i - sqrt(z0i * z0i - 4_mpr * w0 * w0)) / 2_mpr;
            z.push_back(zi2);
            z.push_back(conj(zi2));

            auto pi2 = (p0i - sqrt(p0i * p0i - 4_mpr * w0 * w0)) / 2_mpr;
            p.push_back(pi2);
            p.push_back(conj(pi2));

            B.push_back({1, -2 * (1_mpr / zi2).real(), 1_mpr / abs(zi2) / abs(zi2)});
            A.push_back({1, -2 * (1_mpr / pi2).real(), 1_mpr / abs(pi2) / abs(pi2)});
        }
        return {z, p, H0, B, A};

    }
    return {};
}


