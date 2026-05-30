/** @file
 *****************************************************************************

 Implementation of interfaces for multi-exponentiation routines.

 See multiexp.hpp .

 *****************************************************************************
 * @author     This file is part of libff, developed by SCIPR Lab
 *             and contributors (see AUTHORS).
 * @copyright  MIT license (see LICENSE file)
 *****************************************************************************/

#ifndef MULTIEXP_TCC_
#define MULTIEXP_TCC_

#include <algorithm>
#include <cassert>
#include <type_traits>
#include <typeinfo>
#include <string>    
#include <fstream>  
#include <iostream>
#include <utility>
#include <vector>
#include <cstdlib> // std::system 호출을 위해 반드시 추가!

#include <libff/algebra/fields/bigint.hpp>
#include <libff/algebra/fields/fp_aux.tcc>
#include <libff/algebra/scalar_multiplication/multiexp.hpp>
#include <libff/algebra/scalar_multiplication/wnaf.hpp>
#include <libff/common/profiling.hpp>
#include <libff/common/utils.hpp>

namespace libff {

template<mp_size_t n>
class ordered_exponent {
// to use std::push_heap and friends later
public:
    size_t idx;
    bigint<n> r;

    ordered_exponent(const size_t idx, const bigint<n> &r) : idx(idx), r(r) {};

    bool operator<(const ordered_exponent<n> &other) const
    {
#if defined(__x86_64__) && defined(USE_ASM)
        if (n == 3)
        {
            long res;
            __asm__
                ("// check for overflow           \n\t"
                 "mov $0, %[res]                  \n\t"
                 ADD_CMP(16)
                 ADD_CMP(8)
                 ADD_CMP(0)
                 "jmp done%=                      \n\t"
                 "subtract%=:                     \n\t"
                 "mov $1, %[res]                  \n\t"
                 "done%=:                         \n\t"
                 : [res] "=&r" (res)
                 : [A] "r" (other.r.data), [mod] "r" (this->r.data)
                 : "cc", "%rax");
            return res;
        }
        else if (n == 4)
        {
            long res;
            __asm__
                ("// check for overflow           \n\t"
                 "mov $0, %[res]                  \n\t"
                 ADD_CMP(24)
                 ADD_CMP(16)
                 ADD_CMP(8)
                 ADD_CMP(0)
                 "jmp done%=                      \n\t"
                 "subtract%=:                     \n\t"
                 "mov $1, %[res]                  \n\t"
                 "done%=:                         \n\t"
                 : [res] "=&r" (res)
                 : [A] "r" (other.r.data), [mod] "r" (this->r.data)
                 : "cc", "%rax");
            return res;
        }
        else if (n == 5)
        {
            long res;
            __asm__
                ("// check for overflow           \n\t"
                 "mov $0, %[res]                  \n\t"
                 ADD_CMP(32)
                 ADD_CMP(24)
                 ADD_CMP(16)
                 ADD_CMP(8)
                 ADD_CMP(0)
                 "jmp done%=                      \n\t"
                 "subtract%=:                     \n\t"
                 "mov $1, %[res]                  \n\t"
                 "done%=:                         \n\t"
                 : [res] "=&r" (res)
                 : [A] "r" (other.r.data), [mod] "r" (this->r.data)
                 : "cc", "%rax");
            return res;
        }
        else
#endif
        {
            return (mpn_cmp(this->r.data, other.r.data, n) < 0);
        }
    }
};

/**
 * multi_exp_inner<T, FieldT, Method>() implementes the specified
 * multiexponentiation method.
 * this implementation relies on some rather arcane template magic:
 * function templates cannot be partially specialized, so we cannot just write
 *     template<typename T, typename FieldT>
 *     T multi_exp_inner<T, FieldT, multi_exp_method_naive>
 * thus we resort to using std::enable_if. the basic idea is that *overloading*
 * is what's actually happening here, it's just that, for any given value of
 * Method, only one of the templates will be valid, and thus the correct
 * implementation will be used.
 */

template<typename T, typename FieldT, multi_exp_method Method,
    typename std::enable_if<(Method == multi_exp_method_naive), int>::type = 0>
T multi_exp_inner(
    typename std::vector<T>::const_iterator vec_start,
    typename std::vector<T>::const_iterator vec_end,
    typename std::vector<FieldT>::const_iterator scalar_start,
    typename std::vector<FieldT>::const_iterator scalar_end)
{
    T result(T::zero());

    typename std::vector<T>::const_iterator vec_it;
    typename std::vector<FieldT>::const_iterator scalar_it;

    for (vec_it = vec_start, scalar_it = scalar_start; vec_it != vec_end; ++vec_it, ++scalar_it)
    {
        bigint<FieldT::num_limbs> scalar_bigint = scalar_it->as_bigint();
        result = result + opt_window_wnaf_exp(*vec_it, scalar_bigint, scalar_bigint.num_bits());
    }
    assert(scalar_it == scalar_end);

    return result;
}


template<typename T, typename FieldT, multi_exp_method Method,
    typename std::enable_if<(Method == multi_exp_method_naive_plain), int>::type = 0>
T multi_exp_inner(
    typename std::vector<T>::const_iterator vec_start,
    typename std::vector<T>::const_iterator vec_end,
    typename std::vector<FieldT>::const_iterator scalar_start,
    typename std::vector<FieldT>::const_iterator scalar_end)
{
    T result(T::zero());

    typename std::vector<T>::const_iterator vec_it;
    typename std::vector<FieldT>::const_iterator scalar_it;

    for (vec_it = vec_start, scalar_it = scalar_start; vec_it != vec_end; ++vec_it, ++scalar_it)
    {
        result = result + (*scalar_it) * (*vec_it);
    }
    assert(scalar_it == scalar_end);

    return result;
}


template<typename T, typename FieldT, multi_exp_method Method,
    typename std::enable_if<(Method == multi_exp_method_BDLO12), int>::type = 0>
T multi_exp_inner(
    typename std::vector<T>::const_iterator bases,
    typename std::vector<T>::const_iterator bases_end,
    typename std::vector<FieldT>::const_iterator exponents,
    typename std::vector<FieldT>::const_iterator exponents_end)
{
    UNUSED(exponents_end);
    size_t length = bases_end - bases;
    size_t log2_length = log2(length);
    size_t c = log2_length - (log2_length / 3 - 2);

    const mp_size_t exp_num_limbs =
        std::remove_reference<decltype(*exponents)>::type::num_limbs;
    std::vector<bigint<exp_num_limbs> > bn_exponents(length);
    size_t num_bits = 0;

    for (size_t i = 0; i < length; i++)
    {
        bn_exponents[i] = exponents[i].as_bigint();
        num_bits = std::max(num_bits, bn_exponents[i].num_bits());
    }

    size_t num_groups = (num_bits + c - 1) / c;

    T result;
    bool result_nonzero = false;

    for (size_t k = num_groups - 1; k <= num_groups; k--)
    {
        if (result_nonzero)
        {
            for (size_t i = 0; i < c; i++)
            {
                result = result.dbl();
            }
        }

        std::vector<T> buckets(1 << c);
        std::vector<bool> bucket_nonzero(1 << c);

        for (size_t i = 0; i < length; i++)
        {
            size_t id = 0;
            for (size_t j = 0; j < c; j++)
            {
                if (bn_exponents[i].test_bit(k*c + j))
                {
                    id |= 1 << j;
                }
            }

            if (id == 0)
            {
                continue;
            }

            if (bucket_nonzero[id])
            {
#ifdef USE_MIXED_ADDITION
                buckets[id] = buckets[id].mixed_add(bases[i]);
#else
                buckets[id] = buckets[id] + bases[i];
#endif
            }
            else
            {
                buckets[id] = bases[i];
                bucket_nonzero[id] = true;
            }
        }

#ifdef USE_MIXED_ADDITION
        batch_to_special(buckets);
#endif

        T running_sum;
        bool running_sum_nonzero = false;

        for (size_t i = (1u << c) - 1; i > 0; i--)
        {
            if (bucket_nonzero[i])
            {
                if (running_sum_nonzero)
                {
#ifdef USE_MIXED_ADDITION
                    running_sum = running_sum.mixed_add(buckets[i]);
#else
                    running_sum = running_sum + buckets[i];
#endif
                }
                else
                {
                    running_sum = buckets[i];
                    running_sum_nonzero = true;
                }
            }

            if (running_sum_nonzero)
            {
                if (result_nonzero)
                {
                    result = result + running_sum;
                }
                else
                {
                    result = running_sum;
                    result_nonzero = true;
                }
            }
        }
    }

    return result;
}


template<typename T, typename FieldT, multi_exp_method Method,
    typename std::enable_if<(Method == multi_exp_method_bos_coster), int>::type = 0>
T multi_exp_inner(
    typename std::vector<T>::const_iterator vec_start,
    typename std::vector<T>::const_iterator vec_end,
    typename std::vector<FieldT>::const_iterator scalar_start,
    typename std::vector<FieldT>::const_iterator scalar_end)
{
    const mp_size_t n = std::remove_reference<decltype(*scalar_start)>::type::num_limbs;

    if (vec_start == vec_end)
    {
        return T::zero();
    }

    if (vec_start + 1 == vec_end)
    {
        return (*scalar_start)*(*vec_start);
    }

    std::vector<ordered_exponent<n> > opt_q;
    const size_t vec_len = scalar_end - scalar_start;
    const size_t odd_vec_len = (vec_len % 2 == 1 ? vec_len : vec_len + 1);
    opt_q.reserve(odd_vec_len);
    std::vector<T> g;
    g.reserve(odd_vec_len);

    typename std::vector<T>::const_iterator vec_it;
    typename std::vector<FieldT>::const_iterator scalar_it;
    size_t i;
    for (i=0, vec_it = vec_start, scalar_it = scalar_start; vec_it != vec_end; ++vec_it, ++scalar_it, ++i)
    {
        g.emplace_back(*vec_it);

        opt_q.emplace_back(ordered_exponent<n>(i, scalar_it->as_bigint()));
    }
    std::make_heap(opt_q.begin(),opt_q.end());
    assert(scalar_it == scalar_end);

    if (vec_len != odd_vec_len)
    {
        g.emplace_back(T::zero());
        opt_q.emplace_back(ordered_exponent<n>(odd_vec_len - 1, bigint<n>(0ul)));
    }
    assert(g.size() % 2 == 1);
    assert(opt_q.size() == g.size());

    T opt_result = T::zero();

    while (true)
    {
        ordered_exponent<n> &a = opt_q[0];
        ordered_exponent<n> &b = (opt_q[1] < opt_q[2] ? opt_q[2] : opt_q[1]);

        const size_t abits = a.r.num_bits();

        if (b.r.is_zero())
        {
            // opt_result = opt_result + (a.r * g[a.idx]);
            opt_result = opt_result + opt_window_wnaf_exp(g[a.idx], a.r, abits);
            break;
        }

        const size_t bbits = b.r.num_bits();
        const size_t limit = (abits-bbits >= 20 ? 20 : abits-bbits);

        if (bbits < 1ul<<limit)
        {
            /*
              In this case, exponentiating to the power of a is cheaper than
              subtracting b from a multiple times, so let's do it directly
            */
            // opt_result = opt_result + (a.r * g[a.idx]);
            opt_result = opt_result + opt_window_wnaf_exp(g[a.idx], a.r, abits);
#ifdef DEBUG
            printf("Skipping the following pair (%zu bit number vs %zu bit):\n", abits, bbits);
            a.r.print();
            b.r.print();
#endif
            a.r.clear();
        }
        else
        {
            // x A + y B => (x-y) A + y (B+A)
            mpn_sub_n(a.r.data, a.r.data, b.r.data, n);
            g[b.idx] = g[b.idx] + g[a.idx];
        }

        // regardless of whether a was cleared or subtracted from we push it down, then take back up

        /* heapify A down */
        size_t a_pos = 0;
        while (2*a_pos + 2< odd_vec_len)
        {
            // this is a max-heap so to maintain a heap property we swap with the largest of the two
            if (opt_q[2*a_pos+1] < opt_q[2*a_pos+2])
            {
                std::swap(opt_q[a_pos], opt_q[2*a_pos+2]);
                a_pos = 2*a_pos+2;
            }
            else
            {
                std::swap(opt_q[a_pos], opt_q[2*a_pos+1]);
                a_pos = 2*a_pos+1;
            }
        }

        /* now heapify A up appropriate amount of times */
        while (a_pos > 0 && opt_q[(a_pos-1)/2] < opt_q[a_pos])
        {
            std::swap(opt_q[a_pos], opt_q[(a_pos-1)/2]);
            a_pos = (a_pos-1) / 2;
        }
    }

    return opt_result;
}

template<typename G1_PointT, typename ScalarFieldT>
G1_PointT multi_exp_g1_prove_fpga(
    typename std::vector<G1_PointT>::const_iterator vec_start,
    typename std::vector<G1_PointT>::const_iterator vec_end,
    typename std::vector<ScalarFieldT>::const_iterator scalar_start,
    typename std::vector<ScalarFieldT>::const_iterator scalar_end,
    const size_t chunks)
{
    libff::UNUSED(chunks);

    const size_t length = vec_end - vec_start;
    if (length == 0) {
        return G1_PointT::zero();
    }

    using BaseFieldT = typename std::decay<decltype(vec_start->X)>::type;
    using BaseBigIntT = decltype(vec_start->X.as_bigint());
    using ScalarBigIntT = decltype(scalar_start->as_bigint());

    static const BaseFieldT base_R255 = BaseFieldT(2) ^ 255;

    std::vector<BaseBigIntT> dma_points_X(length);
    std::vector<BaseBigIntT> dma_points_Y(length);
    std::vector<ScalarBigIntT> dma_scalars(length);

    std::vector<G1_PointT> affine_bases(vec_start, vec_end);
    batch_to_special(affine_bases);

    #ifdef MULTICORE
    #pragma omp parallel for
    #endif
    for (size_t i = 0; i < length; ++i) {
        // 🌟 핵심 방어 로직: 해당 점이 항등원(Infinity)이라면 하드웨어가 무시하도록 스칼라를 0으로 강제 세팅
        if (affine_bases[i].is_zero()) {
            // 컴파일러의 타입 추론 모호성을 피하기 위해 명시적으로 0UL(Unsigned Long)을 사용합니다.
            dma_scalars[i]  = 0UL;
            dma_points_X[i] = 0UL;
            dma_points_Y[i] = 0UL;
        } else {
            dma_scalars[i]  = scalar_start[i].as_bigint();
            dma_points_X[i] = (affine_bases[i].X * base_R255).as_bigint();
            dma_points_Y[i] = (affine_bases[i].Y * base_R255).as_bigint();
        }
    }

    // =====================================================================
    // 💾 1. [인풋 평탄화 및 txt 덤프 (Python 입력용)]
    // =====================================================================
    enter_block("Dump FPGA testbench vectors");
    std::ofstream file_scalar("msm_input_scalar.txt");
    std::ofstream file_point_x("msm_input_x.txt");
    std::ofstream file_point_y("msm_input_y.txt");

    if (file_scalar.is_open() && file_point_x.is_open() && file_point_y.is_open()) {
        for (size_t i = 0; i < length; ++i) {
            file_scalar  << dma_scalars[i]  << "\n";
            file_point_x << dma_points_X[i] << "\n";
            file_point_y << dma_points_Y[i] << "\n";
        }
        file_scalar.close(); file_point_x.close(); file_point_y.close();
        print_indent(); printf("💾 [C++] 입력 데이터 %zu개 txt 덤프 완료.\n", length);
    } else {
        print_indent(); printf("⚠️ [C++ 오류] txt 파일을 생성할 수 없습니다.\n");
    }
    leave_block("Dump FPGA testbench vectors");

    // =====================================================================
    // 🚀 2. [진정한 MIL 검증: C++ 내부에서 Python 모델 동기적 실행]
    // =====================================================================
    enter_block("Run Python Golden Model (Co-simulation)");
    print_indent(); printf("🐍 [System] Python 하드웨어 에뮬레이터를 구동하여 최신 버킷을 생성합니다...\n");
    
    int ret = std::system("python3 msm_fpga_emu.py > python_emu.log 2>&1"); 
    
    if (ret != 0) {
        print_indent(); printf("⚠️ [오류] Python 스크립트 실행 실패!\n");
    } else {
        print_indent(); printf("✅ [성공] Python 야코비안 버킷 누적 완료 (결과 덤프됨).\n");
    }
    leave_block("Run Python Golden Model (Co-simulation)");

    // =====================================================================
    // ⚙️ 3. [하드웨어 실행 제어 구조체 선언]
    // =====================================================================
    const size_t WINDOW_COUNT = 24;
    const size_t BITS_PER_WINDOW = 11;
    const size_t BUCKETS_PER_WINDOW = 1 << BITS_PER_WINDOW;

    struct FpgaRawJacobian {
        BaseBigIntT X;
        BaseBigIntT Y;
        BaseBigIntT Z;
    };
    
    std::vector<std::vector<FpgaRawJacobian>> fpga_buckets(
        WINDOW_COUNT, std::vector<FpgaRawJacobian>(BUCKETS_PER_WINDOW)
    );

    // =====================================================================
    // 📥 4. [Python 생성 야코비안 버킷 파일 로드 (DMA 수신 모사)]
    // =====================================================================
    enter_block("Load FPGA buckets from Python model");
    std::ifstream file_bx("buckets_X.txt");
    std::ifstream file_by("buckets_Y.txt");
    std::ifstream file_bz("buckets_Z.txt");
    
    if (file_bx.is_open() && file_by.is_open() && file_bz.is_open()) {
        std::string str_x, str_y, str_z;
        for (size_t w = 0; w < WINDOW_COUNT; ++w) {
            for (size_t b = 0; b < BUCKETS_PER_WINDOW; ++b) {
                file_bx >> str_x;
                file_by >> str_y;
                file_bz >> str_z;
                
                fpga_buckets[w][b].X = BaseBigIntT(str_x.c_str());
                fpga_buckets[w][b].Y = BaseBigIntT(str_y.c_str());
                fpga_buckets[w][b].Z = BaseBigIntT(str_z.c_str());
            }
        }
        file_bx.close(); file_by.close(); file_bz.close();
        print_indent(); printf("📥 [C++] Python 최신 야코비안 버킷 데이터 로드 완료.\n");
    } else {
        print_indent(); printf("⚠️ [C++ 알림] buckets_*.txt 파일 로드 실패.\n");
    }
    leave_block("Load FPGA buckets from Python model");

    // =====================================================================
    // 🔄 5. [도메인 복원 및 대규모 일괄 아핀(Affine) 변환]
    // =====================================================================
    static const BaseFieldT inv_R255 = (BaseFieldT(2) ^ 255).inverse();
    
    std::vector<std::vector<G1_PointT>> window_buckets(
        WINDOW_COUNT, std::vector<G1_PointT>(BUCKETS_PER_WINDOW, G1_PointT::zero())
    );

    std::vector<G1_PointT> active_buckets;
    std::vector<std::pair<size_t, size_t>> active_indices;
    
    active_buckets.reserve(WINDOW_COUNT * BUCKETS_PER_WINDOW);
    active_indices.reserve(WINDOW_COUNT * BUCKETS_PER_WINDOW);

    for (size_t w = 0; w < WINDOW_COUNT; ++w) {
        for (size_t b = 1; b < BUCKETS_PER_WINDOW; ++b) {
            if (fpga_buckets[w][b].Z.is_zero()) continue;

            BaseFieldT restored_X = BaseFieldT(fpga_buckets[w][b].X) * inv_R255;
            BaseFieldT restored_Y = BaseFieldT(fpga_buckets[w][b].Y) * inv_R255;
            BaseFieldT restored_Z = BaseFieldT(fpga_buckets[w][b].Z) * inv_R255;

            active_buckets.emplace_back(G1_PointT(restored_X, restored_Y, restored_Z));
            active_indices.emplace_back(w, b);
        }
    }

    if (!active_buckets.empty()) {
        batch_to_special(active_buckets);
    }

    for (size_t i = 0; i < active_buckets.size(); ++i) {
        size_t w = active_indices[i].first;
        size_t b = active_indices[i].second;
        window_buckets[w][b] = active_buckets[i];
    }

    // =====================================================================
    // 🧮 6. [소프트웨어 리덕션 (Pippenger 누적)]
    // =====================================================================
    std::vector<G1_PointT> window_results(WINDOW_COUNT, G1_PointT::zero());

    #ifdef MULTICORE
    #pragma omp parallel for
    #endif
    for (size_t w = 0; w < WINDOW_COUNT; ++w) {
        G1_PointT running_sum = G1_PointT::zero();
        G1_PointT window_sum = G1_PointT::zero();
        bool running_sum_nonzero = false, window_sum_nonzero = false;

        for (size_t b = BUCKETS_PER_WINDOW - 1; b > 0; --b) {
            if (!window_buckets[w][b].is_zero()) {
                if (running_sum_nonzero) {
                    running_sum = running_sum + window_buckets[w][b];
                } else {
                    running_sum = window_buckets[w][b]; 
                    running_sum_nonzero = true;
                }
            }

            if (running_sum_nonzero) {
                if (window_sum_nonzero) {
                    window_sum = window_sum + running_sum;
                } else {
                    window_sum = running_sum; 
                    window_sum_nonzero = true;
                }
            }
        }
        window_results[w] = window_sum;
    }

    // =====================================================================
    // 🔗 7. [최종 어그리게이션 (Horner's Method)]
    // =====================================================================
    G1_PointT final_result = G1_PointT::zero();
    bool final_nonzero = false;

    for (size_t w = WINDOW_COUNT - 1; w < WINDOW_COUNT; --w) {
        if (final_nonzero) {
            for (size_t shift = 0; shift < BITS_PER_WINDOW; ++shift) {
                final_result = final_result.dbl();
            }
        }
        if (!window_results[w].is_zero()) {
            if (final_nonzero) {
                final_result = final_result + window_results[w];
            } else {
                final_result = window_results[w]; 
                final_nonzero = true;
            }
        }
    }

    // =====================================================================
    // 🏁 8. [최종 검증 및 출력]
    // =====================================================================
    final_result.to_affine_coordinates();
    
    printf("\n==============================================\n");
    printf(" 🏁 [C++ 복원 및 리덕션] 최종 산출 좌표\n");
    printf(" X: "); final_result.X.print(); 
    printf(" Y: "); final_result.Y.print(); 
    printf("==============================================\n\n");

    return final_result;
}

template<typename T, typename FieldT, multi_exp_method Method>
T multi_exp(typename std::vector<T>::const_iterator vec_start,
            typename std::vector<T>::const_iterator vec_end,
            typename std::vector<FieldT>::const_iterator scalar_start,
            typename std::vector<FieldT>::const_iterator scalar_end,
            const size_t chunks)
{
    const size_t total = vec_end - vec_start;
    if ((total < chunks) || (chunks == 1))
    {
        // no need to split into "chunks", can call implementation directly
        return multi_exp_inner<T, FieldT, Method>(
            vec_start, vec_end, scalar_start, scalar_end);
    }

    const size_t one = total/chunks;

    std::vector<T> partial(chunks, T::zero());

#ifdef MULTICORE
#pragma omp parallel for
#endif
    for (size_t i = 0; i < chunks; ++i)
    {
        partial[i] = multi_exp_inner<T, FieldT, Method>(
             vec_start + i*one,
             (i == chunks-1 ? vec_end : vec_start + (i+1)*one),
             scalar_start + i*one,
             (i == chunks-1 ? scalar_end : scalar_start + (i+1)*one));
    }

    T final = T::zero();

    for (size_t i = 0; i < chunks; ++i)
    {
        final = final + partial[i];
    }

    return final;
}


template<typename T, typename FieldT, multi_exp_method Method>
T multi_exp_with_mixed_addition(typename std::vector<T>::const_iterator vec_start,
                                typename std::vector<T>::const_iterator vec_end,
                                typename std::vector<FieldT>::const_iterator scalar_start,
                                typename std::vector<FieldT>::const_iterator scalar_end,
                                const size_t chunks, 
                                bool useFPGA)
{
#ifndef NDEBUG
    assert(std::distance(vec_start, vec_end) == std::distance(scalar_start, scalar_end));
#else
    libff::UNUSED(vec_end);
#endif
    enter_block("Process scalar vector");
    auto value_it = vec_start;
    auto scalar_it = scalar_start;

    const FieldT zero = FieldT::zero();
    const FieldT one = FieldT::one();
    std::vector<FieldT> p;
    std::vector<T> g;

    T acc = T::zero();

    size_t num_skip = 0;
    size_t num_add = 0;
    size_t num_other = 0;

    for (; scalar_it != scalar_end; ++scalar_it, ++value_it)
    {
        if (*scalar_it == zero)
        {
            // do nothing
            ++num_skip;
        }
        else if (*scalar_it == one)
        {
#ifdef USE_MIXED_ADDITION
            acc = acc.mixed_add(*value_it);
#else
            acc = acc + (*value_it);
#endif
            ++num_add;
        }
        else
        {
            p.emplace_back(*scalar_it);
            g.emplace_back(*value_it);
            ++num_other;
        }
    }
    print_indent(); printf("* Elements of w skipped: %zu (%0.2f%%)\n", num_skip, 100.*num_skip/(num_skip+num_add+num_other));
    print_indent(); printf("* Elements of w processed with special addition: %zu (%0.2f%%)\n", num_add, 100.*num_add/(num_skip+num_add+num_other));
    print_indent(); printf("* Elements of w remaining: %zu (%0.2f%%)\n", num_other, 100.*num_other/(num_skip+num_add+num_other));
    leave_block("Process scalar vector");

    // 하드웨어 라우팅 분기 
    if(useFPGA){
        return acc + multi_exp_g1_prove_fpga<T, FieldT>(g.begin(), g.end(), p.begin(), p.end(),chunks);
    }
    else{
        return acc + multi_exp<T, FieldT, Method>(g.begin(), g.end(), p.begin(), p.end(), chunks);
    }
}

template <typename T>
T inner_product(typename std::vector<T>::const_iterator a_start,
                typename std::vector<T>::const_iterator a_end,
                typename std::vector<T>::const_iterator b_start,
                typename std::vector<T>::const_iterator b_end)
{
    return multi_exp<T, T, multi_exp_method_naive_plain>(
        a_start, a_end,
        b_start, b_end, 1);
}

template<typename T>
size_t get_exp_window_size(const size_t num_scalars)
{
    if (T::fixed_base_exp_window_table.empty())
    {
#ifdef LOWMEM
        return 14;
#else
        return 17;
#endif
    }
    size_t window = 1;
    for (long i = T::fixed_base_exp_window_table.size()-1; i >= 0; --i)
    {
#ifdef DEBUG
        if (!inhibit_profiling_info)
        {
            printf("%ld %zu %zu\n", i, num_scalars, T::fixed_base_exp_window_table[i]);
        }
#endif
        if (T::fixed_base_exp_window_table[i] != 0 && num_scalars >= T::fixed_base_exp_window_table[i])
        {
            window = i+1;
            break;
        }
    }

    if (!inhibit_profiling_info)
    {
        print_indent(); printf("Choosing window size %zu for %zu elements\n", window, num_scalars);
    }

#ifdef LOWMEM
    window = std::min((size_t)14, window);
#endif
    return window;
}

template<typename T>
window_table<T> get_window_table(const size_t scalar_size,
                                 const size_t window,
                                 const T &g)
{
    const size_t in_window = 1ul<<window;
    const size_t outerc = (scalar_size+window-1)/window;
    const size_t last_in_window = 1ul<<(scalar_size - (outerc-1)*window);
#ifdef DEBUG
    if (!inhibit_profiling_info)
    {
        print_indent(); printf("* scalar_size=%zu; window=%zu; in_window=%zu; outerc=%zu\n", scalar_size, window, in_window, outerc);
    }
#endif

    window_table<T> powers_of_g(outerc, std::vector<T>(in_window, T::zero()));

    T gouter = g;

    for (size_t outer = 0; outer < outerc; ++outer)
    {
        T ginner = T::zero();
        size_t cur_in_window = outer == outerc-1 ? last_in_window : in_window;
        for (size_t inner = 0; inner < cur_in_window; ++inner)
        {
            powers_of_g[outer][inner] = ginner;
            ginner = ginner + gouter;
        }

        for (size_t i = 0; i < window; ++i)
        {
            gouter = gouter + gouter;
        }
    }

    return powers_of_g;
}

template<typename T, typename FieldT>
T windowed_exp(const size_t scalar_size,
               const size_t window,
               const window_table<T> &powers_of_g,
               const FieldT &pow)
{
    const size_t outerc = (scalar_size+window-1)/window;
    const bigint<FieldT::num_limbs> pow_val = pow.as_bigint();

    /* exp */
    T res = powers_of_g[0][0];

    for (size_t outer = 0; outer < outerc; ++outer)
    {
        size_t inner = 0;
        for (size_t i = 0; i < window; ++i)
        {
            if (pow_val.test_bit(outer*window + i))
            {
                inner |= 1u << i;
            }
        }

        res = res + powers_of_g[outer][inner];
    }

    return res;
}

template<typename T, typename FieldT>
std::vector<T> batch_exp(const size_t scalar_size,
                         const size_t window,
                         const window_table<T> &table,
                         const std::vector<FieldT> &v)
{
    if (!inhibit_profiling_info)
    {
        print_indent();
    }
    std::vector<T> res(v.size(), table[0][0]);

#ifdef MULTICORE
#pragma omp parallel for
#endif
    for (size_t i = 0; i < v.size(); ++i)
    {
        res[i] = windowed_exp(scalar_size, window, table, v[i]);

        if (!inhibit_profiling_info && (i % 10000 == 0))
        {
            printf(".");
            fflush(stdout);
        }
    }

    if (!inhibit_profiling_info)
    {
        printf(" DONE!\n");
    }

    return res;
}

template<typename T, typename FieldT>
std::vector<T> batch_exp_with_coeff(const size_t scalar_size,
                                    const size_t window,
                                    const window_table<T> &table,
                                    const FieldT &coeff,
                                    const std::vector<FieldT> &v)
{
    if (!inhibit_profiling_info)
    {
        print_indent();
    }
    std::vector<T> res(v.size(), table[0][0]);

#ifdef MULTICORE
#pragma omp parallel for
#endif
    for (size_t i = 0; i < v.size(); ++i)
    {
        res[i] = windowed_exp(scalar_size, window, table, coeff * v[i]);

        if (!inhibit_profiling_info && (i % 10000 == 0))
        {
            printf(".");
            fflush(stdout);
        }
    }

    if (!inhibit_profiling_info)
    {
        printf(" DONE!\n");
    }

    return res;
}

template<typename T>
void batch_to_special(std::vector<T> &vec)
{
    enter_block("Batch-convert elements to special form");

    std::vector<T> non_zero_vec;
    for (size_t i = 0; i < vec.size(); ++i)
    {
        if (!vec[i].is_zero())
        {
            non_zero_vec.emplace_back(vec[i]);
        }
    }

    T::batch_to_special_all_non_zeros(non_zero_vec);
    auto it = non_zero_vec.begin();
    T zero_special = T::zero();
    zero_special.to_special();

    for (size_t i = 0; i < vec.size(); ++i)
    {
        if (!vec[i].is_zero())
        {
            vec[i] = *it;
            ++it;
        }
        else
        {
            vec[i] = zero_special;
        }
    }
    leave_block("Batch-convert elements to special form");
}

} // libff

#endif // MULTIEXP_TCC_

