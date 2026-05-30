#include <libff/common/default_types/ec_pp.hpp>
#include <libsnark/zk_proof_systems/ppzksnark/r1cs_gg_ppzksnark/r1cs_gg_ppzksnark.hpp>
#include <libsnark/gadgetlib1/protoboard.hpp>
#include <iostream>
#include <fstream>
#include <chrono>
#include <random>
#include <string>
#include <algorithm>

using namespace std;
using namespace libsnark;

// =================================================================
// 🎲 BN128 모듈러스 한계치를 지키는 완벽한 난수 생성 함수
// =================================================================
std::string generate_valid_scalar_dec(std::mt19937_64& gen) {
    // 0x30644e72e131a029b85045b68181585d2833e84879b9709143e1f593f0000001 의 10진수 값
    const std::string MODULUS = "21888242871839275222246405745257275088548364400416034343698204186575808495617";
    std::uniform_int_distribution<int> dis(0, 9);
    
    while (true) {
        std::string rand_num = "";
        // 77자리의 무작위 숫자 생성
        for (int i = 0; i < 77; ++i) {
            rand_num += std::to_string(dis(gen));
        }

        // 사전순(Lexicographical) 문자열 비교: rand_num이 MODULUS보다 작으면 합격!
        if (rand_num < MODULUS) {
            // 앞부분에 붙은 불필요한 '0'들 제거 (예: 000123 -> 123)
            size_t first_non_zero = rand_num.find_first_not_of('0');
            if (first_non_zero != std::string::npos) {
                rand_num.erase(0, first_non_zero);
            } else {
                rand_num = "0"; // 만약 전부 0이었다면 "0"으로 처리
            }
            return rand_num;
        }
        // 만약 모듈러스보다 크거나 같으면 버리고(Rejection) 다시 생성 (루프 반복)
    }
}

int main() {
    // 🔍 1. 내장 프로파일러 활성화 (MSM 소요 시간 상세 추적)
    libff::inhibit_profiling_info = false;
    libff::inhibit_profiling_counters = false;

    // 타원 곡선 파라미터 초기화
    libff::default_ec_pp::init_public_params();

    // =================================================================
    // 🛠️ 2. 하드웨어 스트레스 테스트용 R1CS 회로 및 데이터 생성
    // =================================================================
    const int NUM_CONSTRAINTS = 65535; // 원하는 데이터 개수만큼 조절 가능
    typedef libff::alt_bn128_Fr FieldT;
    
    protoboard<FieldT> pb;
    pb_variable_array<FieldT> A, B, C;
    A.allocate(pb, NUM_CONSTRAINTS, "A");
    B.allocate(pb, NUM_CONSTRAINTS, "B");
    C.allocate(pb, NUM_CONSTRAINTS, "C");

    // 고성능 난수 엔진 초기화
    std::random_device rd;
    std::mt19937_64 gen(rd()); 

    cout << "🛠️ BN128 모듈러스 이하의 엄격한 무작위 스칼라 생성 및 할당 중 (" << NUM_CONSTRAINTS << "개)..." << endl;

    for (int i = 0; i < NUM_CONSTRAINTS; ++i) {
        // 제약 조건 추가: A * B = C
        pb.add_r1cs_constraint(r1cs_constraint<FieldT>(A[i], B[i], C[i]));

        // 1️⃣ 모듈러스 한계선 아래의 안전한 랜덤 문자열 받아오기
        std::string rand_A = generate_valid_scalar_dec(gen);
        std::string rand_B = generate_valid_scalar_dec(gen);

        // 2️⃣ 스칼라 객체로 파싱하여 할당
        pb.val(A[i]) = FieldT(rand_A.c_str());
        pb.val(B[i]) = FieldT(rand_B.c_str());

        // 3️⃣ 증명이 깨지지 않도록 C 값은 A와 B의 곱으로 강제 재계산
        pb.val(C[i]) = pb.val(A[i]) * pb.val(B[i]);
    }

    // C 배열을 공개 입력(Public Input)으로 설정
    pb.set_input_sizes(NUM_CONSTRAINTS);

    cout << " ✅ 회로 및 스칼라 데이터 세팅 완료!\n" << endl;

    // Protoboard에서 시스템 및 입력 데이터 추출
    r1cs_constraint_system<FieldT> cs = pb.get_constraint_system();
    r1cs_primary_input<FieldT> primary_input = pb.primary_input();
    r1cs_auxiliary_input<FieldT> auxiliary_input = pb.auxiliary_input();

    // =================================================================
    // 🔑 3. 증명 키 생성 및 로드 (Trusted Setup 파일 입출력)
    // =================================================================
    bool GENERATE_NEW_KEYS = true; // 🚨 최초 1회만 true, 파일 생성 후 false로 변경!
    
    r1cs_gg_ppzksnark_proving_key<libff::alt_bn128_pp> pk;
    r1cs_gg_ppzksnark_verification_key<libff::alt_bn128_pp> vk;

    if (GENERATE_NEW_KEYS) {
        cout << "🔑 [Setup] 새로운 키 쌍 생성 및 파일 저장 중..." << endl;
        auto keypair = r1cs_gg_ppzksnark_generator<libff::alt_bn128_pp>(cs);
        pk = keypair.pk;
        vk = keypair.vk;

        std::ofstream pk_out("proving_key.bin", std::ios::binary); pk_out << pk; pk_out.close();
        std::ofstream vk_out("verification_key.bin", std::ios::binary); vk_out << vk; vk_out.close();
    } else {
        cout << "🔑 [Setup] 기존 키 쌍 불러오는 중..." << endl;
        std::ifstream pk_in("proving_key.bin", std::ios::binary); pk_in >> pk; pk_in.close();
        std::ifstream vk_in("verification_key.bin", std::ios::binary); vk_in >> vk; vk_in.close();
    }

    // =================================================================
    // 🚀 4. 증명 생성 (Prover) & 타이머 측정
    // =================================================================
    cout << "\n🚀 증명 생성 (Prover) 시작 - 데이터 덤프 진행됨!" << endl;
    auto prover_start = std::chrono::high_resolution_clock::now();

    auto proof = r1cs_gg_ppzksnark_prover<libff::alt_bn128_pp>(
        pk, 
        primary_input, 
        auxiliary_input
    );

    auto prover_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> prover_time = prover_end - prover_start;

    // =================================================================
    // 🛡️ 5. 최종 검증 (Verifier)
    // =================================================================
    cout << "\n🛡️ 검증 (Verifier) 시작..." << endl;
    bool is_valid = r1cs_gg_ppzksnark_verifier_strong_IC<libff::alt_bn128_pp>(
        vk, 
        primary_input, 
        proof
    );

    cout << "===================================================" << endl;
    cout << "📊 [Performance & Result Report]" << endl;
    cout << "⏱️ 전체 Prover 수행 시간 : " << prover_time.count() << " 초" << endl;
    cout << "✅ 최종 검증 결과        : " << (is_valid ? "True (성공)" : "False (실패)") << endl;
    cout << "===================================================" << endl;

    return 0;
}
