import os

# ==============================================================================
# 1. 시스템 파라미터 (BN254)
# ==============================================================================
P = 21888242871839275222246405745257275088696311157297823662689037894645226208583
R = 21888242871839275222246405745257275088548364400416034343698204186575808495617

R255_P = pow(2, 255, P)
INV_R255_P = pow(R255_P, -1, P)

# 하드웨어 야코비안 항등원 (Z=0 이면 항등원)
HW_INFINITY = (0, 0, 0)

# 일반 대수 아핀 항등원
POINT_AT_INFINITY = (0, 0)

# ==============================================================================
# 2. 하드웨어 내부 로직 모사 (Montgomery Domain Jacobian Mixed Add)
# ==============================================================================
def mon_mul(A, B):
    """FPGA 내부의 몽고메리 모듈러 곱셈기 (A * B * R^-1) 모사"""
    return (A * B * INV_R255_P) % P

def hw_mixed_add(bucket, hw_point):
    """
    [입력]
      - bucket: (X1, Y1, Z1) 버킷의 야코비안 상태 (HW 몽고메리 도메인)
      - hw_point: (X2, Y2) 들어오는 아핀 점 (Z2는 암묵적으로 1, 즉 HW상 R255_P)
    """
    if bucket == HW_INFINITY:
        # 버킷이 비어있다면 아핀 점을 야코비안으로 승격. Z에 몽고메리의 '1'인 R255_P를 부여.
        return (hw_point[0], hw_point[1], R255_P)

    X1, Y1, Z1 = bucket
    X2, Y2 = hw_point

    U1 = X1
    S1 = Y1

    Z1_sq = mon_mul(Z1, Z1)
    U2 = mon_mul(X2, Z1_sq)

    Z1_cu = mon_mul(Z1_sq, Z1)
    S2 = mon_mul(Y2, Z1_cu)

    # Point Doubling
    if U1 == U2:
        if S1 == S2:
            X1_sq = mon_mul(X1, X1)
            M = (3 * X1_sq) % P
            Y1_sq = mon_mul(Y1, Y1)
            S = (4 * mon_mul(X1, Y1_sq)) % P
            
            M_sq = mon_mul(M, M)
            X3 = (M_sq - 2*S) % P
            
            Y1_4 = mon_mul(Y1_sq, Y1_sq)
            Y3 = (mon_mul(M, (S - X3) % P) - 8*Y1_4) % P
            Z3 = (2 * mon_mul(Y1, Z1)) % P
            return (X3, Y3, Z3)
        else:
            return HW_INFINITY

    # Point Addition
    H = (U2 - U1) % P
    R_val = (S2 - S1) % P

    H_sq = mon_mul(H, H)
    H_cu = mon_mul(H_sq, H)
    R_sq = mon_mul(R_val, R_val)
    U1_H_sq = mon_mul(U1, H_sq)

    X3 = (R_sq - H_cu - 2*U1_H_sq) % P
    Y3 = (mon_mul(R_val, (U1_H_sq - X3) % P) - mon_mul(S1, H_cu)) % P
    Z3 = mon_mul(Z1, H)

    return (X3, Y3, Z3)

# (소프트웨어 리덕션 용 일반 대수 아핀 덧셈)
def ec_add_affine(p1, p2):
    if p1 == POINT_AT_INFINITY: return p2
    if p2 == POINT_AT_INFINITY: return p1
    x1, y1 = p1; x2, y2 = p2
    if x1 == x2 and y1 == y2:
        if y1 == 0: return POINT_AT_INFINITY
        lam = (3 * x1 * x1) * pow(2 * y1, -1, P) % P
    elif x1 == x2: return POINT_AT_INFINITY
    else: lam = (y2 - y1) * pow(x2 - x1, -1, P) % P
    x3 = (lam * lam - x1 - x2) % P
    y3 = (lam * (x1 - x3) - y1) % P
    return (x3, y3)

# ==============================================================================
# 메인 검증 파이프라인
# ==============================================================================
print("1. C++ 인풋 파일 로드 중...")
scalars, hw_points = [], []
try:
    with open("msm_input_scalar.txt") as fs, open("msm_input_x.txt") as fx, open("msm_input_y.txt") as fy:
        for s_line, x_line, y_line in zip(fs, fx, fy):
            # 🚨 스칼라는 C++에서 몽고메리 변환 없이 순수 비트로 덤프했으므로 파이썬도 그대로 읽음
            s = int(s_line.strip())
            
            # 하드웨어 포인트는 몽고메리 가중치가 포함된 로우 데이터 그대로 사용
            hw_x = int(x_line.strip())
            hw_y = int(y_line.strip())
            
            scalars.append(s)
            hw_points.append((hw_x, hw_y))
except FileNotFoundError:
    print("[오류] C++ 인풋 파일을 찾을 수 없습니다. C++ 코드를 먼저 실행하세요.")
    exit(1)

print(f"  -> {len(scalars)}개 데이터 인입 완료.")
print("2. [RTL 모사] FPGA 야코비안 버킷 누적(Mixed Add) 시뮬레이션 중...")
WINDOW_COUNT, BITS_PER_WINDOW, BUCKETS_PER_WINDOW = 24, 11, 2048
hw_buckets = [[HW_INFINITY for _ in range(BUCKETS_PER_WINDOW)] for _ in range(WINDOW_COUNT)]

for i in range(len(scalars)):
    s, hw_pt = scalars[i], hw_points[i]
    for w in range(WINDOW_COUNT):
        bucket_idx = (s >> (w * BITS_PER_WINDOW)) & 0x7FF
        if bucket_idx == 0: continue
        # 하드웨어 내부 곱셈기와 100% 동일한 연산 트리거
        hw_buckets[w][bucket_idx] = hw_mixed_add(hw_buckets[w][bucket_idx], hw_pt)

print("3. C++ 전달용 버킷 텍스트 파일(X, Y, Z) 덤프 중...")
with open("buckets_X.txt", "w") as fx, open("buckets_Y.txt", "w") as fy, open("buckets_Z.txt", "w") as fz:
    for w in range(WINDOW_COUNT):
        for b in range(BUCKETS_PER_WINDOW):
            pt = hw_buckets[w][b]
            # hw_buckets 내부의 점은 이미 2^255 도메인 내부에 존재하므로, 변환 없이 그대로 출력!
            fx.write(f"{pt[0]}\n")
            fy.write(f"{pt[1]}\n")
            fz.write(f"{pt[2]}\n")

print("4. Python 소프트웨어 도메인 복원 및 리덕션 연산 중...")
# 하드웨어 버킷을 일반 수학 아핀 도메인으로 복원 (C++의 batch_to_special 과 동일 역할)
affine_buckets = [[POINT_AT_INFINITY for _ in range(BUCKETS_PER_WINDOW)] for _ in range(WINDOW_COUNT)]
for w in range(WINDOW_COUNT):
    for b in range(BUCKETS_PER_WINDOW):
        if hw_buckets[w][b] != HW_INFINITY:
            X_hw, Y_hw, Z_hw = hw_buckets[w][b]
            # 1단계: 몽고메리 가중치 상쇄 (2^255의 역원 곱셈)
            X = (X_hw * INV_R255_P) % P
            Y = (Y_hw * INV_R255_P) % P
            Z = (Z_hw * INV_R255_P) % P
            # 2단계: 야코비안 -> 아핀 좌표 변환 (Z^-2, Z^-3)
            Z_inv = pow(Z, -1, P)
            Z_inv_sq = (Z_inv * Z_inv) % P
            Z_inv_cu = (Z_inv_sq * Z_inv) % P
            affine_buckets[w][b] = ((X * Z_inv_sq) % P, (Y * Z_inv_cu) % P)

window_results = [POINT_AT_INFINITY] * WINDOW_COUNT
for w in range(WINDOW_COUNT):
    running_sum = window_sum = POINT_AT_INFINITY
    for b in range(BUCKETS_PER_WINDOW - 1, 0, -1):
        if affine_buckets[w][b] != POINT_AT_INFINITY:
            running_sum = ec_add_affine(running_sum, affine_buckets[w][b])
        window_sum = ec_add_affine(window_sum, running_sum)
    window_results[w] = window_sum

final_result = POINT_AT_INFINITY
for w in range(WINDOW_COUNT - 1, -1, -1):
    if final_result != POINT_AT_INFINITY:
        for _ in range(BITS_PER_WINDOW):
            final_result = ec_add_affine(final_result, final_result)
    final_result = ec_add_affine(final_result, window_results[w])

print("\n==============================================")
print(" 🏁 [Python 골든 모델] 최종 산출 좌표")
if final_result == POINT_AT_INFINITY:
    print(" X: Point at Infinity")
else:
    print(f" X: {final_result[0]}")
    print(f" Y: {final_result[1]}")
print("==============================================\n")