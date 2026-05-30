import os

# ==============================================================================
# 1. 타원곡선 파라미터 세팅 (BN254 / alt_bn128)
# ==============================================================================
# Fp (베이스 필드, 좌표계) 모듈러스
P = 21888242871839275222246405745257275088696311157297823662689037894645226208583
# Fr (스칼라 필드) 모듈러스
R = 21888242871839275222246405745257275088548364400416034343698204186575808495617

# 2^255 몽고메리 도메인 상쇄를 위한 역원 계산
INV_R255_P = pow(pow(2, 255, P), -1, P)
INV_R255_R = pow(pow(2, 255, R), -1, R)

# 항등원(Point at Infinity) 정의
POINT_AT_INFINITY = (None, None)

# ==============================================================================
# 2. 타원곡선(Elliptic Curve) 연산 코어
# ==============================================================================
def ec_add(p1, p2):
    """표준 아핀(Affine) 좌표계 타원곡선 덧셈 (y^2 = x^3 + 3)"""
    if p1 == POINT_AT_INFINITY: return p2
    if p2 == POINT_AT_INFINITY: return p1
    
    x1, y1 = p1
    x2, y2 = p2
    
    if x1 == x2 and y1 == y2:
        # Point Doubling (점 두배 연산)
        if y1 == 0: return POINT_AT_INFINITY
        lam = (3 * x1 * x1) * pow(2 * y1, -1, P) % P
    elif x1 == x2:
        # x좌표가 같고 부호가 다르면 항등원
        return POINT_AT_INFINITY
    else:
        # Point Addition (서로 다른 두 점 덧셈)
        lam = (y2 - y1) * pow(x2 - x1, -1, P) % P
        
    x3 = (lam * lam - x1 - x2) % P
    y3 = (lam * (x1 - x3) - y1) % P
    return (x3, y3)

# ==============================================================================
# 3. 파일 파싱 및 도메인 복원 (Host -> FPGA 인입 모사)
# ==============================================================================
def load_and_decode_inputs():
    print("[1] 테스트 벡터 파일 로드 및 2^255 도메인 디코딩 중...")
    scalars, points = [], []
    
    try:
        with open("msm_input_scalar.txt", "r") as fs, \
             open("msm_input_x.txt", "r") as fx, \
             open("msm_input_y.txt", "r") as fy:
            
            for s_line, x_line, y_line in zip(fs, fx, fy):
                # 텍스트 라인을 파이썬 초거대 정수로 파싱
                raw_s = int(s_line.strip())
                raw_x = int(x_line.strip())
                raw_y = int(y_line.strip())
                
                # 2^255 몽고메리 가중치 상쇄 (실제 대수적 값으로 복원)
                real_s = (raw_s * INV_R255_R) % R
                real_x = (raw_x * INV_R255_P) % P
                real_y = (raw_y * INV_R255_P) % P
                
                scalars.append(real_s)
                points.append((real_x, real_y))
                
        print(f"  -> 총 {len(scalars)}개의 데이터 로드 완료.\n")
        return scalars, points
    except FileNotFoundError:
        print("  -> [오류] 덤프된 txt 파일을 찾을 수 없습니다.")
        exit(1)

# ==============================================================================
# 4. FPGA 버킷 누적 시뮬레이션 (Hardware Logic Model)
# ==============================================================================
def fpga_bucket_accumulation(scalars, points):
    print("[2] FPGA URAM 버킷 누적(Bucket Accumulation) 에뮬레이션 중...")
    WINDOW_COUNT = 24
    BITS_PER_WINDOW = 11
    BUCKETS_PER_WINDOW = 2048 # 2^11
    
    # 24 x 2048 크기의 2차원 버킷 배열 초기화 (항등원으로 채움)
    buckets = [[POINT_AT_INFINITY for _ in range(BUCKETS_PER_WINDOW)] for _ in range(WINDOW_COUNT)]
    
    for i in range(len(scalars)):
        scalar = scalars[i]
        point = points[i]
        
        for w in range(WINDOW_COUNT):
            # 11비트씩 스칼라 슬라이싱 (비트 마스킹)
            bucket_idx = (scalar >> (w * BITS_PER_WINDOW)) & 0x7FF
            
            # 스칼라 조각이 0이면 버킷에 담지 않음 (Pippenger 최적화)
            if bucket_idx == 0:
                continue
                
            # 해당 윈도우의 타겟 버킷에 점 누적(EC Add)
            buckets[w][bucket_idx] = ec_add(buckets[w][bucket_idx], point)
            
    print("  -> FPGA 버킷 누적 파이프라인 처리 완료.\n")
    return buckets

# ==============================================================================
# 5. 호스트 소프트웨어 리덕션 (CPU Logic Model)
# ==============================================================================
def host_software_reduction(buckets):
    print("[3] 호스트 소프트웨어 리덕션 및 최종 병합(Aggregation) 수행 중...")
    WINDOW_COUNT = 24
    BITS_PER_WINDOW = 11
    BUCKETS_PER_WINDOW = 2048
    
    window_results = [POINT_AT_INFINITY] * WINDOW_COUNT
    
    # 1단계: 각 윈도우별 버킷 누적 (역순 순회)
    for w in range(WINDOW_COUNT):
        running_sum = POINT_AT_INFINITY
        window_sum = POINT_AT_INFINITY
        
        # 2047번 버킷부터 1번 버킷까지 역으로 내려오며 누적
        for b in range(BUCKETS_PER_WINDOW - 1, 0, -1):
            if buckets[w][b] != POINT_AT_INFINITY:
                running_sum = ec_add(running_sum, buckets[w][b])
            window_sum = ec_add(window_sum, running_sum)
            
        window_results[w] = window_sum
        
    # 2단계: 호너의 법칙(Horner's Method)을 이용한 윈도우 시프트 & 결합
    final_result = POINT_AT_INFINITY
    
    for w in range(WINDOW_COUNT - 1, -1, -1):
        # 11비트 시프트 (더블링 11회 반복)
        if final_result != POINT_AT_INFINITY:
            for _ in range(BITS_PER_WINDOW):
                final_result = ec_add(final_result, final_result)
                
        # 현재 윈도우 결과 병합
        final_result = ec_add(final_result, window_results[w])
        
    print("  -> 소프트웨어 리덕션 완료.\n")
    return final_result

# ==============================================================================
# 메인 실행부
# ==============================================================================
if __name__ == "__main__":
    print("=== FPGA MSM 하이브리드 파이프라인 Python 골든 모델 ===\n")
    
    # 1. 파일 읽기 및 복원
    scalars, points = load_and_decode_inputs()
    
    # 2. 하드웨어 에뮬레이션
    fpga_buckets = fpga_bucket_accumulation(scalars, points)
    
    # 3. 소프트웨어 에뮬레이션
    final_point = host_software_reduction(fpga_buckets)
    
    # 4. 최종 결과 출력
    print("========================================================")
    print(" 🎉 최종 산출된 G1 MSM 포인트 (Affine X, Y)")
    print("========================================================")
    if final_point == POINT_AT_INFINITY:
        print("결과: Point at Infinity (항등원)")
    else:
        print(f"X : {final_point[0]}")
        print(f"Y : {final_point[1]}")
    print("========================================================")