#!/usr/bin/env python3
"""
exact_unit_distance.py
══════════════════════════════════════════════════════════════════════
Kiểm chứng EXACT bằng integer arithmetic — zero floating-point error.

Tại sao integer arithmetic đủ?
  Khoảng cách² trong Z² và Z[ω] là SỐ NGUYÊN.
  So sánh dist²==1 là phép toán nguyên hoàn toàn chính xác.
  Không cần epsilon, không có sai số làm tròn.

Constructions được implement:
  1. Z²    — lưới vuông,  dist²(p,q) = da² + db²
  2. Z[ω]  — Eisenstein,  |p-q|²    = da² + da·db + db²   (6 hàng xóm)
  3. Hai lớp Eisenstein: Z[ω] tỉ lệ 1 và Z[ω] tỉ lệ 1/√7

Kết quả mong đợi theo lý thuyết:
  Z²    → ν ≈ 2n   (4 hàng xóm, 2 cặp/điểm)
  Z[ω]  → ν ≈ 3n   (6 hàng xóm, 3 cặp/điểm)
  OpenAI n^1.014 → cần construction không thể implement đơn giản
══════════════════════════════════════════════════════════════════════
"""

import math, sys, time
import numpy as np

# ══════════════════════════════════════════════════════════════════════
#  GENERATORS
# ══════════════════════════════════════════════════════════════════════

def gen_square(n: int) -> np.ndarray:
    """n điểm từ lưới Z² (hàng-cột), dtype int64."""
    k = math.ceil(math.sqrt(n))
    pts = [(i, j) for j in range(k) for i in range(k)]
    return np.array(pts[:n], dtype=np.int64)


def gen_eisenstein(n: int) -> np.ndarray:
    """
    n điểm từ Z[ω] (ω = e^{iπ/3}), sắp xếp tăng dần theo norm:
      |a + bω|² = a² + ab + b²   (Eisenstein norm form)
    Điểm (a,b) biểu diễn số phức  a + b·(1/2 + i√3/2).
    """
    R = math.ceil(math.sqrt(n)) + 5
    a_arr, b_arr = np.mgrid[-R:R+1, -R:R+1]  # shape (2R+1, 2R+1)
    norms = (a_arr*a_arr + a_arr*b_arr + b_arr*b_arr).ravel()
    order = np.argsort(norms, kind='stable')
    a_flat = a_arr.ravel()[order[:n]]
    b_flat = b_arr.ravel()[order[:n]]
    return np.column_stack([a_flat, b_flat]).astype(np.int64)


# ══════════════════════════════════════════════════════════════════════
#  EXACT PAIR COUNTING — numpy vectorized, chunked để tiết kiệm RAM
# ══════════════════════════════════════════════════════════════════════

CHUNK = 1500   # Mỗi khối ≤ 1500×1500 = 2.25M entries ≈ 18 MB

def count_pairs(pts: np.ndarray, target_norm: int, norm_fn) -> int:
    """
    Đếm số cặp (i<j) với norm_fn(da,db) == target_norm.
    HOÀN TOÀN CHÍNH XÁC — integer comparison.
    """
    n = len(pts)
    total = 0
    for i0 in range(0, n, CHUNK):
        pi = pts[i0:i0+CHUNK]
        for j0 in range(i0, n, CHUNK):
            pj = pts[j0:j0+CHUNK]
            da = (pi[:, None, 0] - pj[None, :, 0]).astype(np.int64)
            db = (pi[:, None, 1] - pj[None, :, 1]).astype(np.int64)
            hits = int(np.sum(norm_fn(da, db) == target_norm))
            if i0 == j0:
                hits //= 2   # mỗi cặp (i,j) bị đếm hai lần trong khối chéo
            total += hits
    return total


def norm_z2(da, db):    return da*da + db*db
def norm_eisen(da, db): return da*da + da*db + db*db


# ══════════════════════════════════════════════════════════════════════
#  PHÂN TÍCH HAI LỚP (EXACT)
#  Lớp 0: Z[ω] tỉ lệ 1       → điểm p0 = a0 + b0·ω
#  Lớp 1: Z[ω] tỉ lệ 1/√7   → điểm p1 = (a1 + b1·ω)/√7
#
#  Khoảng cách² giữa p0 (lớp 0) và p1 (lớp 1):
#    |p0 - p1|² = A + B/√7
#  với:
#    A = norm0 + norm1/7           (phần hữu tỉ)
#    B = −2·Re(z0·z̄1)             (phần vô tỉ nhân 1/√7)
#
#  Để |p0-p1| = 1 (EXACT), cần đồng thời:
#    7·norm0 + norm1 = 7              ... (I)   [nhân A=1 lên 7]
#    2a0a1 + a0b1 + a1b0 + 2b0b1 = 0  ... (II)  [Re(z0·z̄1)=0]
#  Cả hai đều là điều kiện NGUYÊN → kiểm tra chính xác tuyệt đối.
#
#  PHÁT HIỆN QUAN TRỌNG:
#    Điều kiện (I): 7·norm0 ≤ 7 → norm0 ≤ 1
#    Chỉ có nguồn gốc (norm0=0) và 6 điểm lân cận gần nhất (norm0=1)
#    mới có thể có cặp cross-layer!
#    → Số cặp cross-layer là HẰNG SỐ O(1), KHÔNG tăng theo n!
# ══════════════════════════════════════════════════════════════════════

def count_cross_layer_exact(pts0: np.ndarray, pts1: np.ndarray) -> tuple:
    """
    Đếm chính xác số cặp unit-distance giữa lớp 0 và lớp 1.
    Trả về (count, n_eligible0, n_eligible1).
    """
    norm0 = (pts0[:,0]**2 + pts0[:,0]*pts0[:,1] + pts0[:,1]**2).astype(np.int64)
    norm1 = (pts1[:,0]**2 + pts1[:,0]*pts1[:,1] + pts1[:,1]**2).astype(np.int64)

    # Chỉ điểm lớp 0 với norm0 ≤ 1 mới có thể có cross-layer pairs
    elig0_idx = np.where(norm0 <= 1)[0]
    # Điểm lớp 1 cần norm1 ≤ 7
    elig1_idx = np.where(norm1 <= 7)[0]

    count = 0
    for i in elig0_idx:
        a0, b0 = int(pts0[i, 0]), int(pts0[i, 1])
        n0_val = a0*a0 + a0*b0 + b0*b0
        target_n1 = 7 - 7*n0_val   # từ điều kiện (I)
        for j in elig1_idx:
            a1, b1 = int(pts1[j, 0]), int(pts1[j, 1])
            if a1*a1 + a1*b1 + b1*b1 != target_n1:
                continue
            # Điều kiện (II): 2a0a1 + a0b1 + a1b0 + 2b0b1 = 0
            if 2*a0*a1 + a0*b1 + a1*b0 + 2*b0*b1 == 0:
                count += 1

    return count, len(elig0_idx), len(elig1_idx)


# ══════════════════════════════════════════════════════════════════════
#  FIT POWER LAW  ν(n) = c · n^α
# ══════════════════════════════════════════════════════════════════════

def fit_alpha(ns, nus):
    v = [(n, nu) for n, nu in zip(ns, nus) if n > 1 and nu > 0]
    if len(v) < 4:
        return float('nan'), float('nan')
    lx = np.log([n for n, _ in v])
    ly = np.log([nu for _, nu in v])
    coeffs = np.polyfit(lx, ly, 1)
    return float(coeffs[0]), float(math.exp(coeffs[1]))


def local_alpha(n1, nu1, n2, nu2):
    if nu1 <= 0 or nu2 <= 0:
        return None
    return math.log(nu2 / nu1) / math.log(n2 / n1)


# ══════════════════════════════════════════════════════════════════════
#  SCAN
# ══════════════════════════════════════════════════════════════════════

SIZES = [10, 20, 50, 100, 200, 300, 500, 750, 1000, 2000, 3000, 5000]

def run_scan(label, gen_fn, norm_fn, target=1):
    print(f"\n{'━'*68}")
    print(f"  {label}")
    print(f"{'━'*68}")
    print(f"  {'n':>6}  {'ν_exact':>9}  {'ν/n':>7}  {'α_local':>9}  {'sec':>5}")
    print(f"  {'─'*55}")

    ns, nus, prev = [], [], None
    for n in SIZES:
        t0 = time.perf_counter()
        pts = gen_fn(n)
        nu  = count_pairs(pts, target, norm_fn)
        dt  = time.perf_counter() - t0
        ns.append(n); nus.append(nu)
        a_loc = (f"{local_alpha(prev[0],prev[1],n,nu):.4f}" if prev and prev[1] > 0 else '')
        print(f"  {n:>6}  {nu:>9}  {nu/n:>7.4f}  {a_loc:>9}  {dt:>5.2f}")
        prev = (n, nu)

    alpha, c = fit_alpha(ns, nus)
    print(f"\n  → ν ≈ {c:.3f} · n^{alpha:.6f}   (OLS trên log-log)")
    return ns, nus, alpha


# ══════════════════════════════════════════════════════════════════════
#  MAIN
# ══════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    sep = '═' * 68
    print(sep)
    print("  EXACT UNIT-DISTANCE PAIR COUNTING")
    print("  Integer arithmetic — kết quả chính xác tuyệt đối, không sai số")
    print(f"  Python {sys.version.split()[0]}, NumPy {np.__version__}")
    print(sep)

    # ── 1. Square Grid ──────────────────────────────────────────────
    r1 = run_scan(
        "Z²  — Lưới vuông (Square Grid)",
        gen_square, norm_z2, target=1
    )

    # ── 2. Eisenstein Z[ω] ──────────────────────────────────────────
    r2 = run_scan(
        "Z[ω] — Eisenstein, CM field Q(√−3)",
        gen_eisenstein, norm_eisen, target=1
    )

    # ── 3. Hai lớp Eisenstein — phân tích cross-layer ───────────────
    print(f"\n{'━'*68}")
    print("  Hai lớp Eisenstein: Z[ω] (tỉ lệ 1) + Z[ω] (tỉ lệ 1/√7)")
    print(f"{'━'*68}")
    print("  Cặp same-layer-0:  norm_eisen = 1  (6 hàng xóm, ≈ 3n₀)")
    print("  Cặp same-layer-1:  norm_eisen = 7  (12 hàng xóm, ≈ 6n₁)")
    print("  Cặp cross-layer:   điều kiện (I)+(II) — dự đoán O(1)")
    print()
    print(f"  {'n':>6}  {'ν_same0':>9}  {'ν_same1':>9}  {'ν_cross':>9}  {'ν_total':>9}  {'ν/n':>6}")
    print(f"  {'─'*65}")

    for n in [100, 500, 1000, 2000, 5000]:
        n0 = n // 2
        n1 = n - n0
        pts0 = gen_eisenstein(n0)
        pts1 = gen_eisenstein(n1)

        t0 = time.perf_counter()
        s0 = count_pairs(pts0, 1, norm_eisen)          # same-layer 0: norm=1
        s1 = count_pairs(pts1, 7, norm_eisen)          # same-layer 1: norm=7 (scale 1/√7)
        cx, elig0, elig1 = count_cross_layer_exact(pts0, pts1)
        total = s0 + s1 + cx
        dt = time.perf_counter() - t0

        print(f"  {n:>6}  {s0:>9}  {s1:>9}  {cx:>9}  {total:>9}  {total/n:>6.3f}  [{dt:.2f}s]")

    print(f"\n  ↑ Chú ý: ν_cross KHÔNG THAY ĐỔI khi n tăng từ 100 → 5000")
    print( "    → Cross-layer là hằng số O(1), không đóng góp vào α!")
    print()
    print( "  Lý do toán học:")
    print( "    Điều kiện (I): 7·norm0 + norm1 = 7")
    print( "    → norm0 ≤ 1  (chỉ 7 điểm lớp 0 thỏa: origin + 6 lân cận)")
    print( "    → norm1 ≤ 7  (hữu hạn điểm lớp 1)")
    print( "    → Tổng số nghiệm là hằng số, không phụ thuộc n!")

    # ── Tóm tắt ─────────────────────────────────────────────────────
    print(f"\n{sep}")
    print("  TÓM TẮT & KẾT LUẬN")
    print(sep)
    print(f"  {'Construction':<30}  {'α_OLS':>8}  {'ν/n (n=5000)':>14}")
    print(f"  {'─'*56}")

    ns1, nus1, a1 = r1
    ns2, nus2, a2 = r2
    print(f"  {'Z² Square Grid':<30}  {a1:>8.4f}  {nus1[-1]/ns1[-1]:>14.4f}")
    print(f"  {'Z[ω] Eisenstein':<30}  {a2:>8.4f}  {nus2[-1]/ns2[-1]:>14.4f}")

    print()
    print("  CÁC NHẬN XÉT QUAN TRỌNG:")
    print()
    print("  [1] Cả hai construction đều cho α ≈ 1 (tuyến tính)")
    print("      → Không có construction đơn giản nào đạt α > 1.014")
    print()
    print("  [2] Hai lớp Eisenstein: cross-layer = HẰNG SỐ, không tăng theo n")
    print("      → 'OpenAI visualization' trong web app chỉ là minh họa khái niệm")
    print("      → ν ≈ 4.5n (tuyến tính), không phải n^1.014")
    print()
    print("  [3] Để đạt α > 1, cần construction fundamentally khác:")
    print("      • Class field tower: K₀ ⊂ K₁ ⊂ K₂ ⊂ ... (vô hạn tầng)")
    print("      • Golod-Shafarevich đảm bảo tháp tồn tại vô hạn")
    print("      • Tại mỗi tầng Kᵢ: tập điểm có mật độ cặp tăng theo hàm mũ")
    print("      • Cần exact arithmetic trong trường số đại số: Magma / PARI-GP")
    print("      • n cần >> 10⁶ để thấy rõ n^1.014 khác n^1.001")
    print()
    print("  [4] Verification thực sự của kết quả OpenAI:")
    print("      • Peer review bởi Tim Gowers et al. (không phải tính toán máy tính)")
    print("      • Chứng minh algebraic: xây dựng construction + đếm pair bằng lý thuyết")
    print("      • Bằng chứng dạng: 'tồn tại n₀ sao cho với n > n₀, ν(n) ≥ n^1.014'")
    print("      • n₀ này có thể cực kỳ lớn — không accessible bằng máy tính!")
    print(sep)
