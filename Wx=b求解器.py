import re
import copy

def safe_eval(expr_str, var_map=None):
    globals_dict = {'__builtins__': None}
    local_map = var_map if var_map is not None else {}
    try:
        result = eval(expr_str, globals_dict, local_map)
        if isinstance(result, (int, float)):
            return float(result)
        raise ValueError(f"非数字结果: {result}")
    except Exception as e:
        raise ValueError(f"求值失败: {expr_str}\n错误: {e}")

def substitute(expr, known):
    if not known:
        return expr
    pattern = r'\b(' + '|'.join(re.escape(v) for v in known) + r')\b'
    return re.sub(pattern, lambda m: str(known[m.group(1)]), expr)

def substitute_eqs(eqs, known):
    return [substitute(eq, known) for eq in eqs]

def extract_vars(expr):
    return set(re.findall(r'\b[a-zA-Z_][a-zA-Z0-9_]*\b', expr))

def filter_identities(eqs, known):
    filtered = []
    for eq in eqs:
        if '=' not in eq:
            filtered.append(eq)
            continue
        left, right = eq.split('=', 1)
        if not (extract_vars(left) | extract_vars(right)) - set(known):
            try:
                val_left = safe_eval(left, known)
                val_right = safe_eval(right, known)
                if abs(val_left - val_right) > 1e-9:
                    raise ValueError(f"矛盾方程: {eq}")
                continue
            except:
                filtered.append(eq)
        else:
            filtered.append(eq)
    return filtered

def parse_linear(eq_str, var_names):
    if '=' not in eq_str:
        raise ValueError("缺少等号")
    left, right = eq_str.split('=', 1)
    expr = f"({left}) - ({right})"
    zero_map = {v: 0.0 for v in var_names}
    expr_zero = safe_eval(expr, zero_map)
    coeffs = []
    for i, var in enumerate(var_names):
        unit_map = {v: 1.0 if v == var else 0.0 for v in var_names}
        expr_one = safe_eval(expr, unit_map)
        coeff = expr_one - expr_zero
        coeffs.append(coeff)
        double_map = {v: 2.0 if v == var else 0.0 for v in var_names}
        expr_two = safe_eval(expr, double_map)
        if abs(expr_two - (2 * coeff + expr_zero)) > 1e-9:
            raise ValueError(f"非线性: {var}")
    for i in range(len(var_names)):
        for j in range(i+1, len(var_names)):
            both_map = {v: 1.0 if v in (var_names[i], var_names[j]) else 0.0 for v in var_names}
            expr_both = safe_eval(expr, both_map)
            expected = expr_zero + coeffs[i] + coeffs[j]
            if abs(expr_both - expected) > 1e-9:
                raise ValueError(f"交叉项: {var_names[i]}, {var_names[j]}")
    const = -expr_zero
    return coeffs, const

def gauss_solve(W, b):
    n = len(W)
    A = copy.deepcopy(W)
    B = copy.deepcopy(b)
    for i in range(n):
        max_row = max(range(i, n), key=lambda r: abs(A[r][i]))
        if abs(A[max_row][i]) < 1e-12:
            raise ValueError("奇异矩阵")
        if max_row != i:
            A[i], A[max_row] = A[max_row], A[i]
            B[i], B[max_row] = B[max_row], B[i]
        for j in range(i+1, n):
            factor = A[j][i] / A[i][i]
            for k in range(i, n):
                A[j][k] -= factor * A[i][k]
            B[j] -= factor * B[i]
    x = [0.0] * n
    for i in range(n-1, -1, -1):
        s = sum(A[i][j] * x[j] for j in range(i+1, n))
        x[i] = (B[i] - s) / A[i][i]
    return x

def solve_direct_assignments(eqs, known):
    new_fixed = {}
    remaining = []
    for eq in eqs:
        if '=' not in eq:
            remaining.append(eq)
            continue
        left, right = eq.split('=', 1)
        left = left.strip()
        if re.fullmatch(r'[a-zA-Z_][a-zA-Z0-9_]*', left):
            vars_in_right = extract_vars(right)
            if not vars_in_right - set(known):
                try:
                    val = safe_eval(right, known)
                    new_fixed[left] = val
                except:
                    remaining.append(eq)
            else:
                remaining.append(eq)
        else:
            remaining.append(eq)
    return new_fixed, remaining


def solve(eqs, init_known=None):
    eqs = list(dict.fromkeys([e.strip() for e in eqs if e.strip()]))
    known = init_known.copy() if init_known else {}
    remaining = eqs

    iteration = 0
    max_iter = 20
    while remaining and iteration < max_iter:
        iteration += 1

        remaining = substitute_eqs(remaining, known)
        remaining = filter_identities(remaining, known)

        # 循环处理直接赋值（右边无未知变量）
        while True:
            new_fixed, remaining = solve_direct_assignments(remaining, known)
            if not new_fixed:
                break
            known.update(new_fixed)
            remaining = substitute_eqs(remaining, new_fixed)
            remaining = filter_identities(remaining, known)

        if not remaining:
            break

        all_vars = set()
        for eq in remaining:
            all_vars |= extract_vars(eq)
        unknowns = sorted(all_vars - set(known))
        if not unknowns:
            if remaining:
                raise ValueError(f"矛盾方程: {remaining}")
            break

        linear_eqs = []
        nonlinear_eqs = []
        for eq in remaining:
            try:
                parse_linear(eq, unknowns)
                linear_eqs.append(eq)
            except ValueError:
                nonlinear_eqs.append(eq)

        linear_vars = set()
        for eq in linear_eqs:
            linear_vars |= extract_vars(eq)
        linear_vars = sorted(linear_vars - set(known))

        # 情况1：方程数等于变量数 -> 直接求解
        if len(linear_eqs) == len(linear_vars) and linear_vars:
            W, b = [], []
            for eq in linear_eqs:
                coeffs, const = parse_linear(eq, linear_vars)
                W.append(coeffs)
                b.append(const)
            try:
                x = gauss_solve(W, b)
            except ValueError as e:
                raise ValueError(f"线性系统求解失败: {e}")

            print("\n=== 线性系统求解 ===")
            print("线性变量:", linear_vars)
            print("系数矩阵 W:")
            for row in W:
                print([round(v, 6) for v in row])
            print("右端向量 b:", [round(v, 6) for v in b])
            print("解向量:")
            for var, val in zip(linear_vars, x):
                print(f"   {var} = {val:.8f}")

            for var, val in zip(linear_vars, x):
                known[var] = val
            continue

        # 情况2：方程数多于变量数 -> 选择一组线性无关的方程求解，并验证其余方程
        elif len(linear_eqs) > len(linear_vars) and linear_vars:
            # 尝试选出 len(linear_vars) 个线性无关的方程
            selected_eqs = []
            selected_indices = []
            # 使用高斯消元思想选取主元
            A_rows = []
            b_vals = []
            for idx, eq in enumerate(linear_eqs):
                coeffs, const = parse_linear(eq, linear_vars)
                row = coeffs[:]
                # 检查当前行是否与已选行线性无关
                # 构造临时增广矩阵，对已选行做消元
                temp_A = [A_rows[i][:] for i in range(len(selected_eqs))]
                temp_b = b_vals[:]
                temp_A.append(row)
                temp_b.append(const)
                # 尝试消元
                nvar = len(linear_vars)
                # 对前 len(selected_eqs) 行已经消元，只需对新行进行消元并检查主元
                # 简化：直接调用一个函数判断加入后是否满秩
                def is_independent(new_row, new_b, selected_rows, selected_bs):
                    m = len(selected_rows)
                    # 构建增广矩阵
                    aug = [selected_rows[i] + [selected_bs[i]] for i in range(m)]
                    aug.append(new_row + [new_b])
                    # 高斯消元求秩
                    mat = [row[:] for row in aug]
                    ncol = len(new_row) + 1
                    rank = 0
                    for col in range(ncol-1):  # 对系数列
                        pivot = None
                        for r in range(rank, len(mat)):
                            if abs(mat[r][col]) > 1e-12:
                                pivot = r
                                break
                        if pivot is None:
                            continue
                        mat[rank], mat[pivot] = mat[pivot], mat[rank]
                        pivot_val = mat[rank][col]
                        for c in range(col, ncol):
                            mat[rank][c] /= pivot_val
                        for r in range(len(mat)):
                            if r != rank and abs(mat[r][col]) > 1e-12:
                                factor = mat[r][col]
                                for c in range(col, ncol):
                                    mat[r][c] -= factor * mat[rank][c]
                        rank += 1
                        if rank == nvar:
                            break
                    # 检查增广矩阵的秩是否等于系数矩阵的秩
                    # 简化：检查最后一行是否全零（系数全零但常数非零则矛盾）
                    # 这里我们只关心是否线性无关，即新行不能由已选行线性表示
                    # 如果加入后秩增加，则无关
                    # 计算当前系数矩阵的秩
                    coeff_mat = [row[:nvar] for row in mat[:rank]]
                    # 实际上，增广矩阵消元后，如果新行消元后系数全零且常数非零，则矛盾，但这里我们只是选方程，先不考虑矛盾
                    # 判定线性无关：消元后新行的系数不全为零
                    if rank > len(selected_rows):
                        return True
                    else:
                        return False
                if is_independent(row, const, [A_rows[i] for i in range(len(selected_eqs))], b_vals):
                    selected_eqs.append(eq)
                    selected_indices.append(idx)
                    A_rows.append(row)
                    b_vals.append(const)
                    if len(selected_eqs) == len(linear_vars):
                        break
            if len(selected_eqs) == len(linear_vars):
                # 求解选出的方程
                W_sel = A_rows
                b_sel = b_vals
                try:
                    x = gauss_solve(W_sel, b_sel)
                except ValueError as e:
                    raise ValueError(f"选取的线性系统求解失败: {e}")
                # 验证剩余线性方程
                all_linear_eqs = linear_eqs
                solved_vars = dict(zip(linear_vars, x))
                all_ok = True
                for eq in all_linear_eqs:
                    # 将解代入方程左-右，应接近0
                    left, right = eq.split('=', 1)
                    expr = f"({left}) - ({right})"
                    # 代入已知变量（包括新解）
                    test_known = {**known, **solved_vars}
                    try:
                        val = safe_eval(expr, test_known)
                        if abs(val) > 1e-9:
                            print(f"验证失败: {eq} 差值 = {val}")
                            all_ok = False
                            break
                    except:
                        all_ok = False
                        break
                if all_ok:
                    print("\n=== 线性系统求解（超定，选择满秩子集） ===")
                    print("线性变量:", linear_vars)
                    print("选用的方程:")
                    for eq in selected_eqs:
                        print(f"   {eq}")
                    print("系数矩阵 W:")
                    for row in W_sel:
                        print([round(v, 6) for v in row])
                    print("右端向量 b:", [round(v, 6) for v in b_sel])
                    print("解向量:")
                    for var, val in zip(linear_vars, x):
                        print(f"   {var} = {val:.8f}")
                    for var, val in zip(linear_vars, x):
                        known[var] = val
                    continue
                else:
                    print("\n⚠️ 超定线性系统矛盾，无法求解")
                    print("当前已知变量:", known)
                    print("剩余方程:", remaining)
                    print("未解出的变量:", unknowns)
                    return
            else:
                # 无法选出满秩子集，系数矩阵秩不足
                print("\n⚠️ 线性系统秩不足，无法求解")
                print("当前已知变量:", known)
                print("剩余方程:", remaining)
                print("未解出的变量:", unknowns)
                return

        else:
            # 方程数小于变量数 或 无非线性方程但线性方程数不等于变量数
            print("\n⚠️ 无法继续求解")
            if nonlinear_eqs:
                print("原因：存在非线性方程，且线性方程数不等于线性变量数")
                print("非线性方程:")
                for eq in nonlinear_eqs:
                    print(f"   {eq}")
            else:
                print(f"原因：线性方程数 ({len(linear_eqs)}) 不等于线性变量数 ({len(linear_vars)})，系统欠定或超定")
                print("线性方程:")
                for eq in linear_eqs:
                    print(f"   {eq}")
            print("当前已知变量:", known)
            print("剩余方程:", remaining)
            print("未解出的变量:", unknowns)
            return

    if not remaining:
        print("\n✅ 所有方程求解完毕")
    else:
        print("\n⚠️ 仍有剩余方程未处理")
        print("剩余方程:", remaining)
        unknowns = sorted(extract_vars(" ".join(remaining)) - set(known))
        if unknowns:
            print("未解出的变量:", unknowns)
        else:
            print("但无未知变量，可能为矛盾恒等式")

    print("\n最终解向量:")
    for var in sorted(known.keys()):
        print(f"   {var} = {known[var]:.8f}")

if __name__ == "__main__":
    # 用户提供的初始已知变量
    init_known = {
        'v1_':15,
        'v3_':-10,
    }
    eqs = [
        "v2_ = 4*i0",
        "(v1_-v1)/20 = (v1-v2)/5+(v1-v3)/10 + 3",
        "i0 = (v1-v3)/10",
        "(v1-v2)/5 = (v2 - v2_)/5 + (v2-v3)/5",
        "(v1-v3)/10 +3+(v2-v3)/5 = (v3-v3_)/15"
    ]
    solve(eqs, init_known)