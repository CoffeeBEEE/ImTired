import re
import copy


# ==================== 安全求值函数（修正：支持纯数值/含变量表达式） ====================
def safe_eval(expr_str, var_map=None):
    """
    安全求值：
    - var_map=None → 仅求值纯数值表达式（无变量）
    - var_map≠None → 代入变量值求值（含变量表达式）
    """
    globals_dict = {'__builtins__': None}
    local_map = var_map if var_map is not None else {}
    try:
        result = eval(expr_str, globals_dict, local_map)
        if isinstance(result, (int, float)):
            return float(result)
        raise ValueError(f"表达式求值结果非数字: {result}")
    except Exception as e:
        raise ValueError(f"求值失败: {expr_str}\n错误: {e}")


# ==================== 提取并链式求解定值变量（核心修正） ====================
def extract_fixed_vars(eqs):
    """
    提取定值变量（支持链式替换，如 i_s=i_1、i_1=0.02 → i_s=0.02）
    返回：{定值变量: 数值}, 剩余非定值方程列表
    """
    fixed_vars = {}
    remaining_eqs = []
    # 匹配「左边单个变量 = 右边表达式」的正则
    single_var_eq_pattern = r'^\s*([a-zA-Z]+(?:_[a-zA-Z0-9]*)?)\s*=\s*([^=]+)\s*$'

    # 第一步：先收集所有「变量=表达式」的方程，尝试迭代求解定值
    var_eq_map = {}  # 暂存「变量: 表达式」
    non_var_eqs = []  # 非「变量=表达式」的方程（如 a + b = c）

    for eq in eqs:
        match = re.match(single_var_eq_pattern, eq)
        if match:
            var_name = match.group(1).strip()
            expr = match.group(2).strip()
            var_eq_map[var_name] = expr
        else:
            non_var_eqs.append(eq)

    # 第二步：迭代求解定值变量（链式替换）
    changed = True
    while changed:
        changed = False
        for var, expr in list(var_eq_map.items()):
            if var in fixed_vars:
                continue  # 已求解，跳过
            try:
                # 尝试代入已有的定值变量，求值表达式
                val = safe_eval(expr, fixed_vars)
                fixed_vars[var] = val
                del var_eq_map[var]  # 从待求解中移除
                changed = True
            except:
                # 表达式含未求解变量，暂存
                continue

    # 第三步：整理剩余方程（未求解的「变量=表达式」 + 非变量方程）
    remaining_eqs = [f"{var}={expr}" for var, expr in var_eq_map.items()] + non_var_eqs
    return fixed_vars, remaining_eqs


# ==================== 代入定值变量（仅替换变量名，不计算表达式） ====================
def substitute_fixed_vars(eqs, fixed_vars):
    """
    仅替换方程中的定值变量名，保留运算结构（如 v_1/R_3 → v_1/200）
    避免错误拼接变量名（如 v_0.005）
    """
    simplified_eqs = []
    # 生成「变量名→数值字符串」的映射（保留浮点数格式）
    var_to_val = {var: str(val) for var, val in fixed_vars.items()}
    # 匹配独立变量名的正则（避免部分匹配，如 R_3 不匹配 R_30）
    var_pattern = re.compile(r'\b(' + '|'.join(re.escape(var) for var in fixed_vars) + r')\b')

    for eq in eqs:
        # 仅替换独立的变量名，不修改运算结构
        simplified_eq = var_pattern.sub(lambda m: var_to_val[m.group(1)], eq)
        simplified_eqs.append(simplified_eq)
    return simplified_eqs


# ==================== 提取所有变量（仅提取非定值变量） ====================
def extract_variables(eqs, fixed_vars):
    """
    提取方程中的变量，排除已定值的变量
    """
    var_set = set()
    var_pattern = r'\b[a-zA-Z]+(?:_[a-zA-Z0-9]*)?\b'
    for eq in eqs:
        eq_stripped = eq.strip()
        if not eq_stripped:
            raise ValueError("检测到空方程，请填写有效的线性方程")
        vars_in_eq = re.findall(var_pattern, eq_stripped)
        if not vars_in_eq:
            raise ValueError(f"方程 {eq} 中未检测到合法变量")
        # 排除定值变量
        var_set.update([v for v in vars_in_eq if v not in fixed_vars])

    if not var_set:
        raise ValueError("未检测到非定值变量")
    return sorted(var_set)


# ==================== 解析单个方程（无修改） ====================
def parse_equation(eq_str, var_names):
    if '=' not in eq_str:
        raise ValueError(f"方程缺少等号: {eq_str}")
    left, right = eq_str.split('=', 1)
    left = left.strip()
    right = right.strip()
    expr_str = f"({left}) - ({right})"

    zero_map = {var: 0.0 for var in var_names}
    expr_zero = safe_eval(expr_str, zero_map)
    const = -expr_zero

    coeffs = []
    for var in var_names:
        perturb_map = {v: 1.0 if v == var else 0.0 for v in var_names}
        expr_one = safe_eval(expr_str, perturb_map)
        coeff = expr_one - expr_zero
        coeffs.append(coeff)

    return coeffs, const


# ==================== 构建线性系统（修正：仅处理非定值变量） ====================
def build_linear_system(eqs, fixed_vars):
    var_names = extract_variables(eqs, fixed_vars)
    n = len(var_names)
    # 过滤空方程 + 确保方程数=变量数
    valid_eqs = [eq.strip() for eq in eqs if eq.strip()]
    if len(valid_eqs) != n:
        raise ValueError(
            f"非定值方程数 ({len(valid_eqs)}) 与非定值变量数 ({n}) 不相等！\n"
            f"非定值变量：{var_names}\n非定值方程：{valid_eqs}"
        )

    W = []
    b = []
    for eq in valid_eqs:
        coeffs, const = parse_equation(eq, var_names)
        W.append(coeffs)
        b.append(const)
    return W, b, var_names


# ==================== 高斯消元求解器（无修改） ====================
def solve_linear(W, b):
    n = len(W)
    A = copy.deepcopy(W)
    B = copy.deepcopy(b)

    for i in range(n):
        # 部分选主元
        max_row = i
        max_val = abs(A[i][i])
        for k in range(i + 1, n):
            if abs(A[k][i]) > max_val:
                max_val = abs(A[k][i])
                max_row = k
        if max_val < 1e-12:
            raise ValueError("矩阵奇异（行列式为0），无法求解")
        # 交换行
        if max_row != i:
            A[i], A[max_row] = A[max_row], A[i]
            B[i], B[max_row] = B[max_row], B[i]

        # 消元
        for j in range(i + 1, n):
            factor = A[j][i] / A[i][i]
            for k in range(i, n):
                A[j][k] -= factor * A[i][k]
            B[j] -= factor * B[i]

    # 回代
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = sum(A[i][j] * x[j] for j in range(i + 1, n))
        x[i] = (B[i] - s) / A[i][i]
    return x


# ==================== 主求解函数（核心重构） ====================
def solve(eqs: list):
    try:
        # 步骤1：提取定值变量（支持链式替换）
        fixed_vars, remaining_eqs = extract_fixed_vars(eqs)

        # 步骤2：代入定值变量，化简方程（仅替换变量名）
        simplified_eqs = substitute_fixed_vars(remaining_eqs, fixed_vars)

        # 步骤3：构建线性系统（仅处理非定值变量）
        W, b, var_names = build_linear_system(simplified_eqs, fixed_vars)
        print("系数矩阵 W:")
        for row in W:
            print([round(val, 6) for val in row])
        print("右端向量 b:\n", [round(val, 6) for val in b])

        # 步骤4：求解
        x = solve_linear(W, b)
        print("\n解向量:\n")
        # 合并定值变量 + 求解变量
        all_solutions = fixed_vars.copy()
        for var, val in zip(var_names, x):
            all_solutions[var] = val
        # 排序输出
        for var in sorted(all_solutions.keys()):
            print(f"{var} = {all_solutions[var]:.8f}")

    except ValueError as e:
        print(f"\n 错误: {e}")


if __name__ == "__main__":
    # 测试案例（原电路方程组）
    eqs = [
        "R_3 = 200",
        "R_4 = 510",
        "R_5 = 300",
        "i_s = i_1",
        "i_1 = 20*10**(-3)",
        "v_1 / R_3 + (v_2 - v_3) / R_5 = i_1",
        "v_2 - v_1 = 3",
        "(v_2 - v_3) / R_5 = v_3 / R_4",
        "i_2 = v_1 / R_3",
        "i_3 = (v_2 - v_3)/R_5",
        "UAB = v_3",
        "UBC = -v_1",
        "UCA = v_1 - v_3",
    ]
    solve(eqs)