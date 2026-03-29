import re
import copy


# ==================== 安全求值函数 ====================
def safe_eval(expr_str, var_map):
    """
    在安全环境中计算表达式，变量值由 var_map 提供。
    支持运算符：+ - * / ** ( )，以及数字、变量名（形如 x_a，a为字母/数字）。
    """
    globals_dict = {'__builtins__': None}
    try:
        result = eval(expr_str, globals_dict, var_map)
        if isinstance(result, (int, float)):
            return float(result)
        else:
            raise ValueError(f"表达式求值结果不是数字: {result}")
    except Exception as e:
        raise ValueError(f"表达式求值失败: {expr_str}\n错误: {e}")


# ==================== 从方程组中提取所有变量（核心修改） ====================
def extract_variables(eqs):
    """
    从所有方程字符串中提取所有形如 x_字母/数字 的变量（如x_a、x_1、x_b5）
    返回按自然顺序排序的变量列表
    """
    var_set = set()
    # 正则规则详解：
    # \b           ：单词边界（避免匹配运算符/数字粘连的部分，如1x、x+中的x）
    # [a-zA-Z]+    ：开头必须是1+个字母（核心规则）
    # (?:_[a-zA-Z0-9]*)? ：可选组，支持：
    #   - 空（纯字母）、_（字母+下划线）、_字母/数字（字母_字母/数字）
    # \b           ：单词边界收尾
    var_pattern = r'\b[a-zA-Z]+(?:_[a-zA-Z0-9]*)?\b'
    for eq in eqs:
        eq_stripped = eq.strip()
        if not eq_stripped:
            raise ValueError("检测到空方程，请填写有效的线性方程")
        # 提取当前方程中的所有合法变量
        vars_in_eq = re.findall(var_pattern, eq_stripped)
        if not vars_in_eq:
            raise ValueError(f"方程 {eq} 中未检测到合法变量（支持：x_1、x_a、x_b5 等）")
        var_set.update(vars_in_eq)

    if not var_set:
        raise ValueError("未检测到任何变量，请确保变量形如 x_1、x_a、x_b5 等")

    # 自然排序（兼顾字母和数字，如x_1 < x_a < x_b < x_b5）
    return sorted(var_set)


# ==================== 解析单个方程（无需修改） ====================
def parse_equation(eq_str, var_names):
    """
    解析一个线性方程，返回系数列表和常数项。
    方程形式：左边 = 右边
    标准化为： Σ coeff_i * x_i = const
    """
    if '=' not in eq_str:
        raise ValueError(f"方程缺少等号: {eq_str}")
    left, right = eq_str.split('=', 1)
    left = left.strip()
    right = right.strip()
    expr_str = f"({left}) - ({right})"

    # 变量值全为0的映射
    zero_map = {var: 0.0 for var in var_names}
    expr_zero = safe_eval(expr_str, zero_map)
    const = -expr_zero

    # 计算每个变量的系数
    coeffs = []
    for var in var_names:
        perturb_map = {v: 1.0 if v == var else 0.0 for v in var_names}
        expr_one = safe_eval(expr_str, perturb_map)
        coeff = expr_one - expr_zero
        coeffs.append(coeff)

    return coeffs, const


# ==================== 构建整个线性系统====================
def build_linear_system(eqs):
    """
    从方程组字符串列表构建 W 矩阵和 b 向量。
    返回 (W, b)，其中 W 为方阵，b 为列向量（列表）。
    """
    var_names = extract_variables(eqs)
    n = len(var_names)
    if len(eqs) != n:
        raise ValueError(f"方程组数量 ({len(eqs)}) 与变量数量 ({n}) 不相等，无法求解方阵系统。")

    W = []
    b = []
    for eq in eqs:
        coeffs, const = parse_equation(eq, var_names)
        W.append(coeffs)
        b.append(const)
    return W, b, var_names  # 新增返回变量名列表


# ==================== 高斯消元求解器====================
def solve_linear(W, b):
    n = len(W)
    A = copy.deepcopy(W)
    B = copy.deepcopy(b)

    for i in range(n):
        # 部分选主元（避免除零/精度问题）
        max_row = i
        max_val = abs(A[i][i])
        for k in range(i + 1, n):
            if abs(A[k][i]) > max_val:
                max_val = abs(A[k][i])
                max_row = k
        if max_val < 1e-12:
            raise ValueError("矩阵奇异（行列式为0），无法求解")
        # 交换主元行
        if max_row != i:
            A[i], A[max_row] = A[max_row], A[i]
            B[i], B[max_row] = B[max_row], B[i]

        # 消去下方行的当前列
        for j in range(i + 1, n):
            factor = A[j][i] / A[i][i]
            for k in range(i, n):
                A[j][k] -= factor * A[i][k]
            B[j] -= factor * B[i]

    # 回代求解
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        s = sum(A[i][j] * x[j] for j in range(i + 1, n))
        x[i] = (B[i] - s) / A[i][i]
    return x

def solve(eqs:list):
    try:
        W, b, var_names = build_linear_system(eqs)
        print("系数矩阵 W:")
        for row in W:
            print([round(val, 4) for val in row])  # 四舍五入方便查看
        print("右端向量 b:", [round(val, 4) for val in b])

        x = solve_linear(W, b)
        print("\n解向量:")
        for var_name, val in zip(var_names, x):
            print(f"{var_name} = {val:.4f}")
    except ValueError as e:
        print(f"错误: {e}")

if __name__ == "__main__":
    #案例
    eqs = [
        "-x_1 + 3**2*(x_1 - x_2) = x_2 + 1",
        "x_2 + x_1 =  5"
    ]
    solve(eqs)