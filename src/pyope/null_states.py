"""
Null States 计算模块

本模块实现 W-algebra null states 的计算功能，包括：
1. 算符展开到自由场表示
2. 系数提取
3. 自由场 Fock 空间基枚举
4. 矩阵构建和秩计算
"""

from typing import Dict, List, Tuple, Set, Any
from fractions import Fraction
from collections import defaultdict
import itertools

from .operators import BasisOperator, d
from .api import NO
from .simplify import simplify
from .constants import One
from .local_operator import OperatorSum, OperatorProduct


def integer_partitions(n: Fraction) -> List[List[Fraction]]:
    """
    生成整数（或半整数）n 的所有分拆

    例如：
    - integer_partitions(3) = [[3], [2,1], [1,1,1]]
    - integer_partitions(Fraction(5,2)) = [[5/2], [2,1/2], [3/2,1], [3/2,1/2,1/2], ...]

    Args:
        n: 要分拆的数（可以是整数或半整数）

    Returns:
        所有可能的分拆列表，每个分拆是一个降序排列的列表
    """
    # 转换为 Fraction 以支持半整数
    n = Fraction(n)

    # 确定最小单位（1 或 1/2）
    if n.denominator == 1:
        min_unit = Fraction(1)
    else:
        min_unit = Fraction(1, 2)

    def partition_helper(target: Fraction, max_value: Fraction) -> List[List[Fraction]]:
        """递归生成分拆"""
        if target == 0:
            return [[]]
        if target < 0 or max_value <= 0:
            return []

        result = []

        # 尝试所有可能的第一个元素（从 max_value 开始递减）
        current = max_value
        while current >= min_unit:
            if current <= target:
                # 使用 current 作为第一个元素
                for rest in partition_helper(target - current, current):
                    result.append([current] + rest)
            current -= min_unit

        return result

    return partition_helper(n, n)


class CoefficientExtractor:
    """
    从算符表达式中提取系数

    用于将算符的线性组合分解为 {算符: 系数} 的形式
    支持 SymPy 表达式系统
    """

    def __init__(self):
        self.coefficients = {}

    def extract(self, expr) -> Dict[Any, Any]:
        """
        提取表达式中的系数

        Args:
            expr: 算符表达式（SymPy 表达式）

        Returns:
            字典 {算符: 系数}
        """
        import sympy as sp

        self.coefficients = {}

        # 如果是加法表达式，分别处理每一项
        if isinstance(expr, sp.Add):
            for term in expr.args:
                self._extract_term(term)
        else:
            # 单个项
            self._extract_term(expr)

        return self.coefficients

    def _extract_term(self, term):
        """提取单个项的系数和算符"""
        import sympy as sp

        # 如果是乘法表达式，分离系数和算符
        if isinstance(term, sp.Mul):
            coeff = 1
            operator_parts = []

            for factor in term.args:
                if self._is_scalar(factor):
                    coeff *= factor
                else:
                    operator_parts.append(factor)

            # 构造算符
            if len(operator_parts) == 0:
                operator = One
            elif len(operator_parts) == 1:
                operator = operator_parts[0]
            else:
                # 多个算符相乘，重新构造乘法表达式
                operator = sp.Mul(*operator_parts)

            # 累加系数
            if operator in self.coefficients:
                self.coefficients[operator] += coeff
            else:
                self.coefficients[operator] = coeff

        else:
            # 单个算符（系数为 1）
            if term in self.coefficients:
                self.coefficients[term] += 1
            else:
                self.coefficients[term] = 1

    def _is_scalar(self, obj) -> bool:
        """判断对象是否为标量（数字）"""
        import sympy as sp
        return isinstance(obj, (int, float, Fraction, complex, sp.Number))


class FockSpaceBasis:
    """
    自由场 Fock 空间基枚举器

    枚举给定 level 的所有自由场算符组合
    """

    def __init__(self, free_fields: List[BasisOperator]):
        """
        Args:
            free_fields: 自由场列表（如 [b, c, beta, gamma]）
        """
        self.free_fields = free_fields
        self.field_weights = {
            field: field.conformal_weight
            for field in free_fields
        }

    def enumerate_basis(self, level: Fraction, max_count: int = 1000) -> List[Any]:
        """
        枚举给定 level 的所有自由场基

        Args:
            level: 目标 level（共形权重）
            max_count: 最大枚举数量

        Returns:
            自由场算符组合列表
        """
        basis = []

        # 使用整数分拆方法
        for partition in self._integer_partitions(level):
            # 对每个分拆，尝试分配自由场
            for assignment in self._assign_fields_to_partition(partition):
                basis.append(self._construct_operator(assignment))

                if len(basis) >= max_count:
                    return basis

        return basis

    def _integer_partitions(self, n: Fraction):
        """
        生成 n 的所有整数分拆（支持半整数）

        为了处理半整数权重，将 n 乘以 2 转换为整数
        """
        n_doubled = int(2 * n)

        def partitions_helper(n, max_val=None):
            if max_val is None:
                max_val = n
            if n == 0:
                yield []
                return
            for i in range(min(n, max_val), 0, -1):
                for p in partitions_helper(n - i, i):
                    yield [i] + p

        # 转换回半整数
        for p in partitions_helper(n_doubled):
            yield [Fraction(x, 2) for x in p]

    def _assign_fields_to_partition(self, partition: List[Fraction]):
        """
        将自由场分配到分拆的各个部分

        对于分拆 [w1, w2, ...], 找到所有可能的自由场组合
        使得每个部分的权重匹配
        """
        # 对每个权重，找到可以达到该权重的自由场（包括导数）
        available_fields = []

        for weight in partition:
            fields_at_weight = []

            for field in self.free_fields:
                base_weight = self.field_weights[field]
                derivative_order = weight - base_weight

                # 检查是否可以通过导数达到目标权重
                if derivative_order >= 0 and derivative_order == int(derivative_order):
                    fields_at_weight.append((field, int(derivative_order)))

            available_fields.append(fields_at_weight)

        # 生成所有组合
        if all(available_fields):
            for combo in itertools.product(*available_fields):
                yield combo

    def _construct_operator(self, assignment: Tuple) -> Any:
        """
        根据分配构造算符

        Args:
            assignment: [(field1, deriv1), (field2, deriv2), ...]

        Returns:
            构造的算符（可能是 NO 的嵌套）
        """
        if len(assignment) == 0:
            return One

        if len(assignment) == 1:
            field, deriv_order = assignment[0]
            if deriv_order == 0:
                return field
            else:
                return d(field, deriv_order)

        # 多个算符，构造正规序乘积
        operators = []
        for field, deriv_order in assignment:
            if deriv_order == 0:
                operators.append(field)
            else:
                operators.append(d(field, deriv_order))

        # 从右到左构造 NO
        result = operators[-1]
        for op in reversed(operators[:-1]):
            result = NO(op, result)

        return result


def extract_coefficients(expr) -> Dict[Any, Any]:
    """
    便捷函数：提取表达式中的系数

    Args:
        expr: 算符表达式

    Returns:
        字典 {算符: 系数}
    """
    extractor = CoefficientExtractor()
    return extractor.extract(expr)


def enumerate_fock_basis(free_fields: List[BasisOperator],
                         level: Fraction,
                         max_count: int = 1000) -> List[Any]:
    """
    便捷函数：枚举自由场 Fock 空间基

    Args:
        free_fields: 自由场列表
        level: 目标 level
        max_count: 最大枚举数量

    Returns:
        自由场算符组合列表
    """
    enumerator = FockSpaceBasis(free_fields)
    return enumerator.enumerate_basis(level, max_count)


class OperatorExpander:
    """
    算符展开器

    将抽象算符（如 W-algebra 生成元）展开为自由场 Fock 空间基的线性组合
    """

    def __init__(self, fock_basis: List[Any]):
        """
        Args:
            fock_basis: 自由场 Fock 空间基列表
        """
        self.fock_basis = fock_basis
        self.basis_index = {str(basis): i for i, basis in enumerate(fock_basis)}

    def expand(self, operator) -> Dict[Any, Any]:
        """
        将算符展开为自由场基的线性组合

        Args:
            operator: 要展开的算符

        Returns:
            字典 {自由场基: 系数}
        """
        # 使用 simplify 简化算符
        from .simplify import simplify
        simplified = simplify(operator)

        # 提取系数
        coeffs = extract_coefficients(simplified)

        # 将系数映射到 Fock 空间基
        expansion = {}
        for op, coeff in coeffs.items():
            # 检查算符是否在 Fock 基中
            op_str = str(op)
            if op_str in self.basis_index:
                expansion[op] = coeff
            else:
                # 如果不在基中，尝试进一步展开
                # 这里可以添加更复杂的展开逻辑
                expansion[op] = coeff

        return expansion


class CoefficientMatrixBuilder:
    """
    系数矩阵构建器

    构建从抽象算符到自由场基的系数矩阵
    """

    def __init__(self, fock_basis: List[Any], operators: List[Any]):
        """
        Args:
            fock_basis: 自由场 Fock 空间基列表
            operators: 抽象算符列表
        """
        self.fock_basis = fock_basis
        self.operators = operators
        self.expander = OperatorExpander(fock_basis)

    def build_matrix(self) -> List[List[Any]]:
        """
        构建系数矩阵 M[i, j] = operator[j] 在 fock_basis[i] 上的系数

        Returns:
            系数矩阵（列表的列表）
        """
        import numpy as np
        from sympy import sympify, simplify as sp_simplify

        n_basis = len(self.fock_basis)
        n_ops = len(self.operators)

        # 初始化矩阵
        matrix = [[0 for _ in range(n_ops)] for _ in range(n_basis)]

        # 为每个算符展开并填充矩阵
        for j, operator in enumerate(self.operators):
            expansion = self.expander.expand(operator)

            # 填充矩阵列
            for basis_op, coeff in expansion.items():
                # 找到基在列表中的索引
                basis_str = str(basis_op)
                if basis_str in self.expander.basis_index:
                    i = self.expander.basis_index[basis_str]
                    matrix[i][j] = coeff

        return matrix

    def compute_rank(self, matrix: List[List[Any]] = None) -> int:
        """
        计算矩阵的秩

        Args:
            matrix: 系数矩阵（如果为 None，则自动构建）

        Returns:
            矩阵的秩
        """
        if matrix is None:
            matrix = self.build_matrix()

        # 处理空矩阵情况
        if len(matrix) == 0 or (len(matrix) > 0 and len(matrix[0]) == 0):
            return 0

        # 转换为 numpy 数组进行秩计算
        import numpy as np
        from fractions import Fraction

        # 将矩阵转换为浮点数（处理 Fraction 和 sympy 表达式）
        numeric_matrix = []
        for row in matrix:
            numeric_row = []
            for val in row:
                if isinstance(val, Fraction):
                    numeric_row.append(float(val))
                elif hasattr(val, 'evalf'):  # sympy 表达式
                    numeric_row.append(float(val.evalf()))
                else:
                    numeric_row.append(float(val))
            numeric_matrix.append(numeric_row)

        # 计算秩
        np_matrix = np.array(numeric_matrix)
        rank = np.linalg.matrix_rank(np_matrix)

        return rank


class NullStatesCalculator:
    """
    Null States 计算器

    完整的 null states 计算流程：
    1. 枚举给定 level 的自由场 Fock 空间基
    2. 枚举给定 level 的抽象算符
    3. 将抽象算符展开为自由场基的线性组合
    4. 构建系数矩阵并计算秩
    5. 计算 null states 数量 = 抽象态数 - 物理态数
    """

    def __init__(self, free_fields: List[BasisOperator]):
        """
        Args:
            free_fields: 自由场列表（如 [b, c, beta, gamma]）
        """
        self.free_fields = free_fields

    def calculate_null_states(
        self,
        level: Fraction,
        abstract_operators: List[Any],
        max_fock_basis: int = 1000
    ) -> Dict[str, Any]:
        """
        计算给定 level 的 null states

        Args:
            level: 目标 level（共形权重）
            abstract_operators: 抽象算符列表（W-algebra 生成元的组合）
            max_fock_basis: 最大 Fock 基数量

        Returns:
            字典包含：
            - 'level': level
            - 'n_abstract': 抽象态数量
            - 'n_fock_basis': Fock 基数量
            - 'rank': 矩阵秩（物理态数量）
            - 'n_null_states': null states 数量
            - 'fock_basis': Fock 空间基列表
            - 'matrix': 系数矩阵
        """
        # 1. 枚举 Fock 空间基
        fock_basis = enumerate_fock_basis(
            self.free_fields,
            level,
            max_count=max_fock_basis
        )

        # 2. 构建系数矩阵
        matrix_builder = CoefficientMatrixBuilder(fock_basis, abstract_operators)
        matrix = matrix_builder.build_matrix()

        # 3. 计算秩
        rank = matrix_builder.compute_rank(matrix)

        # 4. 计算 null states 数量
        n_abstract = len(abstract_operators)
        n_null_states = n_abstract - rank

        return {
            'level': level,
            'n_abstract': n_abstract,
            'n_fock_basis': len(fock_basis),
            'rank': rank,
            'n_null_states': n_null_states,
            'fock_basis': fock_basis,
            'matrix': matrix
        }

    def calculate_character(
        self,
        max_level: Fraction,
        generator_func,
        step: Fraction = Fraction(1, 2)
    ) -> Dict[Fraction, Dict[str, Any]]:
        """
        计算真空特征（vacuum character）到给定 level

        Args:
            max_level: 最大 level
            generator_func: 函数 level -> 抽象算符列表
            step: level 步长（默认 1/2）

        Returns:
            字典 {level: null_states_result}
        """
        results = {}

        current_level = step
        while current_level <= max_level:
            # 生成该 level 的抽象算符
            abstract_ops = generator_func(current_level)

            # 计算 null states
            result = self.calculate_null_states(current_level, abstract_ops)
            results[current_level] = result

            current_level += step

        return results


def calculate_null_states(
    free_fields: List[BasisOperator],
    level: Fraction,
    abstract_operators: List[Any],
    max_fock_basis: int = 1000
) -> Dict[str, Any]:
    """
    便捷函数：计算 null states

    Args:
        free_fields: 自由场列表
        level: 目标 level
        abstract_operators: 抽象算符列表
        max_fock_basis: 最大 Fock 基数量

    Returns:
        null states 计算结果
    """
    calculator = NullStatesCalculator(free_fields)
    return calculator.calculate_null_states(level, abstract_operators, max_fock_basis)


class QuantumNumberCalculator:
    """
    量子数计算器

    计算算符的量子数 (m, r)，用于按量子数分组
    """

    def __init__(self, quantum_number_map: Dict[Any, Tuple[Fraction, Fraction]]):
        """
        Args:
            quantum_number_map: 字典 {算符: (m, r)}
        """
        self.quantum_number_map = quantum_number_map

    def get_quantum_numbers(self, operator) -> Tuple[Fraction, Fraction]:
        """
        获取算符的量子数

        Args:
            operator: 算符

        Returns:
            (m, r) 量子数元组
        """
        from .operators import DerivativeOperator, NormalOrderedOperator
        from .local_operator import OperatorSum

        # 处理导数算符
        if isinstance(operator, DerivativeOperator):
            # 导数不改变量子数
            return self.get_quantum_numbers(operator.base)

        # 处理正规序算符
        if isinstance(operator, NormalOrderedOperator):
            # 正规序：量子数相加
            m1, r1 = self.get_quantum_numbers(operator.left)
            m2, r2 = self.get_quantum_numbers(operator.right)
            return (m1 + m2, r1 + r2)

        # 处理算符和
        if isinstance(operator, OperatorSum):
            # 算符和：所有项应该有相同的量子数
            # 取第一项的量子数
            if operator.terms:
                first_op = list(operator.terms.keys())[0]
                return self.get_quantum_numbers(first_op)
            return (Fraction(0), Fraction(0))

        # 处理 SymPy Add 表达式（算符的线性组合）
        import sympy as sp
        if isinstance(operator, sp.Add):
            # SymPy Add：所有项应该有相同的量子数
            # 取第一项的量子数
            if operator.args:
                # 获取第一项（可能带系数）
                first_term = operator.args[0]
                # 如果第一项是 Mul（系数 * 算符），提取算符部分
                if isinstance(first_term, sp.Mul):
                    # 找到非数值的部分
                    for arg in first_term.args:
                        if not arg.is_number:
                            return self.get_quantum_numbers(arg)
                # 否则直接使用第一项
                return self.get_quantum_numbers(first_term)
            return (Fraction(0), Fraction(0))

        # 处理基本算符
        if operator in self.quantum_number_map:
            return self.quantum_number_map[operator]

        # 尝试字符串匹配（用于处理复合表达式）
        op_str = str(operator)
        for key, value in self.quantum_number_map.items():
            if str(key) == op_str:
                return value

        # 默认返回 (0, 0)
        return (Fraction(0), Fraction(0))


class QuantumNumberGrouper:
    """
    量子数分组器

    按量子数 (m, r) 对算符进行分组，用于块对角矩阵计算
    """

    def __init__(self, quantum_calculator: QuantumNumberCalculator):
        """
        Args:
            quantum_calculator: 量子数计算器
        """
        self.quantum_calculator = quantum_calculator

    def group_operators(
        self,
        operators: List[Any],
        only_non_negative_m: bool = False
    ) -> Dict[Tuple[Fraction, Fraction], List[Any]]:
        """
        按量子数对算符分组

        Args:
            operators: 算符列表
            only_non_negative_m: 是否只保留 m≥0 的算符（利用对称性）

        Returns:
            字典 {(m, r): [算符列表]}
        """
        groups = defaultdict(list)

        for op in operators:
            m, r = self.quantum_calculator.get_quantum_numbers(op)

            # 如果只保留 m≥0 的算符
            if only_non_negative_m and m < 0:
                continue

            groups[(m, r)].append(op)

        return dict(groups)

    def group_fock_basis(
        self,
        fock_basis: List[Any],
        only_non_negative_m: bool = False
    ) -> Dict[Tuple[Fraction, Fraction], List[Any]]:
        """
        按量子数对 Fock 基分组

        Args:
            fock_basis: Fock 空间基列表
            only_non_negative_m: 是否只保留 m≥0 的基

        Returns:
            字典 {(m, r): [Fock 基列表]}
        """
        return self.group_operators(fock_basis, only_non_negative_m)


class GroupedNullStatesCalculator:
    """
    分组 Null States 计算器

    按量子数分组计算 null states，利用矩阵块对角化提高效率
    """

    def __init__(
        self,
        free_fields: List[BasisOperator],
        quantum_number_map: Dict[Any, Tuple[Fraction, Fraction]]
    ):
        """
        Args:
            free_fields: 自由场列表
            quantum_number_map: 量子数映射 {算符: (m, r)}
        """
        self.free_fields = free_fields
        self.quantum_calculator = QuantumNumberCalculator(quantum_number_map)
        self.grouper = QuantumNumberGrouper(self.quantum_calculator)

    def _filter_linearly_independent(
        self,
        operators: List[Any],
        fock_basis: List[Any]
    ) -> List[Any]:
        """
        过滤出线性独立的算符

        通过逐个添加算符并检查矩阵秩是否增加来判断线性独立性

        Args:
            operators: 算符列表
            fock_basis: Fock 空间基列表

        Returns:
            线性独立的算符列表
        """
        if not operators or not fock_basis:
            return []

        independent_ops = []
        current_rank = 0

        for op in operators:
            # 尝试添加这个算符
            test_ops = independent_ops + [op]

            # 构建矩阵并计算秩
            matrix_builder = CoefficientMatrixBuilder(fock_basis, test_ops)
            matrix = matrix_builder.build_matrix()
            new_rank = matrix_builder.compute_rank(matrix)

            # 如果秩增加，说明这个算符是线性独立的
            if new_rank > current_rank:
                independent_ops.append(op)
                current_rank = new_rank

        return independent_ops

    def calculate_null_states_grouped(
        self,
        level: Fraction,
        abstract_operators: List[Any],
        max_fock_basis: int = 1000,
        only_non_negative_m: bool = False,
        filter_linearly_independent: bool = True
    ) -> Dict[str, Any]:
        """
        按量子数分组计算 null states

        Args:
            level: 目标 level
            abstract_operators: 抽象算符列表
            max_fock_basis: 最大 Fock 基数量
            only_non_negative_m: 是否只计算 m≥0 扇区
            filter_linearly_independent: 是否过滤线性独立的算符

        Returns:
            字典包含：
            - 'level': level
            - 'groups': 字典 {(m, r): group_result}
            - 'total_n_abstract': 总抽象态数
            - 'total_rank': 总秩
            - 'total_n_null_states': 总 null states 数
            - 'filtered': 是否进行了线性独立性过滤
        """
        # 1. 枚举 Fock 空间基
        fock_basis = enumerate_fock_basis(
            self.free_fields,
            level,
            max_count=max_fock_basis
        )

        # 2. 按量子数分组
        operator_groups = self.grouper.group_operators(
            abstract_operators,
            only_non_negative_m
        )
        fock_groups = self.grouper.group_fock_basis(
            fock_basis,
            only_non_negative_m
        )

        # 3. 对每个量子数扇区计算
        group_results = {}
        total_n_abstract = 0
        total_rank = 0

        for quantum_numbers in operator_groups.keys():
            m, r = quantum_numbers

            # 获取该量子数的算符和 Fock 基
            ops_in_group = operator_groups.get(quantum_numbers, [])
            fock_in_group = fock_groups.get(quantum_numbers, [])

            if not ops_in_group or not fock_in_group:
                # 如果没有算符或 Fock 基，跳过
                continue

            # 过滤线性独立的算符（如果启用）
            if filter_linearly_independent:
                ops_in_group = self._filter_linearly_independent(ops_in_group, fock_in_group)

            # 构建该扇区的矩阵
            matrix_builder = CoefficientMatrixBuilder(fock_in_group, ops_in_group)
            matrix = matrix_builder.build_matrix()
            rank = matrix_builder.compute_rank(matrix)

            # 记录结果
            n_abstract = len(ops_in_group)
            n_null = n_abstract - rank

            group_results[quantum_numbers] = {
                'quantum_numbers': quantum_numbers,
                'n_abstract': n_abstract,
                'n_fock_basis': len(fock_in_group),
                'rank': rank,
                'n_null_states': n_null,
                'operators': ops_in_group,
                'fock_basis': fock_in_group,
                'matrix': matrix
            }

            total_n_abstract += n_abstract
            total_rank += rank

        total_n_null_states = total_n_abstract - total_rank

        return {
            'level': level,
            'groups': group_results,
            'total_n_abstract': total_n_abstract,
            'total_rank': total_rank,
            'total_n_null_states': total_n_null_states,
            'only_non_negative_m': only_non_negative_m,
            'filtered': filter_linearly_independent
        }


class OperatorEnumerator:
    """
    完整的算符枚举器

    基于整数分拆方法，枚举给定 level 的所有可能算符组合
    包括：
    1. 单个生成元的导数
    2. 两个生成元的正规序乘积及其导数
    3. 更复杂的嵌套组合
    """

    def __init__(self, generators: Dict[str, Dict[str, Any]]):
        """
        Args:
            generators: 字典 {name: {'op': operator, 'weight': weight, ...}}
        """
        self.generators = generators
        self.gen_list = list(generators.items())

    def enumerate_operators(
        self,
        level: Fraction,
        max_derivative_order: int = 10,
        use_partition_method: bool = True
    ) -> List[Any]:
        """
        枚举给定 level 的所有算符

        Args:
            level: 目标 level（共形权重）
            max_derivative_order: 最大导数阶数
            use_partition_method: 是否使用整数分拆方法（推荐，支持任意多生成元）

        Returns:
            算符列表
        """
        if use_partition_method:
            # 使用整数分拆方法（类似 Mathematica）
            operators = self._enumerate_partition_based(level, max_derivative_order)
        else:
            # 旧方法：只支持单个和两个生成元
            operators = []
            operators.extend(self._enumerate_single_generators(level, max_derivative_order))
            operators.extend(self._enumerate_two_generator_products(level, max_derivative_order))

        # 去重（使用字符串表示）
        unique_ops = []
        seen = set()
        for op in operators:
            op_str = str(op)
            if op_str not in seen:
                seen.add(op_str)
                unique_ops.append(op)

        return unique_ops

    def _enumerate_partition_based(
        self,
        level: Fraction,
        max_derivative_order: int
    ) -> List[Any]:
        """
        基于整数分拆的算符枚举（类似 Mathematica 的方法）

        这个方法支持任意数量生成元的乘积，例如：
        - 单个生成元：g
        - 两个生成元：NO(g1, g2)
        - 三个生成元：NO(g1, g2, g3)
        - 四个生成元：NO(g1, g2, g3, g4)

        算法：
        1. 生成 level 的所有整数分拆
        2. 对于每个分拆 [w1, w2, ..., wn]：
           - 为每个权重 wi 生成所有可能的算符（包括导数）
           - 对这些算符列表进行笛卡尔积
           - 将每个组合用 NO 包装

        Args:
            level: 目标 level
            max_derivative_order: 最大导数阶数

        Returns:
            算符列表
        """
        operators = []

        # 生成所有整数分拆
        partitions = integer_partitions(level)

        for partition in partitions:
            # 对于每个分拆，生成所有可能的排列
            # 例如：[2, 1] -> [[2, 1], [1, 2]]
            # 这样可以生成 NO(∂J, J) 和 NO(J, ∂J) 等不同顺序的组合
            unique_permutations = self._get_unique_permutations(partition)

            for perm in unique_permutations:
                # 为排列中的每个权重生成算符列表
                operator_lists = []
                for weight in perm:
                    ops_at_weight = self._generate_operators_at_weight(weight, max_derivative_order)
                    if ops_at_weight:
                        operator_lists.append(ops_at_weight)

                # 如果所有权重都有对应的算符，生成笛卡尔积
                if len(operator_lists) == len(perm):
                    # 使用 itertools.product 生成笛卡尔积
                    for combination in itertools.product(*operator_lists):
                        # 将组合包装成 NO
                        if len(combination) == 1:
                            # 单个算符，直接添加
                            operators.append(combination[0])
                        else:
                            # 多个算符，用嵌套的 NO 包装
                            # NO(a, b, c) -> NO(a, NO(b, c))
                            # NO(a, b, c, d) -> NO(a, NO(b, NO(c, d)))
                            operators.append(self._build_nested_no(combination))

        return operators

    def _get_unique_permutations(self, partition: List[Fraction]) -> List[List[Fraction]]:
        """
        生成分拆的所有唯一排列

        例如：
        - [2, 1] -> [[2, 1], [1, 2]]
        - [1, 1, 1] -> [[1, 1, 1]]  (只有一个唯一排列)
        - [2, 1, 1] -> [[2, 1, 1], [1, 2, 1], [1, 1, 2]]

        Args:
            partition: 整数分拆（降序排列）

        Returns:
            所有唯一排列的列表
        """
        from itertools import permutations

        # 生成所有排列并去重
        unique_perms = set()
        for perm in permutations(partition):
            unique_perms.add(perm)

        # 转换回列表形式
        return [list(perm) for perm in unique_perms]

    def _build_nested_no(self, operators: tuple) -> Any:
        """
        将多个算符构建为嵌套的 NO 形式

        例如：
        - (a, b) -> NO(a, b)
        - (a, b, c) -> NO(a, NO(b, c))
        - (a, b, c, d) -> NO(a, NO(b, NO(c, d)))

        Args:
            operators: 算符元组

        Returns:
            嵌套的 NO 表达式
        """
        if len(operators) == 2:
            return NO(operators[0], operators[1])
        else:
            # 递归构建：NO(first, NO(rest...))
            return NO(operators[0], self._build_nested_no(operators[1:]))

    def _generate_operators_at_weight(
        self,
        weight: Fraction,
        max_derivative_order: int
    ) -> List[Any]:
        """
        生成指定权重的所有可能算符（包括导数）

        对于每个生成元 g，如果 weight(g) + n = target_weight，
        则生成 ∂^n g

        Args:
            weight: 目标权重
            max_derivative_order: 最大导数阶数

        Returns:
            该权重下的所有算符列表
        """
        operators = []

        for name, gen_info in self.generators.items():
            base_op = gen_info['op']
            base_weight = gen_info['weight']

            # 计算需要的导数阶数
            deriv_order = weight - base_weight

            # 检查导数阶数是否有效
            if deriv_order >= 0 and deriv_order <= max_derivative_order:
                # 检查导数阶数是否为整数
                if deriv_order.denominator == 1:
                    deriv_order_int = int(deriv_order)
                    if deriv_order_int == 0:
                        operators.append(base_op)
                    else:
                        operators.append(d(base_op, deriv_order_int))

        return operators

    def _enumerate_single_generators(
        self,
        level: Fraction,
        max_derivative_order: int
    ) -> List[Any]:
        """
        枚举单个生成元的所有可能导数

        对于每个生成元 g，枚举 g, ∂g, ∂²g, ..., ∂ⁿg
        使得 weight(g) + n = level
        """
        operators = []

        for name, gen_info in self.generators.items():
            base_op = gen_info['op']
            base_weight = gen_info['weight']

            # 枚举所有可能的导数阶数
            for deriv_order in range(max_derivative_order + 1):
                total_weight = base_weight + deriv_order

                if total_weight == level:
                    if deriv_order == 0:
                        operators.append(base_op)
                    else:
                        operators.append(d(base_op, deriv_order))

        return operators

    def _enumerate_two_generator_products(
        self,
        level: Fraction,
        max_derivative_order: int
    ) -> List[Any]:
        """
        枚举两个生成元的正规序乘积（优化版）

        使用右结合形式：NO(a1, NO(a2, a3))
        并且只枚举 i <= j 的组合以避免重复

        对于生成元 g1, g2，枚举：
        - NO(∂^n1 g1, ∂^n2 g2)，其中 i <= j
        - ∂^m NO(∂^n1 g1, ∂^n2 g2)
        使得总权重 = level
        """
        operators = []

        for i, (name1, gen1) in enumerate(self.gen_list):
            # 只枚举 j >= i 的组合，避免重复
            for j in range(i, len(self.gen_list)):
                name2, gen2 = self.gen_list[j]

                # 枚举第一个生成元的导数阶数
                for d1 in range(max_derivative_order + 1):
                    # 枚举第二个生成元的导数阶数
                    for d2 in range(max_derivative_order + 1):
                        weight1 = gen1['weight'] + d1
                        weight2 = gen2['weight'] + d2
                        product_weight = weight1 + weight2

                        # 情况 1: NO(∂^d1 g1, ∂^d2 g2) 本身
                        if product_weight == level:
                            op1 = d(gen1['op'], d1) if d1 > 0 else gen1['op']
                            op2 = d(gen2['op'], d2) if d2 > 0 else gen2['op']
                            operators.append(NO(op1, op2))

                        # 情况 2: ∂^m NO(∂^d1 g1, ∂^d2 g2)
                        for outer_deriv in range(1, max_derivative_order + 1):
                            if product_weight + outer_deriv == level:
                                op1 = d(gen1['op'], d1) if d1 > 0 else gen1['op']
                                op2 = d(gen2['op'], d2) if d2 > 0 else gen2['op']
                                no_product = NO(op1, op2)
                                operators.append(d(no_product, outer_deriv))

        return operators



