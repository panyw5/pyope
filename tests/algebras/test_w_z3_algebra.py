import pytest

from pyope import OPE, Zero, simplify, simplify_with_wolfram
from tests.utils.w_z3_algebra import (
    W_Z3_JACOBI_GROUND_TRUTH,
    W_Z3_NULL_GROUND_TRUTH,
    build_null_relations,
    canonicalize_up_to_sign,
    collect_python_jacobi_residuals,
    expected_ope_expr_map,
    expected_jacobi_null_expressions_canonical,
    make_z3_abstract_data,
    make_z3_free_field_data,
    realize_exprs_with_free_fields,
)


class TestWZ3AlgebraComputations:
    def test_pyope_recovers_documented_ope_from_free_field_ground_truth(self):
        """
        `pyope` 从 free-field realization 恢复文档 OPE，并逐项比对 ground truth。

        Reference: `tests/W_Z3-algebra.md`
        """
        ops = make_z3_free_field_data()
        expected = expected_ope_expr_map(ops)
        for key, pole_map in expected.items():
            left_name, right_name = key.split("_", 1)
            result = OPE(ops[left_name], ops[right_name])
            expected_max = max(pole_map, default=0)
            for pole in range(1, result.max_pole + 1):
                diff = simplify(result.pole(pole) - pole_map.get(pole, Zero))
                assert diff == 0 or diff == Zero
            for pole in range(result.max_pole + 1, expected_max + 1):
                diff = simplify(Zero - pole_map.get(pole, Zero))
                assert diff == 0 or diff == Zero

    @pytest.mark.slow
    def test_pyope_checks_abstract_jacobi_against_ground_truth(self):
        """
        `pyope` 收集抽象 Jacobi 残差表达式，并与去重后的 null-expression ground truth 对齐。

        Reference: `tests/W_Z3-algebra.md`, `OPEdefs/OPEdefs.m`
        """
        ops = make_z3_abstract_data("jacobi_expected")
        expected = {
            canonicalize_up_to_sign(expr)
            for expr in expected_jacobi_null_expressions_canonical(ops)
        }
        actual = {
            canonicalize_up_to_sign(expr)
            for expr in collect_python_jacobi_residuals(ops)
        }
        assert actual == expected

    @pytest.mark.slow
    @pytest.mark.mathematica_ref
    def test_pyope_realizes_abstract_jacobi_residuals_to_zero(self):
        """
        抽象 Jacobi 残差在 free-field realization 下全部化为零。

        Reference: user-provided `bc beta-gamma` realization for `W_Z3`
        """
        ops = make_z3_abstract_data("jacobi_realize")
        residuals = collect_python_jacobi_residuals(ops)
        realized = realize_exprs_with_free_fields(residuals, abstract_ops=ops)
        simplified = simplify_with_wolfram(realized)
        assert all(value == 0 or value == Zero for value in simplified)

    @pytest.mark.slow
    @pytest.mark.mathematica_ref
    def test_pyope_checks_selected_weight8_null_states_against_ground_truth(self):
        """
        `pyope` 复现文档列出的 weight-8 null states，并逐条比对 ground truth 为零。

        Reference: `tests/W_Z3-algebra.md`
        """
        ops = make_z3_free_field_data()
        relations = build_null_relations(ops)
        assert len(relations) == len(W_Z3_NULL_GROUND_TRUTH["basis_ids"]) + 1
        simplified = simplify_with_wolfram(list(relations.values()))
        for name, value in zip(relations, simplified):
            assert value == 0 or value == Zero, name
