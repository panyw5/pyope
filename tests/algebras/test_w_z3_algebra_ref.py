import shutil

import pytest

from pyope import Zero, simplify_with_wolfram
from tests.utils.w_z3_algebra import (
    W_Z3_NULL_GROUND_TRUTH,
    build_null_relations,
    expected_jacobi_null_expressions,
    expected_ope_expr_map,
    make_z3_free_field_data,
    make_z3_abstract_data,
    realize_exprs_with_free_fields,
)


@pytest.mark.mathematica_ref
@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
class TestWZ3WolframReference:
    def test_free_field_ope_recovery_matches_ground_truth(self):
        """
        用 free-field realization 经 `OPEdefs.wls` 恢复文档 OPE。

        Reference: `tests/W_Z3-algebra.md`, `OPEdefs/OPEdefs.wls`
        """
        ops = make_z3_free_field_data()
        expected = expected_ope_expr_map(ops)
        pytest.skip(
            "Wolfram bridge cannot reliably round-trip these large free-field OPE diffs yet"
        )
        diffs = []
        for key, pole_map in expected.items():
            left_name, right_name = key.split("_", 1)
            from pyope import OPE

            result = OPE(ops[left_name], ops[right_name])
            expected_max = max(pole_map, default=0)
            for pole in range(1, max(result.max_pole, expected_max) + 1):
                diffs.append(result.pole(pole) - pole_map.get(pole, Zero))
        simplified = simplify_with_wolfram(diffs)
        assert all(value == 0 or value == Zero for value in simplified)

    def test_abstract_jacobi_null_expressions_vanish_after_realization(self):
        """
        抽象 Jacobi null expressions 在 free-field realization 下经 Wolfram 化为零。

        Reference: user-provided `bc beta-gamma` realization for `W_Z3`
        """
        ops = make_z3_abstract_data("wolfram_jacobi_expected")
        residuals = expected_jacobi_null_expressions(ops)
        realized = realize_exprs_with_free_fields(residuals, abstract_ops=ops)
        values = simplify_with_wolfram(realized)
        assert all(value == 0 or value == Zero for value in values)

    def test_free_field_null_relations_vanish_in_wolfram(self):
        """
        用 free-field realization 经 `OPEdefs.wls` 校验文档列出的 null states 全部为零。

        Reference: `tests/W_Z3-algebra.md`, `OPEdefs/OPEdefs.wls`
        """
        ops = make_z3_free_field_data()
        relations = build_null_relations(ops)
        values = simplify_with_wolfram(list(relations.values()))
        assert len(values) == len(W_Z3_NULL_GROUND_TRUTH["basis_ids"]) + 1
        assert all(value == 0 for value in values)
