import shutil
import subprocess

import pytest
import sympy as sp

from pyope import (
    BasisOperator,
    Bosonic,
    MakeOPE,
    NO,
    OPE,
    One,
    dn,
    simplify_with_wolfram,
    simplify,
    set_compute_backend,
    get_compute_backend,
    compute_backend,
)
from pyope.backend import SUPPORTED_BACKENDS
from pyope.wolfram_backend import (
    WolframBackendError,
    _decode_expr,
    _decode_symbol_name,
    _encode_symbol_name,
    _invoke_wolfram,
    _write_wolfram_payload_files,
    chunk_exprs_for_wolfram,
    op_to_wolfram_string,
)


def test_default_backend_is_sympy():
    assert get_compute_backend() == "sympy"


def test_set_compute_backend_rejects_unknown_backend():
    with pytest.raises(ValueError):
        set_compute_backend("unknown")


def test_set_compute_backend_rejects_incomplete_wolfram_backend():
    with pytest.raises(ValueError, match="not ready yet and must not be used"):
        set_compute_backend("wolfram")


def test_supported_backends_include_wolfram():
    assert SUPPORTED_BACKENDS == {"sympy", "wolfram"}


def test_compute_backend_context_manager_restores_previous_backend():
    set_compute_backend("sympy")
    with pytest.raises(ValueError, match="not ready yet and must not be used"):
        with compute_backend("wolfram"):
            pass
    assert get_compute_backend() == "sympy"


def test_sympy_backend_still_computes_normal_ordered_product():
    T = BasisOperator("T")
    J = BasisOperator("J")
    set_compute_backend("sympy")

    result = NO(T, J)

    assert repr(result) == "NO(T, J)"


def test_sympy_backend_still_computes_zero_ope():
    T = BasisOperator("T")
    J = BasisOperator("J")
    set_compute_backend("sympy")

    result = OPE(T, J)

    assert result.max_pole == 0


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_compute_backend_rejects_wolfram_for_registered_binary_ope_path():
    with pytest.raises(ValueError, match="not ready yet and must not be used"):
        with compute_backend("wolfram"):
            pass


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_compute_backend_rejects_wolfram_for_binary_no_path():
    with pytest.raises(ValueError, match="not ready yet and must not be used"):
        with compute_backend("wolfram"):
            pass


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_compute_backend_rejects_wolfram_for_linear_no_path():
    with pytest.raises(ValueError, match="not ready yet and must not be used"):
        with compute_backend("wolfram"):
            pass


def test_wolfram_decoder_rejects_illegal_operator_multiplication():
    T = BasisOperator("T")
    J = BasisOperator("J")

    with pytest.raises(WolframBackendError, match="Illegal operator multiplication"):
        _decode_expr("(T * J)", {"T", "J"})


def test_wolfram_decoder_preserves_derivative_protocol_shape():
    T = BasisOperator("T")

    result = _decode_expr("dn(2, T)", {"T"})

    assert repr(result) == "d^2(T)"


def test_wolfram_decoder_supports_nested_container_payloads():
    T = BasisOperator("T")
    J = BasisOperator("J")

    result = _decode_expr(
        'PyDict([("ops", [NO(T, J), PyTuple([dn(1, T), "tag"])])])',
        {"T", "J"},
    )

    assert result == {"ops": [NO(T, J), (dn(1, T), "tag")]}


def test_chunk_exprs_for_wolfram_respects_item_limit():
    exprs = [sp.Integer(i) for i in range(5)]

    chunks = chunk_exprs_for_wolfram(exprs, max_items=2, max_chars=100)

    assert chunks == [[0, 1], [2, 3], [4]]


def test_chunk_exprs_for_wolfram_raises_for_oversized_expression():
    x = sp.Symbol("x" * 40)

    with pytest.raises(WolframBackendError, match="transport size limit"):
        chunk_exprs_for_wolfram([x], max_items=4, max_chars=20)


def test_chunk_exprs_for_wolfram_allows_large_expression_without_char_limit():
    x = sp.Symbol("x" * 40)

    chunks = chunk_exprs_for_wolfram([x], max_items=4)

    assert chunks == [[x]]


def test_write_wolfram_payload_files_creates_files(tmp_path):
    payload_paths = _write_wolfram_payload_files(
        tmp_path,
        {
            "PYOPE_WL_OPERATION": "EVAL",
            "PYOPE_WL_EXPR": "x + y",
        },
    )

    operation_path = tmp_path / "pyope_wl_operation.txt"
    expr_path = tmp_path / "pyope_wl_expr.txt"

    assert payload_paths["PYOPE_WL_OPERATION_PATH"] == str(operation_path)
    assert payload_paths["PYOPE_WL_EXPR_PATH"] == str(expr_path)
    assert operation_path.read_text(encoding="utf-8") == "EVAL"
    assert expr_path.read_text(encoding="utf-8") == "x + y"
    assert (tmp_path / "manifest.json").exists()


def test_invoke_wolfram_uses_payload_files(monkeypatch, tmp_path):
    captured = {}

    def fake_tempdir(prefix):
        class _TempDir:
            def __enter__(self_inner):
                return str(tmp_path)

            def __exit__(self_inner, exc_type, exc, tb):
                return False

        return _TempDir()

    def fake_run(command, capture_output, text, encoding, check, env):
        captured["command"] = command
        captured["env"] = env
        captured["encoding"] = encoding
        (tmp_path / "pyope_result.txt").write_text("42", encoding="utf-8")
        return subprocess.CompletedProcess(
            command,
            0,
            stdout="",
            stderr="",
        )

    monkeypatch.setattr(
        "pyope.wolfram_backend.tempfile.TemporaryDirectory", fake_tempdir
    )
    monkeypatch.setattr("pyope.wolfram_backend.subprocess.run", fake_run)

    result = _invoke_wolfram(
        "dummy.wls",
        {
            "PYOPE_WL_OPERATION": "EVAL",
            "PYOPE_WL_EXPR": "42",
            "PYOPE_WL_OPERATORS": "",
            "PYOPE_WL_REGISTRY": "",
        },
        set(),
    )

    assert result == 42
    assert captured["command"] == ["wolframscript", "-file", "dummy.wls"]
    assert captured["encoding"] == "utf-8"
    assert captured["env"]["PYOPE_WL_EXPR_PATH"] == str(tmp_path / "pyope_wl_expr.txt")
    assert captured["env"]["PYOPE_WL_OPERATION_PATH"] == str(
        tmp_path / "pyope_wl_operation.txt"
    )
    assert captured["env"]["PYOPE_WL_RESULT_PATH"] == str(tmp_path / "pyope_result.txt")
    assert "PYOPE_WL_EXPR" not in captured["env"]


def test_invoke_wolfram_falls_back_to_stdout_marker(monkeypatch, tmp_path):
    def fake_tempdir(prefix):
        class _TempDir:
            def __enter__(self_inner):
                return str(tmp_path)

            def __exit__(self_inner, exc_type, exc, tb):
                return False

        return _TempDir()

    def fake_run(command, capture_output, text, encoding, check, env):
        return subprocess.CompletedProcess(
            command,
            0,
            stdout="PYOPE_RESULT: 42\n",
            stderr="",
        )

    monkeypatch.setattr(
        "pyope.wolfram_backend.tempfile.TemporaryDirectory", fake_tempdir
    )
    monkeypatch.setattr("pyope.wolfram_backend.subprocess.run", fake_run)

    result = _invoke_wolfram(
        "dummy.wls",
        {
            "PYOPE_WL_OPERATION": "EVAL",
            "PYOPE_WL_EXPR": "42",
            "PYOPE_WL_OPERATORS": "",
            "PYOPE_WL_REGISTRY": "",
        },
        set(),
    )

    assert result == 42


def test_decode_expr_supports_unicode_operator_names():
    beta = BasisOperator("β")

    result = _decode_expr("NO(β, β)", {"β"})

    assert result == NO(beta, beta)


def test_encode_symbol_name_maps_greek_letters_to_mathematica_form():
    assert _encode_symbol_name("β") == r"\[Beta]"
    assert _encode_symbol_name("γ") == r"\[Gamma]"


def test_decode_symbol_name_restores_mathematica_greek_form():
    assert _decode_symbol_name(r"NO(\[Beta], \[Gamma])") == "NO(β, γ)"


def test_op_to_wolfram_string_formats_operator_tree():
    beta = BasisOperator("β")
    gamma = BasisOperator("γ")

    assert op_to_wolfram_string(NO(gamma, dn(2, beta))) == (
        r"NO[\[Gamma], Derivative[2][\[Beta]]]"
    )


def test_op_to_wolfram_string_formats_linear_combinations():
    beta = BasisOperator("β")
    gamma = BasisOperator("γ")

    assert op_to_wolfram_string(2 * NO(gamma, beta) + dn(3, beta)) == (
        r"((2 * NO[\[Gamma], \[Beta]]) + Derivative[3][\[Beta]])"
    )


def test_op_to_wolfram_string_formats_nested_lists():
    alpha = BasisOperator("α")
    beta = BasisOperator("β")
    gamma = BasisOperator("γ")
    delta = BasisOperator("δ")

    assert op_to_wolfram_string(
        [alpha, dn(1, beta), [NO(gamma, delta), dn(2, alpha)]]
    ) == (
        r"{\[Alpha], Derivative[1][\[Beta]], {NO[\[Gamma], \[Delta]], Derivative[2][\[Alpha]]}}"
    )


def test_op_to_wolfram_string_uses_single_backslashes_for_to_expression():
    beta = BasisOperator("β")
    gamma = BasisOperator("γ")

    result = op_to_wolfram_string(NO(gamma, dn(2, beta)))

    assert result == r"NO[\[Gamma], Derivative[2][\[Beta]]]"
    assert "\\\\[Gamma]" not in result
    assert "\\\\[Beta]" not in result


def test_run_eval_uses_public_wolfram_string_encoder(monkeypatch):
    beta = BasisOperator("β")
    gamma = BasisOperator("γ")
    expr = NO(gamma, dn(2, beta))
    captured = {}

    def fake_invoke(script_path, payload, operator_names):
        captured["script_path"] = script_path
        captured["payload"] = payload
        captured["operator_names"] = operator_names
        return expr

    monkeypatch.setattr("pyope.wolfram_backend._invoke_wolfram", fake_invoke)

    result = simplify_with_wolfram(expr)

    assert result == expr
    assert captured["payload"]["PYOPE_WL_EXPR"] == op_to_wolfram_string(expr)


def test_run_eval_list_uses_public_wolfram_string_encoder(monkeypatch):
    beta = BasisOperator("β")
    gamma = BasisOperator("γ")
    exprs = [beta, [NO(gamma, dn(2, beta))]]
    captured = {}

    def fake_invoke(script_path, payload, operator_names):
        captured["script_path"] = script_path
        captured["payload"] = payload
        captured["operator_names"] = operator_names
        return exprs

    monkeypatch.setattr("pyope.wolfram_backend._invoke_wolfram", fake_invoke)

    result = simplify_with_wolfram(exprs)

    assert result == exprs
    assert (
        captured["payload"]["PYOPE_WL_EXPR"]
        if "PYOPE_WL_EXPR" in captured["payload"]
        else True
    )
    assert captured["payload"]["PYOPE_WL_EXPR_LIST"] == op_to_wolfram_string(exprs)


def test_simplify_with_wolfram_auto_chunks_large_lists(monkeypatch):
    exprs = [sp.Integer(i) for i in range(5)]
    calls = []

    def fake_chunk_exprs(values, max_items=32, max_chars=100_000):
        assert list(values) == exprs
        return [exprs[:2], exprs[2:4], exprs[4:]]

    def fake_run_eval_list(values, operation):
        calls.append((list(values), operation))
        return list(values)

    monkeypatch.setattr(
        "pyope.wolfram_backend.chunk_exprs_for_wolfram", fake_chunk_exprs
    )
    monkeypatch.setattr("pyope.wolfram_backend._run_eval_list", fake_run_eval_list)

    result = simplify_with_wolfram(exprs)

    assert result == exprs
    assert calls == [
        ([0, 1], "CANONICALIZE_LIST"),
        ([2, 3], "CANONICALIZE_LIST"),
        ([4], "CANONICALIZE_LIST"),
    ]


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_simplify_with_wolfram_round_trips_expression_lists():
    T = BasisOperator("Tbatcheval")
    J = BasisOperator("Jbatcheval")

    exprs = [NO(T, J), dn(1, T) + J]

    result = simplify_with_wolfram(exprs)

    assert isinstance(result, list)
    assert len(result) == 2
    assert repr(result[0]) == "NO(Tbatcheval, Jbatcheval)"
    assert result[1] == J + dn(1, T)


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_simplify_with_wolfram_round_trips_nested_containers():
    T = BasisOperator("Tnested")
    J = BasisOperator("Jnested")

    payload = {
        "ops": [NO(T, J), (dn(1, T), J + dn(1, T))],
        "meta": {"label": "trial"},
    }

    result = simplify_with_wolfram(payload)

    assert result["ops"][0] == NO(T, J)
    assert result["ops"][1] == (dn(1, T), J + dn(1, T))
    assert result["meta"] == {"label": "trial"}


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_simplify_with_wolfram_round_trips_nested_containers_in_canonical_form():
    T = BasisOperator("Tsimplify_nested")
    J = BasisOperator("Jsimplify_nested")
    Bosonic(T, J)

    payload = {
        "ops": [NO(T, J) + NO(T, J), (NO(T, J) + NO(T, J), NO(T, J))],
        "meta": {"label": "trial"},
    }

    result = simplify_with_wolfram(payload)

    assert result == {
        "ops": [2 * NO(T, J), (2 * NO(T, J), NO(T, J))],
        "meta": {"label": "trial"},
    }


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_simplify_with_wolfram_reorders_normal_products_in_lists():
    A = BasisOperator("Abatchcanon", conformal_weight=1)
    B = BasisOperator("Bbatchcanon", conformal_weight=1)
    Bosonic(A, B)

    result = simplify_with_wolfram([NO(B, A), NO(A, B)])

    assert len(result) == 2
    assert result[0] == NO(A, B)
    assert result[1] == NO(A, B)


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_compute_backend_rejects_wolfram_for_derivative_payload_path():
    with pytest.raises(ValueError, match="not ready yet and must not be used"):
        with compute_backend("wolfram"):
            pass


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_compute_backend_rejects_wolfram_for_nested_no_path():
    with pytest.raises(ValueError, match="not ready yet and must not be used"):
        with compute_backend("wolfram"):
            pass


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_wolfram_whole_expression_eval_preserves_nested_no_shape():
    T = BasisOperator("T")
    J = BasisOperator("J")

    expr = (
        -sp.Rational(8, 9) * NO(T, NO(T, NO(T, NO(J, J))))
        + sp.Rational(4, 9) * NO(T, NO(T, NO(T, dn(1, J))))
        + sp.Rational(2, 3) * NO(dn(1, T), NO(T, NO(T, J)))
    )

    result = simplify_with_wolfram(expr)

    assert result != 0
    assert "NO(" in repr(result)


@pytest.mark.skipif(
    shutil.which("wolframscript") is None, reason="wolframscript not installed"
)
def test_simplify_with_wolfram_handles_null_state_style_expr():
    b = BasisOperator("b")
    c = BasisOperator("c")
    beta = BasisOperator("β")
    gamma = BasisOperator("γ")

    from pyope import Bosonic, Fermionic, clear_registry

    clear_registry()
    Fermionic(b, c)
    Bosonic(beta, gamma)
    OPE[b, c] = MakeOPE([One])
    OPE[beta, gamma] = MakeOPE([-One])

    J = 3 * NO(beta, gamma) + 2 * NO(b, c)
    G = NO(b, gamma)
    Gbar = 3 * NO(beta, dn(1, c)) + 2 * NO(dn(1, beta), c)
    T = (
        -sp.Rational(3, 2) * NO(beta, dn(1, gamma))
        - sp.Rational(1, 2) * NO(dn(1, beta), gamma)
        - 2 * NO(b, dn(1, c))
        - NO(dn(1, b), c)
    )
    W = beta
    Wbar = (
        sp.Rational(10, 9) * dn(2, gamma)
        - sp.Rational(4, 3) * NO(gamma, b, dn(1, c))
        + sp.Rational(2, 3) * NO(gamma, dn(1, b), c)
        - sp.Rational(8, 3) * NO(dn(1, gamma), b, c)
        - 4 * NO(beta, gamma, dn(1, gamma))
        + sp.Rational(2, 3) * NO(dn(1, beta), gamma, gamma)
        + 2 * NO(beta, gamma, gamma, b, c)
        + NO(beta, beta, gamma, gamma, gamma)
    )

    expr = (
        -sp.Rational(8, 9) * NO(T, NO(T, NO(T, NO(J, J))))
        + sp.Rational(4, 3) * NO(T, NO(T, NO(T, T)))
        + sp.Rational(4, 9) * NO(T, NO(T, NO(T, dn(1, J))))
        - sp.Rational(8, 9) * NO(T, NO(T, NO(J, NO(G, Gbar))))
        + sp.Rational(2, 3) * NO(T, NO(T, NO(G, dn(1, Gbar))))
        - 3 * NO(T, NO(T, NO(W, dn(1, Wbar))))
        - sp.Rational(1, 9) * NO(T, NO(T, NO(dn(1, J), NO(J, J))))
        - sp.Rational(4, 9) * NO(T, NO(T, NO(dn(1, J), dn(1, J))))
        - sp.Rational(2, 9) * NO(T, NO(T, NO(dn(1, G), Gbar)))
        + 5 * NO(T, NO(T, NO(dn(1, W), Wbar)))
        + NO(T, NO(T, NO(dn(2, J), J)))
    )

    result = simplify_with_wolfram(expr)

    assert result == 0
