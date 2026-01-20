#!/usr/bin/env python3
"""
Quick summary of generated Mathematica test files.
"""

import os
from pathlib import Path


def count_tests_in_file(filepath: Path) -> int:
    """Count the number of test cases in a Mathematica test file."""
    if not filepath.exists():
        return 0

    with open(filepath, "r") as f:
        content = f.read()

    # Count FormatResult calls
    return content.count("FormatResult[")


def main():
    """Print summary of all test files."""

    test_dir = Path(__file__).parent

    test_files = [
        ("test_virasoro_mathematica.wls", "Virasoro", "T"),
        ("test_sl2_mathematica.wls", "sl(2) Kac-Moody", "J+, J0, J-"),
        ("test_u1_mathematica.wls", "U(1) Current", "T, J"),
        ("test_bcbetagamma_mathematica.wls", "bc-βγ Systems", "b, c, β, γ"),
        ("test_w3_mathematica.wls", "W₃ Algebra", "T, W"),
        ("test_n4_mathematica.wls", "N=4 Small SCFA", "J±, J0, G±, Gt±, T"),
    ]

    print("=" * 80)
    print("MATHEMATICA NESTED NO EXPANSION TEST FILES SUMMARY")
    print("=" * 80)
    print()

    total_tests = 0

    for filename, algebra_name, generators in test_files:
        filepath = test_dir / filename
        num_tests = count_tests_in_file(filepath)
        total_tests += num_tests

        status = "✓" if filepath.exists() else "✗"

        print(f"{status} {filename}")
        print(f"   Algebra: {algebra_name}")
        print(f"   Generators: {generators}")
        print(f"   Test cases: {num_tests}")
        print()

    print("=" * 80)
    print(f"TOTAL TEST CASES: {total_tests}")
    print("=" * 80)
    print()

    # Print usage instructions
    print("USAGE:")
    print("------")
    print("1. Run all tests:")
    print("   ./run_all_mathematica_tests.sh")
    print()
    print("2. Run individual test:")
    print("   wolframscript -file test_virasoro_mathematica.wls")
    print()
    print("3. Parse results and generate Python tests:")
    print("   python3 parse_mathematica_results.py")
    print()
    print("=" * 80)

    # Check if OPEdefs.m exists
    opedefs_path = Path("/Users/lelouch/pyope/OPEdefs/OPEdefs.m")
    if opedefs_path.exists():
        print("✓ OPEdefs.m found at:", opedefs_path)
    else:
        print("✗ WARNING: OPEdefs.m not found at:", opedefs_path)
        print("  Please update the path in all .wls files")
    print()


if __name__ == "__main__":
    main()
