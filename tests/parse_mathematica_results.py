#!/usr/bin/env python3
"""
Parse Mathematica test results and create Python test cases.

This script reads the output from Mathematica nested NO expansion tests
and generates corresponding Python test cases for pyope.
"""

import re
import os
from pathlib import Path
from typing import Dict, List, Tuple


class MathematicaResultParser:
    """Parser for Mathematica OPEdefs output."""

    def __init__(self, result_file: str):
        self.result_file = result_file
        self.tests = []

    def parse(self) -> List[Dict[str, str]]:
        """Parse Mathematica output file and extract test cases."""
        with open(self.result_file, "r") as f:
            content = f.read()

        # Split by test separator
        test_blocks = re.split(r"={40,}", content)

        for block in test_blocks:
            if "Test:" in block and "Output:" in block:
                test_match = re.search(r"Test:\s*(.+?)(?=\n)", block)
                output_match = re.search(
                    r"Output:\s*(.+?)(?=\n\n|\Z)", block, re.DOTALL
                )

                if test_match and output_match:
                    test_name = test_match.group(1).strip()
                    output = output_match.group(1).strip()

                    self.tests.append({"name": test_name, "output": output})

        return self.tests

    def convert_to_python_expr(self, mathematica_expr: str) -> str:
        """Convert Mathematica expression to Python/pyope format."""
        # This is a simplified converter - you may need to enhance it
        expr = mathematica_expr

        # Replace Mathematica function names
        replacements = {
            "Derivative[1]": "D",
            "Derivative[2]": "D2",
            "Derivative[3]": "D3",
            "NO[": "NO(",
            "][": ", ",
            "One": "1",
        }

        for old, new in replacements.items():
            expr = expr.replace(old, new)

        return expr


def generate_python_test_file(
    algebra_name: str, tests: List[Dict[str, str]], output_file: str
):
    """Generate a Python test file from parsed Mathematica results."""

    test_template = '''"""
Test nested NO expansion for {algebra_name} algebra.
Reference results from Mathematica OPEdefs.
"""

import pytest
from pyope import *


class Test{class_name}NestedNO:
    """Test nested normal-ordered products for {algebra_name} algebra."""
    
    def setup_method(self):
        """Set up algebra for each test."""
        # TODO: Initialize algebra and OPEs
        pass
{test_methods}
'''

    method_template = '''
    def test_{method_name}(self):
        """Test: {test_description}"""
        # TODO: Implement test
        # Expected result from Mathematica:
        # {expected_result}
        pass
'''

    # Generate test methods
    test_methods = []
    for i, test in enumerate(tests):
        method_name = f"{algebra_name.lower()}_{i + 1:03d}"
        test_description = test["name"]
        expected_result = test["output"]

        method = method_template.format(
            method_name=method_name,
            test_description=test_description,
            expected_result=expected_result.replace("\n", "\n        # "),
        )
        test_methods.append(method)

    # Generate full test file
    class_name = "".join(word.capitalize() for word in algebra_name.split("_"))
    test_content = test_template.format(
        algebra_name=algebra_name,
        class_name=class_name,
        test_methods="".join(test_methods),
    )

    with open(output_file, "w") as f:
        f.write(test_content)

    print(f"Generated {output_file} with {len(tests)} test cases")


def main():
    """Main function to parse all Mathematica results and generate Python tests."""

    # Define algebra types and their result files
    algebras = {
        "virasoro": "virasoro_results.txt",
        "sl2": "sl2_results.txt",
        "u1": "u1_results.txt",
        "bcbetagamma": "bcbetagamma_results.txt",
        "w3": "w3_results.txt",
        "n4": "n4_results.txt",
    }

    script_dir = Path(__file__).parent
    results_dir = script_dir / "mathematica_results"
    output_dir = script_dir / "generated_tests"

    # Create output directory
    output_dir.mkdir(exist_ok=True)

    print("Parsing Mathematica results and generating Python tests...")
    print("=" * 60)

    for algebra_name, result_file in algebras.items():
        result_path = results_dir / result_file

        if not result_path.exists():
            print(f"⚠ Skipping {algebra_name}: result file not found")
            continue

        print(f"\nProcessing {algebra_name}...")

        # Parse Mathematica results
        parser = MathematicaResultParser(str(result_path))
        tests = parser.parse()

        print(f"  Found {len(tests)} test cases")

        # Generate Python test file
        output_file = output_dir / f"test_{algebra_name}_nested_no_from_mathematica.py"
        generate_python_test_file(algebra_name, tests, str(output_file))

    print("\n" + "=" * 60)
    print("All Python test files generated successfully!")
    print(f"Output directory: {output_dir}")


if __name__ == "__main__":
    main()
