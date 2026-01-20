#!/bin/bash
# Run all Mathematica nested NO expansion tests
# This script executes all algebra test files and saves results

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
OUTPUT_DIR="${SCRIPT_DIR}/mathematica_results"

# Create output directory if it doesn't exist
mkdir -p "${OUTPUT_DIR}"

echo "=========================================="
echo "Running Mathematica Nested NO Tests"
echo "=========================================="
echo ""

# Array of test files
tests=(
    "test_virasoro_mathematica.wls:virasoro"
    "test_sl2_mathematica.wls:sl2"
    "test_u1_mathematica.wls:u1"
    "test_bcbetagamma_mathematica.wls:bcbetagamma"
    "test_w3_mathematica.wls:w3"
    "test_n4_mathematica.wls:n4"
)

# Run each test
for test_info in "${tests[@]}"; do
    IFS=':' read -r test_file test_name <<< "$test_info"
    
    echo "Running ${test_name} tests..."
    echo "  Input: ${test_file}"
    echo "  Output: ${OUTPUT_DIR}/${test_name}_results.txt"
    
    if [ -f "${SCRIPT_DIR}/${test_file}" ]; then
        wolframscript -file "${SCRIPT_DIR}/${test_file}" > "${OUTPUT_DIR}/${test_name}_results.txt" 2>&1
        
        if [ $? -eq 0 ]; then
            echo "  ✓ ${test_name} tests completed successfully"
        else
            echo "  ✗ ${test_name} tests failed"
        fi
    else
        echo "  ✗ Test file not found: ${test_file}"
    fi
    echo ""
done

echo "=========================================="
echo "All tests completed!"
echo "Results saved in: ${OUTPUT_DIR}"
echo "=========================================="
