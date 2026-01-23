# Tasks: Add normal_product() Helper Function

## Implementation Tasks

### Phase 1: Core Implementation
- [x] **Task 1.1**: Implement `normal_product()` function in `src/pyope/api.py`
  - Accept variable number of operators (`*operators`)
  - Handle edge cases: empty list, single operator
  - Build nested NO from right to left
  - Add comprehensive docstring with examples
  - Validation: Function signature and basic logic

- [x] **Task 1.2**: Add type annotations to `normal_product()`
  - Use `*operators: Any` for flexible input
  - Return type `Any` for simplified expressions
  - Ensure compatibility with existing operator types

### Phase 2: API Integration
- [x] **Task 2.1**: Export `normal_product` from `src/pyope/__init__.py`
  - Add to import statement from `api.py`
  - Add to `__all__` list for proper public API
  - Validation: `from pyope import normal_product` works

- [x] **Task 2.2**: Update module-level docstring
  - Mention new function in overview
  - Link to detailed documentation

### Phase 3: Testing
- [x] **Task 3.1**: Create comprehensive test suite in `tests/test_normal_product.py`
  - Test empty parameter list → returns One
  - Test single operator → returns operator itself
  - Test two operators → equivalent to NO(A, B)
  - Test three operators → NO(A, NO(B, C))
  - Test four operators → NO(A, NO(B, NO(C, D)))
  - Test with derivative operators
  - Test with Zero operator
  - Test with One operator
  - Test with scalar coefficients
  - Test with sympy symbols
  - Test with many operators (6+)
  - Validation: All 155 tests pass

- [x] **Task 3.2**: Verify Mathematica equivalence
  - Create Wolfram Script test cases
  - Compare Python output with Mathematica
  - Ensure mathematical correctness
  - Validation: Results match OPEdefs.m

- [x] **Task 3.3**: Run existing test suite
  - Ensure no regressions
  - Verify 226/228 tests still pass
  - Validation: `pytest` output shows expected pass rate

### Phase 4: Documentation
- [x] **Task 4.1**: Create function documentation in `docs/normal_product.md`
  - Explain purpose and motivation
  - Document syntax and conversion rules
  - Provide usage examples for all scenarios
  - List advantages over manual NO nesting
  - Include implementation details
  - Link to related functions
  - Validation: Documentation is clear and complete

- [x] **Task 4.2**: Create demonstration notebook in `demo/normal_product_demo.ipynb`
  - Show basic usage examples
  - Demonstrate edge cases
  - Include practical applications (W-algebra, Sugawara)
  - Add explanatory text in Chinese
  - Make it executable and self-contained
  - Validation: Notebook runs without errors

- [x] **Task 4.3**: Update README (if needed)
  - Add mention of new function
  - Update feature list
  - Validation: README reflects current state

### Phase 5: Code Review and Validation
- [ ] **Task 5.1**: Self-review implementation
  - Check code follows project conventions
  - Ensure type hints are correct
  - Verify docstring format
  - Check for potential bugs
  - Validation: Code review checklist

- [ ] **Task 5.2**: Run static analysis tools
  - Run `mypy` for type checking
  - Run `ruff` for linting
  - Run `black` for formatting
  - Fix any issues found
  - Validation: All tools pass

- [ ] **Task 5.3**: Open pull request
  - Create descriptive PR title
  - Link to this proposal
  - Summarize changes in description
  - Request review from maintainers
  - Validation: PR is created and reviewable

- [ ] **Task 5.4**: Address review feedback
  - Respond to all review comments
  - Make requested changes
  - Push updates to branch
  - Validation: Reviewer approves

### Phase 6: Merge and Release
- [ ] **Task 6.1**: Merge to main branch
  - Ensure CI checks pass
  - Merge using squash or rebase
  - Delete feature branch after merge
  - Validation: Changes appear in main

- [ ] **Task 6.2**: Update version (if needed)
  - Determine if version bump required
  - Update version in `pyproject.toml`
  - Validation: Version reflects change

## Dependencies and Ordering

### Critical Path
1. Task 1.1 → Task 1.2 → Task 2.1 → Task 3.1 → Task 3.3 → Task 5.3 → Task 6.1

### Parallelizable Work
- Tasks 3.2 and 4.1 can run in parallel after Phase 1
- Tasks 4.2 and 4.3 can run in parallel after Phase 2
- Tasks 5.1 and 5.2 can run in parallel after Phase 4

### Blocking Relationships
- Task 3.3 (run existing tests) must complete before Task 5.3 (open PR)
- Task 5.2 (static analysis) must pass before Task 5.3 (open PR)
- All implementation tasks must complete before Task 6.1 (merge)

## Definition of Done

A task is considered **complete** when:
- ✅ Implementation is finished
- ✅ Tests are written and passing
- ✅ Documentation is updated
- ✅ Code follows project conventions
- ✅ No regressions in existing functionality

The entire feature is **done** when:
- ✅ All tasks in Phases 1-4 are complete
- ✅ Code review is approved (Phase 5)
- ✅ Changes are merged to main (Phase 6)
- ✅ Documentation is published (if applicable)

## Estimated Task Complexity

| Task | Complexity | Notes |
|------|------------|-------|
| 1.1  | Low        | Straightforward implementation |
| 1.2  | Low        | Type annotations only |
| 2.1  | Low        | Simple export |
| 2.2  | Low        | Documentation update |
| 3.1  | Medium     | Many test cases to cover |
| 3.2  | Medium     | Cross-language validation |
| 3.3  | Low        | Regression testing |
| 4.1  | Medium     | Comprehensive documentation |
| 4.2  | Medium     | Interactive notebook |
| 4.3  | Low        | Minor update |
| 5.1  | Low        | Self-review |
| 5.2  | Low        | Automated tools |
| 5.3  | Low        | Git operations |
| 5.4  | Variable   | Depends on feedback |
| 6.1  | Low        | Merge operation |
| 6.2  | Low        | Optional version bump |

## Validation Checklist

Before opening PR:
- [ ] All tests pass (`pytest`)
- [ ] Type checking passes (`mypy`)
- [ ] Linting passes (`ruff`)
- [ ] Formatting applied (`black`)
- [ ] Documentation is complete
- [ ] No merge conflicts with main
- [ ] Proposal document is finalized

Before merge:
- [ ] At least one approval from code review
- [ ] All review comments addressed
- [ ] CI checks pass
- [ ] No outstanding issues or concerns
