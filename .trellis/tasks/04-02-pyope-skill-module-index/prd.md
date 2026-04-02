# Build pyope skill module index

## Goal

Create a documentation pack that maps the pyope repository into teachable modules, with example-driven subfiles that can later be turned into a pyope skill.

## Requirements

- Produce a top-level module index covering the public package surface and the major internal module groups.
- Split the material into several topic files instead of a single monolithic note.
- Include representative usage snippets for each major module cluster.
- Distinguish stable user-facing APIs from experimental research-oriented APIs.
- Map source modules to the most useful demos, tests, and reference materials.

## Acceptance Criteria

- [ ] A new documentation folder exists under `docs/` for this skill-prep material.
- [ ] The top-level index explains how the package is structured and which files matter most.
- [ ] The subfiles provide concrete examples for core usage paths.
- [ ] The docs identify experimental areas such as realization, `C_2`, and null-state search.
- [ ] The docs are suitable as a retrieval source when building a later pyope skill.

## Technical Notes

- Prefer the package export surface in `src/pyope/__init__.py` as the primary API boundary.
- Use README examples, demo notebooks, and tests as the main evidence for usage patterns.
- Treat `tmp_*.py` and `tmp_*.wls` as research scratch files rather than primary skill material.