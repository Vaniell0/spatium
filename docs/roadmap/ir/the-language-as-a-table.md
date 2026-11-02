---
date: 2026-10-30
title: the question language as the reader's own tables, and the reference the page reads
files: [include/spatium/io/field_parse.hpp, examples/site_reference.cpp, site/data/reference.json, tests/test_field_parse.cpp]
---

The functions the reader knows, its operators and its named numbers are now tables in `field_parse.hpp` (`functions`: name, arity, a caption, the domain where it has a value, an example; `operators`; `constants`), and the reader takes its names from them: a function added to the table is one the reader reads, and the refusal of an unknown one lists the table. `examples/site_reference.cpp` writes `site/data/reference.json` from the same tables, so the page's reference cannot name a function the reader does not know or miss one it does; every example in it is answered through `respond`, the page's door, and carries the library's LaTeX and value (to twelve digits: the last ones are the libm's). Its `--check` is a ctest. The test holds that every function in the table reads as a call of its arity, that every example reads to a finite number, and, as the case it must fail, that `tan(u)` is refused with the whole table named.
