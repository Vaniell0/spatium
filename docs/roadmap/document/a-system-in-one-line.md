---
date: 2026-11-16
title: a system in one line
files: [include/spatium/io/document.hpp, tests/test_answer_evaluate.cpp]
---

`x + y = 3, x - y = 1` was refused at its comma ("expected 'on', found ','"), while the same equations a line each were solved. `read_system` now splits one line at the commas outside any bracket and reads the pieces as the lines of a system, so the calculator's one line takes a system as typed; `log(2, x) = 3` keeps its comma. Systems of two to four unknowns that are not linear were already answered, every root proved (`x*y = 6`, `x + y = 5`: (2, 3) and (3, 2)).
