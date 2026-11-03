# IR: the field

## Where it stands

A quantity is one `Field` (`io/field.hpp`): numbers, the chart's `u` and `v`, named symbols, and operations -- arithmetic, min/max, less, sin/cos, sqrt, log/exp, abs, any power (`Pow`), noise, gather, instance inputs, and a net of one hidden tanh layer whose weights are data (`Net`). It is written as text and LaTeX (`io/field_text.hpp`; a minus written once) and read back (`io/field_parse.hpp`, the language as tables the page's reference is written from). It lowers to plain data run by one interpreter in any arithmetic (`io/field_pod.hpp`).

## Now

- the language read from MathJSON, the output of an existing LaTeX editor (MathLive), so input is typed as LaTeX and no grammar of our own grows
- refusals as data: what, where, why, and what to write instead -- for a person editing and for automation retrying

## Later

- vectors and matrices as values, once the matrices-and-systems fork is decided
- interval arithmetic, for a certificate that bounds rather than samples

The entries of this area, by date, and what they did not do: `python3 scripts/gen_roadmap.py --area ir`.
