<!--
Delete anything that does not apply. A one-line fix does not need
paragraphs; a new space or algorithm does. Nothing here is a hoop -- each
line is something a reviewer would otherwise have to ask you for, and
asking costs you a round trip.
-->

## What this changes

<!-- What is different afterwards, in a sentence or two. -->

## Why

<!--
The reasoning, not the diff. If this fixes something, what was the wrong
behaviour and what did it cost? If it adds something, what could not be
expressed before?

If you measured anything, put the numbers here, including the ones that
did not go your way -- this repository's own history is mostly plausible
claims that did not survive being run, and a number that surprised you is
the most useful thing in a PR.
-->

## How it was checked

<!--
`ctest` passing is the floor, not the answer. For a bug fix, the useful
sentence is that the new test *fails* without the fix -- a test that
passes either way is not testing the fix.

For anything touching a renderer or a demo, say whether the image
changed, and if it should not have, that it did not.
-->

- [ ] `ctest` is green (`nix develop` then `cmake --preset release`, or see CONTRIBUTING.md)
- [ ] New behaviour has a test; a fix has a test that fails without it
- [ ] Docs touched if the change makes anything in them untrue

## Anything you are unsure about

<!--
Genuinely optional, and genuinely welcome. "I could not work out where
this belongs", "I think this is right but the tolerance is a guess" --
these get a better answer than silence does, and they are not held
against the PR.
-->
