## Test environments

* Local macOS Sonoma 14.2, R 4.4.2

## R CMD check results

`R CMD check --as-cran MALC_0.1.0.tar.gz`

Status: 3 NOTEs

## Notes

This is a new submission.

The local check reports:

* `checking for future file timestamps ... NOTE: unable to verify current time`

This appears to be local to the check environment.

The local check also reports HTML validation notes such as:

* `<main> is not recognized`
* `<script> proprietary attribute "onload"`

These are generated from R's HTML help template by the local HTML validator.
The Rd sources themselves pass the Rd checks, examples, tests, and PDF manual
checks.
