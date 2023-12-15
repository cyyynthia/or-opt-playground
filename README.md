# PuLP playground
Some of my experiments playing with linear programs and PuLP, to make things that are only available in commercial
solvers available in open-source world 🫡

Got inspired from [flopEDT](https://flopedt.org/) as they've been struggling to have a way to compute IIS of their
integer linear programs, so I started toying with PuLP after reading some published paper on the matter.

The code quality is not very good, and I'd recommend to double-check it before using it in any kind of production app.
The snippets in this repository work with CBC MILP. Other solvers are untested, but they *should* work with any
solver.

## Things I've done
- [minimum_constraints_to_remove.py](./minimum_constraints_to_remove.py):
Computes a possible minimal combinaison of constraints to remove in order to make a problem feasible.
Formulas used from this post on Operations Research StackExchange: https://or.stackexchange.com/a/7309

- [compute_iis.py](./compute_iis.py):
Computes an IIS of constraints from an infeasible problem. Uses the *Elastic Filtering* algorithm, from
[*John W. Chinneck, Erik W. Dravnieks, (1991) Locating Minimal Infeasible Constraint Sets in Linear Programs. ORSA Journal on Computing 3(2):157-168.*](https://doi.org/10.1287/ijoc.3.2.157)

## License
Software licensed under the BSD-3-Clause license, please see LICENSE for more details. Each file includes credits
and information about where it's from when appropriate.

---

Fuck Python 🩷
