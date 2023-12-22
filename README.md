# Cynthia's OR optimization playground
My personal playground where I experiment with optimization problems, mostly linear programs.

Got inspired* from [flop!EDT](https://flopedt.org/) as they've been struggling to have a way to compute IIS of their
linear programs, so I started toying with PuLP after reading some published paper on the matter.

*\*by inspired, the author means Pablo Seban from flop!EDT proposed a research question about this as part of a course
at the Toulouse III - Paul Sabatier University, which lead to curiosity taking over and... this. hehe.*

The code quality is not very good, and is certainly not production grade. They're mostly proof of concepts, and full
of misc problems that would need to be addressed in a real solution.

I left a ton of comments in the code along my way, with references, notes, and more. They contain a lot of useful
information about why I did things this way and link to resources which helped me, which of course contains a lot
of further readings.

## Things I've done
- [pulp/compute_iis.py](./pulp/compute_iis.py):
Computes an IIS of constraints of an infeasible problem, using PuLP. Uses a combination of elastic filtering, deletion
filtering, and sensitivity filtering. Not the most efficient solution, but an efficient implementation must come from
something lower level than something like PuLP anyway. This is mostly intended to be a PoC of an implementation.
  - [*John W. Chinneck, Erik W. Dravnieks, (1991) Locating Minimal Infeasible Constraint Sets in Linear Programs. ORSA Journal on Computing 3(2):157-168.*](https://doi.org/10.1287/ijoc.3.2.157)
  - [*John W. Chinneck, (1997) Finding a Useful Subset of Constraints for Analysis in an Infeasible Linear Program. INFORMS Journal on Computing 9(2):164-174.*](https://doi.org/10.1287/ijoc.9.2.164)

- [pulp/minimum_constraints_to_remove.py](./pulp/minimum_constraints_to_remove.py):
Computes a possible minimal combinaison of constraints to remove in order to make a problem feasible. Crude
implementation without many experiments beyond the original PoC.
  - Formulas used from this post on Operations Research StackExchange: https://or.stackexchange.com/a/7309

## Useful resources
- Collection of instances from MIPLIB: https://git.zib.de/miplib2017/revised-submissions-final
- Collection of infeasible instances (netlib): https://github.com/coin-or-tools/Data-Infeas

## License
Software licensed under the BSD-3-Clause license, please see LICENSE for more details.
Each file includes copyright information and the relevant credits.

---

F@#k Python, nya~ 🩷
