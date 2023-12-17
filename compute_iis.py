# -*- coding: utf-8 -*-
#
# Copyright (c) Cynthia Rey, All rights reserved.
# SPDX-License-Identifier: BSD-3-Clause
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice, this
# list of conditions and the following disclaimer.
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution.
# 3. Neither the name of the copyright holder nor the names of its contributors
# may be used to endorse or promote products derived from this software without
# specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
# WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
# ----------
#
# The algorithms implemented in this file are from the following publication:
#     John W. Chinneck, Erik W. Dravnieks, (1991) Locating Minimal Infeasible Constraint Sets in Linear Programs.
#     ORSA Journal on Computing 3(2):157-168.
#     https://doi.org/10.1287/ijoc.3.2.157

from pulp import (
    LpVariable, LpProblem, lpSum,
    LpMinimize, LpConstraintGE, LpConstraintLE, LpStatusInfeasible,
)


# Makes a constraint elastic by adding 2 variables to move the constraint up or down.
# The bounds of the elastic variables are set to only allow the constraint to be relaxed.
def make_constraint_elastic(constraint):
    e_u = LpVariable(f"$elastic_{constraint.name}_up", 0, 0)
    e_d = LpVariable(f"$elastic_{constraint.name}_down", 0, 0)
    elastic_constraint = constraint + e_u - e_d
    elastic_constraint.name = f"{constraint.name}+elastic"

    if elastic_constraint.sense == LpConstraintGE:
        e_u.upBound = None
    elif elastic_constraint.sense == LpConstraintLE:
        e_d.upBound = None
    else:
        e_u.upBound = None
        e_d.upBound = None

    return elastic_constraint, e_u, e_d


# Perform an elastic filtering on the problem, as defined in the linked paper, §4. Elastic filtering.
#
# The output is a set of constraints that contain one or more IIS, but no irrelevant constraints; that is,
# every constraint contained in the output set is part of an IIS. (§4.1., Lemma 3)
#
# From limited testing, it seems doing warm starts (i.e. reusing the results of the previous iteration)
# make zero difference in the performance of the filter, and in some cases worsen it.
def elastic_filter(problem):
    elastic_problem = LpProblem("elastic_filtering", LpMinimize)

    # Turn all the constraints into "elastic" constraints.
    # They then can be "stretched", or moved, to try and make a problem that's solvable.
    elastic_variables = dict()
    for (k, constraint) in problem.constraints.items():
        elastic_constraint, e_u, e_d = make_constraint_elastic(constraint)
        elastic_constraint.__original = constraint

        elastic_variables[elastic_constraint.name] = e_u, e_d
        elastic_problem += elastic_constraint

    # Define objective
    # We want to move the least amount of constraints, so we'll minimize that.
    elastic_problem += lpSum(elastic_variables.values())

    # Scary infinite loop!!!
    # This while true is [theoretically] safe, as each iteration will either:
    # - Have an objective value for the ILP of zero, exiting the loop
    # - Remove constraints from the ILP
    # As the goal is to minimize the sum of all elastic variables, it will reach
    # zero in a finite amount of iterations due to the nature of the ILP we created.
    result = []
    while True:
        status = elastic_problem.solve()
        if status == LpStatusInfeasible or elastic_problem.objective.value() == 0:
            return result

        # For every constraint that has been stretched, turn them back to their original non-elastic form,
        # append it to the result set, and update the elastic variables' upper bound and value to 0.
        for (k, (e_u, e_d)) in elastic_variables.items():
            constraint = elastic_problem.constraints[k]
            if constraint.pi != 0:
                # Enforce the constraint by not allowing it to be stretched
                e_u.upBound = e_d.upBound = 0
                e_u.setInitialValue(0)
                e_d.setInitialValue(0)
                result.append(constraint.__original)


# Perform a deletion filtering on the problem, as defined in the linked paper, §2. Deletion Filtering.
# The deletion filter is enhanced by a sensitivity filter, improving its efficiency (§5. Sensitivity Filtering)
#
# The output is guaranteed to be a single IIS, if an IIS exists in the problem. (§2.1., Theorem 2)
#
# The function accepts a list of filtered constraints, to allow enhancing the speed of the algorithm by only
# operating on the constraints identified by an elastic filter.
def deletion_filter(problem, filtered_constraints=None):
    constraints = list(problem.constraints.values() if filtered_constraints is None else filtered_constraints)

    # If there is a single contraint, we don't need to process any further.
    # Plus, it seems CBC has issues with problems that include zero constraints (who would do that! totally not me-)
    if len(constraints) == 1:
        return constraints

    # For every constraint, check if removing it makes the program feasible.
    # If it does, keep the constraint. Otherwise, drop the constraint.
    # List is copied to make sure deletions are handled correctly.
    to_process = list(constraints)
    while len(to_process) != 0:
        c = to_process.pop()

        sub_problem = LpProblem("feasibility_check", problem.sense)
        sub_problem += problem.objective
        for constraint in constraints:
            if c != constraint:
                sub_problem += constraint

        # Efficiency: this is... not great. We could shortcircuit if we had control over the internals of
        # the solver as soon as we identify if the program will be feasible, or not. No need to go and try to
        # find the optimal solution. For instance, a solver using the two-phase simplex method could stop right
        # after phase 1.
        #
        # That's the downside of doing it from PuLP: we don't have the ability to poke at internal data
        # structures of the solver, or have granular control over what it's doing. It is what it is!
        if sub_problem.solve() == LpStatusInfeasible:
            constraints.remove(c)
            for c in sub_problem.constraints.values():
                if c.pi != 0:
                    constraints.remove(c)
                    to_process.remove(c)

            # Same check as before, for the same reasons.
            if len(constraints) == 1:
                return constraints

    return constraints


# Find a single IIS of an infeasible problem, using an elastic filtering, followed by a deletion/sensitivity filter.
#
# The performance of this approach varies depending on the problem's properties and the size of the IIS.
# It is clear from a limited test that some problems yield an IIS very fast, and some will take their time.
# This is partly due to the inefficient way this is done by interfacing with a solver through a high level modeler,
# which does not let us access the solver internal states nor do partial solves.
#
# In addition, a more sophisticated implementation should use heuristics to determine which algorithm to use, and when.
# For example, a heuristic could be able to tell us running an elastic filtering step is inefficient based on the
# properties of the problem, leading to a faster result while still maintaining the huge speed-up (I measured up
# to 90x speed improvements for some problems such as `gosh.mps`) of having the elastic filtering pass.
#
# Combining the two, we could imagine being able based on a heuristic be able during the solving that the method
# currently picked will be inefficient, and try something else. If we also pour in the ability to optimize the
# algorithms used, there is a huge potential for improvement over this rather crude implementation.
def compute_iis(problem):
    filtered_constraints = elastic_filter(problem)
    return deletion_filter(problem, filtered_constraints)
