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
#
# The algorithms implemented in this file are from the following publication:
# [1]  John W. Chinneck, Erik W. Dravnieks, (1991) Locating Minimal Infeasible Constraint Sets in Linear Programs.
#      ORSA Journal on Computing 3(2):157-168.
#      https://doi.org/10.1287/ijoc.3.2.157
#
# Algorithms used, implementation and behavior were influenced by the following publication:
# [2]  John W. Chinneck, (1997) Finding a Useful Subset of Constraints for Analysis in an Infeasible Linear Program.
#      INFORMS Journal on Computing 9(2):164-174.
#      https://doi.org/10.1287/ijoc.9.2.164

from pulp import (
    LpVariable, LpAffineExpression, LpProblem, lpSum,
    LpMinimize, LpConstraintGE, LpConstraintLE, LpStatusInfeasible,
)


# Creates a deep clone of a problem. Allows having ownership of the problem rather than a borrowed reference.
# The wording above is absolutely a reference to Rust's Borrow Checker. >:D
def deep_clone_problem(problem):
    cloned_vars = {}

    clone = LpProblem(problem.name)  # No need to copy the objective.
    for constraint in problem.constraints.values():
        cloned_constraint = clone_affine_expression_in(constraint, cloned_vars)

        if constraint.sense == LpConstraintGE:
            cloned_constraint = cloned_constraint >= -constraint.constant
        elif constraint.sense == LpConstraintLE:
            cloned_constraint = cloned_constraint <= -constraint.constant
        else:
            cloned_constraint = cloned_constraint == -constraint.constant

        cloned_constraint.name = constraint.name
        clone += cloned_constraint

    return clone


def clone_affine_expression_in(expr, cloned_vars):
    cloned = LpAffineExpression()
    for (x, a) in expr.items():
        if x.name not in cloned_vars:
            cloned_vars[x.name] = LpVariable(x.name, upBound=x.upBound, lowBound=x.lowBound, cat=x.cat)

        y = cloned_vars[x.name]
        cloned += a * y

    cloned.name = expr.name
    return cloned


# Turn a problem into an elastic problem. (Linked paper [1], §3. Elastic Programming and Phase 1 LPs)
#
# This is useful not only for the elastic filter (for quite obvious reasons); but also for the deletion filter.
# The elastic filter *is* the phase 1 of the two-phase simplex method. It therefore gives us insight about the
# feasibility or infeasibility of a problem, the shadow prices gathered tell us, for an infeasible problem,
# which constraints are involved in the IIS that we encountered while solving. Reciprocally, we also know which
# constraints are definitely not involved in *this* specific IIS, and we can therefore eliminate those.
#
# The sensitivity filter therefore becomes almost zero-cost, since it's just observing side products of the
# original goal of identifying feasibility.
def make_problem_elastic(problem):
    elastic_problem = LpProblem(f"{problem.name}+elastic")
    elastic_variables = {}
    for (k, constraint) in problem.constraints.items():
        e_u = LpVariable(f"__elastic_{k}_up", 0, 0)
        e_d = LpVariable(f"__elastic_{k}_down", 0, 0)

        elastic_constraint = constraint + e_u - e_d
        elastic_constraint.name = f"{k}+elastic"

        if elastic_constraint.sense == LpConstraintGE:
            e_u.upBound = None
        elif elastic_constraint.sense == LpConstraintLE:
            e_d.upBound = None
        else:
            e_u.upBound = None
            e_d.upBound = None

        elastic_variables[k] = e_u, e_d
        elastic_problem += elastic_constraint

    # Define objective: we want to move the least amount of constraints, so we'll minimize that.
    # Note: a feasible solution must have a solution where this is zero (no artificial variable is needed).
    elastic_problem += lpSum(elastic_variables.values())
    return elastic_problem, elastic_variables


# Create a PuLP problem from a given set of constraints.
def constraints_to_problem(constraints_set):
    problem = LpProblem()
    for k, constraint in constraints_set.items():
        problem += constraint

    return problem

# Perform an elastic filtering on the problem, as defined in the linked paper [1], §4. Elastic Filtering.
#
# The elastic filter is enhanced with suggestions from the linked paper [2], §6.2.7. Finding Useful Solutions, to
# make the output of the elastic filter as useful as possible, while allowing some optimizations to be performed
# in the deletion filter (such as the use of a sensitivity filer).
#
# The output is a set of constraints that contain one or more IIS, but no irrelevant constraints; that is,
# every constraint contained in the output set is part of an IIS. ([1] §4.1., Lemma 3)
#
# From limited testing, it seems doing warm starts (i.e. reusing the results of the previous iteration)
# make zero difference in the performance of the filter, and in some cases worsen it.
def elastic_filter(problem):
    elastic_problem, elastic_variables = make_problem_elastic(problem)
    result = {}

    # Scary infinite loop!!!
    # This while true is [theoretically] safe, as each iteration will either:
    # - Have an objective value for the LP of zero, exiting the loop
    # - Remove constraints from the LP
    # As the goal is to minimize the sum of all elastic variables, it will reach
    # zero in a finite amount of iterations due to the nature of the LP we created.
    while True:
        status = elastic_problem.solve()
        if status == LpStatusInfeasible or elastic_problem.objective.value() == 0:
            return result

        # For every constraint that has been stretched, turn them back to their original non-elastic form,
        # append it to the result set, and update the elastic variables' upper bound and value to 0.
        for (k, (e_u, e_d)) in elastic_variables.items():
            if e_u.value() != 0 or e_d.value() != 0:
                # Enforce the constraint by forcing the artificial variable to be zero.
                # This is easier than to re-build the problem without the elastic variables.
                #
                # In our case (cannot re-use solver memory structures etc.) it may come at a cost,
                # since the solver will have to go through them again. Let's hope it's negligible! >:3
                e_u.upBound = e_d.upBound = 0
                result[k] = problem.constraints[k]


# Perform a deletion filtering on the problem, as defined in the linked paper [1], §2. Deletion Filtering.
# The deletion filter is enhanced by a sensitivity filter, improving its efficiency. ([1] §5. Sensitivity Filtering)
#
# The output is guaranteed to be a single IIS, if an IIS exists in the problem. ([1] §2.1., Theorem 2)
# Variable bounds are not filtered, although an extension of the deletion algorithm may be used for this
# purpose. ([1] §2.1. Theorems and Discussion)
def deletion_filter(problem):
    # It is impossible a single constraint makes up an IIS
    if len(problem.constraints) == 1:
        return {}

    elastic_problem, _ = make_problem_elastic(problem)

    constraints = problem.constraints.copy()
    elastic_constraints = elastic_problem.constraints.copy()

    # For every constraint, check if removing it makes the program feasible.
    # If it does, keep the constraint. Otherwise, drop the constraint and apply sensitivity filtering.
    to_process = elastic_problem.constraints.copy()
    while len(to_process) != 0:
        key, _ = to_process.popitem()

        sub_problem = LpProblem("feasibility_check", LpMinimize)
        sub_problem.objective = elastic_problem.objective
        for k, constraint in elastic_constraints.items():
            if k != key:
                sub_problem += constraint

        # The problem cannot be infeasible.
        # All constraints are elastic and can be freely violated.
        sub_problem.solve()

        if sub_problem.objective.value() > 0:
            del constraints[key[:-8]]
            del elastic_constraints[key]

            # Sensitivity filter
            for c2 in sub_problem.constraints.values():
                if c2.pi == 0:
                    del constraints[c2.name[:-8]]
                    del elastic_constraints[c2.name]
                    if c2.name in to_process:
                        del to_process[c2.name]

            # Same check as before, for the same reasons.
            if len(constraints) == 1:
                return constraints

    return constraints


# Find a single IIS of an infeasible problem, using an elastic filtering, followed by a deletion/sensitivity filter.
#
# The performance of this approach varies depending on the problem's properties and the size of the IIS.
# It is clear from a limited test that some problems yield an IIS very fast, and some will take their time.
#
# A more sophisticated implementation should use heuristics to determine which algorithm to use, and when.
# For example, a heuristic could be able to tell us running an elastic filtering step is inefficient based on the
# properties of the problem, leading to a faster result while still maintaining the huge speed-up (I measured up
# to 90x speed improvements for some problems such as `gosh.mps`) of having the elastic filtering pass.
#
# Combining the two, we could imagine being able based on a heuristic be able during the solving that the method
# currently picked will be inefficient, and try something else. If we also pour in the ability to optimize the
# algorithms used, there is a huge potential for improvement over this rather crude implementation.
def compute_iis(problem):
    problem = deep_clone_problem(problem)

    # Elastic filter
    filtered_constraints = elastic_filter(problem)
    filtered_problem = constraints_to_problem(filtered_constraints)

    # Deletion filter - produces an IIS (with potentially useless variable bounds)
    iis_constraints = deletion_filter(filtered_problem)
    iis_problem = constraints_to_problem(iis_constraints)

    iis_problem.name = f"iis_of_{problem.name}"
    return iis_problem
