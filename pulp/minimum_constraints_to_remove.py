# -*- coding: utf-8 -*-
#
# Copyright (c) Cynthia Rey, All rights reserved.
# SPDX-License-Identifier: MPL-2.0
#
# This Source Code Form is subject to the terms of the Mozilla Public
# License, v. 2.0. If a copy of the MPL was not distributed with this
# file, You can obtain one at http://mozilla.org/MPL/2.0/.
#
# Formulas used in this file come from this post on Operations Research StackExchange
#     https://or.stackexchange.com/a/7309

from pulp import LpVariable, LpProblem, LpMinimize, LpInteger, lpSum, LpConstraintGE, LpConstraintLE


def extract_constraint_data(constraint, variables):
    return [constraint[v] if v in constraint else 0 for v in variables], constraint.constant


def compute_minimum_set_of_constraints_to_remove(problem):
    prob = LpProblem("min_constraint_deletion", LpMinimize)

    constraint_count = len(problem.constraints)
    variables_count = len(problem.variables())

    x = problem.variables()

    # Transform the problem in something that is of the form `Ax <= b`, and extract A and b.
    # A is our `coefficients_mx`, and b our `constraint_map`.
    #
    # When we have something of the form Ax >= b, we multiply the inequality by -1 to flip the inequality.
    # In case of an equality, the constraint is duplicated, so we have to satisfy both:
    #
    #    Ax <= b
    #    Ax >= b = -Ax <= -b
    #
    # Which is only possible if Ax == b.
    coefficients_mx = []
    constants_list = []
    constraint_map = []
    for (k, constraint) in problem.constraints.items():
        if constraint.sense == LpConstraintLE:
            constraint_map.append(k)
            (a, b) = extract_constraint_data(constraint, x)
            coefficients_mx.append(a)
            constants_list.append(b)
        elif constraint.sense == LpConstraintGE:
            constraint_map.append(k)
            (a, b) = extract_constraint_data(constraint * -1, x)
            coefficients_mx.append(a)
            constants_list.append(b)
        else:
            constraint_count += 1
            constraint_map.append(k)
            constraint_map.append(k)
            (a1, b1) = extract_constraint_data(constraint * -1, x)
            (a2, b2) = extract_constraint_data(constraint, x)
            coefficients_mx.append(a1)
            coefficients_mx.append(a2)
            constants_list.append(b1)
            constants_list.append(b2)

    z = LpVariable.dict("$$z", range(0, constraint_count), 0, 1, LpInteger)
    u = LpVariable.dict("$$u", range(0, variables_count), 0, 1, LpInteger)
    v = LpVariable.dict("$$v", range(0, variables_count))

    # Create the auxiliary problem
    prob += lpSum(z) + lpSum(u)

    for i in range(0, constraint_count):  # Constraint 1
        expression = 0
        for j in range(0, variables_count):
            expression += coefficients_mx[i][j] * x[j]

        # Constraint count is used as Big-M here, for the lack of a better one :shrug:
        prob += expression + constants_list[i] <= constraint_count * z[i]

    for j in range(0, variables_count):  # Constraint 2
        prob += -u[j] <= x[j] - v[j]
        prob += x[j] - v[j] <= u[j]

    # Solve the problem
    if prob.solve() != 1:
        return None

    result = set()
    for i in range(0, constraint_count):
        if z[i].value() == 1:
            result.add(constraint_map[i])

    return list(result)
