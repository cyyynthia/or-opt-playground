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
# The algorithm implemented in this file is from the following publication
#     John W. Chinneck, Erik W. Dravnieks, (1991) Locating Minimal Infeasible Constraint Sets in Linear Programs.
#     ORSA Journal on Computing 3(2):157-168.
#     https://doi.org/10.1287/ijoc.3.2.157

from pulp import LpVariable, LpProblem, LpMinimize, lpSum, LpConstraintGE, LpConstraintLE


def compute_iis(problem):
    iis_problem = LpProblem("iis_calculation_problem", LpMinimize)
    elastic_variables = dict()

    for (k, constraint) in problem.constraints.items():
        e_u = LpVariable(f"$elastic_{k}_up", 0, 0)
        e_d = LpVariable(f"$elastic_{k}_down", 0, 0)
        elastic_constraint = constraint + e_u - e_d
        elastic_constraint.name = f"${k}+elastic"
        elastic_constraint.__original = constraint

        if elastic_constraint.sense == LpConstraintGE:
            e_u.upBound = None
        elif elastic_constraint.sense == LpConstraintLE:
            e_d.upBound = None
        else:
            e_u.upBound = None
            e_d.upBound = None

        elastic_variables[elastic_constraint.name] = [e_u, e_d]
        iis_problem += elastic_constraint

    result = []

    # Scary infinite loop!!!
    # This while true is [theoretically] safe, as each iteration will either:
    # - Have an objective value for the ILP of zero, exiting the loop
    # - Remove constraints from the ILP
    # As the goal is to minimize the sum of all elastic variables, it will reach
    # zero in a finite amount of iterations due to the nature of the ILP we created.
    while True:
        iis_problem += lpSum(elastic_variables.values())
        status = iis_problem.solve()
        if status != 1:
            return None

        if iis_problem.objective.value() == 0:
            return result

        reduced_iis_problem = LpProblem("iis_calculation_problem", LpMinimize)
        for (k, constraint) in iis_problem.constraints.items():
            [e_u, e_d] = elastic_variables[k]
            if e_u.value() > 0 or e_d.value() > 0:
                result.append(constraint)
                del elastic_variables[k]
            else:
                reduced_iis_problem += constraint

        iis_problem = reduced_iis_problem
