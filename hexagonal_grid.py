
# This program is part of the proof for rapidly mixingness of k-heights on hexagonal grids. A block B consist of a
# single hexagon, whose assignment will be encoded as a 6-tuple in clockwise (or, by symmetry, counterclockwise)
# order. Each of these six vertices has exactly one boundary neighbor. In total, these six boundary neighbors can
# impose (k+1)^6 different boundary restrictions to B, and there are k*(k+1)^5 different cover relations that differ on
# a fixed boundary vertex d (by symmetry, this boundary vertex can be assumed to be fixed).

import itertools
import sys


print("PYTHON VERSION:\n" + sys.version)


class NoAdmissibleFilling(Exception):
    pass


def expected_weight(p_block_fillings, p_boundary_constraint):
    # computes the expected weight of a filling of B taken uniformly at random from all possible fillings (
    # block_fillings) that are compatible with the given boundary restriction (boundary_restriction).
    # p_block_fillings:         list containing all possible fillings as 6-tuples, indexed in circular order
    # p_boundary_constraint:    6-tuple of values that represent the boundary constraint, indexed in the same
    #                           circular order, i.e. the i-th value corresponds to the boundary vertex that is adjacent
    #                           to the block vertex to which the i-th value of the filling corresponds.

    # This list will contain the weights of all admissible fillings w.r.t. p_boundary_restriction. The expected weight
    # is then just the average of these values.
    possible_weights = []

    for filling in p_block_fillings:

        # Check whether this filling is admissible
        filling_is_admissible = True
        for i in range(6):
            if abs(filling[i] - p_boundary_constraint[i]) > 1:
                filling_is_admissible = False
        if filling_is_admissible:
            # Yes, it is admissible. Add its weight to possible_weights .
            possible_weights.append(sum(filling))

    if len(possible_weights) == 0:
        # This restriction cannot be fulfilled.
        raise NoAdmissibleFilling()

    # Return the average of the possible_weights .
    return sum(possible_weights) / len(possible_weights)


for k in range(2, 7):
    # Compute the block divergence of k-heights for k = 2, 3, 4, 5, 6

    # Step 1: Iterate over all 6-tuples of {0, ..., k} and keep those that are valid fillings.

    # This list will contain all internally valid fillings, i.e. the values of two successive vertices in circular
    # order differ by at most one.
    block_fillings = []
    for filling in itertools.product(range(k+1), repeat=6):
        is_valid_filling = True
        for i in range(5):
            if abs(filling[i] - filling[i+1]) > 1:
                is_valid_filling = False
        if abs(filling[0] - filling[5]) > 1:
            is_valid_filling = False
        if is_valid_filling:
            block_fillings.append(filling)

    # Step 2: Compute k*(k+1)^5 cover relations, i.e. pairs of boundary constraints that differ on exactly one vertex by one.
    #
    # By symmetry, it is sufficient to only consider cover relations on the boundary of B that differ on the neighbor
    # vertex of the first block vertex (index 0 in the tuple).

    # This list will contain all cover relations as 2-tuples of 6-tuples. It may happen that we use boundary constraints
    # that cannot be satisfied, i.e. for which there is no admissible filling. This happens when the difference of
    # boundary values corresponding to close boundary vertices is too large. However, this will just cause a
    # NoAdmissibleFilling exception later which we can handle, i.e. we simply ignore those.
    cover_relations = []
    for cover_relation_on_other_5 in itertools.product(range(k+1), repeat=5):
        for cover_relation_first_entry in range(k):
            l = tuple([cover_relation_first_entry] + list(cover_relation_on_other_5))
            u = tuple([cover_relation_first_entry + 1] + list(cover_relation_on_other_5))
            cover_relations.append((l, u))

    # Step 3: For each cover relation, compute the expected weight difference on the block and take the maximum.
    maximum_expected_weight_difference = 0
    for cover_relation in cover_relations:
        try:
            (l, u) = cover_relation
            expected_weight_conditioned_on_l = expected_weight(block_fillings, l)
            expected_weight_conditioned_on_u = expected_weight(block_fillings, u)

            expected_weight_difference = expected_weight_conditioned_on_u - expected_weight_conditioned_on_l
            if expected_weight_difference > maximum_expected_weight_difference:
                maximum_expected_weight_difference = expected_weight_difference
        except NoAdmissibleFilling:
            pass

    # Print the computed block divergence.
    print("Block divergence for k=" + str(k) + " : " + str(maximum_expected_weight_difference))
