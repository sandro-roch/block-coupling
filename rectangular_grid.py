
# This program is part of the proof for rapidly mixingness of 3-heights on rectangular grids. A block B is a
# (4x4)-square subgrid. Its boundary consists of 16 vertices that lie on four sides: front, left, back,
# right. Therefore, we iterate over all boundary restrictions on front, left, back and right. When computing the
# expected weight of a filling of B, we count both the number of ways to extend the boundary to a filling of B and
# the sum of the weights of these possible fillings. For both tasks we use a dynamic program which counts the
# admissible partial fillings from the front up to the i-th row; i = 1, 2, 3, 4.


import itertools
import random
import sys


print("PYTHON VERSION:\n" + sys.version)


def expected_weight(p_strips, p_boundary_front, p_boundary_left, p_boundary_right, p_boundary_back):
    # Computes the expected weight of a block filling uniformly chosen at random from all fillings that are compatible
    # with the boundary restriction.
    # p_strips: all k-strips, i.e. k-heights on a path of length 4
    # p_boundary_front: lower side of the boundary restriction, indexed from left to right
    # p_boundary_back: upper side of the boundary restriction, indexed from left to right
    # p_boundary_left: left side of the boundary restriction, index from front to back
    # p_boundary_right: right side of the boundary restriction, indexed from front to back

    # We use a dynamic programming approach and count partial fillings consisting only of the i lowest rows,
    # i.e. fillings of the i rows closest to the front. More precisely, we do not only count their number,
    # but also their total sum of weights. After i iterations, the following dictionary contains number and weight of
    # all partial fillings with given last row, i.e. given i-th row behind the front. For initialization,
    # after running 0 iterations, the last row is the front itself, there is exactly one "empty partial extension"
    # that has no weight yet.

    extensions_by_last_row = {
        p_boundary_front: (1, 0)
    }

    # We run 4 iterations for constructing the first, second, third and fourth row of the filling. For the fourth
    # row, we have to respect the boundary restriction on the back as well.

    for row in range(4):

        # Will be the new version of extension_by_last_row after this iteration
        extensions_by_new_last_row = dict()

        # Iterate over all possibilities for the previous row
        for prev_last_row in extensions_by_last_row:

            # iterate over all (internally valid) possibilities for extending the previous row by a further row
            for new_last_row in p_strips:

                # The new row is internally valid, but we have to ensure that it is also compatible with the previous
                # row and with the boundary restrictions
                is_valid_extension = True

                for j in range(4):
                    if abs(new_last_row[j] - prev_last_row[j]) > 1:
                        # This is not a possible extension because of a conflict with the previous row
                        is_valid_extension = False

                if abs(new_last_row[0] - p_boundary_left[row]) > 1:
                    # This is not a possible extension because of a conflict with the left boundary
                    is_valid_extension = False

                if abs(new_last_row[3] - p_boundary_right[row]) > 1:
                    # This is not a possible extension because of a conflict with the right boundary
                    is_valid_extension = False

                if row == 3:
                    # We are about to construct the last row where we also have to check the backside of the boundary
                    for j in range(3):
                        if abs(new_last_row[j] - p_boundary_back[j]) > 1:
                            # This is not a possible extension because of a conflict with the backside of the boundary
                            is_valid_extension = False

                if is_valid_extension:

                    # It is possible to extend a partial filling whose last row is old_last_row to a partial
                    # filling whose next row is new_last_row. Now update the number and weight of the partial
                    # fillings that have one more row and end in new_last_row!

                    if new_last_row not in extensions_by_new_last_row:
                        extensions_by_new_last_row[new_last_row] = (0, 0)

                    new_number = extensions_by_new_last_row[new_last_row][0] + extensions_by_last_row[prev_last_row][0]
                    new_weight = extensions_by_new_last_row[new_last_row][1] \
                                 + extensions_by_last_row[prev_last_row][1] \
                                 + extensions_by_last_row[prev_last_row][0] * sum(new_last_row)

                    extensions_by_new_last_row[new_last_row] = (new_number, new_weight)

        # ready for the next iteration
        extensions_by_last_row = extensions_by_new_last_row

    # Now compute the expected weight of a filling simply by dividing the total weight of all admissible fillings by
    # their number.

    total_number_of_fillings = 0
    total_weight_of_fillings = 0

    for last_row in extensions_by_last_row:
        total_number_of_fillings += extensions_by_last_row[last_row][0]
        total_weight_of_fillings += extensions_by_last_row[last_row][1]

    return total_weight_of_fillings / total_number_of_fillings


for k in range(2, 4):

    print("CASE k = " + str(k))

    # Step 1: Compute all k-strips, i.e. k-heights on a path of length 4
    print("Computing k-strips...")

    k_strips = []
    for strip in itertools.product(range(k + 1), repeat=4):
        is_valid_k_strip = True
        for i in range(3):
            if abs(strip[i] - strip[i + 1]) > 1:
                is_valid_k_strip = False
        if is_valid_k_strip:
            k_strips.append(strip)

    # Step 2: Compute all possible boundary restrictions, encoded as 16-tuples in counterclockwise order starting with
    # front side
    print("Computing boundaries...")

    boundary_constraints = []
    for sides in itertools.product(k_strips, repeat=4):
        # Two boundary vertices that have distance 2 can only differ by at most two. Hence, we can ignore the boundary
        # constraints in which this property does not hold.
        if abs(sides[0][3] - sides[1][0]) > 2:
            continue
        if abs(sides[1][3] - sides[2][0]) > 2:
            continue
        if abs(sides[2][3] - sides[3][0]) > 2:
            continue
        if abs(sides[3][3] - sides[0][0]) > 2:
            continue
        boundary_constraints.append(sides[0] + sides[1] + sides[2] + sides[3])


    # Step 3: Consider all cover relations where the differing vertex is either the first or the second vertex of the front
    # boundary (other cases are symmetric). Compute the expected weight difference and take the maximum over all of these.
    print("Computing expected weight differences on cover relations...")

    # We will maximize over all boundary_constraints. Shuffling their order does not affect the final result. One could
    # actually remove it. However, if we shuffle their order randomly one can quickly "see" or guess the value that is
    # being approached without the need running the process till it terminates. However, for the proof  of rapid mixing
    # in the case of k = 2, 3 in our article we only use the final result value. However, for lower bounding the block
    # divergence in the case of k = 4, we do interrupt the algorithm.
    random.shuffle(boundary_constraints)

    max_block_divergence = 0

    for r in range(len(boundary_constraints)):

        boundary_constraint = boundary_constraints[r]
        if (r % 1000) == 0:
            print("Processing " + str(r) + "th boundary constraint out of " + str(len(boundary_constraints)) + " in total. Max block divergence found so far: " + str(max_block_divergence))

        # split counterclockwise indexed 16-tuple into front, left, right and back side.

        boundary_front = boundary_constraint[0:4]
        boundary_right = boundary_constraint[4:8]
        boundary_back = boundary_constraint[11:7:-1]
        boundary_left = boundary_constraint[15:11:-1]

        # augment 0-th front vertex if possible

        if boundary_front[0] != k:
            augmented_boundary_front_as_list = list(boundary_front)
            augmented_boundary_front_as_list[0] += 1
            augmented_boundary_front = tuple(augmented_boundary_front_as_list)

            if not (abs(augmented_boundary_front[0] - augmented_boundary_front[1]) > 1) \
                    and not (abs(augmented_boundary_front[0] - boundary_left[0]) > 2):
                # is still a valid boundary after augmentation

                expected_weight_without_augmentation = expected_weight(k_strips, boundary_front, boundary_left,
                                                                       boundary_right, boundary_back)
                expected_weight_with_augmentation = expected_weight(k_strips, augmented_boundary_front, boundary_left,
                                                                    boundary_right, boundary_back)

                expected_weight_difference = expected_weight_with_augmentation - expected_weight_without_augmentation
                if expected_weight_difference > max_block_divergence:
                    max_block_divergence = expected_weight_difference

        # augment 1-th front vertex if possible

        if boundary_front[1] != k:
            augmented_boundary_front_as_list = list(boundary_front)
            augmented_boundary_front_as_list[1] += 1
            augmented_boundary_front = tuple(augmented_boundary_front_as_list)

            if not (abs(augmented_boundary_front[1] - augmented_boundary_front[2]) > 1) \
                    and not (abs(augmented_boundary_front[0] - augmented_boundary_front[1]) > 1):
                # is still a valid boundary after augmentation

                expected_weight_without_augmentation = expected_weight(k_strips, boundary_front, boundary_left,
                                                                       boundary_right, boundary_back)
                expected_weight_with_augmentation = expected_weight(k_strips, augmented_boundary_front, boundary_left,
                                                                    boundary_right, boundary_back)

                expected_weight_difference = expected_weight_with_augmentation - expected_weight_without_augmentation
                if expected_weight_difference > max_block_divergence:
                    max_block_divergence = expected_weight_difference

    print("BLOCK DIVERGENCE IN CASE k = " + str(k) + ": " + str(max_block_divergence))
