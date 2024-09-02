
# This code is a generic brute force attack for computing the block divergence for arbitrary blocks in arbitrary graphs.
# As it iterates over all block constraints and block fillings, it may only be applied for blocks consisting of
# few vertices and few boundary vertices, otherwise the runtime will explode.


import itertools


class NoAdmissibleFilling(Exception):
    pass


class Block:

    # We use 0-based indices for the vertices.
    # k:                                        The k of the k-heights
    # number_vertices:                          Number of vertices within the block. Vertex indices
    #                                           are 0, ..., number_vertices - 1
    # edges:                                    Edges between vertices in the block, given as a list of tuples
    #                                           of vertex indices.
    # number_boundary_vertices:                 Number of vertices on the boundary, boundary vertex indices are
    #                                           0, ..., number_boundary_vertices - 1
    # boundary_edges:                           Edges between vertices in the boundary, given as a list of tuples
    #                                           of boundary vertex indices
    # block_boundary_edges:                     Edges between a block vertex and a boundary vertex, given as a list
    #                                           of tuples of a block vertex index and a boundary vertex index

    def __init__(self, k, number_vertices, edges, number_boundary_vertices, boundary_edges, block_boundary_edges):

        self.k = k
        self.number_vertices = number_vertices
        self.edges = edges
        self.number_boundary_vertices = number_boundary_vertices
        self.boundary_edges = boundary_edges
        self.block_boundary_edges = block_boundary_edges

        self.fillings = self.compute_all_k_heights(self.k, self.number_vertices, self.edges)
        self.boundary_constraints = self.compute_all_k_heights(self.k, self.number_boundary_vertices,
                                                                      self.boundary_edges)

        print("Number of fillings: " + str(len(self.fillings)))
        print("Number of boundary constraints: " + str(len(self.boundary_constraints)))

    @staticmethod
    def compute_all_k_heights(k, n, edges):

        result = []

        for k_height in itertools.product(range(k + 1), repeat=n):
            is_valid_k_height = True
            for edge in edges:
                if abs(k_height[edge[0]] - k_height[edge[1]]) > 1:
                    is_valid_k_height = False
                    break
            if is_valid_k_height:
                result.append(k_height)

        return result

    def compute_expected_weight(self, boundary_constraint):

        total_weight_sum = 0
        number_admissible_fillings = 0

        for filling in self.fillings:
            is_admissible = True
            for block_boundary_edge in self.block_boundary_edges:
                if abs(filling[block_boundary_edge[0]] - boundary_constraint[block_boundary_edge[1]]) > 1:
                    is_admissible = False
                    break
            if is_admissible:
                total_weight_sum += sum(filling)
                number_admissible_fillings += 1

        if number_admissible_fillings == 0:
            # no admissible filling exists with respect to boundary_constraint
            raise NoAdmissibleFilling()

        return total_weight_sum / number_admissible_fillings

    # computes the expected weight difference when doing an augmentation at boundary vertex augmentation_vertex

    def block_divergence(self, augmentation_vertex):

        adjacent_boundary_vertices = []
        for boundary_edge in self.boundary_edges:
            if boundary_edge[0] == augmentation_vertex:
                adjacent_boundary_vertices.append(boundary_edge[1])
            if boundary_edge[1] == augmentation_vertex:
                adjacent_boundary_vertices.append(boundary_edge[0])

        max_expected_difference = 0

        for boundary_constraint in self.boundary_constraints:

            can_be_augmented = True

            if boundary_constraint[augmentation_vertex] == self.k:
                can_be_augmented = False

            for boundary_neighbor in adjacent_boundary_vertices:
                if (boundary_constraint[augmentation_vertex] - boundary_constraint[boundary_neighbor]) == 1:
                    can_be_augmented = False

            if can_be_augmented:

                augmented_boundary_constraint = list(boundary_constraint)
                augmented_boundary_constraint[augmentation_vertex] += 1
                augmented_boundary_constraint = tuple(augmented_boundary_constraint)

                try:
                    expected_weight_without_augmentation = self.compute_expected_weight(boundary_constraint)
                    expected_weight_with_augmentation = self.compute_expected_weight(augmented_boundary_constraint)

                    expected_weight_difference = expected_weight_with_augmentation - expected_weight_without_augmentation
                    if expected_weight_difference > max_expected_difference:
                        max_expected_difference = expected_weight_difference

                except NoAdmissibleFilling:
                    pass

        if max_expected_difference == 0:
            raise Exception("Augmentation does not result in an increase of expected weight.")

        return max_expected_difference