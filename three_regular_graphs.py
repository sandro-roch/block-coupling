
# This program is part of the proof for rapidly mixingness of 3-heights on 3-regular planar graphs. The blocks are
# described in the article. This script goes for k=2 and k=3 through all cases described in the article and bounds
# the block divergence. In the output, each line corresponds to one of these computations. We use a general framework
# for computing the block divergence of small blocks which can be found in generic_block.py . For transparency,
# with each computation, we not just log the results but also the parameters that are used when calling the
# computation. This includes the graph structure of the block.

import math
import sys
import generic_block


print("PYTHON VERSION:\n" + sys.version)

# column description line for reading the output as a CSV
#
print("k; case; number of internally valid fillings; number of valid boundary constraints; result (block divergence); "
      "input number_vertices; input edges; input number_boundary_vertices; input boundary_edges; input "
      "block_boundary_edges")

for k in range(2, 4):

    #
    # TYPE 1.X.1 for 3 <= X <= 10
    #

    for x in range(3, 11):
        number_vertices = x
        edges = [(i, i + 1) for i in range(0, x - 1)] + [(x - 1, 0)]
        number_boundary_vertices = x
        boundary_edges = []
        block_boundary_edges = [(i, i) for i in range(0, x)]
        augmentation_vertex = 0

        block = generic_block.Block(k, number_vertices, edges, number_boundary_vertices, boundary_edges, block_boundary_edges)
        result = block.block_divergence(augmentation_vertex)

        print(f'{k}; Type 1.{x}.1; {len(block.fillings)}; {len(block.boundary_constraints)}; {result}; '
              + f'{number_vertices}; {edges}; {number_boundary_vertices}; {boundary_edges}; {block_boundary_edges}')

    #
    # TYPE 1.X.1Y for 3 <= X <= 10, 2 <= Y <= floor(X / 2)
    #

    for x in range(3, 11):
        number_vertices = x
        edges = [(i, i + 1) for i in range(0, x - 1)] + [(x - 1, 0)]
        number_boundary_vertices = x - 1
        augmentation_vertex = 0
        boundary_edges = []

        for y in range(2, math.floor(x / 2) + 2):
            block_boundary_edges = [(0, 0), (y - 1, 0)]
            for i in range(1, y - 1):
                block_boundary_edges.append((i, i))
            for i in range(y, x):
                block_boundary_edges.append((i, i - 1))

            block = generic_block.Block(k, number_vertices, edges, number_boundary_vertices, boundary_edges, block_boundary_edges)
            result = block.block_divergence(augmentation_vertex)

            print(f'{k}; Type 1.{x}.1{y}; {len(block.fillings)}; {len(block.boundary_constraints)}; {result}; '
                  + f'{number_vertices}; {edges}; {number_boundary_vertices}; {boundary_edges}; {block_boundary_edges}')

    #
    # TYPE 1.L.XYZ WHERE 1 = X < Y < Z < L <= 10 AND (Y - X) <= (Z - Y) <= (L + 1 - Z)
    #

    for L in range(3, 11):

        number_vertices = L
        edges = [(i, i + 1) for i in range(0, L - 1)] + [(L - 1, 0)]
        number_boundary_vertices = L - 2
        boundary_edges = []
        augmentation_vertex = 0

        X = 1

        for Y in range(2, math.floor(L / 2) + 2):
            for Z in range(Y + 1, L + 1):

                if not ((Y - X) <= (Z - Y)):
                    continue

                if not ((Z - Y) <= (L + 1 - Z)):
                    continue

                block_boundary_edges = [(X - 1, 0), (Y - 1, 0), (Z - 1, 0)]

                for i in range(X, Y - 1):
                    block_boundary_edges.append((i, i))

                for i in range(Y, Z - 1):
                    block_boundary_edges.append((i, i - 1))

                for i in range(Z, L):
                    block_boundary_edges.append((i, i - 2))

                block = generic_block.Block(k, number_vertices, edges, number_boundary_vertices, boundary_edges, block_boundary_edges)
                result = block.block_divergence(augmentation_vertex)

                #print("Type 1." + str(L) + "." + str(X) + str(Y) + str(Z) + ": " + str(result))

                print(f'{k}; Type 1.{L}.{X}{Y}{Z}; {len(block.fillings)}; {len(block.boundary_constraints)}; {result}; '
                      + f'{number_vertices}; {edges}; {number_boundary_vertices}; {boundary_edges}; {block_boundary_edges}')

    #
    # TYPE 2.X FOR 1 <= X <= 4
    #

    number_vertices = 8
    number_boundary_vertices = 10
    edges = [(i, i+1) for i in range(7)]
    boundary_edges = []
    block_boundary_edges = [(0, 0), (7, 9)] + [(i - 1, i) for i in range(1, 9)]

    for X in range(1, 5):
        augmentation_vertex = X

        block = generic_block.Block(k, number_vertices, edges, number_boundary_vertices, boundary_edges, block_boundary_edges)
        result = block.block_divergence(augmentation_vertex)

        print(f'{k}; Type 2.{X}; {len(block.fillings)}; {len(block.boundary_constraints)}; {result}; '
              + f'{number_vertices}; {edges}; {number_boundary_vertices}; {boundary_edges}; {block_boundary_edges}')

    #
    # TYPE 2.XY FOR 1 <= X < Y <= 8 WHERE (X - 1) <= (8 - Y)
    #

    number_vertices = 8
    number_boundary_vertices = 9
    edges = [(i, i + 1) for i in range(7)]
    boundary_edges = []

    for X in range(1, 5):
        for Y in range(X + 1, 9):

            if (8 - Y) < (X - 1):
                continue

            augmentation_vertex = X
            block_boundary_edges = [(0, 0), (7, 8), (X - 1, X), (Y - 1, X)]
            for i in range(0, X-1):
                block_boundary_edges.append((i, i + 1))
            for i in range(X, Y - 1):
                block_boundary_edges.append((i, i + 1))
            for i in range(Y, 8):
                block_boundary_edges.append((i, i))

            block = generic_block.Block(k, number_vertices, edges, number_boundary_vertices, boundary_edges,
                                        block_boundary_edges)
            result = block.block_divergence(augmentation_vertex)

            print(f'{k}; Type 2.{X}{Y}; {len(block.fillings)}; {len(block.boundary_constraints)}; {result}; '
                  + f'{number_vertices}; {edges}; {number_boundary_vertices}; {boundary_edges}; {block_boundary_edges}')

    #
    # TYPE 2.XYZ FOR 1 <= X < Y < Z <= 8 WHERE (X - 1) <= (8 - Z)
    #

    number_vertices = 8
    number_boundary_vertices = 8
    edges = [(i, i + 1) for i in range(7)]
    boundary_edges = []

    for X in range(1, 5):
        for Y in range(X + 1, 8):
            for Z in range(Y + 1, 9):

                if (8 - Z) < (X - 1):
                    continue

                augmentation_vertex = X
                block_boundary_edges = [(0, 0), (7, 7), (X - 1, X), (Y - 1, X), (Z - 1, X)]

                for i in range(0, X - 1):
                    block_boundary_edges.append((i, i + 1))
                for i in range(X, Y - 1):
                    block_boundary_edges.append((i, i + 1))
                for i in range(Y, Z - 1):
                    block_boundary_edges.append((i, i))
                for i in range(Z, 8):
                    block_boundary_edges.append((i, i - 1))

                block = generic_block.Block(k, number_vertices, edges, number_boundary_vertices, boundary_edges,
                                            block_boundary_edges)
                result = block.block_divergence(augmentation_vertex)

                print(f'{k}; Type 2.{X}{Y}{Z}; {len(block.fillings)}; {len(block.boundary_constraints)}; {result}; '
                      + f'{number_vertices}; {edges}; {number_boundary_vertices}; {boundary_edges}; {block_boundary_edges}')