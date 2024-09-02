# Block coupling and rapidly mixing k-heights #

This repository contains the supplementary code for the research article *Block coupling and rapidly mixing k-heights* (Felsner, Heldt, Roch, Winkler, 2024); we refer to this article for the theoretical background.

## Requirements

All scripts can be executed using Python 3. The exact version that we used for obtaining our data is Python 3.6.8 . The scripts use only the built-in-modules `math`, `random` and `itertools`, so no additional packages need to be installed.

## Code Files & Execution

* rectangular_grid.py: Computes the block divergence of (4 x 4)-blocks in toroidal hexagonal grid graphs; see Section 4.1 in the article. Run by:

        python3 rectangular_grid.py

* hexagonal_grid.py: Computes the block divergence of 6-blocks in toroidal hexagonal grid graphs; see Section 4.2 in the article. Run by:

        python3 hexagonal_grid.py

* three_regular_graphs.py: Considers the blocks in three regular planar graphs as specified in Section 4.3 in the article. Bounds the block divergence in each of the cases described there.

        python3 three_regular_graphs.py

* generic_block.py: Provides a generic framework for computing the block divergence for a single block with a certain boundary. Can only be used for small blocks; otherwise, the runtime explodes. This file is only used by three_regular_graphs.py and is not intended to be called directly.

For a detailed explanation of the performed computations, see the comments in the code.

## License

This project is licensed under the terms of the CC-BY 4.0 license; see `LICENSE.md`.