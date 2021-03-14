# Optimal Coloring for Cell ID Barcodes (NeuroPAL)

This code was authored by Erdem Varol (Liam Paninski Lab) and Eviatar Yemini (Oliver Hobert Lab) at Columbia University for the NeuroPAL publication,
"NeuroPAL: A Multicolor Atlas for Whole-Brain Neuronal Identification in C. elegans". The publication is available here:

https://www.cell.com/cell/fulltext/S0092-8674(20)31682-2

The algorithmic code herein chooses approximately optimal multicolor solutions for cell identification. Running the algorithm multiple times generates a variety of reporter-fluorophore combinations to test in vivo. The algorithm requires a list of reporters with cell-specific expression (e.g., from WormBase for worms, FlyBase for flies, ZFIN for zebrafish, or MGI for mice), the number of colors available (e.g., 3 colors = RGB), a cell-adjacency matrix specifying which cells must be distinguishable from each other (e.g., neighboring cells), the margin for discriminating colors (e.g., 1 means we can only discriminate whether a color is present or absent whereas 1/3 permits discriminating bright, medium, weak, or no color expression), and a target for reporter sparsity (e.g., are we limited to ~3 transgenic reporters or can we use ~40 as was done in NeuroPAL). Each time the algorithm is run, it chooses a different approximately optimal set of reporter-color assignments. These approximately optimal solutions can then be applied in vivo to implement color barcodes for cell-specific identification.

## Getting Started

The code runs in MATLAB. Please clone the code to a local repository. Everything you need is documented in the file "optimal_color_parameters.m". You can run the NeuroPAL example code by calling the script "find_optimal_colors.m". Thereafter, you should edit "optimal_color_parameters.m" for your specific needs and call "find_optimal_colors.m". Adjust the algorithmic iterations to ensure convergence.

### Prerequisites

This code requires MATLAB. If you receive an error regarding a missing function, please ensure that your MATLAB is up to date and has the appropriate toolboxes.

### Installing

No installation is necessary. Simply clone the repository to your computer.

## Running the code

1. Edit "optimal_color_parameters.m" for your project.
2. Then run "find_optimal_colors.m".

### List of script files

1. optimal_color_parameters.m = edit this file for your project
2. find_optimal_colors.m = run this script to find approximately optimal coloring solutions

### List of function files
3. optimal_color_solver.m = converge on an approximately optimal solution
4. getReporterInfo.m = get reporter info from the NeuroPAL or Brain Atlas tables
5. computeProbAdjacency.m = compute a probabilistic adjacency matrix
6. plotCells.m = plot the cell coloring (e.g., show an approximately optimal worm solution)

### List of data files
7. data.mat = data for running NeuroPAL examples

## Contributing

Please contact the authors to discuss contributing to the work.

## Authors

* **Erdem Varol**
* **Eviatar Yemini**

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details.
