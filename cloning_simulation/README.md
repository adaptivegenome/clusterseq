cloning_simulation
==================

Functionality
-------
Simulates sample of barcoded cells.

We start with 300,000 cells barcoded with one of 14000 unique identifiers. The cells are evenly divided between the identifiers.

After this, we repeatedly perform the following:

* Grow the cell population. The number of cells of each population is increased to reflect exponential growth of the population. We use rate consistent with a doubling every 19 hours- this yields a ~13x growth across 3 days.
* After 3 days, randomly select 300,000 cells from the population (should be around 4,000,000 after 3 days of growth).
* Report the number of unique barcodes present in the 300,000 selected cells.

The above is simulated for 25 cycles (or 75 days).

How to run
------
Simply run the included simulation file with python:
```> python simulation.py```

Several configuration parameters are available at the top of the file, including the number of cells, barcodes, cell doubling frequency, etc. By default, the script will run 200 simulations, and write the number of remaining barcodes out to a file called 'histogram.csv'.

Results
------
The simulation was run for several values of the num_cells_to_keep_passage parameter in order to estimate how many cells must be kept through each passage in order to ensure that all barcodes were present in the sample after 75 days. Results for 300k, 450k, 600k, 750k, 1000k and 2000k are in the results/ directory.