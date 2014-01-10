#!/usr/bin/python
# Author: Lee Baker
# lee@leecbaker.com
#
import math
import random
import array

total_barcodes = 14000
num_cells_to_keep_passage = 1000000
print_messages = False

# all times in hours
cell_doubling_time = 19
passage_interval = 3 * 24
experiment_duration = 24 * 75

#######################

exponential_growth_factor = math.log(2) / cell_doubling_time

# initialize 300k cells, evenly divided between the 14000 barcodes
def initialize_cells():
	global cells
	cells = {}
	for cell in range(num_cells_to_keep_passage):
		barcode = cell * total_barcodes / num_cells_to_keep_passage
		if barcode in cells:
			cells[barcode] += 1
		else:
			cells[barcode] = 1

# grow each population of cells by barcode, independently
def grow_for_hours(hours):
	global cells
	starting_cell_count = sum(cells.values())
	for barcode in cells.keys():
		count = cells[barcode]
		cells[barcode] = int(round(count * math.exp(exponential_growth_factor * hours)))

	final_cell_count = sum(cells.values())
	if print_messages:
		print(" grew for %d hours: from %d to %d cells, with %d barcodes" % (passage_interval, starting_cell_count, final_cell_count, len(cells)))

# pick out 300k random cells to keep, discard the rest
# cells are randomly picked regardless of barcode
def passage():
	global cells
	starting_cell_count = sum(cells.values())
	cell_list = array.array('i')
	for barcode in cells.keys():
		count = cells[barcode]
		for i in range(count):
			cell_list.append(barcode)

	cell_list = random.sample(cell_list, num_cells_to_keep_passage)

	cells = {}
	for c in cell_list:
		if not c in cells:
			cells[c] = 1
		else:
			cells[c] += 1

	if print_messages:
		print(" passaging reduced %d cells to %d cells with %d unique barcodes" % (starting_cell_count, sum(cells.values()), len(cells)))

def run_simulation():
	global cells
	initialize_cells()
	cycles = experiment_duration / passage_interval

	for cycle in range(cycles):
		if print_messages:
			print("Cycle %3d/%d:" % (cycle+1, cycles))
		grow_for_hours(passage_interval)
		passage()
	if print_messages:
		print("Done. Executed for %.1f days, or %.1f population doublings" % (experiment_duration / 24., experiment_duration / cell_doubling_time))
	return len(cells)

def gen_histogram():
	iterations = 10
	barcode_counts = []
	for i in range(iterations):
		barcode_counts.append(run_simulation())

	f = open('histogram.csv','w')
	f.write("\n".join([str(ct) for ct in barcode_counts]))
	f.close()

gen_histogram()
	


