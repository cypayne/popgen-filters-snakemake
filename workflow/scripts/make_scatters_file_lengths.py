#!/usr/bin/python3

## usage: ./make_scatters_file_lengths.py chr_lengths_sorted.txt first_N_scaf_groups max_scatter_length
##  first_N_scaf_groups means that the first N scafs will be 
##  turned into their own individual scaffold group
## ex. ./make_scatters_file_lengths.py ~/refs/Omykiss_Arlee_genome/Arlee-genome.chr-lengths.txt 32 5000000 
import sys

# set global variables
# tab-delimited file with scaffold lengths, sorted in descending order 
input_file = sys.argv[1]
# the first N scafs will get their own groups
first_N_scaf_groups = int(sys.argv[2]) 
# set scatter length 
max_scatter_length = int(sys.argv[3]) 
# set max size for a given scaffold group
scaf_size_cutoff = 5000000  
# header line for output file
print("id\tscatter_idx\tchrom\tstart\tend\tscatter_length")

# initialize counter variables 
scaf_counter = 1
size_sum = 0
scaf_grp_num = first_N_scaf_groups

with open(input_file) as scaf_size_file:
  for line in scaf_size_file:
    current_start = 1  # Start position for the first segment
    scatter_idx = 1  # Index for scatter segments

    values = line.strip().split() 
    scaf_name = values[0]
    scaf_size = int(values[1])
    if scaf_counter <= first_N_scaf_groups:
      # if one of the first N scaffolds, scaf_id should be the
      # scaffold name
      scaf_id = scaf_name
      scaf_counter+=1
      size_sum = scaf_size
    else:
      # if we've passed the first N scaffolds, then start creating
      # multi-scaffold groups
      # if adding this scaffold will push the group size over
      # the size limit, start a new scaffold 
      if size_sum + scaf_size > scaf_size_cutoff:
        # scaf_id will be a group scaffold name
        scaf_grp_num += 1
        group_name = "scaff_grp_"+str(scaf_grp_num)
        scaf_id = group_name
        # reset size sum and add +1 to the scaffold counter
        size_sum = 0
        scaf_counter+=1
      # otherwise add the scaffold length to the size sum
      else:
        size_sum = size_sum + scaf_size

    # Continue outputting segments until the scaffold length is reached
    while current_start <= size_sum:
      current_end = min(current_start + max_scatter_length - 1, size_sum)
      scatter_length = (current_end - current_start + 1)
      print(f"{scaf_id}\tscat_{scatter_idx:04d}\t{scaf_name}\t{current_start}\t{current_end}\t{scatter_length}")
      scatter_idx += 1
      current_start = current_end + 1

