#
import os
import math
import pandas as pd

# pull out scatters
#scatter_groups = pd.read_table(config["scatters"]).set_index("id", drop=False)
#unique_scatters = list(scatter_groups.scatter_idx.unique())
#unique_scatters_table=scatter_groups[['id', 'scatter_idx']].drop_duplicates()
# create matching id-scatter pairs
#unique_scatter_id_matches = list(zip(unique_scatters_table['id'], unique_scatters_table['scatter_idx']))




run_ngsParalog = config['ngsParalog']
run_PCAone = config['PCAone']

if run_ngsParalog or run_PCAone:
    ## count how many scatters there will be
    # Read the input file into a pandas DataFrame
    chrom_lengths = pd.read_table(config['chrom_lengths_list'], header=None, names=['chrom', 'length'])
    # Calculate the number of scatter files (intervals) for each chromosome
    chrom_lengths['num_scatters'] = chrom_lengths['length'].apply(lambda x: math.ceil(x / 5000000))
    # Create a list of (chromosome, scatter_name)
    chromosome_list = []
    scatter_name_list = []
    for _, row in chrom_lengths.iterrows():
        chromosome = row['chrom']
        num_scatters = row['num_scatters']
        # Create scatter names for each chromosome 
        for scatter_number in range(1, num_scatters + 1):
            scatter_name = f"scatter{scatter_number:02d}"
            chromosome_list.append(chromosome)
            scatter_name_list.append(scatter_name)
    # Now zip the two lists (chromosome names and scatter names)
    chrom_scatter_table = pd.DataFrame(list(zip(chromosome_list, scatter_name_list)), columns=['chrom', 'scatter_name'])

    unique_chroms = list(set(chromosome_list))
    chrom_constraint = "|".join(unique_chroms)
else:
    chrom_constraint = "NoChroms"



#if run_ngsParalog: 
#    chroms = pd.read_table(config['chrom_lengths_list'],header=None)
#    chroms = pd.read_table(config['chrom_list'],header=None)
#    unique_chroms = chroms[0].unique().tolist()
    # make constraint 
#    chrom_constraint = "|".join(unique_chroms)
#else:
#    chrom_constraint = "NoChroms"

wildcard_constraints:
    chrom=chrom_constraint,
    #scatter="|".join(unique_scats),



#####

# create scatter bed files
def create_scatter_bed_files(input_file):
    bed_filename_list = []
    with open(input_file, 'r') as infile:
        for line in infile:
            # Split line into chromosome and length
            chrom, length = line.strip().split()
            length = int(length)  # Convert length to integer
            # Create intervals for the current chromosome
            start = 1
            scatter_number = 1
            while start < length:
                # Ensure the end doesn't exceed the chromosome length
                end = min(start + 5000000, length)                  
                scatter_name = f"scatter{scatter_number:02d}"
                # Create output filename based on chromosome and scatter number
                output_filename = f"results/scatter_beds/{chrom}_{scatter_name}.bed"
                bed_filename_list.append(output_filename)
                # Write to the separate BED file for this interval
                with open(output_filename, 'w') as outfile:
                    outfile.write(f"{chrom}\t{start}\t{end}\n")
                    # Increment start and scatter number for the next interval
                start = end
                scatter_number += 1
    return bed_filename_list
