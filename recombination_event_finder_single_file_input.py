"""
A script used to identify recombination breakpoints using parental variant markers of skim-seq sequence data derived from a single inbred line.
The script identifies gross transition points within the data to identify both homozygous and heterozygous regions accurately.
Transitions in data are identified and scored, with the sharpest transitions scoring the highest. 
A logic tree then filters these candidate recombination sites, filtering out lower-scoring duplicates and impossible recombinations based on surrounding data.
If a recombiantion cannot be found between two high-confident regions, a second search ensues that finds the most favorable transition point
using a large sliding window, its size based on the proximity of more confident recombination breakpoints.
It is highly configurable, and works with high accuracy on skim-seq data only 0.028x in coverage.
"""
# Input file example containing the chromosome, position, and if it was the Reference or Alternative allele.
#chr1A_TA299   323292  P2
#chr1A_TA299   323424  P2
#chr1A_TA299   323528  P2
#chr1A_TA299   323867  P1
#...    ... ...

# Output file 1 example of the xy data, containing the chromosome, marker index, and roaming score (useful for plotting the function)
#chr1A_TA299     0       1
#chr1A_TA299     1       2
#chr1A_TA299     2       3
#chr1A_TA299     3       4
#...    ... ...

# Output file 2 example of the detected recombinations containing:
# Chromosome, left marker pos, right marker pos, average pos, transition type, marker index, and breakpoint score (out of 2)
#chr1A_TA299     529372358       529428788       529400573       P1/P2   109847  1.963
#chr2A_TA299     54816056        54820632        54818344        P2/P1   4710    1.927
#chr3A_TA299     42509727        42527879        42518803        P1/P2   5213    1.968
#chr4A_TA299     60813562        60824007        60818784        P2/P1   6995    1.939
#...    ... ...

####################################################################################################################################################
from csv import reader
import csv
import os
from itertools import islice
import numpy as np
import sys

#Read input file from first argument
single_input_file = str(sys.argv[1])

#For chromosome detection purposes, input the exact text that the first chromosome is represented as. Example: "chr1A_TA299"
current_chromosome = "chr1A_TA299"

# The whole size of the sliding window IN BP (physical length, containing both sides). The function within the sliding window will be utilized and used to find the
# Line of best fit for that side. Must be even. Can also be thought of as the smallest recombination site that can reasonabily form.
# The smallest recombination region manually seen at F8 generation is 3119601 bp (so the absolute largest the window should be not much larger than 6000000,
# although it would still likely detect it)
# Recommended: 3000000-8000000
sliding_window_size = 6500000

# Physical bp distance from the first or last marker where the algorithm will begin to look for bends.
# Shoude be equal to or less than half of the sliding window size.
# Keep in mind that if the marker density is too low it will be filtered out anyways, so lean towards the lower end.
telomere_buffer_width = 100000

# Specify the Slope Tolerance for hom and het regions. If the slope of the line of best fit is within the tolerance from 1, -1, or 0, it will be considered a P1, P2,
# or Het region respectively. In an island of matches around a recombination breakpoint, the lines with the closest fit will be considered the correct candidate.
# Err on the higher side-- the algorithm will select the best slope pairing in an island of matches. Must be below 0.5.
slope_tolerance_hom = 0.2
slope_tolerance_het = 0.4

# Score tolerance for something to be considered an anchor site. If something is considered an anchor site, its region type MUST be resolved at both ends.
# A score of 0 would be a perfect anchor site. A tolerance of 0.1 would mean there can be only 0.1 difference in slope between the two sides of the window.
score_tolerance_anchor_hom = 0.08
score_tolerance_anchor_het = 0.1

# If a sliding window only contains alleles of one parent, create an anchor site on that side.
spawn_anchors_next_to_perfect_sites = False

# Correlation Coefficient minimum of homozygous regions for it to be considered P1 or P2. A perfect correlation to a linear function would be +1.
corr_coef_min_hom = 0.6

# Correlation Coefficient minimum of heterozygous regions for it to be considered a Het region. This value should be equal to or less than the homozygous,
# because heterozygous regions are significantly more rough than homozygous regions.
corr_coef_min_het = 0.2

# When desperately searching for a site, if the calculated window size is smaller than this amount of indexes, make it this amount of markers wide.
minimum_desperate_window_size = 20
#The output file is in the following, headerless, format:

# Marker density lower limit. If a window's marker density is below this percentage of the calculated marker density (per bp) of the chromosome,
# the whole jumping window will be passed over.
# Example: if 0.1, a window with a marker density per bp below a tenth of the average of the whole
# chromosome will be passed over and ignored on the first pass.
# The absolute minimum markers required for a window is 3: 2 can mean any errant marker can cause a perfect score on that side
marker_density_lower_limit = 0.25

# Boolean to determine if the marker density filter is passed when calculated all sites or just anchor site candidates.
marker_density_check_on_all_windows = True

# Boolean to output information concerning candidate recombination sites.
verbose_output = False
####################################################################################################################################################

half_window_size = sliding_window_size / 2

def get_window(seq, n):
    # Returns a sliding window of width n over seq data
    iterable = iter(seq)
    result = tuple(islice(iterable, n))
    if len(result) == n:
        yield result
    for element in iterable:
        result = result[1:] + (element,)
        yield result

def find_events(window_list_allele,window_list_pos):
    allele_index_array = []
    allele_roam_y_array = []
    allele_roam_pos = 0

    for n in range(len(window_list_allele)-1):

        first_allele = window_list_allele[n]
        second_allele = window_list_allele[n+1]
        if first_allele != second_allele or (n == len(window_list_allele)-2) or (n == 0):
            allele_index_array.append(n)

        if (first_allele == "P1"):
            allele_roam_pos = allele_roam_pos + 1
            allele_roam_y_array.append(allele_roam_pos)

        elif (first_allele == "P2"):
            allele_roam_pos = allele_roam_pos - 1
            allele_roam_y_array.append(allele_roam_pos)

        elif (first_allele == "P0"):
            allele_roam_y_array.append(allele_roam_pos)
            # Dummy telomere buffer
        else:
            print("ERROR: Not P1 or P2")

        if n == (len(window_list_allele)-2):
            if (second_allele == "P1"):
                allele_roam_pos = allele_roam_pos + 1
                allele_roam_y_array.append(allele_roam_pos)
                #print("P1 found in first half")
            elif (second_allele == "P2"):
                allele_roam_pos = allele_roam_pos - 1
                allele_roam_y_array.append(allele_roam_pos)
                #print("P2 found in first half")
            elif (second_allele == "P0"):
                allele_roam_y_array.append(allele_roam_pos)
                # Dummy telomere buffer
            else:
                print("ERROR: Not P1 or P2")

    
    # If there is not a site to investigate within the range of the sliding window, make a new site just to create an anchor
    allele_index_array_with_anchors = []
    for i in range(len(allele_index_array)-1):
        allele_index_array_with_anchors.append(allele_index_array[i])
        if (window_list_pos[allele_index_array[i+1]] - window_list_pos[allele_index_array[i]]) > half_window_size:
            position_to_add = int(np.round((window_list_pos[allele_index_array[i+1]] + window_list_pos[allele_index_array[i]]) / 2,0))
            closest_index_of_pos_to_add = (min(range(len(window_list_pos)), key=lambda q: abs(window_list_pos[q]-position_to_add)))
            if closest_index_of_pos_to_add != allele_index_array[i] and closest_index_of_pos_to_add != allele_index_array[i+1]:
                allele_index_array_with_anchors.append(closest_index_of_pos_to_add)
    allele_index_array_with_anchors.append(allele_index_array[-1])

    # This array will output the index of the marker that comes just before a P1/P2 or P2/P1 split
    return allele_index_array_with_anchors,allele_roam_y_array

def find_points_in_window(pos_index,x_pos_list,y_pos_list,x_basic_index_list):
    current_physical_pos_down = x_pos_list[pos_index]
    current_physical_pos_up = x_pos_list[pos_index+1]
    current_physical_pos_ave = int(np.round((current_physical_pos_down + current_physical_pos_up) / 2,0))
    downstream_physical_pos = current_physical_pos_ave - (half_window_size)
    upstream_physical_pos = current_physical_pos_ave + (half_window_size)

    x_points_out_Lwindow = []
    y_points_out_Lwindow = []

    x_points_out_Rwindow = []
    y_points_out_Rwindow = []

    for pos in reversed(range(0,pos_index+1)):
        if x_pos_list[pos] > downstream_physical_pos:
            #print("not to downstream pos yet... at " + str(x_pos_list[pos]))
            x_points_out_Lwindow.insert(0,x_basic_index_list[pos])
            y_points_out_Lwindow.insert(0,y_pos_list[pos])
        else:
            x_points_out_Lwindow.insert(0,x_basic_index_list[pos])
            y_points_out_Lwindow.insert(0,y_pos_list[pos])
            break

    for pos in range(pos_index,len(x_pos_list)):
        if x_pos_list[pos] < upstream_physical_pos:
            #print("not to uptream pos yet... at " + str(x_pos_list[pos]))
            x_points_out_Rwindow.append(x_basic_index_list[pos])
            y_points_out_Rwindow.append(y_pos_list[pos])
        else:
            x_points_out_Rwindow.append(x_basic_index_list[pos])
            y_points_out_Rwindow.append(y_pos_list[pos])
            #print("past upstream pos")
            break

    return x_points_out_Lwindow,y_points_out_Lwindow,x_points_out_Rwindow,y_points_out_Rwindow

# Return the slope of the line of best fit
def find_slope_of_LOBF(points_in_x,points_in_y):

    xbar = sum(points_in_x)/len(points_in_x)
    ybar = sum(points_in_y)/len(points_in_y)
    n = len(points_in_x) # or len(Y)

    numer = sum([xi*yi for xi,yi in zip(points_in_x, points_in_y)]) - n * xbar * ybar
    denum = sum([xi**2 for xi in points_in_x]) - n * xbar**2

    if denum != 0:
        b = numer / denum
    else:
        b = 99999

    return b

# Find correlation coefficient
def find_corrcoef(points_in_x,points_in_y):
    return np.corrcoef(points_in_x, points_in_y)[0,1]

def cull_from_list(input_list,list_of_indexes):
    indexes = sorted(list_of_indexes, reverse=True)
    for id in indexes :
        if id < len(input_list):
            input_list.pop(id)

# Finds furthest right anchor point of the same type as left anchor point input.
# Finds furthest left anchor point oof the same type as right anchor point input.
# Grab points on either side of a window NOT based on physical distance-- based on number of alleles (for macro comprehension)
# Start on left new anchor point, work way down alleles grabbing slopes
# The highest scoring breakpoint that matches the expected region layout is considered the winner. Its stats are returned and appended to the island list.
def find_score_desperate(L_anchor_regionT, R_anchor_regionT, L_anchor_index, R_anchor_index, index_archive, regionT_archive, score_archive, tinylist_y, tinylist_basic_x_index, tinylist_bend_indexes, tinylist_bend_indexes_smol, pos_tiny_actual_list):
    # print("looking to find a point with index between:")
    # print(tinylist_bend_indexes_smol[L_anchor_index])
    # print(tinylist_bend_indexes_smol[R_anchor_index])
    if L_anchor_regionT == "P1/P1":
        looking_for_bend_down = True

    elif L_anchor_regionT == "P2/P2":
        looking_for_bend_down = False

    elif L_anchor_regionT == "HET/HET" and R_anchor_regionT == "P1/P1":

        looking_for_bend_down = False
    elif L_anchor_regionT == "HET/HET" and R_anchor_regionT == "P2/P2":

        looking_for_bend_down = True
    else:
        print("ERROR: ANCHOR REGION NOT P1/P1, P2/P2, or HET/HET")
    window_width_units = int(tinylist_bend_indexes_smol[R_anchor_index] - tinylist_bend_indexes_smol[L_anchor_index])
    if window_width_units < minimum_desperate_window_size:
        window_width_units = minimum_desperate_window_size
    window_width_units_L = window_width_units
    window_width_units_R = window_width_units
    tinylist_bend_indexes_smol_desp = []
    all_scores_list = []
    all_indexes_list = []

    for q in tinylist_bend_indexes_smol:
        #real_ind = tinylist_bend_indexes[q]
        if q >= tinylist_bend_indexes_smol[L_anchor_index] and q <= tinylist_bend_indexes_smol[R_anchor_index]:
            tinylist_bend_indexes_smol_desp.append(q)

    for pos_index in tinylist_bend_indexes_smol_desp:
        window_width_units_L = window_width_units
        window_width_units_R = window_width_units

        if pos_index - window_width_units_L < 0:
            window_width_units_L = pos_index
        else:
            window_width_units_L = window_width_units

        if pos_index + window_width_units_R+1 > len(tinylist_basic_x_index):
            window_width_units_R = len(tinylist_basic_x_index) - pos_index - 1
        else:
            window_width_units_R = window_width_units

        if window_width_units_L <= 20 or window_width_units_R <= 20:
            #print("too small")
            pass

        else:
            #print("begin search")
            x_points_out_Lwindow = []
            y_points_out_Lwindow = []
            x_points_out_Rwindow = []
            y_points_out_Rwindow = []

            # print("begin l window")
            for ind in reversed(range(pos_index-window_width_units_L,pos_index+1)):
                # print("l window")
                x_points_out_Lwindow.insert(0,tinylist_basic_x_index[ind])
                y_points_out_Lwindow.insert(0,tinylist_y[ind])

            # print("begin r window")
            for ind in range(pos_index+1,pos_index+window_width_units_R+1):
                # print("r window")
                x_points_out_Rwindow.append(tinylist_basic_x_index[ind])
                y_points_out_Rwindow.append(tinylist_y[ind])

            L_window_slope = find_slope_of_LOBF(x_points_out_Lwindow,y_points_out_Lwindow)
            R_window_slope = find_slope_of_LOBF(x_points_out_Rwindow,y_points_out_Rwindow)

            if looking_for_bend_down:
                score_for_window = L_window_slope - R_window_slope
            else:
                score_for_window = R_window_slope - L_window_slope

            all_scores_list.append(score_for_window)
            all_indexes_list.append(pos_index)

    max_score = max(all_scores_list)
    max_score_list = []
    max_score_index_list = []

    for s in range(len(all_scores_list)):
        if all_scores_list[s] == max_score:
            max_score_list.append(s)
            max_score_index_list.append(all_indexes_list[s])

    if len(max_score_list) > 1:
        # print("WARNING: Somehow multiple high scores when desperately searching for a site. Choosing middle one...")
        middleIndex = int((len(max_score_list) - 1)/2)
        winning_score = max_score_list[middleIndex]
        winning_index = max_score_index_list[middleIndex]
    else:
        winning_score = max_score_list[0]
        winning_index = max_score_index_list[0]

    true_winning_index = tinylist_bend_indexes_smol.index(winning_index)

    return winning_score, true_winning_index

def determine_event_positions(input_file,output_file_1,output_file_2):
    #print(input_file)
    #print(output_file)
    with open(input_file) as csv_file:
        csv_reader = csv.reader(csv_file, delimiter='\t')
        length_of_file = sum(1 for _ in csv_reader)
        csv_file.seek(0)
        csv_reader = csv.reader(csv_file, delimiter='\t')
        s_chrom = []
        s_chrom_temp = []
        s_pos = []
        s_pos_temp = []
        s_allele = []
        s_allele_temp = []
        #If there is a header in the input file, uncomment this
        #header = next(csv_reader)
        current_chrom = current_chromosome
        line_counter = 0
        for lines in csv_reader:
            line_counter = line_counter + 1
            if lines[0] != "chrUn_TA299":
                if (current_chrom != lines[0] or line_counter == length_of_file):
                    current_chrom = lines[0]
                    s_chrom.append(s_chrom_temp)
                    s_pos.append(s_pos_temp)
                    s_allele.append(s_allele_temp)
                    s_chrom_temp = []
                    s_pos_temp = []
                    s_allele_temp = []
                s_chrom_temp.append(lines[0])
                s_pos_temp.append(int(float(lines[1])))
                s_allele_temp.append(lines[2])
            else:
                s_chrom.append(s_chrom_temp)
                s_pos.append(s_pos_temp)
                s_allele.append(s_allele_temp)
                s_chrom_temp = []
                s_pos_temp = []
                s_allele_temp = []
                break

    sequence_counter = 0

    # Creates a file containing data used to represent the function for the chromosomes. Used to easily plot graphs.
    with open(output_file_1, 'w') as csvoutput1:

        writer = csv.writer(csvoutput1, lineterminator='\n', delimiter='\t')
        all=[]
        indexes_to_investigate_gigalist = []
        pos_out_gigalist_y = []
        #Iterates through the list of lists; for each iteration, a whole list of alleles (P1s or P2s) are provided to get_window, the list contained in a chromosome.
        chrom_count = -1
        for allele_tiny_list in s_allele:
            # Iterates through the list of windows
            chrom_count = chrom_count + 1
            pos_tiny_list = s_pos[chrom_count]
            allele_tiny_actual_list = list(allele_tiny_list)
            pos_tiny_actual_list = list(pos_tiny_list)
            index_list_to_investigate,allele_roam_y_value_array = find_events(allele_tiny_actual_list,pos_tiny_actual_list)
            indexes_to_investigate_gigalist.append(index_list_to_investigate)
            pos_out_gigalist_y.append(allele_roam_y_value_array)
            basic_index = -1

            for i in range(len(allele_roam_y_value_array)):
                basic_index = basic_index + 1
                row=[]
                row.append(s_chrom[chrom_count][0])
                row.append(basic_index)
                row.append(allele_roam_y_value_array[i])
                all.append(row)
        writer.writerows(all)

    # Create actual breakpoint data file
    with open(output_file_2, 'w') as csvoutput2:

        writer2 = csv.writer(csvoutput2, lineterminator='\n', delimiter='\t')
        all=[]

        current_window_pos = half_window_size + s_pos[0][0]
        for c in range(len(indexes_to_investigate_gigalist)):
            island_pos_list = []
            island_pos_list_2 = []
            island_region_L = []
            island_region_R = []
            island_region_text = []
            island_score_list = [
            tinylist_bend_indexes = indexes_to_investigate_gigalist[c]
            tinylist_y = pos_out_gigalist_y[c]
            tinylist_pos = s_pos[c]
            tinylist_basic_x_index = list(range(len(tinylist_y)))

            # Here we will find the average marker density, so that we can make sure when checking our windows the marker density
            # does not drop too much.
            tinylist_chrom_width_bp = tinylist_pos[-1] - tinylist_pos[0]
            markers_per_bp = len(tinylist_basic_x_index) / tinylist_chrom_width_bp
            marker_count_minimum  = marker_density_lower_limit * markers_per_bp * half_window_size
            if marker_count_minimum < 3:
                marker_count_minimum = 3
            starting_physical_pos = int(telomere_buffer_width + tinylist_pos[0])
            ending_physical_pos = int(tinylist_pos[-1] - telomere_buffer_width)

            tinylist_bend_indexes_smol = []

            for q in range(len(tinylist_bend_indexes)):
                real_pos = tinylist_pos[tinylist_bend_indexes[q]]
                if real_pos >= starting_physical_pos and real_pos <= ending_physical_pos:
                    tinylist_bend_indexes_smol.append(tinylist_bend_indexes[q])

            island_index_list = []
            island_score_list = []
            island_Pscore_if_het_list = []
            island_pos_list = []
            island_pos_list_2 = []
            island_pos_list_ave = []
            island_region_text = []
            island_region_L = []
            island_region_R = []

            for p in range(len(tinylist_bend_indexes_smol)):
                x_points_out_Lwindow,y_points_out_Lwindow,x_points_out_Rwindow,y_points_out_Rwindow = find_points_in_window(tinylist_bend_indexes_smol[p],tinylist_pos,tinylist_y,tinylist_basic_x_index)

                L_window_slope = find_slope_of_LOBF(x_points_out_Lwindow,y_points_out_Lwindow)
                R_window_slope = find_slope_of_LOBF(x_points_out_Rwindow,y_points_out_Rwindow)
                force_none = False
                if L_window_slope == 99999 or R_window_slope == 99999:
                    force_none = True
                else:
                    corrcoef_L = abs(find_corrcoef(x_points_out_Lwindow,y_points_out_Lwindow))
                    corrcoef_R = abs(find_corrcoef(x_points_out_Rwindow,y_points_out_Rwindow))

                L_window_region = "NONE"
                R_window_region = "NONE"

                #The scores for the window are the difference between the slopes
                score_for_window_L = 1
                score_for_window_R = 1

                if force_none == False:
                    if L_window_slope >=  1 - slope_tolerance_hom:
                        # P1 Region
                        if corrcoef_L >= corr_coef_min_hom:
                            L_window_region = "P1"
                            #score_for_window_L = 1 - L_window_slope
                    elif L_window_slope <= -1 + slope_tolerance_hom:
                        # P2 Region
                        if corrcoef_L >= corr_coef_min_hom:
                            L_window_region = "P2"
                            #score_for_window_L = 1 + L_window_slope
                    elif (L_window_slope <= slope_tolerance_het) and (L_window_slope >= (-1 * slope_tolerance_het)):
                        # Het Region
                        if corrcoef_L >= corr_coef_min_het:
                            L_window_region = "HET"
                            #score_for_window_L = L_window_slope

                    if R_window_slope >=  1 - slope_tolerance_hom:
                        # P1 Region
                        if corrcoef_R >= corr_coef_min_hom:
                            R_window_region = "P1"
                            #score_for_window_R = 1 - R_window_slope
                    elif R_window_slope <= -1 + slope_tolerance_hom:
                        # P2 Region
                        if corrcoef_R >= corr_coef_min_hom:
                            R_window_region = "P2"
                            #score_for_window_R = 1 + R_window_slope
                    elif (R_window_slope <= slope_tolerance_het) and (R_window_slope >= (-1 * slope_tolerance_het)):
                        # Het Region
                        if corrcoef_R >= corr_coef_min_het:
                            R_window_region = "HET"
                            #score_for_window_R = R_window_slope

                    Pslope_if_het = 99999
                    score_for_window_average = abs(L_window_slope - R_window_slope)
                    if L_window_region == "HET" and R_window_region == "HET":
                        Pslope_if_het = 99999
                    # Here, heterozgyous scores are calculated. Heterozgyous scores are scored by counting the slope of the
                    # non-het side twice, giving it twice the importance as the het side.
                    # Ideally, the best score should be around 2 again, but it is possible to be over that.
                    elif L_window_region == "HET" and R_window_region != "NONE":
                        Pslope_if_het = abs(R_window_slope) + score_for_window_average
                    elif R_window_region == "HET" and L_window_region != "NONE":
                        Pslope_if_het = abs(L_window_slope) + score_for_window_average

                if (L_window_region == R_window_region) and L_window_region != "NONE" and R_window_region != "NONE":
                    if ((score_for_window_average <= score_tolerance_anchor_het) and L_window_region == "HET") or score_for_window_average <= score_tolerance_anchor_hom:
                        # See if the marker density is high enough to consider this an anchor site
                        if marker_count_minimum <= len(x_points_out_Lwindow) and marker_count_minimum <= len(x_points_out_Rwindow):
                            # Found anchor site
                            island_index_list.append(p)
                            # The score for the window is going to be very low-- but for an anchor site, thats actually great!
                            # Due to the best score being 2, we take 2 minus the score to get the actual confidence level for the window
                            island_score_list.append(2 - score_for_window_average)
                            island_Pscore_if_het_list.append(Pslope_if_het)
                            island_pos_list.append(tinylist_pos[tinylist_bend_indexes_smol[p]])
                            island_pos_list_2.append(tinylist_pos[tinylist_bend_indexes_smol[p]+1])
                            island_pos_list_ave.append(int(np.round((tinylist_pos[tinylist_bend_indexes_smol[p]] + tinylist_pos[tinylist_bend_indexes_smol[p]+1]) / 2,0)))
                            island_region_text.append(str(L_window_region) + "/" + str(R_window_region))
                            island_region_L.append(L_window_region)
                            island_region_R.append(R_window_region)

                elif (L_window_region != R_window_region) and L_window_region != "NONE" and R_window_region != "NONE":
                    # For anchor sites created from the sides, their scores and positions are made equal to their
                    # breakpoint candidate so they can all be filtered easily in a lump
                    if abs(L_window_slope) >= 1.0 - (score_tolerance_anchor_hom/10) and marker_count_minimum <= len(x_points_out_Lwindow) and marker_count_minimum <= len(x_points_out_Rwindow):
                        if spawn_anchors_next_to_perfect_sites:
                            # Consider this a sequene good enough to make sure there is an anchor site before this as a safeguard...
                            island_index_list.append(p)
                            island_score_list.append(score_for_window_average)
                            island_Pscore_if_het_list.append(Pslope_if_het)
                            island_pos_list.append(tinylist_pos[tinylist_bend_indexes_smol[p]])
                            island_pos_list_2.append(tinylist_pos[tinylist_bend_indexes_smol[p]+1])
                            island_pos_list_ave.append(int(np.round((tinylist_pos[tinylist_bend_indexes_smol[p]] + tinylist_pos[tinylist_bend_indexes_smol[p]+1]) / 2,0)))
                            island_region_text.append(str(L_window_region) + "/" + str(L_window_region))
                            island_region_L.append(L_window_region)
                            island_region_R.append(L_window_region)

                    if marker_density_check_on_all_windows == True:
                        # Need to check this site to make sure it is above the minimum marker density
                        if marker_count_minimum <= len(x_points_out_Lwindow) and marker_count_minimum <= len(x_points_out_Rwindow):
                            island_index_list.append(p)
                            island_score_list.append(score_for_window_average)
                            island_Pscore_if_het_list.append(Pslope_if_het)
                            island_pos_list.append(tinylist_pos[tinylist_bend_indexes_smol[p]])
                            island_pos_list_2.append(tinylist_pos[tinylist_bend_indexes_smol[p]+1])
                            island_pos_list_ave.append(int(np.round((tinylist_pos[tinylist_bend_indexes_smol[p]] + tinylist_pos[tinylist_bend_indexes_smol[p]+1]) / 2,0)))
                            island_region_text.append(str(L_window_region) + "/" + str(R_window_region))
                            island_region_L.append(L_window_region)
                            island_region_R.append(R_window_region)
                    else:
                        # Do not need to check this site to make sure it is above the minimum marker density
                        island_index_list.append(p)
                        island_score_list.append(score_for_window_average)
                        island_Pscore_if_het_list.append(Pslope_if_het)
                        island_pos_list.append(tinylist_pos[tinylist_bend_indexes_smol[p]])
                        island_pos_list_2.append(tinylist_pos[tinylist_bend_indexes_smol[p]+1])
                        island_pos_list_ave.append(int(np.round((tinylist_pos[tinylist_bend_indexes_smol[p]] + tinylist_pos[tinylist_bend_indexes_smol[p]+1]) / 2,0)))
                        island_region_text.append(str(L_window_region) + "/" + str(R_window_region))
                        island_region_L.append(L_window_region)
                        island_region_R.append(R_window_region)

                    if abs(R_window_slope) >= 1.0 - (score_tolerance_anchor_hom/10) and marker_count_minimum <= len(x_points_out_Lwindow) and marker_count_minimum <= len(x_points_out_Rwindow):
                        if spawn_anchors_next_to_perfect_sites:
                            # Consider this a sequene good enough to make sure there is an anchor site after this as a safeguard...
                            island_index_list.append(p)
                            island_score_list.append(score_for_window_average)
                            island_Pscore_if_het_list.append(Pslope_if_het)
                            island_pos_list.append(tinylist_pos[tinylist_bend_indexes_smol[p]])
                            island_pos_list_2.append(tinylist_pos[tinylist_bend_indexes_smol[p]+1])
                            island_pos_list_ave.append(int(np.round((tinylist_pos[tinylist_bend_indexes_smol[p]] + tinylist_pos[tinylist_bend_indexes_smol[p]+1]) / 2,0)))
                            island_region_text.append(str(R_window_region) + "/" + str(R_window_region))
                            island_region_L.append(R_window_region)
                            island_region_R.append(R_window_region)

            indexes_to_cull = []

            # Here regions and anchors with the same region type are compared.
            # IF there is a higher scoring region within half the window size away, remove the lower scoring region
            for center_index in range(len(island_index_list)):
                L_window_winner = True
                R_window_winner = True

                # If the center index involves a heterozgyous site, use that scoring method
                if island_region_L[center_index] == "HET" or island_region_L[center_index] == "HET":
                    if center_index == 0:
                        L_window_winner = True
                    else:
                        for L_search_index in reversed(range(0,center_index)):
                            # If the checked position is within range
                            if island_pos_list_ave[center_index] - island_pos_list_ave[L_search_index] < half_window_size:
                                if island_region_text[center_index] == island_region_text[L_search_index] and island_Pscore_if_het_list[center_index] <= island_Pscore_if_het_list[L_search_index]:
                                    # Looks like you're the loser on this end!
                                    L_window_winner = False
                                    break
                            else:break
                    if len(island_index_list) == center_index+1:
                        R_window_winner = True
                    else:
                        for R_search_index in range(center_index+1,len(island_index_list)):
                            if island_pos_list_ave[R_search_index] - island_pos_list_ave[center_index] < half_window_size:
                                if island_region_text[center_index] == island_region_text[R_search_index] and island_score_list[center_index] < island_score_list[R_search_index]:
                                    # Looks like you're the loser on this end!
                                    R_window_winner = False
                                    break
                            else: break
                else:
                    if center_index == 0:
                        L_window_winner = True
                    else:
                        for L_search_index in reversed(range(0,center_index)):
                            # If the checked position is within range
                            if island_pos_list_ave[center_index] - island_pos_list_ave[L_search_index] < half_window_size:
                                if island_region_text[center_index] == island_region_text[L_search_index] and island_score_list[center_index] <= island_score_list[L_search_index]:
                                    # Looks like you're the loser on this end!
                                    L_window_winner = False
                                    break
                            else:break

                    if len(island_index_list) == center_index+1:
                        R_window_winner = True
                    else:
                        for R_search_index in range(center_index+1,len(island_index_list)):
                            if island_pos_list_ave[R_search_index] - island_pos_list_ave[center_index] < half_window_size:
                                if island_region_text[center_index] == island_region_text[R_search_index] and island_score_list[center_index] < island_score_list[R_search_index]:
                                    # Looks like you're the loser on this end!
                                    R_window_winner = False
                                    break
                            else: break

                if L_window_winner == False or R_window_winner == False:
                    indexes_to_cull.append(center_index)

            cull_from_list(island_index_list,indexes_to_cull)
            cull_from_list(island_score_list,indexes_to_cull)
            cull_from_list(island_Pscore_if_het_list,indexes_to_cull)
            cull_from_list(island_pos_list,indexes_to_cull)
            cull_from_list(island_pos_list_2,indexes_to_cull)
            cull_from_list(island_pos_list_ave,indexes_to_cull)
            cull_from_list(island_region_text,indexes_to_cull)
            cull_from_list(island_region_L,indexes_to_cull)
            cull_from_list(island_region_R,indexes_to_cull)
            indexes_to_cull = []

            # For later, save a list of the original island indexes, scores, and region texts.
            # This is saved after the first step of minute filtering to remove anchors that were erroniously placed
            # from high quality half-sites.
            island_OG_index_archive = island_index_list
            island_OG_score_archive = island_score_list
            island_OG_regionT_archive = island_region_text

            # Cull out duplicate anchor sites
            for center_index in range(len(island_index_list)-1):
                if island_region_L[center_index] == island_region_R[center_index] == island_region_L[center_index+1] == island_region_R[center_index+1]:
                    indexes_to_cull.append(center_index)

            cull_from_list(island_index_list,indexes_to_cull)
            cull_from_list(island_score_list,indexes_to_cull)
            cull_from_list(island_Pscore_if_het_list,indexes_to_cull)
            cull_from_list(island_pos_list,indexes_to_cull)
            cull_from_list(island_pos_list_2,indexes_to_cull)
            cull_from_list(island_pos_list_ave,indexes_to_cull)
            cull_from_list(island_region_text,indexes_to_cull)
            cull_from_list(island_region_L,indexes_to_cull)
            cull_from_list(island_region_R,indexes_to_cull)
            indexes_to_cull = []

            # Check to see if the organization of regions make logical sense: if a region between two anchor sites has a region on the left or right
            # that cannot possibly be paired with anything else in the region, cull it
            seen_Ls = []
            seen_Rs = []
            seen_L_at = []
            seen_R_at = []
            list_of_anchor_indexes = []
            for s in range(len(island_region_L)):
                if island_region_L[s] == island_region_R[s]:
                    list_of_anchor_indexes.append(s)
            if len(list_of_anchor_indexes) < 2:
                pass
            else:
                for a in range(len(list_of_anchor_indexes)-1):
                    seen_Ls = []
                    seen_Rs = []
                    seen_L_at = []
                    seen_R_at = []
                    start_anchor_index = list_of_anchor_indexes[a]
                    end_anchor_index = list_of_anchor_indexes[a+1]
                    seen_Rs.append(island_region_R[start_anchor_index])
                    seen_R_at.append(start_anchor_index)
                    for s in range(start_anchor_index+1,end_anchor_index):
                        seen_Ls.append(island_region_L[s])
                        seen_Rs.append(island_region_R[s])
                        seen_L_at.append(s)
                        seen_R_at.append(s)
                    seen_Ls.append(island_region_L[end_anchor_index])
                    seen_L_at.append(end_anchor_index)

                    for l in range(len(seen_Ls)-1):
                        if seen_Ls[l] not in seen_Rs[:l+1]:
                            indexes_to_cull.append(seen_L_at[l])

                    for r in range(1,len(seen_Rs)):
                        if seen_Rs[r] not in seen_Ls[r:]:
                            indexes_to_cull.append(seen_R_at[r])

            cull_from_list(island_index_list,indexes_to_cull)
            cull_from_list(island_score_list,indexes_to_cull)
            cull_from_list(island_Pscore_if_het_list,indexes_to_cull)
            cull_from_list(island_pos_list,indexes_to_cull)
            cull_from_list(island_pos_list_2,indexes_to_cull)
            cull_from_list(island_pos_list_ave,indexes_to_cull)
            cull_from_list(island_region_text,indexes_to_cull)
            cull_from_list(island_region_L,indexes_to_cull)
            cull_from_list(island_region_R,indexes_to_cull)
            indexes_to_cull = []

            #Merge together anchor sites again after cleaning
            for center_index in range(len(island_index_list)-1):
                if island_region_L[center_index] == island_region_R[center_index] == island_region_L[center_index+1] == island_region_R[center_index+1]:
                    indexes_to_cull.append(center_index)

            cull_from_list(island_index_list,indexes_to_cull)
            cull_from_list(island_score_list,indexes_to_cull)
            cull_from_list(island_Pscore_if_het_list,indexes_to_cull)
            cull_from_list(island_pos_list,indexes_to_cull)
            cull_from_list(island_pos_list_2,indexes_to_cull)
            cull_from_list(island_pos_list_ave,indexes_to_cull)
            cull_from_list(island_region_text,indexes_to_cull)
            cull_from_list(island_region_L,indexes_to_cull)
            cull_from_list(island_region_R,indexes_to_cull)
            indexes_to_cull = []

            go_again = True
            if len(island_score_list) > 2:
                while go_again:
                    #go_again = False
                    go_again_for_duplicate_culling = []
                    go_again_for_noncontiguous_culling = []
                    last_anchor_site = None
                    list_of_anchor_indexes = []
                    for s in range(len(island_region_L)):
                        if island_region_L[s] == island_region_R[s]:
                            list_of_anchor_indexes.append(s)

                    if list_of_anchor_indexes[0] > 0:
                        for l in reversed(range(list_of_anchor_indexes[0])):
                            if island_region_R[l] != island_region_L[l+1]:
                                indexes_to_cull.append(l)
                                break

                    if list_of_anchor_indexes[-1] < len(island_index_list)-1:
                        for l in range(list_of_anchor_indexes[-1]+1,len(island_index_list)):
                            if island_region_R[l-1] != island_region_L[l]:
                                indexes_to_cull.append(l)
                                break

                    if len(list_of_anchor_indexes) < 2:
                        go_again_for_duplicate_culling.append(0)
                        go_again_for_noncontiguous_culling.append(0)
                        go_again == False
                    else:
                        # Outside of anchor work, just in case one doesn't spawn at the ends...
                        break_out_flag = False
                        for a in range(len(list_of_anchor_indexes)-1):
                            start_anchor_index = list_of_anchor_indexes[a]
                            end_anchor_index = list_of_anchor_indexes[a+1]

                            if end_anchor_index - start_anchor_index == 1:
                                # The anchors are right next to each other, but if the anchors are the same, just delete the first duplicate
                                if island_region_R[start_anchor_index] == island_region_L[end_anchor_index]:
                                    go_again_for_duplicate_culling.append(1)
                                    go_again_for_noncontiguous_culling.append(0)
                                    indexes_to_cull.append(start_anchor_index)
                                else:
                                    # This is not good, two anchors that don't match are right next to each other, with nothing in between!
                                    print("WARNING: Desperately search for a site for...")
                                    print(input_file)
                                    print(c+1)
                                    L_anchor_holder = island_region_R[start_anchor_index]
                                    R_anchor_holder = island_region_L[end_anchor_index]
                                    winning_score,winning_index = find_score_desperate(island_region_text[start_anchor_index], island_region_text[end_anchor_index], island_index_list[start_anchor_index], island_index_list[end_anchor_index], island_OG_index_archive, island_OG_regionT_archive, island_OG_score_archive, tinylist_y, tinylist_basic_x_index, tinylist_bend_indexes, tinylist_bend_indexes_smol,pos_tiny_actual_list)
                                    anchor_desp_pos_1 = tinylist_pos[tinylist_bend_indexes_smol[winning_index]]
                                    anchor_desp_pos_2 = tinylist_pos[tinylist_bend_indexes_smol[winning_index]+1]
                                    island_index_list.insert(end_anchor_index,winning_index)
                                    island_score_list.insert(end_anchor_index,winning_score)

                                    # Yes the score won't be accurate if its a het site, max should be 1.0 but it was scored using a different mechanism
                                    island_Pscore_if_het_list.insert(end_anchor_index,winning_score)
                                    island_pos_list.insert(end_anchor_index,anchor_desp_pos_1)
                                    island_pos_list_2.insert(end_anchor_index,anchor_desp_pos_2)
                                    island_pos_list_ave.insert(end_anchor_index,int(np.round((anchor_desp_pos_1 + anchor_desp_pos_2) / 2,0)))
                                    island_region_text.insert(end_anchor_index,str(L_anchor_holder) + "/" + str(R_anchor_holder))
                                    island_region_L.insert(end_anchor_index,L_anchor_holder)
                                    island_region_R.insert(end_anchor_index,R_anchor_holder)
                                    go_again_for_duplicate_culling.append(0)
                                    go_again_for_noncontiguous_culling.append(1)
                                    break

                            elif start_anchor_index+1 == end_anchor_index-1:
                                if island_region_L[start_anchor_index+1] != island_region_R[start_anchor_index] or island_region_R[start_anchor_index+1] != island_region_L[end_anchor_index]:
                                    indexes_to_cull.append(start_anchor_index+1)
                                    go_again_for_duplicate_culling.append(0)
                                    go_again_for_noncontiguous_culling.append(0)
                                    break
                                else:
                                    pass

                            else:
                                # The following indexes work out-- we want to start the index after the first anchor, BUT be 2 behind the last one so we can
                                # still do s+1 to get the site just before the end...
                                for s in range(start_anchor_index+1,end_anchor_index-1):
                                    if verbose_output:
                                        print("general info")
                                        print(c+1)
                                        print(island_region_text)
                                        print(s)
                                        print("Anchors")
                                        print(start_anchor_index)
                                        print(end_anchor_index)
                                        print("island_regions")
                                        print(island_region_L[s])
                                        print(island_region_R[s])
                                        print(island_region_L[s+1])
                                        print(island_region_R[s+1])
                                        print("anchor regions")
                                        print(island_region_L[start_anchor_index])
                                        print(island_region_R[start_anchor_index])
                                        print(island_region_L[end_anchor_index])
                                        print(island_region_R[end_anchor_index])

                                        print("island_region_L[s]")
                                        print(island_region_L[s])
                                        print("island_region_L[s+1]")
                                        print(island_region_L[s+1])
                                        print("island_region_R[s]")
                                        print(island_region_R[s])
                                        print("island_region_R[s+1]")
                                        print(island_region_R[s+1])

                                        print("indexes_to_cull at big")
                                        print(indexes_to_cull)

                                    #see if its an easy fix and if they're the exact same just get the best score
                                    if island_region_L[s] == island_region_L[s+1] and island_region_R[s] == island_region_R[s+1]:
                                        go_again_for_duplicate_culling.append(1)
                                        go_again_for_noncontiguous_culling.append(0)
                                        if island_region_L[s] == "HET" or island_region_R[s] == "HET":
                                            if c == -1: print("just hets!")
                                            if island_Pscore_if_het_list[s] > island_Pscore_if_het_list[s+1]:
                                                indexes_to_cull.append(s+1)
                                                break_out_flag = True
                                                break
                                            else:
                                                indexes_to_cull.append(s)
                                                break_out_flag = True
                                                break
                                        else:
                                            if island_score_list[s] > island_score_list[s+1]:
                                                indexes_to_cull.append(s+1)
                                                break_out_flag = True
                                                break
                                            else:
                                                indexes_to_cull.append(s)
                                                break_out_flag = True
                                                break
                                    # Now to test for NONCONTIC culling
                                    elif s == start_anchor_index+1 and island_region_L[s] != island_region_R[start_anchor_index]:
                                        if c == -1: print("noncontig 1")
                                        indexes_to_cull.append(s)
                                        go_again_for_noncontiguous_culling.append(1)
                                        break_out_flag = True
                                        break
                                    elif s == end_anchor_index-2 and island_region_R[s+1] != island_region_L[end_anchor_index]:
                                        if c == -1: print("noncontig 2")
                                        indexes_to_cull.append(s+1)
                                        go_again_for_noncontiguous_culling.append(1)
                                        break_out_flag = True
                                        break

                                    elif island_region_R[s] != island_region_L[s+1]:
                                        if c == -1: print("noncontig 3")
                                        go_again_for_duplicate_culling.append(0)
                                        go_again_for_noncontiguous_culling.append(1)
                                        if island_region_L[s] == island_region_R[start_anchor_index] and island_region_R[s] == island_region_L[end_anchor_index]:
                                            if c == -1: print("noncontig 4")
                                            indexes_to_cull.append(s+1)
                                            break_out_flag = True
                                            break
                                        elif island_region_L[s+1] == island_region_R[start_anchor_index] and island_region_R[s+1] == island_region_L[end_anchor_index]:
                                            if c == -1: print("noncontig 5")
                                            indexes_to_cull.append(s)
                                            break_out_flag = True
                                            break
                                        # Check to see if it is a messy end to a H site-- P1/H, H/P2, H/P1, P1/P1
                                        elif (island_region_L[s] == "HET" and island_region_L[s+1] == "HET") or (island_region_R[s] == "HET" and island_region_R[s+1] == "HET"):
                                            if c == -1: print("noncontig 6")
                                            if island_Pscore_if_het_list[s] > island_Pscore_if_het_list[s+1]:
                                                indexes_to_cull.append(s+1)
                                                break_out_flag = True
                                                break
                                            else:
                                                indexes_to_cull.append(s)
                                                break_out_flag = True
                                                break
                                        # As a last ditch effort, if the logic still does not work out BUT would if THIS Het recombination site were removed
                                        if s > start_anchor_index+1:
                                            if (island_region_L[s] == "HET" or island_region_R[s] == "HET") and (island_region_R[s-1] == island_region_L[s+1]):
                                                indexes_to_cull.append(s)
                                                break_out_flag = True
                                                break
                                        if s < end_anchor_index-2:
                                            if (island_region_L[s+1] == "HET" or island_region_R[s+1] == "HET") and (island_region_R[s] == island_region_L[s+2]):
                                                indexes_to_cull.append(s+1)
                                                break_out_flag = True
                                                break
                                        if (island_region_R[s] == "HET"):
                                            indexes_to_cull.append(s)
                                            break_out_flag = True
                                            break
                                        if (island_region_L[s+1] == "HET"):
                                            indexes_to_cull.append(s+1)
                                            break_out_flag = True
                                            break

                                        print("ERROR: ALGORITHM CANNOT FIGURE NONCONTIG CULLING OUT")
                                        print(input_file)
                                        print("failed general info")
                                        print("Chromosome "+str(c+1))
                                        print(island_region_text)
                                        print("Anchors")
                                        print(start_anchor_index)
                                        print(end_anchor_index)
                                        print("island_regions")
                                        print(island_region_L[s])
                                        print(island_region_R[s])
                                        print(island_region_L[s+1])
                                        print(island_region_R[s+1])
                                        print("anchor regions")
                                        print(island_region_L[start_anchor_index])
                                        print(island_region_R[start_anchor_index])
                                        print(island_region_L[end_anchor_index])
                                        print(island_region_R[end_anchor_index])

                                        go_again_for_duplicate_culling.append(0)
                                        go_again_for_noncontiguous_culling.append(0)
                                        go_again = False
                                    else:
                                        if c == -1: print("A success")
                                        go_again_for_duplicate_culling.append(0)
                                        go_again_for_noncontiguous_culling.append(0)

                            if break_out_flag == True:
                                go_again = True
                                break

                    if len(indexes_to_cull) == 0 and sum(go_again_for_duplicate_culling) == 0 and sum(go_again_for_noncontiguous_culling) == 0:
                        go_again = False
                        if verbose_output:
                            print("Ah, respite")

                    if verbose_output:
                        print(go_again)
                        print("indexes_to_cull before i actually do it")
                        print(indexes_to_cull)
                    cull_from_list(island_index_list,indexes_to_cull)
                    cull_from_list(island_score_list,indexes_to_cull)
                    cull_from_list(island_Pscore_if_het_list,indexes_to_cull)
                    cull_from_list(island_pos_list,indexes_to_cull)
                    cull_from_list(island_pos_list_2,indexes_to_cull)
                    cull_from_list(island_pos_list_ave,indexes_to_cull)
                    cull_from_list(island_region_text,indexes_to_cull)
                    cull_from_list(island_region_L,indexes_to_cull)
                    cull_from_list(island_region_R,indexes_to_cull)
                    indexes_to_cull = []

            # Cull all anchor regions out
            for a in range(len(island_index_list)):
                if island_region_L[a] == island_region_R[a]:
                    indexes_to_cull.append(a)

            cull_from_list(island_index_list,indexes_to_cull)
            cull_from_list(island_score_list,indexes_to_cull)
            cull_from_list(island_Pscore_if_het_list,indexes_to_cull)
            cull_from_list(island_pos_list,indexes_to_cull)
            cull_from_list(island_pos_list_2,indexes_to_cull)
            cull_from_list(island_pos_list_ave,indexes_to_cull)
            cull_from_list(island_region_text,indexes_to_cull)
            cull_from_list(island_region_L,indexes_to_cull)
            cull_from_list(island_region_R,indexes_to_cull)
            indexes_to_cull = []

            if verbose_output:
                print("after anchor cull")
                print(island_region_text)
                print(island_pos_list)
                print(island_pos_list_2)

            for t in range(len(island_score_list)):
                row = []
                row.append(s_chrom[c][0])
                row.append(island_pos_list[t])
                row.append(island_pos_list_2[t])
                row.append(int(np.round((island_pos_list[t] + island_pos_list_2[t]) / 2,0)))
                row.append(island_region_text[t])
                row.append(tinylist_bend_indexes_smol[island_index_list[t]])
                if island_region_L[t] == "HET" or island_region_R[t] == "HET":
                    row.append(np.round(island_Pscore_if_het_list[t],3))
                else:
                    row.append(np.round(island_score_list[t],3))
                output2_basefile = output_file_2.split("/")[-1]
                row.append(output2_basefile.split(".txt_No.NA.no.H",maxsplit=1)[0])
                all.append(row)
        writer2.writerows(all)

out_file_path_1 = (single_input_file[:-3] + "roaming_indexes_xy.txt")
out_file_path_2 = (single_input_file[:-3] + "recombination_sites.txt")
print("input file:")
print(single_input_file)
determine_event_positions(single_input_file,out_file_path_1,out_file_path_2)