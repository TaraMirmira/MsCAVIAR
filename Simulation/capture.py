import sys
import argparse

def read_set(read_fn):
    f = open(read_fn,'r')
    causal_set = set()
    for line in f:
        line = line.strip()
        causal_set.add(line)
    return causal_set

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='This program takes in the output set of MCaviar.py and the'
    	'true_causal set and calculate the recall rate and configuration size.')
    parser.add_argument('-s1', '--speculated', required=True, dest='speculated_file1',
                        help='speculated set file')
    parser.add_argument('-s2', '--speculated2', required=False, dest='speculated_file2',
                        help='speculated set file 2, for caviar in multiple population')
    parser.add_argument('-s3', '--speculated3', required=False, dest='speculated_file3',
                        help='speculated set file 3, for caviar in multiple population')
    parser.add_argument('-s4', '--speculated4', required=False, dest='speculated_file4',
                        help='speculated set file 4, for caviar in multiple population')
    parser.add_argument('-t', '--true', required=True, dest='true_file',
                        help='true causal set file')
    parser.add_argument('-p', '--parameter', required=True, dest='param',
                        help='the parameter we are manipulating')
    parser.add_argument('-w', '--which', required=True, dest='software',
                        help='the software used (mscaviar, caviar, or paintor')

args = parser.parse_args()
s1_fn = args.speculated_file1
if args.speculated_file2 is not None:
    s2_fn = args.speculated_file2
if args.speculated_file3 is not None:
    s3_fn = args.speculated_file3
if args.speculated_file4 is not None:
    s4_fn = args.speculated_file4
t_fn = args.true_file
parameter = args.param
software = args.software

if software == "mscaviar" or software == "paintor":
    speculated_set = read_set(s1_fn)
    true_set = read_set(t_fn)

    intersect_set = speculated_set.intersection(true_set)

    f1 = open(parameter + "_recall_rate.txt",'a+')
    f1.write(str(float(len(intersect_set)) / float(len(true_set))) + "\n")
    f1.close() 
    f2 = open(parameter + "_config_size.txt",'a+')
    f2.write(str(len(speculated_set)) + "\n")
    f2.close()

elif software == "caviar":
    # intersect and union
    # speculated_set1 = read_set(s1_fn)
    # speculated_set2 = read_set(s2_fn)
    # true_set = read_set(t_fn)

    # intersect_result = speculated_set1.intersection(speculated_set2)
    # union_result = speculated_set1.union(speculated_set2)

    # true_intersect_set = intersect_result.intersection(true_set)
    # true_union_set = union_result.intersection(true_set)

    # f1 = open(parameter + "_intersect_recall_rate.txt",'a+')
    # f1.write(str(float(len(true_intersect_set)) / float(len(true_set))) + "\n")
    # f1.close() 
    # f2 = open(parameter + "_intersect_config_size.txt",'a+')
    # f2.write(str(len(intersect_result)) + "\n")
    # f2.close()

    # f3 = open(parameter + "_union_recall_rate.txt",'a+')
    # f3.write(str(float(len(true_union_set)) / float(len(true_set))) + "\n")
    # f3.close() 
    # f4 = open(parameter + "_union_config_size.txt",'a+')
    # f4.write(str(len(union_result)) + "\n")
    # f4.close()

    # separated as two populations
    speculated_set1 = read_set(s1_fn)
    speculated_set2 = read_set(s2_fn)
    true_set = read_set(t_fn)

    true_pop1_set = speculated_set1.intersection(true_set)
    true_pop2_set = speculated_set2.intersection(true_set)

    f1 = open(parameter + "_pop1_recall_rate.txt",'a+')
    f1.write(str(float(len(true_pop1_set)) / float(len(true_set))) + "\n")
    f1.close() 
    f2 = open(parameter + "_pop1_config_size.txt",'a+')
    f2.write(str(len(speculated_set1)) + "\n")
    f2.close()

    f3 = open(parameter + "_pop2_recall_rate.txt",'a+')
    f3.write(str(float(len(true_pop2_set)) / float(len(true_set))) + "\n")
    f3.close() 
    f4 = open(parameter + "_pop2_config_size.txt",'a+')
    f4.write(str(len(speculated_set2)) + "\n")
    f4.close()
