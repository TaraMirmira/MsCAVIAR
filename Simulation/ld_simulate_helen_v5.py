"""
This script takes LD matrices and simulates summary statistics for multiple populations. Requires 2 populations, can take up to 4.
"""
import argparse
from contextlib import ExitStack
import numpy as np
import scipy.stats as stats
import os

class Main():

  def __init__(self):
    self.ld = []
    self.ld_matrices = []
    self.loci = []
    self.out = []
    self.causal = 1
    self.sims = 100
    self.tau_sqr = 0.5
    self.num_snps = 0 # num of snps in the input LD matrix
    self.upper = 5e-8
    self.lower = 0
    self.sample_size = []
    self.define_parser()  
    self.simulate()

  #define_parser is argparser
  def define_parser(self):
    parser = argparse.ArgumentParser(description = 'This script generates summary statistics using plink LD lists. Assumes LD matrices are whitespace delimited.')
    # should be now matrices not pairwise lists
    parser.add_argument('-l1', '--ld1', dest = 'ld1', required = True, help = 'Location of LD matrix for population 1')
    parser.add_argument('-l2', '--ld2', dest = 'ld2', required = True, help = 'Location of LD matrix for population 2')
    parser.add_argument('-l3', '--ld3', dest = 'ld3', required = False, help = 'Location of LD matrix for population 3')
    parser.add_argument('-l4', '--ld4', dest = 'ld4', required = False, help = 'Location of LD matrix for population 4')
    parser.add_argument('-o','--out_dir', dest = 'out', required = False, default = './', help = 'Directory for output files. Default: ./')
    parser.add_argument('-c', '--causal', dest = 'causal', required = False, default = 1, help = 'Number of causal SNP. Default: 1')
    parser.add_argument('-s', '--sims', dest = 'sims', required = False, default = 100, help = 'Number of simulations for each configuration. Default: 100')
    parser.add_argument('-t', '--tau_sqr', dest = 'tau_sqr', required = False, default = 0.5, help = 'Value for tau squared in simulations: Default: 0.5')
    parser.add_argument('-u', '--upperbound', dest = 'upperbound', required = False, default = 5e-8, help = 'Upperbound of the threshold for the significance of causal snps: Default: 5e-8')
    parser.add_argument('-l', '--lowerbound', dest = 'lowerbound', required = False, default = 0, help = 'Lowerbound of the threshold for the significance of causal snps: Default: 0')
    parser.add_argument('-n', '--sample_size', dest = 'sample_size', required = True, default = 0, help = 'Sample size of the studies')
    args = parser.parse_args()
    self.read_parser(args)

  #make args part of class
  def read_parser(self, args):
    splits = args.out.split('/')
    if len(splits) > 1:  # not in current dir
      output = '/'.join(splits[:]) + '/'
    else:
      output = splits[0]

    # read in LD matrices
    self.ld.append(args.ld1)
    self.ld.append(args.ld2)

    if args.ld3 is not None:
      self.ld.append(args.ld3)
    if args.ld4 is not None:
      self.ld.append(args.ld4)
    
    self.createFolder(output)
    for i in self.ld:
      self.ld_matrices.append(self.read_LD(i))
      self.out.append(output + i.strip().split('/')[-1].split('.')[0])
    self.num_snps = len(self.ld_matrices[0])
    self.loci = np.arange(0, self.num_snps)

    for i in range(len(self.ld_matrices)):
      if len(self.ld_matrices[i]) != len(self.ld_matrices[i][0]):
        raise ValueError('Matrix is not square')
      if len(self.ld_matrices[i]) != len(self.ld_matrices[0]):
        raise ValueError('Matrices are not the same size')

    self.causal = int(args.causal)
    # self.region = int(args.region)
    self.sims = int(args.sims)
    # print(self.sims)
    self.tau_sqr = float(args.tau_sqr)
    self.upper = float(args.upperbound)
    self.lower = float(args.lowerbound)
    self.sample_size = (args.sample_size).split(',')
    self.sample_size = [int(i) for i in self.sample_size]

    if len(self.sample_size) != len(self.ld):
      raise ValueError('number of sample_size is not the same as number of studies!!!!')

  def read_LD(self, read_fn):
    f = open(read_fn,'r')
    SIGMA = []
    for line in f:
        line = line.strip()
        array = line.split()
        SIGMA.append(array)
    return SIGMA

  def createFolder(self, directory):
    """
    this function can ONLY create if there is no folder already existing, if exist, will NOT overwrite
    """
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print('Error: creating directory. ' + directory)

  def test(self, lambda_c, p_list, study):
    count_sig = 0
    # check significant level
    if study == 0:
      # print("study #", study)
      for i in range(self.num_snps):
        if lambda_c[i] != 0: # causal snp
          if p_list[i] <= 5e-8:
            count_sig += 1
            continue
          else:
            # print(p_list[i])
            return False
        else:
          if p_list[i] <= 5e-8:
            count_sig += 1
    else: # other studies
      # print("study #", study)
      for i in range(self.num_snps):
        if lambda_c[i] != 0: # causal snp
          if p_list[i] <= self.upper and p_list[i] > self.lower:
            count_sig += 1
            continue
          else:
            return False
        else: # not causal snp
          if p_list[i] <= 5e-8:
            count_sig += 1
    if (count_sig / self.num_snps) < 0.05 or (count_sig / self.num_snps) > 0.8:
       # non-causal snps should be significant between 5%-50% of the time (not too hard or easy)
      return False
    return True # no False till the end of list

  def draw_causal(self):
    causals = []
    vals = sorted(np.random.choice(self.loci, size = self.causal, replace = False))
    for i in vals:
      causals.append(int(i))

    # NOTE: check for ld between causals for all studies
    for i in range(len(self.ld)):
      for j in range(len(causals)):
        for k in range(len(causals)):
          if j != k:
            if abs(float(self.ld_matrices[i][j][k])) > 0.7:
              causals = self.draw_causal() # redraw until causals have ld less than 70% among them
    return causals

  def simulate(self):
    super_beta = []
    pop_NCP = [[0 for j in range(self.causal)] for i in range(len(self.ld))]
    causals = self.draw_causal()

    sigma_g = 5.2 / np.sqrt(min(self.sample_size))
    super_beta = np.random.normal(loc = sigma_g, scale = 0.0000125, size = self.causal)
    # super_beta = np.random.normal(loc = 0.0, scale = np.sqrt(sigma_g), size = self.causal)
    for i in range(len(pop_NCP)):
      for j in range(len(pop_NCP[i])):
        pop_NCP[i][j] = np.random.normal(loc = super_beta[j] * np.sqrt(self.sample_size[i]), scale = np.sqrt(self.tau_sqr))

    for i in range(len(self.ld)):
      ld_matrix = np.array(self.ld_matrices[i], dtype=np.float64)
      lambda_c = np.zeros(self.num_snps)
      print(i, pop_NCP[i])
      lambda_c[causals] = pop_NCP[i]

      m_name = self.out[i] + '.info'
      m = open(m_name, 'w')
      for j in causals:
        m.write(str(j))
        m.write("\n")
      m.close()
      # np.savetxt(self.out[i] + '.info', np.transpose(np.array([lambda_c], dtype = np.float64)), delimiter = " ", header = "lambda_c")
      count_fail = 0

      l = 1
      while l < (self.sims + 1):
        z = np.random.multivariate_normal(np.dot(ld_matrix,lambda_c),ld_matrix)
        p = stats.norm.pdf(abs(z))*2
        if self.test(lambda_c, p, i):
          # caviar output
          f_name = self.out[i] + '_' + str(l) + '.caviar'
          f = open(f_name,'w')
          for k in range(self.num_snps):
              f.write(str(k) + " ")
              f.write(str(z[k]))
              f.write("\n")
          f.close()
          # paintor output
          p_name = self.out[i] + '_' + str(l) + '.paintor'
          p = open(p_name,'w')
          for k in range(self.num_snps): # for paintor's input file, generate required information
              p.write(str(0) + " ") # fake chr
              p.write(str(0) + " ") # fake pos
              p.write(str(k) + " ") # rsid, just 1-50 for simulations
              p.write("Z" + " ") # fake A0
              p.write("Z" + " ") # fake A1
              p.write(str(z[k]))
              p.write("\n")
          p.close()
          # np.savetxt(, np.transpose(np.array([self.loci, z, p])), delimiter = " ", header = "loci zscore pvalue")
          l += 1
        else:
          count_fail += 1
          if count_fail == 100000: # if failed too many times, redraw the causal variants
            self.simulate()
            exit()
if __name__ == '__main__':
  Main()

