#include <iostream>
#include <fstream>
#include <string>
#include <cmath>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <armadillo>
#include <unistd.h>
#include <vector>
#include "MsUtil.h"
#include "MsPostCal.h"
#include "MsCaviarModel.h"

using namespace std;

/*
 Reads the content of the file, return a vector of paths
 @param fileName the name of file that contains all the paths to z or ld files
 @return vector of paths
 */
vector<string> read_dir(string fileName){
    vector<string> dirs;
    string data;

    ifstream fin(fileName.c_str(), std::ifstream::in);
    if (!fin) {
        cout << "Error: unable to open " << fileName << endl;
        exit(1); // terminate with error
    }

    while(fin.good()){
        getline(fin,data);
        if(data != "") {
            dirs.push_back(data); }}
    fin.close();
    return dirs;
}

vector<int> read_sigma(string sample_size) {
    vector<int> sizes;
    string current_size = "";
    for (int i=0; i < sample_size.size(); i++) {
        if (sample_size[i] == ',') {
            sizes.push_back(stod(current_size)); // push back previous sample size
            current_size = ""; // reset
        }
        else if (isdigit(sample_size[i])) {
            current_size += sample_size[i];
        }
        else {
            cout << "Error: sample size is not in the right format" << endl;
            exit(1);
        }
    }
    if (current_size != "") {
        sizes.push_back(stod(current_size)); // push back last sample size
    }
    return sizes;
}

int main( int argc, char *argv[]  ){
    int totalCausalSNP = 3;
    double gamma = 0.01;
    double rho = 0.95;
    bool histFlag = false;
    int oc = 0;
    double tau_sqr = 0.52;
    double sigma_g_squared = 5.2;

    string ldFile = "";
    string zFile  = "";
    string outputFileName = "";
    string sample_s = "";

    while ((oc = getopt(argc, argv, "vhl:o:z:r:c:g:f:t:s:n:")) != -1) {
        switch (oc) {
            case 'v':
                cout << "Version 0.1\n" << endl;
            case 'h':
                cout << "Options: " << endl;
                cout << "-h, --help                 show this help message and exit " << endl;
                cout << "-o OUTFILE, --out=OUTFILE  specify the output file" << endl;
                cout << "-l LDFile, --ld_file=LDFILE    REQUIRED: the ld input file that contains paths to ld files" << endl;
                cout << "-z ZFile, --z_file=ZFILE   REQUIRED: the z-score and rsID file that contains paths to z files" << endl;
                cout << "-r RHO, --rho-prob=RHO     set $rho$ probability (default 0.95)" << endl;
                cout << "-g GAMMA, --gamma      set $gamma$ the prior of a SNP being causal (default 0.01)" << endl;
                cout << "-c causal          set the maximum number of causal SNPs (default 3)" << endl;
                cout << "-f 1               to out the probaility of different number of causal SNP" << endl;
                cout << "-t TAU_SQR, --tau_sqr=TAU_SQR  set the heterogeneity (t^2) across studies, default is 0.52" << endl;
                cout << "-s SIGMA_G_SQR, --sigma_g_squared=SIGMA_G_SQR    set the NCP variance for the smallest study, default is 5.2" << endl;
                cout << "-n SAMPLE_SIZE, --sample_size    REQUIRED: set the sample sizes (integer) of individual studies, format: enter 50,100 for study 1 with 50 sample size and 2 with 100" << endl;
                exit(0);

            // required options: LD and Z-score filenames, and output file name
            case 'l':
                ldFile = string(optarg);
                break;
            case 'o':
                outputFileName = string(optarg);
                break;
            case 'z':
                zFile = string(optarg);
                break;
            case 'n':
                sample_s = string(optarg);
                break;
            // optional arguments: parameters for fine mapping
            case 'r':
                rho = atof(optarg);
                break;
            case 'c':
                totalCausalSNP = atoi(optarg);
                break;
            case 'g':
                gamma = atof(optarg);
                break;
            case 'f':
                histFlag = true;
                break;
            case 't':
                tau_sqr = atof(optarg);
                break;
            case 's':
                sigma_g_squared = atof(optarg);
                break;
            case ':':
            case '?':
            default:
                break;
        }
    }

    if (ldFile == "" or zFile == "" or outputFileName == "" or sample_s == "") {
        cout << "Error: -l, -z, -o, and -n are required" << endl;
        exit(1);
    }

    vector<string> ldDir = read_dir(ldFile);
    vector<string> zDir = read_dir(zFile);
    vector<int> sample_sizes = read_sigma(sample_s);

    if (ldDir.size() != zDir.size() || ldDir.size() != sample_sizes.size()) {
        cout << "Error: LD files, Z files, and sample sizes do not match in number" << endl;
        exit(1);
    }

    MCaviarModel Mcaviar(ldDir, zDir, sample_sizes, outputFileName, totalCausalSNP, rho, histFlag, gamma, tau_sqr, sigma_g_squared);
    Mcaviar.run();
    Mcaviar.finishUp();
    return 0;
}
