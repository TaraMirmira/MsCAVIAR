#ifndef MCAVIARMODEL_H
#define MCAVIARMODEL_H

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <map>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <armadillo>

#include "MsPostCal.h"

using namespace std;
using namespace arma;

class MCaviarModel{
public:
    double rho;
    double gamma;
    int snpCount;
    int totalCausalSNP;
    vector<mat> * sigma;
    vector< vector<double> > * z_score;
    vector<char> * pcausalSet;
    vector<int> * rank;
    bool histFlag;  // to out the probaility of different number of causal SNP
    MPostCal * post;
    vector< vector<string> > * snpNames;
    vector<string> ldDir;
    vector<string> zDir;
    vector<int> sample_sizes;
    string outputFileName;
    double tau_sqr;
    double sigma_g_squared;
    int num_of_studies;
    vector<double> S_LONG_VEC;
    bool haslowrank = false;

    /*
     consrtuctor for MCaviarModel
     */
    MCaviarModel(vector<string> ldDir, vector<string> zDir, vector<int> sample_sizes, string outputFileName, int totalCausalSNP, double rho, bool histFlag, double gamma=0.01, double tau_sqr = 0.2, double sigma_g_squared = 5.2) {
        this->histFlag = histFlag;
        this->rho = rho;
        this->gamma = gamma;
        this->ldDir = ldDir;
        this->zDir  = zDir;
        this->outputFileName = outputFileName;
        this->totalCausalSNP = totalCausalSNP;
        this->tau_sqr = tau_sqr;
        this->sigma_g_squared = sigma_g_squared;
        this->sample_sizes = sample_sizes;

        //fileSize(ldFile, tmpSize);
        sigma      = new vector<mat>;
        z_score    = new vector<vector<double> >;
        snpNames   = new vector<vector<string> >;

        for(int i = 0; i < ldDir.size(); i++) {
            string ld_file = ldDir[i];
            string z_file = zDir[i];

            vector<double>* temp_LD = new vector<double>;
            vector<string> temp_names;
            vector<double> temp_z;

            importData(ld_file, temp_LD);
            importDataFirstColumn(z_file, temp_names);
            importDataSecondColumn(z_file, temp_z);

            int numSnps = sqrt(temp_LD->size());
            mat temp_sig;
            temp_sig = mat(numSnps, numSnps);
            for (int i = 0; i < numSnps; i++){
                for (int j = 0; j< numSnps; j++){
                    temp_sig(i,j) = temp_LD->at(i * numSnps + j);
                }
            }

            sigma->push_back(temp_sig);
            snpNames->push_back(temp_names);
            z_score->push_back(temp_z);

            delete temp_LD;
        }

        num_of_studies = snpNames->size();
        snpCount = (*snpNames)[0].size();
        pcausalSet = new vector<char>(snpCount);
        rank = new vector<int>(snpCount, 0);

        for (int i = 0; i < z_score->size(); i++){
            for(int j = 0; j < (*z_score)[i].size(); j++){
                S_LONG_VEC.push_back((*z_score)[i][j]);
            }
        }

        /* sigma_g_squared is set to max(5.2, max(abs(z-score)))
        for(int i = 0 ; i < num_of_studies; i++){
            for (int j = 0; j < snpCount; j++){
                if(abs(S_LONG_VEC.at(i*snpCount + j)) > sigma_g_squared){
                    sigma_g_squared = abs(S_LONG_VEC.at(i*snpCount + j));
                }
            }
        }
        */

        //make positive definite
        for (int i = 0; i < sigma->size(); i++){
            //check for low rank
            if(arma::rank(sigma->at(i)) < snpCount){
                haslowrank = true;
                std::cout << "study " << i << " has low rank. Implementing low_rank method.\n";
            }
            
            makeSigmaPositiveSemiDefinite(&(sigma->at(i)), snpCount);
        }

        mat* BIG_SIGMA = new mat(snpCount * num_of_studies, snpCount * num_of_studies, fill::zeros);
        for (int i = 0 ; i < num_of_studies; i++){
            mat temp_sigma = mat(num_of_studies , num_of_studies, fill::zeros);
            temp_sigma(i,i) = 1;
            temp_sigma = kron(temp_sigma, sigma->at(i));
            (*BIG_SIGMA) = (*BIG_SIGMA) + temp_sigma;
        }

        //if low rank, BIG_SIGMA = BIG_B, Stat matrix has new distribution
        if(haslowrank == true){
            //construct big B
            mat* BIG_B = new mat(snpCount*num_of_studies, snpCount*num_of_studies,fill::zeros);
            for(int i = 0; i<num_of_studies; i++){
                mat* tmpmat = new mat(snpCount,snpCount,fill::zeros);
                *tmpmat = BIG_SIGMA->submat(i*snpCount,i*snpCount,(i+1)*snpCount-1,(i+1)*snpCount-1);
                mat* tmpOmega = new mat(snpCount,snpCount,fill::zeros);
                tmpOmega = eigen_decomp(tmpmat,snpCount);

                //construct B
                mat trans_Q = trans(*tmpmat);
                mat sqrt_Omega = sqrt(*tmpOmega);
                mat B_each = sqrt_Omega * trans_Q;

                //merge to Big B
                mat temp_b = mat(num_of_studies , num_of_studies, fill::zeros);
                temp_b(i,i) = 1;
                temp_b = kron(temp_b, B_each);
                (*BIG_B) = (*BIG_B) + temp_b;

                //update S_LONG_VEC
                mat* z_score = new mat(snpCount,1,fill::zeros);
                for(int j = 0; j < snpCount; j++){
                    (*z_score)(j,0) = S_LONG_VEC[i*snpCount+j];
                }

                mat tmpS = inv(sqrt_Omega) * trans_Q;
                mat lowS = tmpS * (*z_score);

                for(int j = 0; j < snpCount; j++){
                    S_LONG_VEC[i*snpCount+j] = lowS(j,0);
                }
                delete(tmpmat);
                delete(tmpOmega);
                delete(z_score);
            }

            *BIG_SIGMA = *BIG_B;
            delete(BIG_B);
        }	

        post = new MPostCal(BIG_SIGMA, &S_LONG_VEC, snpCount, totalCausalSNP, snpNames, gamma, tau_sqr, sigma_g_squared, num_of_studies, sample_sizes, haslowrank);
    }

    /*
     run the greedy algorithm
     @param no param
     @return no return
     */
    void run() {
        post->findOptimalSetGreedy(&S_LONG_VEC, sigma_g_squared, pcausalSet, rank, rho, outputFileName);
    }

    /*
     finish by by printing the set, post and hist file
     @param no param
     @return no return
     */
    void finishUp() {
        ofstream outputFile;
        string outFileNameSet = string(outputFileName)+"_set.txt";
        outputFile.open(outFileNameSet.c_str());
        for(int i = 0; i < snpCount; i++) {
            if((*pcausalSet)[i] == '1')
                outputFile << (*snpNames)[0][i] << endl;
        }
        post->printPost2File(string(outputFileName)+"_post.txt");
        //outputs the histogram data to file
        if(histFlag)
            post->printHist2File(string(outputFileName)+"_hist.txt");
    }

    // destructor
    ~MCaviarModel() {
        delete z_score;
        delete sigma;
        delete snpNames;
        delete pcausalSet;
        delete rank;
    }
};

#endif
