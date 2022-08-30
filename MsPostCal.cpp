#include <vector>
#include <unordered_map>
#include <algorithm>
#include <set>
#include <iostream>
#include <armadillo>
#include <iomanip>
#include <vector>
#include <math.h>
#include "MsUtil.h"
#include "MsPostCal.h"

#include <omp.h>
#include <ctime>

using namespace arma;

vector<vector<int>> MPostCal::get_causal_idxs(vector<int> all_configs) {
    vector<vector<int>> causal_idxs;
    for ( int i = 0; i < num_of_studies; i++ ) {
        vector<int> causal_idxs_studyi;
        int start_loc = std::accumulate(num_snps_all.begin(), num_snps_all.begin()+i, 0);
	for ( int j = 0; j < num_snps_all[i]; j++ ) {
	    if ( all_configs[start_loc + j] == 1 ) {
                causal_idxs_studyi.push_back(j);
	    }
	}
	causal_idxs.push_back(causal_idxs_studyi);
    }
    return causal_idxs;
}

// calibrate for sample size imbalance
mat MPostCal::construct_diagC(vector<int> configure, vector<int> all_configs) {
    mat Identity_M = mat(num_of_studies, num_of_studies, fill::eye);
    mat Matrix_of_sigmaG = mat(num_of_studies, num_of_studies);
    int min_size = * std::min_element(sample_sizes.begin(), sample_sizes.end());

    /*for (int i = 0; i < num_of_studies; i ++) {
        for (int j = 0; j < num_of_studies; j ++) {
            if (i == j) // diagonal: scaled variance
                Matrix_of_sigmaG(i, j) = s_squared * (double(sample_sizes[i]) / min_size);
            else // off-diagonal: covariance
                Matrix_of_sigmaG(i, j) = s_squared * sqrt(long(sample_sizes[i]) * long(sample_sizes[j])) / min_size;
        }
    } */   
    //mat temp1 = t_squared * Identity_M + Matrix_of_sigmaG;
    //mat temp1 = Matrix_of_sigmaG;
    //mat temp2 = mat(snpCount, snpCount, fill::zeros);
    //for(int i = 0; i < snpCount; i++) {
    //    if (configure[i] == 1)
    //       temp2(i, i) = 1;
    //}
    vec diagC_main_diag(all_configs.size(), fill::zeros);
    vector<vector<int>> causal_idxs = get_causal_idxs(all_configs);
    for ( int i = 0; i < causal_idxs.size(); i++ ) {
	int start_loc = std::accumulate(num_snps_all.begin(), num_snps_all.begin()+i, 0);
        for ( int j = 0; j < causal_idxs[i].size(); j++ ) {
            int causal_pos = causal_idxs[i][j];
            diagC_main_diag[start_loc + causal_pos] = s_squared * (double(sample_sizes[i]) / min_size);
	}
    }


    mat diagC = diagmat(diagC_main_diag);
    //mat diagC = kron(temp1, temp2);
    
    //fill in off diagonal
    for ( int i = 0; i < causal_idxs.size(); i++ ) {
	bool causal_in_all = true;
        int causal_in_0 = causal_idxs[0][i];
	vector<int> main_diag_pos;
	main_diag_pos.push_back(causal_in_0);
        int start_loc = std::accumulate(num_snps_all.begin(), num_snps_all.begin()+i, 0);
        for ( int j = 1; j < num_of_studies; j ++ ) {
	    string snp_name = (*SNP_NAME)[0][causal_in_0];
	    int loc_in_j = (*snp_to_idx_all)[j][snp_name];
	    if ( all_configs[start_loc + loc_in_j] == 1 ) {
	       main_diag_pos.push_back(start_loc + loc_in_j);
	    } else {
	        causal_in_all = false;
		continue;
	    }
        }
	if (causal_in_all) {
            for ( int j = 0; j < main_diag_pos.size(); j++ ) {
		int first = main_diag_pos[j];
                for ( int k = 0; k < main_diag_pos.size(); k++ ) {
                   int second = main_diag_pos[k];
		   if ( first == second ) {
                       //main diag position, need to add t_squared
		       diagC(first, first) += t_squared;
		   } else {
		       //off diag pos, need to add s_squared
		       diagC(first, second) = s_squared * sqrt(long(sample_sizes[i]) * long(sample_sizes[j])) / min_size;
		       diagC(second, first) = s_squared * sqrt(long(sample_sizes[i]) * long(sample_sizes[j])) / min_size;
		   }
		}
	    }
	}
    }
    
    return diagC;
}

double MPostCal::likelihood(vector<int> configure, vector<int> all_configs, vector<double> * stat, double sigma_g_squared) {
    //int causalCount = 0;
    int totalCausalCount = 0;
    int totalSnpCount = all_configs.size();
    double matDet   = 0;
    double res      = 0;

    for(int i = 0; i < all_configs.size(); i++)
        totalCausalCount += configure[i];
    if(totalCausalCount == 0){ //TODO how to handle this case for multiple causal vectors
        mat tmpResultMatrixNM = statMatrixtTran * invSigmaMatrix;
        mat tmpResultMatrixNN = tmpResultMatrixNM * statMatrix;

        res = tmpResultMatrixNN(0,0);
        matDet = sigmaDet;
        return (-res/2-sqrt(abs(matDet)));
    }

    mat sigmaC = construct_diagC(configure, all_configs);
    int index_C = 0;
    mat sigmaMatrixTran = sigmaMatrix.t();

    // U is kn by mn matrix of columns corresponding to causal SNP in sigmacC
    // In unequal sample size studies, U is adjusted for the sample sizes
    //mat U(causalCount * num_of_studies, snpCount * num_of_studies, fill::zeros);
    mat U(totalCausalCount, totalSnpCount, fill::zeros);
    for (int i = 0; i < totalSnpCount; i++) {
        if (configure[i] == 1) {
            for (int j = 0; j < totalSnpCount; j++) {
                U(index_C, j) = sigmaC(i, j);
            }
            index_C ++;
        }
    }
    
    index_C = 0;
    
    // V is mn by kn matrix of rows corresponding to causal SNP in sigma
    // In unequal sample size studies, V does not change
    //mat V(causalCount * num_of_studies, snpCount * num_of_studies, fill::zeros);
    mat V(totalCausalCount, totalSnpCount, fill::zeros);
    for (int i = 0; i < totalSnpCount; i++) {
        if (configure[i] == 1) {
            for (int j = 0; j < totalSnpCount; j++) {
                V(index_C, j) = sigmaMatrixTran(i, j);
            }
            index_C ++;
        }
    }
    V = V.t();

    // UV = SigmaC * Sigma (kn by kn)
    mat UV(totalCausalCount, totalCausalCount, fill::zeros);
    UV = U * V;

    //mat I_AA   = mat(snpCount, snpCount, fill::eye); //can ignore, does not get used here
    mat tmp_CC = mat(totalCausalCount, totalCausalCount, fill::eye) + UV;
    matDet = det(tmp_CC) * sigmaDet;

    mat temp1 = invSigmaMatrix * V;
    mat temp2 = mat(totalSnpCount, totalCausalCount, fill::zeros);
    
    #pragma omp critical
    temp2 = temp1 * pinv(tmp_CC);

    mat tmp_AA = invSigmaMatrix - temp2 * U ;

    mat tmpResultMatrix1N = statMatrixtTran * tmp_AA;
    mat tmpResultMatrix11 = tmpResultMatrix1N * statMatrix;
    res = tmpResultMatrix11(0,0);

    if(matDet==0) {
        cout << "Error the matrix is singular and we fail to fix it." << endl;
        exit(0);
    }

    /*
     We compute the log of -res/2-log(det) to see if it is too big or not.
     In the case it is too big we just make it a MAX value.
     */
    double tmplogDet = log(sqrt(abs(matDet)));
    double tmpFinalRes = -res/2 - tmplogDet;

    return tmpFinalRes;
}

//here we still use Woodbury matrix, here sigma_matrix is B, and S is updated already
double MPostCal::lowrank_likelihood(vector<int> configure, vector<int> all_configs, vector<double> * stat, double sigma_g_squared) {
    //int causalCount = 0;
    int totalCausalCount = 0;
    int totalSnpCount = all_configs.size();
    double matDet   = 0;
    double res      = 0;

    for(int i = 0; i < totalSnpCount; i++)
        totalCausalCount += configure[i];
    if(totalCausalCount == 0){ //TODO ok for now, but need to think
        mat tmpResultMatrixNN = statMatrixtTran * statMatrix;
        res = tmpResultMatrixNN(0,0);
        matDet = 1;
        return (-res/2-sqrt(abs(matDet)));
    }

    mat sigmaC = construct_diagC(configure, all_configs);
    int index_C = 0;
    mat sigmaMatrixTran = sigmaMatrix.t();

    // In unequal sample size studies, U is adjusted for the sample sizes
    // here we make U = B * sigmaC, this is still kn by mn
    mat U(totalCausalCount, totalSnpCount, fill::zeros);

    mat small_sigma(totalSnpCount, totalCausalCount, fill::zeros);
    mat small_sigmaC(totalCausalCount, totalCausalCount, fill::zeros);
    for (int i = 0; i < totalSnpCount; i++) {
        if (configure[i] == 1) {
            for (int j = 0; j < totalSnpCount; j++) {
                small_sigma(j, index_C) = sigmaMatrix(j, i);
            }
            small_sigmaC(index_C, index_C) = sigmaC(i, i);
            index_C++;
        }
    }
    U = small_sigma * small_sigmaC;
    U = U.t();

    index_C = 0;

    // here V is B_trans, this is mn by kn
    mat V(totalCausalCount, totalSnpCount, fill::zeros);
    for (int i = 0; i < totalSnpCount; i++) {
        if (configure[i] == 1) {
            for (int j = 0; j < totalSnpCount; j++) {
                V(index_C, j) = sigmaMatrixTran(i, j);
            }
            index_C ++;
        }
    }
    V = V.t();

    // UV = B * SigmaC * Btrans (kn by kn)
    mat UV(totalCausalCount, totalCausalCount, fill::zeros);
    UV = U * V;

    mat I_AA   = mat(totalSnpCount, totalSnpCount, fill::eye);
    mat tmp_CC = mat(totalCausalCount, totalCausalCount, fill::eye) + UV;
    matDet = det(tmp_CC);

    mat temp2 = mat(totalSnpCount, totalCausalCount, fill::zeros);
    #pragma omp critical
    temp2 = V * pinv(tmp_CC);

    mat tmp_AA = I_AA - temp2 * U ;

    mat tmpResultMatrix1N = statMatrixtTran * tmp_AA;
    mat tmpResultMatrix11 = tmpResultMatrix1N * statMatrix;
    res = tmpResultMatrix11(0,0);
    
    if(matDet==0) {
        cout << "Error the matrix is singular and we fail to fix it." << endl;
        exit(0);
    }

    /*
     We compute the log of -res/2-log(det) to see if it is too big or not.
     In the case it is too big we just make it a MAX value.
     */
    double tmplogDet = log(sqrt(abs(matDet)));
    double tmpFinalRes = -res/2 - tmplogDet;

    return tmpFinalRes;
}

int MPostCal::nextBinary(vector<int>& data, int size) {
    int i = 0;
    int total_one = 0;
    int index = size-1;
    int one_countinus_in_end = 0;

    while(index >= 0 && data[index] == 1) {
        index = index - 1;
        one_countinus_in_end = one_countinus_in_end + 1;
    }

    if(index >= 0) {
        while(index >= 0 && data[index] == 0) {
            index = index - 1;
        }
    }
    if(index == -1) {
        while(i <  one_countinus_in_end+1 && i < size) {
            data[i] = 1;
            i=i+1;
        }
        i = 0;
        while(i < size-one_countinus_in_end-1) {
            data[i+one_countinus_in_end+1] = 0;
            i=i+1;
        }
    }
    else if(one_countinus_in_end == 0) {
        data[index] = 0;
        data[index+1] = 1;
    }
    else {
        data[index] = 0;
        while(i < one_countinus_in_end + 1) {
            data[i+index+1] = 1;
            if(i+index+1 >= size)
                printf("ERROR3 %d\n", i+index+1);
            i=i+1;
        }
        i = 0;
        while(i < size - index - one_countinus_in_end - 2) {
            data[i+index+one_countinus_in_end+2] = 0;
            if(i+index+one_countinus_in_end+2 >= size) {
                printf("ERROR4 %d\n", i+index+one_countinus_in_end+2);
            }
            i=i+1;
        }
    }
    i = 0;
    total_one = 0;
    for(i = 0; i < size; i++)
        if(data[i] == 1)
            total_one = total_one + 1;

    return(total_one);
}

vector<int> MPostCal::findConfig(int iter) {
    int numCausal = 0;
    int temp = iter;
    int sum = 0;
    vector<int> config(snpCount, 0);
    int comb = nCr(snpCount,numCausal);
    while(temp > comb) {
        temp = temp - comb;
        numCausal++;
        sum = sum + comb;
        comb = nCr(snpCount,numCausal);
    }

    int times = iter - sum; //this is the number of times we use find_next_binary
    for(int i = 0; i < numCausal; i++){
        config[i] = 1;
    }
    for(int i = 0; i < times; i++){
        temp = nextBinary(config, snpCount);
    }
    printf("num causal in this config is %d\n", temp);
    return config;
}

vector<int> MPostCal::findAllConfigs(vector<int> config) {
    //TODO this is temporary, must fix
    vector<int> all_configs;
    for (int i = 0; i < num_of_studies; i++ ) {
        all_configs.insert(all_configs.end(), config.begin(), config.end());
    }
    return all_configs;
}

double MPostCal::computeTotalLikelihood(vector<double>* stat, double sigma_g_squared) {
    double sumLikelihood = 0;
    long int total_iteration = 0 ;

    for(long int i = 0; i <= maxCausalSNP; i++)
        total_iteration = total_iteration + nCr(snpCount, i);
    cout << "Max Causal = " << maxCausalSNP << endl;

    //clock_t start = clock();
    vector<int> configure;
    vector<int> all_configs;
    int num;

    int chunksize;
    if(total_iteration < 1000){
        chunksize = total_iteration/10;
    }
    else{
        chunksize = total_iteration/1000;
    }
    int curr_iter = 0;

    #pragma omp parallel for schedule(static,chunksize) private(configure,num)
    for(long int i = 0; i < total_iteration; i++) {
        if(i%chunksize == 0){
            configure = findConfig(i);
        }
        else{
            num = nextBinary(configure, snpCount);
        }

	all_configs = findAllConfigs(configure);
        
        vector<int> tempConfigure(snpCount*num_of_studies,0);
        num = 0;
        int sizec = configure.size();
        double tmp_likelihood = 0;

        for(int k = 0; k < sizec; k++){
            if(configure[k] == 1) num += 1;
        }

        for (int m = 0; m < num_of_studies; m++){
            for (int j = 0; j < sizec; j++){
                tempConfigure[snpCount * m + j] = configure[j];
            }
        }

	tempConfigure = all_configs;
        
        if(haslowrank==true){
            tmp_likelihood = lowrank_likelihood(tempConfigure, all_configs, stat, sigma_g_squared) + num * log(gamma) + (snpCount-num) * log(1-gamma);
        }
        else{
            tmp_likelihood = likelihood(tempConfigure, all_configs, stat, sigma_g_squared) + num * log(gamma) + (snpCount-num) * log(1-gamma);
        }
        
        #pragma omp critical
        sumLikelihood = addlogSpace(sumLikelihood, tmp_likelihood);

        for(int j = 0; j < snpCount; j++) {
            for(int k = 0; k < num_of_studies; k++){
                #pragma omp critical
                postValues[j] = addlogSpace(postValues[j], tmp_likelihood * configure[j]);
            }
        }        

        #pragma omp critical
        if(i % 1000 == 0 and i > curr_iter){
            cerr << "\r                                                                 \r" << (double) (i) / (double) total_iteration * 100.0 << "%";
            curr_iter = i;
        }
    }

    //cout << "\ncomputing likelihood of all configurations took  " << (float)(clock()-start)/CLOCKS_PER_SEC << "seconds.\n";

    for(int i = 0; i <= maxCausalSNP; i++)
        histValues[i] = exp(histValues[i]-sumLikelihood);
    
    return(sumLikelihood);
}


vector<char> MPostCal::findOptimalSetGreedy(vector<double> * stat, double sigma_g_squared, vector<int> * rank,  double inputRho, string outputFileName, double cutoff_threshold) {
    int index = 0;
    double rho = double(0);
    double total_post = double(0);

    vector<char> causalSet(snpCount,'0');

    totalLikeLihoodLOG = computeTotalLikelihood(stat, sigma_g_squared);

    export2File(outputFileName+"_log.txt", exp(totalLikeLihoodLOG)); //Output the total likelihood to the log File
    for(int i = 0; i < snpCount; i++)
        total_post = addlogSpace(total_post, postValues[i]);
    printf("\nTotal Likelihood = %e SNP=%d \n", total_post, snpCount);

    std::vector<data> items;
    std::set<int>::iterator it;
    //output the poster to files
    for(int i = 0; i < snpCount; i++) {
        //printf("%d==>%e ",i, postValues[i]/total_likelihood);
        items.push_back(data(exp(postValues[i]-total_post), i, 0));
    }
    printf("\n");
    std::sort(items.begin(), items.end(), by_number());
    for(int i = 0; i < snpCount; i++)
        (*rank)[i] = items[i].index1;
    
    double threshold = cutoff_threshold;
    /*
    if(snpCount > 30){
        threshold = 1/(double)snpCount;
    }
    else{
        threshold = 0.1/(double)snpCount;
    }
    */
    cout << "threshold is " << threshold << "\n";
    
    while(rho < inputRho){
        rho += exp(postValues[(*rank)[index]]-total_post);
        if(exp(postValues[(*rank)[index]]-total_post) > threshold){
            causalSet[(*rank)[index]] = '1';
            printf("%d %e\n", (*rank)[index], rho);
        }
        index++;
    }

    /*
    do{
        rho += exp(postValues[(*rank)[index]]-total_post);
        causalSet[(*rank)[index]] = '1';
        printf("%d %e\n", (*rank)[index], rho);
        index++;
    } while( rho < inputRho);
    */
    printf("\n");
    return(causalSet);
}
