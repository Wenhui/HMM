#include "hmm.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include <fstream>
#include <iostream>
#include <string.h>
#include <vector>

using namespace std;

Model::Model(int stateCount) : N(stateCount) {
    AllocateMemory();
}

Model::Model(char *modelFile) {
    ifstream ifs(modelFile);

    if (ifs.fail()) {   
        cerr << "Failed to open the model file:"<< modelFile << endl;
        exit(1);
    }

    ifs >> N; // the first is the number of states
    AllocateMemory();

    char str[1000];
    int entryCount;

    // load initial state prob
    ifs >> str >> entryCount;
    if (strcmp(str, "InitPr")) {
        cerr << "Error: InitPr expected in model file\n";
    }
    int i;
    for (i=0; i<N; i++) {
        I[i]=0;
    }
    int s;
    double pr;
    for (i=0; i<entryCount; i++) {
        ifs >> s >> pr;
        I[s]=pr;
    }

    // load output prob
    ifs >> str >> entryCount;
    if (strcmp(str, "OutputPr")) {
        cerr << "Error: OutputPr expected in model file\n";
    }
    int j;
    for (i=0; i<N; i++) {
        for (j=0; j<SYMNUM; j++) {
            B[j][i] = 0;
        }
    }
    char sym;
    for (i=0; i<entryCount; i++) {
        ifs >> s >> sym >> pr;
        B[sym-baseChar][s]=pr;
    }

    // load state transition prob
    ifs >> str >> entryCount;
    if (strcmp(str, "TransPr")) {
        cerr << "Error: TransPr expected in model file\n";
    }
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            A[i][j] = 0;
        }
    }
    int s1;
    for (i=0; i<entryCount; i++) {
        ifs >> s >> s1 >> pr;
        A[s][s1]=pr;
    }

}


void Model::AllocateMemory()
{
    I.resize(N);
    A.resize(N);
    B.resize(SYMNUM);
    int i;
    for (i=0; i< N; i++) {
        A[i].resize(N);
    }
    for (i=0; i< SYMNUM; i++) {
        B[i].resize(N);
    }

    ICounter.resize(N);
    ACounter.resize(N);
    BCounter.resize(SYMNUM);
    for (i=0; i< N; i++) {
        ACounter[i].resize(N);
    }
    for (i=0; i< SYMNUM; i++) {
        BCounter[i].resize(N);
    }
    INorm = 0;
    ANorm.resize(N);
    BNorm.resize(N);
}


void Model::ResetCounter()
{
    int i, j;
    for (i=0; i< N; i++) {
        ICounter[i] = 0;

        for (j=0; j<N; j++ ) {
            ACounter[i][j] = 0;
        }
        ANorm[i] =BNorm[i] = 0;
        for (j=0; j<SYMNUM; j++) {
            BCounter[j][i] = 0;
        }
    }
    INorm = 0;
}



void Model::UpdateParameter()
{
    int i, j;
    for (i=0; i< N; i++) {
        I[i] = ICounter[i]/INorm;

        for (j=0; j<N; j++ ) {
            A[i][j] = ACounter[i][j]/ANorm[i];
        }
        for (j=0; j<SYMNUM; j++) {
            B[j][i] = BCounter[j][i]/BNorm[i];
        }
    }
}


Model::~Model()
{
  A.clear();
  B.clear();
  I.clear();
  ACounter.clear();
  BCounter.clear();
  ANorm.clear();
  BNorm.clear();
}


void Model::CountSequence(char *seqFile)
{

    ifstream ifs(seqFile);
    if (ifs.fail()) {
        cerr << "can't open sequence file: "<< seqFile <<endl;
        exit(1);
    }

  // ======== insert your code for counting ============
  // - count how many times it starts with state i
  // - count how many times a particular transition happens
  // - count how many times a particular symbol would be generated from a particular state
  // Increase the corresponding counters as well as the normalizers

  // Use something like the following to iterate over all the characters in the stream
  // It's important that you transform the character to the corresponding true
  // index by subtracting the "basechar" as shown below!
  //
    char c;
    int s;
    int trueIndex;
    int prev_s;

    ifs >> c >> s;
   
    trueIndex = c-baseChar;
    ICounter[s]++;
    BCounter[trueIndex][s]++;
    prev_s = s;

    while (ifs >> c >> s) {
        trueIndex = c-baseChar;
      //    Now you can use something like B[trueIndex][i] to refer
      //    to the entry in B.
        
      // ======= Do your counting here...
      // e.g.,   change ICounter, INorm, ACounter, ANorm, etc


        BCounter[trueIndex][s]++;
        ACounter[prev_s][s]++;

        prev_s = s;
    }

    for (int i = 0; i < N; i++){
        INorm += ICounter[i];
        for (int j = 0; j < N; j++) {
            ANorm[i] += ACounter[i][j];
        }
        for (int j = 0; j < SYMNUM; j++) {
            BNorm[i] += BCounter[j][i];
        }
    }
}

// void Model::UpdateParameter()
// {
//     int i, j;
//     for (i=0; i< N; i++) {
//         I[i] = ICounter[i]/INorm;

//         for (j=0; j<N; j++ ) {
//             A[i][j] = ACounter[i][j]/ANorm[i];
//         }
//         for (j=0; j<SYMNUM; j++) {
//             B[j][i] = BCounter[j][i]/BNorm[i];
//         }
//     }
// }


void Model::RandomInit(int *sym, int seqLen)
{
    /* seed the random number generator */
    srand48(time(NULL));

    int i,j;
    double sum,sumI;
    sumI = 0;

    for (i = 0; i < N; i++) {

        /* initialize the I(state_t) vector */
        I[i] = 0.99/N + 0.01*drand48();
        sumI += I[i];
        /* initialize the A(state_t, state_t+1) matrix */
        sum = 0.0;
        for (j = 0; j < N; j++) {
            A[i][j] = 0.99 / N + 0.01 * drand48();
            sum += A[i][j];
        }
        for (j = 0; j < N; j++)
            A[i][j] /= sum;

        /* initialize the B(output,state) matrix */
        sum = 0.0;
        for (j=0; j<seqLen;j++) {
            B[sym[j]][i] = 0.99 / N + 0.01 * drand48();
        }
        for (j = 0; j<SYMNUM; j++) {
            // B[j][i] = 0.99 / N + 0.01 * drand48();
            sum += B[j][i];
        }
        for (j = 0; j<SYMNUM; j++)
            B[j][i] /= sum;
    }
    for (i = 0; i < N; i++) {
        I[i]  /= sumI;
    }
}


void Model::UniformInit(int *sym, int seqLen)
{

    int i,j;
    double sum,sumI;
    sumI = 0;

    for (i = 0; i < N; i++) {

        /* initialize the I(state_t) vector */
        I[i] = 1;
        sumI += I[i];
        /* initialize the A(state_t, state_t+1) matrix */
        sum = 0.0;
        for (j = 0; j < N; j++) {
            A[i][j] = 1;
            sum += A[i][j];
        }
        for (j = 0; j < N; j++)
            A[i][j] /= sum;

        /* initialize the B(output,state) matrix */
        sum = 0.0;

        for(j=0;j<seqLen; j++) {

            B[sym[j]][i] = 1;
        }

        for (j = 0; j<SYMNUM; j++) {
            sum += B[j][i];
        }
        for (j = 0; j<SYMNUM; j++)
            B[j][i] /= sum;
    }
    for (i = 0; i < N; i++) {
        I[i]  /= sumI;
    }
}


void Model::Save(char *modelFile)
{
    ofstream ofs(modelFile);
    ofs << N << endl;
    int i;
    ofs << "InitPr "<< N << endl;
    for (i=0; i<N; i++) {
        ofs << i << " "<< I[i] << endl;
    }

    int count;
    count = 0;
    int j;
    for (i=0; i<N; i++) {
        for (j=0; j<SYMNUM; j++) {
            if (B[j][i]>0) {
                count++;
            }
        }
    }
    ofs << "OutputPr "<< count << endl;
    for (i=0; i<N; i++) {
        for (j=0; j<SYMNUM; j++) {
            if (B[j][i]>0) {
                ofs << i << " " << (char)(j+baseChar) << " "<< B[j][i] << endl;
            }
        }
    }

    count = 0;
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            if (A[i][j]>0) {
                count++;
            }
        }
    }
    ofs << "TransPr "<< count << endl;
    for (i=0; i<N; i++) {
        for (j=0; j<N; j++) {
            if (A[i][j]>0) {
                ofs << i << " " << j << " "<< A[i][j] << endl;
            }
        }
    }

    ofs.close();
}

void Model::Decode(char *seqFile) {
    vector<TrellisNode> trellis;  // Current trellis information
    trellis.resize(N);

    int i;
    ifstream seq(seqFile);
    if (seq.fail()) {
        cerr << "can't open sequence file : "<< seqFile << endl;
        exit(1);
    }

    char c;
    seq >> c;
    for (i=0; i<N; i++) {
        // probability of starting at state i: p(i) * p(c|s_i)
        // working on the logorithm domain to avoid overflow.
        trellis[i].pr = log(I[i]) + log(B[c-baseChar][i]);
        if(!isnormal(trellis[i].pr)){
            trellis[i].pr = LZERO;
        }
    }

    while (seq >> c) {
	int cIndex = c- baseChar;
    // use cIndex to access matrix B, e.g., B[cIndex][i] gives you
    // the probability of generating symbol c from state i

        double temp[N];
        for (i=0; i<N; i++) {
            temp[i] = 0;
        }

        for(i=0; i<N; i++) {
    		for(int j=0; j<N; j++) {
               	     temp[j] = trellis[j].pr + log(A[j][i]) + log(B[cIndex][i]);
    		}
    	
                
            int max_position = 0;
            for(int k=1; k<N; k++) {
    			if(temp[k] > temp[max_position]) {
                       		max_position = k;
    			}
            }
	
	    trellis[i].newPr = temp[max_position];
        trellis[i].path.push_back(max_position);

        }

        for (i=0; i<N; i++) {
    		trellis[i].pr = trellis[i].newPr;
        }

	

    // =============== insert your code here ========
    // You need to insert the appropriate code to grow
    // the best path for each state according to the Viterbi algorithm
    // Suppose the sequence already processed is <c1, c2, ..., ck>.
    // For the next character c, you generally need to update each of
    // the trellis node i.e., trellis[i], so that it will hold
    // the state transition path ending at state i that most likley
    // "generates" the
    // sequence <c1, ..., ck, c>. Remember trellis[i].path was supposed
    // to hold the most likely path for <c1, ..., ck> that ends at state i.
    // Note that the path does not need to store the state of generating c 
    // since trellis[i] implies the path ends at state i.
    // trellis[i].pr was supposed to hold the probability of data given
    // the path.

    // trellis[i].path is implemented as a STL vector which is very similar to an array.
    // You may use the push_back method to insert a value at the end of the array.
    // You may also use the resize function adjust the size of the array.
    //
    // You may also find the newPr field useful for temporarily storing the
    // probability for the new best path you constructed. (Other states may
    // still need the probability for the old path (stored in the field pr)
    // in order to compute the new probabilities for their new best paths.
    // Remember to work on the logarithm domain to avoid overflow.

    }

    seq.clear();
    seq.seekg(ios::beg);
    vector<int> bestPath;
    int bestTrellis;
    double bestTrellisPr = 0;
    for (int i = 0; i<N; i++ ){
        if (i==0 || (trellis[i].pr > bestTrellisPr)) {
            bestTrellis =i;
            bestTrellisPr = trellis[i].pr;
        }
    }
    int t = trellis[bestTrellis].path.size() - 1;
    int idx = bestTrellis;

    while(t>=0){
      bestPath.push_back(idx);
      idx = trellis[idx].path[t--];
    }
    bestPath.push_back(idx);

    for(int i=bestPath.size()-1;i>=0;i--){
      seq >> c;
      cout << c << " "<< bestPath[i] << endl;
    }
}


void Model::ComputeAlpha(int *seq, int seqLength)
{
  // seq has the sequence of observed symbols, and is an array of
  // length seqLength. It is 0-indexed. seq[0] is the
  // first observed symbol, and seq[seqLength-1] the last.
  // You can think of seqLength as the T in our formulas.
  //
  // seq[t] gives you the *index* of the symbol, which can be used
  // directly as a subscript to access matrix B, e.g., B[seq[t]][i]
  //
  // The alpha array is a seqLength x N 2-dimensional array.
  // alpha[t][i] = prob. of generating observations up to time t
  // and being in state i at time t
  // t goes from 0 to seqLength-1.
  //
  // ======= insert your code for computing the normalized alpha values
  // Store the normalized alpha values in the alpha array and the normalizers
  // in the eta array.
  // The eta array is an array of length seqLength. eta[t] is the normalizer
  // at time t. t goes from 0 to seqLength-1.

    for (int i = 0; i < seqLength; i++) {
    	eta[i] = 0;
    	for (int j = 0; j < N; j++) {
    		alpha[i][j] = 0;
    	}
    }

    for (int i = 0; i < N; i++) {
	    eta[0] += I[i] * B[seq[0]][i];
    }
    eta[0] = 1 / eta[0];


    for (int i = 0; i < N; i++) {
        // probability of starting at state i
        alpha[0][i] = I[i] * B[seq[0]][i] * eta[0];
    }
    

    for (int i = 1; i < seqLength; i++) {
        for (int k = 0; k < N; k++) {
	       double sum = 0;
	       for (int t = 0; t < N; t++) {
		      sum += alpha[i-1][t] * A[t][k];
	       }
	       eta[i] += B[seq[i]][k] * sum;
	    }
        eta[i] = 1 / eta[i];

	    for (int j = 0; j < N; j++) { 
	       for (int k = 0; k < N; k++) {
                alpha[i][j] += alpha[i-1][k] * A[k][j] * B[seq[i]][j] * eta[i];
            }  
        }
    }

}

void Model::ComputeBeta(int *seq, int seqLength)
{
  // seq has the sequence of observed symbols, and is an array of
  // length seqLength. It is 0-indexed. seq[0] is the
  // first observed symbol, and seq[seqLength-1] the last.
  // You can think of seqLength as the T in our formulas.
  //
  // seq[t] gives you the *index* of the symbol, which can be used
  // directly as a subscript to access matrix B, e.g., B[seq[t]][i]

  // The beta array is a seqLength x N 2-dimensional array.
  // beta[t][i] = prob. of generating observations after time t given
  // being in state i at time t
  // t goes from 0 to seqLength-1.


  //======= insert your code for computing the normalized beta values
  // Store the normalized beta values in the beta array.

    for (int i = 0; i < seqLength; i++) {
        for (int j = 0; j < N; j++) {
        	beta[i][j] = 0;
        }
    }	



    for (int i = 0; i < N; i++) {
        // probability of ending at state i is even
        beta[seqLength-1][i] = 1;
    }

    
    
    for (int i = seqLength-2; i >= 0; i--) {
        for (int j = 0; j < N; j++) {
            for (int k = 0; k < N; k++) {
                beta[i][j] += beta[i+1][k] * A[j][k] * B[seq[i+1]][k] * eta[i+1];
            }  
        }
    }

}


void Model::AccumulateCounts(int *seq, int seqLength)
{
    vector<vector<double> > gamma;;

    gamma.resize(seqLength);
    for(int t=0;t<seqLength;t++) {
        gamma[t].resize(N);
    }

    for (int i = 0; i < seqLength; i++) {
        for (int j = 0; j < N; j++) {
     	   gamma[i][j] = 0;
        }
    }

  // ========= Insert your code here for computing the gamma's based on
  // the alpha's and beta's
    for (int i = 0; i < seqLength; i++) {
        double sum = 0;
        for (int j = 0; j < N; j++) {
            sum += alpha[i][j] * beta[i][j];
        }        

        for (int k = 0; k < N; k++) {
            gamma[i][k] = alpha[i][k] * beta[i][k] / sum;
        }
    }



  // ========= Insert your code here for counting for initial state distribution
  // use ICounter and INorm.
    INorm = 0;
    for (int i = 0; i < N; i++) {
        ICounter[i] = gamma[0][i];
        INorm += ICounter[i];
    }



  // ========= Insert your code here for counting for output probabilities
  // use BCounter and BNorm.

    for (int i = 0; i < SYMNUM; i++) {
    	for (int j = 0; j < N; j++) {
    		BCounter[i][j] = 0;
    		BNorm[j] = 0; 
    	}
    }


    for (int k = 0; k < SYMNUM; k++) {
        for (int i = 0; i < N; i++) {
            for (int t = 0; t < seqLength; t++) {
                if (k == seq[t]) {
                    BCounter[k][i] += gamma[t][i];
                }
            }
        }
    }

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < seqLength; j++) {
            BNorm[i] += gamma[j][i];
        }
        
    }




  // ========= Insert your code here for computing the xi's and counting
  //         for state transition probabilities
  //
  // use ACounter and ANorm.
  // You may find the formula (6) in Section 4 is the easiest to use.



  // Note that you do *not* need to pre-compute all the xi_t(i,j)'s before
  // updating the counters. Instead, you can consider a loop like
  // the following and compute the needed xi_t(i,j) on the fly.
  //

    // compute the needed xi_t(i,j) and
    // update ACounter and ANorm

   for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            ACounter[i][j] = 0;
        }
	    ANorm[i] = 0;
   }



    for (int t = 0; t < seqLength-1; t++) {
        // double sum = 0;
        // for (int i = 0; i < N; i++) {
        //     sum += alpha[t][i] * beta[t][i];
        // }        

        // for(int i = 0; i < N; i++) {
        //     for(int j = 0; j < N; j++) {
        //         ACounter[i][j] += alpha[t][i] * A[i][j] * B[seq[t+1]][j] * eta[t+1] * beta[t+1][j] / sum;
        //     } 
        // }

        for(int i = 0; i < N; i++) {
            for(int j = 0; j < N; j++) {
                ACounter[i][j] += gamma[t][i] * A[i][j] * B[seq[t+1]][j] * eta[t+1] * beta[t+1][j] / beta[t][i];
            } 
        }
    }

    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++) {
            ANorm[i] += ACounter[i][j];
        }
    }

  // ========================================

}



void Model::PrintMatrix(vector<vector<double> >& mat, int seqLen)
{
    int i, j;
    cout << "   ";
    for (i=0; i<seqLen; i++ ){
        cout << " "<< i;
    }
    cout << endl;
    for (i=0; i<N; i++) {
        cout << i << " ";
        for (j=0; j<seqLen; j++) {
            cout << mat[j][i] << " ";
        }
        cout << endl;
    }
}


void Model::PrintOutputMatrix()
{
    int i, j;
    cout << "   ";
    for (i=0; i<N; i++ ){
        cout << " "<< i;
    }
    cout << endl;
    for (i=0; i<SYMNUM; i++) {
        bool valid = false;
        for (j=0; j<N;j++) {
            if (B[i][j]>0) {
                valid=true;
            }
        }
        if (valid) {
            char c=i+baseChar;
            cout << c << " ";
            for (j=0; j<N; j++) {
                cout << B[i][j] << " ";
            }
            cout << endl;
        }
    }
}


double Model::UpdateHMM(int *data, int seqLength)
{
    ResetCounter();
    alpha.resize(seqLength);
    beta.resize(seqLength);

    int i;
    for (i=0; i< seqLength; i++ ) {
        alpha[i].resize(N);
        beta[i].resize(N);
    }
    eta.resize(seqLength);

    // cout << "Transition matrix\n";
    // PrintMatrix(A, N);
    // cout << "Output matrix\n";
    // PrintOutputMatrix();

    ComputeAlpha(data, seqLength);
    // cout << "Alpha\n";
    // PrintMatrix(alpha, seqLength);

    ComputeBeta(data, seqLength);
    // cout << "Beta\n";
    // PrintMatrix(beta, seqLength);

    // compute data likelihood
    double prData =0;
    int t;
    for (t=0; t<seqLength; t++) {
        prData += -log(eta[t]);
    }

    ResetCounter();  // Initialize all the counters and the normalizers

    AccumulateCounts(data, seqLength);

    UpdateParameter();


    return (prData);
}


void Model::Train(char *seqFile)
{
    vector<char> buffer(0);
    ifstream ifs(seqFile);
    if (ifs.fail()) {
        cerr << "can't open sequence file for training: "<< seqFile << endl;
    }
    char c;
    while (ifs >> c) {
        buffer.push_back(c);
    }

    /*
       We train for many epochs. As we train, we keep track of how well
       the current HMM models the data ( P(data | hmm) ). Additionally,
       we keep track of a history of how the model has fit the data in
       previous epochs (meanFit). As the HMM settles into a stable
       state, currentFit will asymptote to a particular value. The
       time-averaged meanFit will asymptote to this value also (though
       more slowly). When the currentFit becomes very similar to the
       meanFit, the model isn't changing, so we stop training.
     */

    int maxT = buffer.size();
    int *sym = new int[maxT];
    vector<char>::iterator it;
    int k=0;
    for (it=buffer.begin(); it!=buffer.end(); it++) {
        sym[k++]= (*it)- baseChar;
        //    cerr << "sym: "<< sym[k-1]<<endl;
    }

    // UniformInit(sym,maxT);
    RandomInit(sym,maxT);

    double currentFit, meanFit = -1e-80;
    for(int epoch = 0; fabs( (meanFit - currentFit) / meanFit ) > 1e-6; epoch++) {
        /* Train for many epochs */
        cerr << "\nEpoch " <<  epoch << endl;
        currentFit = UpdateHMM(sym, maxT);
        cerr << " log-likelihood = "<< currentFit << endl;
        if (fabs(currentFit)< 0.00001) {
            break;
        }
        meanFit = 0.2 * currentFit + 0.8 * meanFit;
    }

}
