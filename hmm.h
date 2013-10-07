#ifndef _HMM_HPP 
#define _HMM_HPP

#include <vector>

using namespace std;

// ============== Representation of HMM model =======================
// We assume a fixed vocabulary that consists of all the 94 readable
// ascii characters. That is, the symbol set is fixed as all the 94 
// visible ascii characters (from [33,!] to [126, ~]).
// When a symbol is not observed in the data, its count is automatically
// set to zero, effectively excluded from the model.
// 
// The number of state is stored in a variable N.
// 
// An N-state HMM is represented by the following arrays:
//     1. Initial state probability (I)
//         I is an array of length N (0-indexed), with I[s] 
//         representing the initial probabiltiy of being at state s.
//     2. State transition probability matrix (A)
//         A is a 0-indexed two-dimensional array (NxN), with A[i][j]
//         representing the probability of going to state j from state i.
//     3. Output probability matrix (B)
//         B is a 0-indexed two-dimensional array (MxN), with A[o][s]
//         representing the probability of generating output symbol "o"
//         at state "s". M is the number of unique symbols (currently 94,
//         stored in SYMNUM).
//       
// Note the slight difference between the N used in this program and
// the N used in the assignment text and the course slides.
// There, we have N+1 states, but here, we have N states, so the
// subscript goes from 0 up to N-1, not N!

// HMM Representation on the disk file
// An HMM is stored in a file with the following example format:
// 
// 2
// InitPr 2
// 0 0.5
// 1 0.5
// OutputPr 3
// 0 a 0.5
// 0 b 0.5
// 1 a 1
// TransPr 4
// 0 0 0.5
// 0 1 0.5
// 1 0 0.1
// 1 1 0.9
//
//
// There could be any number of spaces between the numbers or words, but
// the order of keywords and numbers is important. The spelling of the
// keywords "InitPr", "OutputPr", and "TransPr" is strict.
//  -  The first number is the number of states
//  -  The number after each keyword is the number of entries following it
//     that should be interpreted with this keyword
//  -  The entries not specified (e.g., p(b|state_1)) imply that the values
//     are zero.
//  The three keywords themselves suggest what they mean in a straightforward way.
//  

struct TrellisNode {
  vector<int> path; //it stores the best path, in a sequence of state indices, before reaching the current node
  double pr; // the probability of the "path"
  double newPr; // a buffer to store the probability for the new extended path
};

const double LZERO = -1000000000;

class Model {
public:
  Model(int stateCount);
  Model(char *modelFile);

  ~Model();

  void CountSequence(char *seqFile); // estimate an HMM with tagged sequence
  void Save(char *modelFile); // save the model to disk
  void Decode(char *seqFile); // decode a sequence with the Viterbi algorithm
  void Train(char *seqFile); // train an HMM with an untagged sequence using Baum-Welch algorithm
  void ResetCounter(); // initialize counters and the corresponding normalizers
  void UpdateParameter(); // update the parameters based on the counters and the normalizers.



private: 

  void PrintMatrix(vector<vector<double> >& mat, int seqLen);
  void AllocateMemory();
  void RandomInit(int *sym, int seqLen);
  void UniformInit(int *sym, int seqLen);

  void PrintOutputMatrix();


  void ComputeAlpha(int *seq, int seqLength);
  void ComputeBeta(int *seq, int seqLength);
  void AccumulateCounts(int *seq, int seqLength);

  double UpdateHMM(int *data, int seqLength);


  static const char baseChar = '!'; // to shift the index of symbols to
  // start from 0-index.


  //========= these are the parameters of the HMM ===============

  int N; // number of states
  static const int SYMNUM = 94; // number of symbols

  vector<double> I; // initial state probability
  vector<vector<double> > A; // state transition probability
  vector<vector<double> > B; // output probability


  // =========== these are counters for the parameters ===========
  // The counters are useful when estimating the parameters

  vector<double> ICounter; // counter for initial state probability
  double INorm;  // normalizer for initial state probability
  vector<vector<double> > ACounter; // counter for state transition probability
  vector<double> ANorm; // normalizer for transition probability
  vector<vector<double> > BCounter; // counter output probability
  vector<double> BNorm; // normalizer for output probability

  vector<vector<double> > alpha; // alpha array, forward probabilities, 0-indexed
  // alpha[t][i] = prob. of generating observations up to time t and being in state i at time t. If the observations are o_1 ... o_T, then t goes from 0 to T-1.

  vector<vector<double> > beta; // beta array, backward probabilities, 0-indexed
  // beta[t][i] = prob. of generating observations after time t given
  // being in state i at time t
  // If the observations are o_1 ... o_T, then t goes from 0 to T-1.

  vector<double> eta; // alpha normalizer. 0-indexed
  // eta[t] is the normalizer for time t.
  // If the observations are o_1 ... o_T, then t goes from 0 to T-1.
 


};


#endif /* _HMM_HPP */
