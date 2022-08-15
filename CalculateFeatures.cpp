/***********************************************
 * # Copyright 2019-2020. Kuang Mengmeng
 * # GPL version 3.0 applies.
 * ************************************************/

#include <string>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <list>
#include <set>
#include <algorithm>
#include <climits>
#include <cstdio>
#include <cstdlib>
#include <cerrno>
#include <iomanip>
#include "ProbabilisticModel.h"
#include "MultiSequence.h"
#include "CalculateFeatures.h"


CalculateFeatures::CalculateFeatures(MultiSequence *sequences, const ProbabilisticModel &model)
{
    features = Calculation(sequences, model);
}
CalculateFeatures::~CalculateFeatures()
{
    features = "";
}


// string CalculateFeatures::Calculation(MultiSequence *sequences, const ProbabilisticModel &model)
// {
//     assert(sequences);
//     //feature
//     const int numSeqs = sequences->GetNumSequences();

//     const int numPairSeqs = (numSeqs - 1) * numSeqs / 2;
//     int sum_of_pairs = 0;
//     //feature
//     float pid = 0.0;
//     // do all pairwise alignments for family similarity
//     float* PIDs = new float[(numSeqs - 1) * numSeqs / 2 * sizeof(float)];
//     //feature
//     float avg_best_probs = 0.0;
//     float* BestProbs = new float[(numSeqs - 1) * numSeqs / 2 * sizeof(float)];
//     float avg_align_score = 0.0;
//     float* AlignScores = new float[(numSeqs - 1) * numSeqs / 2 * sizeof(float)];

//     float avg_emit_pairs = 0.0;
//     float* EmitPairs = new float[(numSeqs - 1) * numSeqs / 2 * sizeof(float)];
//     //feature
//     float avg_align_length = 0.0;
//     //feature
//     float avg_seq_length = 0.0;
//     for(int i = 0; i < numSeqs; i++)
//     {
//         avg_seq_length += sequences->GetSequenceLength(i);
//     }
//     avg_seq_length /= numSeqs;

//     int pairIdx = -1;
//     for (int a = 0; a < numSeqs - 1; a++)
//     {
//         for (int b = a + 1; b < numSeqs; b++)
//         {
//             pairIdx++;
//             Sequence *seq1 = sequences->GetSequence(a);
//             Sequence *seq2 = sequences->GetSequence(b);
//             pair<SafeVector<char> *, float> alignment = model.ComputeViterbiAlignment(seq1,seq2);
//             SafeVector<char>::iterator iter1 = seq1->GetDataPtr();
//             SafeVector<char>::iterator iter2 = seq2->GetDataPtr();
//             avg_best_probs += (1 + ((float)alignment.second / ( 10 * alignment.first->size())));
//             BestProbs[pairIdx] = 1 + ((float)alignment.second / ( 10 * alignment.first->size() ));
//             float N_correct_match = 0;
//             int i = 1;
//             int j = 1;
//             bool isFirst = true;
//             float tmp_alignment_score = 0.0;
//             float tmp_emit_pairs = 0.0;
//             for (SafeVector<char>::iterator iter = alignment.first->begin();
//                     iter != alignment.first->end(); ++iter)
//             {
//                 if (*iter == 'B')
//                 {
//                     unsigned char c1 = (unsigned char) iter1[i++];
//                     unsigned char c2 = (unsigned char) iter2[j++];
//                     float tmp_num = BLOSUM62[alphabet.find(c1)][alphabet.find(c2)];
//                     if (tmp_num >= -4 && tmp_num <= 11){
//                         tmp_alignment_score += tmp_num;
//                     }else{
//                         tmp_alignment_score += 0;
//                     }
//                     tmp_num =  emitPairs[alphabet.find(c1)][alphabet.find(c2)];
//                     if( tmp_num >= -10 && tmp_num <= 10){
//                         tmp_emit_pairs += tmp_num;
//                     }else{
//                         tmp_emit_pairs += 0;
//                     }
//                     isFirst = true;
//                     if(c1==c2)
//                     {
//                         N_correct_match += 1;
//                     }
//                 }
//                 else if(*iter == 'X')
//                 {
//                     i++;split
//                     {
//                         tmp_alignment_score -= 1;
//                     }
//                 }
//                 else if(*iter == 'Y')
//                 {
//                     j++;
//                     if(isFirst)
//                     {
//                         isFirst = false;
//                         tmp_alignment_score -= 7;
//                     }
//                     else
//                     {
//                         tmp_alignment_score -= 1;
//                     }
//                 }
//             }
//             EmitPairs[pairIdx] = tmp_emit_pairs / alignment.first -> size();
//             avg_emit_pairs += tmp_emit_pairs / alignment.first -> size();
//             AlignScores[pairIdx] = tmp_alignment_score / alignment.first -> size();
//             avg_align_score += tmp_alignment_score / alignment.first -> size();
//             sum_of_pairs += N_correct_match;
//             avg_align_length += (float)alignment.first->size();
//             pid += N_correct_match / alignment.first->size();
//             PIDs[pairIdx] = N_correct_match / alignment.first->size();
//             delete alignment.first;
//         }
//     }
//     pid /= (int)(numSeqs * (numSeqs - 1) / 2);
//     avg_align_length /= (int)(numSeqs * (numSeqs - 1) / 2);splitPID
//     avg_emit_pairs /= (int)(numSeqs * (numSeqs - 1) / 2);
//     //compute the variance of PID and Align Score
//     float variance = 0;
//     for (int k = 0; k < (numSeqs-1)*numSeqs/2; k++) variance += (PIDs[k]-pid)*(PIDs[k]-pid);
//     variance /= (int)(numSeqs * (numSeqs - 1) / 2);
//     float variance_pid = sqrt(variance);

//     variance = 0;
//     for (int k = 0; k < (numSeqs-1)*numSeqs/2; k++) variance += (AlignScores[k]-avg_align_score)*(AlignScores[k]-avg_align_score);
//     variance /= (int)(numSeqs * (numSeqs - 1) / 2);
//     float variance_align_score = sqrt(variance);

//     variance = 0;
//     for (int k = 0; k < (numSeqs-1)*numSeqs/2; k++) variance += (BestProbs[k]-avg_best_probs)*(BestProbs[k]-avg_best_probs);
//     variance /= (int)(numSeqs * (numSeqs - 1) / 2);
//     float variance_best_probs = sqrt(variance);

//     variance = 0;
//     for (int k = 0; k < (numSeqs-1)*numSeqs/2; k++) variance += (EmitPairs[k]-avg_emit_pairs)*(EmitPairs[k]-avg_emit_pairs);
//     variance /= (int)(numSeqs * (numSeqs - 1) / 2);
//     float variance_emit_pairs = sqrt(variance);

//     string ret;
//  //   ret += to_string(avg_seq_length) + ", ";
//     ret += to_string(numSeqs) + ", ";
//  //   ret += to_string(numPairSeqs) + ", ";
//     ret += to_string(avg_align_length) + ", ";
//  //   ret += to_string((float)avg_seq_length / avg_align_length) + ", ";
//     ret += to_string(avg_align_score) + ", ";
//  //   ret += to_string(variance_align_score) + ", ";
//  //   ret += to_string(avg_best_probs) + ", ";
// //    ret += to_string(variance_best_probs) + ", ";
//     ret += to_string(pid) + ", ";
// //    ret += to_string(variance_pid) + ", ";
//     ret += to_string(avg_emit_pairs);
// //    ret += to_string(variance_emit_pairs);
//   //  ret += to_string(sum_of_pairs);

//     return ret;

// }

string CalculateFeatures::Calculation(MultiSequence *sequences, const ProbabilisticModel &model){
	assert(sequences);
	const int numSeqs = sequences->GetNumSequences();
    float identity = 0;

	float* new_final_arr = new float[65535]();
	float* PIDs = new float[ (numSeqs - 1) * numSeqs / 2 * sizeof(float) ];
    float* SOPs = new float[ (numSeqs - 1) * numSeqs / 2 * sizeof(float) ];
    int avg_length = 0;
	string ret_list ;
	for (int a = 0; a < numSeqs - 1; a++) {
		for (int b = a + 1; b < numSeqs; b++) {
				int tmp_sp_idx = 0;
				float tmp_sp = 0;
 			int num_idx = 0;
			 string str_seq1;
			 string str_seq2;
			Sequence *seq1 = sequences->GetSequence(a);
			Sequence *seq2 = sequences->GetSequence(b);
			pair<SafeVector<char> *, float> alignment = model.ComputeViterbiAlignment(seq1,seq2);
			SafeVector<char>::iterator iter1 = seq1->GetDataPtr();
    		SafeVector<char>::iterator iter2 = seq2->GetDataPtr();
			ret_list += to_string(alignment.first -> size()) + ", ";
			string tmp_str;
            float N_correct_match = 0;
			float N_emit = 0;
            int i = 1;int j = 1;
			for (SafeVector<char>::iterator iter = alignment.first->begin();
				iter != alignment.first->end(); ++iter){
					tmp_str += *iter;
				if (*iter == 'B'){
					unsigned char c1 = (unsigned char) iter1[i++];
					unsigned char c2 = (unsigned char) iter2[j++];
   					if(c1==c2) N_correct_match += 1;
						N_emit += emitPairs[alphabet.find(c1)][alphabet.find(c2)];
						if (BLOSUM62[alphabet.find(c1)][alphabet.find(c2)] <= 11 && BLOSUM62[alphabet.find(c1)][alphabet.find(c2)] >= -4  ){
							tmp_sp += BLOSUM62[alphabet.find(c1)][alphabet.find(c2)];
						}
						tmp_sp_idx += 1;
						num_idx ++;
				}
				else if(*iter == 'X'){
                    num_idx ++;
					tmp_sp_idx += 1;
					i++;
        }
 				else if(*iter == 'Y'){
           			num_idx ++;
					 tmp_sp_idx += 1;
					j++;
        }
            }
            if(i!= seq1->GetLength()+1 || j!= seq2->GetLength() + 1 ) cerr << "percent identity error"<< endl;
			ret_list += tmp_str + ", ";
			ret_list += seq1->GetString() + "-";
			ret_list += seq2->GetString() + ", ";
            ret_list += to_string(N_correct_match / alignment.first->size()) + ", ";
			ret_list  += to_string(N_emit / alignment.first->size()) + ", ";
			ret_list += to_string(tmp_sp / tmp_sp_idx) + ";";
			delete alignment.first;
		}
	}

    return ret_list;
}
