
//g++ matchjets.cxx -o short-example `/opt/local/bin/fastjet-config --cxxflags --libs --plugins`

#include "fastjet/PseudoJet.hh"
#include "fastjet/ClusterSequence.hh"
#include "fastjet/ClusterSequenceArea.hh"
#include "fastjet/Selector.hh"
#include "fastjet/tools/Filter.hh"
#include "fastjet/tools/Pruner.hh"

#include <iomanip>

#include <chrono>

using namespace fastjet;
using namespace std;

using std::chrono::high_resolution_clock;
using std::chrono::duration_cast;
using std::chrono::nanoseconds;

const static double R = 0.4;

//This function takes a vector of jets to be geometrically matched to another vector of "candidate" matches. Once matching occurs, a vector of indices is returned allowing one to index the original two vectors to fill responses etc. with the matched jet pairs. Redundantly (for debugging, etc.) the vectors of matches themselves are also updated for later use.
//Note: We need to be able to remove jets from the "candidates" vector after they've been matched, so we make a copy in the function. Also make a copy candidates vector for each iteration on toMatch since this vector has selections applied to it
//Note: In finding which jets were the matches, we know the toMatch jet match will be the 'i'th jet since we are iterating. The candidate_copy jet should be the highest pT match, so the first one in the candidate_copy list. Geometrically match the candidate_copy jet to the nearest candidate jet, since jets have been removed so they don't index to the same jet anymore
std::vector<int> MatchJets(const std::vector<fastjet::PseudoJet> candidates_safe, const std::vector<fastjet::PseudoJet> toMatch, std::vector<fastjet::PseudoJet> & c_matches, std::vector<fastjet::PseudoJet> & t_matches) {
    std::vector<int> match_indices;
    if (candidates_safe.size() == 0 || toMatch.size() == 0) {
        return match_indices; //later, match_indices being empty will tell us there were no matches
    }
    //define candidates outside the loop so list continually dwindles as we remove matched candidates
    std::vector<fastjet::PseudoJet> candidates = candidates_safe;
    for (int i = 0; i < toMatch.size(); ++ i) { //for each jet in toMatch, we try to find a match from a dwindling set of candidates
        fastjet::Selector selectMatchedJets = fastjet::SelectorCircle( R );
        selectMatchedJets.set_reference( toMatch[i] );
        std::vector<fastjet::PseudoJet> matchedToJet = /*sorted_by_pt(*/ selectMatchedJets( candidates )/*)*/; //not sorting by pT because we sort by deltaR later
        if (matchedToJet.size() == 0) { continue; } //means no match to this jet. Remove none from candidates. Continuing on to the next one.
        else { //found at least one match. Need to remove the highest pT one from candidates and add the respective jets to the match vectors.
            //sorting by deltaR vv
            vector<double> dists;
            for (int j = 0; j < matchedToJet.size(); ++ j) {
                dists.push_back(matchedToJet[j].delta_R(toMatch[i]));//fills vector with values of deltaR between matched jets and their match
            }//loop over matched jets to get deltaR values to the reference jet
            matchedToJet = objects_sorted_by_values(matchedToJet, dists);//sorts by deltaR to toMatch[i].
            //sorted by deltaR ^^
      
            match_indices.push_back(i); //push back the toMatch match position
            t_matches.push_back(toMatch[i]);
            c_matches.push_back(matchedToJet[0]); //highest pT match
            for (int j = 0; j < candidates_safe.size(); ++ j) {//finding the match index (needs to be separate from next loop because we're removing jets from candidates but not candidates_safe.)
                if (matchedToJet[0].delta_R(candidates_safe[j]) == 0) {
                    match_indices.push_back(j);
                    break; //should exit only the c_matches loop
                }//found the matched jet in the unsorted non-shrinking list of candidates to add a match index to the list of indices
            }//loop to find the matched jet in the unsorted non-shrinking corresponding list of candidates
            for (int j = 0; j < candidates.size(); ++ j) { //finding which one to delete from candidates before next toMatch iteration.
                if (matchedToJet[0].delta_R(candidates[j]) == 0) { //is probably the same jet
                    candidates.erase(candidates.begin() + j); //removing the jet from the overall list of candidates so it can't be considered next time
                    break; //should exit only the c_matches loop.
                }//found the matched jet in the unsorted shrinking list of candidates and removed it
            }//loop to find the matched jet in the unsorted shrinking corresponding list of candidates
        }//found a match
    }//loop over jets to be matched to candidates
    return match_indices;
}
  
//!this function is similar to the "MatchJets" function, but it instead takes a list of jets which have already been matched ("candidates_safe") and finds the jets to which they correspond in the list of unmatched jets ("toMatch") by (exact) geometrical matching. The remaining jets with no corresponding already-matched jets are either misses or fakes depending on the call.
std::vector<int> FakesandMisses(const std::vector<fastjet::PseudoJet> candidates_safe, const std::vector<fastjet::PseudoJet> toMatch, std::vector<fastjet::PseudoJet> & unmatched) {
    std::vector<int> miss_fake_index;
    std::vector<fastjet::PseudoJet> candidates = candidates_safe;
    for (int i = 0; i < toMatch.size(); ++ i) {
        fastjet::Selector selectMatchedJets = fastjet::SelectorCircle( 0 ); //a "match" is now if we found the same exact jet.
        selectMatchedJets.set_reference( toMatch[i] );
        std::vector<fastjet::PseudoJet> matchedToJet = /*sorted_by_pt(*/ selectMatchedJets( candidates )/*)*/;//not sorting by either pT or deltaR because we should get at most one exact match
        if (matchedToJet.size() == 0) { //means no match to this jet. Remove none from candidates. Add it to unmatched & continue to next one.
            miss_fake_index.push_back(i); //the ith jet in toMatch is a miss or a fake
            unmatched.push_back(toMatch[i]);
            continue;
        }//didn't find a match
        else { //found an exact match
            for (int j = 0; j < candidates.size(); ++ j) { //finding which one to delete from candidates before next toMatch iteration.
                if (matchedToJet[0].delta_R(candidates[j]) == 0) { //is the same jet
                    candidates.erase(candidates.begin() + j); //removing the jet from the overall list of candidates so it can't be considered next time
                    break; //should exit only the candidates loop.
                }//found the match in the corresponding unsorted vector
            }//loop to find the match to remove from consideration in the next go-around
        }//found a match
    }//loop over toMatches
    return miss_fake_index;
}

//tests
int main () {

    const int nJets_ref = 5;
    const int nJets = 3;
    const int nDim = 3;
    //getting some random numbers for tests
    double vals_ref[nJets_ref][nDim];
    double vals[nJets][nDim];
    for (int i = 0; i < nJets_ref; ++ i) {
        for (int j = 0; j < nDim; ++ j) {
            auto start = high_resolution_clock::now();
            auto end = high_resolution_clock::now();
            vals_ref[i][j] = duration_cast<nanoseconds>(end - start).count() / (double) (fmod(rand(), 100.00)); //rand() % 5;
            cout << std::setprecision(5) << (double) vals_ref[i][j] << endl;
        }
    }
    for (int i = 0; i < nJets; ++ i) {
        for (int j = 0; j < nDim; ++ j) {
            auto start = high_resolution_clock::now();
            auto end = high_resolution_clock::now();
            vals[i][j] = duration_cast<nanoseconds>(end - start).count() / (double) (fmod(rand(), 100.00)); //rand() % 5;
            cout << std::setprecision(5) << (double) vals[i][j] << endl;
        }
    }
    
    cout << "~~~NEW TEST~~~" << endl;
    
    //to be matched
    vector<PseudoJet> input1_beforesort;
    vector<PseudoJet> input2_beforesort;
    
    //px,py,pz,E
    for (int i = 0; i < nJets_ref; ++ i) {
        input1_beforesort.push_back(PseudoJet(vals_ref[i][0],vals_ref[i][1],vals_ref[i][2],sqrt(pow(vals_ref[i][0],2) +pow(vals_ref[i][1],2) +pow(vals_ref[i][2],2))));
    }
    for (int i = 0; i < nJets; ++ i) {
        input2_beforesort.push_back(PseudoJet(vals[i][0],vals[i][1],vals[i][2],sqrt(pow(vals[i][0],2) +pow(vals[i][1],2) +pow(vals[i][2],2))));
    }
    
    vector<PseudoJet> input1 = sorted_by_pt(input1_beforesort);
    vector<PseudoJet> input2 = sorted_by_pt(input2_beforesort);
    
    for (int i = 0; i < input1.size(); ++ i) {
        cout << "input1 #" << i << ": pt="<< input1[i].pt() << " rap="
                           << input1[i].rap() << " phi=" << input1[i].phi() << endl;
    }
    
    for (int i = 0; i < input2.size(); ++ i) {
        cout << "input2 #" << i << ": pt="<< input2[i].pt() << " rap="
                           << input2[i].rap() << " phi=" << input2[i].phi() << endl;
    }
    for (int i = 0; i < input1.size(); ++ i) {
        for (int j = 0; j < input2.size(); ++ j) {
        cout << "dist (" << i << ", " << j << "): deltaR="<< input1[i].delta_R(input2[j]) << endl;//" deltaPhi=" << ref.delta_phi_to(a) << endl;
        }
    }
    //final matches
    std::vector<fastjet::PseudoJet> output1;
    std::vector<fastjet::PseudoJet> output2;
    
    std::vector<int> indices;
    std::vector<int> miss_indices;
    std::vector<int> fake_indices;
    
    indices = MatchJets(input2,input1,output2,output1);
    cout << "these jet pairs were matched:" << endl;
    for (int i = 0; i < indices.size(); i+=2) {
         cout << "(" << indices[i] << ", " << indices[i+1] << ")" << endl;
    }
    
    for (int i = 0; i < output1.size(); ++ i) {//output1.size() = output2.size()
        
        cout << "output1 #" << i << ": pt="<< output1[i].pt() << " rap="
                           << output1[i].rap() << " phi=" << output1[i].phi() << endl;
        cout << "output2 #" << i << ": pt="<< output2[i].pt() << " rap="
                           << output2[i].rap() << " phi=" << output2[i].phi() << endl;
    }
    
    vector<PseudoJet> misses; vector<PseudoJet> fakes;
    miss_indices = FakesandMisses(output1,input1,misses);
    fake_indices = FakesandMisses(output2,input2,fakes);
    
    cout << "the kicker: " << endl;
    cout << "misses at: ";
    for (int i = 0; i < miss_indices.size(); ++ i) {
        cout << miss_indices[i] << " ";
    }
    cout << " and fakes at: ";
    for (int i = 0; i < fake_indices.size(); ++ i) {
        cout << fake_indices[i] << " ";
    }
    cout << endl;
    cout << "~~~TEST COMPLETE~~~" << endl;
    return 0;
}
