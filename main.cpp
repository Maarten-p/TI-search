#include <cinttypes>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <future>
#include <map>
#include <tuple>
#include "BooleanFunction.h" 

//The number of threads used
constexpr std::size_t MAX_THREADS = 8;

//The amount of bits in the unshared input/output
//Currently does not support a different number of input and output bits
//The program should find a solution in about 5 hours for 5 bits on a good computer
constexpr std::size_t INPUT_BITS = 5;
constexpr std::size_t OUTPUT_BITS = 5;

//The program currently does NOT work for a different amount of shares
constexpr std::size_t INPUT_SHARES = 3;
constexpr std::size_t OUTPUT_SHARES = 3;

//The amount of bits in the shared input/output
constexpr std::size_t INPUT_SHARES_BITS = INPUT_BITS * INPUT_SHARES;
constexpr std::size_t OUTPUT_SHARES_BITS = OUTPUT_BITS * OUTPUT_SHARES;

//the number of combinations possible of different unshared inputs/outputs
constexpr std::size_t INPUT_SIZE = std::pow(2, INPUT_BITS);
constexpr std::size_t OUTPUT_SIZE = std::pow(2, OUTPUT_BITS);

//the number of combinations possible of different shared inputs/outputs
constexpr std::size_t INPUT_SHARES_SIZE = std::pow(2, INPUT_SHARES_BITS);
constexpr std::size_t OUTPUT_SHARES_SIZE = std::pow(2, OUTPUT_SHARES_BITS);

//the truth table for one share of an output function
typedef BooleanFunction<INPUT_SHARES_SIZE> BlnFunction;

//The correction terms for one shared output function
typedef std::array<std::bitset<INPUT_SHARES_BITS>, OUTPUT_SHARES> CorrectionFunction;

//The correction terms for all shared outputs
typedef std::vector<CorrectionFunction> VecCorrectionFunction;

//One possible combination of unshared inputs
typedef std::bitset<INPUT_BITS> InputBitArray;

//One possible combination of shared inputs
typedef std::bitset<INPUT_SHARES_BITS> SharedInputBitArray;

//Contains the shared truth table for every shared correction term that is used
//The index of every unused shared correction terms stays at 0
//It's faster to have this table precomputed than to compute each element on the fly when it's needed
std::vector<std::bitset<OUTPUT_SHARES_SIZE>> globalTruthTable(INPUT_SHARES_SIZE,0);

//Contains all valid sharings for every possible input function
//It's faster to have this table precomputed than to compute each element on the fly when it's needed
std::vector<std::vector<std::size_t>> globalValidSharesTable;

//Contains a bitset for each share of all output functions that measures on which input shares the function depends
// 110010 -> f1 depends on X1,X2,Y2
//Currently, each output function is so designed that its first share will only depend on the second and third shares of each input
//The second share only depends on the first and third share of each input, and the third share only depends on the first and second share
//This is used to satisfy the non-completeness constraint
std::vector<std::bitset<INPUT_SHARES_BITS>> globalDependence;

//The maximum hamming weight that each individual share of a correction term can have
//This value is slowly raised in an iterative deepening matter
int MAX_HAMMING_WEIGHT = 2;

//The maximum hamming weight that the correction terms of all shares of an output together can have
//This value is slowly raised in an iterative deepening matter
//This limits the search space and gives solutions with a lower hamming weight
//Which is less costly to implement on a chip
int MAX_HAMMING_WEIGHT_TOTAL = 2;

//A struct for the indices of possible combinations of correction terms for different outputs
struct Indices {
    std::vector<std::size_t> indices;

    //Removes one indice from the struct
    Indices without(std::size_t i) const {
        Indices ind;
        ind.indices = indices;
        ind.indices.erase(ind.indices.begin() + i);
        return ind;
    }
};
//redefines the << operator for a cleaner output
std::ostream& operator<<(std::ostream& os, const Indices& ind)
{
    os << '(';
    for(std::size_t i : ind.indices)
        os << i << ' ';
    os << ')';
}


/**
 * Defines a proper comparison function
 * @return true of lhs is smaller than rhs, false otherwise
 */
bool operator<(const Indices& lhs, const Indices& rhs) {
    if(lhs.indices.size() < rhs.indices.size())
        return rhs < lhs; 

    for(std::size_t i = 0; i < lhs.indices.size(); ++i) {
        if(lhs.indices[i] > rhs.indices[i])
            return false;
        else if(lhs.indices[i] < rhs.indices[i])
            return true;
    }
    return false;
}

/**
 * @return true if every (share)th share of every input is zero and if every 
 * (share)th share of the function for which the correction term is used is zero, false otherwise
 */
bool areSharesZero(const std::bitset<INPUT_SHARES_BITS>& bits, std::size_t share, int outputIndex) {
    for(std::size_t i = share; i < bits.size(); i += INPUT_SHARES) {
        if(bits[i]) 
            return false;
        if(globalDependence[outputIndex][i])
            return false;
    }
    return true;
}

//check the MAX_HAMMING_WEIGHT(_TOTAL) comments for more information
bool hammingWeightConstraint(const std::bitset<INPUT_SHARES_BITS> bits) {
    return bits.count() <= MAX_HAMMING_WEIGHT;
}

bool totalHammingWeightConstraint(const std::bitset<INPUT_SHARES_BITS> bits1, const std::bitset<INPUT_SHARES_BITS> bits2, const std::bitset<INPUT_SHARES_BITS> bits3) {
    return (bits1.count()+bits2.count()+bits3.count()) <= MAX_HAMMING_WEIGHT_TOTAL;
}

/**
 * Computes the truth table for the given correction term and add it to the global lookup table
 * @param correction_term
 */
void addToTruthTable(std::bitset<INPUT_SHARES_BITS> correction_term) {
    BlnFunction::BitArray truthTable;
    for (std::size_t j = 0;j<INPUT_SHARES_SIZE;j++) {
        std::bitset<INPUT_SHARES_BITS> bits(j);
        bool result = 0;
        std::bitset<INPUT_SHARES_BITS> im_result = bits & correction_term;
        for (std::size_t k = 0; k<INPUT_SHARES_BITS;k++) {
            result ^= im_result[k];
        }
        truthTable[j] = result;
    }
    globalTruthTable[(int) correction_term.to_ulong()] = truthTable;
}
/**
 * Finds linear correction terms for the sharing (f1, f2, f3).
 * Also adds the correction terms to the global lookup table
 * @note f1, f2 and f3 are assumed to be non-reduced, but nevertheless
 *  noncomplete.
 * @note each VecCorrectionFunction is of size one
 * @warning this doesn't work for more than 8 * sizeof(int) input variables.
 */
std::vector<VecCorrectionFunction> getLinearCorrections(const BlnFunction& f1,
    const BlnFunction& f2, const BlnFunction& f3, int outputIndex) {

    auto futureSpectrum1 = std::async(
        std::launch::async, [&f1] { return f1.walshHadamardTransform(); }
    );
    auto futureSpectrum2 = std::async(
        std::launch::async, [&f2] { return f2.walshHadamardTransform(); }
    );
    auto futureSpectrum3 = std::async(
        std::launch::async, [&f3] { return f3.walshHadamardTransform(); }
    );
    // Collect results
    auto spectrum1 = futureSpectrum1.get(); 
    auto spectrum2 = futureSpectrum2.get(); 
    auto spectrum3 = futureSpectrum3.get(); 

    std::vector<VecCorrectionFunction> solutions;

    for(std::size_t i = 0; i < INPUT_SHARES_SIZE; ++i) {
        std::bitset<INPUT_SHARES_BITS> bits_i(i);
        if(!areSharesZero(bits_i, 0, outputIndex*3) || spectrum1[i] != 0 || !hammingWeightConstraint(bits_i))
            continue;
            
        for(std::size_t j = 0; j < INPUT_SHARES_SIZE; ++j) {
            std::bitset<INPUT_SHARES_BITS> bits_j(j);
            if(!areSharesZero(bits_j, 1, outputIndex*3+1) || spectrum2[j] != 0 || !hammingWeightConstraint(bits_j))
                continue;
            std::bitset<INPUT_SHARES_BITS> bits_l = bits_i ^ bits_j;
            std::size_t l = bits_l.to_ulong();
            if(!areSharesZero(bits_l, 2, outputIndex*3+2) || spectrum3[l] != 0 || !hammingWeightConstraint(bits_l) || !totalHammingWeightConstraint(bits_i,bits_j,bits_l))
                continue;
            CorrectionFunction cf {
                bits_i,
                bits_j,
                bits_l
            };

            // Found a solution, add to list
            solutions.push_back({cf});
            addToTruthTable(bits_i);
            addToTruthTable(bits_j);
            addToTruthTable(bits_l);
        }
    }

    return solutions;
}

/**
 * Finds the collection of possible correction functions for each output function
 * @param realization The truthTable for all shares of all output functions
 */
std::map<Indices, std::vector<VecCorrectionFunction>> buildCorrectionTerms(
    std::vector<std::vector<BlnFunction>>& realization)
{
    std::cout << "Collecting linear correction terms." << std::endl;
    std::map<Indices, std::vector<VecCorrectionFunction>> correctionTerms;
    for(std::size_t i = 0; i < realization.size(); ++i) {
        Indices ind;
        ind.indices.push_back(i);
        correctionTerms[ind] = getLinearCorrections(
            realization[i][0],
            realization[i][1],
            realization[i][2], i
        );
        std::cout << "Collected " << correctionTerms[ind].size()
                  << " linear correction terms for indices " << i << std::endl;
    }

    return correctionTerms;
}

/**
 * Get all possible combinations of indices for a certain level
 * For example level 2 with 4 components -> 0011,0101,0110,1001,1010,1100
 * @param nb_components The number of unshared output functions
 * @param level The number of unshared output functions for which we want to combine the correction terms
 */
std::vector<Indices> getCombinations(int nb_components, int level) {
    std::vector<Indices> result;
    // Generate level-combinations of the components
    // Selector vector, contains true exactly level times
    std::vector<bool> selector(nb_components);
    std::fill(selector.begin() + nb_components - level, selector.end(), true);
    do {
        Indices ind;
        for(int i = 0; i < nb_components; ++i) {
            if(selector[i])
                ind.indices.push_back(i);
        }
        result.push_back(ind);

    } while(std::next_permutation(selector.begin(), selector.end()));
    return result;
}

/**
 * Checks if there exists indices of the level below the current one
 * For example if ind equals 0111 check if there exist correction terms for 0110,0101,0011
 * If one of these terms has no correction terms there is no point in looking for correction terms for 0111
 */
bool canAddIndices(const std::map<Indices, std::vector<VecCorrectionFunction>>& functions,
    const Indices& ind) {
    // Is there at least one correction function for each level - 1 index
    for(std::size_t i = 0; i < ind.indices.size(); ++i) {
        Indices red_ind = ind;
        red_ind.indices.erase(red_ind.indices.begin() + i);
        if(!functions.count(red_ind))
            return false;
    }
    return true;
}

//A comparer for veccorrectionfunctions so that they can be used in a map
struct Comparer {
    bool operator() (const VecCorrectionFunction& b1, const VecCorrectionFunction& b2) const {
        std::size_t size1 = b1.size();
        std::size_t size2 = b2.size();
        if (size1 < size2)
            return true;
        else if (size1 > size2)
            return false;
        for(std::size_t i=0;i<b1.size();i++) {
            std::size_t size11 = b1[i].size();
            std::size_t size22 = b2[i].size();
            if (size11 < size22)
                return true;
            else if (size11 > size22)
                return false;
            for (std::size_t j=0;j<b1[i].size();j++) {
                for (std::size_t k=0;k<b1[i][j].size();k++) {
                    if (b1[i][j][k] < b2[i][j][k])
                        return true;
                    else if (b1[i][j][k] > b2[i][j][k])
                        return false;
                }
            }
        }
        return false;
    }
};

/**
 * Transforms VecCorrectionFunctions of the form (a,b,c) to a map (b,c) -> a
 * Useful for comparing higher levels of indices
 */
std::map<VecCorrectionFunction,std::vector<CorrectionFunction>,Comparer> getMap(std::vector<VecCorrectionFunction> functions) {
    std::map<VecCorrectionFunction,std::vector<CorrectionFunction>,Comparer> results;
    for (size_t i=0;i<functions.size();i++) {
        VecCorrectionFunction key;
        CorrectionFunction value = functions[i][0];
        for (size_t j=1;j<functions[i].size();j++) {
            key.push_back(functions[i][j]);
        }
        if (!results.count(key))
            results[key] = std::vector<CorrectionFunction>();
        results.at(key).push_back(value);
    }
    return results;
}

/**
 * In order to combine (a,c,d) and (b,c,d) two maps {(c,d)} -> {a} and {(c,d)} -> {b} are used
 * For every key (c,d) that exists in both maps, the list of values {a} and {b} are combined with a cross product
 * This allows for a fast, efficient combination of lower levels of indices into higher ones
 */
std::vector<VecCorrectionFunction> combineFunctions(std::map<VecCorrectionFunction,std::vector<CorrectionFunction>,Comparer>& possibilities1, std::map<VecCorrectionFunction,std::vector<CorrectionFunction>,Comparer>& possibilities2) {
    std::vector<VecCorrectionFunction> results;
    for(const auto& pair : possibilities1) {
        if (possibilities2.count(pair.first)) {
            for (const auto& possibility1:pair.second) {
                for(const auto& possibility2:possibilities2.at(pair.first)) {
                    VecCorrectionFunction result;
                    result.push_back(possibility1);
                    result.push_back(possibility2);
                    result.insert(result.end(),pair.first.begin(),pair.first.end());
                    results.push_back(result);
                }
            }
        }
    }
    return results;
}

/**
 * Finds all possible combinations of correction terms for level amount of unshared outputs
 * The correction terms are not yet checked for uniformity, this function only gathers a list of possible candidates
 * @param functions The collection of possible veccorrectionfunctions for each group of indices
 * @param nb_components The number of unshared output functions
 * @param level The number of unshared output functions for which we want to combine the correction terms
 */
std::map<Indices, std::vector<VecCorrectionFunction>> combineCorrectionFunctions(
    const std::map<Indices, std::vector<VecCorrectionFunction>>& functions,
    int nb_components, int level) {

    std::map<Indices, std::vector<VecCorrectionFunction>> result;
    
    for(Indices& ind : getCombinations(nb_components, level)) {
        std::cout << "Combining correction functions for indices ";
        for (auto i: ind.indices)
            std::cout << i << ' ';
        std::cout << std::endl;
        if(!canAddIndices(functions, ind))
            continue; 
        // Remove the first index
        Indices ind_no_first = ind.without(0);
        // Remove the second index
        Indices ind_no_second = ind.without(1);
        std::cout << "\touter loop size is "  << functions.at(ind_no_first).size()
                  << "\n\tinner loop size is " << functions.at(ind_no_second).size()
                  << std::endl;
        std::size_t count = 0;
        if(level==2) {
            for(const auto& correction1 : functions.at(ind_no_first)) {
                for(const auto& correction2 : functions.at(ind_no_second)) {
                    //Useful when you want a quick result for debugging purposes
                    //int v1 = rand() % 300; 
                    //if (v1 == 1) {
                    VecCorrectionFunction correction = correction2;
                    correction.insert(
                        correction.end(), correction1.begin(), correction1.end()
                    );
                    if(!result.count(ind))
                        result[ind] = std::vector<VecCorrectionFunction>();
                    result[ind].push_back(correction);
                    //}
                }
            }
        }
  
        else {
            std::map<VecCorrectionFunction,std::vector<CorrectionFunction>,Comparer> possibilities1 = getMap(functions.at(ind_no_first));
            std::map<VecCorrectionFunction,std::vector<CorrectionFunction>,Comparer> possibilities2 = getMap(functions.at(ind_no_second));
            if(!result.count(ind))
                result[ind] = std::vector<VecCorrectionFunction>();
            result[ind] = combineFunctions(possibilities2,possibilities1);
        }
        
    }
    return result;
}

char nth_letter(int n)
{
    if (n==0) {
        std::cout << "The function nth_letter takes inputs from 1-26, 0 is not a valid input" << std::endl;
        return 'A';
    }
    return "ABCDEFGHIJKLMNOPQRSTUVWXYZ"[n-1];
}

/**
 * Writes the direct sharing and all correction terms to a file
 * @param directSharing first element contains the constant terms, second element contains the linear terms,
 *        third element contains the quadratic terms 
 */
void writeFunctions(const std::string& directory, const std::string& filename, 
    const std::map<Indices, std::vector<VecCorrectionFunction>>& functions,
    std::tuple<std::bitset<OUTPUT_SHARES_BITS>,std::vector<std::bitset<INPUT_SHARES_BITS>>,std::vector<std::vector<std::bitset<INPUT_SHARES_BITS>>>>& directSharing) {
    
    std::ofstream ofs(directory + '/' + filename);
    ofs << "Direct sharing for the function f(...,X,Y,Z) with outputs A,B,C,... :\n\n\n";
    for(std::size_t i=0;i<OUTPUT_BITS;i++) {
        ofs << nth_letter(i+1) << ": ";
        if (std::get<0>(directSharing)[INPUT_SHARES*i])
            ofs << "+1 ";
        for(std::size_t j=0;j<OUTPUT_SHARES_BITS;j++) {
            if (std::get<1>(directSharing)[INPUT_SHARES*i][j])
                ofs << "+ " << nth_letter((26-INPUT_BITS+1) + j/3) << ((j+1)%INPUT_SHARES)+1 << " ";
        }
        for(std::size_t j=0;j<OUTPUT_SHARES_BITS;j++) {
            for(std::size_t k=0;k<OUTPUT_SHARES_BITS;k++) {
                if (std::get<2>(directSharing)[3*i][j][k])
                    ofs << "+ " << nth_letter((26-INPUT_BITS+1) + j/3) << ((j+1)%INPUT_SHARES)+1 << "*" << nth_letter((26-INPUT_BITS+1) + k/3) <<((k+1)%INPUT_SHARES)+1 << " ";
            }
        }
        ofs << "\n";
    }
    ofs << "\n\nCorrection terms:\n\n";
    ofs << "\n------------\n\n";
    for(auto& pair : functions) {
        for(const VecCorrectionFunction& vf : pair.second) {
            for(std::size_t k=0;k<vf.size();k++) {
                for(std::size_t i=0;i<vf[k].size();i++) {
                    ofs << nth_letter(k+1) << (i+1) << ": ";
                    for(std::size_t j=0;j<vf[k][i].size();j++) {
                        if(vf[k][i][j]) 
                            ofs << "+ " << nth_letter((26-INPUT_BITS+1) + j/3) << ((j+1)%INPUT_SHARES)+1 << " ";
                    }
                ofs << "\n";
                }
            }
            ofs << "\n------------\n\n";
        }
    }
    ofs << std::endl;
}

//Returns the total amount of VecCorrectionFunctions in the entire map
std::size_t countCorrectionFunctions(
    const std::map<Indices, std::vector<VecCorrectionFunction>>& functions)
{
    std::size_t total_count = 0;
    for(const auto& pair : functions)
        total_count += pair.second.size();
    return total_count;
}

//Checks if the given sharedInput is a correct sharing of the given input
bool isCorrectSharing(const InputBitArray& input, const SharedInputBitArray& sharedInput) {
    for (std::size_t i = 0; i<INPUT_BITS;i++) {
        bool temp = 0;
        for (std::size_t j = 0; j<INPUT_SHARES;j++) {
            temp = temp ^ sharedInput[i*INPUT_SHARES+j];
        }
        if (temp != input[i]) {
            return false;
        }
    }
    return true;
    
}

//Creates a global variable (sorry) that contains all valid sharings for each possible input
void createGlobalSharingTable() {
    for(std::size_t i = 0; i < INPUT_SIZE; ++i) {
        InputBitArray input(i); 
        std::vector<std::size_t> allShares;
        for(std::size_t j = 0; j < INPUT_SHARES_SIZE; ++j) {
            SharedInputBitArray shared_input(j);
            if (isCorrectSharing(input,shared_input)) {
                allShares.push_back(j);
            }
        }
        globalValidSharesTable.push_back(allShares);
    }       
}

/**
 * Checks whether the given boolean functions are uniform (have the same count of each output when checked over all possible inputs)
 * This is workhorse of the program, about 99% of the time is spend here
 * If you can optimize this function, do so
 * @param components The truthTable for all shares of the output functions which are combined
 * @return true of the given boolean functions are uniform together, false otherwise
 */
bool checkUniformity(const std::vector<std::vector<BlnFunction>>& components) {
    const std::size_t expected_count = std::pow(
        2, 2 * INPUT_BITS - 2 * components.size()
);
    //int test2;
    std::vector<std::vector<std::size_t>> counts(INPUT_SIZE, std::vector<std::size_t>(OUTPUT_SHARES_SIZE,0));
    for(std::size_t i = 0; i < INPUT_SIZE; ++i) {
        InputBitArray input(i); 
        std::bitset<OUTPUT_SHARES_BITS> outputs;
        for(std::size_t validShare : globalValidSharesTable[i]) {
            outputs.reset();
            for(std::size_t k = 0; k < components.size(); ++k) {
                for (std::size_t l = 0; l < OUTPUT_SHARES; ++l) {
                    outputs[k*OUTPUT_SHARES+l] = components[k][l][validShare];
                }
            }
            //test2 = (int) outputs.to_ulong();
            counts[i][(int) outputs.to_ulong()] += 1;
            if(counts[i][(int) outputs.to_ulong()] > expected_count) {
                return false;
            }
        }
    }
//Sanity check, only use for debugging since the performance hit is significant
//    for(std::size_t i = 0; i < INPUT_SIZE; ++i) {
//        for (std::size_t j = 0; j < OUTPUT_SHARES_SIZE; ++j) {
//            if (counts[i][j] != 0 && counts[i][j] != expected_count) {
//                std::cout << "Sanity check failed" << std::endl;
//                std::cout << counts[i][j] << std::endl;
//                std::cout << expected_count << std::endl;
//            } 
//        }
//    }
    return true;
}

typedef std::vector<VecCorrectionFunction>::const_iterator FunctionIter;

/**
 * 
 * @param directory
 * @param realization
 * @param nb_input_variables
 * @param indices
 * @param begin
 * @param end
 * @param batch_nb
 * @param batch_size
 * @return 
 */
std::vector<VecCorrectionFunction> makeBatchUniformWith(
    const std::string& directory,
    const std::vector<std::vector<BlnFunction>>& realization,
    std::size_t nb_input_variables,
    const Indices& indices,
    FunctionIter begin, FunctionIter end,
    std::size_t batch_nb,
        const std::size_t batch_size) {

    std::vector<VecCorrectionFunction> good_correction_functions;

    std::ofstream os(directory + "/batch-" + std::to_string(batch_nb) + ".log");
    
    std::size_t nb_done = 0;
    std::size_t nb_found = 0;

    // Try all possible correction functions
    for(FunctionIter it = begin; it != end; ++it) {
        // Construct the vectorial boolean function, i.e. add vf to realization
        // The result is a flattened vector of vectorial boolean functions
        std::vector<std::vector<BlnFunction>> corrected_realization;
        std::size_t i_corrected = 0;
        std::vector<BlnFunction> temp_realization;
        for(std::size_t i : indices.indices) {
            for(std::size_t s = 0; s < realization[i].size(); ++s) {
                temp_realization.push_back(BlnFunction(globalTruthTable[(int) (*it)[i_corrected][s].to_ulong()])+realization[i][s]);
            }
            corrected_realization.push_back(temp_realization);
            temp_realization.clear();
            ++i_corrected;
        }

        if(checkUniformity(corrected_realization)) {
            ++nb_found;
            good_correction_functions.push_back(*it);
            //std::cout << "Found " << nb_found << " of " << nb_done << " correction functions."
            //   << std::endl;
        }
        ++nb_done;
        if ((nb_done*100)/batch_size != ((nb_done-1)*100)/batch_size && ((nb_done*100)/batch_size % 5 == 0))
            std::cout << (nb_done*100)/batch_size << " percent done" << std::endl;
    }
    return good_correction_functions;

}




std::map<Indices, std::vector<VecCorrectionFunction>> makeUniformWith(
    const std::string& directory,
    const std::vector<std::vector<BlnFunction>>& realization,
    std::size_t nb_input_variables,
    const std::map<Indices, std::vector<VecCorrectionFunction>>& functions,
        const std::size_t level) {

    std::map<Indices, std::vector<VecCorrectionFunction>> filtered_functions;
    
 
    // Produce batches for each set of indices
    for(const auto& pair : functions) {
        std::size_t batch_nb = 0;
        const std::size_t batch_size = pair.second.size() / MAX_THREADS;
        std::vector<std::future<std::vector<VecCorrectionFunction>>> results;
        auto it = pair.second.begin();
        if(pair.second.size() > 2*MAX_THREADS) {
            for(; it < pair.second.end() - batch_size; it += batch_size) {
                results.emplace_back(std::async(std::launch::async,
                    makeBatchUniformWith, directory, realization,
                    nb_input_variables, pair.first, it,
                    it + batch_size, batch_nb, batch_size
                ));
                ++batch_nb;
            }
            // Start batch with remainder
            results.emplace_back(std::async(std::launch::async,
                makeBatchUniformWith, directory, realization,
                nb_input_variables, pair.first, it,
                pair.second.end(), batch_nb, batch_size
            ));

            // Retrieve the results
            int i = 0;
            for(auto& future : results) {
                i +=1;
    //            std::vector<VecCorrectionFunction> corrections = makeBatchUniformWith(directory, realization,
    //                nb_input_variables, pair.first, it,
    //                pair.second.end(), batch_nb
    //            );
                std::vector<VecCorrectionFunction> corrections = future.get();
                if(!corrections.empty() && !filtered_functions.count(pair.first))
                    filtered_functions[pair.first] = std::vector<VecCorrectionFunction>();
                filtered_functions[pair.first].insert(filtered_functions[pair.first].end(),corrections.begin(),corrections.end());
            }
        }
        else {
                std::vector<VecCorrectionFunction> corrections = makeBatchUniformWith(directory, realization,
                    nb_input_variables, pair.first, it,
                    pair.second.end(), batch_nb, pair.second.size()
                );
                if(!corrections.empty() && !filtered_functions.count(pair.first))
                    filtered_functions[pair.first] = std::vector<VecCorrectionFunction>();
                filtered_functions[pair.first].insert(filtered_functions[pair.first].end(),corrections.begin(),corrections.end());
        }
    }
    return filtered_functions;
}

bool makeUniform(
    const std::string& filename,
    std::map<Indices, std::vector<VecCorrectionFunction>> functions,
    const std::vector<std::vector<BlnFunction>> realization,
    std::size_t nb_input_variables, 
    std::tuple<std::bitset<OUTPUT_SHARES_BITS>,std::vector<std::bitset<INPUT_SHARES_BITS>>,std::vector<std::vector<std::bitset<INPUT_SHARES_BITS>>>>& directSharing) {

    std::size_t level = 2;
    if (checkUniformity(realization)) {
        std::cout << "already uniform" << std::endl;
        return true;
    }
       
    while(level <= realization.size() && !functions.empty()) {
        std::cout << "Combining " << countCorrectionFunctions(functions)
                  << " correction functions for level " << level << '.'
                  << std::endl;
        auto candidates = combineCorrectionFunctions(
            functions, realization.size(), level
        );

        std::cout << "Only using indices ";
        for(const auto& pair : candidates)
            for (auto i: pair.first.indices)
                std::cout << i << ' ';
            std::cout << ",";
        std::cout << std::endl;

        std::cout << "Starting search with " << countCorrectionFunctions(candidates)
                  << " candidates." << std::endl;
        functions = makeUniformWith("output", realization, nb_input_variables, candidates, level);

        std::cout << "Found " << countCorrectionFunctions(functions)
                  << " correction functions." << std::endl;
        level++;
    }
    
    std::cout << countCorrectionFunctions(functions) << std::endl;
    writeFunctions("output", filename, functions, directSharing);
    std::cout << "Process stopped at level " << level << '.' << std::endl;
    return (countCorrectionFunctions(functions)!=0);
}

/**
 * Constructs the truth tables for all shares of all outputs
 */
std::vector<std::vector<BlnFunction>> createTruthTable(std::bitset<OUTPUT_SHARES_BITS>& constant_bits, std::vector<std::bitset<INPUT_SHARES_BITS>>& linear_bits,std::vector<std::vector<std::bitset<INPUT_SHARES_BITS>>>& quadratic_bits) {
    std::vector<std::vector<BlnFunction>> truthTables;
    std::vector<BlnFunction> tempTable;
    for (std::size_t i = 0;i<INPUT_SHARES_BITS;i++) {
        std::bitset<INPUT_SHARES_SIZE> truthTable;
        for (std::size_t j = 0;j<INPUT_SHARES_SIZE;j++) {
            std::bitset<INPUT_SHARES_BITS> bits(j);
            bool result = 0;
            std::bitset<INPUT_SHARES_BITS> im_result;
            for (std::size_t k = 0; k<INPUT_SHARES_BITS;k++) {
                im_result[k] = 0;
                for (std::size_t l = 0; l<INPUT_SHARES_BITS;l++) {
                    im_result[k] = im_result[k] ^ (quadratic_bits[i][k][l] & bits[l] & bits[k]);
                }
            }
            for (std::size_t k = 0; k<INPUT_SHARES_BITS;k++) {
                result ^= im_result[k];
            }
            for (std::size_t k = 0; k<INPUT_SHARES_BITS;k++) {
                result = result ^ (linear_bits[i][k] & bits[k]);
            }
            truthTable[j] = result ^ constant_bits[i];
        }
        tempTable.push_back(BlnFunction(truthTable));
        if ((i+1)%INPUT_SHARES == 0) {
            truthTables.push_back(tempTable);
            tempTable.clear();
        }
    }
    return truthTables;
}

/**
 * Creates a global variable that contains information on which share of which output functions depends on which share of which input functions
 * Useful for checking the non-completeness constraint
 */
void getDependence(std::vector<std::bitset<INPUT_SHARES_BITS>>& linear_bits,std::vector<std::vector<std::bitset<INPUT_SHARES_BITS>>>& quadratic_bits) {
    std::vector<std::bitset<INPUT_SHARES_BITS>> dependence;
    for (std::size_t i=0;i<INPUT_SHARES_BITS;i++) {
        std::bitset<INPUT_SHARES_BITS> temp;
        for (std::size_t j=0;j<INPUT_SHARES_BITS;j++) {
            temp[j] = linear_bits[i][j];
        }
        dependence.push_back(temp);
    }
    for (std::size_t i=0;i<INPUT_SHARES_BITS;i++) {
        for (std::size_t j=0;j<INPUT_SHARES_BITS;j++) {
            for (std::size_t k=0;k<INPUT_SHARES_BITS;k++) {
                dependence[i][j] = quadratic_bits[i][j][k] | dependence[i][j];
            }
            for (std::size_t k=0;k<INPUT_SHARES_BITS;k++) {
                dependence[i][j] = quadratic_bits[i][k][j] | dependence[i][j];
            }
        }
    }
    globalDependence = dependence;
}


/**
 * Converts an unshared truthTable, f, into ANF
 * Feel free to bring this function up to date to c++11 using bitsets/std::arrays/vectors
 */
int* getANF(const int f[INPUT_SIZE], int p[INPUT_SIZE]){
    int i, j;
    for (i = 0; i < INPUT_SIZE; ++i) p[i] = f[i];

    for (i = 1; i < INPUT_SIZE; i <<= 1){
        for (j = i; j < INPUT_SIZE; j = (j + 1) | i){
            p[j] ^= p[j ^ i];
        }
    }
    return p;
}

/**
 * Converts ANF to a form with the constant, linear and quadratic terms in different variables
 * Each output function has one constant bit, one vector of linear bits and one matrix of quadratic terms
 */
std::tuple<std::bitset<OUTPUT_BITS>,std::vector<std::bitset<INPUT_BITS>>,std::vector<std::vector<std::bitset<INPUT_BITS>>>> ANFToUnsharedInput(const int p[INPUT_SIZE]) {
    
    std::vector<std::vector<std::bitset<INPUT_BITS>>> quadratic_bits(INPUT_BITS, std::vector<std::bitset<INPUT_BITS>>(INPUT_BITS,std::bitset<INPUT_BITS>()));
    std::bitset<OUTPUT_BITS> constant_bits;
    std::vector<std::bitset<INPUT_BITS>> linear_bits(INPUT_BITS,std::bitset<INPUT_BITS>());
    for (std::size_t i=0;i<OUTPUT_BITS;i++) {
        constant_bits[i] = (p[0] >> i) & 1;
    }
    for (std::size_t i=0;i<INPUT_BITS;i++) {
        for (std::size_t j=0;j<INPUT_BITS;j++) {
            linear_bits[j][i] = (p[1 << i] >> j) & 1;
            }
    }
    for (std::size_t i=1;i<INPUT_BITS;i++) {
        for (std::size_t j=i;j<INPUT_BITS;j++) {
            for (std::size_t k=0;k<INPUT_BITS;k++) {
                quadratic_bits[k][i-1][j] = (p[(1<<j) + (1<<(i-1))] >> k) & 1;
            }
        }
    }
    return std::make_tuple(constant_bits,linear_bits,quadratic_bits); 
}

/**
 * Converts the unshared form to the shared form
 * Each share of each output function has one constant bit, one vector of linear bits and one matrix of quadratic terms
 */
std::tuple<std::bitset<OUTPUT_SHARES_BITS>,std::vector<std::bitset<INPUT_SHARES_BITS>>,std::vector<std::vector<std::bitset<INPUT_SHARES_BITS>>>> unsharedToSharedInput(std::bitset<INPUT_BITS> constant_bits_old, std::vector<std::bitset<INPUT_BITS>> linear_bits_old, std::vector<std::vector<std::bitset<INPUT_BITS>>> quadratic_bits_old) {
    std::bitset<INPUT_SHARES_BITS> constant_bits;
    std::vector<std::bitset<INPUT_SHARES_BITS>> linear_bits(INPUT_SHARES_BITS,std::bitset<INPUT_SHARES_BITS>());
    std::vector<std::vector<std::bitset<INPUT_SHARES_BITS>>> quadratic_bits(INPUT_SHARES_BITS, std::vector<std::bitset<INPUT_SHARES_BITS>>(INPUT_SHARES_BITS,std::bitset<INPUT_SHARES_BITS>()));
    for(std::size_t i=0;i<constant_bits.size();i++) {
        constant_bits[INPUT_SHARES*i] = constant_bits_old[i];
    }
    for(std::size_t i=0;i<linear_bits_old.size();i++) {
        for(std::size_t j=0;j<INPUT_BITS;j++) {
            for(std::size_t k=1;k<INPUT_SHARES;k++) {
                linear_bits[INPUT_SHARES*i+(k-1)][INPUT_SHARES*j+k] = linear_bits_old[i][j];
            }
            linear_bits[INPUT_SHARES*i+INPUT_SHARES-1][INPUT_SHARES*j] = linear_bits_old[i][j];
            //give second share to the first sharing function
            //linear_bits[3*i][3*j+1] = linear_bits_old[i][j];
            //give third share to the second sharing function
            //linear_bits[3*i+1][3*j+2] = linear_bits_old[i][j];
            //give first share to the third sharing function
            //linear_bits[3*i+2][3*j] = linear_bits_old[i][j];
        }
    }
    //hardcoded with 3 shares for now, should be changed to a general form
    for(std::size_t i=0;i<quadratic_bits_old.size();i++) {
        for(std::size_t j=0;j<INPUT_BITS;j++) {
            for(std::size_t k=0;k<INPUT_BITS;k++) {
                quadratic_bits[3*i][3*j+1][3*k+1] = quadratic_bits_old[i][j][k];
                quadratic_bits[3*i][3*j+2][3*k+1] = quadratic_bits_old[i][j][k];
                quadratic_bits[3*i][3*j+1][3*k+2] = quadratic_bits_old[i][j][k];
                quadratic_bits[3*i][3*j+2][3*k+2] = quadratic_bits_old[i][j][k];
                quadratic_bits[3*i+1][3*j][3*k] = quadratic_bits_old[i][j][k];
                quadratic_bits[3*i+1][3*j+2][3*k] = quadratic_bits_old[i][j][k];
                quadratic_bits[3*i+1][3*j][3*k+2] = quadratic_bits_old[i][j][k];
                quadratic_bits[3*i+2][3*j][3*k+1] = quadratic_bits_old[i][j][k];
                quadratic_bits[3*i+2][3*j+1][3*k] = quadratic_bits_old[i][j][k];
            }
        }
    }
    return std::make_tuple(constant_bits,linear_bits,quadratic_bits);
}


/**
 * Reads an unshared binary truth table from a file
 */
int* readTruthTable(const std::string& filename) {
    int truthTable[INPUT_SIZE];
    std::ifstream ifs("input/" + filename);
    if (!ifs.is_open()) {
        std::cout << "failed to open " << filename << '\n';
        return nullptr;
    }
    std::string line;
    int i = 0;
    while(std::getline(ifs, line) && i<INPUT_SIZE) {
        if(line.size() != INPUT_BITS)
            std::cout << "Warning: incorrect input size given, got " << line.size()
                      << " expecting " << INPUT_BITS << '.' << std::endl;
        truthTable[i] = (int) std::bitset<INPUT_BITS>(line).to_ulong();
        i++;
    }
    return truthTable;
}

/**
 * Reads a CorrectionFunction from a file
 * OUTDATED
 */
std::vector<std::bitset<INPUT_SHARES_BITS>> readCandidate(const std::string& filename) {
    
    std::ifstream ifs("control/" + filename);
    std::string line;
    std::vector<std::bitset<INPUT_SHARES_BITS>> correctionFunction;
    int i = 0;
    while(std::getline(ifs, line) && i<INPUT_SHARES_BITS) {
        if(line.size() != INPUT_SHARES_BITS)
            std::cout << "Warning: incorrect input size given, got " << line.size()
                      << " expecting " << INPUT_SHARES_BITS << '.' << std::endl;
        correctionFunction.push_back(std::bitset<INPUT_SHARES_BITS>(line));
        i++;
    }
    return correctionFunction;
}

/**
 * Converts a decimal unshared truth table to a binary unshared truth table
 * also checks if the given truth table is a valid s-box
 */
void decimalToBinaryInputConvertor(const std::string& filename) {
    
    std::ifstream ifs("decimalInput/" + filename);
    std::ofstream ofs("input/" + filename);
    std::string line;
    int i = 0;
    std::bitset<INPUT_SIZE> control;
    while(std::getline(ifs, line) && i<INPUT_SIZE) {
        int number = std::stoi(line);
        if (number >= INPUT_SIZE) {
            std::cout << "bad input given, check original file for errors" << std::endl;
            return;
        }
        std::string binary = std::bitset<INPUT_BITS>(number).to_string();
        control[number] = 1;
        ofs << binary << std::endl;
        i++;
    }
    if (control.count()!= INPUT_SIZE) {
        std::cout << "bad input given, check original file for errors" << std::endl;
    }
}

int main(int argc, char *argv[]) {
    
    if(argc < 2) {
        std::cerr << "Usage: <mode> <filename>. where mode is one of the following:/n 0 : normal operation/n 1 : control/n 2 : decimal to binary convertor"
                  << std::endl;
        return 1;
    }
    int mode = std::stoi(argv[1]);
    if (mode == 2) {
        decimalToBinaryInputConvertor(argv[2]);
        return 0;
    }
    int* truthTable = readTruthTable(argv[2]);
    if (truthTable==nullptr) {
        return 1;
    }
    int p[INPUT_SIZE];
    auto tuple = ANFToUnsharedInput(getANF(truthTable,p));
    auto sharedTuple = unsharedToSharedInput(std::get<0>(tuple),std::get<1>(tuple),std::get<2>(tuple));
    if(mode == 1) {
        std::vector<std::bitset<INPUT_SHARES_BITS>> correctionTerm = readCandidate(argv[2]);
        for(std::size_t i=0;i<INPUT_SHARES_BITS;i++) {
            for(std::size_t j=0;j<INPUT_SHARES_BITS;j++) {
            correctionTerm[i][j] = std::get<1>(sharedTuple)[i][j] + correctionTerm[i][j];
            }
        }
        auto result = createTruthTable(std::get<0>(sharedTuple),correctionTerm,std::get<2>(sharedTuple));
        getDependence(correctionTerm,std::get<2>(sharedTuple));
        for(std::size_t i=0;i<INPUT_SHARES_BITS;i++) {
            bool first_bit = false;
            bool second_bit = false;
            bool third_bit = false;
            for(std::size_t j = 0; j < INPUT_BITS; j += INPUT_SHARES) {
                if(globalDependence[i][j]) {
                    first_bit = true;
                }
                if(globalDependence[i][j+1]) {
                    second_bit = true;
                }
                if(globalDependence[i][j+2]) {
                    third_bit = true;
                }
            }
            if (first_bit && second_bit && third_bit) {
                std::cout << "non-uniformity not met for " << i << std::endl;
            }
            
            std::cout << first_bit << second_bit << third_bit << std::endl;
        }
        createGlobalSharingTable();
        std::cout << checkUniformity(result) << std::endl;
    }
    else {
        auto result = createTruthTable(std::get<0>(sharedTuple),std::get<1>(sharedTuple),std::get<2>(sharedTuple));
        getDependence(std::get<1>(sharedTuple),std::get<2>(sharedTuple));
        std::size_t nb_inputs = INPUT_BITS;
        createGlobalSharingTable();
        bool solutionFound = false;
        while(!solutionFound) {
            solutionFound = makeUniform(
                argv[2], buildCorrectionTerms(result),
                result, nb_inputs, sharedTuple
            ); 
            MAX_HAMMING_WEIGHT_TOTAL *= 1.5;
            MAX_HAMMING_WEIGHT = MAX_HAMMING_WEIGHT_TOTAL/2.5;
            std::cout << MAX_HAMMING_WEIGHT << " " << MAX_HAMMING_WEIGHT_TOTAL << std::endl;
        }
    }
    return 0;
}
