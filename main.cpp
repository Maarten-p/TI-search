#include <cinttypes>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <future>
#include <map>
#include <tuple>
#include "BooleanFunction.h" 

constexpr std::size_t MAX_THREADS = 8;
constexpr std::size_t MAX_HAMMING_WEIGHT = 2;
constexpr std::size_t INPUT_BITS = 4;
constexpr std::size_t INPUT_SHARES = 3;
constexpr std::size_t OUTPUT_SHARES = 3;
constexpr std::size_t OUTPUT_BITS = 4;
constexpr std::size_t INPUT_SHARES_BITS = INPUT_BITS * INPUT_SHARES;
constexpr std::size_t OUTPUT_SHARES_BITS = OUTPUT_BITS * OUTPUT_SHARES;
constexpr std::size_t INPUT_SIZE = std::pow(2, INPUT_BITS * INPUT_SHARES);
constexpr std::size_t OUTPUT_SIZE = std::pow(2, OUTPUT_BITS * OUTPUT_SHARES);

typedef BooleanFunction<INPUT_SIZE> BlnFunction;
typedef std::array<std::bitset<INPUT_SHARES_BITS>, OUTPUT_SHARES> CorrectionFunction;
typedef std::vector<CorrectionFunction> VecCorrectionFunction;
typedef std::bitset<INPUT_BITS> InputBitArray;
typedef std::bitset<INPUT_SHARES * INPUT_BITS> SharedInputBitArray;

std::vector<std::bitset<OUTPUT_SIZE>> globalTruthTable(INPUT_SIZE,0);
std::vector<std::vector<std::size_t>> globalValidSharesTable;

struct Indices {
    std::vector<std::size_t> indices;

    Indices without(std::size_t i) const {
        Indices ind;
        ind.indices = indices;
        ind.indices.erase(ind.indices.begin() + i);
        return ind;
    }
};

std::ostream& operator<<(std::ostream& os, const Indices& ind)
{
    os << '(';
    for(std::size_t i : ind.indices)
        os << i << ' ';
    os << ')';
}

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
 * @return true if every (share)th share of every input is zero, false otherwise
 */
bool areSharesZero(const std::bitset<INPUT_SHARES_BITS>& bits, std::size_t share) {
    for(std::size_t i = share; i < bits.size(); i += INPUT_SHARES) {
        if(bits[i]) 
            return false;
    }
    return true;
}

bool hammingWeightConstraint(const std::bitset<INPUT_SHARES_BITS> bits) {
    return bits.count() <= MAX_HAMMING_WEIGHT;
}

BlnFunction truthTable1(const std::bitset<INPUT_SHARES_BITS>& correction_term, BlnFunction origin) {
    BlnFunction::BitArray truthTable = origin.getTruthTable();
    for (std::size_t j = 0;j<INPUT_SIZE;j++) {
        std::bitset<INPUT_SHARES_BITS> bits(j);
        bool result = 0;
        std::bitset<INPUT_SHARES_BITS> im_result = bits & correction_term;
        for (std::size_t k = 0; k<INPUT_SHARES_BITS;k++) {
            result ^= im_result[k];
        }
        truthTable[j] = truthTable[j] ^ result;
    }
    return BlnFunction(truthTable);
}

void addToTruthTable(std::bitset<INPUT_SHARES_BITS> correction_term) {
    BlnFunction::BitArray truthTable;
    for (std::size_t j = 0;j<INPUT_SIZE;j++) {
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
 * @note f1, f2 and f3 are assumed to be non-reduced, but nevertheless
 *  noncomplete.
 * @note each VecCorrectionFunction is of size one
 * @warning this doesn't work for more than 8 * sizeof(int) input variables.
 */
std::vector<VecCorrectionFunction> getLinearCorrections(const BlnFunction& f1,
    const BlnFunction& f2, const BlnFunction& f3) {

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

    for(std::size_t i = 0; i < spectrum1.size(); ++i) {
        std::bitset<INPUT_SHARES_BITS> bits_i(i);
        if(!areSharesZero(bits_i, 0) || spectrum1[i] != 0 || !hammingWeightConstraint(bits_i))
            continue;

        for(std::size_t j = 0; j < spectrum2.size(); ++j) {
            std::bitset<INPUT_SHARES_BITS> bits_j(j);
            if(!areSharesZero(bits_j, 1) || spectrum2[j] != 0 || !hammingWeightConstraint(bits_j))
                continue;
            
            std::bitset<INPUT_SHARES_BITS> bits_l = bits_i ^ bits_j;
            std::size_t l = bits_l.to_ulong();
            if(!areSharesZero(bits_l, 2) || spectrum3[l] != 0 || !hammingWeightConstraint(bits_l))
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

std::map<Indices, std::vector<VecCorrectionFunction>> buildCorrectionTerms(
    const std::vector<std::vector<BlnFunction>>& realization)
{
    std::cout << "Collecting linear correction terms." << std::endl;
    std::map<Indices, std::vector<VecCorrectionFunction>> correctionTerms;
    for(std::size_t i = 0; i < realization.size(); ++i) {
        Indices ind;
        ind.indices.push_back(i);
        correctionTerms[ind] = getLinearCorrections(
            realization[i][0],
            realization[i][1],
            realization[i][2]
        );
        std::cout << "Collected " << correctionTerms[ind].size()
                  << " linear correction terms for indices " << i << std::endl;
    }

    return correctionTerms;
}

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

bool shouldAddCorrection(const std::map<Indices, std::vector<VecCorrectionFunction>>& functions,
    const Indices& ind, const VecCorrectionFunction& correction) {

    for(std::size_t i = 2; i < ind.indices.size(); ++i) {
        bool is_found = false;
        const Indices ind_part = ind.without(i);
        VecCorrectionFunction correction_part(
            correction.begin(), correction.begin() + i
        );
        correction_part.insert(
            correction_part.end(),
            correction.begin() + i, correction.end()
        );
        for(const auto& vec : functions.at(ind_part)) {
            bool all_equal = true;
            for(std::size_t j : ind_part.indices) {
                if(vec[j] != correction_part[j]) {
                    all_equal = false;
                    break;
                }
            }
            if(all_equal) {
                is_found = true;
                break;
            }
        }

        if(!is_found)
            return false;
    }
    return true;

}

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
                    //int v1 = rand() % 5; 
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
//            std::vector<int> common_indices;
//            std::set_intersection(ind_no_first.indices.begin(), ind_no_first.indices.end(),
//                          ind_no_second.indices.begin(), ind_no_second.indices.end(),
//                          std::back_inserter(common_indices));
//            if (common_indices.size()<level-2)
//                continue;
//            int no_first_indexes = findUncommonIndex(ind_no_first.indices, common_indices);
//            std::vector<int> no_second_indexes = findUncommonIndex(ind_no_second.indices, common_indices);
            std::map<VecCorrectionFunction,std::vector<CorrectionFunction>,Comparer> possibilities1 = getMap(functions.at(ind_no_first));
            std::map<VecCorrectionFunction,std::vector<CorrectionFunction>,Comparer> possibilities2 = getMap(functions.at(ind_no_second));
            if(!result.count(ind))
                result[ind] = std::vector<VecCorrectionFunction>();
            result[ind] = combineFunctions(possibilities2,possibilities1);
        }
        
    }
    return result;
}

void writeFunctions(const std::string& directory, const std::string& file, 
    const std::map<Indices, std::vector<VecCorrectionFunction>>& functions) {

    std::ofstream ofs(directory + '/' + file);
    for(auto& pair : functions) {
        ofs << 1 << '\n';
        for(const VecCorrectionFunction& vf : pair.second) {
            for(const CorrectionFunction& f : vf) {
                for(const std::bitset<INPUT_SHARES_BITS>& share : f) {
                    ofs << '\t' << share << '\n';
                }
                ofs << "\n\n";
            }
            ofs << "\n------------\n";
        }
    }
    ofs << std::endl;
}


std::size_t countCorrectionFunctions(
    const std::map<Indices, std::vector<VecCorrectionFunction>>& functions)
{
    std::size_t total_count = 0;
    for(const auto& pair : functions)
        total_count += pair.second.size();
    return total_count;
}

template <int n>
std::vector<std::bitset<n>> getNSharing(bool bit) {
    std::vector<std::bitset<n>> result;
    for (int i=0; i<pow(2,n);++i) {
        std::bitset<n> input(i);
        bool checksum = false;
        for (bool current_bit : input) {
            checksum ^= current_bit;
        }
        if (bit==checksum) {
            result.push_back(input);
        }
    }
}

std::vector<std::bitset<3>> get3Sharing(bool bit)
{
    std::vector<std::bitset<3>> result;
    std::bitset<3> current;
    for (bool b1 : { false, true }) {
        current[0] = b1;
        for (bool b2 : { false, true }) {
            current[1] = b2;
            current[2] = b1 ^ b2 ^ bit;
            result.push_back(current);
        }
    }
    return result;
}

//template <int shares, int input>
//std::vector<std::bitset<shares*(input+1)> getCombinations(std::vector<std::bitset<shares*input>> oldSharings, std::vector<std::bitset<shares> newSharings) {
//    std::vector<std::bitset<shares*(input+1)> output;
//    for (std::bitset<shares*input> combination:oldSharings) {
//        for (std::bitset<shares> newcombination:newSharings) {
//            output.push_back(bitset<16> result(combination.to_ulong() * 0x100 + newcombination.to_ulong());)
//        }
//    }
//}
//std::vector<SharedInputBitArray> getSharingsForInput(const InputBitArray& input_bits)
//{
//    std::vector<SharedInputBitArray> result(4,0);
//    for(std::size_t i = 0; i < input_bits.size(); ++i) {
//        for (sharedArray : result) {
//            
//        }
//        input_bits.size()/ (2*(i+1))
//        const std::size_t k = INPUT_SHARES * i;
//        std::vector<std::bitset<INPUT_SHARES>> possible_sharings = getNSharing<INPUT_SHARES>(input_bits[i]);
//        for(std::size_t j = 0; j < possible_sharings.size(); ++j)  {
//            result[j][k] = possible_sharings[j][0];
//            result[j][k+1] = possible_sharings[j][1];
//            result[j][k+2] = possible_sharings[j][2];
//            }
//    }
//    return result;
//}

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

void createGlobalSharingTable() {
    const std::size_t input_size_no_shares = std::pow(2, INPUT_BITS);
    for(std::size_t i = 0; i < input_size_no_shares; ++i) {
        InputBitArray input(i); 
        std::vector<std::size_t> allShares;
        for(std::size_t j = 0; j < INPUT_SIZE; ++j) {
            SharedInputBitArray shared_input(j);
            if (isCorrectSharing(input,shared_input)) {
                allShares.push_back(j);
            }
        }
        globalValidSharesTable.push_back(allShares);
    }       
}


bool checkUniformity(const std::vector<std::vector<BlnFunction>>& components,
    std::size_t nb_input_variables, const std::size_t level)
{
    const std::size_t input_size = std::pow(2, INPUT_BITS);
    const std::size_t shared_input_size = std::pow(2, INPUT_SHARES * INPUT_BITS);
    const std::size_t expected_count = std::pow(
        2, 2 * nb_input_variables - 2 * level
);
    //int test2;
    std::vector<std::vector<std::size_t>> counts(input_size, std::vector<std::size_t>(OUTPUT_SIZE,0));
    for(std::size_t i = 0; i < input_size; ++i) {
        InputBitArray input(i); 
        for(std::size_t validShare : globalValidSharesTable[i]) {
            std::bitset<OUTPUT_BITS*OUTPUT_SHARES> outputs;
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
    return true;
}



BlnFunction truthTable2(const std::bitset<INPUT_SHARES_BITS>& correction_term, BlnFunction origin) {
    BlnFunction::BitArray test = origin.getTruthTable();
    for (std::size_t i = 0;i<INPUT_SHARES_BITS;i++) {
        if (correction_term[i]) {
            for (std::size_t j = 0;j<INPUT_SIZE;j+std::pow(2,i)) {
                test[j] = test[j]^1;
            }
        }
    }
    return BlnFunction(test);
}

typedef std::vector<VecCorrectionFunction>::const_iterator FunctionIter;

std::vector<VecCorrectionFunction> makeBatchUniformWith(
    const std::string& directory,
    const std::vector<std::vector<BlnFunction>>& realization,
    std::size_t nb_input_variables,
    const Indices& indices,
    FunctionIter begin, FunctionIter end,
    std::size_t batch_nb,
        const std::size_t level) {

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

        if(checkUniformity(corrected_realization, nb_input_variables, level)) {
            ++nb_found;
            good_correction_functions.push_back(*it);
            //std::cout << "Found " << nb_found << " of " << nb_done << " correction functions."
            //   << std::endl;
        }
        ++nb_done;
        if (nb_done % 10000 == 0)
            std::cout << nb_done << std::endl;
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
    
    if (checkUniformity(realization,nb_input_variables,level)) {
        std::cout << "already uniform" << std::endl;
    }
        
    // Produce batches for each set of indices
    for(const auto& pair : functions) {
        std::size_t batch_nb = 0;
        const std::size_t batch_size = pair.second.size() / MAX_THREADS;
        std::vector<std::future<std::vector<VecCorrectionFunction>>> results;
        auto it = pair.second.begin();
        if(1) {
            for(; it < pair.second.end() - batch_size; it += batch_size) {
                results.emplace_back(std::async(std::launch::async,
                    makeBatchUniformWith, directory, realization,
                    nb_input_variables, pair.first, it,
                    it + batch_size, batch_nb, level
                ));
                ++batch_nb;
            }
            // Start batch with remainder
            results.emplace_back(std::async(std::launch::async,
                makeBatchUniformWith, directory, realization,
                nb_input_variables, pair.first, it,
                pair.second.end(), batch_nb, level
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
                    pair.second.end(), batch_nb, level
                );
                if(!corrections.empty() && !filtered_functions.count(pair.first))
                    filtered_functions[pair.first] = std::vector<VecCorrectionFunction>();
                filtered_functions[pair.first].insert(filtered_functions[pair.first].end(),corrections.begin(),corrections.end());
        }
    }
    return filtered_functions;
}

void makeUniform(
    const std::string& directory,
    std::map<Indices, std::vector<VecCorrectionFunction>> functions,
    const std::vector<std::vector<BlnFunction>> realization,
    std::size_t nb_input_variables) {

    std::size_t level = 2;

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
        functions = makeUniformWith(directory, realization, nb_input_variables, candidates, level);

        std::cout << "Found " << countCorrectionFunctions(functions)
                  << " correction functions." << std::endl;
        level++;
    }
    
    std::cout << functions.size() << std::endl;
    writeFunctions(directory, "end.out", functions);
    std::cout << "Process stopped at level " << level << '.' << std::endl;
}

//std::vector<std::vector<BlnFunction>> readRealization(const std::string& filename)
//{
//    std::vector<std::vector<BlnFunction>> result;
//
//    std::ifstream ifs(filename);
//
//    std::string line;
//    std::vector<BlnFunction> shares;
//    std::size_t nb_so_far = 0;
//    while(std::getline(ifs, line)) {
//        if(line.size() != INPUT_SIZE)
//            std::cout << "Warning: incorrect input size given, got " << line.size()
//                      << " expecting " << INPUT_SIZE << '.' << std::endl;
//        shares.push_back(BlnFunction(BlnFunction::BitArray(line)));
//        ++nb_so_far;
//        if(nb_so_far == OUTPUT_SHARES) {
//            result.push_back(shares);
//            shares.clear();
//            nb_so_far = 0;
//        }
//    }
//    if(nb_so_far)
//        std::cout << "Warning: missing shares for last component." << std::endl;
//    std::cout << "Read realization with " << result.size() << " components."
//              << std::endl;
//    return result;
//}

std::vector<std::vector<BlnFunction>> createTruthTable(std::vector<std::vector<std::bitset<INPUT_SHARES_BITS>>>& functions, std::bitset<OUTPUT_SHARES_BITS>& constant_bits, std::vector<std::bitset<INPUT_SHARES_BITS>>& linear_bits) {
    std::vector<std::vector<BlnFunction>> truthTables;
    std::vector<BlnFunction> tempTable;
    for (std::size_t i = 0;i<functions.size();i++) {
        std::bitset<INPUT_SIZE> truthTable;
        for (std::size_t j = 0;j<INPUT_SIZE;j++) {
            std::bitset<INPUT_SHARES*INPUT_BITS> bits(j);
            bool result = 0;
            std::bitset<INPUT_SHARES*INPUT_BITS> im_result;
            for (std::size_t k = 0; k<INPUT_SHARES*INPUT_BITS;k++) {
                im_result[k] = 0;
                for (std::size_t l = 0; l<INPUT_SHARES*INPUT_BITS;l++) {
                    im_result[k] = im_result[k] ^ (functions[i][k][l] ^ bits[l]);
                }
            }
            for (std::size_t k = 0; k<INPUT_SHARES*INPUT_BITS;k++) {
                result ^= bits[k] ^ im_result[k];
            }
            for (std::size_t k = 0; k<INPUT_SHARES*INPUT_BITS;k++) {
                result = result ^ linear_bits[i][k];
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

std::vector<std::vector<BlnFunction>> readRealization(const std::string& filename)
{
    std::vector<std::vector<BlnFunction>> result;

    std::ifstream ifs(filename);

    std::string line;
    std::vector<std::vector<std::bitset<INPUT_SHARES_BITS>>> functions;
    std::size_t nb_so_far = 0;
    std::getline(ifs, line);
    std::bitset<OUTPUT_SHARES_BITS> constant_bits(line);
    std::vector<std::bitset<INPUT_SHARES_BITS>> linear_bits;
    for (std::size_t j=0;j<OUTPUT_SHARES_BITS;j++) {
        std::getline(ifs,line);
        linear_bits.push_back(std::bitset<INPUT_SHARES_BITS>(line));
    }
    std::vector<std::bitset<INPUT_SHARES_BITS>> function;
    while(std::getline(ifs, line)) {
        if(line.size() != INPUT_SHARES_BITS)
            std::cout << "Warning: incorrect input size given, got " << line.size()
                      << " expecting " << INPUT_SHARES_BITS << '.' << std::endl;
        function.push_back(std::bitset<INPUT_SHARES_BITS>(line));
        nb_so_far++;
        if(nb_so_far == INPUT_SHARES_BITS) {
            functions.push_back(function);
            function.clear();
            nb_so_far = 0;
        }
    }
    if(nb_so_far)
        std::cout << "Warning: missing shares for last component." << std::endl;
    result = createTruthTable(functions, constant_bits, linear_bits);
    std::cout << "Read realization with " << result.size() << " components."
              << std::endl;
    return result;
}

//bool getANF(const int f[N], int p[N]){
//    int i, j;
//    for (i = 0; i < N; ++i) p[i] = f[i];
//
//    for (i = 1; i < N; i <<= 1){
//        for (j = i; j < N; j = (j + 1) | i){
//            p[j] ^= p[j ^ i];
//        }
//    }
//    return true;
//}

int main(int argc, char *argv[])
{
    
//    if(argc < 4) {
//        std::cerr << "Usage: <number of inputs> <output directory> <filename>."
//                  << std::endl;
//        return 1;
//    }
//    std::size_t nb_inputs = std::stoi(argv[1]);
//    if(nb_inputs != INPUT_BITS) {
//        std::cerr << "Expected exactly " << INPUT_BITS << " input bits."
//                  << std::endl;
//        return 1;
//    }
//    
//    std::cout << "Reading realization..." << std::endl;
    auto realization = readRealization("new_format_perm.in");
    std::size_t nb_inputs = INPUT_BITS;
    createGlobalSharingTable();
    makeUniform(
        "gf4", buildCorrectionTerms(realization),
        realization, nb_inputs
    ); 
    
    return 0;
}
