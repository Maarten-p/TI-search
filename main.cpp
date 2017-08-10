#include <cinttypes>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <future>
#include <map>
#include <tuple>
#include "BooleanFunction.h" 

constexpr std::size_t MAX_THREADS = 8;
constexpr std::size_t MAX_HAMMING_WEIGHT = 10000;
constexpr std::size_t MAX_HAMMING_WEIGHT_TOTAL = 1001000;
constexpr std::size_t INPUT_BITS = 4;
constexpr std::size_t OUTPUT_BITS = 4;
constexpr std::size_t INPUT_SHARES = 3;
constexpr std::size_t OUTPUT_SHARES = 3;
constexpr std::size_t INPUT_SHARES_BITS = INPUT_BITS * INPUT_SHARES;
constexpr std::size_t OUTPUT_SHARES_BITS = OUTPUT_BITS * OUTPUT_SHARES;
constexpr std::size_t INPUT_SHARES_SIZE = std::pow(2, INPUT_BITS * INPUT_SHARES);
constexpr std::size_t OUTPUT_SHARES_SIZE = std::pow(2, OUTPUT_BITS * OUTPUT_SHARES);
constexpr std::size_t INPUT_SIZE = std::pow(2, INPUT_BITS);
constexpr std::size_t OUTPUT_SIZE = std::pow(2, OUTPUT_BITS);
constexpr bool CONTROL = false;

typedef BooleanFunction<INPUT_SHARES_SIZE> BlnFunction;
typedef std::array<std::bitset<INPUT_SHARES_BITS>, OUTPUT_SHARES> CorrectionFunction;
typedef std::vector<CorrectionFunction> VecCorrectionFunction;
typedef std::bitset<INPUT_BITS> InputBitArray;
typedef std::bitset<INPUT_SHARES_BITS> SharedInputBitArray;

std::vector<std::bitset<OUTPUT_SHARES_SIZE>> globalTruthTable(INPUT_SHARES_SIZE,0);
std::vector<std::vector<std::size_t>> globalValidSharesTable;
std::vector<std::bitset<INPUT_SHARES_BITS>> globalDependence;

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
bool areSharesZero(const std::bitset<INPUT_SHARES_BITS>& bits, std::size_t share, int outputIndex) {
    for(std::size_t i = share; i < bits.size(); i += INPUT_SHARES) {
        //std::cout<<globalDependence[outputIndex]<<std::endl;
        //std::cout<<i<<std::endl;
        if(bits[i]) 
            return false;
        if(globalDependence[outputIndex][i])
            return false;
    }
    return true;
}

bool hammingWeightConstraint(const std::bitset<INPUT_SHARES_BITS> bits) {
    return bits.count() <= MAX_HAMMING_WEIGHT;
}

bool totalHammingWeightConstraint(const std::bitset<INPUT_SHARES_BITS> bits1, const std::bitset<INPUT_SHARES_BITS> bits2, const std::bitset<INPUT_SHARES_BITS> bits3) {
    return (bits1.count()+bits2.count()+bits3.count()) <= MAX_HAMMING_WEIGHT_TOTAL;
}

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

std::vector<std::vector<int>> constructQuadraticCorrectionTermsSpace() {
    std::vector<std::vector<int>> quadraticSpace;
    for(std::size_t i = 0; i < INPUT_SHARES; ++i) {
        std::vector<int> oneShareSpace;
        for(std::size_t j = 0; j < INPUT_SHARES_BITS; ++j) {
            for(std::size_t k = j; k < INPUT_SHARES_BITS; ++k) {
                if (j/INPUT_SHARES == k/INPUT_SHARES || (j%INPUT_SHARES==i) || (k%INPUT_SHARES==i))
                    continue;
                oneShareSpace.push_back(j*INPUT_SHARES_BITS+k);
            }
        }
        quadraticSpace.push_back(oneShareSpace);
    }
    return quadraticSpace;
}

BlnFunction addQuadraticTerm(BlnFunction& blnFunction,std::size_t j,std::size_t k) {
    std::bitset<INPUT_SHARES_SIZE> truthTable = blnFunction.getTruthTable();
    for(std::size_t i = 0; i < INPUT_SHARES_SIZE; ++i) {
        InputBitArray input(i);
        truthTable[i] = truthTable[i] ^ (input[j] & input[k]);
    }
    return BlnFunction(truthTable);
}

void tryQuadraticCorrectionTerms(std::vector<std::vector<BlnFunction>>& realization,std::size_t outputIndex, std::vector<std::vector<int>> quadraticSpace) {
    std::bitset<INPUT_SHARES> problematicFunctions;
    problematicFunctions.flip();
    for(std::size_t i = 0; i < INPUT_SHARES; ++i) {
        auto spectrum1 = realization[outputIndex][i].walshHadamardTransform(); 
        for(std::size_t j = 0; j < INPUT_SHARES_SIZE; ++j) {
            std::bitset<INPUT_SHARES_BITS> bits_j(j);
            if(!areSharesZero(bits_j, i, outputIndex*3+i) || spectrum1[j] != 0 || !hammingWeightConstraint(bits_j))
                continue;
            problematicFunctions[i] = 0;
            break;
        }
    }
    for(std::size_t i = 0; i < INPUT_SHARES; ++i) {
        if (!problematicFunctions[i])
            continue;
        bool solutionFound = false;
        
        for(std::size_t j = 0; j < INPUT_SHARES_BITS; ++j) {
            for(std::size_t k = j; k < INPUT_SHARES_BITS; ++k) {
                if (j/INPUT_SHARES == k/INPUT_SHARES || (j%INPUT_SHARES==i) || (k%INPUT_SHARES==i))
                    continue;
                realization[outputIndex][i] = addQuadraticTerm(realization[outputIndex][i],j,k);
                std::cout << realization[outputIndex][i] << std::endl;
                auto spectrum1 = realization[outputIndex][i].walshHadamardTransform();
                for(std::size_t l = 0; l < INPUT_SHARES_SIZE; ++l) {
                    std::bitset<INPUT_SHARES_BITS> bits_l(l);
                    if(!areSharesZero(bits_l, i, outputIndex*3+i) || spectrum1[l] != 0 || !hammingWeightConstraint(bits_l))
                        continue;
                    solutionFound = true;
                    std::ofstream ofs("output/quadratic.out");
                    ofs << j << " " << k << std::endl;
                    std::cout << j << " " << k << std::endl;
                    break;
                }
                if (solutionFound) {
                    break;
                }
                else {
                    //adding the same term again removes it
                    addQuadraticTerm(realization[outputIndex][i],j,k);
                }
            }
        if (solutionFound)
            break;
        }
        if (!solutionFound)
            exit(EXIT_FAILURE);
    }
   
}
std::map<Indices, std::vector<VecCorrectionFunction>> buildCorrectionTerms(
    std::vector<std::vector<BlnFunction>>& realization)
{
    std::cout << "Collecting linear correction terms." << std::endl;
    std::map<Indices, std::vector<VecCorrectionFunction>> correctionTerms;
    bool quadraticSpaceExists = false;
    std::vector<std::vector<int>> quadraticSpace;
    for(std::size_t i = 0; i < realization.size(); ++i) {
        Indices ind;
        ind.indices.push_back(i);
        correctionTerms[ind] = getLinearCorrections(
            realization[i][0],
            realization[i][1],
            realization[i][2], i
        );
        if (correctionTerms[ind].size() == 0) {
            if (!quadraticSpaceExists)
                quadraticSpace = constructQuadraticCorrectionTermsSpace();
            std::cout << "test" << std::endl;
            tryQuadraticCorrectionTerms(realization,i, quadraticSpace);
            correctionTerms[ind] = getLinearCorrections(
                realization[i][0],
                realization[i][1],
                realization[i][2], i
            );
        }
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
                    //int v1 = rand() % 30; 
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


bool checkUniformity(const std::vector<std::vector<BlnFunction>>& components,
    std::size_t nb_input_variables) {
    const std::size_t expected_count = std::pow(
        2, 2 * nb_input_variables - 2 * components.size()
);
    //int test2;
    std::vector<std::vector<std::size_t>> counts(INPUT_SIZE, std::vector<std::size_t>(OUTPUT_SHARES_SIZE,0));
    for(std::size_t i = 0; i < INPUT_SIZE; ++i) {
        InputBitArray input(i); 
        for(std::size_t validShare : globalValidSharesTable[i]) {
            std::bitset<OUTPUT_SHARES_BITS> outputs;
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

        if(checkUniformity(corrected_realization, nb_input_variables)) {
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

void makeUniform(
    const std::string& filename,
    std::map<Indices, std::vector<VecCorrectionFunction>> functions,
    const std::vector<std::vector<BlnFunction>> realization,
    std::size_t nb_input_variables) {

    std::size_t level = 2;
    if (checkUniformity(realization,nb_input_variables)) {
        std::cout << "already uniform" << std::endl;
        return;
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
    
    std::cout << functions.size() << std::endl;
    writeFunctions("output", filename, functions);
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

std::vector<std::bitset<INPUT_SHARES_BITS>> readCandidate(const std::string& filename) {
    
    std::ifstream ifs(filename);
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
        std::vector<std::bitset<INPUT_SHARES_BITS>> correctionTerm = readCandidate("control/correction.in");
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
        std::cout << checkUniformity(result,INPUT_BITS) << std::endl;
    }
    else {
    auto result = createTruthTable(std::get<0>(sharedTuple),std::get<1>(sharedTuple),std::get<2>(sharedTuple));
    getDependence(std::get<1>(sharedTuple),std::get<2>(sharedTuple));
    std::size_t nb_inputs = INPUT_BITS;
    createGlobalSharingTable();
    makeUniform(
        "output", buildCorrectionTerms(result),
        result, nb_inputs
    ); 
    }
    return 0;
}
