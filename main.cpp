#include <cinttypes>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <future>
#include <map>
#include <tuple>

#include "BooleanFunction.h" 

constexpr std::size_t MAX_THREADS = 4;

constexpr std::size_t INPUT_BITS = 4;
constexpr std::size_t NB_SHARES = 3;
constexpr std::size_t INPUT_SIZE = std::pow(2, INPUT_BITS * NB_SHARES);

typedef BooleanFunction<INPUT_SIZE> BlnFunction;
typedef std::array<BlnFunction, NB_SHARES> CorrectionFunction;
typedef std::vector<CorrectionFunction> VecCorrectionFunction;
typedef std::bitset<INPUT_BITS> InputBitArray;
typedef std::bitset<NB_SHARES * INPUT_BITS> SharedInputBitArray;

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
bool areSharesZero(const BlnFunction::BitArray& bits, std::size_t share) {
    for(std::size_t i = share; i < bits.size(); i += NB_SHARES) {
        if(bits[i]) 
            return false;
    }
    return true;
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
        BlnFunction::BitArray bits_i(i);
        if(!areSharesZero(bits_i, 0) || spectrum1[i] != 0)
            continue;

        for(std::size_t j = 0; j < spectrum2.size(); ++j) {
            BlnFunction::BitArray bits_j(j);
            if(!areSharesZero(bits_j, 1) || spectrum2[j] != 0)
                continue;

            BlnFunction::BitArray bits_l = bits_i ^ bits_j;
            std::size_t l = bits_l.to_ulong();
            if(!areSharesZero(bits_l, 2) || spectrum3[l] != 0)
                continue;

            CorrectionFunction cf {
                BlnFunction(bits_i),
                BlnFunction(bits_j),
                BlnFunction(bits_l)
            };

            // Found a solution, add to list
            solutions.push_back({cf});
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
        for(const auto& correction1 : functions.at(ind_no_first)) {
            for(const auto& correction2 : functions.at(ind_no_second)) {
                VecCorrectionFunction correction = correction1;
                correction.insert(
                    correction.end(), correction2.begin(), correction2.end()
                );
                if(shouldAddCorrection(functions, ind, correction)) {
                    if(!result.count(ind))
                        result[ind] = std::vector<VecCorrectionFunction>();
                    result[ind].push_back(correction);
                }
            }
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
                for(const BlnFunction& share : f)
                    ofs << '\t' << share << '\n';
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

std::vector<std::bitset<3>> get3Sharing(bool bit)
{
    std::vector<std::bitset<3>> result;
    std::bitset<3> current;
    for(bool b1 = false; !b1; b1 = true) {
        current[0] = b1;
        for(bool b2 = false; !b2; b2 = true) {
            current[1] = b2;
            current[2] = b1 ^ b2;
            result.push_back(current);
        }
    }
}

std::vector<SharedInputBitArray> getSharingsForInput(const InputBitArray& input_bits)
{
    std::vector<SharedInputBitArray> result;
    for(std::size_t i = 0; i < input_bits.size(); ++i) {
        const std::size_t k = 3 * i;
        for(const auto& bits : get3Sharing(input_bits[i])) {
            for(SharedInputBitArray& r : result) {
                r[k] = bits[0];
                r[k + 1] = bits[1];
                r[k + 2] = bits[2];
            }
        }
    }
    return result;
}

bool checkUniformity(const std::vector<BlnFunction>& components,
    std::size_t nb_input_variables, std::size_t expected_count)
{
    const std::size_t input_size = std::pow(2, nb_input_variables);
    const std::size_t shared_input_size = std::pow(2, 3 * nb_input_variables);
    for(std::size_t i = 0; i < input_size; ++i) {
        InputBitArray input(i); 
        std::vector<std::size_t> counts(shared_input_size, 0);
        for(auto& sharing : getSharingsForInput(input)) {
            std::size_t sharingIndex = sharing.to_ulong();
            counts[sharingIndex] += 1; 
            if(counts[sharingIndex] > expected_count)
                return false;
        }
    }
    return true;
}

typedef std::vector<VecCorrectionFunction>::const_iterator FunctionIter;

std::vector<VecCorrectionFunction> makeBatchUniformWith(
    const std::string& directory,
    const std::vector<std::vector<BlnFunction>>& realization,
    std::size_t nb_input_variables,
    const Indices& indices,
    FunctionIter begin, FunctionIter end,
    std::size_t batch_nb) {

    std::vector<VecCorrectionFunction> good_correction_functions;

    std::ofstream os(directory + "/batch-" + std::to_string(batch_nb) + ".log");
    
    std::size_t nb_done = 0;
    std::size_t nb_found = 0;

    const std::size_t expected_count = std::pow(
        2, 2 * nb_input_variables - 2 * realization.size()
    );

    // Try all possible correction functions
    for(FunctionIter it = begin; it != end; ++it) {
        // Construct the vectorial boolean function, i.e. add vf to realization
        // The result is a flattened vector of vectorial boolean functions
        std::vector<BlnFunction> corrected_realization;
        std::size_t i_corrected = 0;
        for(std::size_t i : indices.indices) {
            for(std::size_t s = 0; s < realization.size(); ++s)
                corrected_realization.push_back(realization[i][s] + (*it)[i_corrected][s]);
            ++i_corrected;
        }

        if(checkUniformity(corrected_realization, nb_input_variables, expected_count)) {
            ++nb_found;
            good_correction_functions.push_back(*it);
            os << "Found " << nb_found << " of " << nb_done << " correction functions."
               << std::endl;
        }
        ++nb_done;
    }

}

std::map<Indices, std::vector<VecCorrectionFunction>> makeUniformWith(
    const std::string& directory,
    const std::vector<std::vector<BlnFunction>>& realization,
    std::size_t nb_input_variables,
    const std::map<Indices, std::vector<VecCorrectionFunction>>& functions) {

    std::map<Indices, std::vector<VecCorrectionFunction>> filtered_functions;

    // Produce batches for each set of indices
    for(const auto& pair : functions) {
        std::size_t batch_nb = 0;
        const std::size_t batch_size = pair.second.size() / MAX_THREADS;
        std::vector<std::future<std::vector<VecCorrectionFunction>>> results;
        auto it = pair.second.begin();
        for(; it < pair.second.end() - batch_size; it += batch_size) {
            results.emplace_back(std::async(std::launch::async,
                makeBatchUniformWith, directory, realization,
                nb_input_variables, pair.first, it,
                it + batch_size, batch_nb
            ));
            ++batch_nb;
            std::cout << batch_nb << std::endl;
        }
        // Start batch with remainder
        results.emplace_back(std::async(std::launch::async,
            makeBatchUniformWith, directory, realization,
            nb_input_variables, pair.first, it,
            pair.second.end(), batch_nb
        ));

        // Retrieve the results
        int i = 0;
        for(auto& future : results) {
            i +=1;
            std::cout << i << std::endl;
            std::vector<VecCorrectionFunction> corrections = future.get();
            if(!corrections.empty() && !filtered_functions.count(pair.first))
                filtered_functions[pair.first] = std::vector<VecCorrectionFunction>();
            std::vector<VecCorrectionFunction>& results_so_far = filtered_functions[pair.first];
            results_so_far.insert(results_so_far.begin(), corrections.begin(), corrections.end());
        }
    }
    return filtered_functions;
}

void makeUniform(
    const std::string& directory,
    std::map<Indices, std::vector<VecCorrectionFunction>> functions,
    const std::vector<std::vector<BlnFunction>>& realization,
    std::size_t nb_input_variables) {

    int level = 2;
    writeFunctions(directory, "level1.out", functions);

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
        functions = makeUniformWith(directory, realization, nb_input_variables, candidates);

        std::cout << "Found " << countCorrectionFunctions(candidates)
                  << " correction functions." << std::endl;
    }

    std::cout << "Process stopped at level " << level << '.' << std::endl;
}

std::vector<std::vector<BlnFunction>> readRealization(const std::string& filename)
{
    std::vector<std::vector<BlnFunction>> result;

    std::ifstream ifs(filename);

    std::string line;
    std::vector<BlnFunction> shares;
    std::size_t nb_so_far = 0;
    while(std::getline(ifs, line)) {
        if(line.size() != INPUT_SIZE)
            std::cout << "Warning: incorrect input size given, got " << line.size()
                      << " expecting " << INPUT_SIZE << '.' << std::endl;
        shares.push_back(BlnFunction(BlnFunction::BitArray(line)));
        ++nb_so_far;
        if(nb_so_far == NB_SHARES) {
            result.push_back(shares);
            shares.clear();
            nb_so_far = 0;
        }
    }
    if(nb_so_far)
        std::cout << "Warning: missing shares for last component." << std::endl;
    std::cout << "Read realization with " << result.size() << " components."
              << std::endl;
    return result;
}


int main(int argc, char *argv[])
{
    if(argc < 4) {
        std::cerr << "Usage: <number of inputs> <output directory> <filename>."
                  << std::endl;
        return 1;
    }
    std::size_t nb_inputs = std::stoi(argv[1]);
    if(nb_inputs != INPUT_BITS) {
        std::cerr << "Expected exactly " << INPUT_BITS << " input bits."
                  << std::endl;
        return 1;
    }
    
    std::cout << "Reading realization..." << std::endl;
    auto realization = readRealization(argv[3]);
    makeUniform(
        argv[2], buildCorrectionTerms(realization),
        realization, nb_inputs 
    ); 
    
    return 0;
}
