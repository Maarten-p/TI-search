#include <cinttypes>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <future>
#include <map>
#include <tuple>
#include <bitset>
constexpr int INPUT_BITS = 4;
constexpr int OUTPUT_BITS = 4;
constexpr int INPUT_SHARES = 3;
constexpr int INPUT_SHARED_BITS = INPUT_BITS*INPUT_SHARES;
constexpr int INPUT_SIZE = std::pow(2,INPUT_BITS);
constexpr int OUTPUT_SIZE = std::pow(2,OUTPUT_BITS);

bool getANF(const int f[INPUT_SIZE], int p[INPUT_SIZE]){
    int i, j;
    for (i = 0; i < INPUT_SIZE; ++i) p[i] = f[i];

    for (i = 1; i < INPUT_SIZE; i <<= 1){
        for (j = i; j < INPUT_SIZE; j = (j + 1) | i){
            p[j] ^= p[j ^ i];
        }
    }
    return true;
}

bool unsharedToSharedInput(std::bitset<INPUT_BITS> constant_bits_old, std::vector<std::bitset<INPUT_BITS>> linear_bits_old, std::vector<std::vector<std::bitset<INPUT_BITS>>> quadratic_bits_old) {
    std::bitset<INPUT_SHARED_BITS> constant_bits;
    std::vector<std::bitset<INPUT_SHARED_BITS>> linear_bits;
    std::vector<std::vector<std::bitset<INPUT_SHARED_BITS>>> quadratic_bits;
    for(std::size_t i=0;i<constant_bits.size();i++) {
        constant_bits[3*i] = constant_bits_old[i];
        constant_bits[3*i+1] = constant_bits_old[i];
        constant_bits[3*i+2] = constant_bits_old[i];
    }
    for(std::size_t i=0;i<linear_bits_old.size();i++) {
    }
}

bool ANFToUnsharedInput(const int p[INPUT_SIZE]) {
    
    std::vector<std::vector<std::size_t>> counts(INPUT_SIZE, std::vector<std::size_t>(OUTPUT_SIZE,0));
    std::vector<std::vector<std::bitset<INPUT_BITS>>> quadratic_bits(INPUT_SIZE, std::vector<std::INPUT_SIZE>(INPUT_SIZE,std::bitset<INPUT_BITS));
    std::bitset<OUTPUT_BITS> constant_bits;
    std::vector<std::bitset<INPUT_BITS>> linear_bits(INPUT_BITS,std::bitset<INPUT_BITS>);
    for (int i=0;i<OUTPUT_BITS;i++) {
        constant_bits[i] = p[0][i];
    }
    int next_single_term = 1;
    int number_of_single_terms = 1;
    for (int i=1;i<INPUT_SIZE;i++) {
        if (i==next_single_term) {
            for (int j=0;j<INPUT_BITS;j++) {
                linear_bits[j][i] = p[i][j];
            }
            next_single_term += number_of_single_terms;
            number_of_single_terms++;
        }
        else {
            for (int j=0;j<INPUT_BITS;j++) {
                quadratic_bits[j][i+number_of_single_terms-1-next_single_term][number_of_single_terms-1] = p[i][j];
            }
        }
    }
    unsharedToSharedInput(constant_bits,linear_bits,quadratic_bits); 
}


bool writeToFile(std::bitset<INPUT_SHARED_BITS> constant_bits, std::vector<std::bitset<INPUT_SHARED_BITS>> linear_bits, std::vector<std::vector<std::bitset<INPUT_SHARED_BITS>>> quadratic_bits) {
    return 0;
}

int readTruthTable(const std::string& filename) {
    std::array<int,INPUT_SIZE> truthTable;
    std::ifstream ifs(filename);
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
int main(int argc, char** argv) {
    std::array<int,INPUT_SIZE> ANFForm;
    int p[INPUT_SIZE] = {};
    getANF(readTruthTable("raw_input/test.in"),p);
    return 0;
}

