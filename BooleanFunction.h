#ifndef _BOOLEANFUNCTION_H_INCLUDE_GUARD
#define _BOOLEANFUNCTION_H_INCLUDE_GUARD

#include <vector>
#include <bitset>
#include <cmath>
#include <ostream>

template<std::size_t NIn>
class BooleanFunction {
public:
    typedef std::bitset<NIn> BitArray;
private:

    BitArray truthTable;

    typedef std::vector<int>::iterator Iterator;
    void discreteFourierTransform(Iterator begin, Iterator end) const {
        const std::size_t middle = std::distance(begin, end) / 2;
        if(middle == 0)
            return;
        std::vector<int> upperHalf(
            begin,
            begin + middle 
        );
        auto itUpper = begin;
        auto itLower = begin + middle;
        for(std::size_t i = 0; i < middle; ++i) {
            *itUpper += *itLower;
            *itLower = upperHalf[i] - *itLower; 
            ++itUpper;
            ++itLower;
        }

        discreteFourierTransform(begin, begin + middle);
        discreteFourierTransform(begin + middle, end);
    }

public:
    explicit BooleanFunction(const BitArray& truthTable)
        : truthTable(truthTable)
    {
    }

    bool operator[](const std::size_t& index) const
    {
        return truthTable[index];
    }

    BooleanFunction& operator+=(const BooleanFunction& rhs)
    {
        truthTable ^= rhs.truthTable;
    }

    std::vector<int> walshHadamardTransform() const
    {
        std::vector<int> expandedTruthTable;
        expandedTruthTable.resize(NIn);
        for(std::size_t i = 0; i < truthTable.size(); ++i)
            expandedTruthTable[i] = truthTable[i];
        discreteFourierTransform(expandedTruthTable.begin(), expandedTruthTable.end());

        for(int& i : expandedTruthTable)
            i *= -2;

        expandedTruthTable[0] += NIn;

        return expandedTruthTable;
    }

    bool operator==(const BooleanFunction<NIn>& rhs) const {
        return truthTable == rhs.truthTable;
    }

    BitArray getTruthTable() const {
        return truthTable;
    }
};

template<std::size_t NIn>
std::ostream& operator<<(std::ostream& os, const BooleanFunction<NIn>& rhs)
{
    typename BooleanFunction<NIn>::BitArray tt = rhs.getTruthTable();
    for(std::size_t i = 0; i < tt.size(); ++i)
        os << tt[i] << ' ';
    return os;
}

template<std::size_t NIn>
BooleanFunction<NIn> operator+(BooleanFunction<NIn> lhs,
    const BooleanFunction<NIn>& rhs)
{
    lhs += rhs;
    return lhs;
}

template<std::size_t NIn>
bool operator<(const BooleanFunction<NIn>& b1, const BooleanFunction<NIn>& b2) {
    return b1.BitArray.to_ulong() < b2.BitArray.to_ulong();
}

#endif // _BOOLEANFUNCTION_H_INCLUDE_GUARD
