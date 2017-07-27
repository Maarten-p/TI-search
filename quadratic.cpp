#include <NTL/mat_GF2.h>
#include <iostream>
#include <vector>
#include <set>
#include <utility>
#include <thread>

typedef std::set<std::pair<long, long>> PairSet;

struct Position {
    Position(short i, short j)
        : i(i), j(j) {}
    short i;
    short j;
};

typedef std::vector<Position> SparseMatrix;

void setBlock(NTL::mat_GF2& matrix, long nb_shares, long i, long j)
{
    for(long s = 0; s < nb_shares - 1; ++s) {
        for(long t = 0; t < nb_shares - 1; ++t) {
            if(s != t) {
                matrix.put(i + s, j + t, NTL::GF2(1)); 
                matrix.put(j + t, i + s, NTL::GF2(1)); 
            }
        }
    }
}

void setBlockFull(NTL::mat_GF2& matrix, long nb_shares, long i, long j)
{
    for(long s = 0; s < nb_shares; ++s) {
        for(long t = 0; t < nb_shares; ++t) {
            matrix.put(i + s, j + t, NTL::GF2(1)); 
            matrix.put(j + t, i + s, NTL::GF2(1)); 
        }
    }
}

NTL::mat_GF2 getInitialMatrix(long nb_inputs, long nb_shares, const PairSet& terms)
{
    const long size = nb_inputs * (nb_shares - 1);
    NTL::mat_GF2 matrix;
    matrix.SetDims(size, size);
    for(long k = 0; k < nb_inputs; ++k) {
        const long i = (nb_shares - 1) * k;
        for(long l = 0; l < k; ++l) {
            const long j = (nb_shares - 1) * l;
            if(terms.count(std::make_pair(k, l)) || terms.count(std::make_pair(l, k)))
                setBlock(matrix, nb_shares, i, j);
        }
    }
    return matrix;
}

NTL::mat_GF2 getSumMatrix(long nb_inputs, long nb_shares, const PairSet& terms)
{
    const long size = nb_inputs * nb_shares;
    NTL::mat_GF2 matrix;
    matrix.SetDims(size, size);
    for(long k = 0; k < nb_inputs; ++k) {
        const long i = nb_shares * k;
        for(long l = 0; l < k; ++l) {
            const long j = nb_shares * l;
            if(terms.count(std::make_pair(k, l)) || terms.count(std::make_pair(l, k)))
                setBlockFull(matrix, nb_shares, i, j);
        }
    }
    return matrix;
}

void getLowRankMatricesImpl(const NTL::mat_GF2& base, const SparseMatrix& current,
    long max_rank, long max_terms, std::vector<SparseMatrix>& matrices)
{
    if(max_terms == 0)
        return;

    for(long i = 0; i < base.NumRows(); ++i) {
        for(long j = i; j < base.NumCols(); ++j) {
            if(i == j || base.get(i, j) == NTL::GF2(1))
                continue;

            NTL::mat_GF2 copy(base);
            copy.put(i, j, NTL::GF2(1));
            copy.put(j, i, NTL::GF2(1));
            SparseMatrix copy_sparse = current;
            copy_sparse.emplace_back(i, j);

            NTL::mat_GF2 temp(copy);
            if(NTL::gauss(temp/*, max_rank + 1*/) <= max_rank)
                matrices.push_back(copy_sparse); 
            getLowRankMatricesImpl(copy, copy_sparse, max_rank, max_terms - 1, matrices);
        }
    }
}

/**
 * @pre base is a symmetric matrix
 */
void getLowRankMatrices(const NTL::mat_GF2& base, long max_rank, long max_terms,
    std::vector<SparseMatrix>& matrices)
{
    const SparseMatrix sparse;
    getLowRankMatricesImpl(base, sparse, max_rank, max_terms, matrices);
}

void copyBlock(const NTL::mat_GF2& matrix, NTL::mat_GF2& result, long k, long l,
    long nb_inputs, long nb_shares, long exclude_share)
{
    const long i = nb_shares * k;
    const long j = nb_shares * l;

    const long a = (nb_shares - 1) * k;
    const long b = (nb_shares - 1) * l;

    long offset_s = 0;
    for(long s = 0; s < nb_shares; ++s) {
        if(s == exclude_share) {
            offset_s = 1;
            continue;
        }
        long offset_t = 0;
        for(long t = 0; t < nb_shares; ++t) {
            if(t == exclude_share) {
                offset_t = 1;
                continue;
            }
            NTL::GF2 value = matrix.get(a + s - offset_s, b + t - offset_t);
            result.put(i + s, j + t, value); 
            result.put(j + t, i + s, value); 
        }
    }
}

NTL::mat_GF2 expandMatrix(const NTL::mat_GF2& matrix, long nb_inputs, long nb_shares,
    long exclude_share)
{
    const long size = nb_inputs * nb_shares;
    NTL::mat_GF2 result;
    result.SetDims(size, size);
    for(long k = 0; k < nb_inputs; ++k) {
        for(long l = 0; l < k; ++l)
            copyBlock(matrix, result, k, l, nb_inputs, nb_shares, exclude_share);
    }
    return result;
}

bool validateMatrix(const NTL::mat_GF2& matrix, long nb_shares, long exclude_share)
{
    for(long i = 0; i < matrix.NumRows(); ++i) {
        for(long j = 0; j < matrix.NumCols(); ++j) {
            if(matrix.get(i, j) == NTL::GF2(1) &&
              (i % nb_shares == exclude_share || j % nb_shares == exclude_share))
                return false;
        }
    }
    return true;
}

NTL::mat_GF2 transform(NTL::mat_GF2 base, const SparseMatrix& mod)
{
    for(const Position& pos : mod) {
        base.put(pos.i, pos.j, NTL::GF2(1));
        base.put(pos.j, pos.i, NTL::GF2(1));
    }
    return base;
}

typedef std::vector<SparseMatrix>::const_iterator MatIter;

void printSolution(MatIter begin, MatIter end, const std::vector<SparseMatrix>& matrices,
    long nb_inputs, long nb_shares, long max_rank, const NTL::mat_GF2& base, const NTL::mat_GF2& sum)
{
    for(auto it = begin; it != end; ++it) {
        NTL::mat_GF2 m1_exp = expandMatrix(transform(base, *it), nb_inputs, nb_shares, 0); // For share 0
        for(const SparseMatrix& m2_sparse : matrices) {
            NTL::mat_GF2 m2 = transform(base, m2_sparse);
            NTL::mat_GF2 m2_exp = expandMatrix(m2, nb_inputs, nb_shares, 1); // For share 1
            NTL::mat_GF2 m3 = m1_exp + m2_exp + sum;
            if(validateMatrix(m3, nb_shares, 2) && NTL::gauss(m3/*, max_rank + 1*/) <= max_rank) {
                std::cout << "Found matrices: " << std::endl;
                std::cout << m1_exp << std::endl;
                std::cout << m2_exp << std::endl;
                std::cout << m1_exp + m2_exp + sum << std::endl;
                return;
            }
        }
    }
} 

int main()
{
    PairSet pairs;
    pairs.insert(std::make_pair(0l, 2l));
    pairs.insert(std::make_pair(0l, 3l));
    pairs.insert(std::make_pair(1l, 2l));

    const long max_rank = 6;
    const long max_rank_third = 6;
    const long max_terms = 3;
    const long nb_inputs = 4;
    const long nb_shares = 3;
    const long nb_threads = 4;
    
    // PairSet pairs;
    // pairs.insert(std::make_pair(0l, 5l));
    // pairs.insert(std::make_pair(0l, 7l));
    // pairs.insert(std::make_pair(1l, 4l));
    // pairs.insert(std::make_pair(1l, 5l));
    // pairs.insert(std::make_pair(1l, 6l));
    // pairs.insert(std::make_pair(1l, 7l));
    // pairs.insert(std::make_pair(2l, 5l));
    // pairs.insert(std::make_pair(2l, 6l));
    // pairs.insert(std::make_pair(3l, 4l));
    // pairs.insert(std::make_pair(3l, 5l));
    // pairs.insert(std::make_pair(3l, 7l));

    // const long max_rank = 10;
    // const long max_rank_third = 14;
    // const long max_terms = 4;
    // const long nb_inputs = 8;
    // const long nb_shares = 3;
    // const long nb_threads = 24;
    

    const NTL::mat_GF2 base = getInitialMatrix(nb_inputs, nb_shares, pairs);
    const NTL::mat_GF2 sum = getSumMatrix(nb_inputs, nb_shares, pairs);
    std::vector<SparseMatrix> matrices;
    getLowRankMatrices(base, max_rank, max_terms, matrices);

    const size_t nb_per_thread = matrices.size() / nb_threads;

    std::cout << "Found " << matrices.size() << " low rank matrices." << std::endl;
    std::cout << "Starting " << nb_threads << " threads working on "
              << nb_per_thread << " matrices each." << std::endl;

    std::vector<std::thread> threads;
    auto it = matrices.begin();
    const auto end = matrices.end() - nb_per_thread;
    for(; it < end; it += nb_per_thread) {
        // Note: thread safety is ok since printSolutions doesn't modify the
        //       matrices or arrays
        threads.emplace_back(
            printSolution, it, it + nb_per_thread, std::cref(matrices),
            nb_inputs, nb_shares, max_rank_third, std::cref(base), std::cref(sum)
        );
    }
    // Add the remaining matrices
    threads.emplace_back(
        printSolution, it, matrices.end(), std::cref(matrices),
        nb_inputs, nb_shares, max_rank, std::cref(base), std::cref(sum)
    );

    for(auto& t : threads)
        t.join();
     
}
