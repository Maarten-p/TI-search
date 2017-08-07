#include <cstdlib>

bool getANF(const int f[N], int p[N]){
    int i, j;
    for (i = 0; i < N; ++i) p[i] = f[i];

    for (i = 1; i < N; i <<= 1){
        for (j = i; j < N; j = (j + 1) | i){
            p[j] ^= p[j ^ i];
        }
    }
    return true;
}

bool ANFToProperInput(const int p[N], const std::string& filename) {
    
}

int main(int argc, char** argv) {
    int N = argv[0]
    int p[N] = {};
    getANF(argv[1],p);
    return 0;
}

