#include <libff/common/profiling.hpp>

#include <libff/algebra/curves/bn128/bn128_pp.hpp>
#include <libff/algebra/curves/mnt/mnt4/mnt4_pp.hpp>
#include <libff/algebra/curves/mnt/mnt6/mnt6_pp.hpp>

#include <vector>
#include <chrono>
#include <fstream>

#include <gmp.h>
#include <gmpxx.h>


using namespace libff;
using namespace std;

template <typename ppT>
void expo_test2(const vector<unsigned long> delta_vec){

    #ifdef CURVE_MNT6
    cout << "MNT6 checked" << endl;
    #endif
    cout << "Size of the vector: " << delta_vec.size() << endl;
    auto t2 = chrono::high_resolution_clock::now();
    auto t3 = chrono::high_resolution_clock::now();
    auto t4 = t2 - t2;

    Fr<ppT> temp = Fr<ppT>::random_element();
    G1<ppT> g1 = G1<ppT>::one();
    G1<ppT> temp_g1 = (Fr<ppT>::random_element()) * G1<ppT>::one();

    t4 = t2 - t2;
    for (size_t i = 0; i < delta_vec.size(); ++i){

        temp = delta_vec[i];
        t2 = chrono::high_resolution_clock::now();
        temp_g1 = temp_g1 + (temp * g1);
        t3 = chrono::high_resolution_clock::now();
        t4 += t3 - t2;
        // temp_g2 = temp * g2;

    }
    cout << "EC Expo Test: " << chrono::duration<double, micro>(t4).count() / delta_vec.size() << " us" << endl;


    const string MODULUS = "25195908475657893494027183240048398571429282126204032027777137836043662020707595556264018525880784406918290641249515082189298559149176184502808489120072844992687392807287776735971418347270261896375014971824691165077613379859095700097330459748808428401797429100642458691817195118746121515172654632282216869987549182422433637259085141865462043576798423387184774447920739934236584823824281198163815010674810451660377306056201619676256133844143603833904414952634432190114657544454178424020924616515723350778707749817125772467962926386356373289912154831438167899885040445364023527381951378636564391212010397122822120720357";
    const string BASE = "2988348162058574136915891421498819466320163312926952423791023078876139";
    mpz_class g, delta_mpz, n, r;

    g.set_str(BASE.c_str(), 10);
    n.set_str(MODULUS.c_str(), 10);
    t4 = t2 - t2;
    for (size_t i = 0; i < delta_vec.size(); i++){
        mpz_set_ui(delta_mpz.get_mpz_t(), delta_vec[i]);
        t2 = chrono::high_resolution_clock::now();
        mpz_powm(r.get_mpz_t(), g.get_mpz_t(), delta_mpz.get_mpz_t(), n.get_mpz_t());
        mpz_mul(r.get_mpz_t(), r.get_mpz_t(), n.get_mpz_t());
        mpz_mod(r.get_mpz_t(), r.get_mpz_t(), n.get_mpz_t());
        t3 = chrono::high_resolution_clock::now();
        t4 += t3 - t2;
    }
    cout << "GMP exp: " << chrono::duration<double, micro>(t4).count() / delta_vec.size() << " us " << endl;

    cout << "============== EXPO COMPARE 2 ==================" << endl;

}



template<typename ppT>
bool expo_test(const size_t &limit = 10000, const int &txn_count = 1000) {

    auto t1 = chrono::high_resolution_clock::now(), t2 = chrono::high_resolution_clock::now();
    auto t3 = t2 - t2;
    srand(time(NULL));

    vector<unsigned long> exponents;
    for (size_t i = 0; i < txn_count; i++){
        unsigned long temp = rand() % limit;
        exponents.push_back(temp);
    }

    expo_test2<ppT>(exponents);
    return true;
}
int main(int argc, char **argv){

    libff::inhibit_profiling_info = true;
    libff::inhibit_profiling_counters = true;

    int L = atoi(argv[1]);
    size_t limit = (int)pow(10, L);
    cout << "=============================================================" << endl;
    cout << "Bounds of TX value: " << limit << endl;
    cout << "=============================================================" << endl;

    #ifdef CURVE_MNT6
        cout << "MNT6" << endl;
        mnt6_pp::init_public_params();
        expo_test<mnt6_pp>(limit);
        cout << "=============================================================" << endl;

    #endif

    #ifdef CURVE_MNT4
        // cout << "MNT4" << endl;
        // mnt4_pp::init_public_params();
        // expo_test<mnt4_pp>(limit);
        // cout << "=============================================================" << endl;
    #endif
    return 0;
}
