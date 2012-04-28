#define _CRT_RAND_S

#define NODE_ELS 8
#include "Bzet.h"

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <bitset>
#include <assert.h>

using namespace std;

#define NBITS 10000
#define SIZE 80000
#define MAX_DENSITY 20
#define COMPLETE

int ntests = 5;
int rawsize = (int) floor(SIZE / 8.0);

clock_t gstart, gend;

inline void glob_start() {
    gstart = clock();
}

inline void glob_end() {
    gend = clock();
    printf("Time taken: %f\n", (float) (gend - gstart) / CLOCKS_PER_SEC);
}

unsigned int rands() {
    unsigned int s;
    //rand_s(&s);
    s = rand();
    return s;
}

void gen(vector<unsigned int>& result, unsigned int nbits, unsigned int space) {
    vector<unsigned int> nums;
    for (unsigned int i = 0; i < space; i++)
        nums.push_back(i);
       
    cout << "generating" << endl;
    for (unsigned int j = 0; j < nbits; j++) {
        //cout << j << endl;
        unsigned int n = rands() % nums.size();

        unsigned int bit = nums[n];
        nums.erase(nums.begin() + n);

        result.push_back(bit);
    }
}

int main() {
    //cout << "START\n";
    //freopen ("correctnesstest.txt", "w", stdout);
    int density;

    for (int i = 0; i < ntests; i++) {
        srand((unsigned)time(0)); 
        density = (rands() % MAX_DENSITY) + 1;

        clock_t start = clock();

#ifndef NBITS
        int nbits = (int) floor(SIZE * (density / 100.0));
#else
        int nbits = NBITS;
#endif

        cout << "Doing round " << (i + 1) << " of " << ntests << endl;
        cout << "Density is " << density << ", " << nbits << " to be set" << endl;

        vector<unsigned int> bits1;
        gen(bits1, nbits, SIZE);
        
        cout << "Done generating bit indices" << endl;

        bitset<SIZE> bitset1;
        BZET bzet1;

#ifdef COMPLETE
        for (int x = 0; x < SIZE; x++)
            assert(bzet1.at(x) == 0);
#endif

        cout << "Building first bitset" << endl;

        glob_start();

        int nset = 0;
        for (int j = 0; j < nbits; j++) {
            nset++;
            /*BZET::align(bzet1, BZET(bits1[j]));
            cout << j << " " << bits1[j] << endl;
            BZET(bits1[j]).printBzet();
            bzet1.printBzet();*/
            bitset1.set(bits1[j]);
            bzet1.set(bits1[j]);
#ifdef COMPLETE
            bool test = (bzet1.count() == nset);
            if (!test) {
                //bzet1.dump();
                bzet1.printBzet();
                cout << "count(): " << bzet1.count() << endl;
                cout << "expected: " << nset << endl;
                assert(bzet1.count() == nset);
            }

            assert(bzet1.at(bits1[j]));
#endif
        }

        glob_end();

        cout << "Done setting" << endl;

        cout << "First bzet is size " << bzet1.size() << ", ratio is " << (1.0 * bzet1.size() / rawsize * 100.0) << endl;

        cout << "Building second bitset" << endl;

        vector<unsigned int> bits2;
        gen(bits2, nbits, SIZE);

        glob_start();

        bitset<SIZE> bitset2;
        BZET bzet2;

        for (int j = 0; j < nbits; j++) {
            int bit = bits2[j];

            bitset2.set(bit);
            bzet2.set(bit);
        }

#ifdef COMPLETE
        assert(bzet2.count() == nbits);
#endif

        glob_end();

        cout << "Second bzet is size " << bzet2.size() << ", ratio is " << (1.0 * bzet2.size() / rawsize * 100) << endl;

        //and
        cout << "Testing AND" << endl;

        glob_start();
        BZET bzetAND = bzet1 & bzet2;
        bitset<SIZE> bitsetAND = bitset1 & bitset2;

        if (bzetAND.count() != bitsetAND.count()) {
            cout << "AND count mismatch " << bzetAND.count() << " " << bitsetAND.count() << endl;
            exit(1);
        }

        cout << "Bit counts match for AND: " << bzetAND.count() << endl;

        int64_t* bits = new int64_t[bzetAND.count()];
        bzetAND.getBits(bits);

        for (int i = 0; i < bzetAND.count(); i++) {
            int bit = bits[i];
            if (!bitsetAND.test(bit)) {
                cout << "bit " << bit << " not set in bitset" << endl;
                exit(1);
            }
        }

        delete bits;

        glob_end();

        cout << "There are " << bzetAND.count() << " bits in common" << endl;

        glob_start();

        //or
        cout << "Testing OR" << endl;
        BZET bzetOR = bzet1 | bzet2;
        bitset<SIZE> bitsetOR = bitset1 | bitset2;

        if (bzetOR.count() != bitsetOR.count()) {
            cout << "OR count mismatch" << endl;
            exit(1);
        }

#ifdef COMPLETE
        assert(bzetOR.count() == nbits * 2 - bzetAND.count());

        bits = new int64_t[bzetOR.count()];
        bzetOR.getBits(bits);

        for (int i = 0; i < bzetOR.count(); i++) {
            int bit = bits[i];
            if (!bitsetOR.test(bit)) {
                cout << "bit " << bit << " not set in bitset" << endl;
                exit(1);
            }
        }

        delete bits;
#endif
        
        glob_end();

        glob_start();

        //xor
        cout << "Testing XOR" << endl;
        BZET bzetXOR = bzet1 ^ bzet2;
        bitset<SIZE> bitsetXOR = bitset1 ^ bitset2;

        if (bzetOR.count() != bitsetOR.count()) {
            cout << "XOR count mismatch" << endl;
            exit(1);
        }

#ifdef COMPLETE

        bits = new int64_t[bzetXOR.count()];
        bzetXOR.getBits(bits);

        for (int i = 0; i < bzetXOR.count(); i++) {
            int bit = bits[i];
            if (!bitsetXOR.test(bit)) {
                cout << "bit " << bit << " not set in bitset" << endl;
                exit(1);
            }
        }

        delete bits;
#endif

        glob_end();

        glob_start();

        /*//not
        cout << "Testing NOT" << endl;
        BZET mask(0, SIZE);
        BZET::align(mask, bzet1);
        BZET::align(mask, bzet2);
        BZET bzetNOT1 = ~bzet1 & mask;
        BZET bzetNOT2 = ~bzet2 & mask;
        bitset<SIZE> bitsetNOT1 = ~bitset1;
        bitset<SIZE> bitsetNOT2 = ~bitset2;

        if (bzetNOT1.count() != bitsetNOT1.count()) {
            cout << "NOT count mismatch1" << bzetNOT1.count() << "-" << bitsetNOT1.count() << endl;
            exit(1);
        }
        if (bzetNOT2.count() != bitsetNOT2.count()) {
            cout << "NOT count mismatch2" << endl;
            exit(1);
        }*/
        /*
        cout << "bzet1 NOT " << bzetNOT1.size() << ", ratio is " << (1.0 * bzetNOT1.size() / rawsize * 100) << endl;
        cout << "bzet2 NOT " << bzetNOT2.size() << ", ratio is " << (1.0 * bzetNOT2.size() / rawsize * 100) << endl;

        bits = new int64_t[bzetNOT1.count()];
        bzetNOT1.getBits(bits);

        for (int i = 0; i < bzetNOT1.count(); i++) {
            int bit = bits[i];
            if (!bitsetNOT1.test(bit)) {
                cout << "bit " << bit << " not set in bitset" << endl;
                exit(1);
            }
        }

        delete bits;

        bits = new int64_t[bzetNOT2.count()];
        bzetNOT2.getBits(bits);

        for (int i = 0; i < bzetNOT2.count(); i++) {
            int bit = bits[i];
            if (!bitsetNOT2.test(bit)) {
                cout << "bit " << bit << " not set in bitset" << endl;
                exit(1);
            }
        }

        delete bits;
#endif
        */
        glob_end();

        clock_t end = clock();

        cout << "Time taken: " << (1.0 * (end - start) / CLOCKS_PER_SEC) << endl;
        cout << endl;
    }
}