#define _CRT_RAND_S

#define NODE_ELS 8
#include "Bzet.h"

#include <stdlib.h>
#include <time.h>
#include <bitset>
#include <stdio.h>
#include <iostream>
#include <assert.h>
#include <vector>

using namespace std;

#define SIZE 100000
#define NBITS 10000
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
    rand_s(&s);
    return s;
}

void gen(vector<unsigned int>& result, unsigned int nbits, unsigned int space) {
    vector<unsigned int> nums;
    for (unsigned int i = 0; i < space; i++)
        nums.push_back(i);
       
    for (unsigned int j = 0; j < nbits; j++) {
        unsigned int n = rands() % nums.size();

        unsigned int bit = nums[n];
        nums.erase(nums.begin() + n);

        result.push_back(bit);
    }
}

void dumpstep(BZET_PTR b) {
    int n = b->nhalfnodes;

    cout << "bufn = " << b->nbufhalfnodes << endl;

    for (int i = 0; i < n; i++) {
        if ((i % 10) == 0)
            printf("\n");
        printf("0x%.2X ", b->step[i]);
    }
    printf("\n");
}

int main() {
    /*Bzet8 *b1 = Bzet8_new();
    Bzet8 *b2 = Bzet8_new();

    Bzet8_SET(b1, 1);
    Bzet8_SET(b1, 1000);
    Bzet8_SET(b2, 1);
    Bzet8_SET(b2, 1000);

    cout << Bzet8_TEST(b1, 1000) << endl;

    Bzet8 *result = Bzet8_AND(b1, b2);
    Bzet8_HEX(result);

    return 1;*/

    freopen ("correctnesstest.txt", "w", stdout);
    int density;

    for (int i = 0; i < ntests; i++) {
        //srand((unsigned)time(0)); 
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
        
        cout << "Done generating bit indices." << endl;

        bitset<SIZE> bitset1;
        BZET_PTR bzet1 = BZET_FUNC(new)();

        cout << "Building first bitset." << endl;

        glob_start();

        int nset = 0;
        for (int j = 0; j < nbits; j++) {
            nset++;

            cout << "bzet1 nbuf = " << bzet1->nbufhalfnodes << endl;

            BZET_PTR old = BZET_FUNC(clone)(bzet1);

            bitset1.set(bits1[j]);
            BZET_FUNC(SET)(bzet1, bits1[j]);

#ifdef COMPLETE
            cout << j << ": " << bits1[j] << endl;
            assert(BZET_FUNC(TEST)(bzet1, bits1[j]) == 1);

            bool test = (BZET_FUNC(COUNT)(bzet1) == nset);
            /*BZET_PTR x = BZET_FUNC(new)(bits1[j]);
            align(x, bzet1);
            BZET_FUNC(HEX)(x);
            BZET_FUNC(destroy)(x);
            BZET_FUNC(HEX)(bzet1);*/
            if (!test) {
                BZET_FUNC(HEX)(old);
                dumpstep(bzet1);
                BZET_FUNC(HEX)(bzet1);
                //dumpstep(bzet1);
                cout << "count(): " << BZET_FUNC(COUNT)(bzet1) << endl;
                cout << "expected: " << nset << endl;
                exit(1);
            }
            BZET_FUNC(destroy)(old);
#endif
        }

        glob_end();

        return 0;

#ifdef COMPLETE
        cout << "Done building first bitset" << endl;
#endif

        cout << "First bzet is size " << BZET_FUNC(size)(bzet1) << ", ratio is " << (1.0 * BZET_FUNC(size)(bzet1) / rawsize * 100.0) << endl;

        cout << "Building second bitset" << endl;

        bits1.clear();
        gen(bits1, nbits, SIZE);
        
        cout << "Done generating bit indices." << endl;

        bitset<SIZE> bitset2;
        BZET_PTR bzet2 = BZET_FUNC(new)();

        cout << "Building second bitset." << endl;

        glob_start();

        nset = 0;
        for (int j = 0; j < nbits; j++) {
            nset++;

            bitset2.set(bits1[j]);
            BZET_FUNC(SET)(bzet2, bits1[j]);

#ifdef COMPLETE
            assert(BZET_FUNC(TEST)(bzet2, bits1[j]) == 1);

            bool test = (BZET_FUNC(COUNT)(bzet2) == nset);
            if (!test) {
                cout << "count(): " << BZET_FUNC(COUNT)(bzet2) << endl;
                cout << "expected: " << nset << endl;
                exit(1);
            }
#endif
        }

        glob_end();

        cout << "Second bzet is size " << BZET_FUNC(size)(bzet2) << ", ratio is " << (1.0 * BZET_FUNC(size)(bzet2) / rawsize * 100) << endl;

        //and
        cout << "Testing AND" << endl;

        glob_start();
        BZET_PTR bzetAND = BZET_FUNC(AND)(bzet1, bzet2);
        bitset<SIZE> bitsetAND = bitset1 & bitset2;

        if (BZET_FUNC(COUNT)(bzetAND) != bitsetAND.count()) {
            cout << "AND count mismatch " << BZET_FUNC(COUNT)(bzetAND) << " " << bitsetAND.count() << endl;
            exit(1);
        }

#ifdef COMPLETE

        int64_t* bits = new int64_t[BZET_FUNC(COUNT)(bzetAND)];
        BZET_FUNC(getBits)(bzetAND, bits);

        for (int i = 0; i < BZET_FUNC(COUNT)(bzetAND); i++) {
            int bit = bits[i];
            if (!bitsetAND.test(bit)) {
                cout << "bit " << bit << " not set in bitset" << endl;
                exit(1);
            }
        }

        delete bits;
#endif

        glob_end();

        cout << "There are " << BZET_FUNC(COUNT)(bzetAND) << " bits in common" << endl;

       /* glob_start();

        //or
        cout << "Testing OR" << endl;
        Bzet4 bzetOR = bzet1 | bzet2;
        bitset<SIZE> bitsetOR = bitset1 | bitset2;

        if (bzetOR.count() != bitsetOR.count()) {
            cout << "OR count mismatch" << endl;
            exit(1);
        }

#ifdef COMPLETE
        cout << "Bitcount: ";
        count_bits(bzetOR);

        assert(bzetOR.count() == nbits * 2 - bzetAND.count());

        bits = new long long[bzetOR.count()];
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
        Bzet4 bzetXOR = bzet1 ^ bzet2;
        bitset<SIZE> bitsetXOR = bitset1 ^ bitset2;

        if (bzetOR.count() != bitsetOR.count()) {
            cout << "XOR count mismatch" << endl;
            exit(1);
        }

#ifdef COMPLETE
        cout << "Bitcount: ";
        count_bits(bzetXOR);

        bits = new long long[bzetXOR.count()];
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

        //not
        cout << "Testing NOT" << endl;
        Bzet4 bzetNOT1 = ~bzet1;
        Bzet4 bzetNOT2 = ~bzet2;
        bitset<SIZE> bitsetNOT1 = ~bitset1;
        bitset<SIZE> bitsetNOT2 = ~bitset2;

#ifdef COMPLETE
        if (bzetNOT1.count() != bitsetNOT1.count()) {
            cout << "NOT count mismatch" << endl;
            exit(1);
        }
        if (bzetNOT2.count() != bitsetNOT2.count()) {
            cout << "NOT count mismatch" << endl;
            exit(1);
        }
        //cout << "Bitcount: ";
        //count_bits(bzetOR);

        cout << "bzet1 NOT " << bzetNOT1.size() << ", ratio is " << (1.0 * bzetNOT1.size() / rawsize * 100) << endl;
        cout << "bzet2 NOT " << bzetNOT2.size() << ", ratio is " << (1.0 * bzetNOT2.size() / rawsize * 100) << endl;

        bits = new long long[bzetNOT1.count()];
        bzetNOT1.getBits(bits);

        for (int i = 0; i < bzetNOT1.count(); i++) {
            int bit = bits[i];
            if (!bitsetNOT1.test(bit)) {
                cout << "bit " << bit << " not set in bitset" << endl;
                exit(1);
            }
        }

        delete bits;

        bits = new long long[bzetNOT2.count()];
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
        glob_end();
        */
        clock_t end = clock();

        cout << "Time taken: " << (1.0 * (end - start) / CLOCKS_PER_SEC) << endl;
        cout << endl;
    }
}