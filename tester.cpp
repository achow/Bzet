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
//#define NBITS 10000
#define MAX_DENSITY 20
#define COMPLETE

int ntests = 50;
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
    //freopen ("correctnesstest.txt", "w", stdout);
    int density;

    for (int i = 0; i < ntests; i++) {
        srand(time(NULL)); 
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

            //cout << "bzet1 nbuf = " << bzet1->nbufhalfnodes << endl;

            //BZET_PTR old = BZET_FUNC(clone)(bzet1);

            bitset1.set(bits1[j]);
            BZET_FUNC(SET)(bzet1, bits1[j]);

#ifdef COMPLETE
            //cout << j << ": " << bits1[j] << endl;
            assert(BZET_FUNC(TEST)(bzet1, bits1[j]) == 1);
            //BZET_FUNC(HEX)(bzet1);
            bool test = (BZET_FUNC(COUNT)(bzet1) == nset);
            /*BZET_PTR x = BZET_FUNC(new)(bits1[j]);
            align(x, bzet1);
            BZET_FUNC(HEX)(x);
            BZET_FUNC(destroy)(x);
            BZET_FUNC(HEX)(bzet1);*/
            if (!test) {
                //BZET_FUNC(HEX)(old);
                dumpstep(bzet1);
                BZET_FUNC(HEX)(bzet1);
                //dumpstep(bzet1);
                cout << "count(): " << BZET_FUNC(COUNT)(bzet1) << endl;
                cout << "expected: " << nset << endl;
                exit(1);
            }
            //BZET_FUNC(destroy)(old);
#endif
        }

        glob_end();

#ifdef COMPLETE
        cout << "Done building first bitset" << endl;
#endif

        cout << "First bzet is size " << BZET_FUNC(size)(bzet1) << ", ratio is " << (1.0 * BZET_FUNC(size)(bzet1) / rawsize * 100.0) << endl;

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

        cout << "AND counts match" << endl;
        cout << "There are " << BZET_FUNC(COUNT)(bzetAND) << " bits in common" << endl;
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

        cout << "Bits verified, AND succeeded" << endl;
#endif

        glob_end();

        glob_start();

        //or
        cout << "Testing OR" << endl;
        BZET_PTR bzetOR = BZET_FUNC(OR)(bzet1, bzet2);
        bitset<SIZE> bitsetOR = bitset1 | bitset2;

        if (BZET_FUNC(COUNT)(bzetOR) != bitsetOR.count()) {
            cout << "OR count mismatch bzet=" << BZET_FUNC(COUNT)(bzetOR) << "bitset=" << bitsetOR.count() << endl;
            BZET_FUNC(HEX)(bzetOR);
            exit(1);
        }

        cout << "OR counts match" << endl;

#ifdef COMPLETE
        bits = new int64_t[BZET_FUNC(COUNT)(bzetOR)];
        BZET_FUNC(getBits)(bzetOR, bits);

        for (int i = 0; i < BZET_FUNC(COUNT)(bzetOR); i++) {
            int bit = bits[i];
            if (!bitsetOR.test(bit)) {
                cout << "bit " << bit << " not set in bitset" << endl;
                exit(1);
            }
        }

        delete bits;
#endif
        
        glob_end();

        cout << "OR bits verified, OR succeeded" << endl;

        glob_start();

        //xor
        cout << "Testing XOR" << endl;
        BZET_PTR bzetXOR = BZET_FUNC(XOR)(bzet1, bzet2);
        bitset<SIZE> bitsetXOR = bitset1 ^ bitset2;

        if (BZET_FUNC(COUNT)(bzetXOR) != bitsetXOR.count()) {
            cout << "XOR count mismatch" << endl;
            exit(1);
        }

#ifdef COMPLETE
        bits = new int64_t[BZET_FUNC(COUNT)(bzetXOR)];
        BZET_FUNC(getBits)(bzetXOR, bits);

        for (int i = 0; i < BZET_FUNC(COUNT)(bzetXOR); i++) {
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
        /*cout << "Testing NOT" << endl;
        BZET_PTR bzetNOT1 = BZET_FUNC(NOT)(bzet1);
        BZET_PTR bzetNOT2 = BZET_FUNC(NOT)(bzet2);
        bitset<SIZE> bitsetNOT1 = ~bitset1;
        bitset<SIZE> bitsetNOT2 = ~bitset2;

#ifdef COMPLETE
        if (BZET_FUNC(COUNT)(bzetNOT1) != bitsetNOT1.count()) {
            cout << "NOT count mismatch" << endl;
            exit(1);
        }
        if (BZET_FUNC(COUNT)(bzetNOT2) != bitsetNOT2.count()) {
            cout << "NOT count mismatch" << endl;
            exit(1);
        }*/

        /*
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

        delete bits;*/

        glob_end();
        
        clock_t end = clock();

        cout << "Time taken: " << (1.0 * (end - start) / CLOCKS_PER_SEC) << endl;

        BZET_FUNC(destroy)(bzet1);
        BZET_FUNC(destroy)(bzet2);
        BZET_FUNC(destroy)(bzetAND);
        BZET_FUNC(destroy)(bzetOR);
        BZET_FUNC(destroy)(bzetXOR);

        cout << endl;
    }
}