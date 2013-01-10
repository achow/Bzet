#include "Bzet.h"

#include <math.h>
#include <stdlib.h>
#include <time.h>
#include <iostream>
#include <vector>
#include <bitset>
#include <assert.h>

using namespace std;

// Define NBITS if a specific density is to be tested.

// Maximum number of Bzet-encoded bits
#ifndef SIZE
#define SIZE 10000
#endif

#define MAX_DENSITY 50

// Check if every bit set matches in std::bitset vs Bzet
#define COMPLETE

#ifndef NTESTS
#define NTESTS 50
#endif

// Which tests to run
#define TESTAND
#define TESTOR
#define TESTXOR
#define TESTNOT // Very slow!

int rawsize = (int) floor(SIZE / 8.0);

// Global timers
clock_t gstart, gend;

inline void glob_start() {
    gstart = clock();
}

inline void glob_end() {
    gend = clock();
    printf("Time taken: %f\n", (float) (gend - gstart) / CLOCKS_PER_SEC);
}

// Get a random integer. This is a separate function to allow for easy
// modification of number generation.
unsigned int rands() {
    unsigned int s;
    s = rand();
    return s;
}

// Generate (into result) nbits random integers in range [0, space]
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

void show_raw(BZET& b) {
    uint8_t *raw = new uint8_t[b.size()];
    b.hex(raw);
    for (int i = 0; i < b.size(); i++) {
        printf("0x%x ", raw[i]);
    }
    printf("\n");
    delete raw;
}

int main() {
    //cout << "START\n";
    //freopen ("correctnesstest.txt", "w", stdout);
    int density, r;
    int64_t *bits;

    for (int i = 0; i < NTESTS; i++) {
        // Seed rand()
        //srand((unsigned)time(0)); 

        // Randomly generate a density
        density = (rands() % MAX_DENSITY) + 1;

        // Start timer for this run.
        clock_t start = clock();

        // BEGIN Generating bitsets.

        // Get number of bits to set
#ifndef NBITS
        int nbits = (int) floor(SIZE * (density / 100.0));
#else
        int nbits = NBITS;
#endif

        cout << "Doing round " << (i + 1) << " of " << NTESTS << endl;
        cout << "Density is " << density << ", " << nbits << " to be set" << endl;

        cout << "Building first bitset" << endl;

        glob_start();

        // Generate first set of bit indices
        vector<unsigned int> bits1;
        gen(bits1, nbits, SIZE);
        
        bitset<SIZE> bitset1;
        BZET bzet1;

        // Bzet must be empty.
        // Note: Size of bzet is 1 byte when empty.
        assert(bzet1.empty() && bzet1.size() == 1);

        // Set bits in Bzet and std::bitset
        for (int j = 0; j < nbits; j++) {
            /*BZET::align(bzet1, BZET(bits1[j]));
            cout << j << " " << bits1[j] << endl;
            BZET(bits1[j]).printBzet();
            bzet1.printBzet();*/
            bitset1.set(bits1[j]);
            bzet1.set(bits1[j]);
#ifdef COMPLETE
            // Verify that we actually set the bit in the Bzet
            bool test = (bzet1.count() == j + 1);
            if (!test) {
                //bzet1.dump();
                bzet1.printBzet();
                cout << "count(): " << bzet1.count() << endl;
                cout << "expected: " << (j + 1) << endl;
                assert(bzet1.count() == j + 1);
            }

            assert(bzet1.at(bits1[j]));
#endif
        }

        glob_end();

        cout << "Done setting" << endl;

        cout << "First bzet is size " << bzet1.size() << ", ratio is " << (1.0 * bzet1.size() / rawsize * 100.0) << endl;

        cout << "Building second bitset" << endl;

        glob_start();

        // Generate second set of bit indices
        vector<unsigned int> bits2;
        gen(bits2, nbits, SIZE);

        bitset<SIZE> bitset2;
        BZET bzet2;

        for (int j = 0; j < nbits; j++) {
            int bit = bits2[j];

            bitset2.set(bit);
            bzet2.set(bit);
#ifdef COMPLETE
            // Verify that we actually set the bit in the Bzet
            bool test = (bzet2.count() == j + 1);
            if (!test) {
                bzet2.printBzet();
                cout << "count(): " << bzet2.count() << endl;
                cout << "expected: " << (j + 1) << endl;
                assert(bzet2.count() == j + 1);
            }

            assert(bzet2.at(bits2[j]));
#endif
        }

        glob_end();

        cout << "Second bzet is size " << bzet2.size() << ", ratio is " << (1.0 * bzet2.size() / rawsize * 100) << endl;

        // END Generating bitsets

#ifdef TESTAND
        cout << "Testing AND" << endl;

        glob_start();
        BZET bzetAND = bzet1 & bzet2;
        bitset<SIZE> bitsetAND = bitset1 & bitset2;

        if ((size_t) bzetAND.count() != (size_t) bitsetAND.count()) {
            cout << "AND count mismatch " << bzetAND.count() << " " << bitsetAND.count() << endl;
            exit(1);
        }

        cout << "Bit counts match for AND: " << bzetAND.count() << endl;

        bits = new int64_t[bzetAND.count()];
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
#endif // TESTAND

#ifdef TESTOR
        cout << "Testing OR" << endl;
        glob_start();

        BZET bzetOR = bzet1 | bzet2;
        bitset<SIZE> bitsetOR = bitset1 | bitset2;

        if ((size_t) bzetOR.count() != (size_t) bitsetOR.count()) {
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
#endif // TESTOR

#ifdef TESTXOR
        cout << "Testing XOR" << endl;

        glob_start();

        BZET bzetXOR = bzet1 ^ bzet2;
        bitset<SIZE> bitsetXOR = bitset1 ^ bitset2;

        if ((size_t) bzetOR.count() != (size_t) bitsetOR.count()) {
            cout << "XOR count mismatch" << endl;
            exit(1);
        }

#ifdef COMPLETE

        bits = new int64_t[bzetXOR.count()];
        bzetXOR.getBits(bits);

        for (int i = 0; i < bzetXOR.count(); i++) {
            int bit = bits[i];
            if (!bitsetXOR.test(bit)) {
                bzet1.printBzet();
                bzet2.printBzet();
                bzetXOR.printBzet();
                cout << "bit " << bit << " not set in bitset" << endl;
                exit(1);
            }
        }

        delete bits;
#endif

        glob_end();
#endif // TESTXOR

#ifdef TESTNOT
        cout << "Testing NOT" << endl;

        glob_start();

        BZET mask;
        mask.setRange(0, SIZE);
        BZET::align(mask, bzet1);
        BZET::align(mask, bzet2);
        BZET bzetNOT1 = ~bzet1 & mask;
        BZET bzetNOT2 = ~bzet2 & mask;
        bitset<SIZE> bitsetNOT1 = ~bitset1;
        bitset<SIZE> bitsetNOT2 = ~bitset2;

        cout << "bzet1 NOT " << bzetNOT1.size() << ", ratio is " << (1.0 * bzetNOT1.size() / rawsize * 100) << endl;
        cout << "bzet2 NOT " << bzetNOT2.size() << ", ratio is " << (1.0 * bzetNOT2.size() / rawsize * 100) << endl;

        cout << "Verifying NOT operation" << endl;

#ifdef COMPLETE
        for (int i = 0; i < SIZE; i++) {
            if (bitsetNOT1.test(i) != bzetNOT1.at(i)) {
                cout << "bit " << i << " mismatch: should be " << bitsetNOT1.test(i) << endl;
                exit(1);
            }
            if (bitsetNOT2.test(i) != bzetNOT2.at(i)) {
                cout << "bit " << i << " mismatch: should be " << bitsetNOT2.test(i) << endl;
                exit(1);
            }
        }
#else
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

        cout << "Verifying NOT counts" << endl;

        if ((size_t) bzetNOT1.count() != (size_t) bitsetNOT1.count()) {
            cout << "NOT count mismatch1" << bzetNOT1.count() << " - >" << bitsetNOT1.count() << endl;
            cout << "LAST is " << bzet1.lastBit() << endl;
            cout << "NOT LAST is " << bzetNOT1.lastBit() << endl;
            exit(1);
        }
        if ((size_t) bzetNOT2.count() != (size_t) bitsetNOT2.count()) {
            cout << "NOT count mismatch2" << endl;
            exit(1);
        }

        glob_end();
#endif // TESTNOT

        clock_t end = clock();

        cout << "Time taken: " << (1.0 * (end - start) / CLOCKS_PER_SEC) << endl;
        cout << endl;
    }
}
