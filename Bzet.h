/*****************************************************************************
 * Generic code for Bzet8/Bzet16/Bzet32
 *
 * Things to do:
 *  - Good error messages, better error handling in general (?)
 *  - Currently C-style C++, changes needed to be C compatible
 *  - Performance can be improved by using popcount compiler intrinsics
 *    instead of manually counting bits
 *
 * Author: Alex Chow
 *****************************************************************************/

#ifndef BZET_H_ 
#define BZET_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <assert.h>

#ifndef NODE_ELS
#define NODE_ELS 8
#endif

// TODO: bool typedef?

#define PASTE_(x,y) x ## y
#define PASTE(x,y) PASTE_(x,y)

#define PASTE_UNDER_(x,y) x ## _ ## y
#define PASTE_UNDER(x,y) PASTE_UNDER_(x,y)

#define STRINGIFY_(x) #x
#define STRINGIFY(x) STRINGIFY_(x)

#define BZET PASTE(Bzet, NODE_ELS)
#define BZET_PTR BZET *
#define BZET_FUNC(x) PASTE_UNDER(BZET, x)
#define POW PASTE(pow, NODE_ELS)

#if NODE_ELS == 32
typedef uint32_t halfnode_t;
#define NPOWERS 7
static const unsigned int PASTE(powersof, NODE_ELS)[NPOWERS] = 
    { 1, 32, 1024, 32768, 1048576, 33554432, 1073741824 };
#elif NODE_ELS == 16
typedef uint16_t halfnode_t
#define NPOWERS 9
static const unsigned int PASTE(powersof, NODE_ELS)[NPOWERS] = 
    { 1, 16, 256, 4096, 65536, 1048576, 16777216, 268435456, 4294967296 };
#elif NODE_ELS == 8
typedef uint8_t halfnode_t;
#define NPOWERS 10
static const unsigned int PASTE(powersof, NODE_ELS)[NPOWERS] = 
    { 1, 8, 64, 512, 4096, 32768, 262144, 2097152, 16777216, 134217728 };
#else
#error "Invalid NODE_ELS provided"
#endif

#define INITIAL_ALLOC 1024
#define RESIZE_SCALE 2

#define EXPORT_TAGS

enum OP { 
    OP_0000, OP_0001, OP_0010, OP_0011, OP_0100, OP_0101, OP_0110, OP_0111,
    OP_1000, OP_1001, OP_1010, OP_1011, OP_1100, OP_1101, OP_1110, OP_1111,
    AND = 1, XOR = 6, OR = 7, NOR = 8, NOT = 10, NAND = 14 };
enum ACTION { DA0, DA1, DB0, DB1, CA, CB, NA, NB };
enum NODETYPE { SATURATED, EMPTY, NORMAL };

static const ACTION optable[64] = {       
        DB0, DA0, DB0, DA0 , //00 0000 FALSE   Result
        DB0, DA0, CB,  CA  , //01 0001 AND        |
        CB,  CA,  NB,  DA0 , //02 0010 A<-B       |
        DB0, CA,  DB0, CA  , //03 0011 A          V
        DB0, DA0, DB0, CA  , //04 0100 ????
        CB,  DA0, CB,  DA0 , //05 0101 B
        CB,  CA,  NB,  NA  , //06 0110 XOR
        CB,  CA,  DB1, DA1 , //07 0111 OR
        NB,  NA,  DB0, DA0 , //08 1000 NOR
        NB,  NA,  CB,  CA  , //09 1001 EQ
        NB,  DA0, NB,  DA0 , //10 1010 ~B
        CB,  CA,  CB,  CA  , //11 1011 ????
        DB0, NA,  DB0, NA  , //12 1100 ~A
        DB0, NA,  CB,  CA  , //13 1101 A->B
        DB1, DA1, NB,  NA  , //14 1110 NAND
        DB1, DA1, DB1, DA1 }; //15 1111 TRUE
     // 0T   T0   1T   T1 
     // 0    1    2    3   

// Bzet struct
typedef struct {
    size_t nbufhalfnodes;
    size_t nhalfnodes;
    halfnode_t* bzet; //points to the bzet
    unsigned char* step; //points to an array that holds stepThrough values
    unsigned char depth;
} BZET;


// FORWARD DECLARATIONS


// Constructors and Destructors

// Bzet_new()
// Creates an empty Bzet
EXPORT_TAGS BZET_PTR BZET_FUNC(new)();

// Bzet_new(bit)
// Creates a Bzet with bit bit set
EXPORT_TAGS BZET_PTR BZET_FUNC(new)(int64_t bit);

// Bzet_new(startbit, len)
// Range constructor. Creates a bzet with len bits starting from startbit set.
EXPORT_TAGS BZET_PTR BZET_FUNC(new)(int64_t startbit, int64_t len);

// Bzet_destroy(b)
// Frees all memory associated with b
EXPORT_TAGS void BZET_FUNC(destroy)(BZET_PTR b);

// Bzet_clone(b)
// Creates a clone of b
EXPORT_TAGS BZET_PTR BZET_FUNC(clone)(BZET_PTR b);

// Bzet4(void* data, int size);


// Operators

// Bzet_setequal(left, right)
// Does a deep copy of the contents of right to left
EXPORT_TAGS void BZET_FUNC(setequal)(BZET_PTR left, BZET_PTR right);

// Bzet_NOT(b)
// Returns a new bzet which is the bitwise NOT of b
EXPORT_TAGS BZET_PTR BZET_FUNC(NOT)(BZET_PTR b);

// Bzet_INVERT(b)
// In place NOT of b
EXPORT_TAGS void BZET_FUNC(INVERT)(BZET_PTR b);

// Bzet_OR(left, right)
// Bitwise OR of left and right
EXPORT_TAGS BZET_PTR BZET_FUNC(OR)(BZET_PTR left, BZET_PTR right);

// Bzet_AND(left, right)
// Bitwise AND of left and right
EXPORT_TAGS BZET_PTR BZET_FUNC(AND)(BZET_PTR left, BZET_PTR right);

// Bzet_XOR(left, right)
// Bitwise XOR of left and right
EXPORT_TAGS BZET_PTR BZET_FUNC(XOR)(BZET_PTR left, BZET_PTR right);

// Bzet_COMPARE(left, right)
// Test for equality between two bzets
EXPORT_TAGS BZET_PTR BZET_FUNC(COMPARE)(BZET_PTR left, BZET_PTR right);

// Bzet_binop(left, right, op)
// Generic binary operations
EXPORT_TAGS BZET_PTR BZET_FUNC(binop)(BZET_PTR left, BZET_PTR right, OP op);


// Bit setting, unsetting, and getting

// Bzet_TEST(b, bit)
// Checks if bit bit is set
EXPORT_TAGS bool BZET_FUNC(TEST)(BZET_PTR b, int64_t bit);

// Bzet_RANGE(b, start, len)
// Sets len bits starting from start
EXPORT_TAGS void BZET_FUNC(RANGE)(BZET_PTR b, int64_t start, int64_t len);

// Bzet_SET(bit)
// Sets bit bit
EXPORT_TAGS void BZET_FUNC(SET)(BZET_PTR b, int64_t bit);

// Bzet_UNSET(bit)
// Unsets bit bit
EXPORT_TAGS void BZET_FUNC(UNSET)(BZET_PTR b, int64_t bit);

// Bzet_FIRST(b)
// Gets the first bit set
EXPORT_TAGS int64_t BZET_FUNC(FIRST)(BZET_PTR b);

// Bzet_LAST(b)
// Gets the last bit set
EXPORT_TAGS int64_t BZET_FUNC(LAST)(BZET_PTR b);

// Bzet_COUNT(b)
// Counts the number of bits set
EXPORT_TAGS int64_t BZET_FUNC(COUNT)(BZET_PTR b);


// Other exported functions

// Bzet_LEV(b)
// Returns the depth of the bzet
EXPORT_TAGS int BZET_FUNC(LEV)(BZET_PTR b);

// Bzet_size(b)
// Returns the actual size of the bzet
EXPORT_TAGS size_t BZET_FUNC(size)(BZET_PTR b);

// Bzet_HEX(b)
// "Pretty print"s the bzet to stdout
EXPORT_TAGS void BZET_FUNC(HEX)(BZET_PTR b);

// Bzet_repr(b, target)
// Copies the canonical bzet representation to target
EXPORT_TAGS void BZET_FUNC(repr)(BZET_PTR b, void *target);

// Bzet_CLEAN(b)
// Empties the bzet
EXPORT_TAGS void BZET_FUNC(CLEAN)(BZET_PTR b);

// Bzet_EMPTY(b)
// Checks if the bzet is empty
EXPORT_TAGS bool BZET_FUNC(EMPTY)(BZET_PTR b);

// Bzet_getBits(b, bits, limit, start)
// Retrieve the indices of bits set
EXPORT_TAGS int64_t BZET_FUNC(getBits)(BZET_PTR b, int64_t* bits, int64_t limit = 0, int64_t start = 0);


// Auxiliary functions

NODETYPE _binop(BZET_PTR left, BZET_PTR right, OP op, int lev, size_t left_loc = 1, size_t right_loc = 1, size_t loc = 1);

// Recursively print the bzet in "pretty print"
void _printBzet(BZET_PTR b, int stdOffset, FILE* target, int depth, size_t loc = 0, int offset = 0, bool pad = 0);

// Align two Bzets to the same level
void align(BZET_PTR b1, BZET_PTR b2);

// Strip leading unnecessary nodes
void normalize(BZET_PTR b);

// Return the location of the next node after traversing the subtree starting at loc
size_t stepThrough(BZET_PTR b, size_t loc);

// In-place bitwise NOT of the subtree whose root is at loc
void subtreeNot(BZET_PTR b, size_t loc, int depth);


// Inline auxiliary functions

// Error message printing and optional exiting
inline
void display_error(char* message, bool fatal = false, FILE* output = stderr) {
    fprintf(output, "%s\n", message);
#if (defined _DEBUG || defined DEBUG)
    dump();
#endif
    if (fatal)
        exit(1);
}

// Common Bzet constructor initialization
inline 
BZET_PTR init() {
    // Allocate bzet struct
    BZET_PTR b = (BZET_PTR) malloc(sizeof(BZET));
    if (!b)
        return NULL;

    // Allocate bzet node array
    b->bzet = (halfnode_t*) malloc(INITIAL_ALLOC * sizeof(halfnode_t));
    if (!b->bzet) {
        free(b);
        return NULL;
    }

    // Allocate step array
    b->step = (unsigned char*) malloc(INITIAL_ALLOC * sizeof(unsigned char));
    if (!b->step) {
        free(b);
        free(b->bzet);
        return NULL;
    }

    b->depth = 0;
    b->nhalfnodes = 0;
    b->nbufhalfnodes = INITIAL_ALLOC;

    return b;
}

// Resizes the buffers in the Bzet if necessary
inline
void resize(BZET_PTR b, size_t nhalfnodes) {
    // If reallocation is required
    if (nhalfnodes > b->nbufhalfnodes) {
        // Get new size required by growing it by RESIZE_SCALE repeatedly
        while (nhalfnodes > b->nbufhalfnodes)
            b->nbufhalfnodes *= RESIZE_SCALE;

        // Reallocate bzet
        halfnode_t *bzet_temp = (halfnode_t *) realloc(b->bzet, b->nbufhalfnodes * sizeof(halfnode_t));
        /*if (!bzet_temp) {
            // TODO: Add some error message
            return;
        }*/

        // Reallocate step
        unsigned char *step_temp = (unsigned char*) realloc(b->step, b->nbufhalfnodes * sizeof(unsigned char));
        /*if (!step_temp) {
            // TODO: Add some error message
            return;
        }*/

        // Check that realloc succeeded
        // TODO: Replace with better checking above
        if (!b->bzet || !b->step) {
            fprintf(stderr, "Fatal error: Resizing bzet failed attempting to allocate %l bytes\n", b->nbufhalfnodes * sizeof(halfnode_t));
            display_error("", true);
        }

        b->bzet = bzet_temp;
        b->step = step_temp;
    }

    // Update internal halfnode counter
    b->nhalfnodes = nhalfnodes;
}

// Trusty fast pow8/pow16/pow32 (assumes n >= 0)
inline
size_t POW(int n) {
    // Lookup power in table if possible
    if (n < NPOWERS)
        return PASTE(powersof, NODE_ELS)[n];

    return (size_t) pow((double) NODE_ELS, n);
}

// Does operations for two data bits, used in Bzet_binop
inline
int do_data_op(OP op, int left_data_bit, int right_data_bit) {
    // Use op directly to build result bit
    return (op >> (3 - (((int) left_data_bit << 1) | (int) right_data_bit))) & 0x1;
}

/*
#if (defined _DEBUG || defined DEBUG)
        void printBytes(FILE* target = stdout) const;
        void validateBzet(size_t loc = 1, int lev = 0);
        void dump() const;
#endif

        static int buildDepth(int64_t size);
        static unsigned char dust(unsigned char x);
        int depthAt(size_t loc) const;

        void appendSubtree(const Bzet4& src, size_t loc);
        void dropNodes(size_t loc, int n);
*/

#if 0
#if (defined _DEBUG || defined DEBUG)
inline
void Bzet4::printBytes(FILE* target) const {
    for (int i = 0; i < b->size; ++i) {
        if (i % 10 == 0)
            printf("\n");
        fprintf(target, " 0x%.2X", b->bzet[i]); 
    } 
    fprintf(target, "\n");
}
#endif

#if (defined _DEBUG || defined DEBUG)
void Bzet4::validateBzet(size_t loc, int lev) {
    if (!lev)
        lev = b->bzet[0];

    if (lev == 1)
        return;

    size_t nextLoc = loc;
    for (int i = 3; i >= 0; --i) {
        if ((b->bzet[loc] >> i) & 1) {
            bool check = nextLoc + 1 >= 0 && nextLoc < b->size;
            if (!check) {
                printf("assert failed at nextLoc = %d, size = %d\n", (int) nextLoc, (int) b->size);
                printf("bytes are: ");
                printBytes();
            }
            assert(nextLoc + 1 >= 0 && nextLoc < b->size);
            validateBzet(nextLoc + 1, lev - 1);
            nextLoc = stepThrough(nextLoc + 1) - 1;
        }
    }

    bool check = (dust(b->bzet[loc]) == b->bzet[loc]);
    if (!check)
        printf("assert failed at loc = %d, lev = %d\n", (int) loc, lev);
    assert(check);
    /*
    for (int i = 1; i < b->size; ++i) {
        bool test = (b->step[i] > 0 && b->step[i] <= b->size);
        if (!test) {
            printf("assertion failed, i = %d, b->step[i] = %d\n", i, b->step[i]);
            dump();
        }
        assert(test);
    }*/
}

void Bzet4::dump() const {
    printf("DUMP\n");
    printBzet();
//#if (defined _DEBUG || defined DEBUG)
    printf("BYTES: ");
    printBytes();
//#endif
    printf("\nSize is %d\n", (int) b->size);
    printf("STEP: ");
    int x = 0;
    for (int i = 1; i < b->size; ++i) {
        if (x % 10 == 0)
            printf("\n");
        printf("%5d", (int) b->step[i]);
        if (b->step[i] == 2 && b->step[i + 1] != 1) {
            printf("\n | invalid step at %d | \n", i);
        }
        x++;
    }
    printf("\n");
}
#endif
#endif 

/*****************************************************************************
 * Function Implementations
 *****************************************************************************/

// Bzet_new()
BZET_PTR BZET_FUNC(new)() { 
    BZET_PTR b = init();
    if (!b)
        return NULL;

    BZET_FUNC(CLEAN)(b);

    return b;
}

// Bzet_new(bit)
BZET_PTR BZET_FUNC(new)(int64_t bit) {
    if (bit < 0)
        return NULL;

    BZET_PTR b = init();
    if (!b)
        return NULL;

    // Build depth
    int depth = 0;
    while (POW(depth + 1) < (size_t) bit)
        depth++;

    // Resize to accomodate full bzet
    resize(b, 2*depth + 1);

    //set depth byte
    b->depth = (unsigned char) depth;

    // Set level 0 data node
    b->bzet[2*depth] = (halfnode_t) 0x1 << (NODE_ELS - 1 - (bit % NODE_ELS));

    // Zero out other nodes, since we will be only filling the tree portions next
    memset(b->bzet, 0x00, 2*depth * sizeof(halfnode_t));

    // Set tree bits, working linearly from the node at the head of the bzet
    for (int i = depth; i >= 1; i--) {
        // Get weight of bits at this level
        size_t cpow = POW(i);

        // Calculate tree bit to set
        int setbit = (int) (bit / cpow);

        // Set the tree bit
        b->bzet[2*(depth - i) + 1] = (halfnode_t) 0x1 << (NODE_ELS - 1 - setbit);
    }

    // Set b->step
    // It is only necessary to set steps corresponding to data halfnodes, since steps
    // corresponding to tree halfnodes are only there to simplify things and are 
    // never used
    // No need to check b->size <= 255 to make sure b->step values doesn't overflow
    // since it would never happen (4^255 ~ 10^153)
    for (int i = 0; i < b->nhalfnodes; i += 2)
        b->step[i] = (unsigned char) b->nhalfnodes - i;

#if (defined _DEBUG || defined DEBUG)
    validateBzet();
#endif
}

// Bzet_new(startbit, len)
// Pretty inefficient, but it does the job
BZET_PTR BZET_FUNC(new)(int64_t startbit, int64_t len) {
    if (len <= 0 || startbit < 0)
        return NULL;

    // Create Bzet with startbit set
    BZET_PTR b = BZET_FUNC(new)(startbit);
    if (!b)
        return NULL;

    // OR in the other bits
    for (int64_t i = startbit + 1; i < startbit + len; i++) {
        // Create new bzet with bit i set
        BZET_PTR next = BZET_FUNC(new)(i);
        if (!next) {
            //TODO: some error message
            BZET_FUNC(destroy)(b);
            return NULL;
        }

        // Get a new bzet with the new bit ORed in
        BZET_PTR temp = BZET_FUNC(OR)(b, next);
        if (!next) {
            //TODO: some error message
            BZET_FUNC(destroy)(next);
            BZET_FUNC(destroy)(b);
            return NULL;
        }

        // Free unnecessary bzets
        BZET_FUNC(destroy)(next);
        BZET_FUNC(destroy)(b);

        // Set b to the correct bzet
        b = temp;
    }

    return b;
}

// Bzet_destroy(b)
void BZET_FUNC(destroy)(BZET_PTR b) {
    if (!b)
        return;

    free(b->bzet);
    free(b->step);
    free(b);
}

// Bzet_clone(b)
BZET_PTR BZET_FUNC(clone)(BZET_PTR b) {
    if (!b)
        return NULL;

    // Allocate Bzet
    BZET_PTR copy = (BZET_PTR) malloc(sizeof(BZET));
    if (!copy)
        return NULL;

    // Copy contents over
    copy->depth = b->depth;
    copy->nbufhalfnodes = b->nbufhalfnodes;
    copy->nhalfnodes = b->nhalfnodes;

    // Deep copy of bzet and step
    copy->bzet = (halfnode_t *) malloc(copy->nbufhalfnodes * sizeof(halfnode_t));
    if (!copy->bzet) {
        free(copy);
        return NULL;
    }

    copy->step = (unsigned char *) malloc(copy->nbufhalfnodes * sizeof(unsigned char));
    if (!copy->step) {
        free(copy->bzet);
        free(copy);
        return NULL;
    }

    memcpy(copy->bzet, b->bzet, copy->nhalfnodes * sizeof(halfnode_t));
    memcpy(copy->step, b->step, copy->nhalfnodes * sizeof(unsigned char));

    return copy;
}


// Operators

// Bzet_setequal(left, right)
void BZET_FUNC(setequal)(BZET_PTR left, BZET_PTR right) {
    if (!left || !right)
        return;

    // Allocate new bzet and step
    halfnode_t *bzet_temp = (halfnode_t *) malloc(right->nbufhalfnodes * sizeof(halfnode_t));
    if (!bzet_temp)
        return;

    unsigned char *step_temp = (unsigned char *) malloc(right->nbufhalfnodes * sizeof(unsigned char));
    if (!step_temp) {
        free(bzet_temp);
        return;
    }

    // Free existing bzet and step
    free(left->bzet);
    free(left->step);

    // Copy contents of bzet and step over
    memcpy(left->bzet, right->bzet, left->nhalfnodes * sizeof(halfnode_t));
    memcpy(left->step, right->step, left->nhalfnodes * sizeof(unsigned char));

    // Copy other fields over
    left->depth = right->depth;
    left->nbufhalfnodes = right->nbufhalfnodes;
    left->nhalfnodes = right->nhalfnodes;
}

// Bzet_NOT(b)
BZET_PTR BZET_FUNC(NOT)(BZET_PTR b) {
    if (!b)
        return NULL;

    // Create a clone of b
    BZET_PTR not = BZET_FUNC(clone)(b);
    if (!not)
        return NULL;

    // NOT it in place
    BZET_FUNC(INVERT)(not);

    return not;
}

// Bzet_INVERT(b)
void BZET_FUNC(INVERT)(BZET_PTR b) {
    if (!b)
        return;

    subtreeNot(b, 0, b->depth);
}

// Bzet_OR(left, right)
BZET_PTR BZET_FUNC(OR)(BZET_PTR left, BZET_PTR right);

// Bzet_AND(left, right)
BZET_PTR BZET_FUNC(AND)(BZET_PTR left, BZET_PTR right);

// Bzet_XOR(left, right)
BZET_PTR BZET_FUNC(XOR)(BZET_PTR left, BZET_PTR right);

// Bzet_COMPARE(left, right)
BZET_PTR BZET_FUNC(COMPARE)(BZET_PTR left, BZET_PTR right);

// Bzet_binop(left, right, op)
BZET_PTR BZET_FUNC(binop)(BZET_PTR left, BZET_PTR right);


// Bzet_TEST(b, bit)
bool BZET_FUNC(TEST)(BZET_PTR b, int64_t bit) {
    // Test if a bit is set by ANDing the bzet with a bzet with only bit set
    // b & Bzet(bit)

    // Create temporary bzets
    BZET_PTR test_bzet = BZET_FUNC(new)(bit);
    if (!test_bzet) {
        // TODO: Some error message
        return false;
    }

    BZET_PTR result_bzet = BZET_FUNC(AND)(b, test_bzet);
    if (!result_bzet) {
        // TODO: Some error message
        BZET_FUNC(destroy)(test_bzet);
        return false;
    }
    
    // Cache result
    bool ret = BZET_FUNC(EMPTY)(result_bzet);

    // Free temporary bzets
    BZET_FUNC(destroy)(test_bzet);
    BZET_FUNC(destroy)(result_bzet);

    return ret;
}

// Bzet_RANGE(b, start, len)
void BZET_FUNC(RANGE)(BZET_PTR b, int64_t start, int64_t len);

// Bzet_SET(bit)
void BZET_FUNC(SET)(BZET_PTR b, int64_t bit) {
    // Set a bit by ORing a bzet with bit set into b
    // b | Bzet(bit)

    // Create temporary bzets
    BZET_PTR temp_bzet = BZET_FUNC(new)(bit);
    if (!temp_bzet) {
        // TODO: Some error message
        return;
    }

    BZET_PTR result_bzet = BZET_FUNC(OR)(b, temp_bzet);
    if (!result_bzet) {
        // TODO: Some error message
        BZET_FUNC(destroy)(temp_bzet);
        return;
    }

    // Shallow copy result to b
    free(b->bzet);
    free(b->step);

    memcpy(result_bzet, b, sizeof(*result_bzet));

    free(result_bzet);

    // Free temporary bzet
    BZET_FUNC(destroy)(temp_bzet);
}

// Bzet_UNSET(bit)
void BZET_FUNC(UNSET)(BZET_PTR b, int64_t bit) {
    // Set a bit by ANDing b with the negation of a bzet with bit set
    // b & ~Bzet(bit)

    // Create temporary bzets
    BZET_PTR temp_bzet = BZET_FUNC(new)(bit);
    if (!temp_bzet) {
        // TODO: Some error message
        return;
    }

    BZET_PTR notb = BZET_FUNC(NOT)(b);
    if (!notb) {
        // TODO: Some error message
        BZET_FUNC(destroy)(temp_bzet);
        return;
    }

    BZET_PTR result_bzet = BZET_FUNC(AND)(temp_bzet, notb);
    if (!result_bzet) {
        // TODO: Some error message
        BZET_FUNC(destroy)(temp_bzet);
        BZET_FUNC(destroy)(notb);
        return;
    }

    // Shallow copy result to b
    free(b->bzet);
    free(b->step);

    memcpy(result_bzet, b, sizeof(*result_bzet));

    free(result_bzet);

    // Free temporary bzets
    BZET_FUNC(destroy)(temp_bzet);
    BZET_FUNC(destroy)(notb);
}

// Bzet_FIRST(b)
int64_t BZET_FUNC(FIRST)(BZET_PTR b);

// Bzet_LAST(b)
int64_t BZET_FUNC(LAST)(BZET_PTR b);

// Bzet_COUNT(b)
int64_t BZET_FUNC(COUNT)(BZET_PTR b);


// Other exported functions

// Bzet_LEV(b)
int BZET_FUNC(LEV)(BZET_PTR b) {
    return (int) b->depth;
}

// Bzet_size(b)
size_t BZET_FUNC(size)(BZET_PTR b) {
    return (b->nhalfnodes * sizeof(b->nhalfnodes) + 1);
}

// Bzet_HEX(b)
void BZET_FUNC(HEX)(BZET_PTR b) {
    _printBzet(b, 0, stdout, b->depth);
}

// Bzet_repr(b, target)
void BZET_FUNC(repr)(BZET_PTR b, void *target) {
    unsigned char *chartarget = (unsigned char *) target;
    chartarget[0] = b->depth;
    memcpy(chartarget + 1, b->bzet, b->nhalfnodes * sizeof(b->nhalfnodes));
}

// Bzet_CLEAN(b)
void BZET_FUNC(CLEAN)(BZET_PTR b) {
    resize(b, 0);
    b->depth = 0;
}

// Bzet_EMPTY(b)
bool BZET_FUNC(EMPTY)(BZET_PTR b) {
    return (b->nhalfnodes == 0);
}

// Bzet_getBits(b, bits, limit, start)
int64_t BZET_FUNC(getBits)(BZET_PTR b, int64_t* bits, int64_t limit, int64_t start);


// Auxiliary functions

NODETYPE _binop(BZET_PTR left, BZET_PTR right, OP op, int lev, size_t left_loc, size_t right_loc, size_t loc);

void _printBzet(BZET_PTR b, int stdOffset, FILE* target, int depth, size_t loc, int offset, bool pad) {
    // Print level info
    if (loc == 0) {
        fprintf(target, "%.2XL", b->depth);

        // If empty, we're done.
        if (BZET_FUNC(EMPTY)(b))
            fprintf(target, "\n");
    }

    // No reading past the bzet!
    if (loc >= b->nhalfnodes)
        return;

    halfnode_t data_bits = b->bzet[loc];
    halfnode_t tree_bits = b->bzet[loc + 1];

    // Print offset if any
    if (pad) {
        // First 3 spaces to align with XXL level byte
        // (sizeof(halfnode_t) + 3)*offset for spaces an internal node takes up
        //     sizeof(halfnode_t) - number of bytes printed out
        //     3 - "pretty" part: [?-?]
        int blanks = 3 + (sizeof(halfnode_t) + 3) * offset + stdOffset;
        for (int j = 0; j < blanks; j++)
            fprintf(target, " ");
    }

    // If depth is 0 or no tree bits are set, this is a data node
    if (depth == 0 || !tree_bits) {
        fprintf(target, "D(%.*X)\n", sizeof(halfnode_t), data_bits);
    }
    // Otherwise this is a tree node
    else {     
        // Print the current node
        fprintf(target, "[%.*X-%.*X]", sizeof(halfnode_t), data_bits, 
            sizeof(halfnode_t), tree_bits);

        // Recursively print subtrees
        bool firstNode = true;
        depth--;
        for (int i = NODE_ELS - 1; i >= 0; i--) {
            // TODO: Use popcount
            // If tree bit set
            if ((tree_bits >> i) & 0x1) {
                // Print first node without offset
                if (firstNode) {
                    _printBzet(b, stdOffset, target, depth, loc + 1, offset + 1);
                    firstNode = false;
                } 
                else {
                    _printBzet(b, stdOffset, target, depth, stepThrough(b, loc + 1), offset + 1, true);
                    loc = stepThrough(b, loc + 1) - 1;
                }
            }
        }
    }
}

void align(BZET_PTR b1, BZET_PTR b2);
void normalize(BZET_PTR b);

size_t stepThrough(BZET_PTR b, size_t loc) {
    // Make sure loc is in range
    if (loc >= b->nhalfnodes)
        return -1;

    // Get step offset
    unsigned char step_offset = b->step[loc];

    // If the step offset is 0, the offset is too large to be actually stored
    // Compute it by examining the step of subtrees stemming from this node
    if (step_offset == 0) {
        // Get node tree bits
        unsigned char tree_bits = b->bzet[loc + 1];

        // Advance to location of first subtree node
        loc++;

        // Step through each subtree
        // TODO: Performance gain by using popcount
        for (int i = NODE_ELS - 1; i >= 0; i--) {
            if ((tree_bits >> i) & 1)
                loc = stepThrough(b, loc);
        }

        return loc;
    }

    // Step offset is nonzero, just add it loc and return it
    return loc + step_offset;
}

void subtreeNot(BZET_PTR b, size_t loc, int depth) {
    // Check if loc is in range
    if (loc >= b->nhalfnodes)
        return;

    // For level 0 nodes, simply NOT it
    if (depth == 0) {
        b->bzet[loc] = ~b->bzet[loc];
        return;
    }

    // Handling all other nodes

    // Get data and tree bits
    halfnode_t data_bits = b->bzet[loc];
    halfnode_t tree_bits = b->bzet[loc + 1];

    // Make result data bits by bitwise NOT and dusting with tree bits
    // and replace old data bits with it
    b->bzet[loc] = ~data_bits & ~tree_bits;

    // Repeat recursively for all subtrees
    loc++;
    // TODO: Performance gain by using popcount for tree bits
    for (int i = 0; i < NODE_ELS; i++) {
        if ((tree_bits >> i) & 1)
            subtreeNot(b, loc, depth - 1);

        loc = stepThrough(b, loc);
    }
}

#endif