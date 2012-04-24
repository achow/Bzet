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

#ifndef BZET
#define BZET PASTE(Bzet, NODE_ELS)
#endif

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
enum NODETYPE { SATURATED, EMPTY, NORMAL, LITERAL };

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
    unsigned char* step; //points to an array that holds step_through values
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
EXPORT_TAGS bool BZET_FUNC(COMPARE)(BZET_PTR left, BZET_PTR right);

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

// Creates a bzet from specified data
BZET_PTR bitstobzet(void *data, size_t len);
void treetobits(unsigned char *buf, halfnode_t *node, int depth);

NODETYPE _binop(BZET_PTR result, BZET_PTR left, BZET_PTR right, OP op, int lev, size_t left_loc = 1, size_t right_loc = 1);

// Recursively print the bzet in "pretty print"
void _printBzet(BZET_PTR b, int stdOffset, FILE* target, int depth, size_t loc = 0, int offset = 0, bool pad = 0);

// Align two Bzets to the same level
void align(BZET_PTR b1, BZET_PTR b2);

// Strip leading unnecessary nodes
void normalize(BZET_PTR b);

// Return the location of the next node after traversing the subtree starting at loc
size_t step_through(BZET_PTR b, size_t loc);

// In-place bitwise NOT of the subtree whose root is at loc
void subtree_not(BZET_PTR b, size_t loc, int depth);


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
// TODO: Add error messages
inline 
BZET_PTR init(size_t initial_alloc = INITIAL_ALLOC) {
    // Allocate bzet struct
    BZET_PTR b = (BZET_PTR) malloc(sizeof(BZET));
    if (!b)
        return NULL;

    // Allocate bzet node array
    b->bzet = (halfnode_t*) malloc(initial_alloc * sizeof(halfnode_t));
    if (!b->bzet) {
        free(b);
        return NULL;
    }

    // Allocate step array
    b->step = (unsigned char*) malloc(initial_alloc * sizeof(unsigned char));
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

// Append a subtree from src to dest starting at loc in src
inline
void append_subtree(BZET_PTR dst, BZET_PTR src, size_t loc) {
    if (!dst || !src)
        return;

    // Calculate copy size and cache copy destination
    size_t copy_size = step_through(src, loc) - loc;
    size_t dst_loc = dst->nhalfnodes;

    // Resize dst to accomodate copy_size new elements
    resize(dst, dst->nhalfnodes + copy_size);

    // Do copy
    memcpy(dst + dst_loc, src + loc, copy_size * sizeof(halfnode_t));
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

        void append_subtree(const Bzet4& src, size_t loc);
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
            nextLoc = step_through(nextLoc + 1) - 1;
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
    while (POW(depth + 1) <= (size_t) bit)
        depth++;

    // Resize to accomodate full bzet
    resize(b, 2*depth + 1);

    //set depth byte
    b->depth = (unsigned char) depth;

    // Set level 0 data node
    // Each depth contains 2 halfnodes, b->bzet[2*depth] is the data node
    b->bzet[2*depth] = (halfnode_t) 0x1 << (NODE_ELS - 1 - (bit % NODE_ELS));

    // Zero out other nodes, since we will be only filling the tree portions next
    memset(b->bzet, 0x00, 2*depth * sizeof(halfnode_t));

    // Set tree bits, working linearly from the node at the head of the bzet
    for (int i = depth; i >= 1; i--) {
        // Get weight of bits at this level
        size_t cpow = POW(i);

        // Calculate tree bit to set
        int setbit = (int) (bit / cpow) % NODE_ELS;

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

    return b;
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

    // Resize left's bzet and step buffers
    resize(left, right->nhalfnodes);

    // Copy contents of bzet over
    left->depth = right->depth;
    left->nbufhalfnodes = right->nbufhalfnodes;
    left->nhalfnodes = right->nhalfnodes;
    memcpy(left->bzet, right->bzet, left->nhalfnodes * sizeof(halfnode_t));
    memcpy(left->step, right->step, left->nhalfnodes * sizeof(unsigned char));
}

// Bzet_NOT(b)
BZET_PTR BZET_FUNC(NOT)(BZET_PTR b) {
    if (!b)
        return NULL;

    // Create a clone of b
    BZET_PTR n = BZET_FUNC(clone)(b);
    if (!n)
        return NULL;

    // NOT it in place
    BZET_FUNC(INVERT)(n);

    return n;
}

// Bzet_INVERT(b)
void BZET_FUNC(INVERT)(BZET_PTR b) {
    if (!b)
        return;

    subtree_not(b, 0, b->depth);
}

// Bzet_OR(left, right)
BZET_PTR BZET_FUNC(OR)(BZET_PTR left, BZET_PTR right) {
    if (!left || !right)
        return NULL;

    // If left bzet is empty, the bitwise OR will be equal to the right bzet
    if (left->nhalfnodes == 0)
        return BZET_FUNC(clone)(right);
    // If right bzet is empty, the bitwise OR will be equal to the left bzet
    else if (right->nhalfnodes == 0)
        return BZET_FUNC(clone)(left);
    // Otherwise just operate on them
    else {
        align(left, right);
        BZET_PTR result = BZET_FUNC(binop)(left, right, OR);
        normalize(left);
        normalize(right);
        return result;
    }
}

// Bzet_AND(left, right)
BZET_PTR BZET_FUNC(AND)(BZET_PTR left, BZET_PTR right) {
    if (!left || !right)
        return NULL;

    // If either bzet is empty, the bitwise AND will be an empty bzet
    if (left->nhalfnodes == 0 || right->nhalfnodes == 0)
        return BZET_FUNC(new)();
    // Otherwise just operate on them
    else {
        align(left, right);
        BZET_PTR result = BZET_FUNC(binop)(left, right, AND);
        normalize(left);
        normalize(right);
        return result;
    }
}

// Bzet_XOR(left, right)
BZET_PTR BZET_FUNC(XOR)(BZET_PTR left, BZET_PTR right) {
    if (!left || !right)
        return NULL;

    // If left bzet is empty, the bitwise XOR will be equal to the right bzet
    if (left->nhalfnodes == 0)
        return BZET_FUNC(clone)(right);
    // If right bzet is empty, the bitwise XOR will be equal to the left bzet
    else if (right->nhalfnodes == 0)
        return BZET_FUNC(clone)(left);
    // Otherwise just operate on them
    else {
        align(left, right);
        BZET_PTR result = BZET_FUNC(binop)(left, right, XOR);
        normalize(left);
        normalize(right);
        return result;
    }
}

// Bzet_COMPARE(left, right)
bool BZET_FUNC(COMPARE)(BZET_PTR left, BZET_PTR right) {
    if (left->depth != right->depth || left->nhalfnodes != right->nhalfnodes ||
        memcmp(left->bzet, right->bzet, left->nhalfnodes * sizeof(halfnode_t)))
        return false;

    return true;
}

// Bzet_binop(left, right, op)
BZET_PTR BZET_FUNC(binop)(BZET_PTR left, BZET_PTR right, OP op) {
    if (!left || !right)
        return NULL;

    BZET_PTR result = BZET_FUNC(new)();
    _binop(result, left, right, op, left->depth);
    return result;
}


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
void BZET_FUNC(RANGE)(BZET_PTR b, int64_t start, int64_t len) {
    if (!b)
        return;

    // Create bzet mask
    BZET_PTR mask = BZET_FUNC(new)();
    if (!mask)
        // TODO: Add error message
        return;

    for (int64_t i = start; i < start + len; i++)
        BZET_FUNC(SET)(mask, start);

    // OR the mask into the original bzet
    BZET_PTR result = BZET_FUNC(OR)(b, mask);
    if (!result) {
        // TODO: Add error message
        BZET_FUNC(destroy)(mask);
        return;
    }
    
    // Free the mask
    BZET_FUNC(destroy)(mask);

    // Gut result and give its internals to b
    free(b->bzet);
    free(b->step);
    memcpy(b, result, sizeof(BZET));
    free(result);
}

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

/*BZET_PTR bitstobzet(void *data, size_t len) {
    unsigned char *bytes = (unsigned char *) data;

    // Find depth required to store len bytes
    int depth = 0;
    while (len * 8 < PASTE(pow, NODE_ELS)(depth + 1))
        depth++;

    // Create a bzet with initial alloc of maximum required buffers
    // b->nhalfnodes at this point is 0
    BZET_PTR b = init(depth / sizeof(halfnode_t) + 1);

    // Set actual depth
    b->depth = depth;


}*/

NODETYPE bitstotree(BZET_PTR b, int depth, unsigned char *data, size_t len) {
    if (len < PASTE(pow, NODE_ELS)(depth) / 8)
        return EMPTY;

    // If a level 0 node is to be created
    if (depth == 0) {
        // Build a level 0 data node
        halfnode_t data_node = 0;
        for (int i = 0; i < sizeof(halfnode_t); i++)
            data_node = (data_node << 8) | data[i];

        // If data node is saturated
        if (data_node == (halfnode_t) -1)
            return SATURATED;
        // If data node is empty
        else if (data_node == 0)
            return EMPTY;
        // Otherwise commit to bzet
        else {
            b->bzet[b->nhalfnodes] = data_node;
            b->step[b->nhalfnodes] = 1;
            b->nhalfnodes++;
            return NORMAL;
        }
    }

    // Create nodes
    halfnode_t data_bits = 0, tree_bits = 0;
    
    // Reserve space for nodes and cache location
    size_t loc = b->nhalfnodes;
    resize(b, b->nhalfnodes + 2);

    // Build subtrees
    size_t subtree_size = PASTE(pow, NODE_ELS)(depth - 1) / 8;
    for (int i = 0; i < NODE_ELS; i++) {
        // Recurse
        NODETYPE ret = bitstotree(b, depth - 1, data, len);

        // Modify data and len as necessary
        data += subtree_size;
        len -= subtree_size;

        // If subtree built was empty
        if (ret == EMPTY) {
            // Nothing to do, bits are already set to 0
        }
        // If subtree was saturated
        else if (ret == SATURATED) {
            // Set data bit
            data_bits = (data_bits << 1) | 0x1;
        }
        // If subtree was a data literal
        else if (ret == LITERAL) {
            // Set data and tree bit
            data_bits = (data_bits << 1) | 0x1;
            tree_bits = (tree_bits << 1) | 0x1;
        }
        // Otherwise it is a normal tree node
        else {
            // Set tree bit
            tree_bits = (tree_bits << 1) | 0x1;
        }
    }

    // If a literal subtree is smaller than the formed tree
    if ((b->nhalfnodes - loc) * sizeof(halfnode_t) >= PASTE(pow, NODE_ELS)(depth + 1)) {
        size_t locstart = loc;
        unsigned char *buf = (unsigned char *) malloc(PASTE(pow, NODE_ELS)(depth + 1) * sizeof(halfnode_t));
        // TODO: More elegant error handling
        if (!buf)
            display_error("malloc failed in bitstotree", true);

        // Convert tree to bytes
        treetobits(buf, b->bzet + loc, depth);

        unsigned char *bufptr = buf;

        // Write back as necessary
        // For each halfnode's worth
        for (int i = 0; i < PASTE(pow, NODE_ELS)(depth + 1) / sizeof(halfnode_t); i++) {
            // Copy to satisfy byte ordering
            halfnode_t node = 0;
            for (int j = 0; j < sizeof(halfnode_t); j++) {
                node |= bufptr[0];
                bufptr++;
            }
            b->bzet[loc] = node;
            loc++;
        }

        // Shrink size
        resize(b, loc);

        // Set step
        if (b->nhalfnodes - locstart > 255)
            b->step[locstart] = 0;
        else
            b->step[locstart] = (unsigned char) (b->nhalfnodes - locstart);

        return LITERAL;
    }

    // Write nodes
    b->bzet[loc] = data_bits;
    b->bzet[loc + 1] = tree_bits;

    // Set step
    if (b->nhalfnodes - loc > 255)
        b->step[loc] = 0;
    else
        b->step[loc] = (unsigned char) (b->nhalfnodes - loc);

    return NORMAL;
}

void treetobits(unsigned char *buf, halfnode_t *node, int depth) {
    // If depth is 0, we already have a data literal
    if (depth == 0) {
        // Work byte by byte to write to buffer in big endian
        halfnode_t data_node = node[0];
        for (int i = 0; i < sizeof(halfnode_t); i++) {
            buf[i] = data_node & ((halfnode_t) -1);
            data_node <<= 8;
        }

        return;
    }

    // Otherwise we have a subtree to traverse
    halfnode_t data_bits = node[0];
    halfnode_t tree_bits = node[1];
    node += 2;

    // Work through each subtree
    for (int i = 0; i < sizeof(halfnode_t) * 8; i++) {
        bool data_bit = (data_bits << i) & 0x1;
        bool tree_bit = (tree_bits << i) & 0x1;

        // If there is a subtree
        if (tree_bit && !data_bit)
            treetobits(buf, node, depth - 1);
        // If there is a data literal subtree
        else if (tree_bit && data_bit) {
            // Work byte by byte to write to buffer in big endian
            size_t nodes = PASTE(pow, NODE_ELS)(depth - 1) / 8 / sizeof(halfnode_t);
            // For each node
            for (int j = 0; j < nodes; j++) {
                halfnode_t data = node[j];
                // For each byte
                for (int k = 0; k < sizeof(halfnode_t); k++) {
                    buf[0] = data & 0xFF;
                    buf++;
                    data >>= 8;
                }
            }
        }
        // If this element is saturated
        else if (data_bit) {
            memset(buf, 0xFF, PASTE(pow, NODE_ELS)(depth - 1) / 8);
            buf += PASTE(pow, NODE_ELS)(depth - 1) / 8;
        }
        // Otherwise this element is empty
        else {
            memset(buf, 0x00, PASTE(pow, NODE_ELS)(depth - 1) / 8);
            buf += PASTE(pow, NODE_ELS)(depth - 1) / 8;
        }
    }
}

NODETYPE _binop(BZET_PTR result, BZET_PTR left, BZET_PTR right, OP op, 
                int lev, size_t left_loc, size_t right_loc) {
    // Get corresponding nodes of the left and right tree
    halfnode_t c_left_data = left->bzet[left_loc];
    halfnode_t c_left_tree = left->bzet[left_loc + 1];
    halfnode_t c_right_data = right->bzet[right_loc];
    halfnode_t c_right_tree = right->bzet[right_loc];

    // Create new node to hold data
    halfnode_t node_data = 0;
    halfnode_t node_tree = 0;

    // Reserve space for this node
    size_t loc = result->nhalfnodes;
    resize(result, result->nhalfnodes + 2);

    // For each element in the node
    for (int i = NODE_ELS - 1; i >= 0; --i) {
        int cur_left_tree_bit = (c_left_tree >> i) & 0x1;
        int cur_left_data_bit = (c_left_data >> (i + NODE_ELS)) & 0x1;
        int cur_right_tree_bit = (c_right_tree >> i) & 0x1;
        int cur_right_data_bit = (c_right_data >> (i + NODE_ELS)) & 0x1;
 
        // TT: If both tree bits are on
        if (cur_left_tree_bit && cur_right_tree_bit) {
            // TODO: Handling for data literal operations
            // If both corresponding subtrees are data literals
            if (cur_left_data_bit && cur_left_tree_bit) {
            }
            // If left subtree is a data literal and right is not
            else if (cur_left_data_bit) {
            }
            // If right subtree is a data literal and left is not
            else if (cur_right_data_bit) {
            }
            // Otherwise both subtrees are actual trees
            else {
                // Recurse
                NODETYPE cn = _binop(result, left, right, op, lev - 1, left_loc + 2, right_loc + 2);

                // Saturated subtree
                if (cn == SATURATED) {
                    // Turn on data bit
                    node_data = (node_data << 1) | 0x1;
                    node_tree <<= 1;
                }
                // Empty subtree
                else if (cn == EMPTY) {
                    // Shift data and tree nodes over
                    node_data <<= 1;
                    node_tree <<= 1;
                }
                // "Tree" subtree
                else if (cn == NORMAL) {
                    // Turn on tree bit
                    node_data <<= 1;
                    node_tree = (node_tree << 1) | 0x1;
                }
                // Otherwise data literal subtree
                else {
                    // TODO: Handle data literal subtree
                }

                // Advance location pointers
                left_loc = step_through(left, left_loc + 2) - 2;
                right_loc = step_through(right, right_loc + 2) - 2;
            }
        } 
        // ?T or T?: If only one of the tree bits are on
        else if (cur_left_tree_bit || cur_right_tree_bit) {
            // Look up action from optable
            ACTION action;
            // ?T
            if (cur_right_tree_bit) {
                // 1T
                if (cur_left_data_bit)
                    action = optable[(op << 2) + 2];
                // 0T
                else
                    action = optable[op << 2];
            }
            // T?
            else {
                // T1
                if (cur_right_data_bit)
                    action = optable[(op << 2) + 3];
                // T0
                else
                    action = optable[(op << 2) + 1];
            }

            // Used for NA and NB actions
            size_t end;

            // Execute action
            switch (action) {
                // Delete left subtree, set data bit off
                case DA0:
                    // Skip the left subtree
                    left_loc = step_through(left, left_loc + 2) - 2;

                    // Data bit is already 0
                    node_data <<= 1;
                    node_tree <<= 1;

                    break;

                // Delete right subtree, set data bit off
                case DB0:
                    // Skip the right subtree
                    right_loc = step_through(right, right_loc + 2) - 2;

                    // Data bit is already 0
                    node_data <<= 1;
                    node_tree <<= 1;

                    break;

                // Delete left subtree, set data bit on
                case DA1: 
                    // Skip the left subtree
                    left_loc = step_through(left, left_loc + 2) - 2;

                    // Turn on data bit
                    node_data = (node_data << 1) | 0x1;
                    node_tree <<= 1;

                    break;

                // Delete right subtree, set data bit on
                case DB1:
                    // Skip the right subtree
                    right_loc = step_through(right, right_loc + 2) - 2;

                    // Turn on data bit
                    node_data = (node_data << 1) | 0x1;
                    node_tree <<= 1;

                    break;

                // Copy left subtree into result
                case CA: 
                    // Append left subtree
                    append_subtree(result, left, left_loc + 2);

                    // Move through left subtree
                    left_loc = step_through(left, left_loc + 2) - 2;

                    // Turn on tree bit
                    node_data <<= 1;
                    node_tree = (node_tree << 1) | 0x1;

                    break;

                // Copy right subtree into result
                case CB:
                    // Append left subtree
                    append_subtree(result, left, left_loc + 2);

                    // Move through left subtree
                    left_loc = step_through(left, left_loc + 2) - 2;

                    // Turn on tree bit
                    node_data <<= 1;
                    node_tree = (node_tree << 1) | 0x1;

                    break;

                // Copy left subtree into result and negate
                case NA:
                    end = result->nhalfnodes;

                    // Append left subtree
                    append_subtree(result, left, left_loc + 2);

                    // Negate
                    subtree_not(result, end, lev - 1);

                    // Move through left subtree
                    left_loc = step_through(left, left_loc + 2) - 2;

                    // Turn on tree bit
                    node_data <<= 1;
                    node_tree = (node_tree << 1) | 0x1;

                    break;

                // Copy right subtree into result and negate
                case NB:
                    end = result->nhalfnodes;

                    // Append right subtree
                    append_subtree(result, right, right_loc + 2);

                    // Negate
                    subtree_not(result, end, lev - 1);

                    // Move through left subtree
                    left_loc = step_through(right, right_loc + 2) - 2;

                    // Turn on tree bit
                    node_data <<= 1;
                    node_tree = (node_tree << 1) | 0x1;
                    break;

                default:
                    // Should be impossible to get here, but just in case?
                    display_error("Bzet4::_binop: Something went terribly, terribly wrong", true);
                    break;
            }
        } 
        // Only data bits
        else /*if (!cur_left_tree_bit && !cur_right_tree_bit)*/ {
            // Shift in data bit
            node_data = (node_data << 1) | do_data_op(op, cur_left_data_bit, cur_right_data_bit);
            node_tree <<= 1;
        }
    }

    // Write nodes
    result->bzet[loc] = node_data;
    result->bzet[loc + 1] = node_tree;

    // Set step
    if (result->nhalfnodes - loc > 255)
        result->step[loc] = 0;
    else 
        result->step[loc] = (unsigned char) (result->nhalfnodes - loc);

    // If resulting node is empty
    if (node_tree == 0 && node_data == 0) {
        // Drop newly committed node if not root node
        if (result->nhalfnodes > 2)
            resize(result, result->nhalfnodes - 2);
        return EMPTY;
    }
    // If resulting node is saturated
    else if (node_tree == 0 && node_data == (halfnode_t) -1) {
        // Drop newly committed node if not root node
        if (result->nhalfnodes > 2)
            resize(result, result->nhalfnodes - 2);
        return SATURATED;
    }

    // TODO: See if collapsible to bit literal

#if (defined _DEBUG || defined DEBUG)
    assert(dust(m_bzet[loc]) == m_bzet[loc]);
#endif

    return NORMAL;
}

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
        int blanks = 3 + (sizeof(halfnode_t) * 2 + 3) * offset + stdOffset;
        for (int j = 0; j < blanks; j++)
            fprintf(target, " ");
    }

    // If depth is 0 or no tree bits are set, this is a data node
    if (depth == 0 || !tree_bits) {
        fprintf(target, "D(%.*X)\n", sizeof(halfnode_t) * 2, data_bits);
    }
    // Otherwise this is a tree node
    else {     
        // Print the current node
        fprintf(target, "[%.*X-%.*X]", sizeof(halfnode_t) * 2, data_bits, 
            sizeof(halfnode_t) * 2, tree_bits);

        // Recursively print subtrees
        bool firstNode = true;
        depth--;
        for (int i = NODE_ELS - 1; i >= 0; i--) {
            // TODO: Use popcount
            // If tree bit set
            if ((tree_bits >> i) & 0x1) {
                // Print first node without offset
                if (firstNode) {
                    _printBzet(b, stdOffset, target, depth, loc + 2, offset + 1);
                    firstNode = false;
                } 
                else {
                    _printBzet(b, stdOffset, target, depth, step_through(b, loc + 2), offset + 1, true);
                    loc = step_through(b, loc + 2) - 1;
                }
            }
        }
    }
}

void align(BZET_PTR b1, BZET_PTR b2) {
    if (!b1 || !b2)
        return;

    // Get difference in depths
    int diffdepth = abs(b1->depth - b2->depth);

    // Nothing to be done if they are the same depth
    if (diffdepth == 0)
        return;

    // If b2 needs to be grown
    if (b1->depth > b2->depth) {
        b2->depth = b1->depth;
        size_t old_size = b2->nhalfnodes;

        // Resize b2 to accommodate new heading nodes
        resize(b2, b2->nhalfnodes + diffdepth * 2);

        // Move bzet and step in b2 to accommodate new heading nodes
        memmove(b2->bzet + diffdepth * 2, b2->bzet, old_size * sizeof(halfnode_t));
        memmove(b2->step + diffdepth * 2, b2->step, old_size * sizeof(halfnode_t));

        // Add new nodes and new step
        size_t loc = 0;
        for (int i = 0; i < diffdepth; i++) {
            b2->bzet[loc] = 0;
            b2->bzet[loc + 1] = (halfnode_t) (0x1 << (NODE_ELS - 1));

            if (b2->nhalfnodes - loc > 255)
                b2->step[loc] = 0;
            else
                b2->step[loc] = (unsigned char) (b2->nhalfnodes - loc);

            loc += 2;
        }
    }
    // Otherwise b1 needs to be grown
    else {
        b1->depth = b2->depth;
        size_t old_size = b1->nhalfnodes;

        // Resize b1 to accommodate new heading nodes
        resize(b1, b1->nhalfnodes + diffdepth * 2);

        // Move bzet and step in b1 to accommodate new heading nodes
        memmove(b1->bzet + diffdepth * 2, b1->bzet, old_size * sizeof(halfnode_t));
        memmove(b1->step + diffdepth * 2, b1->step, old_size * sizeof(halfnode_t));

        // Add new nodes and new step
        size_t loc = 0;
        for (int i = 0; i < diffdepth; i++) {
            b1->bzet[loc] = 0;
            b1->bzet[loc + 1] = (halfnode_t) (0x1 << (NODE_ELS - 1));

            if (b1->nhalfnodes - loc > 255)
                b1->step[loc] = 0;
            else
                b1->step[loc] = (unsigned char) (b1->nhalfnodes - loc);

            loc += 2;
        }
    }
}

void normalize(BZET_PTR b) {
    if (!b)
        return;

    int leading = 0;
    int loc = 0;
    int depth = b->depth;
    // Count leading nodes that can be stripped.
    // This occurs when the data node is all zero and the leftmost tree bit is set
    while (depth > 0 && b->bzet[loc] == 0 && 
        b->bzet[loc + 1] == (halfnode_t) (0x1 << (NODE_ELS - 1))) {
        leading++;
        loc += 2;
    }

    // If there are leading nodes that can be stripped
    if (leading) {
        // Modify depth
        b->depth -= leading;
        // Move bzet
        memmove(b->bzet, b->bzet + loc, loc * sizeof(b->bzet[0]));
        // Move step
        memmove(b->step, b->step + loc, loc * sizeof(b->step[0]));
        // Resize
        resize(b, b->nhalfnodes - loc);
    }
}

size_t step_through(BZET_PTR b, size_t loc) {
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
                loc = step_through(b, loc);
        }

        return loc;
    }

    // Step offset is nonzero, just add it loc and return it
    return loc + step_offset;
}

void subtree_not(BZET_PTR b, size_t loc, int depth) {
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

    // Make result data bits by bitwise NOT
    b->bzet[loc] = (~data_bits & ~tree_bits) | data_bits & tree_bits;

    // Repeat recursively for all subtrees
    loc += 2;
    // TODO: Performance gain by using popcount for tree bits
    for (int i = 0; i < NODE_ELS; i++) {
        // If tree bit set
        if ((tree_bits >> i) & 1) {
            // If subtree is in tree form
            if (((data_bits >> i) & 1) == 0) {
                subtree_not(b, loc, depth - 1);
                loc = step_through(b, loc);
            }
            // If subtree is in literal form
            else {
                // Get number of halfnodes 
                size_t nnodes = PASTE(pow, NODE_ELS)(depth - 1) / 8 / sizeof(halfnode_t);

                // NOT each halfnode
                for (size_t i = 0; i < nnodes; i++) {
                    b->bzet[loc] = ~b->bzet[loc];
                    loc++;
                }
            }
        }
    }
}

#endif
