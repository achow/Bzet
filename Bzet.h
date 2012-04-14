/*****************************************************************************
 * Generic code for Bzet8/Bzet16/Bzet32
 *
 * Things to do:
 *  - Good error messages, better error handling in general (?)
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
void _printBzet(BZET_PTR b, int stdOffset = 0, FILE* target = stdout, size_t loc = 1, int depth = 0, int offset = 0, bool pad = 0);
void align(BZET_PTR b1, BZET_PTR b2);
void normalize(BZET_PTR b);
size_t stepThrough(BZET_PTR b, size_t loc);
void subtreeNot(BZET_PTR b, size_t loc, int depth = 0);


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
        if (!bzet_temp)
            return;

        // Reallocate step
        unsigned char *step_temp = (unsigned char*) realloc(b->step, b->nbufhalfnodes * sizeof(unsigned char));
        if (!step_temp)
            return;

        b->bzet = bzet_temp;
        b->step = step_temp;

        /*
        // Check that realloc succeeded
        if (!b->bzet || !b->step) {
            fprintf(stderr, "Fatal error: Resizing bzet failed attempting to allocate %d bytes\n", (int) (b->nbufhalfnodes * sizeof(halfnode_t)));
            display_error("", true);
        }*/
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
EXPORT_TAGS BZET_PTR BZET_FUNC(OR)(BZET_PTR left, BZET_PTR right);

// Bzet_AND(left, right)
EXPORT_TAGS BZET_PTR BZET_FUNC(AND)(BZET_PTR left, BZET_PTR right);

// Bzet_XOR(left, right)
EXPORT_TAGS BZET_PTR BZET_FUNC(XOR)(BZET_PTR left, BZET_PTR right);

// Bzet_COMPARE(left, right)
EXPORT_TAGS BZET_PTR BZET_FUNC(COMPARE)(BZET_PTR left, BZET_PTR right);

// Bzet_binop(left, right, op)
EXPORT_TAGS BZET_PTR BZET_FUNC(binop)(BZET_PTR left, BZET_PTR right);


/*****************************************************************************
 * 
 *  Function name: Bzet4
 *
 *  Purpose:       Constructor for Bzet4. Allows for Bzet compression of a
 *                 specified file.
 *
 *  Inputs:        const char* filename: file to compress
 *  Return values: None
 * 
 *  Author:        Alex Chow
 *  Date:          4/10/2012
 *
 ****************************************************************************
/*Bzet4::Bzet4(const char* filename) {
    //initialize Bzet
    init();
    clear();

    //open file
    FILE* f = fopen(filename, "r");

    if (f) {
        //retrieve file size
        fseek(f, 0, SEEK_END);
        int64_t filesize = ftell(f);
        fseek(f, 0, SEEK_SET);

        //aligned size
        int64_t alignedSize = ((filesize & 1) ? filesize : filesize + 1);
        int64_t alignedSizeBits = alignedSize * 8;
        int depth = buildDepth(alignedSizeBits - 1);

        //initial size, counting level 1 nodes
        int64_t bytesreqd = alignedSize / 2;
        alignedSize = alignedSize / 2;

        //each level up, there will be a factor of 4 fewer bytes required
        while (alignedSize > NODE_ELS) {

        }
    }

#if (defined _DEBUG || defined DEBUG)
    validateBzet();
#endif
}

/*****************************************************************************
 * 
 *  Function name: operator~
 *
 *  Purpose:       Bitwise negation operator
 *
 *  Inputs:        None
 *  Return values: Bitwise negation of *this
 * 
 *  Author:        Alex Chow
 *  Date:          11/1/2011
 *
 ****************************************************************************
Bzet4 Bzet4::operator~() const {
    //clone bzet and NOT in place
    Bzet4 not(*this);
    not.subtreeNot(1, (int) b->bzet[0]);

#if (defined _DEBUG || defined DEBUG)
    not.validateBzet();
#endif

    return not;
}

/*****************************************************************************
 * 
 *  Function name:  appendSubtree
 *
 *  Purpose:        Appends subtree starting at loc 
 *
 *  Inputs:         Bzet4& src: source to copy from
 *                  size_t loc: start of the subtree to append
 *  Return values:  None
 * 
 *  Author:         Alex Chow
 *  Date:           12/12/2011
 *
 ****************************************************************************
 void Bzet4::appendSubtree(const Bzet4& src, size_t loc) {
     size_t step = src.stepThrough(loc);

     //get length to copy
     size_t len = src.stepThrough(loc) - loc;
     size_t old_size = b->size;

     //resize bzet to accommodate new nodes
     resize(b->size + len);

     //copy bzet
     memcpy(b->bzet + old_size, src.b->bzet + loc, len);
     
     //copy b->step
     memcpy(b->step + old_size, src.b->step + loc, len);
}

/*****************************************************************************
 * 
 *  Function name:  dropNodes
 *
 *  Purpose:        Drops the specified number of nodes starting from loc
 *
 *  Inputs:         size_t loc: location of first node to drop
 *                  int n: number of nodes to drop
 *                
 *  Return values:  None
 * 
 *  Author:         Alex Chow
 *  Date:           12/12/2011
 *
 ****************************************************************************
void Bzet4::dropNodes(size_t loc, int n) {
#if (defined _DEBUG || defined DEBUG)
    assert(loc > 0 && loc < b->size && n >= 0);
#endif

    memmove(b->bzet + loc, b->bzet + loc + n, b->size - (loc + n));
    memmove(b->step + loc, b->step + loc + n, b->size - (loc + n));

    resize(b->size - n);

#if (defined _DEBUG || defined DEBUG)
    validateBzet();
#endif
}

/*****************************************************************************
 * 
 *  Function name:  do_data_op
 *
 *  Purpose:        Does the specified operation on provided bits
 *
 *  Inputs:         OP op: operation to perform
 *                  int left_data_bit: the left data bit
 *                  int right_data_bit: the right data bit
 *                
 *  Return values:  result of doing op on the provided bits
 * 
 *  Author:         Alex Chow
 *  Date:           12/11/2011
 *
 ****************************************************************************
int Bzet4::do_data_op(OP op, int left_data_bit, int right_data_bit) {
    //use op directly to build result bit
    return (op >> (3 - (((int) left_data_bit << 1) | (int) right_data_bit))) & 0x1;
}

/*****************************************************************************
 * 
 *  Function name: _binop
 *
 *  Purpose:       Actual implementation of binop
 *
 *  Inputs:        Bzet4 left: the left operand
 *                 Bzet4 right: the right operand
 *                 OP op: operation to do
 *                 size_t left_loc: current location in left
 *                 size_t right_loc: current locatin in right
 *                 size_t loc: current location in this bzet
 *  Return values: NODETYPE: type of root node of built subtree (SATURATED, 
 *                           EMPTY, NORMAL)
 * 
 *  Author:        Alex Chow
 *  Date:          12/11/2011
 *
 ****************************************************************************
NODETYPE Bzet4::_binop(const Bzet4& left, const Bzet4& right, OP op, int lev, size_t left_loc, size_t right_loc, size_t loc) {
    //left_loc and right_loc are unmodified until an operation is done
    //so they both point to the corresponding nodes in left and right used to build the current node
    unsigned char c_left = left.b->bzet[left_loc];
    unsigned char c_right = right.b->bzet[right_loc];

    //special handling for level 1
    if (lev == 1) {
        //build data nodes 0 and 1
        unsigned char result_c1 = 0x00;
        for (int i = 7; i >= 0; --i) {
            int bita = (c_left >> i) & 1;
            int bitb = (c_right >> i) & 1;
            result_c1 |= (do_data_op(op, bita, bitb) << i);
        }

        //get corresponding data nodes 2 and 3
        c_left = left.b->bzet[left_loc + 1];
        c_right = right.b->bzet[right_loc + 1];

        //build data nodes 2 and 3
        unsigned char result_c2 = 0x00;
        for (int i = 7; i >= 0; --i) {
            int bita = (c_left >> i) & 1;
            int bitb = (c_right >> i) & 1;
            result_c2 |= (do_data_op(op, bita, bitb) << i);
        }

        //saturated node
        if (result_c1 == 0xFF && result_c2 == 0xFF) {
            //build bzet manually if bzet is saturated
            if (b->size == 1) {
                resize(2);
                b->bzet[0]++;
                b->bzet[1] = 0x80;
                b->step[1] = 1;
                return NORMAL;
            }

            return SATURATED;
        }
        //empty node
        else if (!result_c1 && !result_c2) {
            //build bzet manually if bzet is empty
            if (b->size == 1) {
                clear();
                return NORMAL;
            }

            return EMPTY;
        }

        //otherwise commit data
        //create two nodes to accommodate data nodes
        resize(b->size + 2);

        //write result
        b->bzet[loc] = result_c1;
        b->bzet[loc + 1] = result_c2;

        //set steps
        b->step[loc] = 2;
        b->step[loc + 1] = 1;
        
        return NORMAL;
    }

    //create new node to hold data
    resize(b->size + 1);
    b->bzet[loc] = 0x00;

    //go through each element in the current node
    for (int i = NODE_ELS - 1; i >= 0; --i) {
        int cur_left_tree_bit = (c_left >> i) & 0x1;
        int cur_left_data_bit = (c_left >> (i + NODE_ELS)) & 0x1;
        int cur_right_tree_bit = (c_right >> i) & 0x1;
        int cur_right_data_bit = (c_right >> (i + NODE_ELS)) & 0x1;
 
        //TT: if both tree bits are on
        if (cur_left_tree_bit && cur_right_tree_bit) {
            //turn on tree bit
            //b->bzet[loc] |= 0x01 << i;

            //recurse
            NODETYPE cn = _binop(left, right, op, lev - 1, left_loc + 1, right_loc + 1, b->size);

            if (cn == SATURATED) {
                //saturated subtree, turn on data bit
                b->bzet[loc] |= 0x10 << i;
            }
            else if (cn == EMPTY) {
                //do nothing, elements are already set to 0
            }
            else {
                //subtree exists, turn on tree bit
                b->bzet[loc] |= 0x01 << i;
            }

            //advance location pointers
            left_loc = left.stepThrough(left_loc + 1) - 1;
            right_loc = right.stepThrough(right_loc + 1) - 1;
        } 
        //?T or T?: if only one of the tree bits are on
        else if (cur_left_tree_bit || cur_right_tree_bit) {
            //look up action from optable
            ACTION action;
            //?T
            if (cur_right_tree_bit) {
                //1T
                if (cur_left_data_bit)
                    action = optable[(op << 2) + 2];
                //0T
                else
                    action = optable[op << 2];
            }
            //T?
            else {
                //T1
                if (cur_right_data_bit)
                    action = optable[(op << 2) + 3];
                //T0
                else
                    action = optable[(op << 2) + 1];
            }

            //used for NA and NB actions
            size_t end;

            //execute action
            switch (action) {
                //delete left subtree, set data bit off
                case DA0:
                    //skip the left subtree
                    left_loc = left.stepThrough(left_loc + 1) - 1;

                    //nothing to be done, data bit is already 0
                    break;

                //delete right subtree, set data bit off
                case DB0:
                    //skip the right subtree
                    right_loc = right.stepThrough(right_loc + 1) - 1;

                    //nothing to be done, data bit is already 0
                    break;

                //delete left subtree, set data bit on
                case DA1: 
                    //skip the left subtree
                    left_loc = left.stepThrough(left_loc + 1) - 1;
                    //turn on data bit in the result node
                    b->bzet[loc] |= 0x80 >> ((NODE_ELS - 1) - i);
                    break;

                //delete right subtree, set data bit on
                case DB1:
                    //skip the right subtree
                    right_loc = right.stepThrough(right_loc + 1) - 1;
                    //turn on data bit in the result node
                    b->bzet[loc] |= 0x80 >> ((NODE_ELS - 1) - i);
                    break;

                //copy left subtree into result
                case CA: 
                    //turn on tree bit
                    b->bzet[loc] |= 0x08 >> ((NODE_ELS - 1) - i);

                    //append left subtree
                    appendSubtree(left, left_loc + 1);

                    //move through left subtree
                    left_loc = left.stepThrough(left_loc + 1) - 1;
                    break;

                //copy right subtree into result
                case CB:
                    //turn on tree bit
                    b->bzet[loc] |= 0x08 >> ((NODE_ELS - 1) - i);

                    //append right subtree
                    appendSubtree(right, right_loc + 1);

                    //move through right subtree
                    right_loc = right.stepThrough(right_loc + 1) - 1;
                    break;

                //copy left subtree into result and negate
                case NA:
                    //turn on tree bit
                    b->bzet[loc] |= 0x08 >> ((NODE_ELS - 1) - i);

                    end = b->size;

                    //append left subtree
                    appendSubtree(left, left_loc + 1);

                    //negate appended subtree
                    subtreeNot(end, left.depthAt(left_loc + 1));

                    //move through left subtree
                    left_loc = left.stepThrough(left_loc + 1) - 1;
                    break;

                //copy right subtree into result and negate
                case NB:
                    //turn on tree bit
                    b->bzet[loc] |= 0x08 >> ((NODE_ELS - 1) - i);

                    end = b->size;

                    //append right subtree
                    appendSubtree(right, right_loc + 1);

                    //negate appended subtree
                    subtreeNot(end, right.depthAt(right_loc + 1));

                    //move through right subtree
                    right_loc = right.stepThrough(right_loc + 1) - 1;
                    break;

                default:
                    //should be impossible to get here, but just in case...
                    display_error("Bzet4::_binop: Something went terribly, terribly wrong", true);
                    break;
            }
        } 
        //only data bits
        else if (!cur_left_tree_bit && !cur_right_tree_bit) {
            if (lev != 1) {
                //process data bits
                int result_bit = do_data_op(op, cur_left_data_bit, cur_right_data_bit);

                //if result_bit is on, turn on the corresponding data bit in the result node
                if (result_bit) {
                    b->bzet[loc] |= 0x10 << i;
                }
            }
        }
        //malformed bzet
        else {
            display_error("Bzet4::_binop: Malformed bzet (data bit and tree bit are on simultaneously)", true);
        }
    }

    //set step
    if (b->size - loc > 255)
        b->step[loc] = 0;
    else 
        b->step[loc] = (unsigned char) (b->size - loc);

    //resulting node is empty
    if (b->bzet[loc] == 0x00) {
        //build bzet manually if bzet is empty
        if (b->size == 2) {
            clear();
            return NORMAL;
        }

        resize(b->size - 1);
        return EMPTY;
    }
    //resulting node is saturated
    else if (b->bzet[loc] == 0xF0) {
        if (b->size == 2) {
            resize(2);
            b->bzet[0]++;
            b->bzet[1] = 0x80;
            b->step[1] = 1;
            return NORMAL;
        }

        resize(b->size - 1);
        return SATURATED;
    }

#if (defined _DEBUG || defined DEBUG)
    assert(dust(b->bzet[loc]) == b->bzet[loc]);
#endif

    return NORMAL;
}
    
 /*****************************************************************************
 * 
 *  Function name: binop
 *
 *  Purpose:       Does a binary operation specified by OP using left and right
 *                 as operands
 *
 *  Inputs:        Bzet4 left: left size of operation
 *                 Bzet4 right: right hand side of operation
 *                 OP op: operation to do
 *
 *  Return values: A Bzet4 holding the result of the operation
 *
 *  Notes:         Implementation delegated to _binop
 * 
 *  Author:        Alex Chow
 *  Date:          12/10/2011
 *
 ****************************************************************************
Bzet4 Bzet4::binop(Bzet4& left, Bzet4& right, OP op) {
    //align bzets
    align(left, right);

    Bzet4 result;

    //set level byte
    result.b->bzet[0] = left.b->bzet[0];

    //operate
    result._binop(left, right, op, left.b->bzet[0]);

#if (defined _DEBUG || defined DEBUG)
    result.validateBzet();
#endif

    result.normalize();

#if (defined _DEBUG || defined DEBUG)
    result.validateBzet();
#endif

    return result;
}

/*****************************************************************************
 * 
 *  Function name: operator|
 *
 *  Purpose:       Bitwise OR operator
 *
 *  Inputs:        Bzet4& right: Right hand side of OR operation
 *  Return values: Bitwise OR of this Bzet and right
 * 
 *  Author:        Alex Chow
 *  Date:          4/10/2012
 *
 ****************************************************************************
Bzet4 Bzet4::operator|(const Bzet4& right) const {
    //this Bzet is empty
    if (empty()) {
        return right;
    }
    
    //the other Bzet is empty
    if (right.empty()) {
        return *this;
    }

    return binop(Bzet4(*this), Bzet4(right), OR);
}

/*****************************************************************************
 * 
 *  Function name: operator&
 *
 *  Purpose:       Bitwise AND operator
 *
 *  Inputs:        Bzet4& right: Right hand side of AND operation
 *  Return values: Bitwise AND of this Bzet and right
 * 
 *  Author:        Alex Chow
 *  Date:          4/10/2012
 *
 ****************************************************************************
Bzet4 Bzet4::operator&(const Bzet4& right) const {
    //this Bzet is empty
    if (empty()) {
        return *this;
    }
    
    //the other Bzet is empty
    if (right.empty()) {
        return right;
    }

    return binop(Bzet4(*this), Bzet4(right), AND);
}

/*****************************************************************************
 * 
 *  Function name: operator^
 *
 *  Purpose:       Bitwise XOR operator
 *
 *  Inputs:        Bzet4& right: Right hand side of XOR operation
 *  Return values: Bitwise XOR of this Bzet and right
 * 
 *  Author:        Alex Chow
 *  Date:          4/10/2012
 *
 ****************************************************************************
Bzet4 Bzet4::operator^(const Bzet4& right) const {
    //if the left Bzet is empty
    if (empty()) {
        return right;
    }

    //if the right Bzet is empty
    if (right.empty()) {
        return *this;
    }

    return binop(Bzet4(*this), Bzet4(right), XOR);
}

/*****************************************************************************
 * 
 *  Function name: operator==
 *
 *  Purpose:       Test equality between two Bzets
 *
 *  Inputs:        Bzet4 right: Right hand side of == operation
 *  Return values: true if this and right have the same bzet representation
 *                 false otherwise
 * 
 *  Author:        Alex Chow
 *  Date:          11/06/2011
 *
 ****************************************************************************
bool Bzet4::operator==(const Bzet4& right) const {
    //if sizes are the same, check content
    if (b->size == right.b->size) {
        if (!memcmp(b->bzet, right.b->bzet, b->size))
            return true;
    }

    //otherwise not equal
    return false;
}


/*****************************************************************************
 * 
 *  Function name: firstBit
 *
 *  Purpose:       Returns the location of the first bit on
 *
 *  Inputs:        None
 *  Return values: Location of the first bit on (e.g. bit 0)
 *                 -1 if no bits are on
 *
 *  Author:        Alex Chow
 *  Date:          10/28/2011
 *
 ****************************************************************************
int64_t Bzet4::firstBit() const {
    //empty Bzet
    if (b->size == 1)
        return -1;

    int level = b->bzet[0]; //current level in tree
    int64_t bitBase = 0; //number of bits prior to the current node
    size_t loc = 1; //location in Bzet

    //just traverse the leftmost branch
    while (level) {
        //special calculation at level 1
        if (level == 1) {
            //check data nodes 0 and 1
            unsigned char c = b->bzet[loc];

            //if any bit is set in node 0 and 1, first bit found
            if (c)
                for (int i = 7; i >= 0; --i)
                    if ((c >> i) & 1)
                        return bitBase + 7 - (int64_t) i;

            //if no data bit set in first 2 nodes, check next 2 nodes
            c = b->bzet[loc + 1];
            for (int i = 7; i >= 0; --i)
                if ((c >> i) & 1)
                    return bitBase + 15 - i;
#ifdef DEBUG
            assert(c);
#endif
        }

        //calculation for all other nodes
        unsigned char data_bits = (b->bzet[loc] >> NODE_ELS) & 0xF;
        unsigned char tree_bits = b->bzet[loc] & 0xF;

        //check each data bit then the corresponding tree bit
        for (int i = NODE_ELS - 1; i >= 0; --i) {
            //check data bit
            //if on, return the lowest bit
            if ((data_bits >> i) & 1)
                return bitBase + (NODE_ELS - 1 - i) * pow4(level);

            //check tree bit
            if ((tree_bits >> i) & 1) {
                //if on, advance bitBase to the correct prefix and check that node
                bitBase += (NODE_ELS - 1 - i) * pow4(level);
                break;
            }
        }

        --level;
        ++loc;
    }

    //it shouldn't get here unless the Bzet was malformed
    display_error("Bzet4::firstBit(): Warning: Bzet is malformed", true);
    return -1;
}

/*****************************************************************************
 * 
 *  Function name: lastBit
 *
 *  Purpose:       Returns the location of the last bit on
 *
 *  Inputs:        None
 *  Return values: Location of the last bit on (e.g. bit 0)
 *                 -1 if no bits are on
 *
 *  Author:        Alex Chow
 *  Date:          10/28/2011
 *
 ****************************************************************************
int64_t Bzet4::lastBit() const {
    //empty Bzet
    if (size() == 1)
        return -1;

    size_t loc = 1; //location in b->bzet (node)
    int level = b->bzet[0]; //current level in tree, initialize to bzet depth
    int64_t bitBase = 0; //number of bits prior to the current node

    do {
        //special calculation at level 1
        if (level == 1) {
            //check data nodes 2 and 3
            unsigned char c = b->bzet[loc + 1];

            if (c)
                for (int i = 0; i < 8; ++i)
                    if ((c >> i) & 1)
                        return bitBase + 15 - i;

            //if no data bit set in last 2 nodes, check first 2
            c = b->bzet[loc];
            for (int i = 0; i < 8; ++i)
                if ((c >> i) & 1)
                    return bitBase + 7 - i;
#if (defined _DEBUG || defined DEBUG)
            assert(c);
#endif
        }

        unsigned char c = b->bzet[loc];
        unsigned char data_bits = (c >> NODE_ELS) & 0xF;
        unsigned char tree_bits = c & 0xF;

        //get last tree and data bit
        int last_tree_bit = 0;
        int last_data_bit = 0;
        for (int i = NODE_ELS - 1; i >= 0; --i) {
            if ((tree_bits >> i) & 1)
                last_tree_bit = NODE_ELS - 1 - i;
            if ((data_bits >> i) & 1)
                last_data_bit = NODE_ELS - 1 - i;
        }

        //if this node has the farthest data bit set
        //if so, last bit found
        if (last_data_bit > last_tree_bit) {
            return bitBase + (last_data_bit + 1) * pow4(level) - 1;
        }

#if (defined _DEBUG || defined DEBUG)
        assert(last_data_bit != last_tree_bit);
#endif

        //if no tree bits are on, we're at the end
        if (!tree_bits) {
            //get the location of the last data bit
            int last_data_bit = 0;
            for (int i = NODE_ELS - 1; i >= 0; --i)
                if ((data_bits >> i) & 1)
                    last_data_bit = NODE_ELS - 1 - i;

            return bitBase + (last_data_bit + 1) * pow4(level) - 1;
        }

        //otherwise step through the nodes until you get to the node
        //corresponding to the last tree bit on in the current node
        ++loc;
        for (int i = NODE_ELS - 1; i > NODE_ELS - 1 - last_tree_bit; --i) 
            if ((tree_bits >> i) & 1) 
                loc = stepThrough(loc);

        bitBase += last_tree_bit * pow4(level);

        //loc now points to location of the node corresponding to the last
        //tree bit, so repeat this process
        --level; //moved down a level
    } while (loc < b->size);

    //if it gets here, then the Bzet was malformed...
    display_error("Bzet4::lastBit(): Warning: Bzet is malformed");
    return -1;
}

/*****************************************************************************
 * 
 *  Function name: count
 *
 *  Purpose:       Returns the number of bits on in the Bzet
 *
 *  Inputs:        None
 *  Return values: The number of bits on in the Bzet
 * 
 *  Author:        Alex Chow
 *  Date:          10/29/2011
 *
 ****************************************************************************
int64_t Bzet4::count() const {
    //empty bzet has 0
    if (empty())
        return 0;

    int64_t bitCount = 0;

    //go through each node
    for (size_t i = HEADER_SIZE; i < b->size; ++i) {
        int lev = depthAt(i);

        //if lev is 1, all data bits have weight 1, do bit count on all 8 bits
        if (lev == 1) {
            unsigned char c = b->bzet[i];
            int nub->bits = 0;
            for (int j = 0; j < 8; ++j) 
                nub->bits += ((c >> j) & 1);
            bitCount += nub->bits;
        }
        //otherwise do weighted bit count on data bits only
        else {
            unsigned char data_bits = (b->bzet[i] >> NODE_ELS) & 0xF;

            //count the number of data bits on
            int nub->bits = 0;
            for (int j = 0; j < NODE_ELS; ++j) 
                nub->bits += (data_bits >> j) & 1;

            //add weighted bit count to running total
            bitCount += nub->bits * pow4(lev);
        }
    }

    return bitCount;
}

/*****************************************************************************
 * 
 *  Function name: getBits
 *
 *  Purpose:       Fills an array pointed to by bits with the locations of
 *                 each bit on in the Bzet
 *
 *  Inputs:        int64_t* bits: pointer to an array of bits
 *                 int64_t limit: max number of bits to record
 *                                default 0 means no limit
 *                 int64_t start: look for bits starting from this bit
 *                                default 0
 *  Return values: Number of elements of bits filled in
 *                 -1 if bzet is empty
 * 
 *  Author:        Alex Chow
 *  Date:          10/31/2011
 *
 ****************************************************************************
int64_t Bzet4::getBits(int64_t* bits, int64_t limit, int64_t start) {
    if (empty())
        return -1;

    //counter for number of elements of bits filled in
    int64_t bitcount = 0;

    //set limit to number of bits on if unspecified
    if (!limit)
        limit = count();

    //create a disposable copy
    Bzet4 copy(*this);

    //if start is defined, move to the first bit that is after or equal to start
    if (start) {
        int64_t bit = copy.firstBit();
        while (bit < start && bit != -1) {
            copy.unset(bit);
            bit = copy.firstBit();
        }
    }
    
    //keep getting the location of the first bit
    //unset it after to move location of the first bit
    for (int64_t i = 0; i < limit; ++i) {
        bits[i] = copy.firstBit();
        copy.unset((int64_t) bits[i]);
        ++bitcount;
    }

    //return number of bit indices filled in
    return (limit < bitcount) ? limit : bitcount;
}

/*****************************************************************************
 * 
 *  Function name: empty
 *
 *  Purpose:       Returns whether or not the Bzet is empty
 *
 *  Inputs:        None
 *  Return values: true  if Bzet is empty (no bits are on)
 *                 false if Bzet is non-empty (some bit is on)
 * 
 *  Author:        Alex Chow
 *  Date:          10/28/2011
 *
 ****************************************************************************
bool Bzet4::empty() const {
    return (b->size == 1);
}

/*****************************************************************************
 * 
 *  Function name: at
 *
 *  Purpose:       Returns the value at the specified bit
 *
 *  Inputs:        int64_t bit: Location of target bit in the Bzet
 *  Return values: 0 if specified bit is out of range
 *                 Otherwise the value (1 or 0) of the specified bit
 * 
 *  Author:        Alex Chow
 *  Date:          10/23/2011
 *
 ****************************************************************************
bool Bzet4::at(int64_t bit) const {
    //get number of bits this Bzet stores
    int64_t totalBits = pow4(b->bzet[0] + 1);

    //if bit is out of range
    if (bit >= totalBits)
        return 0;

    return !((*this) & Bzet4(bit)).empty();
}

/*****************************************************************************
 * 
 *  Function name: setRange
 *
 *  Purpose:       Sets a range of bits from start to (start + len - 1)
 *
 *  Inputs:        int64_t start: First bit to set
 *                 int64_t len: Number of bits to set, starting from start
 *  Return values: None
 * 
 *  Author:        Alex Chow
 *  Date:          11/3/2011
 *
 ****************************************************************************
void Bzet4::setRange(int64_t start, int64_t len) {
    //validate parameters
    if (start < 0 || len <= 0) {
        display_error("Bzet4::setRange: invalid start or len, doing nothing");
        return;
    }

    //set each bit individually
    for (int64_t i = start; i < start + len; ++i)
        set(i);
}

/*****************************************************************************
 * 
 *  Function name: set
 *
 *  Purpose:       Sets specified bit within the Bzet to 1
 *
 *  Inputs:        int64_t bit: Location of target bit in the Bzet
 *  Return values: None
 * 
 *  Author:        Alex Chow
 *  Date:          10/23/2011
 *
 ****************************************************************************
void Bzet4::set(int64_t bit) {
    //verify valid bit to set
    if (bit < 0) {
        display_error("Bzet4::set: invalid bit to set, doing nothing");
        return;
    }

    //OR the bit in
    *this = *this | Bzet4(bit);
}

/*****************************************************************************
 * 
 *  Function name: unset
 *
 *  Purpose:       Sets specified bit within the Bzet to 0
 *
 *  Inputs:        int64_t bit: Location of target bit in the Bzet
 *  Return values: None
 * 
 *  Author:        Alex Chow
 *  Date:          10/25/2011
 *
 ****************************************************************************
void Bzet4::unset(int64_t bit) {
    //verify valid bit to unset
    if (bit < 0) {
        display_error("Bzet4::unset: invalid bit to unset, doing nothing");
        return;
    }

    Bzet4 b(bit);
    align(*this, b);
    *this = *this & ~b;
}

/*****************************************************************************
 * 
 *  Function name: depth
 *
 *  Purpose:       Returns depth of Bzet
 *
 *  Inputs:        None
 *  Return values: Depth of Bzet
 * 
 *  Author:        Alex Chow
 *  Date:          10/24/2011
 *
 ****************************************************************************
int Bzet4::depth() const {
    return b->bzet[0];
}

/*****************************************************************************
 * 
 *  Function name: size
 *
 *  Purpose:       Returns size of Bzet in bytes
 *
 *  Inputs:        None
 *  Return values: Size of Bzet in bytes
 * 
 *  Author:        Alex Chow
 *  Date:          10/23/2011
 *
 ****************************************************************************
size_t Bzet4::size() const {
    return b->size;
}

/*****************************************************************************
 * 
 *  Function name: clear
 *
 *  Purpose:       Resets Bzet4 to a Bzet of length 0
 *
 *  Inputs:        None
 *  Return values: None
 * 
 *  Author:        Alex Chow
 *  Date:          10/23/2011
 *
 ****************************************************************************
void Bzet4::clear() {
    resize(1);
    b->bzet[0] = 0x00; //set level byte to 0
}

/*****************************************************************************
 * 
 *  Function name: hex
 *
 *  Purpose:       Copies the hex representation of Bzet to str
 *
 *  Inputs:        void* str: Copy target
 *  Return values: None
 * 
 *  Notes:         Allocation and deallocation of buffers is handled by caller
 *                 Use as follows
 *                     Bzet4 b();
 *                     void* buf = malloc(b.size());
 *                     b.hex(buf);
 *                     free(buf);
 *
 *  Author:        Alex Chow
 *  Date:          10/23/2011
 *
 ****************************************************************************
void Bzet4::hex(void* str) const {
    memcpy(str, b->bzet, b->size); 
}

/*****************************************************************************
 * 
 *  Function name: printBzet
 *
 *  Purpose:       Prints human-readable representation of the Bzet to
 *                 target file descriptor
 *                 Implementation delegated to _printBzet
 *
 *  Inputs:        int stdOffset: user can specify a standard offset of all lines
 *                                default is 0
 *                 FILE* target: file descriptor to print to
 *                               default is stdout
 *  Return values: None
 * 
 *  Author:        Alex Chow
 *  Date:          10/26/2011
 *
 ****************************************************************************
void Bzet4::printBzet(int stdOffset, FILE* target) const {
    _printBzet(stdOffset, target);
}

/*****************************************************************************
 * 
 *  Function name: pow4
 *
 *  Purpose:       Returns 4**x
 *
 *  Inputs:        int x: Exponent for 4**x
 *  Return values: 4**x
 * 
 *  Author:        Alex Chow
 *  Date:          10/26/2011
 *
 ****************************************************************************
size_t Bzet4::pow4(int x) {
    //reference table for powers of less than 10
    if (x < 10)
        return powersof4[x];

    //otherwise calculate power
    return (size_t) pow((double) 4, x);
}

/*****************************************************************************
 * 
 *  Function name: buildDepth
 *
 *  Purpose:       Returns the depth of the Bzet tree needed to hold bit
 *
 *  Inputs:        int64_t bit: The bit that the Bzet needs to hold
 *  Return values: Tree depth needed to hold bit
 * 
 *  Author:        Alex Chow
 *  Date:          10/23/2011
 *
 ****************************************************************************
int Bzet4::buildDepth(int64_t bit) {
    //invalid bit to build depth
    if (bit < 0)
        return -1;

    //if bit index < 16, then we're at level 1
    int newDepth = 1;

    //otherwise, for each power of 4 we go up one level
    while (bit > 15) {
        bit /= (int64_t) NODE_ELS;
        ++newDepth;
    }

    return newDepth;
}

/*****************************************************************************
 * 
 *  Function name: dust
 *
 *  Purpose:       Turns off a tree bit if its corresponding data bit is on
 *
 *  Inputs:        unsigned char x: byte to dust
 *  Return values: The dusted byte
 * 
 *  Author:        Alex Chow
 *  Date:          10/25/2011
 *
 ****************************************************************************
unsigned char Bzet4::dust(unsigned char x) {
    return ~(x >> NODE_ELS) & x;
}

/*****************************************************************************
 * 
 *  Function name: align
 *
 *  Purpose:       Aligns the two Bzets by inserting 0x08 nodes to extend the 
 *                 shorter Bzet to the same depth
 *
 *  Inputs:        Bzet4& b1: First Bzet
 *                 Bzet4& b2: Second Bzet
 *  Return values: None
 * 
 *  Author:        Alex Chow
 *  Date:          10/25/2011
 *
 ****************************************************************************
void Bzet4::align(Bzet4& b1, Bzet4& b2) {
    int diffDepth = b1.depth() - b2.depth();

    //already aligned, do nothing
    if (diffDepth == 0)
        return;

    //b1 is bigger
    if (diffDepth > 0) {
        //resize b2
        b2.resize(b2.size() + diffDepth);

        //modify b2 depth
        b2.b->bzet[0] = b2.depth() + diffDepth;

        //shift Bzet to make room for extra nodes
        for (size_t i = b2.size() - 1; i >= HEADER_SIZE + diffDepth; --i) {
            b2.b->bzet[i] = b2.b->bzet[i - diffDepth];
            b2.b->step[i] = b2.b->step[i - diffDepth];
        }

        //fill in extra nodes
        for (int i = 0; i < diffDepth; ++i) {
            b2.b->bzet[HEADER_SIZE + i] = 0x08;

            //set b->step
            if ((size_t) b2.b->step[diffDepth + 1] + (diffDepth - i) >= 255)
                b2.b->step[HEADER_SIZE + i] = 0;
            else
                b2.b->step[HEADER_SIZE + i] = b2.b->step[diffDepth + 1] + (diffDepth - i);
        }
    } 
    //b2 is bigger
    else {
        //change sign of diffDepth
        diffDepth = -diffDepth;

        //resize b1
        b1.resize(b1.size() + diffDepth);

        //modify b1 depth
        b1.b->bzet[0] = b1.depth() + diffDepth;

        //shift Bzet to make room for extra nodes
        for (size_t i = b1.size() - 1; i >= HEADER_SIZE + diffDepth; --i) {
            b1.b->bzet[i] = b1.b->bzet[i - diffDepth];
            b1.b->step[i] = b1.b->step[i - diffDepth];
        }

        //fill in extra nodes
        for (int i = 0; i < diffDepth; ++i) {
            b1.b->bzet[HEADER_SIZE + i] = 0x08;

            //set b->step
            if ((size_t) b1.b->step[diffDepth + 1] + (diffDepth - i) >= 255)
                b1.b->step[HEADER_SIZE + i] = 0;
            else
                b1.b->step[HEADER_SIZE + i] = b1.b->step[diffDepth + 1] + (diffDepth - i);
        }
    }

#if (defined _DEBUG || defined DEBUG)
    b1.validateBzet();
    b2.validateBzet();
#endif
}

/*****************************************************************************
 * 
 *  Function name: depthAt
 *
 *  Purpose:       Returns what depth the byte (node) at loc is
 *
 *  Inputs:        size_t tLoc: target byte
 *  Return values: -1 if tLoc or curLoc are out of range
 *                 otherwise depth of node at loc
 * 
 *  Author:        Alex Chow
 *  Date:          10/28/2011
 *
 ****************************************************************************
int Bzet4::depthAt(size_t tLoc) const {
    //check if out of range
    if (tLoc > size())
        return -1;

    //if loc is 1, just return bzet depth
    if (tLoc == 1)
        return depth();

    //the node at b->bzet[1] is at depth b->bzet[0]
    int cdepth = depth();
    size_t nextLoc = 1;

    while (true) {
        //found target
        if (nextLoc == tLoc)
            return cdepth;

        //special handling for level 1 data
        if (cdepth == 1 && (int) nextLoc - (int) tLoc <= 1)
            return 1;

        size_t next = stepThrough(nextLoc);
        //check leftmost nodes
        if (next > tLoc) {
            nextLoc++;
            cdepth--;
        } 
        else {
            //otherwise move to the next node horizontally
            nextLoc = next;
        }
    }
  
    return cdepth;
}

/*****************************************************************************
 * 
 *  Function name: stepThrough
 *
 *  Purpose:       Returns the location of the next node after walking through
 *                 the subtree starting at loc
 *
 *  Inputs:        size_t loc: location of starting node of target subtree
 *  Return values: -1 if loc is out of valid range (1 to b->size - 1)
 *                 otherwise the location of the next node after walking
 *                 through the subtree starting at loc
 * 
 *  Author:        Alex Chow
 *  Date:          10/28/2011
 *
 ****************************************************************************
size_t Bzet4::stepThrough(size_t loc) const {
    //make sure loc is in range
    if (loc <= 0 || loc >= b->size)
        return -1;

    //retrieve the offset of the end of the subtree at loc
    unsigned char loc_offset = b->step[loc];

    //if the offset is 0, the offset is too large to be actually stored
    //compute it by using the subtrees stemming from this node
    if (loc_offset == 0) {
        //get node tree bits
        unsigned char tree_bits = b->bzet[loc] & 0xF;

        //advance to location of first subtree node
        loc++;

        //step through each subtree
        for (int i = NODE_ELS - 1; i >= 0; i--) {
            if ((tree_bits >> i) & 1)
                loc = stepThrough(loc);
        }

        return loc;
    }

    return loc + b->step[loc];
}

/*****************************************************************************
 * 
 *  Function name: subtreeNot
 *
 *  Purpose:       Does a bitwise not of the subtree starting at loc
 *
 *  Inputs:        size_t loc: location of starting node of target subtree
 *                 int depth: depth of the node at loc
 *  Return values: None
 * 
 *  Author:        Alex Chow
 *  Date:          11/21/2011
 *
 ****************************************************************************
void Bzet4::subtreeNot(size_t loc, int depth) {
    if (loc == 0 || loc >= b->size)
        return;

    if (!depth)
        depth = depthAt(loc);

    //if no tree bits are on, just not the data bits
    if (!(b->bzet[loc] & 0xF) && depth > 1) {
        b->bzet[loc] = ~(b->bzet[loc] | 0xF);
        return;
    }

    //if depth is 1, we just not this and the next node
    if (depth == 1) {
        b->bzet[loc] = ~b->bzet[loc];
        b->bzet[loc + 1] = ~b->bzet[loc + 1];
        return;
    }

    size_t nextLoc = loc + 1;

    //decrement depth
    --depth;

    for (int i = NODE_ELS - 1; i >= 0; --i) {
        //if tree bit is on, recursively not each node
        if ((b->bzet[loc] >> i) & 1) {
            subtreeNot(nextLoc, depth);
            nextLoc = stepThrough(nextLoc);
        }
    }

    //get the current byte
    unsigned char c = b->bzet[loc]; 

    //new data bits
    unsigned char data_bits = (~(c >> NODE_ELS)) << NODE_ELS;

    //add tree bits
    unsigned char newc = data_bits | (c & 0xF);

    //set up for dusting with tree bits priority
    unsigned char tree_bits = c & 0xF;

    //dust using tree bits and write back
    b->bzet[loc] = newc & (~(tree_bits << NODE_ELS) | tree_bits);
}

/*****************************************************************************
 * 
 *  Function name: _printBzet
 *
 *  Purpose:       Prints human-readable representation of the Bzet to
 *                 target file descriptor
 *
 *  Inputs:        int stdOffset: user can specify a standard offset of all lines
 *                                default is 0
 *                 FILE* target: file descriptor to print to
 *                               default is stdout
 *                 size_t loc: location within the Bzet to print, used during recursion
 *                          default is 1 (first node)
 *                 int offset: offset in tree levels, used during recursion
 *                             default is 0
 *  Return values: None
 *
 *  Notes:         Public interface through printBzet(int stdOffset, FILE* target)
 * 
 *  Author:        Alex Chow
 *  Date:          10/26/2011
 *
 ****************************************************************************
void Bzet4::_printBzet(int stdOffset, FILE* target, size_t loc, int depth, int offset, bool pad) const {
    //print level info
    if (loc == 1) {
        fprintf(target, "%.2XL", b->bzet[0]);

        //quick fix for empty bzet
        if (empty())
            fprintf(target, "\n");
    }

    if (!depth)
        depth = depthAt(loc);

    //no reading past the bzet!
    if (loc >= size())
        return;

    unsigned char c = b->bzet[loc];
    unsigned char data_bits = (c >> NODE_ELS) & 0xF;
    unsigned char tree_bits = c & 0xF;

    //print offset if any
    if (pad) {
        int blanks = 5*offset + 3 + stdOffset;
        for (int j = 0; j < blanks; j++)
            fprintf(target, " ");
    }

    //if depth is 1, print special level 1 interpretation
    if (depth == 1) {   
        fprintf(target, "D(%.2X%.2X)\n", b->bzet[loc], b->bzet[loc + 1]);
    } 
    //if no tree bits on, then this is a data node
    else if ((c & 0xF0) == c) {
        fprintf(target, "D(%X)\n", data_bits);
    }
    else {     
        //print the current node
        fprintf(target, "[%X-%X]", data_bits, tree_bits);

        //recursively print child nodes
        bool firstNode = true;
        --depth;
        for (int i = NODE_ELS - 1; i >= 0; --i) {
            //if tree bit set
            if ((tree_bits >> i) & 0x1) {
                //print first node without offset
                if (firstNode) {
                    _printBzet(stdOffset, target, loc + 1, depth, offset + 1);
                    firstNode = false;
                } else {
                    _printBzet(stdOffset, target, stepThrough(loc + 1), depth, offset + 1, true);
                    loc = stepThrough(loc + 1) - 1;
                }
            }
        }
    }
}

/*****************************************************************************
 * 
 *  Function name: resize
 *
 *  Purpose:       Resizes b->bzet to the specified size in bytes
 *
 *  Inputs:        int size: New size of b->bzet
 *  Return values: None
 * 
 *  Author:        Alex Chow
 *  Date:          10/23/2011
 *
 ****************************************************************************
void Bzet4::resize(size_t size) {
    //if reallocation is required
    if (size > b->bufsize) {
        //get new size required by growing it by RESIZE_SCALE repeatedly
        while (size > b->bufsize)
            b->bufsize *= RESIZE_SCALE;

        //reallocate
        b->bzet = (unsigned char*) realloc(b->bzet, b->bufsize * sizeof(unsigned char));
        b->step = (unsigned char*) realloc(b->step, b->bufsize * sizeof(unsigned char));

        //check that realloc succeeded
        if (!b->bzet || !b->step) {
            fprintf(stderr, "Fatal error: Resizing bzet failed attempting to allocate %d bytes\n", (int) size);
            display_error("", true);
        }
    }

    //update internal size
    b->size = size;
}

/*****************************************************************************
 * 
 *  Function name: normalize
 *
 *  Purpose:       Normalizes b->bzet to get rid of extra nodes
 *
 *  Inputs:        None
 *  Return values: None
 * 
 *  Author:        Alex Chow
 *  Date:          10/26/2011
 *
 ****************************************************************************
void Bzet4::normalize() {
    //no need to normalize bzet with only 1 level
    if (b->bzet[0] == 0x01)
        return;

    //count leading 0x08 nodes
    int i = 1;
    while (b->bzet[i] == 0x08)
        ++i;
    --i; //i is now the number of leading 0x08 nodes

    //if there are leading 0x08 nodes, get rid of them
    if (i) {
        //drop i nodes from the beginning
        dropNodes(1, i);

        //change depth accordingly
        b->bzet[0] -= i;
    }
}

/*****************************************************************************
 * 
 *  Function name: loadBzet
 *
 *  Purpose:       Overwrites the current Bzet with a literal Bzet (byte form)
 *
 *  Inputs:        void* bzet_literal: Pointer to a Bzet literal
 *                 int size: Size of the Bzet literal
 *  Return values: None
 * 
 *  Author:        Alex Chow
 *  Date:          11/4/2011
 *
 ****************************************************************************
void Bzet4::loadBzet(void* bzet_literal, int size) {
    if (bzet_literal) {
        resize(size);
        memcpy(b->bzet, bzet_literal, size);
    } else {
        display_error("Bzet4::loadBzet(void*, int64_t): null pointer, doing nothing");
    }
}

*/
#endif
