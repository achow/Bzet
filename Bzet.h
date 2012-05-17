/*****************************************************************************
 * Generic code for Bzet8/Bzet16/Bzet32/Bzet64
 *
 * Things to do:
 *  - Good error messages, better error handling in general (?)
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

//#define USE_LITERAL

#define PASTE_(x,y) x ## y
#define PASTE(x,y) PASTE_(x,y)

#define PASTE_UNDER_(x,y) x ## _ ## y
#define PASTE_UNDER(x,y) PASTE_UNDER_(x,y)

#define STRINGIFY_(x) #x
#define STRINGIFY(x) STRINGIFY_(x)

#ifndef BZET
#define BZET PASTE(Bzet, NODE_ELS)
#endif

#define POW PASTE(pow, NODE_ELS)

#if NODE_ELS == 64
typedef uint64_t halfnode_t;
#define NPOWERS 6
size_t const PASTE(powersof, NODE_ELS)[NPOWERS] = 
    { 1, 64, 4096, 262144, 16777216, 1073741824 };
#elif NODE_ELS == 32
typedef uint32_t halfnode_t;
#define NPOWERS 7
size_t const PASTE(powersof, NODE_ELS)[NPOWERS] = 
    { 1, 32, 1024, 32768, 1048576, 33554432, 1073741824 };
#elif NODE_ELS == 16
typedef uint16_t halfnode_t;
#define NPOWERS 8
size_t const PASTE(powersof, NODE_ELS)[NPOWERS] = 
    { 1, 16, 256, 4096, 65536, 1048576, 16777216, 268435456 };
#elif NODE_ELS == 8
typedef uint8_t halfnode_t;
#define NPOWERS 11
size_t const PASTE(powersof, NODE_ELS)[NPOWERS] = 
    { 1, 8, 64, 512, 4096, 32768, 262144, 2097152, 16777216, 134217728, 1073741824 };
#elif NODE_ELS == 4
typedef uint8_t node_t;
#define NPOWERS 17
size_t const PASTE(powersof, NODE_ELS)[NPOWERS] = 
    { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144, 1048576, 4194304, 16777216, 67108864, 268435456, 1073741824};
#else
#error "Invalid NODE_ELS provided"
#endif

#if STEP_BYTES == 2
typedef unsigned short step_t;
#else
#define STEP_BYTES 1
typedef unsigned char step_t;
#endif

#define STEP_T_MAX ((size_t) 1 << (STEP_BYTES * 8))

#define INITIAL_ALLOC 1024
#define RESIZE_SCALE 2

enum OP { 
    OP_0000, OP_0001, OP_0010, OP_0011, OP_0100, OP_0101, OP_0110, OP_0111,
    OP_1000, OP_1001, OP_1010, OP_1011, OP_1100, OP_1101, OP_1110, OP_1111,
    OP_AND = 1, OP_XOR = 6, OP_OR = 7, OP_NOR = 8, OP_NOT = 10, OP_NAND = 14 };
enum ACTION { DA0, DA1, DB0, DB1, CA, CB, NA, NB };

#ifdef USE_LITERAL
enum NODETYPE { SATURATED, EMPTY, NORMAL, LITERAL };
#else
enum NODETYPE { SATURATED, EMPTY, NORMAL };
#endif

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

// Bzet class
class BZET {
    public:
        BZET();
        BZET(int64_t bit);
        BZET(int64_t startbit, int64_t len);
        BZET(const BZET& b);
        ~BZET();

        BZET& operator=(const BZET& right);
        BZET operator~() const;
        BZET operator|(BZET& right);
		BZET operator&(BZET& right);
		BZET operator^(BZET& right);
		bool operator==(const BZET& right) const;

        static BZET binop(BZET& left, BZET& right, OP op);

        bool at(int64_t bit) const;
        void setRange(int64_t start, int64_t len);
        void set(int64_t bit);
        void seqset(int64_t bit);
        void unset(int64_t bit);

        int depth() const;
        size_t size() const;

        static void align(BZET& b1, BZET& b2);
        void normalize();

        void clear();

        void hex(void* str) const;
        void printBzet(int stdOffset = 0, FILE* target = stdout) const;
        void exportTo(FILE* stream) const;
        void importFrom(FILE* stream, size_t size);

		int64_t firstBit() const;
        int64_t lastBit() const;
        int64_t count() const;
        int64_t getBits(int64_t* bits, int64_t limit = 0, int64_t start = 0);
		bool empty() const;

        void printStep() {
#if NODE_ELS == 4
            for (size_t i = 0; i < m_nnodes; i++) {
#else
            for (size_t i = 0; i < m_nhalfnodes; i++) {
#endif
                if ((i % 10) == 0)
                    printf("\n");
                printf("0x%.*X ", sizeof(step_t) * 2, (step_t) m_step[i]);
                //printf("%d ", m_step[i]);
            }
            printf("\n");
        }

    private:
#ifdef USE_LITERAL
        //BZET_PTR bitstobzet(void *data, size_t len);
        //void treetobits(unsigned char *buf, halfnode_t *node, int depth);
#endif

        NODETYPE _binop(BZET& left, BZET& right, OP op, int lev, size_t left_loc = 0, size_t right_loc = 0);
        void _printBzet(int stdOffset, FILE* target, int depth, size_t loc = 0, int offset = 0, bool pad = 0) const;
        int64_t _count(size_t loc, int depth) const;
        void subtree_not(size_t loc, int depth);
        static NODETYPE _seqset(BZET& b, size_t locb, BZET& right, size_t locright, int depth);
        static bool _at(BZET& b, size_t locb, BZET& right, size_t locright, int depth);
        size_t build_step(size_t loc, int depth);

        size_t step_through(size_t loc, int depth) const;

        // Inline auxiliary functions
        static void display_error(const char* message, bool fatal = false, FILE* output = stderr);
        void init(size_t initial_alloc = INITIAL_ALLOC);
#if NODE_ELS == 4
        void resize(size_t nnodes);
#else
        void resize(size_t nhalfnodes);
#endif
        static size_t POW(int n);
        static int do_data_op(OP op, int left_data_bit, int right_data_bit);
        void append_subtree(BZET& src, size_t loc, int depth);
        void set_step(size_t loc, int depth);
        //void mod_step(size_t loc, int depth, int diff);

#if NODE_ELS == 4
        size_t m_nbufnodes;
        size_t m_nnodes;
        node_t* m_bzet;
#else
        size_t m_nbufhalfnodes;
        size_t m_nhalfnodes;
        halfnode_t* m_bzet; //points to the bzet
#endif
        step_t* m_step; //points to an array that holds step_through values
        unsigned char m_depth;
};

// Inline auxiliary functions

// Error message printing and optional exiting
inline
void BZET::display_error(const char* message, bool fatal, FILE* output) {
    fprintf(output, "%s\n", message);
    if (fatal)
        exit(1);
}

// Common Bzet constructor initialization
// TODO: Add error messages
inline 
void BZET::init(size_t initial_alloc) {
    // Allocate bzet node array
#if NODE_ELS == 4
    m_bzet = (node_t *) malloc(initial_alloc * sizeof(node_t));
#else
    m_bzet = (halfnode_t *) malloc(initial_alloc * sizeof(halfnode_t));
#endif
    if (!m_bzet) {
        display_error("init malloc failed\n", true);
    }

    // Allocate step array
    m_step = (step_t *) malloc(initial_alloc * sizeof(step_t));
    if (!m_step) {
        display_error("init malloc failed\n", true);
        free(m_bzet);
    }

    m_depth = 0;
#if NODE_ELS == 4
    m_nnodes = 0;
    m_nbufnodes = INITIAL_ALLOC;
#else
    m_nhalfnodes = 0;
    m_nbufhalfnodes = INITIAL_ALLOC;
#endif
}

// Resizes the buffers in the Bzet if necessary
inline
#if NODE_ELS == 4
void BZET::resize(size_t nnodes) {
    // If reallocation is required
    if (nnodes > m_nbufnodes) {
        // Get new size required by growing it by RESIZE_SCALE repeatedly
        while (nnodes > m_nbufnodes)
            m_nbufnodes *= RESIZE_SCALE;

        // Reallocate bzet
        node_t *bzet_temp = (node_t *) realloc(m_bzet, m_nbufnodes * sizeof(node_t));

        // Reallocate step
        step_t *step_temp = (step_t *) realloc(m_step, m_nbufnodes * sizeof(step_t));

        // Check that realloc succeeded
        // TODO: Replace with better checking above
        if (!m_bzet || !m_step) {
            fprintf(stderr, "Fatal error: Resizing bzet failed attempting to allocate %l bytes\n", m_nbufnodes * sizeof(node_t));
            display_error("", true);
        }

        m_bzet = bzet_temp;
        m_step = step_temp;
    }

    // Update internal halfnode counter
    m_nnodes = nnodes;
}
#else
void BZET::resize(size_t nhalfnodes) {
    // If reallocation is required
    if (nhalfnodes > m_nbufhalfnodes) {
        // Get new size required by growing it by RESIZE_SCALE repeatedly
        while (nhalfnodes > m_nbufhalfnodes)
            m_nbufhalfnodes *= RESIZE_SCALE;

        // Reallocate bzet
        halfnode_t *bzet_temp = (halfnode_t *) realloc(m_bzet, m_nbufhalfnodes * sizeof(halfnode_t));

        // Reallocate step
        step_t *step_temp = (step_t *) realloc(m_step, m_nbufhalfnodes * sizeof(step_t));

        // Check that realloc succeeded
        // TODO: Replace with better checking above
        if (!m_bzet || !m_step) {
            fprintf(stderr, "Fatal error: Resizing bzet failed attempting to allocate %l bytes\n", m_nbufhalfnodes * sizeof(halfnode_t));
            display_error("", true);
        }

        m_bzet = bzet_temp;
        m_step = step_temp;
    }

    // Update internal halfnode counter
    m_nhalfnodes = nhalfnodes;
}
#endif

// Trusty fast pow8/pow16/pow32 (assumes n >= 0)
inline
size_t BZET::POW(int n) {
    // Lookup power in table if possible
    if (n < NPOWERS)
        return PASTE(powersof, NODE_ELS)[n];

    return (size_t) pow((double) NODE_ELS, n);
}

// Does operations for two data bits, used in Bzet_binop
inline
int BZET::do_data_op(OP op, int left_data_bit, int right_data_bit) {
    // Use op directly to build result bit
    return (op >> (3 - (((int) left_data_bit << 1) | (int) right_data_bit))) & 0x1;
}

// Append a subtree from src to dest starting at loc in src
inline
void BZET::append_subtree(BZET& src, size_t loc, int depth) {
    // Calculate copy size and cache copy destination
    size_t copy_size = src.step_through(loc, depth) - loc;

#if NODE_ELS == 4
    size_t dst_loc = m_nnodes;

    // Resize to accomodate copy_size new elements
    resize(m_nnodes + copy_size);

    // Do copy
    memcpy(m_bzet + dst_loc, src.m_bzet + loc, copy_size * sizeof(node_t));
    memcpy(m_step + dst_loc, src.m_step + loc, copy_size * sizeof(step_t));
#else
    size_t dst_loc = m_nhalfnodes;

    // Resize to accomodate copy_size new elements
    resize(m_nhalfnodes + copy_size);

    // Do copy
    memcpy(m_bzet + dst_loc, src.m_bzet + loc, copy_size * sizeof(halfnode_t));
    memcpy(m_step + dst_loc, src.m_step + loc, copy_size * sizeof(step_t));
#endif
}

// Build step at loc
inline
void BZET::set_step(size_t loc, int depth) {
#if NODE_ELS == 4
    if (depth == 1) {
        m_step[loc] = 2;
        return;
    }

    size_t locdiff = m_nnodes - loc;

    if (locdiff >= STEP_T_MAX)
        m_step[loc] = 0;
    else
        m_step[loc] = (step_t) locdiff;
#else
    if (depth == 0) {
        m_step[loc] = 1;
        return;
    }

    size_t locdiff = m_nhalfnodes - loc;

    if (locdiff >= STEP_T_MAX) {
        if (locdiff <= STEP_T_MAX * STEP_T_MAX - 1) {
            // MSB first
            m_step[loc] = (step_t) (locdiff / STEP_T_MAX);
            m_step[loc + 1] = (step_t) (locdiff % STEP_T_MAX);
        }
        else {
            m_step[loc] = 0;
            m_step[loc + 1] = 0;
        }
    }
    else {
        m_step[loc] = 0;
        m_step[loc + 1] = (step_t) locdiff;
    }
#endif
}

/*inline
void BZET::mod_step(size_t loc, int depth, int diff) {
    if (depth == 0)
        return;

    size_t step = m_step[loc] * STEP_T_MAX + m_step[loc + 1];
    step -= diff;

    m_step[loc] = (step_t) (step / STEP_T_MAX);
    m_step[loc + 1] = (step_t) (step % STEP_T_MAX);
}*/

#ifdef BZET_IMPL_

// -- CONSTRUCTORS / DESTRUCTORS --

// Default constructor. Creates an empty bzet.
BZET::BZET() { 
    init();
}

// Creates a bzet with bit bit set
BZET::BZET(int64_t bit) {
    if (bit < 0) {
        display_error("bit < 0 in Bzet_new(bit)\n", true);
    }

    init();

#if NODE_ELS == 4
    // Build depth
    int depth = 1;

    int64_t tempbit = bit;
    // Level 1 holds 16 bits, numbers 0 to 15
    while (tempbit > 15) {
        tempbit /= (int64_t) NODE_ELS;
        ++depth;
    }

    // Resize to accomodate full bzet
    // +1 since level 1 requires 2 bytes to store data literals
    resize(depth + 1);

    //set depth bit
    m_depth = (unsigned char) depth;

    // Set data nodes at level 1
    // Set bits containing level 0 nodes 0 and 1 (bits 0-7)
    m_bzet[depth - 1] = (node_t) 0x80 >> (bit & 0xF);
    // Set bits containing level 0 nodes 2 and 3 (bits 8-15) if previous node not set
    m_bzet[depth] = (m_bzet[depth - 1]) ? 0x00 : (node_t) 0x80 >> (bit & 0x7);

    // Since level 1 set, divide by 16 to see if tree nodes needed
    bit >>= 4; 

    // Set tree nodes if needed
    int loc = depth - 1;
    while (bit > 0) {
        m_bzet[--loc] = (node_t) 0x08 >> (bit & 0x3);
        bit >>= 2;
    }

    // Set m_step
    // No need to check m_size <= 255 to make sure m_step values doesn't overflow
    // Since it would never happen (4^255 ~ 10^153)
    for (int i = 0; i < m_nnodes; ++i)
        m_step[i] = (step_t) (m_nnodes - i);
#else
    // Build depth
    int depth = 0;
    while (POW(depth + 1) <= (size_t) bit)
        depth++;

    // Resize to accomodate full bzet
    resize(2*depth + 1);

    //set depth byte
    m_depth = (unsigned char) depth;

    // Set level 0 data node
    // Each depth contains 2 halfnodes, b->bzet[2*depth] is the data node
    m_bzet[2*depth] = (halfnode_t) 0x1 << (NODE_ELS - 1 - (bit % NODE_ELS));

    // Zero out other nodes, since we will be only filling the tree portions next
    memset(m_bzet, 0x00, 2*depth * sizeof(halfnode_t));

    // Set tree bits, working linearly from the node at the head of the bzet
    for (int i = depth; i >= 1; i--) {
        // Get weight of bits at this level
        size_t cpow = POW(i);

        // Calculate tree bit to set
        int setbit = (int) (bit / cpow) % NODE_ELS;

        // Set the tree bit
        m_bzet[2*(depth - i) + 1] = (halfnode_t) 0x1 << (NODE_ELS - 1 - setbit);
    }

    // Set b->step
    // No need to check b->size <= 255 to make sure b->step values doesn't
    // overflow since it would never happen (4^255 ~ 10^153)
    for (int i = 0; i < m_nhalfnodes; i += 2) {
        m_step[i] = 0;
        m_step[i + 1] = (step_t) (m_nhalfnodes - i);
    }
    m_step[m_nhalfnodes - 1] = 1;
#endif
}

// Range constructor
BZET::BZET(int64_t startbit, int64_t len) {
    if (len <= 0 || startbit < 0)
        display_error("len <= 0 || startbit < 0", true);

    init();

    setRange(startbit, len);
}

// Copy constructor
BZET::BZET(const BZET& b) {
    // Copy contents over
    m_depth = b.m_depth;
#if NODE_ELS == 4
    m_nbufnodes = b.m_nbufnodes;
    m_nnodes = b.m_nnodes;

    // Deep copy of bzet and step
    m_bzet = (node_t *) malloc(m_nbufnodes * sizeof(node_t));
    if (!m_bzet) {
        display_error("copy malloc failed\n", true);
    }

    m_step = (step_t *) malloc(m_nbufnodes * sizeof(step_t));
    if (!m_step) {
        display_error("copy malloc failed\n", true);
        free(m_bzet);
    }

    memcpy(m_bzet, b.m_bzet, m_nnodes * sizeof(node_t));
    memcpy(m_step, b.m_step, m_nnodes * sizeof(step_t));
#else
    m_nbufhalfnodes = b.m_nbufhalfnodes;
    m_nhalfnodes = b.m_nhalfnodes;

    // Deep copy of bzet and step
    m_bzet = (halfnode_t *) malloc(m_nbufhalfnodes * sizeof(halfnode_t));
    if (!m_bzet) {
        display_error("copy malloc failed\n", true);
    }

    m_step = (step_t *) malloc(m_nbufhalfnodes * sizeof(step_t));
    if (!m_step) {
        display_error("copy malloc failed\n", true);
        free(m_bzet);
    }

    memcpy(m_bzet, b.m_bzet, m_nhalfnodes * sizeof(halfnode_t));
    memcpy(m_step, b.m_step, m_nhalfnodes * sizeof(step_t));
#endif
}

// Destructor
BZET::~BZET() {
    free(m_bzet);
    free(m_step);
}

// -- OPERATIONS ---

// Assignment operator
BZET& BZET::operator=(const BZET& right) {
#if NODE_ELS == 4
    // Resize left's bzet and step buffers
    resize(right.m_nnodes);

    // Copy contents of bzet over
    m_depth = right.m_depth;

    memcpy(m_bzet, right.m_bzet, m_nnodes * sizeof(node_t));
    memcpy(m_step, right.m_step, m_nnodes * sizeof(step_t));
#else
    // Resize left's bzet and step buffers
    resize(right.m_nhalfnodes);

    // Copy contents of bzet over
    m_depth = right.m_depth;

    memcpy(m_bzet, right.m_bzet, m_nhalfnodes * sizeof(halfnode_t));
    memcpy(m_step, right.m_step, m_nhalfnodes * sizeof(step_t));
#endif

    return *this;
}

// Bitwise NOT
BZET BZET::operator~() const {
    // Clone this bzet
    BZET copy(*this);

    // NOT it in place
    copy.subtree_not(0, copy.m_depth);

    return copy;
}

// Bitwise OR
BZET BZET::operator|(BZET& right) {
    // If left bzet is empty, the bitwise OR will be equal to the right bzet
    if (empty())
        return BZET(right);
    // If right bzet is empty, the bitwise OR will be equal to the left bzet
    else if (right.empty())
        return BZET(*this);
    // Otherwise just operate on them
    else {
        return binop(*((BZET *) this), right, OP_OR);
    }
}

// Bitwise AND
BZET BZET::operator&(BZET& right) {
    // If either bzet is empty, the bitwise AND will be an empty bzet
    if (empty() || right.empty())
        return BZET();
    // Otherwise just operate on them
    else {
        return binop(*((BZET *) this), right, OP_AND);
    }
}

// Bitwise XOR
BZET BZET::operator^(BZET& right) {
    // If left bzet is empty, the bitwise XOR will be equal to the right bzet
    if (empty())
        return BZET(right);
    // If right bzet is empty, the bitwise XOR will be equal to the left bzet
    else if (right.empty())
        return BZET(*this);
    // Otherwise just operate on them
    else {
        return binop(*((BZET *) this), right, OP_XOR);
    }
}

// Comparison operator
bool BZET::operator==(const BZET& right) const {
    if (m_depth != right.m_depth || 
#if NODE_ELS == 4
        m_nnodes != right.m_nnodes ||
        memcmp(m_bzet, right.m_bzet, m_nnodes * sizeof(node_t)))
#else
        m_nhalfnodes != right.m_nhalfnodes ||
        memcmp(m_bzet, right.m_bzet, m_nhalfnodes * sizeof(halfnode_t)))
#endif
        return false;

    return true;
}

// Generic binary operation
BZET BZET::binop(BZET& left, BZET& right, OP op) {
    // Make result bzet
    BZET result;

    // Align both bzets
    align(left, right);

    // Set depth
    result.m_depth = left.m_depth;

    // Operate
    result._binop(left, right, op, left.m_depth);

    // Normalize all bzets
    result.normalize();
    left.normalize();
    right.normalize();
    return result;
}

// -- BIT OPERATIONS --

// Test if a bit is set
bool BZET::at(int64_t bit) const {
    // Test if a bit is set by ANDing the bzet with a bzet with only bit set
    // and checking if the result is empty
    // *this & Bzet(bit)
    if (empty())
        return false;

    BZET temp(bit);
    align(*((BZET *) this), temp);
    return _at(*((BZET *) this), 0, temp, 0, m_depth);
    //return !(*((BZET*) this) & temp).empty();
}

// Set a range of bits
void BZET::setRange(int64_t start, int64_t len) {
    // Create bzet mask
    BZET mask;
    for (int64_t i = start; i < start + len; i++)
        mask.set(i);

    // OR the mask into the original bzet
    *this = *this | mask;
}

// Set a bit
void BZET::set(int64_t bit) {
    // Set a bit by ORing a bzet with bit set into b
    // b | Bzet(bit)
    if (empty())
        *this = BZET(bit);

    BZET temp(bit);
    // Make result bzet
    BZET result;

    // Align both bzets
    align(*this, temp);

    // Operate
    result._binop(*this, temp, OP_OR, m_depth);

    // Normalize
    result.normalize();

    // Shallow copy result to *this
#if NODE_ELS == 4
    m_nnodes = result.m_nnodes;
    m_nbufnodes = result.m_nbufnodes;
#else
    m_nhalfnodes = result.m_nhalfnodes;
    m_nbufhalfnodes = result.m_nbufhalfnodes;
#endif
    free(m_bzet);
    free(m_step);
    m_bzet = result.m_bzet;
    m_step = result.m_step;
    result.m_bzet = NULL;
    result.m_step = NULL;
}

// Sequential set, bit must be greater than the last bit set in the bzet
void BZET::seqset(int64_t bit) {
    if (empty())
        *this = BZET(bit);
    else {
        BZET temp(bit);
        align(*this, temp);
        _seqset(*this, 0, temp, 0, m_depth);
        normalize();
    }
}

// Unset a bit
void BZET::unset(int64_t bit) {
    // Set a bit by ANDing b with the negation of a bzet with bit set
    // b & ~Bzet(bit)
    BZET nb = BZET(bit);
    align(*this, nb);
    nb = ~nb;
    *this = *this & nb;
}

// Get the first bit set
// TODO: Add support for bit literal subtrees
int64_t BZET::firstBit() const {
    if (empty())
        return -1;

#if NODE_ELS == 4
    int level = m_depth;
    int64_t bit = 0;
    size_t loc = 0;

    // Move through the tree to find the first bit
    while (true) {
        // Special calculation at level 1
        if (level == 1) {
            // Check data nodes 0 and 1
            node_t c = m_bzet[loc];

            // If any bit is set in node 0 and 1, first bit found
            if (c)
                for (int i = 7; i >= 0; i--)
                    if ((c >> i) & 1)
                        return bit + 7 - (int64_t) i;

            // If no data bit set in first 2 nodes, check next 2 nodes
            c = m_bzet[loc + 1];
            for (int i = 7; i >= 0; i--)
                if ((c >> i) & 1)
                    return bit + 15 - i;
        }

        // Calculation for all other nodes
        unsigned char data_bits = (m_bzet[loc] >> NODE_ELS) & 0xF;
        unsigned char tree_bits = m_bzet[loc] & 0xF;

        // Check each data bit then the corresponding tree bit
        for (int i = NODE_ELS - 1; i >= 0; i--) {
            // Check data bit
            // If on, return the lowest bit
            if ((data_bits >> i) & 1)
                return bit + (NODE_ELS - 1 - i) * pow4(level);

            // Check tree bit
            if ((tree_bits >> i) & 1) {
                // If on, advance bitBase to the correct prefix and check that node
                bit += (NODE_ELS - 1 - i) * pow4(level);
                break;
            }
        }

        level--;
        loc++;
    }
#else
    int64_t bit = 0;
    int depth = m_depth;
    size_t loc = 0;
    // Move through the tree to find the first bit
    while (true) {
        halfnode_t node_data = m_bzet[loc];
        // At level 0 nodes, first bit will be here somewhere
        if (depth == 0) {
            for (int i = NODE_ELS - 1; i >= 0; i--)
                if ((node_data >> i) & 1)
                    return bit + ((NODE_ELS - 1) - i);
        }

        // For all other nodes, look for first data bit
        halfnode_t node_tree = m_bzet[loc + 1];
        for (int i = NODE_ELS - 1; i >= 0; i--) {
            int data_bit = (node_data >> i) & 1;
            int tree_bit = (node_tree >> i) & 1;

            // If data bit set, first bit found
            if (data_bit) {
                return bit + ((NODE_ELS - 1) - i) * POW(depth);
            }
            // If tree bit is on, break out and examine subtree for first bit
            else if (tree_bit) {
                // Add weight of empty subtrees before this element as bit offset
                bit += ((NODE_ELS - 1) - i) * POW(depth);
                break;
            }
        }

        // Move on to subtree
        depth--;
        loc += 2;
    }
#endif
}

// Bzet_LAST(b)
// TODO: Add support for bit literal subtrees
int64_t BZET::lastBit() const {
    if (empty())
        return -1;

#if NODE_ELS == 4
    // TODO: Fix lastBit() for Bzet4
    size_t loc = 0; 
    int level = m_depth; 
    int64_t bit = 0;

    do {
        // Cpecial calculation at level 1
        if (level == 1) {
            // Check data nodes 2 and 3
            node_t c = m_bzet[loc + 1];

            if (c)
                for (int i = 0; i < 8; ++i)
                    if ((c >> i) & 1)
                        return bit + 15 - i;

            // If no data bit set in last 2 nodes, check first 2
            c = m_bzet[loc];
            for (int i = 0; i < 8; ++i)
                if ((c >> i) & 1)
                    return bit + 7 - i;
        }

        node_t c = m_bzet[loc];
        node_t data_bits = (c >> NODE_ELS) & 0xF;
        node_t tree_bits = c & 0xF;

        //get last tree and data bit
        int last_tree_bit = 0;
        int last_data_bit = 0;
        for (int i = NODE_ELS - 1; i >= 0; i--) {
            if ((tree_bits >> i) & 1)
                last_tree_bit = NODE_ELS - 1 - i;
            if ((data_bits >> i) & 1)
                last_data_bit = NODE_ELS - 1 - i;
        }

        // If this node has the farthest data bit set
        // If so, last bit found
        if (last_data_bit > last_tree_bit) {
            return bit + (last_data_bit + 1) * pow4(level) - 1;
        }

        // If no tree bits are on, we're at the end
        if (!tree_bits) {
            //get the location of the last data bit
            int last_data_bit = 0;
            for (int i = NODE_ELS - 1; i >= 0; i--)
                if ((data_bits >> i) & 1)
                    last_data_bit = NODE_ELS - 1 - i;

            return bit + (last_data_bit + 1) * pow4(level) - 1;
        }

        // Otherwise step through the nodes until you get to the node
        // corresponding to the last tree bit on in the current node
        ++loc;
        for (int i = NODE_ELS - 1; i > NODE_ELS - 1 - last_tree_bit; i--) 
            if ((tree_bits >> i) & 1) 
                loc = step_through(loc, level - 1);

        bit += last_tree_bit * pow4(level);

        // Loc now points to location of the node corresponding to the last
        // tree bit, so repeat this process
        level--;
    } while (loc < m_nnodes);
#else
    int64_t bit = 0;
    int depth = m_depth;
    size_t loc = 0;
    // Move through the tree to find the first bit
    while (true) {
        halfnode_t node_data = m_bzet[loc];
        // At level 0 nodes, last bit will be here somewhere
        if (depth == 0) {
            for (int i = 0; i < NODE_ELS; i++)
                if ((node_data >> i) & 1)
                    return bit + ((NODE_ELS - 1) - i);
        }

        // For all other nodes, look for first data bit
        halfnode_t node_tree = m_bzet[loc + 1];

        // Get position of last tree and data bits
        int last_tree_bit = 0;
        int last_data_bit = 0;
        int skipped_subtrees = -1;
        for (int i = NODE_ELS - 1; i >= 0; i--) {
            if ((node_data >> i) & 1)
                last_data_bit = (NODE_ELS - 1) - i;
            if ((node_tree >> i) & 1) {
                last_tree_bit = (NODE_ELS - 1) - i;
                skipped_subtrees++;
            }
        }

        // If last data bit occurs after last tree bit, last bit found
        if (last_data_bit > last_tree_bit)
            return bit + (last_data_bit + 1) * POW(depth) - 1;

        // Otherwise last data bit is in the final subtree
        // Skip subtrees to get to last subtree
        for (int i = 0; i < skipped_subtrees; i++)
            loc = step_through(loc + 2, depth - 1) - 2;
        
        // Add bits skipped to bit offset
        bit += (last_tree_bit) * POW(depth);

        // Move on to last subtree
        depth--;
        loc += 2;
    }
#endif
}

// Count the number of bits set
int64_t BZET::count() const {
    return _count(0, m_depth);
}

// -- OTHER BZET OPERATIONS --

// Get depth of bzet
int BZET::depth() const {
    return (int) m_depth;
}

// Get size of bzet
size_t BZET::size() const {
#if NODE_ELS == 4
    return (m_nnodes + 1);
#else
    return (m_nhalfnodes * sizeof(halfnode_t) + 1);
#endif
}

// Prints out Bzet tree in human-readable form
void BZET::printBzet(int stdOffset, FILE* target) const {
    _printBzet(stdOffset, target, m_depth);
}

// Copy canonical form of bzet to supplied buffer
void BZET::hex(void *buf) const {
    unsigned char *charbuf = (unsigned char *) buf;
    charbuf[0] = m_depth;
#if NODE_ELS == 4
    memcpy(charbuf + 1, m_bzet, m_nnodes * sizeof(node_t));
#else
    memcpy(charbuf + 1, m_bzet, m_nhalfnodes * sizeof(halfnode_t));
#endif
}

void BZET::exportTo(FILE* stream) const {
    // Write depth out
    fwrite((void *) &m_depth, 1, 1, stream);

#if NODE_ELS == 4
    fwrite(m_bzet, sizeof(node_t), m_nnodes, stream);
#else
    // Write out bzet in little endian
    int endian_test = 1;
    // Little endian
    if (((unsigned char *) &endian_test)[0] == 1) {
        fwrite(m_bzet, sizeof(halfnode_t), m_nhalfnodes, stream);
    }
    // Big endian
    else {
        for (size_t i = 0; i < m_nhalfnodes; i++) {
            // Write out current halfnode in little endian
            halfnode_t curhalf = m_bzet[i];
            for (size_t byte = 0; byte < sizeof(halfnode_t); byte++) {
                fwrite((unsigned char *) &curhalf + byte, 1, 1, stream);
            }
        }
    }
#endif
}

void BZET::importFrom(FILE* stream, size_t size) {
    // Read depth byte
    fread(&m_depth, 1, 1, stream);

#if NODE_ELS == 4
    resize((size - 1) / sizeof(node_t));

    // Read in bzet
    fread(m_bzet, 1, size - 1, stream);
#else
    // Resize to accommodate bzet
    resize((size - 1) / sizeof(halfnode_t));

    // Read in bzet
    int endian_test = 1;
    // Little endian
    if (((unsigned char *) &endian_test)[0] == 1) {
        fread(m_bzet, 1, size - 1, stream);
    }
    // Big endian
    else {
        for (size_t i = 0; i < m_nhalfnodes; i++) {
            // Read in halfnodes as big endian
            halfnode_t curhalf = 0;
            unsigned char bytestage = 0x00;
            for (size_t byte = 0; byte < sizeof(halfnode_t); byte++) {
                fread(&bytestage, 1, 1, stream);
                curhalf |= (halfnode_t) bytestage << (byte * 8);
            }
            m_bzet[i] = curhalf;
        }
    }
#endif

    // Build m_step
    build_step(0, m_depth);
}

// Clears the bzet
void BZET::clear() {
    resize(0);
    m_depth = 0;
}

// Checks whether or not the bzet is empty
bool BZET::empty() const {
#if NODE_ELS == 4
    return (m_nnodes == 0);
#else
    return (m_nhalfnodes == 0);
#endif
}

// Get indices of bits set in the bzet
int64_t BZET::getBits(int64_t* bits, int64_t limit, int64_t start) {
    int64_t bit;
    size_t loc = 0;

    // Set number of bits to get, which is the smaller of the
    // number of bits set in b and limit (if set)
    int64_t bitcount = count();
    limit = limit ? ((limit > bitcount) ? bitcount : limit) : bitcount;

    // Clone this bzet to use to get bits
    BZET b_copy(*this);

    for (int64_t i = 0; i < limit; i++) {
        // Get the first bit set and commit to bits
        bit = b_copy.firstBit();
#ifdef DEBUG
        if (bit < 0) {
            printf("getbits: bit < 0 (%d of limit %d)", i, limit);
            display_error("", true);
        }
#endif
        // Unset first bit
        b_copy.unset(bit);

        if (bit >= start) {
            bits[loc] = bit;
            loc++;
        }
    }

    return limit;
}


// Auxiliary functions

#ifdef USE_LITERAL
/*
NODETYPE bitstotree(BZET& b, int depth, unsigned char *data, size_t len) {
    if (len < POW(depth) / 8)
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
    size_t subtree_size = POW(depth - 1) / 8;
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
    if ((b->nhalfnodes - loc) * sizeof(halfnode_t) >= POW(depth + 1)) {
        size_t locstart = loc;
        unsigned char *buf = (unsigned char *) malloc(POW(depth + 1) * sizeof(halfnode_t));
        // TODO: More elegant error handling
        if (!buf)
            display_error("malloc failed in bitstotree", true);

        // Convert tree to bytes
        treetobits(buf, b->bzet + loc, depth);

        unsigned char *bufptr = buf;

        // Write back as necessary
        // For each halfnode's worth
        for (int i = 0; i < POW(depth + 1) / sizeof(halfnode_t); i++) {
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

void treetobits(halfnode_t *dst, halfnode_t *node, int depth) {
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
            size_t nodes = POW(depth - 1) / 8 / sizeof(halfnode_t);
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
            memset(buf, 0xFF, POW(depth - 1) / 8);
            buf += POW(depth - 1) / 8;
        }
        // Otherwise this element is empty
        else {
            memset(buf, 0x00, POW(depth - 1) / 8);
            buf += POW(depth - 1) / 8;
        }
    }
}
*/
#endif

NODETYPE BZET::_binop(BZET& left, BZET& right, OP op, int lev, size_t left_loc, size_t right_loc) {
#if NODE_ELS == 4
    // Handle level 1 operations
    if (lev == 1) {
        node_t c_left = left.m_bzet[left_loc];
        node_t c_right = right.m_bzet[right_loc];

        // Build data nodes 0 and 1
        node_t result_c1 = 0x00;
        for (int i = 7; i >= 0; --i) {
            int bita = (c_left >> i) & 1;
            int bitb = (c_right >> i) & 1;
            result_c1 |= (do_data_op(op, bita, bitb) << i);
        }

        // Get corresponding data nodes 2 and 3
        c_left = left.m_bzet[left_loc + 1];
        c_right = right.m_bzet[right_loc + 1];

        // Build data nodes 2 and 3
        node_t result_c2 = 0x00;
        for (int i = 7; i >= 0; --i) {
            int bita = (c_left >> i) & 1;
            int bitb = (c_right >> i) & 1;
            result_c2 |= (do_data_op(op, bita, bitb) << i);
        }

        // Saturated node
        if (result_c1 == 0xFF && result_c2 == 0xFF) {
            // Build bzet manually if bzet is saturated
            if (empty()) {
                resize(2);
                m_bzet[0]++;
                m_bzet[1] = 0x80;
                m_step[1] = 1;
                return NORMAL;
            }

            return SATURATED;
        }
        // Empty node
        else if (!result_c1 && !result_c2) {
            // Build bzet manually if bzet is empty
            if (empty()) {
                clear();
                return NORMAL;
            }

            return EMPTY;
        }

        size_t loc = m_nnodes;

        // Otherwise commit data
        // Create two nodes to accommodate data nodes
        resize(m_nnodes + 2);

        // Write result
        m_bzet[loc] = result_c1;
        m_bzet[loc + 1] = result_c2;

        /*// Set steps
        m_step[loc] = 2;
        m_step[loc + 1] = 1;*/
        
        return NORMAL;
    }
#else
    // Handle level 0 bit operations
    if (lev == 0) {
        halfnode_t node_data = 0;
        halfnode_t left_data = left.m_bzet[left_loc];
        halfnode_t right_data = right.m_bzet[right_loc];
        // Compute data node
        for (int i = NODE_ELS - 1; i >= 0; i--) {
            int left_data_bit = (left_data >> i) & 0x1;
            int right_data_bit = (right_data >> i) & 0x1;
            node_data = (node_data << 1) | (halfnode_t) do_data_op(op, left_data_bit, right_data_bit);
        }

        // Empty node
        if (node_data == 0) {
            return EMPTY;
        }
        // Saturated node
        else if (node_data == (halfnode_t) -1) {
            // If bzet is empty, commit the node
            if (m_nhalfnodes == 0) {
                resize(1);
                m_bzet[0] = node_data;
                m_step[0] = 0x1;
            }
            return SATURATED;
        }
        // Normal node, commit
        else {
            size_t loc = m_nhalfnodes;
            resize(loc + 1);
            m_bzet[loc] = node_data;
            m_step[loc] = 0x1;
            return NORMAL;
        }
    }
#endif

    // Get corresponding nodes of the left and right tree
#if NODE_ELS == 4
    node_t c_left = left.m_bzet[left_loc];
    node_t c_right = right.m_bzet[right_loc];
    left_loc++;
    right_loc++;
#else 
    halfnode_t c_left_data = left.m_bzet[left_loc];
    halfnode_t c_left_tree = left.m_bzet[left_loc + 1];
    halfnode_t c_right_data = right.m_bzet[right_loc];
    halfnode_t c_right_tree = right.m_bzet[right_loc + 1];
    left_loc += 2;
    right_loc += 2;
#endif

    // Create new node to hold data
#if NODE_ELS == 4
    node_t new_node = 0;
#else
    halfnode_t node_data = 0;
    halfnode_t node_tree = 0;
#endif

    // Reserve space for this node
#if NODE_ELS == 4
    size_t loc = m_nnodes;
    resize(m_nnodes + 1);
#else
    size_t loc = m_nhalfnodes;
    resize(m_nhalfnodes + 2);
#endif

    // For each element in the node
    for (int i = NODE_ELS - 1; i >= 0; i--) {
#if NODE_ELS == 4
        int cur_left_tree_bit = (c_left >> i) & 0x1;
        int cur_left_data_bit = (c_left >> (i + NODE_ELS)) & 0x1;
        int cur_right_tree_bit = (c_right >> i) & 0x1;
        int cur_right_data_bit = (c_right >> (i + NODE_ELS)) & 0x1;
#else
        int cur_left_tree_bit = (c_left_tree >> i) & 0x1;
        int cur_left_data_bit = (c_left_data >> i) & 0x1;
        int cur_right_tree_bit = (c_right_tree >> i) & 0x1;
        int cur_right_data_bit = (c_right_data >> i) & 0x1;
#endif
        
        // TT: If both tree bits are on
        if (cur_left_tree_bit && cur_right_tree_bit) {
#ifdef USE_LITERAL
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
#endif
                // Recurse
                NODETYPE cn = _binop(left, right, op, lev - 1, left_loc, right_loc);

                // Saturated subtree
                if (cn == SATURATED) {
                    // Turn on data bit
#if NODE_ELS == 4
                    new_node |= 0x10 << i;
#else
                    node_data = (node_data << 1) | 0x1;
                    node_tree <<= 1;
#endif
                }
                // Empty subtree
                else if (cn == EMPTY) {
#if NODE_ELS == 4
                    // Do nothing, element is already set to 0
#else
                    // Shift data and tree nodes over
                    node_data <<= 1;
                    node_tree <<= 1;
#endif
                }
                // "Tree" subtree
                else if (cn == NORMAL) {
                    // Turn on tree bit
#if NODE_ELS == 4
                    // Subtree exists, turn on tree bit
                    new_node |= 0x01 << i;
#else
                    node_data <<= 1;
                    node_tree = (node_tree << 1) | 0x1;
#endif
                }
#ifdef USE_LITERAL
                // Otherwise data literal subtree
                else {
                    // TODO: Handle data literal subtree
                }
#endif

                // Advance location pointers
                left_loc = left.step_through(left_loc, lev - 1);
                right_loc = right.step_through(right_loc, lev - 1);
#ifdef USE_LITERAL
            }
#endif
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
                    left_loc = left.step_through(left_loc, lev - 1);

                    // Data bit is already 0
#if NODE_ELS == 4
#else
                    node_data <<= 1;
                    node_tree <<= 1;
#endif

                    break;

                // Delete right subtree, set data bit off
                case DB0:
                    // Skip the right subtree
                    right_loc = right.step_through(right_loc, lev - 1);

                    // Data bit is already 0
#if NODE_ELS == 4
#else
                    node_data <<= 1;
                    node_tree <<= 1;
#endif

                    break;

                // Delete left subtree, set data bit on
                case DA1: 
                    // Skip the left subtree
                    left_loc = left.step_through(left_loc, lev - 1);

                    // Turn on data bit
#if NODE_ELS == 4
                    new_node |= 0x80 >> ((NODE_ELS - 1) - i);
#else
                    node_data = (node_data << 1) | 0x1;
                    node_tree <<= 1;
#endif

                    break;

                // Delete right subtree, set data bit on
                case DB1:
                    // Skip the right subtree
                    right_loc = right.step_through(right_loc, lev - 1);

                    // Turn on data bit
#if NODE_ELS == 4
                    new_node |= 0x80 >> ((NODE_ELS - 1) - i);
#else
                    node_data = (node_data << 1) | 0x1;
                    node_tree <<= 1;
#endif

                    break;

                // Copy left subtree into result
                case CA: 
                    // Append left subtree
                    append_subtree(left, left_loc, lev - 1);

                    // Move through left subtree
                    left_loc = left.step_through(left_loc, lev - 1);

                    // Turn on tree bit
#if NODE_ELS == 4
                    new_node |= 0x08 >> ((NODE_ELS - 1) - i);
#else
                    node_data <<= 1;
                    node_tree = (node_tree << 1) | 0x1;
#endif

                    break;

                // Copy right subtree into result
                case CB:
                    // Append right subtree
                    append_subtree(right, right_loc, lev - 1);

                    // Move through right subtree
                    right_loc = right.step_through(right_loc, lev - 1);

                    // Turn on tree bit
#if NODE_ELS == 4
                    new_node |= 0x08 >> ((NODE_ELS - 1) - i);
#else
                    node_data <<= 1;
                    node_tree = (node_tree << 1) | 0x1;
#endif

                    break;

                // Copy left subtree into result and negate
                case NA:
#if NODE_ELS == 4
                    end = m_nnodes;
#else
                    end = m_nhalfnodes;
#endif

                    // Append left subtree
                    append_subtree(left, left_loc, lev - 1);

                    // Negate
                    subtree_not(end, lev - 1);

                    // Move through left subtree
                    left_loc = left.step_through(left_loc, lev - 1);

                    // Turn on tree bit
#if NODE_ELS == 4
                    new_node |= 0x08 >> ((NODE_ELS - 1) - i);
#else
                    node_data <<= 1;
                    node_tree = (node_tree << 1) | 0x1;
#endif

                    break;

                // Copy right subtree into result and negate
                case NB:
#if NODE_ELS == 4
                    end = m_nnodes;
#else
                    end = m_nhalfnodes;
#endif

                    // Append right subtree
                    append_subtree(right, right_loc, lev - 1);

                    // Negate
                    subtree_not(end, lev - 1);

                    // Move through right subtree
                    right_loc = right.step_through(right_loc, lev - 1);

                    // Turn on tree bit
#if NODE_ELS == 4
                    new_node |= 0x08 >> ((NODE_ELS - 1) - i);
#else
                    node_data <<= 1;
                    node_tree = (node_tree << 1) | 0x1;
#endif
                    break;

                default:
                    // Should be impossible to get here, but just in case?
                    display_error("_binop: Something went terribly, terribly wrong", true);
                    break;
            }
        } 
        // Only data bits
#ifdef USE_LITERAL
        else {
#else
        else if (!cur_left_tree_bit && !cur_right_tree_bit) {
#endif
            // Shift in data bit
#if NODE_ELS == 4
            if (do_data_op(op, cur_left_data_bit, cur_right_data_bit)) {
                new_node |= 0x10 << i;
            }
#else
            node_data = (node_data << 1) | (halfnode_t) do_data_op(op, cur_left_data_bit, cur_right_data_bit);
            node_tree <<= 1;
#endif
        }
#ifndef USE_LITERAL
        else {
            display_error("_binop: Invalid bit state (data bit and tree bit set simultaneously)", true);
        }
#endif
    }

    // Write nodes
#if NODE_ELS == 4
    m_bzet[loc] = new_node;
#else
    m_bzet[loc] = node_data;
    m_bzet[loc + 1] = node_tree;
#endif

    // Set step
    set_step(loc, lev);

#if NODE_ELS == 4
    if (m_bzet[loc] == 0x00) {
        // Nuild bzet manually if bzet is empty
        if (m_nnodes == 1) {
            clear();
            return NORMAL;
        }

        resize(m_nnodes - 1);
        return EMPTY;
    }
    // Resulting node is saturated
    else if (m_bzet[loc] == 0xF0) {
        if (m_nnodes == 1) {
            resize(2);
            m_bzet[0]++;
            m_bzet[1] = 0x80;
            m_step[1] = 1;
            return NORMAL;
        }

        resize(m_nnodes - 1);
        return SATURATED;
    }
#else
    // If resulting node is empty
    if (node_tree == 0 && node_data == 0) {
        // Drop newly committed node if not root node
        if (m_nhalfnodes > 2) {
            resize(m_nhalfnodes - 2);
        }
        // If root node, the bzet is empty
        else {
            clear();
            return NORMAL;
        }

        return EMPTY;
    }
    // If resulting node is saturated
    else if (node_tree == 0 && node_data == (halfnode_t) -1) {
        // Drop newly committed node if not root node
        if (m_nhalfnodes > 2) {
            resize(m_nhalfnodes - 2);
        }

        return SATURATED;
    }
#endif

    // TODO: See if collapsible to bit literal

    return NORMAL;
}

void BZET::_printBzet(int stdOffset, FILE* target, int depth, size_t loc, int offset, bool pad) const {
    // Print level info
    if (loc == 0) {
        fprintf(target, "%.2XL", m_depth);

        // If empty, we're done.
        if (empty())
            fprintf(target, "\n");
    }

#if NODE_ELS == 4
    // No reading past the bzet!
    if (loc >= size())
        return;

    unsigned char c = m_bzet[loc];
    unsigned char data_bits = (c >> NODE_ELS) & 0xF;
    unsigned char tree_bits = c & 0xF;

    // Print offset if any
    if (pad) {
        int blanks = 5*offset + 3 + stdOffset;
        for (int j = 0; j < blanks; j++)
            fprintf(target, " ");
    }

    // If depth is 1, print special level 1 interpretation
    if (depth == 1) {   
        fprintf(target, "D(%.2X%.2X)\n", m_bzet[loc], m_bzet[loc + 1]);
    } 
    // If no tree bits on, then this is a data node
    else if ((c & 0xF0) == c) {
        fprintf(target, "D(%X)\n", data_bits);
    }
    else {     
        // Print the current node
        fprintf(target, "[%X-%X]", data_bits, tree_bits);

        // Recursively print child nodes
        bool firstNode = true;
        depth--;
        for (int i = NODE_ELS - 1; i >= 0; --i) {
            // If tree bit set
            if ((tree_bits >> i) & 0x1) {
                // Print first node without offset
                if (firstNode) {
                    _printBzet(stdOffset, target, depth, loc + 1, offset + 1);
                    firstNode = false;
                } else {
                    _printBzet(stdOffset, target, depth, step_through(loc + 1, depth), offset + 1, true);
                    loc = step_through(loc + 1, depth) - 1;
                }
            }
        }
    }
#else
    // No reading past the bzet!
    if (loc >= m_nhalfnodes)
        return;

    halfnode_t data_bits = m_bzet[loc];
    halfnode_t tree_bits = m_bzet[loc + 1];

    // Print offset if any
    if (pad) {
        // First 3 spaces to align with XXL level byte
        // (sizeof(halfnode_t) * 4 + 3)*offset for spaces an internal node takes up
        //     sizeof(halfnode_t) - number of bytes printed out
        //     3 - "pretty" part: [?-?]
        int blanks = 3 + (sizeof(halfnode_t) * 4 + 3) * offset + stdOffset;
        for (int j = 0; j < blanks; j++)
            fprintf(target, " ");
    }

    // If depth is 0 or no tree bits are set, this is a data node
    if (depth == 0 || !tree_bits) {
        fprintf(target, "D(%.*X)\n", sizeof(halfnode_t) * 2, data_bits);
        /*if (m_step[loc] != 1) {
            printf("data node with step != 1 at %d, nhalf=%d\n", loc, m_nhalfnodes);
            exit(1);
        }*/
    }
    // Otherwise this is a tree node
    else {     
        // Print the current node
        fprintf(target, "[%.*X-%.*X]", sizeof(halfnode_t) * 2, data_bits, 
            sizeof(halfnode_t) * 2, tree_bits);

        // Recursively print subtrees
        bool firstNode = true;
        depth--;
        loc += 2;
        for (int i = NODE_ELS - 1; i >= 0; i--) {
            // TODO: Use popcount
            // If tree bit set
            if ((tree_bits >> i) & 0x1) {
                // Print first node without offset
                if (firstNode) {
                    _printBzet(stdOffset, target, depth, loc, offset + 1);
                    firstNode = false;
                } 
                else {
                    _printBzet(stdOffset, target, depth, loc, offset + 1, true);
                }
                loc = step_through(loc, depth);
            }
        }
    }
#endif
}

void BZET::align(BZET& b1, BZET& b2) {
    // Get difference in depths
    int diffdepth = abs(b1.m_depth - b2.m_depth);
    // Nothing to be done if they are the same depth
    if (diffdepth == 0)
        return;

    // If either is empty, there's no need to align
    if (b1.empty() || b2.empty())
        return;

    // If b2 needs to be grown
    if (b1.m_depth > b2.m_depth) {
        b2.m_depth = b1.m_depth;
#if NODE_ELS == 4
        size_t old_size = b2.m_nnodes;

        // Resize b2 to accommodate new heading nodes
        b2.resize(b2.m_nnodes + diffdepth);

        // Move bzet and step in b2 to accommodate new heading nodes
        memmove(b2.m_bzet + diffdepth, b2.m_bzet, old_size * sizeof(node_t));
        memmove(b2.m_step + diffdepth, b2.m_step, old_size * sizeof(step_t));

        // Add new nodes and new step
        size_t loc = 0;
        for (int i = 0; i < diffdepth; i++) {
            b2.m_bzet[loc] = (node_t) ((node_t) 0x1 << (NODE_ELS - 1));
            b2.set_step(loc, 2); // Cheat with this parameter since depth != 1

            loc++;
        }
#else
        size_t old_size = b2.m_nhalfnodes;

        // Resize b2 to accommodate new heading nodes
        b2.resize(b2.m_nhalfnodes + diffdepth * 2);

        // Move bzet and step in b2 to accommodate new heading nodes
        memmove(b2.m_bzet + diffdepth * 2, b2.m_bzet, old_size * sizeof(halfnode_t));
        memmove(b2.m_step + diffdepth * 2, b2.m_step, old_size * sizeof(step_t));

        // Add new nodes and new step
        size_t loc = 0;
        for (int i = 0; i < diffdepth; i++) {
            b2.m_bzet[loc] = 0x80;
            b2.set_step(loc, 2); // Cheat with this parameter since depth != 0

            loc += 2;
        }
#endif
    }
    // Otherwise b1 needs to be grown
    else {
        b1.m_depth = b2.m_depth;
#if NODE_ELS == 4
        size_t old_size = b1.m_nnodes;

        // Resize b1 to accommodate new heading nodes
        b1.resize(b1.m_nnodes + diffdepth);

        // Move bzet and step in b1 to accommodate new heading nodes
        memmove(b1.m_bzet + diffdepth, b1.m_bzet, old_size * sizeof(node_t));
        memmove(b1.m_step + diffdepth, b1.m_step, old_size * sizeof(step_t));

        // Add new nodes and new step
        size_t loc = 0;
        for (int i = 0; i < diffdepth; i++) {
            b1.m_bzet[loc] = (node_t) ((node_t) 0x1 << (NODE_ELS - 1));
            b1.set_step(loc, 2); // Cheat with this parameter since depth != 1

            loc++;
        }
#else
        size_t old_size = b1.m_nhalfnodes;

        // Resize b1 to accommodate new heading nodes
        b1.resize(b1.m_nhalfnodes + diffdepth * 2);

        // Move bzet and step in b1 to accommodate new heading nodes
        memmove(b1.m_bzet + diffdepth * 2, b1.m_bzet, old_size * sizeof(halfnode_t));
        memmove(b1.m_step + diffdepth * 2, b1.m_step, old_size * sizeof(step_t));

        // Add new nodes and new step
        size_t loc = 0;
        for (int i = 0; i < diffdepth; i++) {
            b1.m_bzet[loc] = 0;
            b1.m_bzet[loc + 1] = (halfnode_t) ((halfnode_t) 0x1 << (NODE_ELS - 1));
            b1.set_step(loc, 2); // Cheat with this paramter since depth != 0

            loc += 2;
        }
#endif
    }
}

void BZET::normalize() {
#if NODE_ELS == 4
    //no need to normalize bzet with only 1 level
    if (m_bzet[0] == 0x01)
        return;

    //count leading 0x08 nodes
    int leading = 0;
    while (m_bzet[leading] == 0x08)
        leading++;

    // If there are leading 0x08 nodes, get rid of them
    if (leading) {
        // Modify depth
        m_depth -= (unsigned char) leading;

        // Move bzet and step
        memmove(m_bzet, m_bzet + leading, (m_nnodes - leading) * sizeof(node_t));
        memmove(m_step, m_step + leading, (m_nnodes - leading) * sizeof(step_t));

        // Resize
        resize(m_nnodes - leading);
    }
#else
    int leading = 0;
    int loc = 0;
    int depth = m_depth;
    // Count leading nodes that can be stripped.
    // This occurs when the data node is all zero and only the leftmost tree bit is set
    while (depth > 0 && m_bzet[loc] == 0 && 
           m_bzet[loc + 1] == (halfnode_t) ((halfnode_t) 0x1 << (NODE_ELS - 1))) {
        leading++;
        loc += 2;
    }

    // If there are leading nodes that can be stripped
    if (leading) {
        // Modify depth
        m_depth -= (unsigned char) leading;

        // Move bzet and step
        memmove(m_bzet, m_bzet + loc, (m_nhalfnodes - loc) * sizeof(halfnode_t));
        memmove(m_step, m_step + loc, (m_nhalfnodes - loc) * sizeof(step_t));

        // Resize
        resize(m_nhalfnodes - loc);
    }
#endif
}

size_t BZET::step_through(size_t loc, int depth) const {
    // Make sure loc is in range
#if NODE_ELS == 4
    if (loc >= m_nnodes) {
        printf("stepthrough fail trying %d, nnodes=%d\n", loc, m_nnodes);
#else
    if (loc >= m_nhalfnodes) {
        printf("stepthrough fail trying %d, nhalf=%d\n", loc, m_nhalfnodes);
#endif
        display_error("", true);
        //return -1;
    }

#if NODE_ELS == 4
    if (depth == 1)
        return loc + 2;
#else
    if (depth == 0)
        return loc + 1;
#endif

    // Get step offset
#if NODE_ELS == 4
    size_t step_offset = m_step[loc];
#else
    size_t step_offset = ((size_t) m_step[loc]) * STEP_T_MAX + (size_t) m_step[loc + 1];
#endif

    // If the step offset is 0, the offset is too large to be actually stored
    // Compute it by examining the step of subtrees stemming from this node
    if (step_offset == 0) {
#if NODE_ELS == 4
        // Get node tree bits
        node_t tree_bits = m_bzet[loc] & 0xF;

        // Advance to location of first subtree node
        loc++;
#else
        // Get node tree bits
        halfnode_t tree_bits = m_bzet[loc + 1];

        // Advance to location of first subtree node
        loc += 2;
#endif

        // Step through each subtree
        // TODO: Performance gain by using popcount
        for (int i = NODE_ELS - 1; i >= 0; i--) {
            if ((tree_bits >> i) & 1)
                loc = step_through(loc, depth - 1);
        }
        return loc;
    }

    // Step offset is nonzero, just add it loc and return it
    return loc + step_offset;
}

void BZET::subtree_not(size_t loc, int depth) {
#if NODE_ELS == 4
    if (loc >= m_nnodes)
        return;

    // If no tree bits are on, just not the data bits
    if (!(m_bzet[loc] & 0xF) && depth > 1) {
        m_bzet[loc] = ~(m_bzet[loc] | 0xF);
        return;
    }

    // If depth is 1, we just not this and the next node
    if (depth == 1) {
        m_bzet[loc] = ~m_bzet[loc];
        m_bzet[loc + 1] = ~m_bzet[loc + 1];
        return;
    }

    size_t nextLoc = loc + 1;

    // Decrement depth
    --depth;

    for (int i = NODE_ELS - 1; i >= 0; --i) {
        // If tree bit is on, recursively not each node
        if ((m_bzet[loc] >> i) & 1) {
            subtree_not(nextLoc, depth);
            nextLoc = step_through(nextLoc, depth - 1);
        }
    }

    // Get the current byte
    node_t c = m_bzet[loc]; 

    // New data bits
    node_t data_bits = (~(c >> NODE_ELS)) << NODE_ELS;

    // Add tree bits
    node_t newc = data_bits | (c & 0xF);

    // Set up for dusting with tree bits priority
    node_t tree_bits = c & 0xF;

    // Dust using tree bits and write back
    m_bzet[loc] = newc & (~(tree_bits << NODE_ELS) | tree_bits);
#else
    // Check if loc is in range
    if (loc >= m_nhalfnodes) {
        printf("subtree_not: loc >= m_halfnodes %d %d\n", loc, m_nhalfnodes);
        exit(1);
        return;
    }

    // For level 0 nodes, simply NOT it
    if (depth == 0) {
        m_bzet[loc] = ~m_bzet[loc];
        return;
    }

    // Handling all other nodes

    // Get data and tree bits
    halfnode_t data_bits = m_bzet[loc];
    halfnode_t tree_bits = m_bzet[loc + 1];

    // Make result data bits by bitwise NOT
    m_bzet[loc] = (~data_bits & ~tree_bits) | (data_bits & tree_bits);

    // Repeat recursively for all subtrees
    loc += 2;
    // TODO: Performance gain by using popcount for tree bits
    for (int i = 0; i < NODE_ELS; i++) {
        // If tree bit set
        if ((tree_bits >> i) & 1) {
            // If subtree is in tree form
            if (((data_bits >> i) & 1) == 0) {
                subtree_not(loc, depth - 1);
                loc = step_through(loc, depth - 1);
            }
            // If subtree is in literal form
            else {
                // Get number of halfnodes 
                size_t nnodes = POW(depth - 1) / 8 / sizeof(halfnode_t);

                // NOT each halfnode
                for (size_t i = 0; i < nnodes; i++) {
                    m_bzet[loc] = ~m_bzet[loc];
                    loc++;
                }
            }
        }
    }
#endif
}

// TODO: use popcount, add support for bit literal subtree
int64_t BZET::_count(size_t loc, int depth) const {
    if (empty())
        return 0;

    int64_t count = 0;

#if NODE_ELS == 4
    // Just return bit count with weight 1
    if (depth == 1) {
        node_t data = m_bzet[loc];
        for (int i = 0; i < 8; i++)
            count += ((data >> i) & 1);

        data = m_bzet[loc + 1];
        for (int i = 0; i < 8; i++)
            count += ((data >> i) & 1);
#else
    // Just return bit count with weight 1
    if (depth == 0) {
        halfnode_t data = m_bzet[loc];
        for (int i = 0; i < NODE_ELS; i++)
            count += ((data >> i) & 1);
#endif

        return count;
    }

#if NODE_ELS == 4
    // Handle all other levels
    node_t data_node = m_bzet[loc] >> NODE_ELS;
    node_t tree_node = m_bzet[loc] & 0xF;
    for (int i = NODE_ELS - 1; i >= 0; i--) {
        // TODO: Add handling for variable data literal
        int data_bit = (data_node >> i) & 1;
        int tree_bit = (tree_node >> i) & 1;

        // If data bit set, add its weighted count to running count
        if (data_bit) {
            count += POW(depth);
        }
        // If tree bit set, recurse
        else if (tree_bit) {
            count += _count(loc + 1, depth - 1);
            loc = step_through(loc + 1, depth - 1) - 1;
        }
    }
#else
    // Handle all other levels
    halfnode_t data_node = m_bzet[loc];
    halfnode_t tree_node = m_bzet[loc + 1];
    for (int i = NODE_ELS - 1; i >= 0; i--) {
        // TODO: Add handling for variable data literal
        int data_bit = (data_node >> i) & 1;
        int tree_bit = (tree_node >> i) & 1;

        // If data bit set, add its weighted count to running count
        if (data_bit) {
            count += POW(depth);
        }
        // If tree bit set, recurse
        else if (tree_bit) {
            count += _count(loc + 2, depth - 1);
            loc = step_through(loc + 2, depth - 1) - 2;
        }
    }
#endif

    return count;
}

// Implementation of seqset
NODETYPE BZET::_seqset(BZET& b, size_t locb, BZET& right, size_t locright, int depth) {
    // At depth 0, we are at the lowest level, set the bit and finish
#if NODE_ELS == 4
    if (depth == 1) {
#else
    if (depth == 0) {
#endif
        b.m_bzet[locb] |= right.m_bzet[locright];
#if NODE_ELS == 4
        b.m_bzet[locb + 1] |= right.m_bzet[locright];
        // If node is saturated, drop these node and signal it was dropped
        if (b.m_bzet[locb] == 0xFF && b.m_bzet[locb + 1] == 0xFF) {
            // If these are the only nodes
            if (b.m_nnodes == 2) {
                b.resize(1);
                b.m_bzet[0] = 0x80;
                b.m_step[0] = 1;
                b.m_depth = 2;
                return NORMAL;
            }

            b.resize(b.m_nnodes - 1);
#else
        // If node is saturated, drop this node and signal it was dropped
        if (b.m_bzet[locb] == (halfnode_t) -1) {
            // If this is the only node
            if (b.m_nhalfnodes == 1) {
                b.resize(2);
                b.m_bzet[0] = ((halfnode_t) 1 << (NODE_ELS - 1));
                b.m_bzet[1] = 0;
                b.m_step[0] = 2;
                b.m_depth = 1;
                return NORMAL;
            }

            b.resize(b.m_nhalfnodes - 1);
#endif
            return SATURATED;
        }

        return NORMAL;
    }

#if NODE_ELS == 4
    node_t bdata_bits = b.m_bzet[locb] >> NODE_ELS;
    node_t btree_bits = b.m_bzet[locb] & 0xF;
    node_t righttree_bits = right.m_bzet[locright] & 0xF;
#else
    halfnode_t& bdata_bits = b.m_bzet[locb];
    halfnode_t& btree_bits = b.m_bzet[locb + 1];
    //halfnode_t& rightdata_bits = right.m_bzet[locright];
    halfnode_t& righttree_bits = right.m_bzet[locright + 1];
#endif
    size_t loc = locb;

    // Data bit higher up is set, we're done
    if ((bdata_bits | righttree_bits) == bdata_bits) {
        return NORMAL;
    }
    // Same tree bit is set, recurse
    else if ((btree_bits | righttree_bits) == btree_bits) {
        // Skip all subtrees in b before the subtree we want
        // TODO: use popcount to speed up?
#if NODE_ELS == 4
        locb++;
#else
        locb += 2;
#endif
        bool first = true;
        for (int i = 0; i < NODE_ELS; i++)
            if ((btree_bits >> i) & 1) {
                if (!first)
                    locb = b.step_through(locb, depth - 1);
                else
                    first = false;
            }

        // Recurse
#if NODE_ELS == 4
        NODETYPE ret = _seqset(b, locb, right, locright + 1, depth - 1);
#else
        NODETYPE ret = _seqset(b, locb, right, locright + 2, depth - 1);
#endif
        if (ret == SATURATED) {
#if NODE_ELS == 4
            // Flip tree and data bits
            b.m_bzet[loc] ^= ((righttree_bits << NODE_ELS) | righttree_bits);

            // If this node is saturated
            if (b.m_bzet[loc] == 0xFF) {
                // If this is the only node
                if (b.m_nnodes == 1) {
                    b.resize(1);
                    b.m_bzet[0] = 0x80;
                    b.m_step[0] = 1;
                    b.m_depth++;
                }

                b.resize(b.m_nnodes - 1);
#else
            // Flip tree and data bits
            bdata_bits ^= righttree_bits;
            btree_bits ^= righttree_bits;

            // If this node is saturated
            if (b.m_bzet[loc] == (halfnode_t) -1) {
                // If this is the only node
                if (b.m_nhalfnodes == 2) {
                    b.resize(2);
                    b.m_bzet[0] = ((halfnode_t) 1 << (NODE_ELS - 1));
                    b.m_bzet[1] = 0;
                    b.m_step[0] = 2;
                    b.m_depth++;
                    return NORMAL;
                }

                b.resize(b.m_nhalfnodes - 2);
#endif
                return SATURATED;
            }
        }

        // Update step
        // Cheat with depth since depth != 0
        b.set_step(loc, 2);        
        return NORMAL;
    }
    // Tree bit is not on, set tree bit and append the subtree
    else {
#if NODE_ELS == 4
        // Set tree bit
        b.m_bzet[loc] |= righttree_bits;
        // Append subtree
        b.append_subtree(right, locright + 1, depth - 1);
#else
        // Set tree bit
        btree_bits |= righttree_bits;
        // Append subtree
        b.append_subtree(right, locright + 2, depth - 1);
#endif
        // Update step
        // Cheat with depth since depth != 0
        b.set_step(locb, 1);
        return NORMAL;
    }
}

// Implementation of at
bool BZET::_at(BZET& b, size_t locb, BZET& right, size_t locright, int depth) {
#if NODE_ELS == 4
    // At depth 1, we are at the lowest level, check if the bit is set
    if (depth == 1) {
        if ((b.m_bzet[locb] & right.m_bzet[locright]) ||
            (b.m_bzet[locb + 1] & right.m_bzet[locright + 1]))
            return true;
        
        return false;
    }

    node_t bdata_bits = b.m_bzet[locb] >> NODE_ELS;
    node_t btree_bits = b.m_bzet[locb] & 0xF;
    node_t righttree_bits = right.m_bzet[locright] & 0xF;
#else
    // At depth 0, we are at the lowest level, check if the bit is set
    if (depth == 0) {
        if (b.m_bzet[locb] & right.m_bzet[locright])
            return true;

        return false;
    }

    halfnode_t& bdata_bits = b.m_bzet[locb];
    halfnode_t& btree_bits = b.m_bzet[locb + 1];
    //halfnode_t& rightdata_bits = right.m_bzet[locright];
    halfnode_t& righttree_bits = right.m_bzet[locright + 1];
#endif

    // Data bit higher up is set, we're done. Bit found.
    if (bdata_bits & righttree_bits) {
        return true;
    }
    // Same tree bit is set, recurse
    else if (btree_bits & righttree_bits) {
        // Skip all subtrees in b before the subtree we want
#if NODE_ELS == 4
        locb++;
#else
        locb += 2;
#endif
        for (int i = NODE_ELS - 1; i >= 0; i--)
            if ((btree_bits >> i) & (righttree_bits >> i)) {
                break;
            }
            else if ((btree_bits >> i) & 1) {
                locb = b.step_through(locb, depth - 1);
            }

        // Recurse
#if NODE_ELS == 4
        return _at(b, locb, right, locright + 1, depth - 1);
#else
        return _at(b, locb, right, locright + 2, depth - 1);
#endif
    }
    // Tree bit is not on, so bit is not set.
    else {
        return false;
    }
}

size_t BZET::build_step(size_t loc, int depth) {
#if NODE_ELS == 4
    if (depth == 1) {
        m_step[loc] = 2;
        return 2;
    }
#else
    if (depth == 0) {
        m_step[loc] = 1;
        return 1;
    }
#endif

    bool overflow = false;
    size_t curloc = loc;

#if NODE_ELS == 4
    size_t step = 1;
    node_t tree_bits = m_bzet[loc] & 0xF;
    loc++;
#else
    size_t step = 2;
    halfnode_t tree_bits = m_bzet[loc + 1];
    loc += 2;
#endif
    for (size_t i = 0; i < NODE_ELS; i++) {
        if ((tree_bits >> i) & 1) {
            size_t s = build_step(loc, depth - 1);
            if (overflow || s == 0) {
                step = 0;
                overflow = true;
            }
            else
                step += s;
                
            loc = step_through(loc, depth - 1);
        }
    }

#if NODE_ELS == 4
    if (m_step[curloc] >= STEP_T_MAX)
        m_step[curloc] = 0;
    else
        m_step[curloc] = (step_t) STEP_T_MAX;
#else
    m_step[curloc] = (step_t) (step / STEP_T_MAX);
    m_step[curloc + 1] = (step_t) (step % STEP_T_MAX);
#endif

    return step;
}

#endif // BZET_IMPL_

#endif // BZET_H_
