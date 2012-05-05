/*****************************************************************************
 * Generic code for Bzet8/Bzet16/Bzet32
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

		int64_t firstBit() const;
        int64_t lastBit() const;
        int64_t count() const;
        int64_t getBits(int64_t* bits, int64_t limit = 0, int64_t start = 0);
		bool empty() const;

    private:
#ifdef USE_LITERAL
        //BZET_PTR bitstobzet(void *data, size_t len);
        //void treetobits(unsigned char *buf, halfnode_t *node, int depth);
#endif

        NODETYPE _binop(BZET& left, BZET& right, OP op, int lev, size_t left_loc = 0, size_t right_loc = 0);
        void _printBzet(int stdOffset, FILE* target, int depth, size_t loc = 0, int offset = 0, bool pad = 0) const;
        int64_t _count(size_t loc, int depth) const;
        void subtree_not(size_t loc, int depth);
        static bool _seqset(BZET& b, size_t locb, BZET& right, size_t locright, int depth);

        size_t step_through(size_t loc, int depth) const;

        // Inline auxiliary functions
        static void display_error(const char* message, bool fatal = false, FILE* output = stderr);
        void init(size_t initial_alloc = INITIAL_ALLOC);
        void resize(size_t nhalfnodes);
        static size_t POW(int n);
        static int do_data_op(OP op, int left_data_bit, int right_data_bit);
        void append_subtree(BZET& src, size_t loc, int depth);
        void set_step(size_t loc, int depth);
        void mod_step(size_t loc, int depth, int diff);

        size_t m_nbufhalfnodes;
        size_t m_nhalfnodes;
        halfnode_t* m_bzet; //points to the bzet
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
    m_bzet = (halfnode_t *) malloc(initial_alloc * sizeof(halfnode_t));
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
    m_nhalfnodes = 0;
    m_nbufhalfnodes = INITIAL_ALLOC;
}

// Resizes the buffers in the Bzet if necessary
inline
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
    size_t dst_loc = m_nhalfnodes;

    // Resize to accomodate copy_size new elements
    resize(m_nhalfnodes + copy_size);

    // Do copy
    memcpy(m_bzet + dst_loc, src.m_bzet + loc, copy_size * sizeof(halfnode_t));
    memcpy(m_step + dst_loc, src.m_step + loc, copy_size * sizeof(step_t));
}

// Build step at loc
inline
void BZET::set_step(size_t loc, int depth) {
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
}

inline
void BZET::mod_step(size_t loc, int depth, int diff) {
    if (depth == 0)
        return;

    size_t step = m_step[loc] * STEP_T_MAX + m_step[loc + 1];
    step -= diff;

    m_step[loc] = (step_t) (step / STEP_T_MAX);
    m_step[loc + 1] = (step_t) (step % STEP_T_MAX);
}

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
}

// Destructor
BZET::~BZET() {
    free(m_bzet);
    free(m_step);
}

// -- OPERATIONS ---

// Assignment operator
BZET& BZET::operator=(const BZET& right) {
    // Resize left's bzet and step buffers
    resize(right.m_nhalfnodes);

    // Copy contents of bzet over
    m_depth = right.m_depth;

    memcpy(m_bzet, right.m_bzet, m_nhalfnodes * sizeof(halfnode_t));
    memcpy(m_step, right.m_step, m_nhalfnodes * sizeof(step_t));

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
    if (m_nhalfnodes == 0)
        return BZET(right);
    // If right bzet is empty, the bitwise OR will be equal to the left bzet
    else if (right.m_nhalfnodes == 0)
        return BZET(*this);
    // Otherwise just operate on them
    else {
        return binop(*((BZET *) this), right, OP_OR);
    }
}

// Bitwise AND
BZET BZET::operator&(BZET& right) {
    // If either bzet is empty, the bitwise AND will be an empty bzet
    if (m_nhalfnodes == 0 || right.m_nhalfnodes == 0)
        return BZET();
    // Otherwise just operate on them
    else {
        return binop(*((BZET *) this), right, OP_AND);
    }
}

// Bitwise XOR
BZET BZET::operator^(BZET& right) {
    // If left bzet is empty, the bitwise XOR will be equal to the right bzet
    if (m_nhalfnodes == 0)
        return BZET(right);
    // If right bzet is empty, the bitwise XOR will be equal to the left bzet
    else if (right.m_nhalfnodes == 0)
        return BZET(*this);
    // Otherwise just operate on them
    else {
        return binop(*((BZET *) this), right, OP_XOR);
    }
}

// Comparison operator
bool BZET::operator==(const BZET& right) const {
    if (m_depth != right.m_depth || m_nhalfnodes != right.m_nhalfnodes ||
        memcmp(m_bzet, right.m_bzet, m_nhalfnodes * sizeof(halfnode_t)))
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

    BZET temp(bit);
    return !(*((BZET*) this) & temp).empty();
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
    m_nhalfnodes = result.m_nhalfnodes;
    m_nbufhalfnodes = result.m_nbufhalfnodes;
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
}

// Bzet_LAST(b)
// TODO: Add support for bit literal subtrees
int64_t BZET::lastBit() const {
    if (empty())
        return -1;

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
    return (m_nhalfnodes * sizeof(halfnode_t) + 1);
}

// Prints out Bzet tree in human-readable form
void BZET::printBzet(int stdOffset, FILE* target) const {
    _printBzet(stdOffset, target, m_depth);
}

// Copy canonical form of bzet to supplied buffer
void BZET::hex(void *buf) const {
    unsigned char *charbuf = (unsigned char *) buf;
    charbuf[0] = m_depth;
    memcpy(charbuf + 1, m_bzet, m_nhalfnodes * sizeof(halfnode_t));
}

// Clears the bzet
void BZET::clear() {
    resize(0);
    m_depth = 0;
}

// Checks whether or not the bzet is empty
bool BZET::empty() const {
    return (m_nhalfnodes == 0);
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

    // Get corresponding nodes of the left and right tree
    halfnode_t c_left_data = left.m_bzet[left_loc];
    halfnode_t c_left_tree = left.m_bzet[left_loc + 1];
    halfnode_t c_right_data = right.m_bzet[right_loc];
    halfnode_t c_right_tree = right.m_bzet[right_loc + 1];
    left_loc += 2;
    right_loc += 2;

    // Create new node to hold data
    halfnode_t node_data = 0;
    halfnode_t node_tree = 0;

    // Reserve space for this node
    size_t loc = m_nhalfnodes;
    resize(m_nhalfnodes + 2);

    // For each element in the node
    for (int i = NODE_ELS - 1; i >= 0; i--) {
        int cur_left_tree_bit = (c_left_tree >> i) & 0x1;
        int cur_left_data_bit = (c_left_data >> i) & 0x1;
        int cur_right_tree_bit = (c_right_tree >> i) & 0x1;
        int cur_right_data_bit = (c_right_data >> i) & 0x1;
        
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
                    node_data <<= 1;
                    node_tree <<= 1;

                    break;

                // Delete right subtree, set data bit off
                case DB0:
                    // Skip the right subtree
                    right_loc = right.step_through(right_loc, lev - 1);

                    // Data bit is already 0
                    node_data <<= 1;
                    node_tree <<= 1;

                    break;

                // Delete left subtree, set data bit on
                case DA1: 
                    // Skip the left subtree
                    left_loc = left.step_through(left_loc, lev - 1);

                    // Turn on data bit
                    node_data = (node_data << 1) | 0x1;
                    node_tree <<= 1;

                    break;

                // Delete right subtree, set data bit on
                case DB1:
                    // Skip the right subtree
                    right_loc = right.step_through(right_loc, lev - 1);

                    // Turn on data bit
                    node_data = (node_data << 1) | 0x1;
                    node_tree <<= 1;

                    break;

                // Copy left subtree into result
                case CA: 
                    // Append left subtree
                    append_subtree(left, left_loc, lev - 1);

                    // Move through left subtree
                    left_loc = left.step_through(left_loc, lev - 1);

                    // Turn on tree bit
                    node_data <<= 1;
                    node_tree = (node_tree << 1) | 0x1;

                    break;

                // Copy right subtree into result
                case CB:
                    // Append right subtree
                    append_subtree(right, right_loc, lev - 1);

                    // Move through right subtree
                    right_loc = right.step_through(right_loc, lev - 1);

                    // Turn on tree bit
                    node_data <<= 1;
                    node_tree = (node_tree << 1) | 0x1;

                    break;

                // Copy left subtree into result and negate
                case NA:
                    end = m_nhalfnodes;

                    // Append left subtree
                    append_subtree(left, left_loc, lev - 1);

                    // Negate
                    subtree_not(end, lev - 1);

                    // Move through left subtree
                    left_loc = left.step_through(left_loc, lev - 1);

                    // Turn on tree bit
                    node_data <<= 1;
                    node_tree = (node_tree << 1) | 0x1;

                    break;

                // Copy right subtree into result and negate
                case NB:
                    end = m_nhalfnodes;

                    // Append right subtree
                    append_subtree(right, right_loc, lev - 1);

                    // Negate
                    subtree_not(end, lev - 1);

                    // Move through right subtree
                    right_loc = right.step_through(right_loc, lev - 1);

                    // Turn on tree bit
                    node_data <<= 1;
                    node_tree = (node_tree << 1) | 0x1;
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
            node_data = (node_data << 1) | (halfnode_t) do_data_op(op, cur_left_data_bit, cur_right_data_bit);
            node_tree <<= 1;
        }
#ifndef USE_LITERAL
        else {
            display_error("_binop: Invalid bit state (data bit and tree bit set simultaneously)", true);
        }
#endif
    }

    // Write nodes
    m_bzet[loc] = node_data;
    m_bzet[loc + 1] = node_tree;

    // Set step
    set_step(loc, lev);

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
        if (m_step[loc] != 1) {
            printf("data node with step != 1 at %d, nhalf=%d\n", loc, m_nhalfnodes);
            exit(1);
        }
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
        size_t old_size = b2.m_nhalfnodes;

        // Resize b2 to accommodate new heading nodes
        b2.resize(b2.m_nhalfnodes + diffdepth * 2);

        // Move bzet and step in b2 to accommodate new heading nodes
        memmove(b2.m_bzet + diffdepth * 2, b2.m_bzet, old_size * sizeof(halfnode_t));
        memmove(b2.m_step + diffdepth * 2, b2.m_step, old_size * sizeof(step_t));

        // Add new nodes and new step
        size_t loc = 0;
        for (int i = 0; i < diffdepth; i++) {
            b2.m_bzet[loc] = 0;
            b2.m_bzet[loc + 1] = (halfnode_t) ((halfnode_t) 0x1 << (NODE_ELS - 1));
            b2.set_step(loc, 1); // Cheat with this parameter since depth != 0

            loc += 2;
        }
    }
    // Otherwise b1 needs to be grown
    else {
        b1.m_depth = b2.m_depth;
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
            b1.set_step(loc, 1); // Cheat with this paramter since depth != 0

            loc += 2;
        }
    }
}

void BZET::normalize() {
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
        m_depth -= (step_t) leading;

        // Move bzet and step
        memmove(m_bzet, m_bzet + loc, (m_nhalfnodes - loc) * sizeof(halfnode_t));
        memmove(m_step, m_step + loc, (m_nhalfnodes - loc) * sizeof(step_t));

        // Resize
        resize(m_nhalfnodes - loc);
    }
}

size_t BZET::step_through(size_t loc, int depth) const {
    // Make sure loc is in range
    if (loc >= m_nhalfnodes) {
        printf("stepthrough fail trying %d, nhalf=%d\n", loc, m_nhalfnodes);
        display_error("", true);
        //return -1;
    }

    if (depth == 0)
        return loc + 1;

    // Get step offset
    size_t step_offset = ((size_t) m_step[loc]) * STEP_T_MAX + (size_t) m_step[loc + 1];

    // If the step offset is 0, the offset is too large to be actually stored
    // Compute it by examining the step of subtrees stemming from this node
    if (step_offset == 0) {
        // Get node tree bits
        halfnode_t tree_bits = m_bzet[loc + 1];

        // Advance to location of first subtree node
        loc += 2;

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
}

// TODO: use popcount, add support for bit literal subtree
int64_t BZET::_count(size_t loc, int depth) const {
    // Just return bit count with weight 1
    if (depth == 0) {
        halfnode_t data = m_bzet[loc];
        int64_t count = 0;
        for (int i = 0; i < NODE_ELS; i++)
            count += ((data >> i) & 1);
        return count;
    }

    // Handle all other levels
    halfnode_t data_node = m_bzet[loc];
    halfnode_t tree_node = m_bzet[loc + 1];
    int64_t count = 0;
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

    return count;
}

// Implementation of seqset
bool BZET::_seqset(BZET& b, size_t locb, BZET& right, size_t locright, int depth) {
    // At depth 0, we are at the lowest level, set the bit and finish
    if (depth == 0) {
        b.m_bzet[locb] |= right.m_bzet[locright];
        // If node is saturated, signal seqset to "handle" it
        if (b.m_bzet[locb] == (halfnode_t) -1) {
            return true;
        }

        return false;
    }

    halfnode_t& bdata_bits = b.m_bzet[locb];
    halfnode_t& btree_bits = b.m_bzet[locb + 1];
    //halfnode_t& rightdata_bits = right.m_bzet[locright];
    halfnode_t& righttree_bits = right.m_bzet[locright + 1];

    // Data bit higher up is set, we're done
    if ((bdata_bits | righttree_bits) == bdata_bits) {
        return false;
    }
    // Same tree bit is set, recurse
    else if ((btree_bits | righttree_bits) == btree_bits) {
        // Skip all subtrees in b before the subtree we want
        // TODO: use popcount to speed up?
        size_t loc = locb;
        locb += 2;
        bool first = true;
        for (int i = 0; i < NODE_ELS; i++)
            if ((btree_bits >> i) & 1) {
                if (!first)
                    locb = b.step_through(locb, depth - 1);
                else
                    first = false;
            }

        // Recurse
        bool ret = _seqset(b, locb, right, locright + 2, depth - 1);

        // Update step
        // Cheat with depth since depth != 0
        b.set_step(loc, 1);
        return ret;
    }
    // Tree bit is not on, set tree bit and append the subtree
    else {
        // Set tree bit
        btree_bits |= righttree_bits;
        // Append subtree
        b.append_subtree(right, locright + 2, depth - 1);
        // Update step
        // Cheat with depth since depth != 0
        b.set_step(locb, 1);
        return false;
    }
}

#endif // BZET_IMPL_

#endif // BZET_H_
