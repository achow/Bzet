#include "Bzet4.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>

#define NODE_ELS 4
#define HEADER_SIZE 1
#ifndef INITIAL_ALLOC
#define INITIAL_ALLOC 1024
#endif
#define RESIZE_SCALE 2

const int Bzet4::powersof4[10] = { 1, 4, 16, 64, 256, 1024, 4096, 16384, 65536, 262144 };

const ACTION Bzet4::optable[64] = {       
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

#if (defined _DEBUG || defined DEBUG)
void Bzet4::validateBzet(size_t loc, int lev) {
    if (!lev)
        lev = m_bzet[0];

    if (lev == 1)
        return;

    size_t nextLoc = loc;
    for (int i = 3; i >= 0; --i) {
        if ((m_bzet[loc] >> i) & 1) {
            bool check = nextLoc + 1 >= 0 && nextLoc < m_size;
            if (!check) {
                printf("assert failed at nextLoc = %d, size = %d\n", (int) nextLoc, (int) m_size);
                printf("bytes are: ");
                printBytes();
            }
            assert(nextLoc + 1 >= 0 && nextLoc < m_size);
            validateBzet(nextLoc + 1, lev - 1);
            nextLoc = stepThrough(nextLoc + 1) - 1;
        }
    }

    bool check = (dust(m_bzet[loc]) == m_bzet[loc]);
    if (!check)
        printf("assert failed at loc = %d, lev = %d\n", (int) loc, lev);
    assert(check);
    /*
    for (int i = 1; i < m_size; ++i) {
        bool test = (m_step[i] > 0 && m_step[i] <= m_size);
        if (!test) {
            printf("assertion failed, i = %d, m_step[i] = %d\n", i, m_step[i]);
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
    printf("\nSize is %d\n", (int) m_size);
    printf("STEP: ");
    int x = 0;
    for (int i = 1; i < m_size; ++i) {
        if (x % 10 == 0)
            printf("\n");
        printf("%5d", (int) m_step[i]);
        if (m_step[i] == 2 && m_step[i + 1] != 1) {
            printf("\n | invalid step at %d | \n", i);
        }
        x++;
    }
    printf("\n");
}
#endif

/*****************************************************************************
 * 
 *  Function name: Bzet4
 *
 *  Purpose:       Default constructor for Bzet4
 *                 Creates a Bzet4 with depth 0
 *
 *  Inputs:        None
 *  Return values: None
 * 
 *  Author:        Alex Chow
 *  Date:          10/23/2011
 *
 *****************************************************************************/
Bzet4::Bzet4() { 
    init();
    clear();
}

/*****************************************************************************
 * 
 *  Function name: Bzet4
 *
 *  Purpose:       Copy constructor for Bzet4
 *
 *  Inputs:        Bzet4& src: source Bzet4 to copy from
 *  Return values: None
 * 
 *  Author:        Alex Chow
 *  Date:          10/25/2011
 *
 *****************************************************************************/
Bzet4::Bzet4(const Bzet4& src) { 
    init();

    //resize and copy from src
    resize(src.size());
    src.hex(m_bzet);

    memcpy(m_step, src.m_step, m_size * sizeof(*m_step));
}

/*****************************************************************************
 * 
 *  Function name: Bzet4
 *
 *  Purpose:       Constructor for Bzet4. Initializes a Bzet4 with bit bit 
 *                 turned on
 *
 *  Inputs:        int64_t bit: What bit to turn on
 *  Return values: None
 * 
 *  Author:        Alex Chow
 *  Date:          10/23/2011
 *
 *****************************************************************************/
Bzet4::Bzet4(int64_t bit) {
    init();

    if (bit < 0)
        display_error("Bzet4::Bzet4(bit): invalid bit to set", true);

    //build depth
    int depth = 1;

    int64_t tempbit = bit;
    //level 1 holds 16 bits, numbers 0 to 15
    while (tempbit > 15) {
        tempbit /= (int64_t) NODE_ELS;
        ++depth;
    }

    //resize to accomodate full bzet
    //+1 since level 1 requires 2 bytes to store data literals
    resize(HEADER_SIZE + depth + 1);

    //set depth bit
    m_bzet[0] = (unsigned char) depth;

    //set data nodes at level 1
    //set bits containing level 0 nodes 0 and 1 (bits 0-7)
    m_bzet[depth] = (unsigned char) 0x80 >> (bit & 0xF);
    //set bits containing level 0 nodes 2 and 3 (bits 8-15) if previous node not set
    m_bzet[depth + 1] = (m_bzet[depth]) ? 0x00 : (unsigned char) 0x80 >> (bit & 0x7);

    //since level 1 set, divide by 16 to see if tree nodes needed
    bit >>= 4; 

    //set tree nodes if needed
    int loc = depth;
    while (bit > 0 && loc > 1) {
        m_bzet[--loc] = (unsigned char) 0x08 >> (bit & 0x3);
        bit >>= 2; //divide by 4
    }

    //set m_step
    //no need to check m_size <= 255 to make sure m_step values doesn't overflow
    //since it would never happen (4^255 ~ 10^153)
    for (int i = HEADER_SIZE; i < m_size; ++i)
        m_step[i] = (unsigned char) m_size - i;

#if (defined _DEBUG || defined DEBUG)
    validateBzet();
#endif
}

/*****************************************************************************
 * 
 *  Function name: Bzet4
 *
 *  Purpose:       Constructor for Bzet4. Initializes a Bzet4 with len bits 
 *                 turned on, starting from startbit.
 *
 *  Inputs:        int64_t startbit: Starting bit to turn on
 *                 int64_t len: Number of bits to turn on, starting from startbit
 *  Return values: None
 * 
 *  Author:        Alex Chow
 *  Date:          11/3/2011
 *
 *****************************************************************************/
Bzet4::Bzet4(int64_t startbit, int64_t len) {
    //intialize Bzet
    init();
    clear();

    //set the range
    setRange(startbit, len);

#if (defined _DEBUG || defined DEBUG)
    validateBzet();
#endif
}

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
 *****************************************************************************/
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
}*/

/*****************************************************************************
 * 
 *  Function name: Bzet4
 *
 *  Purpose:       Constructor for Bzet4. Allows for creation of a Bzet from
 *                 arbitrary input
 *
 *  Inputs:        void* data: Pointer to data to compress
 *                 int64_t size: Size of data in bytes
 *  Return values: None
 * 
 *  Author:        Alex Chow
 *  Date:          11/4/2011
 *
 *****************************************************************************/
Bzet4::Bzet4(void* data, int size) {
    //initialize Bzet
    init();
    clear();

    //make sure pointer is non-null
    if (data) {
        //retype data to pointer to bytes
        unsigned char* pData = (unsigned char*) data;
        int64_t bitno = 0;
        for (int i = 0; i < size; ++i) {
            unsigned char c = pData[i];
            //set each bit if it's on
            for (int j = 7; j >= 0; --j) {
                if ((c >> j) & 1)
                    set(bitno);
                ++bitno;
            }
        }
    } else {
        display_error("Bzet4::Bzet4(void*, int): null pointer, empty Bzet4 created instead");
    }

#if (defined _DEBUG || defined DEBUG)
    validateBzet();
#endif
}

/*****************************************************************************
 * 
 *  Function name: ~Bzet4
 *
 *  Purpose:       Destructor
 *
 *  Inputs:        None
 *  Return values: None
 * 
 *  Author:        Alex Chow
 *  Date:          10/24/2011
 *
 *****************************************************************************/
Bzet4::~Bzet4() {
    free(m_bzet);
    free(m_step);
}

/*****************************************************************************
 * 
 *  Function name: operator=
 *
 *  Purpose:       Assignment operator
 *
 *  Inputs:        Bzet4& right: Right hand side of = operation
 *  Return values: *this
 * 
 *  Author:        Alex Chow
 *  Date:          10/24/2011
 *
 *****************************************************************************/
Bzet4& Bzet4::operator=(const Bzet4& right) {
    //resize and copy from right
    resize(right.size());
    right.hex(m_bzet);
    memcpy(m_step, right.m_step, m_size * sizeof(*m_step));

    return *this;
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
 *****************************************************************************/
Bzet4 Bzet4::operator~() const {
    //clone bzet and NOT in place
    Bzet4 not(*this);
    not.subtreeNot(1, (int) m_bzet[0]);

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
 *****************************************************************************/
 void Bzet4::appendSubtree(const Bzet4& src, size_t loc) {
     size_t step = src.stepThrough(loc);

     //get length to copy
     size_t len = src.stepThrough(loc) - loc;
     size_t old_size = m_size;

     //resize bzet to accommodate new nodes
     resize(m_size + len);

     //copy bzet
     memcpy(m_bzet + old_size, src.m_bzet + loc, len);
     
     //copy m_step
     memcpy(m_step + old_size, src.m_step + loc, len);
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
 *****************************************************************************/
void Bzet4::dropNodes(size_t loc, int n) {
#if (defined _DEBUG || defined DEBUG)
    assert(loc > 0 && loc < m_size && n >= 0);
#endif

    memmove(m_bzet + loc, m_bzet + loc + n, m_size - (loc + n));
    memmove(m_step + loc, m_step + loc + n, m_size - (loc + n));

    resize(m_size - n);

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
 *****************************************************************************/
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
 *****************************************************************************/
NODETYPE Bzet4::_binop(const Bzet4& left, const Bzet4& right, OP op, int lev, size_t left_loc, size_t right_loc, size_t loc) {
    //left_loc and right_loc are unmodified until an operation is done
    //so they both point to the corresponding nodes in left and right used to build the current node
    unsigned char c_left = left.m_bzet[left_loc];
    unsigned char c_right = right.m_bzet[right_loc];

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
        c_left = left.m_bzet[left_loc + 1];
        c_right = right.m_bzet[right_loc + 1];

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
            if (m_size == 1) {
                resize(2);
                m_bzet[0]++;
                m_bzet[1] = 0x80;
                m_step[1] = 1;
                return NORMAL;
            }

            return SATURATED;
        }
        //empty node
        else if (!result_c1 && !result_c2) {
            //build bzet manually if bzet is empty
            if (m_size == 1) {
                clear();
                return NORMAL;
            }

            return EMPTY;
        }

        //otherwise commit data
        //create two nodes to accommodate data nodes
        resize(m_size + 2);

        //write result
        m_bzet[loc] = result_c1;
        m_bzet[loc + 1] = result_c2;

        //set steps
        m_step[loc] = 2;
        m_step[loc + 1] = 1;
        
        return NORMAL;
    }

    //create new node to hold data
    resize(m_size + 1);
    m_bzet[loc] = 0x00;

    //go through each element in the current node
    for (int i = NODE_ELS - 1; i >= 0; --i) {
        int cur_left_tree_bit = (c_left >> i) & 0x1;
        int cur_left_data_bit = (c_left >> (i + NODE_ELS)) & 0x1;
        int cur_right_tree_bit = (c_right >> i) & 0x1;
        int cur_right_data_bit = (c_right >> (i + NODE_ELS)) & 0x1;
 
        //TT: if both tree bits are on
        if (cur_left_tree_bit && cur_right_tree_bit) {
            //turn on tree bit
            //m_bzet[loc] |= 0x01 << i;

            //recurse
            NODETYPE cn = _binop(left, right, op, lev - 1, left_loc + 1, right_loc + 1, m_size);

            if (cn == SATURATED) {
                //saturated subtree, turn on data bit
                m_bzet[loc] |= 0x10 << i;
            }
            else if (cn == EMPTY) {
                //do nothing, elements are already set to 0
            }
            else {
                //subtree exists, turn on tree bit
                m_bzet[loc] |= 0x01 << i;
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
                    m_bzet[loc] |= 0x80 >> ((NODE_ELS - 1) - i);
                    break;

                //delete right subtree, set data bit on
                case DB1:
                    //skip the right subtree
                    right_loc = right.stepThrough(right_loc + 1) - 1;
                    //turn on data bit in the result node
                    m_bzet[loc] |= 0x80 >> ((NODE_ELS - 1) - i);
                    break;

                //copy left subtree into result
                case CA: 
                    //turn on tree bit
                    m_bzet[loc] |= 0x08 >> ((NODE_ELS - 1) - i);

                    //append left subtree
                    appendSubtree(left, left_loc + 1);

                    //move through left subtree
                    left_loc = left.stepThrough(left_loc + 1) - 1;
                    break;

                //copy right subtree into result
                case CB:
                    //turn on tree bit
                    m_bzet[loc] |= 0x08 >> ((NODE_ELS - 1) - i);

                    //append right subtree
                    appendSubtree(right, right_loc + 1);

                    //move through right subtree
                    right_loc = right.stepThrough(right_loc + 1) - 1;
                    break;

                //copy left subtree into result and negate
                case NA:
                    //turn on tree bit
                    m_bzet[loc] |= 0x08 >> ((NODE_ELS - 1) - i);

                    end = m_size;

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
                    m_bzet[loc] |= 0x08 >> ((NODE_ELS - 1) - i);

                    end = m_size;

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
                    m_bzet[loc] |= 0x10 << i;
                }
            }
        }
        //malformed bzet
        else {
            display_error("Bzet4::_binop: Malformed bzet (data bit and tree bit are on simultaneously)", true);
        }
    }

    //set step
    if (m_size - loc > 255)
        m_step[loc] = 0;
    else 
        m_step[loc] = (unsigned char) (m_size - loc);

    //resulting node is empty
    if (m_bzet[loc] == 0x00) {
        //build bzet manually if bzet is empty
        if (m_size == 2) {
            clear();
            return NORMAL;
        }

        resize(m_size - 1);
        return EMPTY;
    }
    //resulting node is saturated
    else if (m_bzet[loc] == 0xF0) {
        if (m_size == 2) {
            resize(2);
            m_bzet[0]++;
            m_bzet[1] = 0x80;
            m_step[1] = 1;
            return NORMAL;
        }

        resize(m_size - 1);
        return SATURATED;
    }

#if (defined _DEBUG || defined DEBUG)
    assert(dust(m_bzet[loc]) == m_bzet[loc]);
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
 *****************************************************************************/
Bzet4 Bzet4::binop(Bzet4& left, Bzet4& right, OP op) {
    //align bzets
    align(left, right);

    Bzet4 result;

    //set level byte
    result.m_bzet[0] = left.m_bzet[0];

    //operate
    result._binop(left, right, op, left.m_bzet[0]);

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
 *****************************************************************************/
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
 *****************************************************************************/
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
 *****************************************************************************/
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
 *****************************************************************************/
bool Bzet4::operator==(const Bzet4& right) const {
    //if sizes are the same, check content
    if (m_size == right.m_size) {
        if (!memcmp(m_bzet, right.m_bzet, m_size))
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
 *****************************************************************************/
int64_t Bzet4::firstBit() const {
    //empty Bzet
    if (m_size == 1)
        return -1;

    int level = m_bzet[0]; //current level in tree
    int64_t bitBase = 0; //number of bits prior to the current node
    size_t loc = 1; //location in Bzet

    //just traverse the leftmost branch
    while (level) {
        //special calculation at level 1
        if (level == 1) {
            //check data nodes 0 and 1
            unsigned char c = m_bzet[loc];

            //if any bit is set in node 0 and 1, first bit found
            if (c)
                for (int i = 7; i >= 0; --i)
                    if ((c >> i) & 1)
                        return bitBase + 7 - (int64_t) i;

            //if no data bit set in first 2 nodes, check next 2 nodes
            c = m_bzet[loc + 1];
            for (int i = 7; i >= 0; --i)
                if ((c >> i) & 1)
                    return bitBase + 15 - i;
#ifdef DEBUG
            assert(c);
#endif
        }

        //calculation for all other nodes
        unsigned char data_bits = (m_bzet[loc] >> NODE_ELS) & 0xF;
        unsigned char tree_bits = m_bzet[loc] & 0xF;

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
 *****************************************************************************/
int64_t Bzet4::lastBit() const {
    //empty Bzet
    if (size() == 1)
        return -1;

    size_t loc = 1; //location in m_bzet (node)
    int level = m_bzet[0]; //current level in tree, initialize to bzet depth
    int64_t bitBase = 0; //number of bits prior to the current node

    do {
        //special calculation at level 1
        if (level == 1) {
            //check data nodes 2 and 3
            unsigned char c = m_bzet[loc + 1];

            if (c)
                for (int i = 0; i < 8; ++i)
                    if ((c >> i) & 1)
                        return bitBase + 15 - i;

            //if no data bit set in last 2 nodes, check first 2
            c = m_bzet[loc];
            for (int i = 0; i < 8; ++i)
                if ((c >> i) & 1)
                    return bitBase + 7 - i;
#if (defined _DEBUG || defined DEBUG)
            assert(c);
#endif
        }

        unsigned char c = m_bzet[loc];
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
    } while (loc < m_size);

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
 *****************************************************************************/
int64_t Bzet4::count() const {
    //empty bzet has 0
    if (empty())
        return 0;

    int64_t bitCount = 0;

    //go through each node
    for (size_t i = HEADER_SIZE; i < m_size; ++i) {
        int lev = depthAt(i);

        //if lev is 1, all data bits have weight 1, do bit count on all 8 bits
        if (lev == 1) {
            unsigned char c = m_bzet[i];
            int num_bits = 0;
            for (int j = 0; j < 8; ++j) 
                num_bits += ((c >> j) & 1);
            bitCount += num_bits;
        }
        //otherwise do weighted bit count on data bits only
        else {
            unsigned char data_bits = (m_bzet[i] >> NODE_ELS) & 0xF;

            //count the number of data bits on
            int num_bits = 0;
            for (int j = 0; j < NODE_ELS; ++j) 
                num_bits += (data_bits >> j) & 1;

            //add weighted bit count to running total
            bitCount += num_bits * pow4(lev);
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
 *****************************************************************************/
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
 *****************************************************************************/
bool Bzet4::empty() const {
    return (m_size == 1);
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
 *****************************************************************************/
bool Bzet4::at(int64_t bit) const {
    //get number of bits this Bzet stores
    int64_t totalBits = pow4(m_bzet[0] + 1);

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
 *****************************************************************************/
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
 *****************************************************************************/
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
 *****************************************************************************/
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
 *****************************************************************************/
int Bzet4::depth() const {
    return m_bzet[0];
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
 *****************************************************************************/
size_t Bzet4::size() const {
    return m_size;
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
 *****************************************************************************/
void Bzet4::clear() {
    resize(1);
    m_bzet[0] = 0x00; //set level byte to 0
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
 *****************************************************************************/
void Bzet4::hex(void* str) const {
    memcpy(str, m_bzet, m_size); 
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
 *****************************************************************************/
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
 *****************************************************************************/
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
 *****************************************************************************/
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
 *****************************************************************************/
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
 *****************************************************************************/
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
        b2.m_bzet[0] = b2.depth() + diffDepth;

        //shift Bzet to make room for extra nodes
        for (size_t i = b2.size() - 1; i >= HEADER_SIZE + diffDepth; --i) {
            b2.m_bzet[i] = b2.m_bzet[i - diffDepth];
            b2.m_step[i] = b2.m_step[i - diffDepth];
        }

        //fill in extra nodes
        for (int i = 0; i < diffDepth; ++i) {
            b2.m_bzet[HEADER_SIZE + i] = 0x08;

            //set m_step
            if ((size_t) b2.m_step[diffDepth + 1] + (diffDepth - i) >= 255)
                b2.m_step[HEADER_SIZE + i] = 0;
            else
                b2.m_step[HEADER_SIZE + i] = b2.m_step[diffDepth + 1] + (diffDepth - i);
        }
    } 
    //b2 is bigger
    else {
        //change sign of diffDepth
        diffDepth = -diffDepth;

        //resize b1
        b1.resize(b1.size() + diffDepth);

        //modify b1 depth
        b1.m_bzet[0] = b1.depth() + diffDepth;

        //shift Bzet to make room for extra nodes
        for (size_t i = b1.size() - 1; i >= HEADER_SIZE + diffDepth; --i) {
            b1.m_bzet[i] = b1.m_bzet[i - diffDepth];
            b1.m_step[i] = b1.m_step[i - diffDepth];
        }

        //fill in extra nodes
        for (int i = 0; i < diffDepth; ++i) {
            b1.m_bzet[HEADER_SIZE + i] = 0x08;

            //set m_step
            if ((size_t) b1.m_step[diffDepth + 1] + (diffDepth - i) >= 255)
                b1.m_step[HEADER_SIZE + i] = 0;
            else
                b1.m_step[HEADER_SIZE + i] = b1.m_step[diffDepth + 1] + (diffDepth - i);
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
 *****************************************************************************/
int Bzet4::depthAt(size_t tLoc) const {
    //check if out of range
    if (tLoc > size())
        return -1;

    //if loc is 1, just return bzet depth
    if (tLoc == 1)
        return depth();

    //the node at m_bzet[1] is at depth m_bzet[0]
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
 *  Return values: -1 if loc is out of valid range (1 to m_size - 1)
 *                 otherwise the location of the next node after walking
 *                 through the subtree starting at loc
 * 
 *  Author:        Alex Chow
 *  Date:          10/28/2011
 *
 *****************************************************************************/
size_t Bzet4::stepThrough(size_t loc) const {
    //make sure loc is in range
    if (loc <= 0 || loc >= m_size)
        return -1;

    //retrieve the offset of the end of the subtree at loc
    unsigned char loc_offset = m_step[loc];

    //if the offset is 0, the offset is too large to be actually stored
    //compute it by using the subtrees stemming from this node
    if (loc_offset == 0) {
        //get node tree bits
        unsigned char tree_bits = m_bzet[loc] & 0xF;

        //advance to location of first subtree node
        loc++;

        //step through each subtree
        for (int i = NODE_ELS - 1; i >= 0; i--) {
            if ((tree_bits >> i) & 1)
                loc = stepThrough(loc);
        }

        return loc;
    }

    return loc + m_step[loc];
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
 *****************************************************************************/
void Bzet4::subtreeNot(size_t loc, int depth) {
    if (loc == 0 || loc >= m_size)
        return;

    if (!depth)
        depth = depthAt(loc);

    //if no tree bits are on, just not the data bits
    if (!(m_bzet[loc] & 0xF) && depth > 1) {
        m_bzet[loc] = ~(m_bzet[loc] | 0xF);
        return;
    }

    //if depth is 1, we just not this and the next node
    if (depth == 1) {
        m_bzet[loc] = ~m_bzet[loc];
        m_bzet[loc + 1] = ~m_bzet[loc + 1];
        return;
    }

    size_t nextLoc = loc + 1;

    //decrement depth
    --depth;

    for (int i = NODE_ELS - 1; i >= 0; --i) {
        //if tree bit is on, recursively not each node
        if ((m_bzet[loc] >> i) & 1) {
            subtreeNot(nextLoc, depth);
            nextLoc = stepThrough(nextLoc);
        }
    }

    //get the current byte
    unsigned char c = m_bzet[loc]; 

    //new data bits
    unsigned char data_bits = (~(c >> NODE_ELS)) << NODE_ELS;

    //add tree bits
    unsigned char newc = data_bits | (c & 0xF);

    //set up for dusting with tree bits priority
    unsigned char tree_bits = c & 0xF;

    //dust using tree bits and write back
    m_bzet[loc] = newc & (~(tree_bits << NODE_ELS) | tree_bits);
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
 *****************************************************************************/
void Bzet4::_printBzet(int stdOffset, FILE* target, size_t loc, int depth, int offset, bool pad) const {
    //print level info
    if (loc == 1) {
        fprintf(target, "%.2XL", m_bzet[0]);

        //quick fix for empty bzet
        if (empty())
            fprintf(target, "\n");
    }

    if (!depth)
        depth = depthAt(loc);

    //no reading past the bzet!
    if (loc >= size())
        return;

    unsigned char c = m_bzet[loc];
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
        fprintf(target, "D(%.2X%.2X)\n", m_bzet[loc], m_bzet[loc + 1]);
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
 *  Purpose:       Resizes m_bzet to the specified size in bytes
 *
 *  Inputs:        int size: New size of m_bzet
 *  Return values: None
 * 
 *  Author:        Alex Chow
 *  Date:          10/23/2011
 *
 *****************************************************************************/
void Bzet4::resize(size_t size) {
    //if reallocation is required
    if (size > m_bufsize) {
        //get new size required by growing it by RESIZE_SCALE repeatedly
        while (size > m_bufsize)
            m_bufsize *= RESIZE_SCALE;

        //reallocate
        m_bzet = (unsigned char*) realloc(m_bzet, m_bufsize * sizeof(unsigned char));
        m_step = (unsigned char*) realloc(m_step, m_bufsize * sizeof(unsigned char));

        //check that realloc succeeded
        if (!m_bzet || !m_step) {
            fprintf(stderr, "Fatal error: Resizing bzet failed attempting to allocate %d bytes\n", (int) size);
            display_error("", true);
        }
    }

    //update internal size
    m_size = size;
}

/*****************************************************************************
 * 
 *  Function name: normalize
 *
 *  Purpose:       Normalizes m_bzet to get rid of extra nodes
 *
 *  Inputs:        None
 *  Return values: None
 * 
 *  Author:        Alex Chow
 *  Date:          10/26/2011
 *
 *****************************************************************************/
void Bzet4::normalize() {
    //no need to normalize bzet with only 1 level
    if (m_bzet[0] == 0x01)
        return;

    //count leading 0x08 nodes
    int i = 1;
    while (m_bzet[i] == 0x08)
        ++i;
    --i; //i is now the number of leading 0x08 nodes

    //if there are leading 0x08 nodes, get rid of them
    if (i) {
        //drop i nodes from the beginning
        dropNodes(1, i);

        //change depth accordingly
        m_bzet[0] -= i;
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
 *****************************************************************************/
void Bzet4::loadBzet(void* bzet_literal, int size) {
    if (bzet_literal) {
        resize(size);
        memcpy(m_bzet, bzet_literal, size);
    } else {
        display_error("Bzet4::loadBzet(void*, int64_t): null pointer, doing nothing");
    }
}
