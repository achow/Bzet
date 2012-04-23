#define NODE_ELS 8
#define BZET _Bzet8

#include "Bzet.h"

class Bzet8 {
    public:
        Bzet8() { m_bzet = BZET_FUNC(new)(); }
        Bzet8(uint64_t bit) { m_bzet = BZET_FUNC(new)(bit); }
        Bzet8(Bzet8& b) { m_bzet = BZET_FUNC(clone)(b.m_bzet); }
        ~Bzet8() { BZET_FUNC(destroy)(m_bzet); }

        bool operator==(Bzet8& left, Bzet8& right) {
            return BZET_FUNC(COMPARE)(left.m_bzet, right.m_bzet);
        }
        void operator=(Bzet8& left, Bzet8& right) {
            BZET_FUNC(setequal(left.m_bzet, right.m_bzet));
        }
        Bzet8 operator~(Bzet8& b) {
            Bzet8 n(b);
            BZET_FUNC(INVERT)(n.m_bzet);
            return n;
        }
        Bzet8 operator&(Bzet8& left, Bzet8& right) {
            Bzet8 result;
            BZET_FUNC(destroy)(result.m_bzet);
            result.m_bzet = BZET_FUNC(AND)(left.m_bzet, right.m_bzet);
        }
        Bzet8 operator^(Bzet8& left, Bzet8& right) {
            Bzet8 result;
            BZET_FUNC(destroy)(result.m_bzet);
            result.m_bzet = BZET_FUNC(XOR)(left.m_bzet, right.m_bzet);
        }

        Bzet8 operator|(Bzet8& left, Bzet8& right) {
            Bzet8 result;
            BZET_FUNC(destroy)(result.m_bzet);
            result.m_bzet = BZET_FUNC(OR)(left.m_bzet, right.m_bzet);
        }
        bool at(int64_t bit) { return BZET_FUNC(TEST)(m_bzet, bit); }
        void set(int64_t bit) { BZET_FUNC(SET)(m_bzet, bit); }
        void unset(int64_t bit) { BZET_FUNC(UNSET)(m_bzet, bit); }
        void setRange(int64_t start, int64_t len) { BZET_FUNC(RANGE)(m_bzet, start, len) }
        int64_t firstBit() { return BZET_FUNC(FIRST)(m_bzet); }
        int64_t lastBit() { return BZET_FUNC(LAST)(m_bzet); }
        int64_t count() { return BZET_FUNC(COUNT)(m_bzet); }
        int depth() { return BZET_FUNC(LEV)(m_bzet); }
        size_t size() { return BZET_FUNC(size)(m_bzet); }
        void printBzet() { BZET_FUNC(HEX)(); }
        bool empty() { BZET_FUNC(EMPTY)(m_bzet); }
        void clear() { BZET_FUNC(CLEAN)(m_bzet); }
        
        
    private:
        BZET_PTR m_bzet;
}
