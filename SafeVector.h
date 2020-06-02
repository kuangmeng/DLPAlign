/////////////////////////////////////////////////////////////////
// SafeVector.h
//
// STL vector with array bounds checking.  To enable bounds
// checking, #define ENABLE_CHECKS.
/////////////////////////////////////////////////////////////////

#ifndef SAFEVECTOR_H
#define SAFEVECTOR_H

#include <cassert>
#include <vector>
using namespace std;

/////////////////////////////////////////////////////////////////
// SafeVector
//
// Class derived from the STL std::vector for bounds checking.
/////////////////////////////////////////////////////////////////

template<class TYPE>
class SafeVector: public std::vector<TYPE>
{
public:

    // miscellaneous constructors
    SafeVector() :
        std::vector<TYPE>()
    {
    }
    SafeVector(size_t size) :
        std::vector<TYPE>(size)
    {
    }
    SafeVector(size_t size, const TYPE &value) :
        std::vector<TYPE>(size, value)
    {
    }
    SafeVector(const SafeVector &source) :
        std::vector<TYPE>(source)
    {
    }

};

// some commonly used vector types
typedef SafeVector<int> VI;
typedef SafeVector<VI> VVI;
typedef SafeVector<VVI> VVVI;
typedef SafeVector<float> VF;
typedef SafeVector<VF> VVF;
typedef SafeVector<VVF> VVVF;

#endif
