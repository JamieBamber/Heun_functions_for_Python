/*
C++ sgn function
*/

#ifndef SGN_HPP_
#define SGN_HPP_

namespace MyMath{
template <typename T> int sgn(T val) {
    return (T(0) < val) - (val < T(0));
}
}

#endif /* SGN_HPP_ */
