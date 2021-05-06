#pragma once

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <vector>

namespace gtl {
#define SQR(x) pow((double)x, 2.0)

#if defined(min) || defined(max)
#error Error: min or max are defined as preprocessor macros, probably in <windows.h>.  Define NOMINMAX macro before including any system headers!
#if !defined(NOMINMAX)
#define NOMINMAX
#endif
#endif

// useful constants
#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif

#define EPS std::numeric_limits<Type>::epsilon()

    //! Convert a_value from degrees to radians.
    template <typename Type>
    Type DegToRad(Type a)
    {
        return (Type)(a * 0.01745329251994);
    }

    //!  Convert a_value from radians to degrees.
    template <typename Type>
    Type RadToDeg(Type a)
    {
        return (Type)(a * 57.2957795130823);
    }

    template <typename Type>
    inline const Type& min3(const Type& a, const Type& b, const Type& c)
    {
        return std::min(std::min(a, b), c);
    }

    template <typename Type>
    inline const Type& max3(const Type& a, const Type& b, const Type& c)
    {
        return std::max(std::max(a, b), c);
    }

    template <typename Type>
    constexpr Type clamp(const Type& value, const Type& min, const Type& max)
    {
        return (value < min) ? min : (value > max) ? max : value;
    }

    template <typename Type>
    inline Type rand(const Type min, const Type max)
    {
        return (std::rand() / (Type)RAND_MAX) * (max - min) + min;
    }

    template <typename Type>
    inline Type sign(Type a)
    {
        return (a >= (Type)0.0) ? (Type)1.0 : (Type)-1.0;
    }

    template <typename Type>
    inline bool equals(const Type& a, const Type& b, Type eps = std::numeric_limits<Type>::epsilon())
    {
        return (std::abs(a - b) < eps);
    }

    inline bool isPowerOfTwo(int value)
    {
        return (value & (value - 1)) == 0;
    }

    inline int nearestPowerOfTwo(int value)
    {
        int vlog = (int)(log((float)value) * 1.4427f);

        return (int)pow(2.f, (float)vlog);
    }

} // namespace gtl
