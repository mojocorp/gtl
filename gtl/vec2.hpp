#pragma once

#include <gtl/gtl.hpp>

namespace gtl {
    /*!
    \class vec2 vec2.hpp geometry/vec2.hpp
    \brief 2 dimensional vector.
    \ingroup base

    This class is used by many other classes.

    \sa vec3
    */
    template <typename Type>
    class vec2 {
    public:
        Type x;
        Type y;

        //! The default constructor.The vector will be null.
        vec2()
        {
            x = 0;
            y = 0;
        }

        //! Constructs an instance with initial values from \a v.
        vec2(const Type v[2])
        {
            x = v[0];
            y = v[1];
        }

        //! Constructs an instance with the initial values from \a a_x and \a a_y.
        vec2(Type a_x, Type a_y)
        {
            x = a_x;
            y = a_y;
        }

        //! Constructs an instance with initial values from \a a_vec.
        vec2(const vec2<Type>& a_vec)
        {
            x = a_vec.x;
            y = a_vec.y;
        }

        //! Set new x and y coordinates for the vector. Returns reference to self.
        vec2<Type>& setValue(const Type v[2])
        {
            x = v[0];
            y = v[1];

            return *this;
        }

        //! Set new x and y coordinates for the vector. Returns reference to self.
        vec2<Type>& setValue(Type a_x, Type a_y)
        {
            x = a_x;
            y = a_y;

            return *this;
        }

        //! Returns a pointer to an array containing the coordinates of the vector.
        const Type* getValue() const
        {
            return &x;
        }

        //! Calculates and returns the dot product of this vector with \a a_vec.
        Type dot(const vec2<Type>& a_vec) const
        {
            return (x * a_vec.x + y * a_vec.y);
        }

        //! Return length of vector.
        Type length() const
        {
            return (Type)std::sqrt(x * x + y * y);
        }

        //! Return squared length of vector.
        Type sqrLength() const
        {
            return (x * x + y * y);
        }

        //! Normalize the vector to unit length. Return value is the original length of the vector before normalization.
        Type normalize()
        {
            const Type magnitude = length();

            if (magnitude != 0.0)
                (*this) *= (Type)(1.0 / magnitude);
            else
                setValue(0.0, 0.0);

            return magnitude;
        }

        //! Returns the normalized unit vector form of this vector.
        vec2<Type> normalized() const
        {
            vec2<Type> v(*this);
            v.normalize();

            return v;
        }

        //! Returns the cross product of this vector with \a a_vec.
        Type cross(const vec2<Type>& a_vec)
        {
            return (x * a_vec.y - y * a_vec.x);
        }

        //! Negate the vector (i.e. point it in the opposite direction).
        void negate()
        {
            x = -x;
            y = -y;
        }

        //! Index operator. Returns modifiable x or y value.
        Type& operator[](int i)
        {
            return (&x)[i];
        }

        //! Index operator. Returns x or y value.
        const Type& operator[](int i) const
        {
            return (&x)[i];
        }

        //! Multiply components of vector with value \a d. Returns reference to self.
        vec2<Type>& operator*=(const Type d)
        {
            x *= d;
            y *= d;

            return *this;
        }

        //! Divides components of vector with value \a d. Returns reference to self.
        vec2<Type>& operator/=(const Type d)
        {
            const Type inv = 1.0f / d;

            x *= inv;
            y *= inv;

            return *this;
        }

        //! Multiply components of vector with value \a a_vec.
        vec2<Type>& operator*=(const vec2<Type>& a_vec)
        {
            x *= a_vec.x;
            y *= a_vec.y;

            return *this;
        }

        //! Adds this vector and vector \a a_vec. Returns reference to self.
        vec2<Type>& operator+=(const vec2<Type>& a_vec)
        {
            x += a_vec.x;
            y += a_vec.y;

            return *this;
        }

        //! Subtracts vector \a a_vec from this vector. Returns reference to self.
        vec2<Type>& operator-=(const vec2<Type>& a_vec)
        {
            x -= a_vec.x;
            y -= a_vec.y;

            return *this;
        }

        //! Non-destructive negation operator.
        vec2<Type> operator-() const
        {
            return vec2<Type>(-x, -y);
        }

        friend vec2<Type> operator*(const vec2<Type>& a_vec, const Type d)
        {
            return vec2<Type>(a_vec.x * d, a_vec.y * d);
        }

        friend vec2<Type> operator*(const Type d, const vec2<Type>& a_vec)
        {
            return a_vec * d;
        }

        friend vec2<Type> operator/(const vec2<Type>& a_vec, const Type d)
        {
            return vec2<Type>(a_vec.x / d, a_vec.y / d);
        }

        friend vec2<Type> operator*(const vec2<Type>& v1, const vec2<Type>& v2)
        {
            return vec2<Type>(v1.x * v2.x, v1.y * v2.y);
        }

        friend vec2<Type> operator+(const vec2<Type>& v1, const vec2<Type>& v2)
        {
            return vec2<Type>(v1.x + v2.x, v1.y + v2.y);
        }

        friend vec2<Type> operator-(const vec2<Type>& v1, const vec2<Type>& v2)
        {
            return vec2<Type>(v1.x - v2.x, v1.y - v2.y);
        }

        //! Check the two given vector for equality.
        friend bool operator==(const vec2<Type>& v1, const vec2<Type>& v2)
        {
            return v1.x == v2.x && v1.y == v2.y;
        }

        //! Check the two given vector for inequality.
        friend bool operator!=(const vec2<Type>& v1, const vec2<Type>& v2)
        {
            return !(v1 == v2);
        }

        //! Check for equality with given tolerance.
        bool equals(const vec2<Type>& a_vec, const Type a_tolerance = 1E-2) const
        {
            return ((*this - a_vec).sqrLength() <= a_tolerance * a_tolerance);
        }

        friend std::ostream& operator<<(std::ostream& os, const vec2<Type>& vect)
        {
            return os << vect.x << " " << vect.y;
        }

        //! Largest representable vector
        static vec2<Type> max()
        {
            return vec2<Type>(std::numeric_limits<Type>::max(), std::numeric_limits<Type>::max());
        }
    };

    typedef vec2<int> vec2i;
    typedef vec2<float> vec2f;
    typedef vec2<double> vec2d;
} // namespace gtl
