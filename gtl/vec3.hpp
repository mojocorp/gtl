#pragma once

#include <gtl/gtl.hpp>

namespace gtl {
    /*!
    \class vec3 vec3.hpp geometry/vec3.hpp
    \brief 3 dimensional vector.
    \ingroup base

    This class is used by many other classes.

    \sa vec2
    */
    template <typename Type>
    class vec3 {
    public:
        Type x;
        Type y;
        Type z;

        //! The default constructor.The vector will be null.
        vec3()
        {
            x = 0;
            y = 0;
            z = 0;
        }

        //! Constructs an instance with initial values from \a v.
        vec3(const Type v[3])
        {
            x = v[0];
            y = v[1];
            z = v[2];
        }

        //! Constructs an instance with the initial values from \a a_x, \a a_y and \a a_z.
        vec3(Type a_x, Type a_y, Type a_z)
        {
            x = a_x;
            y = a_y;
            z = a_z;
        }

        //! Constructs an instance with initial values from \a a_vec.
        vec3(const vec3<Type>& a_vec)
        {
            x = a_vec.x;
            y = a_vec.y;
            z = a_vec.z;
        }

        //! Set new x, y and z values for the vector. Returns reference to self.
        vec3<Type>& setValue(const Type v[3])
        {
            x = v[0];
            y = v[1];
            z = v[2];

            return *this;
        }

        //! Set new x, y and z values for the vector. Returns reference to self.
        vec3<Type>& setValue(Type a_x, Type a_y, Type a_z)
        {
            x = a_x;
            y = a_y;
            z = a_z;

            return *this;
        }

        //! Returns a pointer to an array containing the coordinates of the vector.
        const Type* getValue() const
        {
            return &x;
        }

        //! Calculates and returns the dot product of this vector with \a a_vec.
        Type dot(const vec3<Type>& a_vec) const
        {
            return (x * a_vec.x + y * a_vec.y + z * a_vec.z);
        }

        //! Return length of vector.
        Type length() const
        {
            return (Type)std::sqrt(x * x + y * y + z * z);
        }

        //! Return squared length of vector.
        Type sqrLength() const
        {
            return (x * x + y * y + z * z);
        }

        //! Normalize the vector to unit length. Return value is the original length of the vector before normalization.
        Type normalize()
        {
            const Type magnitude = length();

            if (magnitude != 0.0)
                (*this) *= (Type)(1.0 / magnitude);
            else
                setValue(0.0, 0.0, 0.0);

            return magnitude;
        }

        //! Returns the normalized unit vector form of this vector.
        vec3<Type> normalized() const
        {
            vec3<Type> v(*this);
            v.normalize();

            return v;
        }

        //! Returns the cross product of this vector with \a a_vec.
        vec3<Type> cross(const vec3<Type>& a_vec) const
        {
            return vec3<Type>(y * a_vec.z - a_vec.y * z,
                              z * a_vec.x - a_vec.z * x,
                              x * a_vec.y - a_vec.x * y);
        }

        //! Negate the vector (i.e. point it in the opposite direction).
        void negate()
        {
            x = -x;
            y = -y;
            z = -z;
        }

        //! Index operator. Returns modifiable x, y or z value.
        Type& operator[](int i) { return (&x)[i]; }

        //! Index operator. Returns x, y or z value.
        const Type& operator[](int i) const { return (&x)[i]; }

        //! Multiply components of vector with value \a d. Returns reference to self.
        vec3<Type>& operator*=(const Type d)
        {
            x *= d;
            y *= d;
            z *= d;

            return *this;
        }

        //! Divides components of vector with value \a d. Returns reference to self.
        vec3<Type>& operator/=(const Type d)
        {
            *this *= (1.0f / d);

            return *this;
        }

        //! Multiply components of vector with value \a a_vec.
        vec3<Type>& operator*=(const vec3<Type>& a_vec)
        {
            x *= a_vec.x;
            y *= a_vec.y;
            z *= a_vec.z;

            return *this;
        }

        //! Adds this vector and vector \a a_vec. Returns reference to self.
        vec3<Type>& operator+=(const vec3<Type>& a_vec)
        {
            x += a_vec.x;
            y += a_vec.y;
            z += a_vec.z;

            return *this;
        }

        //! Subtracts vector \a a_vec from this vector. Returns reference to self.
        vec3<Type>& operator-=(const vec3<Type>& a_vec)
        {
            x -= a_vec.x;
            y -= a_vec.y;
            z -= a_vec.z;

            return *this;
        }

        //! Non-destructive negation operator.
        vec3<Type> operator-() const
        {
            return vec3<Type>(-x, -y, -z);
        }

        friend vec3<Type> operator*(const vec3<Type>& a_vec, const Type d)
        {
            return vec3<Type>(a_vec.x * d, a_vec.y * d, a_vec.z * d);
        }

        friend vec3<Type> operator*(const Type d, const vec3<Type>& a_vec)
        {
            return a_vec * d;
        }

        friend vec3<Type> operator/(const vec3<Type>& a_vec, const Type d)
        {
            return vec3<Type>(a_vec.x / d, a_vec.y / d, a_vec.z / d);
        }

        friend vec3<Type> operator*(const vec3<Type>& v1, const vec3<Type>& v2)
        {
            return vec3<Type>(v1.x * v2.x, v1.y * v2.y, v1.z * v2.z);
        }

        friend inline vec3<Type> operator+(const vec3<Type>& v1, const vec3<Type>& v2)
        {
            return vec3<Type>(v1.x + v2.x, v1.y + v2.y, v1.z + v2.z);
        }

        friend inline vec3<Type> operator-(const vec3<Type>& v1, const vec3<Type>& v2)
        {
            return vec3<Type>(v1.x - v2.x, v1.y - v2.y, v1.z - v2.z);
        }

        //! Check the two given vector for equality.
        friend bool operator==(const vec3<Type>& v1, const vec3<Type>& v2)
        {
            return (v1.x == v2.x && v1.y == v2.y && v1.z == v2.z);
        }

        //! Check the two given vector for inequality.
        friend bool operator!=(const vec3<Type>& v1, const vec3<Type>& v2)
        {
            return !(v1 == v2);
        }

        //! Check for equality with given tolerance.
        bool equals(const vec3<Type>& a_vec, const Type a_tolerance = 1E-2) const
        {
            return ((*this - a_vec).sqrLength() <= a_tolerance * a_tolerance);
        }

        friend std::ostream& operator<<(std::ostream& os, const vec3<Type>& vect)
        {
            return os << vect.x << " " << vect.y << " " << vect.z;
        }

        //! Largest representable vector
        static vec3<Type> max()
        {
            return vec3<Type>(std::numeric_limits<Type>::max(),
                              std::numeric_limits<Type>::max(),
                              std::numeric_limits<Type>::max());
        }
    };

    typedef vec3<int> vec3i;
    typedef vec3<float> vec3f;
    typedef vec3<double> vec3d;
} // namespace gtl
