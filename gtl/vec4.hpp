#pragma once

#include <gtl/gtl.hpp>

namespace gtl {
    /*!
    \class vec4 vec4.hpp geometry/vec4.hpp
    \brief 4 dimensional vector.
    \ingroup base

    This class is used by many other classes.

    \sa vec3
    */
    template <typename Type>
    class vec4 {
    public:
        Type x;
        Type y;
        Type z;
        Type w;

        //! The default constructor.The vector will be null.
        vec4()
        {
            x = 0;
            y = 0;
            z = 0;
            w = 0;
        }

        //! Constructs an instance with initial values from \a v.
        vec4(const Type v[4])
        {
            x = v[0];
            y = v[1];
            z = v[2];
            w = v[3];
        }

        //! Constructs an instance with the initial values from \a a_x, \a a_y, \a a_z and \a a_w.
        vec4(Type a_x, Type a_y, Type a_z, Type a_w)
        {
            x = a_x;
            y = a_y;
            z = a_z;
            w = a_w;
        }

        //! Constructs an instance with initial values from \a a_vec.
        vec4(const vec4<Type>& a_vec)
        {
            x = a_vec.x;
            y = a_vec.y;
            z = a_vec.z;
            w = a_vec.w;
        }

        //! Set new x, y, z and w values for the vector. Returns reference to self.
        vec4<Type>& setValue(const Type v[4])
        {
            x = v[0];
            y = v[1];
            z = v[2];
            w = v[3];

            return *this;
        }

        //! Set new x, y, z and w values for the vector. Returns reference to self.
        vec4<Type>& setValue(Type a_x, Type a_y, Type a_z, Type a_w)
        {
            x = a_x;
            y = a_y;
            z = a_z;
            w = a_w;

            return *this;
        }

        //! Returns a pointer to an array containing the coordinates of the vector.
        const Type* getValue() const
        {
            return &x;
        }

        //! Calculates and returns the dot product of this vector with \a a_vec.
        Type dot(const vec4<Type>& a_vec) const
        {
            return (x * a_vec.x + y * a_vec.y + z * a_vec.z + w * a_vec.w);
        }

        //! Return length of vector.
        Type length() const
        {
            return (Type)std::sqrt(x * x + y * y + z * z + w * w);
        }

        //! Return squared length of vector.
        Type sqrLength() const
        {
            return (x * x + y * y + z * z + w * w);
        }

        //! Normalize the vector to unit length. Return value is the original length of the vector before normalization.
        Type normalize()
        {
            const Type magnitude = length();

            if (magnitude != 0.0)
                (*this) *= (Type)(1.0 / magnitude);
            else
                setValue(0.0, 0.0, 0.0, 0.0);

            return magnitude;
        }

        //! Returns the normalized unit vector form of this vector.
        vec4<Type> normalized() const
        {
            vec4<Type> v(*this);
            v.normalize();

            return v;
        }

        //! Negate the vector (i.e. point it in the opposite direction).
        void negate()
        {
            x = -x;
            y = -y;
            z = -z;
            w = -w;
        }

        //! Return this vector reflected off the surface with the given normal \a N. N should be normalized.
        vec4<Type> reflect(const vec4<Type>& N) const
        {
            const vec4<Type>& I(*this);

            return I - 2 * N.dot(I) * N;
        }

        //! Refract this vector through a surface with the given normal \a N and ratio of indices of refraction \a eta.
        vec4<Type> refract(const vec4<Type>& N, Type eta) const
        {
            const vec4<Type>& I(*this);
            const Type k = 1.0 - eta * eta * (1.0 - N.dot(I) * N.dot(I));

            return (k < 0.0) ? 0 : eta * I - (eta * N.dot(I) + std::sqrt(k)) * N;
        }

        //! Index operator. Returns modifiable x, y, z or w value.
        Type& operator[](int i)
        {
            return (&x)[i];
        }

        //! Index operator. Returns x, y, z or w value.
        const Type& operator[](int i) const
        {
            return (&x)[i];
        }

        //! Multiply components of vector with value \a d. Returns reference to self.
        vec4<Type>& operator*=(const Type d)
        {
            x *= d;
            y *= d;
            z *= d;
            w *= d;

            return *this;
        }

        //! Divides components of vector with value \a d. Returns reference to self.
        vec4<Type>& operator/=(const Type d)
        {
            const Type inv = 1.0f / d;

            x *= inv;
            y *= inv;
            z *= inv;
            w *= inv;

            return *this;
        }

        //! Multiply components of vector with value \a a_vec.
        vec4<Type>& operator*=(const vec4<Type>& a_vec)
        {
            x *= a_vec.x;
            y *= a_vec.y;
            z *= a_vec.z;
            w *= a_vec.w;

            return *this;
        }

        //! Adds this vector and vector \a a_vec. Returns reference to self.
        vec4<Type>& operator+=(const vec4<Type>& a_vec)
        {
            x += a_vec.x;
            y += a_vec.y;
            z += a_vec.z;
            w += a_vec.w;

            return *this;
        }

        //! Subtracts vector \a a_vec from this vector. Returns reference to self.
        vec4<Type>& operator-=(const vec4<Type>& a_vec)
        {
            x -= a_vec.x;
            y -= a_vec.y;
            z -= a_vec.z;
            w -= a_vec.w;

            return *this;
        }

        //! Non-destructive negation operator.
        vec4<Type> operator-() const
        {
            return vec4<Type>(-x, -y, -z, -w);
        }

        friend vec4<Type> operator*(const vec4<Type>& a_vec, const Type d)
        {
            return vec4<Type>(a_vec.x * d,
                              a_vec.y * d,
                              a_vec.z * d,
                              a_vec.w * d);
        }
        friend vec4<Type> operator*(const Type d, const vec4<Type>& a_vec)
        {
            return a_vec * d;
        }
        friend vec4<Type> operator/(const vec4<Type>& a_vec, const Type d)
        {
            return vec4<Type>(a_vec.x / d,
                              a_vec.y / d,
                              a_vec.z / d,
                              a_vec.w / d);
        }
        friend vec4<Type> operator*(const vec4<Type>& v1, const vec4<Type>& v2)
        {
            return vec4<Type>(v1.x * v2.x,
                              v1.y * v2.y,
                              v1.z * v2.z,
                              v1.w * v2.w);
        }
        friend vec4<Type> operator+(const vec4<Type>& v1, const vec4<Type>& v2)
        {
            return vec4<Type>(v1.x + v2.x,
                              v1.y + v2.y,
                              v1.z + v2.z,
                              v1.w + v2.w);
        }
        friend vec4<Type> operator-(const vec4<Type>& v1, const vec4<Type>& v2)
        {
            return vec4<Type>(v1.x - v2.x,
                              v1.y - v2.y,
                              v1.z - v2.z,
                              v1.w - v2.w);
        }

        //! Check the two given vector for equality.
        friend bool operator==(const vec4<Type>& v1, const vec4<Type>& v2)
        {
            return (v1.x == v2.x && v1.y == v2.y && v1.z == v2.z && v1.w == v2.w);
        }

        //! Check the two given vector for inequality.
        friend bool operator!=(const vec4<Type>& v1, const vec4<Type>& v2)
        {
            return !(v1 == v2);
        }

        //! Check for equality with given tolerance.
        bool equals(const vec4<Type>& a_vec, const Type a_tolerance = 1E-2) const
        {
            return ((*this - a_vec).sqrLength() <= a_tolerance * a_tolerance);
        }

        friend std::ostream& operator<<(std::ostream& os, const vec4<Type>& vect)
        {
            return os << vect.x << '\t' << vect.y << '\t' << vect.z << '\t' << vect.w;
        }

        //! Largest representable vector
        static vec4<Type> max()
        {
            return vec4<Type>(std::numeric_limits<Type>::max(),
                              std::numeric_limits<Type>::max(),
                              std::numeric_limits<Type>::max(),
                              std::numeric_limits<Type>::max());
        }
    };

    typedef vec4<int> vec4i;
    typedef vec4<float> vec4f;
    typedef vec4<double> vec4d;
} // namespace gtl
