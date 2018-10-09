#pragma once

#include <gtl/gtl.hpp>

namespace gtl {
    /*!
    \class complex complex.hpp geometry/complex.hpp
    \brief Represents a complex number.
    \ingroup base
    */
    template <typename Type>
    class complex {
    public:
        //! Constructs an instance using given real and imaginary values.
        complex(Type a_real, Type a_imaginary)
        {
            m_real = a_real;
            m_imag = a_imaginary;
        }

        //! Constructs an instance using values from a given complex instance
        complex(const complex<Type>& c)
        {
            m_real = c.m_real;
            m_imag = c.m_imag;
        }

        //! Set new real and imaginary values for the complex. Returns reference to self.
        complex<Type>& setValue(Type a_real, Type a_imaginary)
        {
            m_real = a_real;
            m_imag = a_imaginary;
            return *this;
        }

        //! Gets the real value of the complex number.
        const Type& getReal() const
        {
            return m_real;
        }

        //! Sets the real value of the complex number.
        void setReal(Type real)
        {
            m_real = real;
        }

        //! Gets the imaginary value of the complex number.
        const Type& getImaginary() const
        {
            return m_imag;
        }

        //! Sets the imaginary value of the complex number.
        void setImaginary(Type imaginary)
        {
            m_imag = imaginary;
        }

        //! Negates the complex number.
        void negate()
        {
            m_real = -m_real;
            m_imag = -m_imag;
        }

        //! Check for equality with given tolerance.
        bool equals(const complex<Type>& a_complex, const Type a_tolerance = 1E-2) const
        {
            return ((std::abs(m_real - a_complex.m_real) <= a_tolerance) && (std::abs(m_imag - a_complex.m_imag) <= a_tolerance));
        }

        //! Gets the modulus of the complex number.
        Type modulus() const
        {
            return (Type)std::sqrt(m_real * m_real + m_imag * m_imag);
        }

        //! Gets the squared modulus of the complex number.
        Type sqrModulus() const
        {
            return (m_real * m_real + m_imag * m_imag);
        }

        //! Returns the complex conjugate.
        complex<Type> getConjugate() const
        {
            return complex<Type>(m_real, -m_imag);
        }

        //! Normalize the complex number.
        void normalize()
        {
            const Type inv_modulus = (Type)(1.0 / modulus());

            m_real = m_real * inv_modulus;
            m_imag = m_imag * inv_modulus;
        }

        //! Check the two given complex for equality.
        friend bool operator==(const complex<Type>& c1, const complex<Type>& c2)
        {
            return c1.m_real == c2.m_real && c1.m_imag == c2.m_imag;
        }

        //! Check the two given complex for inequality.
        friend bool operator!=(const complex<Type>& c1, const complex<Type>& c2)
        {
            return !(c1 == c2);
        }

        //! Non-destructive negation operator.
        complex<Type> operator-() const
        {
            return complex<Type>(-m_real, -m_imag);
        }

        complex<Type>& operator*=(const complex<Type>& c)
        {
            m_real += m_real * c.m_real - m_imag * c.m_imag;
            m_imag += m_real * c.m_imag + m_imag * c.m_real;
            return *this;
        }

        complex<Type>& operator*=(const Type& d)
        {
            m_real *= d;
            m_imag *= d;
            return *this;
        }

        complex<Type>& operator/=(const complex<Type>& c)
        {
            const Type inv_sqr_modulus = (Type)(1.0 / sqrModulus());

            m_real = (m_real * c.m_real + m_imag * c.m_imag) * inv_sqr_modulus;
            m_imag = (m_imag * c.m_real - m_real * c.m_imag) * inv_sqr_modulus;

            return *this;
        }

        complex<Type>& operator/=(const Type& d)
        {
            m_real /= d;
            m_imag /= d;
            return *this;
        }

        complex<Type>& operator+=(const complex<Type>& c)
        {
            m_real += c.m_real;
            m_imag += c.m_imag;
            return *this;
        }

        complex<Type>& operator+=(const Type& d)
        {
            m_real += d;
            return *this;
        }

        complex<Type>& operator-=(const complex<Type>& c)
        {
            m_real -= c.m_real;
            m_imag -= c.m_imag;
            return *this;
        }

        complex<Type>& operator-=(const Type& d)
        {
            m_real -= d;
            return *this;
        }

        friend complex<Type> operator*(const complex<Type>& c1, const complex<Type>& c2)
        {
            return complex<Type>(c1.m_real * c2.m_real - c1.m_imag * c2.m_imag, c1.m_real * c2.m_imag + c1.m_imag * c2.m_real);
        }
        friend complex<Type> operator*(const complex<Type>& c, const Type d)
        {
            return complex<Type>(c.m_real * d, c.m_imag * d);
        }
        friend complex<Type> operator*(const Type d, const complex<Type>& c)
        {
            return c * d;
        }

        friend complex<Type> operator/(const complex<Type>& c, const Type d)
        {
            return complex<Type>(c.m_real / d, c.m_imag / d);
        }
        friend complex<Type> operator/(const complex<Type>& c1, const complex<Type>& c2)
        {
            return complex<Type>(c1 * c2.getConjugate() / c2.sqrModulus());
        }

        friend complex<Type> operator+(const complex<Type>& c1, const complex<Type>& c2)
        {
            return complex<Type>(c1.m_real + c2.m_real, c1.m_imag + c2.m_imag);
        }
        friend complex<Type> operator+(const complex<Type>& c, const Type d)
        {
            return complex<Type>(c.m_real + d, c.m_imag);
        }
        friend complex<Type> operator+(const Type d, const complex<Type>& c)
        {
            return c + d;
        }

        friend complex<Type> operator-(const complex<Type>& c1, const complex<Type>& c2)
        {
            return complex<Type>(c1.m_real - c2.m_real, c1.m_imag - c2.m_imag);
        }
        friend complex<Type> operator-(const complex<Type>& c, const Type d)
        {
            return complex<Type>(c.m_real - d, c.m_imag);
        }
        friend complex<Type> operator-(const Type d, const complex<Type>& c)
        {
            return c - d;
        }

    private:
        Type m_real;
        Type m_imag;
    };

    typedef complex<int> complexi;
    typedef complex<float> complexf;
    typedef complex<double> complexd;
}
