#pragma once

#include <gtl/gtl.hpp>
#include <gtl/quaternion.hpp>
#include <gtl/vec3.hpp>

namespace gtl {
    // forward declaration
    template <typename Type>
    class quaternion;

    /*!
    \class matrix4 matrix4.hpp geometry/matrix4.hpp
    \brief Row-major 4x4 matrix class.
    \ingroup base

    This class is used by many other classes.

    \sa matrix3
    */
    template <typename Type>
    class matrix4 {
    public:
        //! The default constructor. The matrix will be identity.
        matrix4()
        {
            makeIdentity();
        }

        //! Constructs a matrix instance with the given initial elements.
        matrix4(const Type a11, const Type a12, const Type a13, const Type a14,
                const Type a21, const Type a22, const Type a23, const Type a24,
                const Type a31, const Type a32, const Type a33, const Type a34,
                const Type a41, const Type a42, const Type a43, const Type a44)
        {
            m_data[0][0] = a11;
            m_data[0][1] = a12;
            m_data[0][2] = a13;
            m_data[0][3] = a14;
            m_data[1][0] = a21;
            m_data[1][1] = a22;
            m_data[1][2] = a23;
            m_data[1][3] = a24;
            m_data[2][0] = a31;
            m_data[2][1] = a32;
            m_data[2][2] = a33;
            m_data[2][3] = a34;
            m_data[3][0] = a41;
            m_data[3][1] = a42;
            m_data[3][2] = a43;
            m_data[3][3] = a44;
        }

        //! Constructs a matrix instance with the initial elements from the \a matrix argument.
        matrix4(const matrix4<Type>& a_matrix)
        {
            *this = a_matrix;
        }

        //! Set the matrix to be the identity matrix. \sa isIdentity().
        void makeIdentity()
        {
            m_data[0][0] = m_data[1][1] = m_data[2][2] = m_data[3][3] = (Type)1.0;
            m_data[0][1] = m_data[0][2] = m_data[0][3] = (Type)0.0;
            m_data[1][0] = m_data[1][2] = m_data[1][3] = (Type)0.0;
            m_data[2][0] = m_data[2][1] = m_data[2][3] = (Type)0.0;
            m_data[3][0] = m_data[3][1] = m_data[3][2] = (Type)0.0;
        }

        //! Check if matrix is identity.
        bool isIdentity() const
        {
            return ((m_data[0][0] == 1.0) && (m_data[0][1] == (Type)0.0) && (m_data[0][2] == (Type)0.0) && (m_data[0][3] == (Type)0.0) && (m_data[1][0] == 0.0) && (m_data[1][1] == (Type)1.0) && (m_data[1][2] == (Type)0.0) && (m_data[1][3] == (Type)0.0) && (m_data[2][0] == 0.0) && (m_data[2][1] == (Type)0.0) && (m_data[2][2] == (Type)1.0) && (m_data[2][3] == (Type)0.0) && (m_data[3][0] == 0.0) && (m_data[3][1] == (Type)0.0) && (m_data[3][2] == (Type)0.0) && (m_data[3][3] == (Type)1.0));
        }

        //! Inverse the current matrix.
        matrix4<Type>& invert()
        {
            // inversion with Cramer's Rule.
            Type tmp[12];
            Type dst[16];

            transpose();

            //calculate pairs for first 8 elements (cofactor)
            tmp[0] = m_data[2][2] * m_data[3][3];
            tmp[1] = m_data[3][2] * m_data[2][3];
            tmp[2] = m_data[1][2] * m_data[3][3];
            tmp[3] = m_data[3][2] * m_data[1][3];
            tmp[4] = m_data[1][2] * m_data[2][3];
            tmp[5] = m_data[2][2] * m_data[1][3];
            tmp[6] = m_data[0][2] * m_data[3][3];
            tmp[7] = m_data[3][2] * m_data[0][3];
            tmp[8] = m_data[0][2] * m_data[2][3];
            tmp[9] = m_data[2][2] * m_data[0][3];
            tmp[10] = m_data[0][2] * m_data[1][3];
            tmp[11] = m_data[1][2] * m_data[0][3];

            //calculate first 8 elements (cofactor)
            dst[0] = tmp[0] * m_data[1][1] + tmp[3] * m_data[2][1] + tmp[4] * m_data[3][1];
            dst[0] -= tmp[1] * m_data[1][1] + tmp[2] * m_data[2][1] + tmp[5] * m_data[3][1];
            dst[1] = tmp[1] * m_data[0][1] + tmp[6] * m_data[2][1] + tmp[9] * m_data[3][1];
            dst[1] -= tmp[0] * m_data[0][1] + tmp[7] * m_data[2][1] + tmp[8] * m_data[3][1];
            dst[2] = tmp[2] * m_data[0][1] + tmp[7] * m_data[1][1] + tmp[10] * m_data[3][1];
            dst[2] -= tmp[3] * m_data[0][1] + tmp[6] * m_data[1][1] + tmp[11] * m_data[3][1];
            dst[3] = tmp[5] * m_data[0][1] + tmp[8] * m_data[1][1] + tmp[11] * m_data[2][1];
            dst[3] -= tmp[4] * m_data[0][1] + tmp[9] * m_data[1][1] + tmp[10] * m_data[2][1];
            dst[4] = tmp[1] * m_data[1][0] + tmp[2] * m_data[2][0] + tmp[5] * m_data[3][0];
            dst[4] -= tmp[0] * m_data[1][0] + tmp[3] * m_data[2][0] + tmp[4] * m_data[3][0];
            dst[5] = tmp[0] * m_data[0][0] + tmp[7] * m_data[2][0] + tmp[8] * m_data[3][0];
            dst[5] -= tmp[1] * m_data[0][0] + tmp[6] * m_data[2][0] + tmp[9] * m_data[3][0];
            dst[6] = tmp[3] * m_data[0][0] + tmp[6] * m_data[1][0] + tmp[11] * m_data[3][0];
            dst[6] -= tmp[2] * m_data[0][0] + tmp[7] * m_data[1][0] + tmp[10] * m_data[3][0];
            dst[7] = tmp[4] * m_data[0][0] + tmp[9] * m_data[1][0] + tmp[10] * m_data[2][0];
            dst[7] -= tmp[5] * m_data[0][0] + tmp[8] * m_data[1][0] + tmp[11] * m_data[2][0];

            //calculate pairs for second 8 elements (cofactors)
            tmp[0] = m_data[2][0] * m_data[3][1];
            tmp[1] = m_data[3][0] * m_data[2][1];
            tmp[2] = m_data[1][0] * m_data[3][1];
            tmp[3] = m_data[3][0] * m_data[1][1];
            tmp[4] = m_data[1][0] * m_data[2][1];
            tmp[5] = m_data[2][0] * m_data[1][1];
            tmp[6] = m_data[0][0] * m_data[3][1];
            tmp[7] = m_data[3][0] * m_data[0][1];
            tmp[8] = m_data[0][0] * m_data[2][1];
            tmp[9] = m_data[2][0] * m_data[0][1];
            tmp[10] = m_data[0][0] * m_data[1][1];
            tmp[11] = m_data[1][0] * m_data[0][1];

            //calculate second 8 elements (cofactors)
            dst[8] = tmp[0] * m_data[1][3] + tmp[3] * m_data[2][3] + tmp[4] * m_data[3][3];
            dst[8] -= tmp[1] * m_data[1][3] + tmp[2] * m_data[2][3] + tmp[5] * m_data[3][3];
            dst[9] = tmp[1] * m_data[0][3] + tmp[6] * m_data[2][3] + tmp[9] * m_data[3][3];
            dst[9] -= tmp[0] * m_data[0][3] + tmp[7] * m_data[2][3] + tmp[8] * m_data[3][3];
            dst[10] = tmp[2] * m_data[0][3] + tmp[7] * m_data[1][3] + tmp[10] * m_data[3][3];
            dst[10] -= tmp[3] * m_data[0][3] + tmp[6] * m_data[1][3] + tmp[11] * m_data[3][3];
            dst[11] = tmp[5] * m_data[0][3] + tmp[8] * m_data[1][3] + tmp[11] * m_data[2][3];
            dst[11] -= tmp[4] * m_data[0][3] + tmp[9] * m_data[1][3] + tmp[10] * m_data[2][3];
            dst[12] = tmp[2] * m_data[2][2] + tmp[5] * m_data[3][2] + tmp[1] * m_data[1][2];
            dst[12] -= tmp[4] * m_data[3][2] + tmp[0] * m_data[1][2] + tmp[3] * m_data[2][2];
            dst[13] = tmp[8] * m_data[3][2] + tmp[0] * m_data[0][2] + tmp[7] * m_data[2][2];
            dst[13] -= tmp[6] * m_data[2][2] + tmp[9] * m_data[3][2] + tmp[1] * m_data[0][2];
            dst[14] = tmp[6] * m_data[1][2] + tmp[11] * m_data[3][2] + tmp[3] * m_data[0][2];
            dst[14] -= tmp[10] * m_data[3][2] + tmp[2] * m_data[0][2] + tmp[7] * m_data[1][2];
            dst[15] = tmp[10] * m_data[2][2] + tmp[4] * m_data[0][2] + tmp[9] * m_data[1][2];
            dst[15] -= tmp[8] * m_data[1][2] + tmp[11] * m_data[2][2] + tmp[5] * m_data[0][2];

            //calculate determinant
            const Type det = (Type)1.0 / (m_data[0][0] * dst[0] + m_data[1][0] * dst[1] + m_data[2][0] * dst[2] + m_data[3][0] * dst[3]);

            //calculate matrix inverse
            m_data[0][0] = dst[0] * det;
            m_data[0][1] = dst[4] * det;
            m_data[0][2] = dst[8] * det;
            m_data[0][3] = dst[12] * det;
            m_data[1][0] = dst[1] * det;
            m_data[1][1] = dst[5] * det;
            m_data[1][2] = dst[9] * det;
            m_data[1][3] = dst[13] * det;
            m_data[2][0] = dst[2] * det;
            m_data[2][1] = dst[6] * det;
            m_data[2][2] = dst[10] * det;
            m_data[2][3] = dst[14] * det;
            m_data[3][0] = dst[3] * det;
            m_data[3][1] = dst[7] * det;
            m_data[3][2] = dst[11] * det;
            m_data[3][3] = dst[15] * det;

            return *this;
        }

        //! Return a new matrix which is the inverse matrix of this.
        matrix4<Type> inverse() const
        {
            matrix4<Type> mat(*this);

            return mat.invert();
        }

        //! Transpose the current matrix.
        matrix4<Type>& transpose()
        {
            *this = matrix4<Type>(m_data[0][0], m_data[1][0], m_data[2][0], m_data[3][0],
                                  m_data[0][1], m_data[1][1], m_data[2][1], m_data[3][1],
                                  m_data[0][2], m_data[1][2], m_data[2][2], m_data[3][2],
                                  m_data[0][3], m_data[1][3], m_data[2][3], m_data[3][3]);
            return *this;
        }

        //! Returns the determinant of the 3x3 submatrix specified by the row and column indices.
        Type det3(int r1, int r2, int r3, int c1, int c2, int c3) const
        {
            // More or less directly from "Advanced Engineering Mathematics"
            // (E. Kreyszig), 6th edition.
            Type a11 = m_data[r1][c1];
            Type a12 = m_data[r1][c2];
            Type a13 = m_data[r1][c3];
            Type a21 = m_data[r2][c1];
            Type a22 = m_data[r2][c2];
            Type a23 = m_data[r2][c3];
            Type a31 = m_data[r3][c1];
            Type a32 = m_data[r3][c2];
            Type a33 = m_data[r3][c3];

            Type M11 = a22 * a33 - a32 * a23;
            Type M21 = -(a12 * a33 - a32 * a13);
            Type M31 = a12 * a23 - a22 * a13;

            return (a11 * M11 + a21 * M21 + a31 * M31);
        }

        //! Returns the determinant of the upper left 3x3 submatrix.
        Type det3() const
        {
            return this->det3(0, 1, 2, 0, 1, 2);
        }

        //! Returns the determinant of the matrix.
        Type det4() const
        {
            Type det = 0.0;
            det += m_data[0][0] * det3(1, 2, 3, 1, 2, 3);
            det -= m_data[1][0] * det3(0, 2, 3, 1, 2, 3);
            det += m_data[2][0] * det3(0, 1, 3, 1, 2, 3);
            det -= m_data[3][0] * det3(0, 1, 2, 1, 2, 3);
            return det;
        }

        //! Check if the \a a_matrix matrix is equal to this one, within the given tolerance value.
        bool equals(const matrix4<Type>& a_matrix, Type a_tolerance = 1E-2) const
        {
            for (int i = 0; i < 4; i++) {
                for (int j = 0; j < 4; j++) {
                    if (std::abs(m_data[i][j] - a_matrix.m_data[i][j]) > a_tolerance)
                        return false;
                }
            }
            return true;
        }

        //! Return pointer to the matrix' 4x4 array.
        operator Type*()
        {
            return &(m_data[0][0]);
        }

        //! Returns pointer to the 4 element array representing a matrix row. \a i should be within [0, 3].
        Type* operator[](int i)
        {
            return m_data[i];
        }

        //! Returns pointer to the 4 element array representing a matrix row. \a i should be within [0, 3].
        const Type* operator[](int i) const
        {
            return m_data[i];
        }

        //! Assignment operator. Copies the elements from \a a_matrix to the matrix.
        matrix4<Type>& operator=(const matrix4<Type>& m)
        {
            m_data[0][0] = m.m_data[0][0];
            m_data[0][1] = m.m_data[0][1];
            m_data[0][2] = m.m_data[0][2];
            m_data[0][3] = m.m_data[0][3];
            m_data[1][0] = m.m_data[1][0];
            m_data[1][1] = m.m_data[1][1];
            m_data[1][2] = m.m_data[1][2];
            m_data[1][3] = m.m_data[1][3];
            m_data[2][0] = m.m_data[2][0];
            m_data[2][1] = m.m_data[2][1];
            m_data[2][2] = m.m_data[2][2];
            m_data[2][3] = m.m_data[2][3];
            m_data[3][0] = m.m_data[3][0];
            m_data[3][1] = m.m_data[3][1];
            m_data[3][2] = m.m_data[3][2];
            m_data[3][3] = m.m_data[3][3];

            return *this;
        }

        matrix4<Type>& operator=(const quaternion<Type>& a_quat)
        {
            setRotate(a_quat);

            return *this;
        }

        //! Right-multiply with the \a a_matrix matrix. \sa multRight().
        matrix4<Type>& operator*=(const matrix4<Type>& a_matrix)
        {
            return multRight(a_matrix);
        }

        //! Multiplies matrix \a m1 with matrix \a m2 and returns the resultant matrix.
        friend matrix4<Type> operator*(const matrix4<Type>& m1, const matrix4<Type>& m2)
        {
            matrix4<Type> m = m1;
            m *= m2;
            return m;
        }

        /*! Compare matrices to see if they are equal. For two matrices to be equal,
        all their individual elements must be equal.

        \sa equals().
        */
        friend bool operator==(const matrix4<Type>& m1, const matrix4<Type>& m2)
        {
            return (
                m1.m_data[0][0] == m2.m_data[0][0] && m1.m_data[0][1] == m2.m_data[0][1] && m1.m_data[0][2] == m2.m_data[0][2] && m1.m_data[0][3] == m2.m_data[0][3] &&

                m1.m_data[1][0] == m2.m_data[1][0] && m1.m_data[1][1] == m2.m_data[1][1] && m1.m_data[1][2] == m2.m_data[1][2] && m1.m_data[1][3] == m2.m_data[1][3] &&

                m1.m_data[2][0] == m2.m_data[2][0] && m1.m_data[2][1] == m2.m_data[2][1] && m1.m_data[2][2] == m2.m_data[2][2] && m1.m_data[2][3] == m2.m_data[2][3] &&

                m1.m_data[3][0] == m2.m_data[3][0] && m1.m_data[3][1] == m2.m_data[3][1] && m1.m_data[3][2] == m2.m_data[3][2] && m1.m_data[3][3] == m2.m_data[3][3]);
        }

        friend bool operator!=(const matrix4<Type>& m1, const matrix4<Type>& m2)
        {
            return !(m1 == m2);
        }

        void setRotate(const quaternion<Type>& a_quat)
        {
            *this = a_quat.getMatrix();
        }

        //! Set matrix to be a pure scaling matrix. (no translation or rotation components)
        void setScale(const Type a_scale)
        {
            setScale(vec3<Type>(a_scale, a_scale, a_scale));
        }

        //! Set matrix to be a pure scaling matrix. (no translation or rotation components)
        void setScale(const vec3<Type>& a_scale)
        {
            makeIdentity();
            m_data[0][0] = a_scale[0];
            m_data[1][1] = a_scale[1];
            m_data[2][2] = a_scale[2];
        }

        //! Make this matrix into a pure translation matrix. (no scale or rotation components)
        void setTranslate(const vec3<Type>& t)
        {
            makeIdentity();
            m_data[3][0] = t[0];
            m_data[3][1] = t[1];
            m_data[3][2] = t[2];
        }

        //! Let this matrix be right-multiplied by m. Returns reference to self.
        matrix4<Type>& multRight(const matrix4<Type>& m)
        {
            // Trivial cases
            if (m.isIdentity())
                return *this;
            else if (isIdentity())
                return (*this = m);

            matrix4<Type> tmp;

            tmp[0][0] = m_data[0][0] * m.m_data[0][0] + m_data[0][1] * m.m_data[1][0] + m_data[0][2] * m.m_data[2][0] + m_data[0][3] * m.m_data[3][0];
            tmp[0][1] = m_data[0][0] * m.m_data[0][1] + m_data[0][1] * m.m_data[1][1] + m_data[0][2] * m.m_data[2][1] + m_data[0][3] * m.m_data[3][1];
            tmp[0][2] = m_data[0][0] * m.m_data[0][2] + m_data[0][1] * m.m_data[1][2] + m_data[0][2] * m.m_data[2][2] + m_data[0][3] * m.m_data[3][2];
            tmp[0][3] = m_data[0][0] * m.m_data[0][3] + m_data[0][1] * m.m_data[1][3] + m_data[0][2] * m.m_data[2][3] + m_data[0][3] * m.m_data[3][3];
            tmp[1][0] = m_data[1][0] * m.m_data[0][0] + m_data[1][1] * m.m_data[1][0] + m_data[1][2] * m.m_data[2][0] + m_data[1][3] * m.m_data[3][0];
            tmp[1][1] = m_data[1][0] * m.m_data[0][1] + m_data[1][1] * m.m_data[1][1] + m_data[1][2] * m.m_data[2][1] + m_data[1][3] * m.m_data[3][1];
            tmp[1][2] = m_data[1][0] * m.m_data[0][2] + m_data[1][1] * m.m_data[1][2] + m_data[1][2] * m.m_data[2][2] + m_data[1][3] * m.m_data[3][2];
            tmp[1][3] = m_data[1][0] * m.m_data[0][3] + m_data[1][1] * m.m_data[1][3] + m_data[1][2] * m.m_data[2][3] + m_data[1][3] * m.m_data[3][3];
            tmp[2][0] = m_data[2][0] * m.m_data[0][0] + m_data[2][1] * m.m_data[1][0] + m_data[2][2] * m.m_data[2][0] + m_data[2][3] * m.m_data[3][0];
            tmp[2][1] = m_data[2][0] * m.m_data[0][1] + m_data[2][1] * m.m_data[1][1] + m_data[2][2] * m.m_data[2][1] + m_data[2][3] * m.m_data[3][1];
            tmp[2][2] = m_data[2][0] * m.m_data[0][2] + m_data[2][1] * m.m_data[1][2] + m_data[2][2] * m.m_data[2][2] + m_data[2][3] * m.m_data[3][2];
            tmp[2][3] = m_data[2][0] * m.m_data[0][3] + m_data[2][1] * m.m_data[1][3] + m_data[2][2] * m.m_data[2][3] + m_data[2][3] * m.m_data[3][3];
            tmp[3][0] = m_data[3][0] * m.m_data[0][0] + m_data[3][1] * m.m_data[1][0] + m_data[3][2] * m.m_data[2][0] + m_data[3][3] * m.m_data[3][0];
            tmp[3][1] = m_data[3][0] * m.m_data[0][1] + m_data[3][1] * m.m_data[1][1] + m_data[3][2] * m.m_data[2][1] + m_data[3][3] * m.m_data[3][1];
            tmp[3][2] = m_data[3][0] * m.m_data[0][2] + m_data[3][1] * m.m_data[1][2] + m_data[3][2] * m.m_data[2][2] + m_data[3][3] * m.m_data[3][2];
            tmp[3][3] = m_data[3][0] * m.m_data[0][3] + m_data[3][1] * m.m_data[1][3] + m_data[3][2] * m.m_data[2][3] + m_data[3][3] * m.m_data[3][3];

            return (*this = tmp);
        }

        //! Let this matrix be left-multiplied by m. Returns reference to self.
        matrix4<Type>& multLeft(const matrix4<Type>& m)
        {
            // Trivial cases
            if (m.isIdentity())
                return *this;
            else if (isIdentity())
                return (*this = m);

            matrix4<Type> tmp;

            tmp[0][0] = m.m_data[0][0] * m_data[0][0] + m.m_data[0][1] * m_data[1][0] + m.m_data[0][2] * m_data[2][0] + m.m_data[0][3] * m_data[3][0];
            tmp[0][1] = m.m_data[0][0] * m_data[0][1] + m.m_data[0][1] * m_data[1][1] + m.m_data[0][2] * m_data[2][1] + m.m_data[0][3] * m_data[3][1];
            tmp[0][2] = m.m_data[0][0] * m_data[0][2] + m.m_data[0][1] * m_data[1][2] + m.m_data[0][2] * m_data[2][2] + m.m_data[0][3] * m_data[3][2];
            tmp[0][3] = m.m_data[0][0] * m_data[0][3] + m.m_data[0][1] * m_data[1][3] + m.m_data[0][2] * m_data[2][3] + m.m_data[0][3] * m_data[3][3];
            tmp[1][0] = m.m_data[1][0] * m_data[0][0] + m.m_data[1][1] * m_data[1][0] + m.m_data[1][2] * m_data[2][0] + m.m_data[1][3] * m_data[3][0];
            tmp[1][1] = m.m_data[1][0] * m_data[0][1] + m.m_data[1][1] * m_data[1][1] + m.m_data[1][2] * m_data[2][1] + m.m_data[1][3] * m_data[3][1];
            tmp[1][2] = m.m_data[1][0] * m_data[0][2] + m.m_data[1][1] * m_data[1][2] + m.m_data[1][2] * m_data[2][2] + m.m_data[1][3] * m_data[3][2];
            tmp[1][3] = m.m_data[1][0] * m_data[0][3] + m.m_data[1][1] * m_data[1][3] + m.m_data[1][2] * m_data[2][3] + m.m_data[1][3] * m_data[3][3];
            tmp[2][0] = m.m_data[2][0] * m_data[0][0] + m.m_data[2][1] * m_data[1][0] + m.m_data[2][2] * m_data[2][0] + m.m_data[2][3] * m_data[3][0];
            tmp[2][1] = m.m_data[2][0] * m_data[0][1] + m.m_data[2][1] * m_data[1][1] + m.m_data[2][2] * m_data[2][1] + m.m_data[2][3] * m_data[3][1];
            tmp[2][2] = m.m_data[2][0] * m_data[0][2] + m.m_data[2][1] * m_data[1][2] + m.m_data[2][2] * m_data[2][2] + m.m_data[2][3] * m_data[3][2];
            tmp[2][3] = m.m_data[2][0] * m_data[0][3] + m.m_data[2][1] * m_data[1][3] + m.m_data[2][2] * m_data[2][3] + m.m_data[2][3] * m_data[3][3];
            tmp[3][0] = m.m_data[3][0] * m_data[0][0] + m.m_data[3][1] * m_data[1][0] + m.m_data[3][2] * m_data[2][0] + m.m_data[3][3] * m_data[3][0];
            tmp[3][1] = m.m_data[3][0] * m_data[0][1] + m.m_data[3][1] * m_data[1][1] + m.m_data[3][2] * m_data[2][1] + m.m_data[3][3] * m_data[3][1];
            tmp[3][2] = m.m_data[3][0] * m_data[0][2] + m.m_data[3][1] * m_data[1][2] + m.m_data[3][2] * m_data[2][2] + m.m_data[3][3] * m_data[3][2];
            tmp[3][3] = m.m_data[3][0] * m_data[0][3] + m.m_data[3][1] * m_data[1][3] + m.m_data[3][2] * m_data[2][3] + m.m_data[3][3] * m_data[3][3];

            return (*this = tmp);
        }

        //! Multiplies matrix by given column vector, giving vector result.
        void multMatrixVec(const vec3<Type>& src, vec3<Type>& dst) const
        {
            const Type x = m_data[0][0] * src[0] + m_data[0][1] * src[1] + m_data[0][2] * src[2] + m_data[0][3];
            const Type y = m_data[1][0] * src[0] + m_data[1][1] * src[1] + m_data[1][2] * src[2] + m_data[1][3];
            const Type z = m_data[2][0] * src[0] + m_data[2][1] * src[1] + m_data[2][2] * src[2] + m_data[2][3];
            const Type w = m_data[3][0] * src[0] + m_data[3][1] * src[1] + m_data[3][2] * src[2] + m_data[3][3];

            dst.setValue(x / w, y / w, z / w);
        }

        //! Multiplies given row vector by matrix, giving vector result.
        void multVecMatrix(const vec3<Type>& src, vec3<Type>& dst) const
        {
            const Type x = src[0] * m_data[0][0] + src[1] * m_data[1][0] + src[2] * m_data[2][0] + m_data[3][0];
            const Type y = src[0] * m_data[0][1] + src[1] * m_data[1][1] + src[2] * m_data[2][1] + m_data[3][1];
            const Type z = src[0] * m_data[0][2] + src[1] * m_data[1][2] + src[2] * m_data[2][2] + m_data[3][2];
            const Type w = src[0] * m_data[0][3] + src[1] * m_data[1][3] + src[2] * m_data[2][3] + m_data[3][3];

            dst.setValue(x / w, y / w, z / w);
        }

        /*!Multiplies given row vector by matrix, giving vector result. 
        src is assumed to be a direction vector, so translation part of matrix is ignored. 
        Note: if you wish to transform surface points and normals by a matrix, call multVecMatrix() for the points and call multDirMatrix() on the inverse transpose of the matrix for the normals.
        */
        void multDirMatrix(const vec3<Type>& src, vec3<Type>& dst) const
        {
            dst.setValue(src[0] * m_data[0][0] + src[1] * m_data[1][0] + src[2] * m_data[2][0],
                         src[0] * m_data[0][1] + src[1] * m_data[1][1] + src[2] * m_data[2][1],
                         src[0] * m_data[0][2] + src[1] * m_data[1][2] + src[2] * m_data[2][2]);
        }

        friend std::ostream& operator<<(std::ostream& os, const matrix4<Type>& mat)
        {
            for (unsigned int i = 0; i < 4; i++) {
                os << mat[i][0] << '\t' << mat[i][1] << '\t' << mat[i][2] << '\t' << mat[i][3] << std::endl;
            }
            return os;
        }

        //! Returns a null rotation.
        static matrix4<Type> identity()
        {
            return matrix4<Type>(1, 0, 0, 0,
                                 0, 1, 0, 0,
                                 0, 0, 1, 0,
                                 0, 0, 0, 1);
        }

        /*! Multiplies this matrix by another that applies an orthographic
            projection for a window with lower-left corner (\a left, \a bottom),
            upper-right corner (\a right, \a top), and the specified \a nearPlane
            and \a farPlane clipping planes.

            \sa frustum(), perspective()
        */
        void ortho(Type left, Type right, Type bottom, Type top, Type nearPlane, Type farPlane)
        {
            if (left == right || bottom == top || nearPlane == farPlane)
                return;

            const Type width = right - left;
            const Type invheight = top - bottom;
            const Type clip = farPlane - nearPlane;
            matrix4<Type> m;
            m.m_data[0][0] = Type(2.0) / width;
            m.m_data[1][0] = 0;
            m.m_data[2][0] = 0;
            m.m_data[3][0] = -(left + right) / width;
            m.m_data[0][1] = 0;
            m.m_data[1][1] = Type(2.0) / invheight;
            m.m_data[2][1] = 0;
            m.m_data[3][1] = -(top + bottom) / invheight;
            m.m_data[0][2] = 0;
            m.m_data[1][2] = 0;
            m.m_data[2][2] = Type(-2.0) / clip;
            m.m_data[3][2] = -(nearPlane + farPlane) / clip;
            m.m_data[0][3] = 0;
            m.m_data[1][3] = 0;
            m.m_data[2][3] = 0;
            m.m_data[3][3] = 1;

            *this *= m;
        }

        /*! Multiplies this matrix by another that applies a perspective
            frustum projection for a window with lower-left corner (\a left, \a bottom),
            upper-right corner (\a right, \a top), and the specified \a nearPlane
            and \a farPlane clipping planes.

            \sa ortho(), perspective()
        */
        void frustum(Type left, Type right, Type bottom, Type top, Type nearPlane, Type farPlane)
        {
            if (left == right || bottom == top || nearPlane == farPlane)
                return;

            const Type width = right - left;
            const Type invheight = top - bottom;
            const Type clip = farPlane - nearPlane;
            matrix4<Type> m;
            m.m_data[0][0] = Type(2.0) * nearPlane / width;
            m.m_data[1][0] = 0;
            m.m_data[2][0] = (left + right) / width;
            m.m_data[3][0] = 0;
            m.m_data[0][1] = 0;
            m.m_data[1][1] = Type(2.0) * nearPlane / invheight;
            m.m_data[2][1] = (top + bottom) / invheight;
            m.m_data[3][1] = 0;
            m.m_data[0][2] = 0;
            m.m_data[1][2] = 0;
            m.m_data[2][2] = -(nearPlane + farPlane) / clip;
            m.m_data[3][2] = Type(-2.0) * nearPlane * farPlane / clip;
            m.m_data[0][3] = 0;
            m.m_data[1][3] = 0;
            m.m_data[2][3] = -1;
            m.m_data[3][3] = 0;

            *this *= m;
        }

        /*! Multiplies this matrix by another that applies a perspective
            projection. The vertical field of view will be \a verticalAngle degrees
            within a window with a given \a aspectRatio that determines the horizontal
            field of view.
            The projection will have the specified \a nearPlane and \a farPlane clipping
            planes which are the distances from the viewer to the corresponding planes.

            \sa ortho(), frustum()
        */
        void perspective(Type verticalAngle, Type aspectRatio, Type nearPlane, Type farPlane)
        {
            if (nearPlane == farPlane || aspectRatio == 0.0f)
                return;

            const Type radians = DegToRad(verticalAngle / 2.0f);
            const Type sine = std::sin(radians);
            if (sine == 0)
                return;
            const Type cotan = std::cos(radians) / sine;
            const Type clip = farPlane - nearPlane;
            matrix4<Type> m;
            m.m_data[0][0] = cotan / aspectRatio;
            m.m_data[1][0] = 0;
            m.m_data[2][0] = 0;
            m.m_data[3][0] = 0;
            m.m_data[0][1] = 0;
            m.m_data[1][1] = cotan;
            m.m_data[2][1] = 0;
            m.m_data[3][1] = 0;
            m.m_data[0][2] = 0;
            m.m_data[1][2] = 0;
            m.m_data[2][2] = -(nearPlane + farPlane) / clip;
            m.m_data[3][2] = -(Type(2.0) * nearPlane * farPlane) / clip;
            m.m_data[0][3] = 0;
            m.m_data[1][3] = 0;
            m.m_data[2][3] = -1;
            m.m_data[3][3] = 0;

            *this *= m;
        }

        /*! Multiplies this matrix by a viewing matrix derived from an eye
            point. The \a center value indicates the center of the view that
            the \a eye is looking at.  The \a up value indicates which direction
            should be considered up with respect to the \a eye.
        
            \note The \a up vector must not be parallel to the line of sight
            from \a eye to \a center.
        */
        void lookAt(const vec3<Type>& eye, const vec3<Type>& center, const vec3<Type>& up)
        {
            const vec3<Type> forward = (center - eye).normalized();
            const vec3<Type> side = forward.cross(up).normalized();
            const vec3<Type> upVector = side.cross(forward);

            matrix4<Type> m;
            m.m_data[0][0] = side.x();
            m.m_data[1][0] = side.y();
            m.m_data[2][0] = side.z();
            m.m_data[3][0] = -side.dot(eye);
            m.m_data[0][1] = upVector.x();
            m.m_data[1][1] = upVector.y();
            m.m_data[2][1] = upVector.z();
            m.m_data[3][1] = -upVector.dot(eye);
            m.m_data[0][2] = -forward.x();
            m.m_data[1][2] = -forward.y();
            m.m_data[2][2] = -forward.z();
            m.m_data[3][2] = -forward.dot(eye);
            m.m_data[0][3] = 0;
            m.m_data[1][3] = 0;
            m.m_data[2][3] = 0;
            m.m_data[3][3] = 1;

            *this *= m;
        }

        /*! Multiplies this matrix by another that performs the scale and bias
            transformation used by OpenGL to transform from normalized device
            coordinates (NDC) to viewport (window) coordinates. That is it maps
            points from the cube ranging over [-1, 1] in each dimension to the
            viewport with it's near-lower-left corner at (\a left, \a bottom, \a nearPlane)
            and with size (\a width, \a height, \a farPlane - \a nearPlane).
 
            This matches the transform used by the fixed function OpenGL viewport
            transform controlled by the functions glViewport() and glDepthRange().
        */
        void viewport(Type left, Type bottom, Type width, Type height, Type nearPlane, Type farPlane)
        {
            const Type w2 = width / Type(2.0);
            const Type h2 = height / Type(2.0);

            matrix4<Type> m;
            m.m_data[0][0] = w2;
            m.m_data[1][0] = 0;
            m.m_data[2][0] = 0;
            m.m_data[3][0] = left + w2;
            m.m_data[0][1] = 0;
            m.m_data[1][1] = h2;
            m.m_data[2][1] = 0;
            m.m_data[3][1] = bottom + h2;
            m.m_data[0][2] = 0;
            m.m_data[1][2] = 0;
            m.m_data[2][2] = (farPlane - nearPlane) / 2.0f;
            m.m_data[3][2] = (nearPlane + farPlane) / 2.0f;
            m.m_data[0][3] = 0;
            m.m_data[1][3] = 0;
            m.m_data[2][3] = 0;
            m.m_data[3][3] = 1;

            *this *= m;
        }

    private:
        Type m_data[4][4];
    };

    typedef matrix4<int> matrix4i;
    typedef matrix4<float> matrix4f;
    typedef matrix4<double> matrix4d;
} // namespace gtl
