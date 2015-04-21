#pragma once

#include <gtl/gtl.hpp>
#include <gtl/vec3.hpp>
#include <gtl/quaternion.hpp>

namespace gtl
{
    /*!
    \class matrix3 matrix3.hpp geometry/matrix3.hpp
    \brief Row-major 3x3 matrix class.
    \ingroup base

    This class is used by many other classes.

    \sa matrix4
    */
    template<typename Type>
    class matrix3
    {
    public:
        //! The default constructor. The matrix will be identity.
        matrix3()
        {
            makeIdentity();
        }

        //! Constructs a matrix instance with the given initial elements.
        matrix3(const Type a11, const Type a12, const Type a13,
                const Type a21, const Type a22, const Type a23,
                const Type a31, const Type a32, const Type a33)
        {
            m_data[0][0] = a11; m_data[0][1] = a12; m_data[0][2] = a13;
            m_data[1][0] = a21; m_data[1][1] = a22; m_data[1][2] = a23;
            m_data[2][0] = a31; m_data[2][1] = a32; m_data[2][2] = a33;
        }

        //! Constructs a matrix instance with the initial elements from the \a matrix argument.
        matrix3(const matrix3<Type> & a_matrix)
        {
            *this = a_matrix;
        }

        //! Create from matrix4 excluding row and column k
        matrix3(const matrix4<Type> & a_matrix, int k)
        {
            for (int i = 0; i< 3; i++) {
                for (int j = 0; j< 3; j++) {
                    m_data[i][j] = a_matrix[(i>=k)?i+1:i][(j>=k)?j+1:j];
                }
            }
        }

        //! Set the matrix to be the identity matrix. \sa isIdentity().
        void makeIdentity()
        {
            m_data[0][0]=m_data[1][1]=m_data[2][2]= (Type)1.0;
            m_data[0][1]=m_data[0][2]=m_data[1][0]= (Type)0.0;
            m_data[1][2]=m_data[2][0]=m_data[2][1]= (Type)0.0;
        }

        //! Check if matrix is identity. 
        bool isIdentity() const
        {
            return ((m_data[0][0] == 1.0) && (m_data[0][1] == (Type)0.0) && (m_data[0][2] == (Type)0.0) &&
                    (m_data[1][0] == 0.0) && (m_data[1][1] == (Type)1.0) && (m_data[1][2] == (Type)0.0) &&
                    (m_data[2][0] == 0.0) && (m_data[2][1] == (Type)0.0) && (m_data[2][2] == (Type)1.0));

        }

        //! Inverse the current matrix.	
        matrix3<Type> & invert()
        {
            Type tmp[3][3];
            // Invert using cofactors.
            tmp[0][0] = m_data[1][1]*m_data[2][2] - m_data[1][2]*m_data[2][1];
            tmp[0][1] = m_data[0][2]*m_data[2][1] - m_data[0][1]*m_data[2][2];
            tmp[0][2] = m_data[0][1]*m_data[1][2] - m_data[0][2]*m_data[1][1];
            tmp[1][0] = m_data[1][2]*m_data[2][0] - m_data[1][0]*m_data[2][2];
            tmp[1][1] = m_data[0][0]*m_data[2][2] - m_data[0][2]*m_data[2][0];
            tmp[1][2] = m_data[0][2]*m_data[1][0] - m_data[0][0]*m_data[1][2];
            tmp[2][0] = m_data[1][0]*m_data[2][1] - m_data[1][1]*m_data[2][0];
            tmp[2][1] = m_data[0][1]*m_data[2][0] - m_data[0][0]*m_data[2][1];
            tmp[2][2] = m_data[0][0]*m_data[1][1] - m_data[0][1]*m_data[1][0];

            const Type mInvDet = (Type)1.0 / (m_data[0][0]*tmp[0][0]+m_data[0][1]*tmp[1][0]+m_data[0][2]*tmp[2][0]);

            m_data[0][0] = tmp[0][0] * mInvDet; m_data[0][1] = tmp[0][1] * mInvDet; m_data[0][2] = tmp[0][2] * mInvDet;
            m_data[1][0] = tmp[1][0] * mInvDet; m_data[1][1] = tmp[1][1] * mInvDet; m_data[1][2] = tmp[1][2] * mInvDet;
            m_data[2][0] = tmp[2][0] * mInvDet; m_data[2][1] = tmp[2][1] * mInvDet; m_data[2][2] = tmp[2][2] * mInvDet;

            return *this;
        }


        //! Return a new matrix which is the inverse matrix of this.	
        matrix3<Type> inverse() const
        {
            matrix3<Type> mat(*this);

            return mat.invert();
        }

        //! Transpose the current matrix.
        matrix3<Type> & transpose()
        {
            *this = matrix3<Type>(m_data[0][0], m_data[1][0], m_data[2][0],
                                  m_data[0][1], m_data[1][1], m_data[2][1],
                                  m_data[0][2], m_data[1][2], m_data[2][2]);
            return *this;
        }

        //! Returns the determinant of the matrix.
        Type det() const
        {
            return (m_data[0][0]*m_data[1][1]*m_data[2][2] +
                    m_data[0][1]*m_data[1][2]*m_data[2][0] +
                    m_data[0][2]*m_data[1][0]*m_data[2][1] -
                    m_data[0][2]*m_data[1][1]*m_data[2][0] -
                    m_data[0][1]*m_data[1][0]*m_data[2][2] -
                    m_data[0][0]*m_data[1][2]*m_data[2][1]);
        }

        //! Check if the \a a_matrix matrix is equal to this one, within the given tolerance value.
        bool equals(const matrix3<Type> & a_matrix, Type a_tolerance=1E-2) const
        {
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    if (std::abs(m_data[i][j] - a_matrix.m_data[i][j]) > a_tolerance)
                        return false;
                }
            }
            return true;
        }

        //! Return pointer to the matrix' 3x3 array.
        operator Type * ()
        { 
            return &(m_data[0][0]); 
        }

        //! Returns pointer to the 3 element array representing a matrix row. \a i should be within [0, 2].
        Type * operator [](int i)
        { 
            return m_data[i]; 
        }

        //! Returns pointer to the 3 element array representing a matrix row. \a i should be within [0, 2].
        const Type * operator [](int i) const
        { 
            return m_data[i]; 
        }

        //! Assignment operator. Copies the elements from \a a_matrix to the matrix.
        matrix3<Type> & operator =(const matrix3<Type> & m)
        {
            m_data[0][0] = m.m_data[0][0];
            m_data[0][1] = m.m_data[0][1];
            m_data[0][2] = m.m_data[0][2];
            m_data[1][0] = m.m_data[1][0];
            m_data[1][1] = m.m_data[1][1];
            m_data[1][2] = m.m_data[1][2];
            m_data[2][0] = m.m_data[2][0];
            m_data[2][1] = m.m_data[2][1];
            m_data[2][2] = m.m_data[2][2];

            return *this;
        }

        matrix3<Type> & operator =(const quaternion<Type> & a_quat)
        {
            setRotate(a_quat);

            return *this;
        }

        //! Right-multiply with the \a a_matrix matrix. \sa multRight().
        matrix3<Type> & operator *=(const matrix3<Type> & a_matrix)
        {
            return multRight(a_matrix);
        }

        //! Multiplies matrix \a m1 with matrix \a m2 and returns the resultant matrix.
        friend matrix3<Type> operator *(const matrix3<Type> & m1, const matrix3<Type> & m2)
        { 
            matrix3<Type> m = m1;
            m *= m2; 
            return m; 
        }

        /*! Compare matrices to see if they are equal. For two matrices to be equal,
        all their individual elements must be equal.

        \sa equals().
        */
        friend bool operator ==(const matrix3<Type> & m1, const matrix3<Type> & m2)
        {
            return (
                m1.m_data[0][0] == m2.m_data[0][0] &&
                m1.m_data[0][1] == m2.m_data[0][1] &&
                m1.m_data[0][2] == m2.m_data[0][2] &&

                m1.m_data[1][0] == m2.m_data[1][0] &&
                m1.m_data[1][1] == m2.m_data[1][1] &&
                m1.m_data[1][2] == m2.m_data[1][2] &&

                m1.m_data[2][0] == m2.m_data[2][0] &&
                m1.m_data[2][1] == m2.m_data[2][1] &&
                m1.m_data[2][2] == m2.m_data[2][2]
            );
        }


        friend bool operator !=(const matrix3<Type> & m1, const matrix3<Type> & m2)
        { 
            return !(m1 == m2);
        }

        void setRotate(const quaternion<Type> & a_quat)
        {
            const vec4<Type> & values = a_quat.getValue();

            m_data[0][0] = (Type)(1 - 2.0 * (values[1] * values[1] + values[2] * values[2]));
            m_data[0][1] = (Type)(    2.0 * (values[0] * values[1] + values[2] * values[3]));
            m_data[0][2] = (Type)(    2.0 * (values[2] * values[0] - values[1] * values[3]));

            m_data[1][0] = (Type)(    2.0 * (values[0] * values[1] - values[2] * values[3]));
            m_data[1][1] = (Type)(1 - 2.0 * (values[2] * values[2] + values[0] * values[0]));
            m_data[1][2] = (Type)(    2.0 * (values[1] * values[2] + values[0] * values[3]));

            m_data[2][0] = (Type)(    2.0 * (values[2] * values[0] + values[1] * values[3]));
            m_data[2][1] = (Type)(    2.0 * (values[1] * values[2] - values[0] * values[3]));
            m_data[2][2] = (Type)(1 - 2.0 * (values[1] * values[1] + values[0] * values[0]));
        }

        //! Set matrix to be a pure scaling matrix. (no translation or rotation components)
        void setScale(const Type s)
        {
            setScale(vec3<Type>(s,s,s));
        }

        //! Set matrix to be a pure scaling matrix. (no translation or rotation components)
        void setScale(const vec3<Type> & s)
        {
            makeIdentity();
            m_data[0][0] = s[0];
            m_data[1][1] = s[1];
            m_data[2][2] = s[2];
        }

        //! Let this matrix be right-multiplied by m. Returns reference to self.
        matrix3<Type> & multRight(const matrix3<Type> & m)
        {
            // Trivial cases
            if (m.isIdentity())
                return *this;
            else if (isIdentity())
                return (*this = m);

            matrix3<Type> tmp(
                m_data[0][0]*m[0][0]+m_data[0][1]*m[1][0]+m_data[0][2]*m[2][0],
                m_data[0][0]*m[0][1]+m_data[0][1]*m[1][1]+m_data[0][2]*m[2][1],
                m_data[0][0]*m[0][2]+m_data[0][1]*m[1][2]+m_data[0][2]*m[2][2],
                m_data[1][0]*m[0][0]+m_data[1][1]*m[1][0]+m_data[1][2]*m[2][0],
                m_data[1][0]*m[0][1]+m_data[1][1]*m[1][1]+m_data[1][2]*m[2][1],
                m_data[1][0]*m[0][2]+m_data[1][1]*m[1][2]+m_data[1][2]*m[2][2],
                m_data[2][0]*m[0][0]+m_data[2][1]*m[1][0]+m_data[2][2]*m[2][0],
                m_data[2][0]*m[0][1]+m_data[2][1]*m[1][1]+m_data[2][2]*m[2][1],
                m_data[2][0]*m[0][2]+m_data[2][1]*m[1][2]+m_data[2][2]*m[2][2]);

            return (*this = tmp);
        }

        //! Multiplies matrix by given column vector, giving vector result. 
        void multMatrixVec(const vec3<Type> & src, vec3<Type> & dst) const
        {
            Type x = m_data[0][0]*src[0] + m_data[0][1]*src[1] + m_data[0][2]*src[2];
            Type y = m_data[1][0]*src[0] + m_data[1][1]*src[1] + m_data[1][2]*src[2];
            Type z = m_data[2][0]*src[0] + m_data[2][1]*src[1] + m_data[2][2]*src[2];

            dst.setValue(x, y, z);
        }

        //! Multiplies given row vector by matrix, giving vector result. 
        void multVecMatrix(const vec3<Type> & src, vec3<Type> & dst) const
        {
            Type x = src[0]*m_data[0][0] + src[1]*m_data[1][0] + src[2]*m_data[2][0];
            Type y = src[0]*m_data[0][1] + src[1]*m_data[1][1] + src[2]*m_data[2][1];
            Type z = src[0]*m_data[0][2] + src[1]*m_data[1][2] + src[2]*m_data[2][2];

            dst.setValue(x, y, z);
        }

        friend std::ostream & operator<<(std::ostream & os, const matrix3<Type> & mat)
        { 
            for (unsigned int i=0; i<3; i++) {
                os << mat[i][0] <<'\t'<< mat[i][1] <<'\t'<< mat[i][2] <<std::endl;
            }
            return os;
        }
    private:
        Type m_data[3][3];
    };

    typedef matrix3<int>    matrix3i;
    typedef matrix3<float>  matrix3f;
    typedef matrix3<double> matrix3d;
} // namespace gtl
