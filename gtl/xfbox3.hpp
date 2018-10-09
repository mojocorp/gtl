#pragma once

#include <gtl/box3.hpp>
#include <gtl/gtl.hpp>
#include <gtl/matrix4.hpp>
#include <gtl/vec3.hpp>

namespace gtl {
    /*!
    \class xfbox3 xfbox3.hpp geometry/xfbox3.hpp
    \brief The xfbox3 class is an object oriented box.
    \ingroup base

    This box class is used by many other classes.

    \sa box3
    */
    template <typename Type>
    class xfbox3 : public box3<Type> {
    public:
        //! The default constructor makes an empty box.
        xfbox3()
            : box3<Type>()
        {
            m_matrix.makeIdentity();
            m_invertedMatrix.makeIdentity();
        }

        //!	Constructs a box with the given corners.
        xfbox3(const vec3<Type>& a_min, const vec3<Type>& a_max)
            : box3<Type>(a_min, a_max)
        {
            m_matrix.makeIdentity();
            m_invertedMatrix.makeIdentity();
        }

        xfbox3(const box3<Type>& a_box)
            : box3<Type>(a_box.m_min, a_box.m_max)
        {
            m_matrix.makeIdentity();
            m_invertedMatrix.makeIdentity();
        }

        void transform(const matrix4<Type>& m)
        {
            setTransform(m_matrix * m);
        }

        void setTransform(const matrix4<Type>& m)
        {
            m_matrix = m;
            m_invertedMatrix = m.inverse();
        }

        const matrix4<Type>& getTransform() const
        {
            return m_matrix;
        }

        const matrix4<Type>& getInverse() const
        {
            return m_invertedMatrix;
        }

        vec3<Type> getCenter() const
        {
            vec3<Type> transcenter;
            m_matrix.multVecMatrix(box3<Type>::getCenter(), transcenter);
            return transcenter;
        }

        void extendBy(const vec3<Type>& pt)
        {
            vec3<Type> trans;
            m_invertedMatrix.multVecMatrix(pt, trans);
            box3<Type>::extendBy(trans);
        }

        void extendBy(const box3<Type>& bb)
        {
            if (bb.isEmpty())
                return;

            if (box3<Type>::isEmpty()) {
                *this = bb;
                m_matrix.makeIdentity();
                m_invertedMatrix.makeIdentity();
                return;
            }
        }
        /*
        void extendBy(const xfbox3<Type> &bb)
        {

        }

        bool intersect(const vec3<Type> &pt) const
        {
        vec3<Type> p;
        m_invertedMatrix.multVecMatrix(pt, p);
        return box3<Type>::intersect(p);
        }

        bool intersect(const box3<Type> &bb) const
        {
        return false;
        }

        bool intersect (const xfbox3<Type> &bb) const
        {
        return false;	
        }
        */
        box3<Type> project() const
        {
            box3<Type> box = (box3<Type>)*this;

            if (!box.isEmpty())
                box.transform(m_matrix);

            return box;
        }

        Type getVolume() const
        {
            if (box3<Type>::isEmpty())
                return 0.0;

            return std::abs(box3<Type>::getVolume() * m_matrix.det3());
        }

        //! Check \a b1 and \a b2 for equality.
        friend bool operator==(const xfbox3<Type>& b1, const xfbox3<Type>& b2)
        {
            return (b1.getMin() == b2.getMin() && b1.getMax() == b2.getMax() && b1.m_matrix == b2.m_matrix);
        }

        //! Check \a b1 and \a b2 for inequality.
        friend bool operator!=(const xfbox3<Type>& b1, const xfbox3<Type>& b2)
        {
            return !(b1 == b2);
        }

    private:
        matrix4<Type> m_matrix;
        matrix4<Type> m_invertedMatrix;
    };

    typedef xfbox3<int> xfbox3i;
    typedef xfbox3<float> xfbox3f;
    typedef xfbox3<double> xfbox3d;
} // namespace gtl
