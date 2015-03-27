#pragma once

#include <gtl/gtl.hpp>
#include <gtl/vec2.hpp>

namespace gtl
{
    /*!
    \class box2 box2.hpp geometry/box2.hpp
    \brief Axis-Aligned 2D Bounding Box Class..
    \ingroup base

    This box class is used by many other classes.

    \sa box3
    */
    template<typename Type>
    class box2
    {
    public:
        //! The default constructor makes an empty box.
        box2()
        {
            makeEmpty();
        }

        //!	Constructs a box with the given corners.
        box2(const vec2<Type> & a_min, const vec2<Type> & a_max)
        {
            m_min = a_min;
            m_max = a_max;
        }

        //! Reset the boundaries of the box with the given corners.
        void setBounds(const vec2<Type> & a_min, const vec2<Type> & a_max)
        {
            m_min = a_min;
            m_max = a_max;
        }

        //! Check if this has been marked as an empty box. \sa makeEmpty().
        bool isEmpty() const
        {
            return (m_max[0] < m_min[0] || m_max[1] < m_min[1]);
        }

        //! Marks this as an empty box.	\sa isEmpty().
        void makeEmpty()
        {
            m_min =  vec2<Type>::max();
            m_max = -vec2<Type>::max();
        }

        //! Returns the lower left corner of the box. \sa getCenter(), getMax().
        const vec2<Type> & getMin() const
        { 
            return m_min; 
        }

        //! Returns the upper right corner of the box. \sa getMin().
        const vec2<Type> & getMax() const
        { 
            return m_max; 
        }

        //! Returns width and height of box.
        vec2<Type> getSize() const
        {
            return m_max - m_min;
        }

        //! Returns the center point of the box.
        vec2<Type> getCenter() const
        {
            return vec2<Type>((m_max[0] + m_min[0]) * 0.5f,
                              (m_max[1] + m_min[1]) * 0.5f);
        }

        //! Extend the boundaries of the box by the given point.
        void extendBy(const vec2<Type> & a_point)
        {
            if (isEmpty()) {
                setBounds(a_point, a_point);
            } else {
                if (a_point[0] < m_min[0]) m_min[0] = a_point[0];
                if (a_point[1] < m_min[1]) m_min[1] = a_point[1];

                if (a_point[0] > m_max[0]) m_max[0] = a_point[0];
                if (a_point[1] > m_max[1]) m_max[1] = a_point[1];
            }
        }

        //! Extend the boundaries of the box by the given \a a_box parameter.
        void extendBy(const box2<Type> & a_box)
        {
            if (isEmpty()) {
                *this = a_box;
            } else {
                extendBy(a_box.getMin());
                extendBy(a_box.getMax());
            }
        }

        //! Give the surface of the box
        Type getSurface() const
        {
            const vec2<Type> size = m_max - m_min;
            return size[0] * size[1];
        }

        //! Check if \a a_point lies within the boundaries of this box.
        bool intersect(const vec2<Type> & a_point) const
        {
            return !(a_point[0] < m_min[0] || a_point[0] > m_max[0] ||
                     a_point[1] < m_min[1] || a_point[1] > m_max[1]);
        }

        //! Check if the given box lies wholly or partly within the boundaries of this box.
        bool intersect(const box2<Type> & a_box) const
        {
            if ((m_max[0] < a_box.m_min[0]) || (m_min[0] > a_box.m_max[0]) ||
                (m_max[1] < a_box.m_min[1]) || (m_min[1] > a_box.m_max[1])){
                    return false;
            }
            return true;
        }

        //! Check \a b1 and \a b2 for equality.
        friend bool operator ==(const box2<Type> & b1, const box2<Type> & b2)
        { 
            return b1.getMin() == b2.getMin() && b1.getMax() == b2.getMax(); 
        }

        //! Check \a b1 and \a b2 for inequality.
        friend bool operator !=(const box2<Type> & b1, const box2<Type> & b2)
        { 
            return !(b1 == b2); 
        }

    private:
        vec2<Type> m_min;
        vec2<Type> m_max;
    };

    typedef box2<int>    box2i;
    typedef box2<float>  box2f;
    typedef box2<double> box2d;
} // namespace gtl
