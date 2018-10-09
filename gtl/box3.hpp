#pragma once

#include <gtl/gtl.hpp>
#include <gtl/matrix4.hpp>
#include <gtl/plane.hpp>
#include <gtl/ray.hpp>
#include <gtl/vec3.hpp>

namespace gtl {
    // forward declaration
    template <typename Type>
    class plane;
    template <typename Type>
    class ray;

    /*!
    \class box3 box3.hpp geometry/box3.hpp
    \brief Axis-Aligned 3D Bounding Box Class..
    \ingroup base

    This box class is used by many other classes.

    \sa xfbox3
    */
    template <typename Type>
    class box3 {
    public:
        //! The default constructor makes an empty box.
        box3()
        {
            makeEmpty();
        }

        //! Constructs a box with the given corners.
        box3(const vec3<Type>& a_min, const vec3<Type>& a_max)
        {
            m_min = a_min;
            m_max = a_max;
        }

        //! Reset the boundaries of the box with the given corners.
        void setBounds(const vec3<Type>& a_min, const vec3<Type>& a_max)
        {
            m_min = a_min;
            m_max = a_max;
        }

        //! Check if this has been marked as an empty box. \sa makeEmpty().
        bool isEmpty() const
        {
            return (m_max[0] < m_min[0] || m_max[1] < m_min[1] || m_max[2] < m_min[2]);
        }

        //! Marks this as an empty box.	\sa isEmpty().
        void makeEmpty()
        {
            m_min = vec3<Type>::max();
            m_max = -vec3<Type>::max();
        }

        //! Returns the lower left corner of the box. \sa getCenter(), getMax().
        const vec3<Type>& getMin() const
        {
            return m_min;
        }

        //! Returns the upper right corner of the box. \sa getMin().
        const vec3<Type>& getMax() const
        {
            return m_max;
        }

        //! Returns width, height and depth of box.
        vec3<Type> getSize() const
        {
            return m_max - m_min;
        }

        //! Returns the center point of the box.
        vec3<Type> getCenter() const
        {
            return vec3<Type>((m_max[0] + m_min[0]) * 0.5f,
                              (m_max[1] + m_min[1]) * 0.5f,
                              (m_max[2] + m_min[2]) * 0.5f);
        }

        //! Extend the boundaries of the box by the given point.
        void extendBy(const vec3<Type>& a_point)
        {
            if (isEmpty()) {
                setBounds(a_point, a_point);
            } else {
                if (a_point[0] < m_min[0])
                    m_min[0] = a_point[0];
                if (a_point[1] < m_min[1])
                    m_min[1] = a_point[1];
                if (a_point[2] < m_min[2])
                    m_min[2] = a_point[2];

                if (a_point[0] > m_max[0])
                    m_max[0] = a_point[0];
                if (a_point[1] > m_max[1])
                    m_max[1] = a_point[1];
                if (a_point[2] > m_max[2])
                    m_max[2] = a_point[2];
            }
        }

        //! Extend the boundaries of the box by the given \a a_box parameter.
        void extendBy(const box3<Type>& a_box)
        {
            if (isEmpty()) {
                *this = a_box;
            } else {
                extendBy(a_box.getMin());
                extendBy(a_box.getMax());
            }
        }

        //! Give the volume of the box (0 for an empty box)
        Type getVolume() const
        {
            if (isEmpty())
                return 0.0;

            return (m_max[0] - m_min[0]) * (m_max[1] - m_min[1]) * (m_max[2] - m_min[2]);
        }

        //! Transforms box3 by matrix, enlarging box3 to contain result.
        void transform(const matrix4<Type>& m)
        {
            // a transformed empty box is still empty
            if (isEmpty())
                return;

            vec3<Type> corners[8];
            corners[0] = m_min;
            corners[1][0] = m_min[0];
            corners[1][1] = m_max[1];
            corners[1][2] = m_min[2];
            corners[2][0] = m_max[0];
            corners[2][1] = m_max[1];
            corners[2][2] = m_min[2];
            corners[3][0] = m_max[0];
            corners[3][1] = m_min[1];
            corners[3][2] = m_min[2];
            corners[4] = m_max;
            corners[5][0] = m_min[0];
            corners[5][1] = m_max[1];
            corners[5][2] = m_max[2];
            corners[6][0] = m_min[0];
            corners[6][1] = m_min[1];
            corners[6][2] = m_max[2];
            corners[7][0] = m_max[0];
            corners[7][1] = m_min[1];
            corners[7][2] = m_max[2];

            box3<Type> newbox;
            for (int i = 0; i < 8; ++i) {
                m.multVecMatrix(corners[i], corners[i]);
                newbox.extendBy(corners[i]);
            }

            setBounds(newbox.m_min, newbox.m_max);
        }

        //! Check if \a a_point lies within the boundaries of this box.
        bool intersect(const vec3<Type>& a_point) const
        {
            return !(a_point[0] < m_min[0] || a_point[0] > m_max[0] || a_point[1] < m_min[1] || a_point[1] > m_max[1] || a_point[2] < m_min[2] || a_point[2] > m_max[2]);
        }

        //! Check if the given box lies wholly or partly within the boundaries of this box.
        bool intersect(const box3<Type>& a_box) const
        {
            if ((m_max[0] < a_box.m_min[0]) || (m_min[0] > a_box.m_max[0]) || (m_max[1] < a_box.m_min[1]) || (m_min[1] > a_box.m_max[1]) || (m_max[2] < a_box.m_min[2]) || (m_min[2] > a_box.m_max[2])) {
                return false;
            }
            return true;
        }

        //! Check if the given ray intersect the box.
        bool intersect(const ray<Type>& a_ray, Type& tmin, Type& tmax) const
        {
            // Amy Williams, Steve Barrus, R. Keith Morley, and Peter Shirley
            // "An Efficient and Robust ray-Box Intersection Algorithm"
            //  Journal of graphics tools, 10(1):49-54, 2005
            if (isEmpty())
                return false;

            const vec3<Type> inv_direction(1 / a_ray.getDirection()[0], 1 / a_ray.getDirection()[1], 1 / a_ray.getDirection()[2]);

            const vec3<int> sign(inv_direction[0] < 0, inv_direction[1] < 0, inv_direction[2] < 0);

            tmin = ((sign[0] ? m_max : m_min).x() - a_ray.getOrigin().x()) * inv_direction.x();
            tmax = (((1 - sign[0]) ? m_max : m_min).x() - a_ray.getOrigin().x()) * inv_direction.x();

            const Type tymin = ((sign[1] ? m_max : m_min).y() - a_ray.getOrigin().y()) * inv_direction.y();
            const Type tymax = (((1 - sign[1]) ? m_max : m_min).y() - a_ray.getOrigin().y()) * inv_direction.y();

            if ((tmin > tymax) || (tymin > tmax))
                return false;

            if (tymin > tmin)
                tmin = tymin;
            if (tymax < tmax)
                tmax = tymax;

            const Type tzmin = ((sign[2] ? m_max : m_min).z() - a_ray.getOrigin().z()) * inv_direction.z();
            const Type tzmax = (((1 - sign[2]) ? m_max : m_min).z() - a_ray.getOrigin().z()) * inv_direction.z();

            if ((tmin > tzmax) || (tzmin > tmax))
                return false;

            if (tzmin > tmin)
                tmin = tzmin;
            if (tzmax < tmax)
                tmax = tzmax;

            return (tmin >= 0 && tmax >= 1);
        }

        //! Check if the given plane intersect the box.
        bool intersect(const plane<Type>& a_plane) const
        {
            // Empty boxes can cause problems.
            if (isEmpty())
                return false;

            const vec3<Type>& pnorm = a_plane.getNormal();

            // Use separating axis theorem to test overlap.
            const vec3<Type> vmin(pnorm[0] > 0.0 ? m_min[0] : m_max[0],
                                  pnorm[1] > 0.0 ? m_min[1] : m_max[1],
                                  pnorm[2] > 0.0 ? m_min[2] : m_max[2]);
            if (a_plane.isInHalfSpace(vmin))
                return false;

            const vec3<Type> vmax(pnorm[0] > 0.0 ? m_max[0] : m_min[0],
                                  pnorm[1] > 0.0 ? m_max[1] : m_min[1],
                                  pnorm[2] > 0.0 ? m_max[2] : m_min[2]);
            if (a_plane.isInHalfSpace(vmax))
                return true;

            return false;
        }

        //! Check if the given triangle intersect the box.
        bool intersect(const vec3<Type>& a_p0, const vec3<Type>& a_p1, const vec3<Type>& a_p2) const
        {
            // "Fast 3D Triangle-Box Overlap Testing"
            // Tomas Akenine-Moller
            // Journal of Graphics Tools

            // Bullet 1:
            box3<Type> tribox;
            tribox.extendBy(a_p0);
            tribox.extendBy(a_p1);
            tribox.extendBy(a_p2);

            if (!this->intersect(tribox))
                return false;

            const vec3<Type> boxcenter = this->getCenter();

            // move everything so that the boxcenter is in (0,0,0)
            const vec3<Type> v[3] = { a_p0 - boxcenter,
                                      a_p1 - boxcenter,
                                      a_p2 - boxcenter };

            // compute triangle edges
            const vec3<Type> e[3] = { v[1] - v[0],
                                      v[2] - v[1],
                                      v[0] - v[2] };

            // Bullet 2:
            const vec3<Type> normal = e[0].cross(e[1]);

            if (!this->intersect(plane<Type>(normal, a_p0)))
                return false;

            // Bullet 3:
            const vec3<Type> boxhalfsize = this->getSize() * 0.5f;
            const Type f[3][3] = { { 1.0f, 0.0f, 0.0f },
                                   { 0.0f, 1.0f, 0.0f },
                                   { 0.0f, 0.0f, 1.0f } };

            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    const vec3<Type> a = e[i].cross(f[j]);
                    const float p0 = a.dot(v[0]);
                    const float p1 = a.dot(v[1]);
                    const float p2 = a.dot(v[2]);
                    const float min = min3(p0, p1, p2);
                    const float max = max3(p0, p1, p2);
                    float radius = std::abs(a[0]) * boxhalfsize[0] + std::abs(a[1]) * boxhalfsize[1] + std::abs(a[2]) * boxhalfsize[2];
                    if (min > radius || max < -radius)
                        return false;
                }
            }

            return true;
        }

        //! Check \a b1 and \a b2 for equality.
        friend bool operator==(const box3<Type>& b1, const box3<Type>& b2)
        {
            return b1.getMin() == b2.getMin() && b1.getMax() == b2.getMax();
        }

        //! Check \a b1 and \a b2 for inequality.
        friend bool operator!=(const box3<Type>& b1, const box3<Type>& b2)
        {
            return !(b1 == b2);
        }

    private:
        vec3<Type> m_min;
        vec3<Type> m_max;
    };

    typedef box3<int> box3i;
    typedef box3<float> box3f;
    typedef box3<double> box3d;
} // namespace gtl
