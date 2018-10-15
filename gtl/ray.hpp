#pragma once

#include <gtl/box3.hpp>
#include <gtl/gtl.hpp>
#include <gtl/plane.hpp>
#include <gtl/vec3.hpp>

namespace gtl {
    /*!
    \class ray ray.hpp geometry/ray.hpp
    \brief Represents a directed ray in 3D.
    \ingroup base

    This class is used by many other classes.

    \sa vec3
    */
    template <typename Type>
    class ray {
    public:
        //! Default constructor. Does nothing.
        ray() {}

        //! Constructs an instance with initial values from \a a_ray.
        ray(const ray<Type>& a_ray)
            : m_origin(a_ray.m_origin)
            , m_direction(a_ray.m_direction)
        {
        }

        //! Create a ray from a_origin with the direction a_direction. If \a normalize is true, \a a_direction will be normalized.
        ray(const vec3<Type>& a_origin, const vec3<Type>& a_direction, bool normalize = true)
        {
            setValue(a_origin, a_direction, normalize);
        }

        //! Set position and direction of the ray. If \a normalize is true, \a a_direction will be normalized.
        void setValue(const vec3<Type>& a_origin, const vec3<Type>& a_direction, bool normalize = true)
        {
            m_origin = a_origin;
            m_direction = normalize ? a_direction.normalized() : a_direction;
        }

        //! Return the ray origin.
        const vec3<Type>& getOrigin() const
        {
            return m_origin;
        }

        //! Return the ray normalized direction.
        const vec3<Type>& getDirection() const
        {
            return m_direction;
        }

        //! Get point on ray at t
        vec3<Type> getValue(Type t) const
        {
            return vec3<Type>(m_origin + t * m_direction);
        }

        //! Distance point line
        Type getDistance(const vec3<Type>& pt) const
        {
            return std::sqrt(getSqrDistance(pt));
        }

        //! Squared distance point line
        Type getSqrDistance(const vec3<Type>& pt) const
        {
            const Type t = (pt - m_origin).dot(m_direction);

            return (getValue(t) - pt).sqrLength();
        }

        //! Reflects about \a normal and return a new ray with origin at \a newOrigin.
        ray<Type> reflect(const vec3<Type>& newOrigin, const vec3<Type>& normal) const
        {
            return ray<Type>(newOrigin, m_direction - 2.0f * (m_direction.dot(normal)) * normal);
        }

        //! Project the given point on the ray.
        vec3<Type> project(const vec3<Type>& pt) const
        {
            const Type numerator = (pt - m_origin).dot(m_direction);
            const Type denumerator = m_direction.length();

            return (m_origin + m_direction * (numerator / denumerator));
        }

        /*! Intersect the ray with the given triangle defined by vert0,vert1,vert2.
        Return true if there is an intersection.
        If there is an intersection, a vector a_tuv is returned, where t is the
        distance to the plane in which the triangle lies and (u,v) represents the
        coordinates inside the triangle.
        */
        bool intersect(const vec3<Type>& vert0, const vec3<Type>& vert1, const vec3<Type>& vert2, vec3<Type>& a_tuv) const
        {
            // Tomas Moller and Ben Trumbore.
            // Fast, minimum storage ray-triangle intersection.
            // Journal of graphics tools, 2(1):21-28, 1997

            // find vectors for two edges sharing vert0
            const vec3<Type> edge1 = vert1 - vert0;
            const vec3<Type> edge2 = vert2 - vert0;

            // begin calculating determinant - also used to calculate U parameter
            const vec3<Type> pvec = m_direction.cross(edge2);

            // if determinant is near zero, ray lies in plane of triangle
            const Type det = edge1.dot(pvec);

            if (det < EPS)
                return false;

            // calculate distance from vert0 to ray origin
            const vec3<Type> tvec = m_origin - vert0;

            // calculate U parameter and test bounds
            a_tuv[1] = tvec.dot(pvec);

            if (a_tuv[1] < 0.0 || a_tuv[1] > det)
                return false;

            // prepare to test V parameter
            const vec3<Type> qvec = tvec.cross(edge1);

            // calculate V parameter and test bounds
            a_tuv[2] = m_origin.dot(qvec);

            if (a_tuv[2] < 0.0 || a_tuv[1] + a_tuv[2] > det)
                return false;

            // calculate t, scale parameters, ray intersects triangle
            a_tuv[0] = edge2.dot(qvec);

            const Type inv_det = (Type)1.0 / det;

            a_tuv[0] *= inv_det;
            a_tuv[1] *= inv_det;
            a_tuv[2] *= inv_det;

            return true;
        }

        /*! Calculate the shortest line between two lines in 3D
          Calculate also the values of mua and mub where
            Pa = P1 + mua (P2 - P1)
            Pb = P3 + mub (P4 - P3)
          Return false if no solution exists.

          Two lines in 3 dimensions generally don't intersect at a point, they may
          be parallel (no intersections) or they may be coincident (infinite intersections)
          but most often only their projection onto a plane intersect..
          When they don't exactly intersect at a point they can be connected by a line segment,
          the shortest line segment is unique and is often considered to be their intersection in 3D.
          */
        bool intersect(const ray<Type>& a_ray, Type& mua, Type& mub) const
        {
            // Based on code from Paul Bourke
            // http://astronomy.swin.edu.au/~pbourke
            const vec3<Type> p1 = m_origin;
            const vec3<Type> p2 = m_origin + 1.0f * m_direction;
            const vec3<Type> p3 = a_ray.m_origin;
            const vec3<Type> p4 = a_ray.m_origin + 1.0f * a_ray.m_direction;

            const vec3<Type> p13 = p1 - p3;
            const vec3<Type> p43 = p4 - p3;

            if (std::abs(p43[0]) < EPS && std::abs(p43[1]) < EPS && std::abs(p43[2]) < EPS)
                return false;

            const vec3<Type> p21 = p2 - p1;

            if (std::abs(p21[0]) < EPS && std::abs(p21[1]) < EPS && std::abs(p21[2]) < EPS)
                return false;

            const Type d1343 = p13.dot(p43);
            const Type d4321 = p43.dot(p21);
            const Type d1321 = p13.dot(p21);
            const Type d4343 = p43.dot(p43);
            const Type d2121 = p21.dot(p21);

            const Type denom = d2121 * d4343 - d4321 * d4321;

            if (std::abs(denom) < EPS)
                return false;

            const Type numer = d1343 * d4321 - d1321 * d4343;

            mua = numer / denom;
            mub = (d1343 + d4321 * mua) / d4343;

            return true;
        }

        //! Check the two given ray for equality.
        friend bool operator==(const ray<Type>& r1, const ray<Type>& r2)
        {
            return (r1.m_origin == r2.m_origin && r1.m_direction == r2.m_direction);
        }

        //! Check the two given ray for inequality.
        friend bool operator!=(const ray<Type>& r1, const ray<Type>& r2)
        {
            return !(r1 == r2);
        }

    private:
        vec3<Type> m_origin; //!< ray origin
        vec3<Type> m_direction; //!< Normalized direction
    };

    typedef ray<int> rayi;
    typedef ray<float> rayf;
    typedef ray<double> rayd;
} // namespace gtl
