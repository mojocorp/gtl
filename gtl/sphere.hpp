#pragma once

#include <gtl/box3.hpp>
#include <gtl/gtl.hpp>
#include <gtl/ray.hpp>
#include <gtl/vec3.hpp>

namespace gtl {
    /*!
    \class sphere sphere.hpp geometry/sphere.hpp
    \brief Represents a sphere in 3D.
    \ingroup base

    This class is used by many other classes.

    \sa   
    */
    template <typename Type>
    class sphere {
    public:
        //! The default constructor.
        sphere()
        {
            m_center.setValue(0.0, 0.0, 0.0);
            m_radius = 1;
        }

        //! Constructs an instance with initial values from \a a_sphere.
        sphere(const sphere<Type>& a_sphere)
        {
            m_center = a_sphere.m_center;
            m_radius = a_sphere.m_radius;
        }

        //! Construct a sphere given center and radius
        sphere(const vec3<Type>& c, Type r)
        {
            m_center = c;
            m_radius = r;
        }

        //! Construct a sphere given the poles
        sphere(const vec3<Type>& p1, const vec3<Type>& p2)
        {
            setPoles(p1, p2);
        }

        //! Compute a fast approximation of the bounding ball for the given point set.
        sphere(const std::vector<vec3<Type>>& points)
        {
            //based on the algorithm given by [Jack Ritter, 1990]
            if (points.empty())
                return;

            // find a large diameter to start with
            // first get the bounding box
            box3<Type> bbox;
            for (size_t i = 0; i < points.size(); i++) {
                bbox.extendBy(points[i]);
            }

            setPoles(bbox.getMin(), bbox.getMax());

            m_radius /= (Type)2.0;

            for (size_t i = 0; i < points.size(); i++) {
                extendBy(points[i]);
            }
        }

        //! Change the center and radius
        void setValue(const vec3<Type>& c, Type r)
        {
            m_center = c;
            m_radius = r;
        }

        //! Specify pair of antipodal points
        void setPoles(const vec3<Type>& p1, const vec3<Type>& p2)
        {
            m_center = (Type)0.5 * (p1 + p2);
            m_radius = (Type)0.5 * (p1 - p2).length();
        }

        //! Set the center
        void setCenter(const vec3<Type>& c)
        {
            m_center = c;
        }

        //! Set the radius
        void setRadius(Type r)
        {
            m_radius = r;
        }

        //! Return the center
        const vec3<Type>& getCenter() const
        {
            return m_center;
        }

        //! Return the radius
        Type getRadius() const
        {
            return m_radius;
        }

        //! Volume of the sphere
        Type getVolume() const
        {
            return (Type)(M_PI * (4.0 / 3.0) * std::pow(m_radius, 3));
        }

        //! Surface of the sphere
        Type getSurface() const
        {
            return (Type)(M_PI * 4.0 * std::pow(m_radius, 2));
        }

        //! Make the sphere containing a given box
        void circumscribe(const box3<Type>& box)
        {
            m_center = box.getCenter();
            m_radius = (box.getMax() - m_center).length();
        }

        //! Extend the boundaries of the sphere by the given point.
        void extendBy(const vec3<Type>& a_point)
        {
            if (intersect(a_point))
                return;

            vec3<Type> dir = m_center - a_point;

            dir.normalize();

            const vec3<Type> p1 = m_center + m_radius * dir;

            setPoles(p1, a_point);
        }

        //! Extend the boundaries of the sphere by the given sphere.
        void extendBy(const sphere<Type>& sphere)
        {
            if (intersect(sphere))
                return;

            vec3<Type> dir = m_center - sphere.getCenter();

            dir.normalize();

            const vec3<Type> p1 = m_center + m_radius * dir;
            const vec3<Type> p2 = sphere.getCenter() - sphere.getRadius() * dir;

            setPoles(p1, p2);
        }

        //! Returns true if the given point p lies within the sphere.
        bool intersect(const vec3<Type>& p) const
        {
            return (p - m_center).sqrLength() < m_radius * m_radius;
        }

        //! Intersect with a sphere, returning true if there is an intersection.
        bool intersect(const sphere<Type>& s) const
        {
            const Type d1 = (s.getCenter() - m_center).sqrLength();
            const Type d2 = m_radius + s.getRadius();

            return (d1 < d2 * d2);
        }

        //! Intersect ray and sphere, returning true if there is an intersection.
        bool intersect(const ray<Type>& r, Type& t0, Type& t1) const
        {
            const vec3<Type> r_to_s = r.getOrigin() - m_center;

            //Compute A, B and C coefficients
            const Type A = r.getDirection().sqrLength();
            const Type B = 2.0f * r_to_s.dot(r.getDirection());
            const Type C = r_to_s.sqrLength() - m_radius * m_radius;

            //Find discriminant
            Type disc = B * B - 4.0 * A * C;

            // if discriminant is negative there are no real roots
            if (disc < 0.0)
                return false;

            disc = (Type)std::sqrt((double)disc);

            t0 = (-B + disc) / (2.0 * A);
            t1 = (-B - disc) / (2.0 * A);

            // check if we're inside it
            if ((t0 < 0.0 && t1 > 0) || (t0 > 0 && t1 < 0))
                return false;

            if (t0 > t1)
                std::swap(t0, t1);

            return (t0 > 0);
        }

        //! Intersect with an axis aligned box, returning true if there is an intersection.
        bool intersect(const box3<Type>& b) const
        {
            // Arvo's algorithm.
            Type d = 0;

            //find the square of the distance from the sphere to the box
            for (unsigned int i = 0; i < 3; i++) {
                if (m_center[i] < b.getMin()[i]) {
                    Type s = m_center[i] - b.getMin()[i];
                    d += s * s;
                } else if (m_center[i] > b.getMax()[i]) {
                    Type s = m_center[i] - b.getMax()[i];
                    d += s * s;
                }
            }
            return d <= m_radius * m_radius;
        }

        //! Distribute N points on the sphere (uniform).
        void getUniformSurfacePoints(std::vector<vec3<Type>>& points, size_t N)
        {
            // "generalized spiral set" [Saff and Kuijlaars, 1997]
            points.resize(N);

            Type theta = 0.0;
            for (size_t k = 1; k <= N; k++) {
                const Type cosphi = (Type)(-1.0 + 2.0 * (k - 1) / (Type)(N - 1));
                const Type sinphi = (Type)std::sqrt(1.0 - cosphi * cosphi);

                if (k == 1 || k == N)
                    theta = 0.0;
                else
                    theta += (Type)3.6 / (sinphi * std::sqrt((Type)N));

                points[k - 1] = vec3<Type>(m_radius * sinphi * std::cos(theta),
                                           m_radius * sinphi * std::sin(theta),
                                           m_radius * cosphi);
            }
        }

        //! Check the two given sphere for equality.
        friend bool operator==(const sphere<Type>& s1, const sphere<Type>& s2)
        {
            return (s1.m_center == s2.m_center && s1.m_radius == s2.m_radius);
        }

        //! Check the two given sphere for inequality.
        friend bool operator!=(const sphere<Type>& s1, const sphere<Type>& s2)
        {
            return !(s1 == s2);
        }

    private:
        vec3<Type> m_center; //!< sphere center
        Type m_radius; //!< sphere radius
    };

    typedef sphere<int> spherei;
    typedef sphere<float> spheref;
    typedef sphere<double> sphered;
} // namespace gtl
