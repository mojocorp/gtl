#pragma once

#include <gtl/gtl.hpp>
#include <gtl/matrix4.hpp>
#include <gtl/vec3.hpp>
#include <gtl/vec4.hpp>

namespace gtl {
    // forward declaration
    template <typename Type>
    class matrix4;

    /*!
    \class quat quat.hpp geometry/quat.hpp
    \brief Object that stores a rotation.
    \ingroup base

    All angles are in degree and all rotations are right-handed.
    The rotation is internally stored as a quaternion.
    This class is used by many other classes.

    \sa matrix4
    */
    template <typename Type>
    class quaternion {
    public:
        //! The default constructor just initializes a valid rotation.
        quaternion()
        {
            makeIdentity();
        }

        //! Construct a quaternion initialized with the given axis-of-rotation and rotation angle.
        quaternion(const vec3<Type>& a_axis, Type a_degree)
        {
            setValue(a_axis, a_degree);
        }

        //! Construct a quaternion initialized with the given quaternion components.
        quaternion(Type q0, Type q1, Type q2, Type q3)
        {
            setValue(q0, q1, q2, q3);
        }

        //! Construct a quaternion initialized with the given rotation matrix.
        quaternion(const matrix4<Type>& a_matrix)
        {
            setValue(a_matrix);
        }

        //! Construct a quaternion which is the minimum rotation necessary to make vector \a rotate_from point in the direction of vector \a rotate_to.
        quaternion(const vec3<Type>& a_rotate_from, const vec3<Type>& a_rotate_to)
        {
            setValue(a_rotate_from, a_rotate_to);
        }

        //! Set the quaternion to be the identity. \sa isIdentity().
        void makeIdentity()
        {
            setValue(0.0, 0.0, 0.0, 1.0);
        }

        //! Check if quaternion is identity.
        bool isIdentity() const
        {
            return ((m_data[0] == 0.0) && (m_data[1] == 0.0) && (m_data[2] == 0.0) && (m_data[3] == 1.0));
        }

        //! Normalizes a rotation quaternion to unit 4D length. Return value is the original length of the vector before normalization.
        Type normalize()
        {
            return m_data.normalize();
        }

        //! Returns the conjugated quaternion */
        quaternion<Type> conjugate() const
        {
            return quaternion<Type>(-m_data[0], -m_data[1], -m_data[2], m_data[3]);
        }

        //! Set the rotation.
        void setValue(Type q0, Type q1, Type q2, Type q3)
        {
            m_data[0] = q0;
            m_data[1] = q1;
            m_data[2] = q2;
            m_data[3] = q3;
            m_data.normalize();
        }

        //! Reset the rotation by the four quaternions in the array.
        void setValue(const Type q[4])
        {
            setValue(q[0], q[1], q[2], q[3]);
        }

        //! Set the rotation from the components of the given matrix.
        void setValue(const matrix4<Type>& a_matrix)
        {
            const Type scalerow = a_matrix[0][0] + a_matrix[1][1] + a_matrix[2][2];

            if (scalerow > 0.0f) {
                Type s = std::sqrt(scalerow + a_matrix[3][3]);

                m_data[3] = s * (Type)0.5;

                s = (Type)0.5 / s;

                m_data[0] = (a_matrix[1][2] - a_matrix[2][1]) * s;
                m_data[1] = (a_matrix[2][0] - a_matrix[0][2]) * s;
                m_data[2] = (a_matrix[0][1] - a_matrix[1][0]) * s;
            } else {
                int i = 0;
                if (a_matrix[1][1] > a_matrix[0][0])
                    i = 1;
                if (a_matrix[2][2] > a_matrix[i][i])
                    i = 2;

                const int j = (i + 1) % 3;
                const int k = (j + 1) % 3;

                float s = std::sqrt((a_matrix[i][i] - (a_matrix[j][j] + a_matrix[k][k])) + a_matrix[3][3]);

                m_data[i] = s * 0.5f;

                s = 0.5f / s;

                m_data[3] = (a_matrix[j][k] - a_matrix[k][j]) * s;
                m_data[j] = (a_matrix[i][j] + a_matrix[j][i]) * s;
                m_data[k] = (a_matrix[i][k] + a_matrix[k][i]) * s;
            }

            if (a_matrix[3][3] != 1.0f) {
                operator*=((Type)1.0 / std::sqrt(a_matrix[3][3]));
            }
        }

        //! Reset rotation with the given axis-of-rotation and rotation angle.
        void setValue(const vec3<Type>& a_axis, Type a_degree)
        {
            m_data[3] = (Type)std::cos(gtl::DegToRad(a_degree) / (Type)2.0);

            const Type sineval = (Type)std::sin(gtl::DegToRad(a_degree) / (Type)2.0);

            const vec3<Type> a = a_axis.normalized();

            m_data[0] = a[0] * sineval;
            m_data[1] = a[1] * sineval;
            m_data[2] = a[2] * sineval;
        }

        //! Construct a rotation which is the minimum rotation necessary to make vector \a rotate_from point in the direction of vector \a rotate_to.
        void setValue(const vec3<Type>& a_rotate_from, const vec3<Type>& a_rotate_to)
        {
            const vec3<Type> from = a_rotate_from.normalized();
            const vec3<Type> to = a_rotate_to.normalized();

            const Type dot = from.dot(to);

            vec3<Type> crossvec = from.cross(to);

            const Type crosslen = crossvec.length();

            if (crosslen == 0.0f) { // Parallel vectors
                // Check if they are pointing in the same direction.
                if (dot > 0.0) {
                    setValue(0.0f, 0.0f, 0.0f, 1.0f);
                } else {
                    // Ok, so they are parallel and pointing in the opposite direction
                    // of each other.
                    // Try crossing with x axis.
                    vec3<Type> t = from.cross(vec3<Type>(1.0f, 0.0f, 0.0f));
                    // If not ok, cross with y axis.
                    if (t.length() == 0.0f)
                        t = from.cross(vec3<Type>(0.0f, 1.0f, 0.0f));

                    t.normalize();
                    setValue(t[0], t[1], t[2], 0.0);
                }
            } else { // Vectors are not parallel
                crossvec.normalize();
                // The fabs() wrapping is to avoid problems when `dot' "overflows"
                // a tiny wee bit, which can lead to sqrt() returning NaN.
                crossvec *= (Type)std::sqrt(0.5 * std::abs(1.0 - dot));

                // The fabs() wrapping is to avoid problems when `dot' "underflows"
                // a tiny wee bit, which can lead to sqrt() returning NaN.
                setValue(crossvec[0], crossvec[1], crossvec[2], (Type)std::sqrt(0.5 * std::abs(1.0 + dot)));
            }
        }

        //! Euler-to-quaternion conversion
        void setValue(Type yaw, Type pitch, Type roll)
        {
            // Gainer, 1972
            // calculate trig identities
            yaw = gtl::DegToRad(yaw) * (Type)0.5;
            pitch = gtl::DegToRad(pitch) * (Type)0.5;
            roll = gtl::DegToRad(roll) * (Type)0.5;

            const Type c1 = std::cos(roll);
            const Type c2 = std::cos(yaw);
            const Type c3 = std::cos(pitch);

            const Type s1 = std::sin(roll);
            const Type s2 = std::sin(yaw);
            const Type s3 = std::sin(pitch);

            const Type c2c1 = c2 * c1;
            const Type s2s1 = s2 * s1;

            m_data[0] = s3 * c2c1 - c3 * s2s1;
            m_data[1] = c3 * s2 * c1 + s3 * c2 * s1;
            m_data[2] = c3 * c2 * s1 - s3 * s2 * c1;
            m_data[3] = c3 * c2c1 + s3 * s2s1;
        }

        //! Return the four quaternion components representing the rotation.
        const vec4<Type>& getValue() const
        {
            return m_data;
        }

        //! Return the rotation in the form of an axis-of-rotation and a rotation angle.
        void getValue(vec3<Type>& a_axis, Type& a_degree) const
        {
            if ((m_data[3] >= -1.0) && (m_data[3] <= 1.0)) {
                a_degree = (Type)(std::acos(m_data[3])) * 2.0f;

                const Type scale = (Type)std::sin(gtl::DegToRad(a_degree) / 2.0f);

                if (scale != 0.0f) {
                    a_axis[0] = m_data[0] / scale;
                    a_axis[1] = m_data[1] / scale;
                    a_axis[2] = m_data[2] / scale;

                    return;
                }
            }

            // quat can't be converted to axis and rotation angle, so we just
            // set up values to indicate this.
            a_axis.setValue(0.0, 0.0, 1.0);
            a_degree = 0.0f;
        }

        void getValue(Type& yaw, Type& pitch, Type& roll) const
        {
            const Type q12 = m_data[0] * m_data[0];
            const Type q22 = m_data[1] * m_data[1];
            const Type q32 = m_data[2] * m_data[2];
            const Type q42 = m_data[3] * m_data[3];

            yaw = gtl::RadToDeg(std::asin(-2.0 * (m_data[0] * m_data[2] - m_data[3] * m_data[1])));
            pitch = gtl::RadToDeg(std::atan(2.0 * (m_data[3] * m_data[0] + m_data[1] * m_data[2]) / (q42 - q12 - q22 + q32)));
            roll = gtl::RadToDeg(std::atan(2.0 * (m_data[0] * m_data[1] + m_data[3] * m_data[2]) / (q42 + q12 - q22 - q32)));
        }

        //! Return this rotation in the form of a matrix.
        matrix4<Type> getMatrix() const
        {
            matrix4<Type> matrix;

            // calculate coefficients
            const Type x2 = m_data[0] + m_data[0];
            const Type y2 = m_data[1] + m_data[1];
            const Type z2 = m_data[2] + m_data[2];
            const Type xx = m_data[0] * x2;
            const Type xy = m_data[0] * y2;
            const Type xz = m_data[0] * z2;
            const Type yy = m_data[1] * y2;
            const Type yz = m_data[1] * z2;
            const Type zz = m_data[2] * z2;
            const Type wx = m_data[3] * x2;
            const Type wy = m_data[3] * y2;
            const Type wz = m_data[3] * z2;

            matrix[0][0] = (Type)1.0 - (yy + zz);
            matrix[1][0] = xy - wz;
            matrix[2][0] = xz + wy;
            matrix[3][0] = (Type)0.0;

            matrix[0][1] = xy + wz;
            matrix[1][1] = (Type)1.0 - (xx + zz);
            matrix[2][1] = yz - wx;
            matrix[3][1] = (Type)0.0;

            matrix[0][2] = xz - wy;
            matrix[1][2] = yz + wx;
            matrix[2][2] = (Type)1.0 - (xx + yy);
            matrix[3][2] = (Type)0.0;

            matrix[0][3] = (Type)0.0;
            matrix[1][3] = (Type)0.0;
            matrix[2][3] = (Type)0.0;
            matrix[3][3] = (Type)1.0;

            return matrix;
        }

        //! Invert the rotation. Returns reference to self.
        quaternion<Type>& invert()
        {
            // Optimize by doing 1 div and 4 muls instead of 4 divs.
            const Type inv = 1.0 / m_data.sqrLength();

            m_data[0] = -m_data[0] * inv;
            m_data[1] = -m_data[1] * inv;
            m_data[2] = -m_data[2] * inv;
            m_data[3] = m_data[3] * inv;

            return *this;
        }

        //! Non-destructively inverses the rotation and returns the result.
        quaternion<Type> inverse() const
        {
            quaternion<Type> quat(*this);

            return quat.invert();
        }

        //! Multiplies the quaternions.
        quaternion<Type>& operator*=(const quaternion<Type>& a_quat)
        {
            const Type tx = m_data[0];
            const Type ty = m_data[1];
            const Type tz = m_data[2];
            const Type tw = m_data[3];

            const Type qx = a_quat.m_data[0];
            const Type qy = a_quat.m_data[1];
            const Type qz = a_quat.m_data[2];
            const Type qw = a_quat.m_data[3];

            setValue(qw * tx + qx * tw + qy * tz - qz * ty,
                     qw * ty - qx * tz + qy * tw + qz * tx,
                     qw * tz + qx * ty - qy * tx + qz * tw,
                     qw * tw - qx * tx - qy * ty - qz * tz);

            return *this;
        }

        //! Multiplies components of quaternion with scalar value \a a_scale.
        quaternion<Type>& operator*=(Type a_scale)
        {
            m_data *= a_scale;

            return *this;
        }

        //! Check the two given quaternion for equality.
        friend bool operator==(const quaternion<Type>& q1, const quaternion<Type>& q2)
        {
            return ((q1.m_data[0] == q2.m_data[0]) && (q1.m_data[1] == q2.m_data[1]) && (q1.m_data[2] == q2.m_data[2]) && (q1.m_data[3] == q2.m_data[3]));
        }

        //! Check the two given quaternion for inequality.
        friend bool operator!=(const quaternion<Type>& q1, const quaternion<Type>& q2)
        {
            return !(q1 == q2);
        }

        //! Check for equality with given tolerance.
        bool equals(const quaternion<Type>& a_quat, Type a_tolerance = 1E-2) const
        {
            return m_data.equals(a_quat.m_data, a_tolerance);
        }

        //! Multiplies the two quaternions and returns the result.
        friend quaternion<Type> operator*(const quaternion<Type>& q1, const quaternion<Type>& q2)
        {
            quaternion<Type> q(q1);
            q *= q2;
            return q;
        }

        friend vec3<Type> operator*(const quaternion<Type>& q1, const vec3<Type>& v1)
        {
            vec3<Type> vec;
            q1.multVec(v1, vec);
            return vec;
        }

        //! Rotate the \a src vector and put the result in \a dst.
        void multVec(const vec3<Type>& src, vec3<Type>& dst) const
        {
            const Type QwQx = m_data[3] * m_data[0];
            const Type QwQy = m_data[3] * m_data[1];
            const Type QwQz = m_data[3] * m_data[2];
            const Type QxQy = m_data[0] * m_data[1];
            const Type QxQz = m_data[0] * m_data[2];
            const Type QyQz = m_data[1] * m_data[2];

            const Type Vx = src[0], Vy = src[1], Vz = src[2];

            const Type rx = 2 * (Vy * (-QwQz + QxQy) + Vz * (QwQy + QxQz));
            const Type ry = 2 * (Vx * (QwQz + QxQy) + Vz * (-QwQx + QyQz));
            const Type rz = 2 * (Vx * (-QwQy + QxQz) + Vy * (QwQx + QyQz));

            const Type QwQw = m_data[3] * m_data[3];
            const Type QxQx = m_data[0] * m_data[0];
            const Type QyQy = m_data[1] * m_data[1];
            const Type QzQz = m_data[2] * m_data[2];

            rx += Vx * (QwQw + QxQx - QyQy - QzQz);
            ry += Vy * (QwQw - QxQx + QyQy - QzQz);
            rz += Vz * (QwQw - QxQx - QyQy + QzQz);

            dst.setValue(rx, ry, rz);
        }

        //! Scale the angle of rotation by \a a_scale_factor.
        void scaleAngle(Type a_scale_factor)
        {
            vec3<Type> axis;
            Type deg;

            getValue(axis, deg);
            setValue(axis, deg * a_scale_factor);
        }

        //! Interpolates along the shortest path between the two rotation positions (from \a rot0 to \a rot1, 0<= t <= 1).
        static quaternion<Type> slerp(const quaternion<Type>& rot0, const quaternion<Type>& rot1, Type t)
        {
            quaternion<Type> to = rot1;

            Type dot = rot0.m_data.dot(to.m_data);

            // Find the correct direction of the interpolation.
            if (dot < 0.0f) {
                dot = -dot;
                to.m_data = -to.m_data;
            }

            // fallback to linear interpolation, in case we run out of floating
            // point precision
            Type scale0 = (Type)(1.0 - t);
            Type scale1 = t;

            if ((1.0f - dot) > 1E-3) {
                Type angle = (Type)std::acos(dot);
                Type sinangle = (Type)sin(angle);
                if (sinangle > 1E-3) {
                    // calculate spherical interpolation
                    scale0 = (Type)(std::sin((1.0 - t) * angle)) / sinangle;
                    scale1 = (Type)(std::sin(t * angle)) / sinangle;
                }
            }

            const vec4<Type> vec = (scale0 * rot0.m_data) + (scale1 * to.m_data);

            return quaternion<Type>(vec[0], vec[1], vec[2], vec[3]);
        }

        //! Returns a null rotation.
        static quaternion<Type> identity()
        {
            return quaternion<Type>(0.0, 0.0, 0.0, 1.0);
        }

    private:
        vec4<Type> m_data;
    };

    typedef quaternion<float> quaternionf;
    typedef quaternion<double> quaterniond;
} // namespace gtl
