#include <gtest/gtest.h>
#include <gtl/matrix3.hpp>
#include <gtl/quaternion.hpp>

using namespace gtl;

TEST(Matrix3Test, contructors)
{
    matrix3f mat3f1(1.0f, 0.0f, 0.0f,
                    0.0f, 1.0f, 0.0f,
                    0.0f, 0.0f, 1.0f);

    EXPECT_TRUE(mat3f1.isIdentity());

    matrix3f mat3f2;
    mat3f2.makeIdentity();

    EXPECT_TRUE(mat3f1 == mat3f2);
    EXPECT_TRUE(mat3f1.equals(mat3f2));

    mat3f1.inverse();

    quaternionf rot;
    mat3f1.setRotate(rot);
}
