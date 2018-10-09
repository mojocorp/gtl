#include <gtest/gtest.h>
#include <gtl/vec3.hpp>

using namespace gtl;

TEST(Vec3Test, sanity)
{
    EXPECT_EQ(3 * sizeof(int), sizeof(vec3i));
    EXPECT_EQ(3 * sizeof(float), sizeof(vec3f));
    EXPECT_EQ(3 * sizeof(double), sizeof(vec3d));
}

TEST(Vec3Test, contructors)
{
    vec3f v3f1(1.0f, 0.0f, 0.0f);
    vec3f v3f2(v3f1);

    EXPECT_TRUE(v3f1 == v3f2);

    EXPECT_TRUE(v3f1 + v3f2 == 2.0f * v3f1);

    EXPECT_TRUE(v3f1 - v3f2 == vec3f(0.0f, 0.0f, 0.0f));

    v3f1.setValue(1.0f, 0.0f, 0.0f);

    EXPECT_FLOAT_EQ(1.0f, v3f1.length());

    v3f1.normalize();
    v3f2.normalize();

    EXPECT_FLOAT_EQ(1.0f, v3f1.length());

    EXPECT_TRUE(v3f1.cross(v3f2) == vec3f(0.0f, 0.0f, 0.0f));

    EXPECT_EQ(v3f1.x, v3f1[0]);
    EXPECT_EQ(v3f1.y, v3f1[1]);
    EXPECT_EQ(v3f1.z, v3f1[2]);

    v3f1 /= 4.0f;

    v3f1.negate();
    v3f1.getValue();
}
