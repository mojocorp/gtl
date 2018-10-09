#include <gtest/gtest.h>
#include <gtl/vec4.hpp>

using namespace gtl;

TEST(Vec4Test, sanity)
{
    EXPECT_EQ(4 * sizeof(int), sizeof(vec4i));
    EXPECT_EQ(4 * sizeof(float), sizeof(vec4f));
    EXPECT_EQ(4 * sizeof(double), sizeof(vec4d));
}

TEST(Vec4Test, contructors)
{
    vec4f v4f1(1.0f, 0.0f, 0.0f, 0.0f);
    vec4f v4f2(v4f1);

    EXPECT_TRUE(v4f1 == v4f2);

    EXPECT_TRUE(v4f1 + v4f2 == 2.0f * v4f1);

    EXPECT_TRUE(v4f1 - v4f2 == vec4f(0.0f, 0.0f, 0.0f, 0.0f));

    v4f1.setValue(1.0f, 0.0f, 0.0f, 0.0f);

    EXPECT_FLOAT_EQ(1.0f, v4f1.length());

    v4f1.normalize();
    v4f2.normalize();

    EXPECT_FLOAT_EQ(1.0f, v4f1.length());

    EXPECT_EQ(v4f1.x, v4f1[0]);
    EXPECT_EQ(v4f1.y, v4f1[1]);
    EXPECT_EQ(v4f1.z, v4f1[2]);

    v4f2.negate();
    v4f2.getValue();

    float v[4];
    vec4f v4f22(v);
}
