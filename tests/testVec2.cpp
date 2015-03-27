#include <gtest/gtest.h>
#include <gtl/vec2.hpp>

using namespace gtl;

TEST(Vec2Test, sanity)
{
    EXPECT_EQ(2*sizeof(int),    sizeof(vec2i));
    EXPECT_EQ(2*sizeof(float),  sizeof(vec2f));
    EXPECT_EQ(2*sizeof(double), sizeof(vec2d));
}

TEST(Vec2Test, contructors)
{
    vec2f v2f1(1.0f, 0.0f);
    vec2f v2f2(v2f1);

    EXPECT_TRUE(v2f1 == v2f2);

    EXPECT_TRUE(v2f1 + v2f2 == 2.0f*v2f1);

    EXPECT_TRUE(v2f1 - v2f2 == vec2f(0.0f,0.0f));

    v2f1.setValue(1.0f,0.0f);

    EXPECT_FLOAT_EQ(1.0f, v2f1.length());

    v2f1.normalize();
    v2f2.normalize();

    EXPECT_FLOAT_EQ(1.0f, v2f1.length());

    EXPECT_FLOAT_EQ(0.0f, v2f1.cross(v2f2));

    EXPECT_EQ(v2f1.x, v2f1[0]);
    EXPECT_EQ(v2f1.y, v2f1[1]);
}

