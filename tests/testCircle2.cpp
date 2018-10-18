#include <gtest/gtest.h>
#include <gtl/circle2.hpp>

using namespace gtl;

TEST(Circle2Test, contructors)
{
    circle2f c1(vec2f(-1.0f, -1.0f), 1.5f);
    circle2f c2(vec2f(1.0f, 1.0f), 2.0f);

    vec2f p1, p2;
    EXPECT_TRUE(c1.intersect(c2, p1, p2));

    EXPECT_FLOAT_EQ(7.0685835f, c1.getArea());
    EXPECT_FLOAT_EQ(9.424778f, c1.getCircumference());
}
