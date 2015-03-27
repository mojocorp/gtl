#include <gtest/gtest.h>
#include <gtl/box2.hpp>

using namespace gtl;

TEST(Box2Test, contructors)
{
    box2f box2f1, box2f2;

    EXPECT_TRUE(box2f1.isEmpty());

    vec2f min(-1.0f,-1.0f);
    vec2f max( 1.0f, 1.0f);

    box2f1.setBounds(min, max);

    EXPECT_TRUE(box2f1.getMin() == min);
    EXPECT_TRUE(box2f1.getMax() == max);
 
    EXPECT_TRUE(box2f1.getSize() == vec2f(2.0f,2.0f));
    EXPECT_TRUE(box2f1.getCenter() == vec2f(0.0f,0.0f));

    box2f2.extendBy(min);
    box2f2.extendBy(max);

    EXPECT_TRUE(box2f2.getMin() == min);
    EXPECT_TRUE(box2f2.getMax() == max);

    EXPECT_TRUE(box2f1 == box2f2);

    box2f1.makeEmpty();

    EXPECT_TRUE(box2f1.isEmpty());

    box2f1.setBounds(vec2f(-0.5f,-0.5f), vec2f(0.5f,0.5f));

    EXPECT_TRUE(box2f2.intersect(box2f1));

    EXPECT_TRUE(box2f2.intersect(min));
    EXPECT_FALSE(box2f1.intersect(min));

    EXPECT_TRUE(box2f2.intersect(box2f1));

    box2f1.setBounds(vec2f(2.0f,2.0f), vec2f(3.0f,3.0f));

    EXPECT_FALSE(box2f2.intersect(box2f1));
}
