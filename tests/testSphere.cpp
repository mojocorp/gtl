#include <gtest/gtest.h>
#include <gtl/sphere.hpp>

using namespace gtl;

TEST(SphereTest, contructors)
{
    vec3f center(0.0f, 0.0f, 0.0f);
    float radius = 1.0f;

    spheref sphere(center, radius);

    EXPECT_TRUE(sphere.getCenter() == center);
    EXPECT_FLOAT_EQ(radius, sphere.getRadius());

    sphere.setValue(center, radius);

    EXPECT_TRUE(sphere.getCenter() == center);
    EXPECT_FLOAT_EQ(radius, sphere.getRadius());

    sphere.setCenter(center);
    sphere.setRadius(radius);

    EXPECT_TRUE(sphere.getCenter() == center);
    EXPECT_FLOAT_EQ(radius, sphere.getRadius());

    sphere.setPoles(vec3f(0.0f, 1.0f, 0.0f), vec3f(0.0f, -1.0f, 0.0f));

    EXPECT_TRUE(sphere.getCenter() == center);
    EXPECT_FLOAT_EQ(radius, sphere.getRadius());

    EXPECT_TRUE(sphere.intersect(center));
    EXPECT_FALSE(sphere.intersect(vec3f(0.0f, 2.0f, 0.0f)));

    spheref sphere2(sphere);

    EXPECT_TRUE(sphere == sphere2);

    float c = 2.0f * sqrt(1.0f / 3.0f);

    box3f box(vec3f(0.0f, 0.0f, 0.0f), vec3f(c, c, c));

    sphere2.circumscribe(box);

    EXPECT_FLOAT_EQ(1.0f, sphere2.getRadius());
}
