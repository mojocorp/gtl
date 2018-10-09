#include <gtest/gtest.h>
#include <gtl/ray.hpp>

using namespace gtl;

TEST(RayTest, contructors)
{
    vec3f origin(0.0f, 0.0f, 0.0f);
    vec3f direction(1.0f, 0.0f, 0.0f);
    vec3f middle(0.5f, 0.0f, 0.0f);

    rayf ray(origin, direction);

    EXPECT_TRUE(ray.getDirection() == direction);
    EXPECT_TRUE(ray.getOrigin() == origin);

    ray.setValue(origin, direction);

    EXPECT_TRUE(ray.getDirection() == direction);
    EXPECT_TRUE(ray.getOrigin() == origin);

    vec3f point(0.0f, 1.0f, 0.0f);

    EXPECT_FLOAT_EQ(1.0f, ray.getDistance(point));

    point.setValue(0.5f, 1.0f, 0.0f);

    EXPECT_TRUE(ray.project(point) == middle);

    EXPECT_TRUE(ray.getValue(0.5f) == middle);

    rayf rayf2(ray);

    EXPECT_TRUE(ray == rayf2);

    vec3f origin2(0.5f, 1.0f, -0.5f);
    vec3f direction2(0.0f, 0.0f, 1.0f);
    vec3f middle2(0.5f, 1.0f, 0.0f);

    rayf2.setValue(origin2, direction2);

    float mua, mub;
    EXPECT_TRUE(rayf2.intersect(ray, mua, mub));

    EXPECT_FLOAT_EQ(0.5f, mua);
    EXPECT_FLOAT_EQ(0.5f, mub);
}
