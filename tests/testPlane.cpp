#include <gtest/gtest.h>
#include <gtl/plane.hpp>

using namespace gtl;

TEST(PlaneTest, contructors)
{
    vec3f normal(0.0f,1.0f,0.0f);
    float distance = 1.0f;

    planef plane(normal, distance);

    EXPECT_TRUE(plane.getNormal() == normal);
    EXPECT_FLOAT_EQ(distance, plane.getDistanceFromOrigin());
 
    plane.setNormal(normal);
    plane.setDistance(distance);

    EXPECT_TRUE(plane.getNormal() == normal);
    EXPECT_TRUE(plane.getDistanceFromOrigin() == distance);
 
    EXPECT_FLOAT_EQ(-1.0f, plane.getDistance(vec3f(0.0f,0.0f,0.0f)) );

    EXPECT_FALSE(plane.isInHalfSpace(vec3f(0.0f,0.0f,0.0f)));

    plane = planef(vec3f(0.0f,1.0f,0.0f),vec3f(2.0f,1.0f,0.0f),vec3f(1.0f,1.0f,-1.0f));

    EXPECT_TRUE(plane.getNormal() == normal);
    EXPECT_FLOAT_EQ(distance, plane.getDistanceFromOrigin());
 
    planef plane2 = plane;

    EXPECT_TRUE(plane2 == plane);

    rayf ray(vec3f(0.0f,0.0f,0.0f), vec3f(0.0f,1.0f,0.0));

    float mu;
    EXPECT_TRUE(plane.intersect(ray, mu));

    EXPECT_FLOAT_EQ(0.0f, plane.getDistance(ray.getValue(mu)));


    ////////////
    std::vector< vec3f > points;

    points.push_back(vec3f(0,0,0));
    points.push_back(vec3f(1,1,1));
    points.push_back(vec3f(2,2,2));

    planef plane3(points);

    vec3f p = plane3.project(vec3f(0.0f,1.0f,0.0));
}
