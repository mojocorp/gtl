#include <gtest/gtest.h>
#include <gtl/matrix4.hpp>

using namespace gtl;

TEST(Matrix4Test, contructors)
{
    matrix4f mat4f1(1.0f, 0.0f, 0.0f, 0.0f,
                    0.0f, 1.0f, 0.0f, 0.0f,
                    0.0f, 0.0f, 1.0f, 0.0f,
                    0.0f, 0.0f, 0.0f, 1.0f);

    EXPECT_TRUE(mat4f1.isIdentity());

    matrix4f mat4f2;
    mat4f2.makeIdentity();

    EXPECT_TRUE(mat4f1 == mat4f2);
    EXPECT_TRUE(mat4f1.equals(mat4f2));
}
