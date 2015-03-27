#include <gtest/gtest.h>
#include <gtl/complex.hpp>

using namespace gtl;

TEST(ComplexTest, contructors)
{
    complexf c1(1,1);
    complexf c2(1,1);

    EXPECT_TRUE(c1 == c2);

    EXPECT_FLOAT_EQ(1.0f, c1.getReal());
    EXPECT_FLOAT_EQ(1.0f, c1.getImaginary());

    c1.normalize();
    
    EXPECT_FLOAT_EQ(1.0f, c1.modulus());

    c2 = c1;
    c1.negate();
    
    EXPECT_TRUE(c1.equals(-c2));

    c1 = c1 + c2;
    c1 += c2;
    c1 = c1 - c2;
    c1 -= c2;
    c1 = c1 / c2;
    c1 /= c2;
    c1 = c1 * c2;
    c1 *= c2;
    
    c2 = c1;
    c1 /= c1;
    c2 = c2 / c2;

    EXPECT_TRUE(c1.equals(c2));

    c1 += 10.0f;
    c1 -= 10.0f;
    c1 *= 10.0f;
    c1 /= 10.0f;

    c1.modulus();
    c1.sqrModulus();

    complexf c3 = c1.getConjugate();

    c3.normalize();
    c3 = c3 * 2.0f;
    c3 = c3 / 2.0f;
}
