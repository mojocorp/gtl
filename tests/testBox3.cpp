#include <gtest/gtest.h>
#include <gtl/box3.hpp>
#include <gtl/xfbox3.hpp>

#include <cstdlib>
#include <ctime>

using namespace gtl;

TEST(Box3Test, contructors)
{
    box3f box3f1, box3f2;

    EXPECT_TRUE(box3f1.isEmpty());

    vec3f min(-1.0f, -1.0f, -1.0f);
    vec3f max(1.0f, 1.0f, 1.0f);

    box3f1.setBounds(min, max);

    EXPECT_TRUE(box3f1.getMin() == min);
    EXPECT_TRUE(box3f1.getMax() == max);

    EXPECT_TRUE(box3f1.getSize() == vec3f(2.0f, 2.0f, 2.0f));
    EXPECT_TRUE(box3f1.getCenter() == vec3f(0.0f, 0.0f, 0.0f));

    box3f2.extendBy(min);
    box3f2.extendBy(max);

    EXPECT_TRUE(box3f2.getMin() == min);
    EXPECT_TRUE(box3f2.getMax() == max);

    EXPECT_TRUE(box3f1 == box3f2);

    box3f1.makeEmpty();

    EXPECT_TRUE(box3f1.isEmpty());

    box3f1.setBounds(vec3f(-0.5f, -0.5f, -0.5f), vec3f(0.5f, 0.5f, 0.5f));

    EXPECT_TRUE(box3f2.intersect(box3f1));

    EXPECT_TRUE(box3f2.intersect(min));
    EXPECT_FALSE(box3f1.intersect(min));

    EXPECT_TRUE(box3f2.intersect(box3f1));

    box3f1.setBounds(vec3f(2.0f, 2.0f, 2.0f), vec3f(3.0f, 3.0f, 3.0f));

    EXPECT_FALSE(box3f2.intersect(box3f1));
}

TEST(XfBox3Test, constructors)
{
    xfbox3f xfbox;
}

bool triBoxOverlapMoller(const gtl::box3f& box, const gtl::vec3f triverts[3]);

static gtl::vec3f getRandomVector(float min, float max)
{
    return gtl::vec3f(min + (std::rand() / (float)RAND_MAX) * (max - min),
                      min + (std::rand() / (float)RAND_MAX) * (max - min),
                      min + (std::rand() / (float)RAND_MAX) * (max - min));
}

TEST(Box3Test, triBoxOverlapCompare)
{
    gtl::vec3f boxcenter(0.0f, 0.0f, 0.0f);
    gtl::vec3f boxhalfsize(0.5f, 0.5f, 0.5f);

    gtl::box3f box(gtl::vec3f(boxcenter - boxhalfsize),
                   gtl::vec3f(boxcenter + boxhalfsize));
    gtl::vec3f triverts[3];

    std::srand(std::time(nullptr));
    for (int i = 0; i < 100000; i++) {
        triverts[0] = getRandomVector(-4.0f, 4.0f);
        triverts[1] = getRandomVector(-4.0f, 4.0f);
        triverts[2] = getRandomVector(-4.0f, 4.0f);

        bool r1 = triBoxOverlapMoller(box, triverts);
        bool r2 = box.intersect(triverts[0], triverts[1], triverts[2]);
        EXPECT_EQ(r1, r2);
    }
}

TEST(Box3Test, triBoxOverlap_moller)
{
    gtl::vec3f boxcenter(0.0f, 0.0f, 0.0f);
    gtl::vec3f boxhalfsize(0.5f, 0.5f, 0.5f);

    gtl::box3f box(gtl::vec3f(boxcenter - boxhalfsize),
                   gtl::vec3f(boxcenter + boxhalfsize));
    gtl::vec3f triverts[3];

    std::srand(0);
    for (int i = 0; i < 100000; i++) {
        triverts[0] = getRandomVector(-4.0f, 4.0f);
        triverts[1] = getRandomVector(-4.0f, 4.0f);
        triverts[2] = getRandomVector(-4.0f, 4.0f);

        triBoxOverlapMoller(box, triverts);
    }
}

TEST(Box3Test, triBoxOverlap_gtl)
{
    gtl::vec3f boxcenter(0.0f, 0.0f, 0.0f);
    gtl::vec3f boxhalfsize(0.5f, 0.5f, 0.5f);

    gtl::box3f box(gtl::vec3f(boxcenter - boxhalfsize),
                   gtl::vec3f(boxcenter + boxhalfsize));
    gtl::vec3f triverts[3];

    std::srand(0);
    for (int i = 0; i < 100000; i++) {
        triverts[0] = getRandomVector(-4.0f, 4.0f);
        triverts[1] = getRandomVector(-4.0f, 4.0f);
        triverts[2] = getRandomVector(-4.0f, 4.0f);

        box.intersect(triverts[0], triverts[1], triverts[2]);
    }
}

/********************************************************/
/* AABB-triangle overlap test code                      */
/* by Tomas Akenine-Möller                              */
/* Function: int triBoxOverlap(float boxcenter[3],      */
/*          float boxhalfsize[3],float triverts[3][3]); */
/* History:                                             */
/*   2001-03-05: released the code in its first version */
/*   2001-06-18: changed the order of the tests, faster */
/*                                                      */
/* Acknowledgement: Many thanks to Pierre Terdiman for  */
/* suggestions and discussions on how to optimize code. */
/* Thanks to David Hunt for finding a ">="-bug!         */
/********************************************************/
#include <cmath>
#include <cstdio>

#define X 0
#define Y 1
#define Z 2

#define CROSS(dest, v1, v2)                  \
    dest[0] = v1[1] * v2[2] - v1[2] * v2[1]; \
    dest[1] = v1[2] * v2[0] - v1[0] * v2[2]; \
    dest[2] = v1[0] * v2[1] - v1[1] * v2[0];

#define DOT(v1, v2) (v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2])

#define SUB(dest, v1, v2)    \
    dest[0] = v1[0] - v2[0]; \
    dest[1] = v1[1] - v2[1]; \
    dest[2] = v1[2] - v2[2];

#define FINDMINMAX(x0, x1, x2, min, max) \
    min = max = x0;                      \
    if (x1 < min)                        \
        min = x1;                        \
    if (x1 > max)                        \
        max = x1;                        \
    if (x2 < min)                        \
        min = x2;                        \
    if (x2 > max)                        \
        max = x2;

static int planeBoxOverlap(float normal[3], float vert[3], float maxbox[3]) // -NJMP-
{
    int q;
    float vmin[3], vmax[3], v;
    for (q = X; q <= Z; q++) {
        v = vert[q]; // -NJMP-
        if (normal[q] > 0.0f) {
            vmin[q] = -maxbox[q] - v; // -NJMP-
            vmax[q] = maxbox[q] - v; // -NJMP-
        } else {
            vmin[q] = maxbox[q] - v; // -NJMP-
            vmax[q] = -maxbox[q] - v; // -NJMP-
        }
    }
    if (DOT(normal, vmin) > 0.0f)
        return 0; // -NJMP-
    if (DOT(normal, vmax) >= 0.0f)
        return 1; // -NJMP-
    return 0;
}

/*======================== X-tests ========================*/
#define AXISTEST_X01(a, b, fa, fb)                   \
    p0 = a * v0[Y] - b * v0[Z];                      \
    p2 = a * v2[Y] - b * v2[Z];                      \
    if (p0 < p2) {                                   \
        min = p0;                                    \
        max = p2;                                    \
    } else {                                         \
        min = p2;                                    \
        max = p0;                                    \
    }                                                \
    rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z]; \
    if (min > rad || max < -rad)                     \
        return 0;

#define AXISTEST_X2(a, b, fa, fb)                    \
    p0 = a * v0[Y] - b * v0[Z];                      \
    p1 = a * v1[Y] - b * v1[Z];                      \
    if (p0 < p1) {                                   \
        min = p0;                                    \
        max = p1;                                    \
    } else {                                         \
        min = p1;                                    \
        max = p0;                                    \
    }                                                \
    rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z]; \
    if (min > rad || max < -rad)                     \
        return 0;

/*======================== Y-tests ========================*/
#define AXISTEST_Y02(a, b, fa, fb)                   \
    p0 = -a * v0[X] + b * v0[Z];                     \
    p2 = -a * v2[X] + b * v2[Z];                     \
    if (p0 < p2) {                                   \
        min = p0;                                    \
        max = p2;                                    \
    } else {                                         \
        min = p2;                                    \
        max = p0;                                    \
    }                                                \
    rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z]; \
    if (min > rad || max < -rad)                     \
        return 0;

#define AXISTEST_Y1(a, b, fa, fb)                    \
    p0 = -a * v0[X] + b * v0[Z];                     \
    p1 = -a * v1[X] + b * v1[Z];                     \
    if (p0 < p1) {                                   \
        min = p0;                                    \
        max = p1;                                    \
    } else {                                         \
        min = p1;                                    \
        max = p0;                                    \
    }                                                \
    rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z]; \
    if (min > rad || max < -rad)                     \
        return 0;

/*======================== Z-tests ========================*/
#define AXISTEST_Z12(a, b, fa, fb)                   \
    p1 = a * v1[X] - b * v1[Y];                      \
    p2 = a * v2[X] - b * v2[Y];                      \
    if (p2 < p1) {                                   \
        min = p2;                                    \
        max = p1;                                    \
    } else {                                         \
        min = p1;                                    \
        max = p2;                                    \
    }                                                \
    rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y]; \
    if (min > rad || max < -rad)                     \
        return 0;

#define AXISTEST_Z0(a, b, fa, fb)                    \
    p0 = a * v0[X] - b * v0[Y];                      \
    p1 = a * v1[X] - b * v1[Y];                      \
    if (p0 < p1) {                                   \
        min = p0;                                    \
        max = p1;                                    \
    } else {                                         \
        min = p1;                                    \
        max = p0;                                    \
    }                                                \
    rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y]; \
    if (min > rad || max < -rad)                     \
        return 0;

bool triBoxOverlapMoller(const gtl::box3f& box, const gtl::vec3f triverts[3])
{
    /*    use separating axis theorem to test overlap between triangle and box */
    /*    need to test for overlap in these directions: */
    /*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
    /*       we do not even need to test these) */
    /*    2) normal of the triangle */
    /*    3) crossproduct(edge from tri, {x,y,z}-directin) */
    /*       this gives 3x3=9 more tests */

    float min, max, p0, p1, p2, rad, fex, fey, fez; // -NJMP- "d" local variable removed
    float normal[3];

    const gtl::vec3f boxcenter = box.getCenter();
    const gtl::vec3f boxsize = box.getSize() * 0.5f;
    float boxhalfsize[3] = { boxsize[0], boxsize[1], boxsize[2] };

    float v0[3], v1[3], v2[3];
    float e0[3], e1[3], e2[3];
    /* This is the fastest branch on Sun */
    /* move everything so that the boxcenter is in (0,0,0) */
    SUB(v0, triverts[0], boxcenter);
    SUB(v1, triverts[1], boxcenter);
    SUB(v2, triverts[2], boxcenter);

    /* compute triangle edges */
    SUB(e0, v1, v0); /* tri edge 0 */
    SUB(e1, v2, v1); /* tri edge 1 */
    SUB(e2, v0, v2); /* tri edge 2 */

    /* Bullet 3:  */

    /*  test the 9 tests first (this was faster) */
    fex = fabsf(e0[X]);
    fey = fabsf(e0[Y]);
    fez = fabsf(e0[Z]);

    AXISTEST_X01(e0[Z], e0[Y], fez, fey);
    AXISTEST_Y02(e0[Z], e0[X], fez, fex);
    AXISTEST_Z12(e0[Y], e0[X], fey, fex);

    fex = fabsf(e1[X]);
    fey = fabsf(e1[Y]);
    fez = fabsf(e1[Z]);

    AXISTEST_X01(e1[Z], e1[Y], fez, fey);
    AXISTEST_Y02(e1[Z], e1[X], fez, fex);
    AXISTEST_Z0(e1[Y], e1[X], fey, fex);

    fex = fabsf(e2[X]);
    fey = fabsf(e2[Y]);
    fez = fabsf(e2[Z]);

    AXISTEST_X2(e2[Z], e2[Y], fez, fey);
    AXISTEST_Y1(e2[Z], e2[X], fez, fex);
    AXISTEST_Z12(e2[Y], e2[X], fey, fex);

    /* Bullet 1: */

    /*  first test overlap in the {x,y,z}-directions */
    /*  find min, max of the triangle each direction, and test for overlap in */
    /*  that direction -- this is equivalent to testing a minimal AABB around */
    /*  the triangle against the AABB */

    /* test in X-direction */
    FINDMINMAX(v0[X], v1[X], v2[X], min, max);

    if (min > boxhalfsize[X] || max < -boxhalfsize[X])
        return 0;

    /* test in Y-direction */
    FINDMINMAX(v0[Y], v1[Y], v2[Y], min, max);

    if (min > boxhalfsize[Y] || max < -boxhalfsize[Y])
        return 0;

    /* test in Z-direction */
    FINDMINMAX(v0[Z], v1[Z], v2[Z], min, max);

    if (min > boxhalfsize[Z] || max < -boxhalfsize[Z])
        return 0;

    /* Bullet 2: */
    /*  test if the box intersects the plane of the triangle */
    /*  compute plane equation of triangle: normal*x+d=0 */
    CROSS(normal, e0, e1);

    if (!planeBoxOverlap(normal, v0, boxhalfsize))
        return 0;

    return 1; /* box and triangle overlaps */
}
