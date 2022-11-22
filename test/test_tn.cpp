#include <tinynurbs/tinynurbs.h>
#include <iostream>

int main()
{
    tinynurbs::RationalSurface<float> srf;
    srf.degree_u = 3;
    srf.degree_v = 3;
    srf.knots_u = {0, 0, 0, 0, 1, 1, 1, 1};
    srf.knots_v = {0, 0, 0, 0, 1, 1, 1, 1};

    // 2D array of control points using tinynurbs::array2<T> container
    // Example from geometrictools.com/Documentation/NURBSCircleSphere.pdf
    srf.control_points = {4, 4, 
                        {glm::vec3(0, 0, 1), glm::vec3(0, 0, 1), glm::vec3(0, 0, 1), glm::vec3(0, 0, 1),
                        glm::vec3(2, 0, 1), glm::vec3(2, 4, 1),  glm::vec3(-2, 4, 1),  glm::vec3(-2, 0, 1),
                        glm::vec3(2, 0, -1), glm::vec3(2, 4, -1), glm::vec3(-2, 4, -1), glm::vec3(-2, 0, -1),
                        glm::vec3(0, 0, -1), glm::vec3(0, 0, -1), glm::vec3(0, 0, -1), glm::vec3(0, 0, -1)
                        }
                        };
    srf.weights = {4, 4,
                {1,       1.f/3.f, 1.f/3.f, 1,
                1.f/3.f, 1.f/9.f, 1.f/9.f, 1.f/3.f,
                1.f/3.f, 1.f/9.f, 1.f/9.f, 1.f/3.f,
                1,       1.f/3.f, 1.f/3.f, 1
                }
                };

    tinynurbs::RationalSurface<float> left, right;
    std::tie(left, right) = tinynurbs::surfaceSplitV<float>(srf, 0.25);

    tinynurbs::surfaceSaveOBJ("output_surface.obj", srf);

    return 0;
}