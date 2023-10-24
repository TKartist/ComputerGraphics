//
//  Textures.h
//  Raytracer
//
//  Created by Piotr Didyk on 19.08.21.
//

#ifndef Textures_h
#define Textures_h

#include "glm/glm.hpp"

glm::vec3 checkerboardTexture(glm::vec2 uv)
{
    int n = 20;
    int f = int(floor(n * uv.x) + floor(n * uv.y)) % 2;
    if (f == 0)
    {
        std::cout << "1" << std::endl;
        return glm::vec3(1.0);
    }
    else
    {
        std::cout << "0" << std::endl;
        return glm::vec3(0.0);
    }
}
glm::vec3 rainbowTexture(glm::vec2 uv)
{
    int n = 20;
    int f = int(int(n * (uv.x + uv.y))) % 3;
    if (f == 0)
    {
        return glm::vec3(1, 0, 0);
    }
    if (f == 1)
    {
        return glm::vec3(0, 1, 0);
    }
    return glm::vec3(0, 0, 1);
}

#endif /* Textures_h */
