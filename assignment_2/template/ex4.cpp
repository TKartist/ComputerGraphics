/**
@file main.cpp
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include "glm/glm.hpp"

#include "Image.h"
#include "Material.h"

using namespace std;

/**
 Class representing a single ray.
 */
class Ray
{
public:
    glm::vec3 origin;    ///< Origin of the ray
    glm::vec3 direction; ///< Direction of the ray
                         /**
                          Contructor of the ray
                          @param origin Origin of the ray
                          @param direction Direction of the ray
                          */
    Ray(glm::vec3 origin, glm::vec3 direction) : origin(origin), direction(direction)
    {
    }
};

class Object;

/**
 Structure representing the even of hitting an object
 */
struct Hit
{
    bool hit;               ///< Boolean indicating whether there was or there was no intersection with an object
    glm::vec3 normal;       ///< Normal vector of the intersected object at the intersection point
    glm::vec3 intersection; ///< Point of Intersection
    float distance;         ///< Distance from the origin of the ray to the intersection point
    Object *object;         ///< A pointer to the intersected object
};

/**
 General class for the object
 */
class Object
{
public:
    glm::vec3 color;   ///< Color of the object
    Material material; ///< Structure describing the material of the object
                       /** A function computing an intersection, which returns the structure Hit */
    virtual Hit intersect(Ray ray) = 0;

    /** Function that returns the material struct of the object*/
    Material getMaterial()
    {
        return material;
    }
    /** Function that set the material
     @param material A structure desribing the material of the object
    */
    void setMaterial(Material material)
    {
        this->material = material;
    }
};

/**
 Implementation of the class Object for sphere shape.
 */
class Sphere : public Object
{
private:
    float radius;     ///< Radius of the sphere
    glm::vec3 center; ///< Center of the sphere

public:
    /**
     The constructor of the sphere
     @param radius Radius of the sphere
     @param center Center of the sphere
     @param color Color of the sphere
     */
    Sphere(float radius, glm::vec3 center, glm::vec3 color) : radius(radius), center(center)
    {
        this->color = color;
    }
    Sphere(float radius, glm::vec3 center, Material material) : radius(radius), center(center)
    {
        this->material = material;
    }
    /** Implementation of the intersection function*/
    Hit intersect(Ray ray)
    {

        glm::vec3 c = center - ray.origin;

        float cdotc = glm::dot(c, c);
        float cdotd = glm::dot(c, ray.direction);

        Hit hit;

        float D = 0;
        if (cdotc > cdotd * cdotd)
        {
            D = sqrt(cdotc - cdotd * cdotd);
        }
        if (D <= radius)
        {
            hit.hit = true;
            float t1 = cdotd - sqrt(radius * radius - D * D);
            float t2 = cdotd + sqrt(radius * radius - D * D);

            float t = t1;
            if (t < 0)
                t = t2;
            if (t < 0)
            {
                hit.hit = false;
                return hit;
            }

            hit.intersection = ray.origin + t * ray.direction;
            hit.normal = glm::normalize(hit.intersection - center);
            hit.distance = glm::distance(ray.origin, hit.intersection);
            hit.object = this;
        }
        else
        {
            hit.hit = false;
        }
        return hit;
    }
};

/**
 Light class
 */
class Light
{
public:
    glm::vec3 position; ///< Position of the light source
    glm::vec3 color;    ///< Color/intentisty of the light source
    Light(glm::vec3 position) : position(position)
    {
        color = glm::vec3(1.0);
    }
    Light(glm::vec3 position, glm::vec3 color) : position(position), color(color)
    {
    }
};

vector<Light *> lights; ///< A list of lights in the scene
glm::vec3 ambient_light(1.0, 1.0, 1.0);
vector<Object *> objects; ///< A list of all objects in the scene

/** Function for computing color of an object according to the Phong Model
 @param point A point belonging to the object for which the color is computed
 @param normal A normal vector the the point
 @param view_direction A normalized direction from the point to the viewer/camera
 @param material A material structure representing the material of the object
*/
glm::vec3 PhongModel(glm::vec3 point, glm::vec3 normal, glm::vec3 view_direction, Material material)
{

    glm::vec3 color(0.0);

    /*

     Assignment 2

     Phong model.
     Your code should implement a loop over all the lightsources in the array lights and agredate the contribution of each of them to the final color of the object.
     Outside of the loop add also the ambient component from ambient_light.

    */
    for (Light *light : lights)
    {
        // getting diffuse
        glm::vec3 lightSource = glm::normalize(light->position - point);
        float diffuse = max(glm::dot(lightSource, normal), 0.0f);

        // getting specular
        glm::vec3 halfVector = glm::normalize(lightSource + view_direction);
        float specular = glm::pow(max(glm::dot(normal, halfVector), 0.0f), 4 * material.shininess);

        color += light->color * (material.diffuse * diffuse + material.specular * specular);
    }

    color += ambient_light * material.ambient;

    // The final color has to be clamped so the values do not go beyond 0 and 1.
    color = glm::clamp(color, glm::vec3(0.0), glm::vec3(1.0));
    return color;
}

/**
 Functions that computes a color along the ray
 @param ray Ray that should be traced through the scene
 @return Color at the intersection point
 */
glm::vec3 trace_ray(Ray ray)
{

    Hit closest_hit;

    closest_hit.hit = false;
    closest_hit.distance = INFINITY;

    for (int k = 0; k < objects.size(); k++)
    {
        Hit hit = objects[k]->intersect(ray);
        if (hit.hit == true && hit.distance < closest_hit.distance)
            closest_hit = hit;
    }

    glm::vec3 color(0.0);

    if (closest_hit.hit)
    {
        color = PhongModel(closest_hit.intersection, closest_hit.normal, glm::normalize(-ray.direction), closest_hit.object->getMaterial());
        /*

         Assignment 2

         Replace the above line of the code with the call of the function for computing Phong model below.

        */
        // color = PhongModel(closest_hit.intersection, closest_hit.normal, glm::normalize(-ray.direction), closest_hit.object->getMaterial());
    }
    else
    {
        color = glm::vec3(0.0, 0.0, 0.0);
    }
    return color;
}
/**
 Function defining the scene
 */
void sceneDefinition()
{
    Material red_specular;
    red_specular.diffuse = glm::vec3(1.0f, 0.3f, 0.3f);
    red_specular.ambient = glm::vec3(0.01f, 0.03f, 0.03f);
    red_specular.specular = glm::vec3(0.5);
    red_specular.shininess = 10.0;

    Material blue_specular;
    blue_specular.diffuse = glm::vec3(0.7f, 0.7f, 1.0f);
    blue_specular.ambient = glm::vec3(0.07f, 0.07f, 0.1f);
    blue_specular.specular = glm::vec3(0.6);
    blue_specular.shininess = 100.0;

    Material green_specular;
    green_specular.diffuse = glm::vec3(0.7f, 0.9f, 0.7f);
    green_specular.ambient = glm::vec3(0.07f, 0.09f, 0.07f);
    green_specular.specular = glm::vec3(0.0);
    green_specular.shininess = 0.0;

    objects.push_back(new Sphere(0.5, glm::vec3(-1.0, -2.5, 6.0), red_specular));
    objects.push_back(new Sphere(1.0, glm::vec3(1.0, -2.0, 8.0), blue_specular));
    objects.push_back(new Sphere(1.0, glm::vec3(3.0, -2.0, 6.0), green_specular));
    /*

     Assignment 2

     Add here all the objects to the scene. Remember to add them using the new constructor for the sphere with material structure.
     You will also need to define the materials.
     Example of adding one sphere:

     Material red_specular;
     red_specular.diffuse = glm::vec3(1.0f, 0.3f, 0.3f);
     red_specular.ambient = glm::vec3(0.01f, 0.03f, 0.03f);
     red_specular.specular = glm::vec3(0.5);
     red_specular.shininess = 10.0;

     objects.push_back(new Sphere(0.5, glm::vec3(-1,-2.5,6), red_specular));


     Remember also about adding some lights. For example a white light of intensity 0.4 and position in (0,26,5):

     lights.push_back(new Light(glm::vec3(0, 26, 5), glm::vec3(0.4)));

    */

    lights.push_back(new Light(glm::vec3(0.0, 26.0, 5.0), glm::vec3(0.4)));
    lights.push_back(new Light(glm::vec3(0.0, 1.0, 12.0), glm::vec3(0.4)));
    lights.push_back(new Light(glm::vec3(0.0, 5.0, 1.0), glm::vec3(0.4)));
}

int main(int argc, const char *argv[])
{

    clock_t t = clock(); // variable for keeping the time of the rendering

    int width = 1024; // width of the image
    int height = 768; // height of the image
    float fov = 90;   // field of view

    sceneDefinition(); // Let's define a scene

    if (argc > 2)
    {
        float angle = 2 * M_PI / 60 * stoi(argv[2]);
        float radius = sqrt(pow(lights[1]->position.x, 2) + pow(lights[1]->position.y, 2) + pow(lights[1]->position.z, 2));
        // Assuming x is constant, y and z coordinates are equal to R(radius) * cos(angle) and R * cos(angle) if we consider the coordinates are moving in a circle.
        float y = radius * cos(angle);
        float z = radius * sin(angle);
        lights[1]->position = glm::vec3(lights[1]->position.x, y, z);
    }

    Image image(width, height); // Create an image where we will store the result

    float s = 2 * tan(0.5 * fov / 180 * M_PI) / width;
    float X = -s * width / 2;
    float Y = s * height / 2;

    for (int i = 0; i < width; i++)
        for (int j = 0; j < height; j++)
        {

            float dx = X + i * s + s / 2;
            float dy = Y - j * s - s / 2;
            float dz = 1;

            glm::vec3 origin(0, 0, 0);
            glm::vec3 direction(dx, dy, dz);
            direction = glm::normalize(direction);

            Ray ray(origin, direction);

            image.setPixel(i, j, trace_ray(ray));
        }

    t = clock() - t;
    cout << "It took " << ((float)t) / CLOCKS_PER_SEC << " seconds to render the image." << endl;
    cout << "I could render at " << (float)CLOCKS_PER_SEC / ((float)t) << " frames per second." << endl;

    // Writing the final results of the rendering
    if (argc > 2)
    {
        image.writeImage(argv[1]);
    }
    else
    {
        image.writeImage("./result.ppm");
    }

    return 0;
}
