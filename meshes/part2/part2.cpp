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
#include "OBJloader.h"

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

class Triangle : public Object
{
private:
    Face face;
    glm::vec3 p_max;
    glm::vec3 p_min;

public:
    /**
     The constructor of the sphere
     @param face Face of the triangle
     @param color Color of the sphere
     */
    Triangle(Face face, glm::vec3 color) : face(face)
    {
        this->color = color;
    }
    Triangle(Face face, Material material) : face(face)
    {
        this->material = material;
    }
    /** Implementation of the intersection function*/
    Hit intersect(Ray ray)
    {
        float n_normal_d = glm::dot(face.triangleNormal, ray.direction);
        Hit hit;
        hit.hit = false;
        if (n_normal_d < 0.0001f)
        {
            return hit;
        }
        glm::vec3 distanceVector = face.p1 - ray.origin;
        float D = glm::dot(face.triangleNormal, face.p1) * -1;
        float t = glm::dot(face.triangleNormal, distanceVector) / n_normal_d;
        glm::vec3 onPlanePoint = ray.origin + ray.direction * t;

        glm::vec3 edge0 = face.p2 - face.p1;
        glm::vec3 edge1 = face.p3 - face.p2;
        glm::vec3 edge2 = face.p1 - face.p3;
        glm::vec3 c0 = onPlanePoint - face.p1;
        glm::vec3 c1 = onPlanePoint - face.p2;
        glm::vec3 c2 = onPlanePoint - face.p3;

        if (glm::dot(face.triangleNormal, glm::cross(edge0, c0)) >= 0 && glm::dot(face.triangleNormal, glm::cross(edge1, c1)) >= 0 && glm::dot(face.triangleNormal, glm::cross(edge2, c2)) >= 0)
        {
            hit.hit = true;
            hit.intersection = onPlanePoint;

            // float alpha = glm::determinant(glm::mat3x3(onPlanePoint, face.p2, face.p3)) / face.det;
            // float beta = glm::determinant(glm::mat3x3(face.p1, onPlanePoint, face.p3)) / face.det;
            // float gamma = glm::determinant(glm::mat3x3(face.p1, face.p2, onPlanePoint)) / face.det;

            // hit.normal = alpha * face.n1 + beta * face.n2 + gamma * face.n3;
            hit.normal = glm::normalize(face.triangleNormal);
            hit.distance = glm::distance(ray.origin, hit.intersection);
            hit.object = this;
        }
        return hit;
    }

    Face getFace()
    {
        return face;
    }
};

class Plane : public Object
{

private:
    glm::vec3 normal;
    glm::vec3 point;

public:
    Plane(glm::vec3 point, glm::vec3 normal) : point(point), normal(normal)
    {
    }
    Plane(glm::vec3 point, glm::vec3 normal, Material material) : point(point), normal(normal)
    {
        this->material = material;
    }
    Hit intersect(Ray ray)
    {

        Hit hit;
        hit.hit = false;

        float a = glm::dot(point - ray.origin, normal);
        float b = glm::dot(ray.direction, normal);

        float t = a / b;

        if (b != 0 && t > 0)
        {
            hit.hit = true;
            hit.intersection = ray.origin + ray.direction * t;
            hit.normal = glm::normalize(normal);
            hit.distance = glm::distance(ray.origin, hit.intersection);
            hit.object = this;
        }

        return hit;
    }
};

class Box : public Object
{

private:
    glm::vec3 pmin;
    glm::vec3 pmax;

public:
    Box(glm::vec3 pmin, glm::vec3 pmax) : pmin(pmin), pmax(pmax)
    {
    }
    Hit intersect(Ray ray)
    {
        Hit hit;
        hit.hit = false;

        // lb is the corner of AABB with minimal coordinates - left bottom, rt is maximal corner
        // r.org is origin of ray
        float t1 = (pmin.x - ray.origin.x) / ray.direction.x;
        float t2 = (pmax.x - ray.origin.x) / ray.direction.x;
        float t3 = (pmin.y - ray.origin.y) / ray.direction.y;
        float t4 = (pmax.y - ray.origin.y) / ray.direction.y;
        float t5 = (pmin.z - ray.origin.z) / ray.direction.z;
        float t6 = (pmax.z - ray.origin.z) / ray.direction.z;

        float tmin = max(max(min(t1, t2), min(t3, t4)), min(t5, t6));
        float tmax = min(min(max(t1, t2), max(t3, t4)), max(t5, t6));

        // Check for intersection
        if (tmin <= tmax)
        {
            hit.hit = true;
            hit.object = this;
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

struct bvhStruct
{
    glm::vec3 pmin;
    glm::vec3 pmax;
    vector<int> objectIndices;
    Box *box;
    struct bvhStruct *rightTree;
    struct bvhStruct *leftTree;
};

vector<Light *> lights; ///< A list of lights in the scene
glm::vec3 ambient_light(1.0, 1.0, 1.0);
vector<Object *> objects; ///< A list of all objects in the scene
vector<Face> tris;
bvhStruct *tree;

/** Function for computing color of an object according to the Phong Model
 @param point A point belonging to the object for which the color is computed
 @param normal A normal vector the the point
 @param view_direction A normalized direction from the point to the viewer/camera
 @param material A material structure representing the material of the object
*/
glm::vec3 PhongModel(glm::vec3 point, glm::vec3 normal, glm::vec3 view_direction, Material material)
{

    glm::vec3 color(0.0);

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

vector<int> traverseTree(bvhStruct *branch, Ray ray)
{
    if (branch->leftTree == nullptr && branch->rightTree == nullptr)
    {
        return branch->objectIndices;
    }
    bool leftHit = branch->leftTree != nullptr ? branch->leftTree->box->intersect(ray).hit : false;
    bool rightHit = branch->rightTree != nullptr ? branch->rightTree->box->intersect(ray).hit : false;
    vector<int> leftResult, rightResult;
    if (leftHit)
    {
        leftResult = traverseTree(branch->leftTree, ray);
    }
    if (rightHit)
    {
        rightResult = traverseTree(branch->rightTree, ray);
    }
    if (!leftHit && !rightHit)
    {
        return vector<int>();
    }
    leftResult.insert(leftResult.end(), rightResult.begin(), rightResult.end());
    return leftResult;
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
    bvhStruct *traverser = new bvhStruct();
    traverser->pmin = tree->pmin;
    traverser->pmax = tree->pmax;
    traverser->box = tree->box;
    traverser->leftTree = tree->leftTree;
    traverser->rightTree = tree->rightTree;
    traverser->objectIndices = {};
    for (int i : tree->objectIndices)
    {
        traverser->objectIndices.push_back(i);
    }
    vector<int> target = traverseTree(traverser, ray);
    for (int k : target)
    {
        Hit hit = objects[k]->intersect(ray);
        if (hit.hit == true && hit.distance < closest_hit.distance)
            closest_hit = hit;
    }

    for (int k = 0; k < 6; k++)
    {
        Hit hit = objects[k]->intersect(ray);
        if (hit.hit == true && hit.distance < closest_hit.distance)
            closest_hit = hit;
    }

    glm::vec3 color(0.0);

    if (closest_hit.hit)
    {
        color = PhongModel(closest_hit.intersection, closest_hit.normal, glm::normalize(-ray.direction), closest_hit.object->getMaterial());
        // color = PhongModel(closest_hit.intersection, closest_hit.normal, glm::normalize(-ray.direction), closest_hit.object->getMaterial());
    }
    else
    {
        color = glm::vec3(0.0, 0.0, 0.0);
    }
    return color;
}

glm::mat3x3 getTranslationMatrix(float xRad, float yRad, float zRad)
{
    float cX = cos(xRad);
    float cY = cos(yRad);
    float cZ = cos(zRad);

    float sX = sin(xRad);
    float sY = sin(yRad);
    float sZ = sin(zRad);

    return glm::mat3x3(
        cY * cZ, sX * sY * cZ - cX * sZ, cX * sY * cZ + sX * sZ,
        cY * sZ, sX * sY * sZ + cX * cZ, cX * sY * sZ - sX * cZ,
        -sY, sX * cY, cX * cY);
}

vector<glm::vec3> getBoundingBox(vector<int> points)
{
    glm::vec3 pmin = glm::vec3(INT_MAX);
    glm::vec3 pmax = glm::vec3(INT_MIN);
    for (int p : points)
    {
        pmin.x = min({pmin.x, tris[p].p1.x, tris[p].p2.x, tris[p].p3.x}, less<float>());
        pmin.y = min({pmin.y, tris[p].p1.y, tris[p].p2.y, tris[p].p3.y}, less<float>());
        pmin.z = min({pmin.z, tris[p].p1.z, tris[p].p2.z, tris[p].p3.z}, less<float>());
        pmax.x = max({pmax.x, tris[p].p1.x, tris[p].p2.x, tris[p].p3.x});
        pmax.y = max({pmax.y, tris[p].p1.y, tris[p].p2.y, tris[p].p3.y});
        pmax.z = max({pmax.z, tris[p].p1.z, tris[p].p2.z, tris[p].p3.z});
    }
    return {pmin, pmax};
}

bool inRange(glm::vec3 min, glm::vec3 max, glm::vec3 point)
{
    if (min.x > point.x || max.x < point.x)
    {
        return false;
    }
    if (min.y > point.y || max.y < point.y)
    {
        return false;
    }
    if (min.z > point.z || max.z < point.z)
    {
        return false;
    }
    return true;
}

bvhStruct *bvh_node(vector<int> points, int direction)
{
    vector<glm::vec3> coords = getBoundingBox(points);
    bvhStruct *current = new bvhStruct();
    current->pmin = coords[0];
    current->pmax = coords[1];
    current->objectIndices = {};
    current->box = new Box(current->pmin, current->pmax);
    if (points.size() < 500)
    {
        for (int l : points)
        {
            current->objectIndices.push_back(l);
        }
        return current;
    }
    glm::vec3 newMax, newMin;
    if (direction % 3 == 0)
    { // cut y-val
        float yVal = (coords[1].y - coords[0].y) / 2 + coords[0].y;
        newMax = glm::vec3(coords[1].x, yVal, coords[1].z);
        newMin = glm::vec3(coords[0].x, yVal, coords[0].z);
    }
    else if (direction % 3 == 1)
    { // cut x-val
        float xVal = (coords[1].x - coords[0].x) / 2 + coords[0].x;
        newMax = glm::vec3(xVal, coords[1].y, coords[1].z);
        newMin = glm::vec3(xVal, coords[0].y, coords[0].z);
    }
    else
    { // cut z-val
        float zVal = (coords[1].z - coords[0].z) / 2 + coords[0].z;
        newMax = glm::vec3(coords[1].x, coords[1].y, zVal);
        newMin = glm::vec3(coords[0].x, coords[0].y, zVal);
    }
    vector<int> leftPoints;
    vector<int> rightPoints;
    for (int i : points)
    {
        if (inRange(current->pmin, newMax, tris[i].p1) || inRange(current->pmin, newMax, tris[i].p2) || inRange(current->pmin, newMax, tris[i].p3))
        {
            leftPoints.push_back(i);
        }
        if (inRange(newMin, current->pmax, tris[i].p1) || inRange(newMin, current->pmax, tris[i].p2) || inRange(newMin, current->pmax, tris[i].p3))
        {
            rightPoints.push_back(i);
        }
    }
    current->leftTree = leftPoints.empty() ? nullptr : bvh_node(leftPoints, direction + 1);
    current->rightTree = rightPoints.empty() ? nullptr : bvh_node(rightPoints, direction + 1);
    return current;
}

/**
 Function defining the scene
 */
void sceneDefinition()
{
    glm::vec3 bunnyStartingPos = glm::vec3(0.0f, -3.0f, 9.0f);
    glm::vec3 armaStartingPos = glm::vec3(-3.0f, -2.0f, 9.0f);
    glm::vec3 lucyStartingPos = glm::vec3(3.0f, -2.0f, 9.0f);
    glm::mat3x3 armaRotate = getTranslationMatrix(glm::radians(-15.0f), glm::radians(150.0f), 0.0f);
    glm::mat3x3 noRotate = getTranslationMatrix(0.0f, 0.0f, 0.0f);

    // passing the filepath and 3d object position and rotation
    Triangles bunnyTriangles = loadOBJ("../meshes/bunny_small.obj", bunnyStartingPos, noRotate);
    tris.insert(tris.end(), bunnyTriangles.faces.begin(), bunnyTriangles.faces.end());
    glm::vec3 bunny_min = bunnyTriangles.p_min;
    glm::vec3 bunny_max = bunnyTriangles.p_max;
    Triangles armaTriangles = loadOBJ("../meshes/armadillo_small.obj", armaStartingPos, armaRotate);
    tris.insert(tris.end(), armaTriangles.faces.begin(), armaTriangles.faces.end());
    glm::vec3 arma_min = armaTriangles.p_min;
    glm::vec3 arma_max = armaTriangles.p_max;
    Triangles lucyTriangles = loadOBJ("../meshes/lucy_small.obj", lucyStartingPos, noRotate);
    tris.insert(tris.end(), lucyTriangles.faces.begin(), lucyTriangles.faces.end());
    glm::vec3 lucy_min = lucyTriangles.p_min;
    glm::vec3 lucy_max = lucyTriangles.p_max;

    Material red_specular;
    red_specular.diffuse = glm::vec3(1.0f, 0.3f, 0.3f);
    red_specular.ambient = glm::vec3(0.01f, 0.03f, 0.03f);
    red_specular.specular = glm::vec3(0.5);
    red_specular.shininess = 10.0;

    Material green_specular;
    green_specular.diffuse = glm::vec3(0.7f, 0.9f, 0.7f);
    green_specular.ambient = glm::vec3(0.07f, 0.09f, 0.07f);
    green_specular.specular = glm::vec3(0.0);
    green_specular.shininess = 0.0;

    Material white_wall;
    white_wall.ambient = glm::vec3(0.2f);
    white_wall.diffuse = glm::vec3(1.0f);
    white_wall.specular = glm::vec3(1.0);
    white_wall.shininess = 100.0;

    Material purple_wall;
    purple_wall.ambient = glm::vec3(0.1f, 0.01f, 0.0f);
    purple_wall.diffuse = glm::vec3(0.7f, 0.7f, 1.0f);

    Material pink_wall;
    pink_wall.diffuse = glm::vec3(1.3f, 0.5f, 0.8f);
    pink_wall.ambient = glm::vec3(0.03f, 0.01f, 0.0f);

    // extending x-axis
    objects.push_back(new Plane(glm::vec3(-15, 0, 0), glm::vec3(1, 0, 0), pink_wall));
    objects.push_back(new Plane(glm::vec3(15, 0, 0), glm::vec3(-1, 0, 0), purple_wall));

    // extending y-axis
    objects.push_back(new Plane(glm::vec3(0, -3, 0), glm::vec3(0, 1, 0), white_wall));
    objects.push_back(new Plane(glm::vec3(0, 27, 0), glm::vec3(0, -1, 0), white_wall));

    // extending z-axis
    objects.push_back(new Plane(glm::vec3(0, 0, -0.01), glm::vec3(0, 0, 1), green_specular));
    objects.push_back(new Plane(glm::vec3(0, 0, 30), glm::vec3(0, 0, -1), green_specular));
    vector<int> idk;
    for (int i = 0; i < tris.size(); i++)
    {
        objects.push_back(new Triangle(tris[i], red_specular));
        idk.push_back(i);
    }
    tree = bvh_node(idk, 0);
    lights.push_back(new Light(glm::vec3(0, 26, 5), glm::vec3(0.3)));
    lights.push_back(new Light(glm::vec3(0, 1, 12), glm::vec3(0.3)));
    lights.push_back(new Light(glm::vec3(0, 5, 1), glm::vec3(0.3)));
}

int main(int argc, const char *argv[])
{

    clock_t t = clock(); // variable for keeping the time of the rendering

    int width = 1024; // width of the image
    int height = 768; // height of the image
    float fov = 90;   // field of view

    sceneDefinition(); // Let's define a scene

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

            glm::vec3 origin(0, 0, 0); // z-axis:front and back y-axis: left and right x-axis: up and down
            glm::vec3 direction(dx, dy, dz);
            direction = glm::normalize(direction);

            Ray ray(origin, direction);

            image.setPixel(i, j, trace_ray(ray));
        }

    t = clock() - t;
    cout << "It took " << ((float)t) / CLOCKS_PER_SEC << " seconds to render the image." << endl;
    cout << "I could render at " << (float)CLOCKS_PER_SEC / ((float)t) << " frames per second." << endl;

    // Writing the final results of the rendering
    if (argc == 2)
    {
        image.writeImage(argv[1]);
    }
    else
    {
        image.writeImage("./result.ppm");
    }

    return 0;
}
