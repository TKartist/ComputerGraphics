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
	glm::vec3 origin;	 ///< Origin of the ray
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
	bool hit;				///< Boolean indicating whether there was or there was no intersection with an object
	glm::vec3 normal;		///< Normal vector of the intersected object at the intersection point
	glm::vec3 intersection; ///< Point of Intersection
	float distance;			///< Distance from the origin of the ray to the intersection point
	Object *object;			///< A pointer to the intersected object
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
	float radius;	  ///< Radius of the sphere
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

class Triangle : public Object
{
private:
	Face face;

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

	Hit intersect(Ray ray)
	{
		Hit hit;
		hit.hit = false;

		float n_normal_d = glm::dot(face.triangleNormal, ray.direction);

		if (n_normal_d < 0.001f)
		{
			return hit;
		}

		float t = glm::dot(face.triangleNormal, face.p1 - ray.origin) / n_normal_d;
		glm::vec3 onPlanePoint = ray.origin + ray.direction * t;

		// Barycentric coordinates calculation
		glm::vec3 edge0 = face.p2 - face.p1;
		glm::vec3 edge1 = face.p3 - face.p1;
		glm::vec3 c0 = onPlanePoint - face.p1;

		float denom = glm::dot(edge0, edge0) * glm::dot(edge1, edge1) - glm::dot(edge0, edge1) * glm::dot(edge0, edge1);
		float alpha = (glm::dot(edge1, edge1) * glm::dot(c0, edge0) - glm::dot(edge0, edge1) * glm::dot(c0, edge1)) / denom;
		float beta = (glm::dot(edge0, edge0) * glm::dot(c0, edge1) - glm::dot(edge0, edge1) * glm::dot(c0, edge0)) / denom;
		float gamma = 1.0f - alpha - beta;

		if (alpha >= 0 && beta >= 0 && gamma >= 0)
		{
			hit.hit = true;
			hit.intersection = onPlanePoint;
			hit.normal = glm::normalize(alpha * face.n1 + beta * face.n2 + gamma * face.n3);
			hit.distance = glm::distance(ray.origin, hit.intersection);
			hit.object = this;
		}

		return hit;
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

		/*
		 Excercise 1 - Plane-ray intersection
		 */

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

/**
 Light class
 */
class Light
{
public:
	glm::vec3 position; ///< Position of the light source
	glm::vec3 color;	///< Color/intentisty of the light source
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
/**
 Function defining the scene
 */
void sceneDefinition()
{
	glm::vec3 bunnyStartingPos = glm::vec3(0.0f, -3.0f, 9.0f);
	glm::vec3 armaStartingPos = glm::vec3(-5.0f, -2.0f, 9.0f);
	glm::vec3 lucyStartingPos = glm::vec3(6.0f, -2.0f, 9.0f);
	glm::mat3x3 armaRotate = getTranslationMatrix(glm::radians(-15.0f), glm::radians(150.0f), 0.0f);
	glm::mat3x3 noRotate = getTranslationMatrix(0.0f, 0.0f, 0.0f);
	// glm::vec3 lucyStartingPos = glm::vec3(0.0f, -3.0f, 9.0f);

	// passing the filepath and 3d object position
	vector<Face> bunny = loadOBJ("./meshes/bunny.obj", bunnyStartingPos, noRotate);
	vector<Face> arma = loadOBJ("./meshes/armadillo.obj", armaStartingPos, armaRotate);
	vector<Face> lucy = loadOBJ("./meshes/lucy.obj", lucyStartingPos, noRotate);

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

	Material white_specular;
	white_specular.diffuse = glm::vec3(1.0f, 0.5f, 0.31f);
	white_specular.ambient = glm::vec3(1.0f, 0.5f, 0.31f);
	white_specular.specular = glm::vec3(0.5f, 0.5f, 0.5f);
	white_specular.shininess = 32.0f;

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

	Material turquoise_wall;
	turquoise_wall.ambient = glm::vec3(0.01f, 0.0f, 0.1f);
	turquoise_wall.diffuse = glm::vec3(0.4f, 0.7f, 0.7f);

	// objects.push_back(new Sphere(0.5, glm::vec3(-1.0, -2.5, 6.0), red_specular));
	// objects.push_back(new Sphere(1.0, glm::vec3(1.0, -2.0, 8.0), blue_specular));
	// objects.push_back(new Sphere(1.0, glm::vec3(3.0, -2.0, 6.0), green_specular));

	// for (int i = 0; i < bunny.size(); i++)
	// {
	// 	objects.push_back(new Triangle(bunny[i], red_specular));
	// }
	for (int i = 0; i < arma.size(); i++)
	{
		objects.push_back(new Triangle(arma[i], blue_specular));
	}
	// for (int i = 0; i < lucy.size(); i++)
	// {
	// 	objects.push_back(new Triangle(lucy[i], green_specular));
	// }

	// extending x-axis
	objects.push_back(new Plane(glm::vec3(-15, 0, 0), glm::vec3(1, 0, 0), pink_wall));
	objects.push_back(new Plane(glm::vec3(15, 0, 0), glm::vec3(-1, 0, 0), purple_wall));

	// extending y-axis
	objects.push_back(new Plane(glm::vec3(0, -3, 0), glm::vec3(0, 1, 0), white_wall));
	objects.push_back(new Plane(glm::vec3(0, 27, 0), glm::vec3(0, -1, 0), white_wall));

	// extending z-axis
	objects.push_back(new Plane(glm::vec3(0, 0, -0.01), glm::vec3(0, 0, 1), green_specular));
	objects.push_back(new Plane(glm::vec3(0, 0, 30), glm::vec3(0, 0, -1), green_specular));

	lights.push_back(new Light(glm::vec3(0, 26, 5), glm::vec3(0.5)));
	lights.push_back(new Light(glm::vec3(0, 1, 12), glm::vec3(0.5)));
	lights.push_back(new Light(glm::vec3(0, 5, 1), glm::vec3(0.5)));
}

int main(int argc, const char *argv[])
{

	clock_t t = clock(); // variable for keeping the time of the rendering

	int width = 1024; // width of the image
	int height = 768; // height of the image
	float fov = 90;	  // field of view

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
