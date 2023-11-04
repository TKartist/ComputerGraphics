/**
@file main.cpp
*/

#include <iostream>
#include <fstream>
#include <cmath>
#include <ctime>
#include <vector>
#include "glm/glm.hpp"
#include "glm/gtx/transform.hpp"

#include "Image.h"
#include "Material.h"

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
	glm::vec2 uv;			///< Coordinates for computing the texture (texture coordinates)
};

/**
 General class for the object
 */
class Object
{

protected:
	glm::mat4 transformationMatrix;		   ///< Matrix representing the transformation from the local to the global coordinate system
	glm::mat4 inverseTransformationMatrix; ///< Matrix representing the transformation from the global to the local coordinate system
	glm::mat4 normalMatrix;				   ///< Matrix for transforming normal vectors from the local to the global coordinate system

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
	 @param material A structure describing the material of the object
	*/
	void setMaterial(Material material)
	{
		this->material = material;
	}
	/** Functions for setting up all the transformation matrices
	 @param matrix The matrix representing the transformation of the object in the global coordinates */
	void setTransformation(glm::mat4 matrix)
	{

		transformationMatrix = matrix;
		inverseTransformationMatrix = glm::inverse(transformationMatrix);
		normalMatrix = glm::transpose(inverseTransformationMatrix);

		/* ----- Assignment 5 ---------
		 Set the two remaining matrices

		inverseTransformationMatrix =
		normalMatrix =

		 */
	}

	Ray localRay(Ray ray)
	{
		glm::vec4 transformedOrigin = inverseTransformationMatrix * glm::vec4(ray.origin, 1.0f);
		glm::vec4 transformedDirection = inverseTransformationMatrix * glm::vec4(ray.direction, 0.0f);

		ray.origin = glm::vec3(transformedOrigin);
		ray.direction = glm::normalize(glm::vec3(transformedDirection));
		return ray;
	}

	Hit globalHit(Hit hit, glm::vec3 origin)
	{
		if (!hit.hit)
		{
			return hit;
		}
		Hit global = hit;
		glm::vec3 intersects = hit.intersection;
		global.intersection = glm::vec3(transformationMatrix * glm::vec4(intersects, 1.0));
		global.normal = glm::normalize(glm::vec3(transpose(inverseTransformationMatrix) * glm::vec4(intersects, 0.0f)));
		global.distance = glm::distance(origin, global.intersection);
		return global;
	}
};

/**
 Implementation of the class Object for sphere shape.
 */
class Sphere : public Object
{
private:
	float radius = 1.0f;				///< Radius of the sphere
	glm::vec3 center = glm::vec3(0.0f); ///< Center of the sphere

public:
	/**
	 The constructor of the sphere
	 @param radius Radius of the sphere
	 @param center Center of the sphere
	 @param color Color of the sphere
	 */
	Sphere(glm::vec3 color)
	{
		this->color = color;
	}
	Sphere(Material material)
	{
		this->material = material;
	}
	/** Implementation of the intersection function*/
	Hit intersect(Ray ray)
	{

		Hit hit;
		hit.hit = false;
		glm::vec3 o = ray.origin;
		Ray local = localRay(ray);

		glm::vec3 c = center - local.origin;

		float a = glm::dot(c, local.direction);

		float D_squared = glm::dot(c, c) - a * a;
		float D = sqrt(D_squared);

		float closest_t;
		if (D > radius)
		{
			return hit;
		}
		else if (D < radius)
		{
			float t1 = a - sqrt(radius * radius - D_squared);
			closest_t = t1 >= 0 ? t1 : a + sqrt(radius * radius - D_squared);
		}
		else
		{
			closest_t = a + sqrt(radius * radius - D_squared);
		}
		if (closest_t < 0)
		{
			return hit;
		}

		hit.hit = true;
		hit.distance = closest_t;
		hit.intersection = local.origin + local.direction * hit.distance;
		hit.normal = glm::normalize(hit.intersection - center);
		hit.object = this;

		hit.uv.x = 0.5 - asin(hit.normal.y) / M_PI;
		hit.uv.y = 2 * (0.5 + atan2(hit.normal.z, hit.normal.x) / (2 * M_PI));

		return globalHit(hit, o);
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
		float DdotN = glm::dot(ray.direction, normal);
		if (DdotN < 0)
		{

			float PdotN = glm::dot(point - ray.origin, normal);
			float t = PdotN / DdotN;

			if (t > 0)
			{
				hit.hit = true;
				hit.normal = normal;
				hit.distance = t;
				hit.object = this;
				hit.intersection = t * ray.direction + ray.origin;
			}
		}
		return hit;
	}
};

class Cone : public Object
{
public:
	Cone(Material material)
	{
		this->material = material;
	}
	Hit intersect(Ray ray)
	{

		Hit hit;
		hit.hit = false;

		/*  ---- Assignemnt 5 -----

		 Implement the ray-cone intersection. Before intersecting the ray with the cone,
		 make sure that you transform the ray into the local coordinate system.
		 Remember about normalizing all the directions after transformations.
		*/
		glm::vec3 o = ray.origin;
		Ray local = localRay(ray);

		float a = pow(local.direction.x, 2) + pow(local.direction.z, 2) - pow(local.direction.y, 2);
		float b = 2.0f * (local.direction.x * local.origin.x + local.direction.z * local.origin.z - local.direction.y * local.origin.y);
		float c = pow(local.origin.x, 2) + pow(local.origin.z, 2) - pow(local.origin.y, 2);

		float disc = pow(b, 2) - 4 * a * c;

		if (disc < 0)
		{
			return hit;
		}

		float t1 = (-b - sqrt(disc)) / (2.0f * a);
		float t2 = (-b + sqrt(disc)) / (2.0f * a);

		float t = min(t1, t2);

		if (t > 0)
		{
			hit.intersection = local.origin + t * local.direction;
			float y = hit.intersection.y;

			Plane *base = new Plane(glm::vec3(0, 1, 0), glm::vec3(0, 1, 0));
			Hit baseHit = base->intersect(local);
			if (baseHit.hit)
			{
				auto i = baseHit.intersection - glm::vec3(0, 1, 0);
				float n = glm::dot(i, i);
				if (sqrt(n) <= 1)
				{
					hit.hit = true;
					hit.object = this;
					hit.intersection = glm::vec3((transformationMatrix * glm::vec4(baseHit.intersection, 1)));
					hit.normal = glm::normalize(glm::vec3(normalMatrix * glm::vec4(baseHit.normal, 0)));
					hit.distance = glm::distance(hit.intersection, o);
					return hit;
				}
			}
			if (y >= 0 && y <= 1)
			{
				hit.hit = true;
				hit.object = this;
				return globalHit(hit, o);
			}
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
glm::vec3 ambient_light(0.001, 0.001, 0.001);
vector<Object *> objects; ///< A list of all objects in the scene

/** Function for computing color of an object according to the Phong Model
 @param point A point belonging to the object for which the color is computer
 @param normal A normal vector the the point
 @param uv Texture coordinates
 @param view_direction A normalized direction from the point to the viewer/camera
 @param material A material structure representing the material of the object
*/
glm::vec3 PhongModel(glm::vec3 point, glm::vec3 normal, glm::vec2 uv, glm::vec3 view_direction, Material material)
{

	glm::vec3 color(0.0);
	for (int light_num = 0; light_num < lights.size(); light_num++)
	{

		glm::vec3 light_direction = glm::normalize(lights[light_num]->position - point);
		glm::vec3 reflected_direction = glm::reflect(-light_direction, normal);

		float NdotL = glm::clamp(glm::dot(normal, light_direction), 0.0f, 1.0f);
		float VdotR = glm::clamp(glm::dot(view_direction, reflected_direction), 0.0f, 1.0f);

		glm::vec3 diffuse_color = material.diffuse;
		if (material.texture)
		{
			diffuse_color = material.texture(uv);
		}

		glm::vec3 diffuse = diffuse_color * glm::vec3(NdotL);
		glm::vec3 specular = material.specular * glm::vec3(pow(VdotR, material.shininess));

		// distance to the light
		float r = glm::distance(point, lights[light_num]->position);
		r = max(r, 0.1f);

		color += lights[light_num]->color * (diffuse + specular) / r / r;
	}
	color += ambient_light * material.ambient;

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
		color = PhongModel(closest_hit.intersection, closest_hit.normal, closest_hit.uv, glm::normalize(-ray.direction), closest_hit.object->getMaterial());
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
void sceneDefinition(int frameIndex)
{

	Material green_diffuse;
	green_diffuse.ambient = glm::vec3(0.03f, 0.1f, 0.03f);
	green_diffuse.diffuse = glm::vec3(0.3f, 1.0f, 0.3f);

	Material red_specular;
	red_specular.diffuse = glm::vec3(1.0f, 0.2f, 0.2f);
	red_specular.ambient = glm::vec3(0.01f, 0.02f, 0.02f);
	red_specular.specular = glm::vec3(0.5);
	red_specular.shininess = 10.0;

	Material blue_specular;
	blue_specular.ambient = glm::vec3(0.02f, 0.02f, 0.1f);
	blue_specular.diffuse = glm::vec3(0.2f, 0.2f, 1.0f);
	blue_specular.specular = glm::vec3(0.6);
	blue_specular.shininess = 100.0;

	glm::mat4 s1m = glm::translate(glm::mat4(1.0), glm::vec3(1, -2, 8)) * glm::scale(glm::mat4(1.0), glm::vec3(1)) * glm::rotate(glm::radians(0.0f), glm::vec3(1.0f, 1.0f, 1.0f));
	Sphere *s1 = new Sphere(blue_specular);
	s1->setTransformation(s1m);
	objects.push_back(s1);
	glm::mat4 s2m = glm::translate(glm::mat4(1.0), glm::vec3(-1, -2.5, 6)) * glm::scale(glm::mat4(1.0), glm::vec3(0.5)) * glm::rotate(glm::radians(0.0f), glm::vec3(1.0f, 1.0f, 1.0f));
	Sphere *s2 = new Sphere(red_specular);
	s2->setTransformation(s2m);
	objects.push_back(s2);

	// ------ Assignment 5 -------

	// You can remove the green sphere as it should be replaced with a green cone
	// objects.push_back(new Sphere(1.0, glm::vec3(3, -2, 6), green_diffuse));

	// Textured sphere
	Material textured;
	textured.texture = &rainbowTexture;
	float calc = frameIndex * 12 / 120;
	if (calc > 6)
	{
		calc = 6 - remainder(calc, 6);
	}
	glm::mat4 s3m = glm::translate(glm::mat4(1.0), glm::vec3(calc, 4, 23)) * glm::scale(glm::mat4(1.0), glm::vec3(7.0)) * glm::rotate(glm::radians(0.0f), glm::vec3(1.0f, 1.0f, 1.0f));
	Sphere *s3 = new Sphere(textured);
	s3->setTransformation(s3m);
	objects.push_back(s3);

	// Planes
	Material red_diffuse;
	red_diffuse.ambient = glm::vec3(0.09f, 0.06f, 0.06f);
	red_diffuse.diffuse = glm::vec3(0.9f, 0.6f, 0.6f);

	Material blue_diffuse;
	blue_diffuse.ambient = glm::vec3(0.06f, 0.06f, 0.09f);
	blue_diffuse.diffuse = glm::vec3(0.6f, 0.6f, 0.9f);
	objects.push_back(new Plane(glm::vec3(0, -3, 0), glm::vec3(0.0, 1, 0)));
	objects.push_back(new Plane(glm::vec3(0, 1, 30), glm::vec3(0.0, 0.0, -1.0), green_diffuse));
	objects.push_back(new Plane(glm::vec3(-15, 1, 0), glm::vec3(1.0, 0.0, 0.0), red_diffuse));
	objects.push_back(new Plane(glm::vec3(15, 1, 0), glm::vec3(-1.0, 0.0, 0.0), blue_diffuse));
	objects.push_back(new Plane(glm::vec3(0, 27, 0), glm::vec3(0.0, -1, 0)));
	objects.push_back(new Plane(glm::vec3(0, 1, -0.01), glm::vec3(0.0, 0.0, 1.0), green_diffuse));

	/* ----- Assignment 5 -------
	Create two conse and add them to the collection of our objects.
	Remember to create them with corresponding materials and transformation matrices


	Cone *cone1 = new Cone(...);
	cone1->setTransformation(...);
	objects.push_back(cone1);

	Cone *cone2 = new Cone(...);
	cone2->setTransformation(...);
	objects.push_back(cone2);

	*/

	Material yellow_specular;
	yellow_specular.ambient = glm::vec3(1.0f);
	yellow_specular.diffuse = glm::vec3(1.0f, 1.0f, 0.0f);
	yellow_specular.specular = glm::vec3(0.6);
	yellow_specular.shininess = 100.0;

	glm::mat4 transformationMatrix1 = glm::translate(glm::mat4(1.0), glm::vec3(5, 9, 14)) * glm::scale(glm::mat4(1.0), glm::vec3(3, 12, 3)) * glm::rotate(glm::radians(180.0f), glm::vec3(1.0f, 0.0f, 1.0f));
	Cone *cone1 = new Cone(yellow_specular);
	cone1->setTransformation(transformationMatrix1);
	objects.push_back(cone1);

	glm::mat4 cone2Transformation = glm::translate(glm::mat4(1.0), glm::vec3(6, -3, 7)) * glm::rotate(glm::radians(70.0f), glm::vec3(0, 0, 1)) * glm::scale(glm::mat4(1.0), glm::vec3(1, 3, 1));

	Cone *cone2 = new Cone(green_diffuse);
	cone2->setTransformation(cone2Transformation);
	objects.push_back(cone2);

	lights.push_back(new Light(glm::vec3(0, 26, 5), glm::vec3(1.0, 1.0, 1.0)));
	lights.push_back(new Light(glm::vec3(0, 1, 12), glm::vec3(0.1)));
	lights.push_back(new Light(glm::vec3(0, 5, 1), glm::vec3(0.4)));
}

/**
 Function performing tonemapping of the intensities computed using the raytracer
 @param intensity Input intensity
 @return Tonemapped intensity in range (0,1)
 */
glm::vec3 toneMapping(glm::vec3 intensity)
{
	float gamma = 1.0 / 2.0;
	float alpha = 12.0f;
	return glm::clamp(alpha * glm::pow(intensity, glm::vec3(gamma)), glm::vec3(0.0), glm::vec3(1.0));
}

int main(int argc, const char *argv[])
{

	clock_t t = clock(); // variable for keeping the time of the rendering

	int width = 1024; // width of the image
	int height = 768; // height of the image
	float fov = 90;	  // field of view

	int frameIndex = 0;
	if (argc > 2)
	{
		frameIndex = atoi(argv[2]);
	}
	sceneDefinition(frameIndex); // Let's define a scene

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

			image.setPixel(i, j, toneMapping(trace_ray(ray)));
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
