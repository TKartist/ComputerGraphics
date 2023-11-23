
#ifndef OBJloader_h
#define OBJloader_h

#include <iostream>
#include <string>
#include <fstream>
#include <vector>
#include <sstream>
#include <algorithm>
#include "glm/glm.hpp"
#include <cmath>
#include <ctime>
#include <GLFW/glfw3.h>

using namespace std;

struct Face
{
    glm::vec3 p1;
    glm::vec3 p2;
    glm::vec3 p3;
    glm::vec3 n1; // vector normal of p1
    glm::vec3 n2; // vecotr normal of p2
    glm::vec3 n3; // vector normal of p3
    glm::vec3 triangleNormal;
    float det;
};
struct Triangles
{
    vector<Face> faces;
};

static Triangles loadOBJ(const char *file_name, const glm::vec3 startingPos)
{
    vector<glm::vec3> vectors;
    vector<glm::vec3> vectorNormals;
    stringstream ss;
    ifstream in_file(file_name);

    // file not open error handler
    if (!in_file.is_open())
    {
        throw "Error: OBJLOADER: could not open file";
    }

    string line = "";
    glm::vec3 val;
    string prefix = "";
    int faceVectorIndex = 0;
    vector<GLint> faceVectors;
    vector<GLint> faceNormals;
    int crawler = 0;
    glm::vec3 p_min = glm::vec3(INT_MAX);
    glm::vec3 p_max = glm::vec3(INT_MIN);
    Triangles triangles;

    // read one line at a time
    while (getline(in_file, line))
    {
        ss.clear();
        ss.str(line);
        ss >> prefix;

        if (prefix == "v")
        { // vertex position
            ss >> val.x >> val.y >> val.z;
            vectors.push_back(val);
        }
        else if (prefix == "vn")
        {
            ss >> val.x >> val.y >> val.z;
            vectorNormals.push_back(val);
        }
        else if (prefix == "f")
        {
            int counter = 0;
            GLint tmp_glint;
            while (ss >> tmp_glint)
            {
                if (counter == 0)
                {
                    faceVectors.push_back(tmp_glint);
                }
                else if (counter == 1)
                {
                    faceNormals.push_back(tmp_glint);
                }
                if (ss.peek() == 47)
                {
                    counter++;
                    ss.ignore(2);
                }
                else if (ss.peek() == 32)
                {
                    counter = 0;
                    ss.ignore(1);
                }
                if (counter > 1)
                {
                    counter = 0;
                }
            }
            Face face;
            int size = faceVectors.size();
            face.p1 = vectors[faceVectors[size - 3] - 1] + startingPos;
            face.p2 = vectors[faceVectors[size - 2] - 1] + startingPos;
            face.p3 = vectors[faceVectors[size - 1] - 1] + startingPos;
            if (vectorNormals.size() > 0)
            {
                face.n1 = vectorNormals[faceNormals[size - 3] - 1];
                face.n2 = vectorNormals[faceNormals[size - 3] - 2];
                face.n3 = vectorNormals[faceNormals[size - 3] - 3];
            }
            glm::vec3 PQ = face.p2 - face.p1;
            glm::vec3 PR = face.p3 - face.p1;
            face.triangleNormal = glm::cross(PQ, PR);
            glm::mat3x3 triangle = glm::mat3x3(face.p1, face.p2, face.p3);
            face.det = glm::determinant(triangle);
            triangles.faces.push_back(face);
        }
    }
    return triangles;
}

#endif