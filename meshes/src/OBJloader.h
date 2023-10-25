
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
    glm::vec3 n1;
    glm::vec3 n2;
    glm::vec3 n3;
    glm::vec3 triangleNormal;
    float area;
};

static vector<Face>
loadOBJ(const char *file_name)
{
    vector<glm::vec3> vectors;
    vector<glm::vec3> vectorNormals;
    vector<Face> faces;
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
    vector<GLint> faceVectors;
    vector<GLint> faceNormals;

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
                else if (counter == 2)
                {
                    faceNormals.push_back(tmp_glint);
                }

                if (tmp_glint == '/')
                {
                    ++counter;
                    ss.ignore(1, '/');
                }
                else if (tmp_glint == ' ')
                {
                    ++counter;
                    ss.ignore(1, ' ');
                }
                if (counter > 2)
                {
                    counter = 0;
                }
            }
            Face face;
            face.p1 = vectors[faceVectors[0]];
            face.p2 = vectors[faceVectors[1]];
            face.p3 = vectors[faceVectors[2]];
            if (vectorNormals.size() > 0)
            {
                face.n1 = vectors[faceNormals[0]];
                face.n2 = vectors[faceNormals[1]];
                face.n3 = vectors[faceNormals[2]];
            }
            glm::vec3 PQ = face.p2 - face.p1;
            glm::vec3 PR = face.p3 - face.p1;
            face.triangleNormal = glm::cross(PQ, PR);
            face.area = glm::length(face.triangleNormal) / 2;
            faces.push_back(face);
        }
    }

    return faces;
}

#endif