#define _USE_MATH_DEFINES

#include <cstdlib>
#include <iostream>
#include <memory>
#include <vector>
#include <cmath>

#include "glut37/glut.h"

using namespace std;

struct Vec4 {

    float x, y, z, w;
    Vec4(float x = 0, float y = 0, float z = 0, float w = 1) : x(x), y(y), z(z), w(w) {}

    Vec4 operator*(float matrix[4][4]) const {
        Vec4 result;
        result.x = x * matrix[0][0] + y * matrix[0][1] + z * matrix[0][2] + w * matrix[0][3];
        result.y = x * matrix[1][0] + y * matrix[1][1] + z * matrix[1][2] + w * matrix[1][3];
        result.z = x * matrix[2][0] + y * matrix[2][1] + z * matrix[2][2] + w * matrix[2][3];
        result.w = x * matrix[3][0] + y * matrix[3][1] + z * matrix[3][2] + w * matrix[3][3];
        return result;
    }

};


struct Vec3 {

    float x, y, z;
    Vec3(float x = 0, float y = 0, float z = 0) : x(x), y(y), z(z) {}

    Vec3 operator-(const Vec3& v) const { return Vec3(x - v.x, y - v.y, z - v.z); }
    Vec3 operator+(const Vec3& v) const { return Vec3(x + v.x, y + v.y, z + v.z); }

    Vec3 operator*(const Vec3& v) const { return Vec3(x * v.x, y * v.y, z * v.z); }
    Vec3 operator/(const Vec3& v) const { return Vec3(x / v.x, y / v.y, z / v.z); }

    Vec3 operator*(float d) const { return Vec3(x * d, y * d, z * d); }
    Vec3 operator/(float d) const { return Vec3(x / d, y / d, z / d); }

    Vec3 operator*(float matrix[4][4]) const {

        Vec3 result;
        result.x = x * matrix[0][0] + y * matrix[0][1] + z * matrix[0][2] + matrix[0][3];
        result.y = x * matrix[1][0] + y * matrix[1][1] + z * matrix[1][2] + matrix[1][3];
        result.z = x * matrix[2][0] + y * matrix[2][1] + z * matrix[2][2] + matrix[2][3];
        return result;
    }

    float dot(const Vec3& v) const { return x * v.x + y * v.y + z * v.z; }

    Vec3 normalize() const {
        float mag = sqrt(x * x + y * y + z * z);
        return Vec3(x / mag, y / mag, z / mag);
    }

    float norm() const {
        return std::sqrt(x * x + y * y + z * z);
    }

    Vec3 cross(const Vec3& v) const {
        return Vec3(
            y * v.z - z * v.y,
            z * v.x - x * v.z,
            x * v.y - y * v.x
        );
    }

    Vec3 pow(float p) const {

        return Vec3(std::pow(x, p), std::pow(y, p), std::pow(z, p));
    }
};

Vec3 operator*(float d, const Vec3& v) {
    return Vec3(d * v.x, d * v.y, d * v.z);
}

Vec3 operator/(float d, const Vec3& v) {
    return Vec3(d / v.x, d / v.y, d / v.z);
}

struct HitInfo {

    bool hit;
    float distance;
    Vec3 point;

    HitInfo(bool hit = false, float distance = -1, Vec3 point = NULL) : hit(hit), distance(distance), point(point) {}
};

float EPSILON = 1e-10;

int gWidth = 512;
int gHeight = 512;

float gL = -0.1;
float gR = 0.1;
float gB = -0.1;
float gT = 0.1;
float gN = -0.1;
float gF = -1000;

int gNumVertices = 0;    // Number of 3D vertices.
int gNumTriangles = 0;    // Number of triangles.
int* gIndexBuffer = NULL; // Vertex indices for the triangles.
Vec3* gVertexBuffer = NULL; // Array for vertices

Vec3* gPixels = new Vec3[gWidth * gHeight];

float gTransformMatrix[4][4] = {
    {2, 0, 0, 0},
    {0, 2, 0, 0},
    {0, 0, 2, -7},
    {0, 0, 0, 1}
};

float gPerspectiveProjectionMatrix[4][4] = {
   {2.0f * gN / (gR - gL), 0, (gR + gL) / (gR - gL), 0},
    {0, 2.0f * gN / (gT - gB), (gT + gB) / (gT - gB), 0},
    {0, 0, -(gF + gN) / (gF - gN), -2.0f * gF * gN / (gF - gN)},
    {0, 0, -1.0f, 0}
};

float gViewportTransformMatrix[4][4] = {
    {gWidth/ 2.0f, 0, 0, (gWidth - 1) / 2.0f},
    {0, gHeight / 2.0f, 0, (gHeight - 1) / 2.0f},
    {0, 0, (gF - gN) / 2.0f, (gF + gN) / 2.0f},
    {0, 0, 0, 1}
};

void createScene()
{
    int width = 32;
    int height = 16;

    float theta, phi;
    int t;

    gNumVertices = (height - 2) * width + 2;
    gNumTriangles = (height - 2) * (width - 1) * 2;

    // TODO: Allocate an array for gNumVertices vertices.
    gVertexBuffer = new Vec3[gNumVertices];
    gIndexBuffer = new int[3 * gNumTriangles];


    t = 0;
    for (int j = 1; j < height - 1; ++j)
    {
        for (int i = 0; i < width; ++i)
        {
            theta = (float)j / (height - 1) * M_PI;
            phi = (float)i / (width - 1) * M_PI * 2;

            float   x = sinf(theta) * cosf(phi);
            float   y = cosf(theta);
            float   z = -sinf(theta) * sinf(phi);

            // TODO: Set vertex t in the vertex array to {x, y, z}.
            // Applying the actual position of vertices using matrix multiplication
            gVertexBuffer[t] = Vec3(x, y, z) * gTransformMatrix;

            t++;
        }
    }

    // TODO: Set vertex t in the vertex array to {0, 1, 0}.
    // Applying the actual position of vertices using matrix multiplication
    gVertexBuffer[t] = Vec3(0, 1, 0) * gTransformMatrix;

    t++;

    // TODO: Set vertex t in the vertex array to {0, -1, 0}.
    // Applying the actual position of vertices using matrix multiplication
    gVertexBuffer[t] = Vec3(0, -1, 0) * gTransformMatrix;

    t++;

    t = 0;
    for (int j = 0; j < height - 3; ++j)
    {
        for (int i = 0; i < width - 1; ++i)
        {
            gIndexBuffer[t++] = j * width + i;
            gIndexBuffer[t++] = (j + 1) * width + (i + 1);
            gIndexBuffer[t++] = j * width + (i + 1);
            gIndexBuffer[t++] = j * width + i;
            gIndexBuffer[t++] = (j + 1) * width + i;
            gIndexBuffer[t++] = (j + 1) * width + (i + 1);
        }
    }
    for (int i = 0; i < width - 1; ++i)
    {
        gIndexBuffer[t++] = (height - 2) * width;
        gIndexBuffer[t++] = i;
        gIndexBuffer[t++] = i + 1;
        gIndexBuffer[t++] = (height - 2) * width + 1;
        gIndexBuffer[t++] = (height - 3) * width + (i + 1);
        gIndexBuffer[t++] = (height - 3) * width + i;
    }

    // The index buffer has now been generated. Here's how to use to determine
    // the vertices of a triangle. Suppose you want to determine the vertices
    // of triangle i, with 0 <= i < gNumTriangles. Define:
    //
    // k0 = gIndexBuffer[3*i + 0]
    // k1 = gIndexBuffer[3*i + 1]
    // k2 = gIndexBuffer[3*i + 2]
    //
    // Now, the vertices of triangle i are at positions k0, k1, and k2 (in that
    // order) in the vertex array (which you should allocate yourself at line
    // 27).
    //
    // Note that this assumes 0-based indexing of arrays (as used in C/C++,
    // Java, etc.) If your language uses 1-based indexing, you will have to
    // add 1 to k0, k1, and k2.


    for (int i = 0, l = gNumVertices; i < l; i++) {

        Vec4 vec4 = Vec4(gVertexBuffer[i].x, gVertexBuffer[i].y, gVertexBuffer[i].z) * gPerspectiveProjectionMatrix;

        gVertexBuffer[i].x = vec4.x / vec4.w;
        gVertexBuffer[i].y = vec4.y / vec4.w;
        gVertexBuffer[i].z = vec4.z / vec4.w;
    }

    for (int i = 0, l = gNumVertices; i < l; i++) {

        Vec4 vec4 = Vec4(gVertexBuffer[i].x, gVertexBuffer[i].y, gVertexBuffer[i].z) * gViewportTransformMatrix;

        gVertexBuffer[i].x = vec4.x / vec4.w;
        gVertexBuffer[i].y = vec4.y / vec4.w;
        gVertexBuffer[i].z = vec4.z / vec4.w;
    }
}


Vec3 getBarycentricCoords(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c) {

    Vec3 v0 = b - a, v1 = c - a, v2 = p - a;

    float d00 = v0.dot(v0);
    float d01 = v0.dot(v1);
    float d11 = v1.dot(v1);
    float d20 = v2.dot(v0);
    float d21 = v2.dot(v1);
    float denom = d00 * d11 - d01 * d01;

    if (abs(denom) < EPSILON) {
        return { -1, -1, -1 };
    }

    float v = (d11 * d20 - d01 * d21) / denom;
    float w = (d00 * d21 - d01 * d20) / denom;
    float u = 1.0f - v - w;

    return { u, v, w };
}


void render() {

    for (int i = 0; i < gNumTriangles; i++) {

    	const Vec3& v0 = gVertexBuffer[gIndexBuffer[i * 3 + 0]];
        const Vec3& v1 = gVertexBuffer[gIndexBuffer[i * 3 + 1]];
        const Vec3& v2 = gVertexBuffer[gIndexBuffer[i * 3 + 2]];

        int minX = max(0, (int)min({ v0.x, v1.x, v2.x }));
        int maxX = min(gWidth - 1, (int)max({ v0.x, v1.x, v2.x }));

        int minY = max(0, (int)min({ v0.y, v1.y, v2.y }));
        int maxY = min(gHeight - 1, (int)max({ v0.y, v1.y, v2.y }));

        for (int x = minX; x <= maxX; x++) {
            for (int y = minY; y <= maxY; y++) {

            	Vec3 p(x, y, 0);
                Vec3 v0xy(v0.x, v0.y, 0);
                Vec3 v1xy(v1.x, v1.y, 0);
                Vec3 v2xy(v2.x, v2.y, 0);

                Vec3 barycentricCoords = getBarycentricCoords(p, v0xy, v1xy, v2xy);

                if (barycentricCoords.x >= 0 && barycentricCoords.y >= 0 && barycentricCoords.z >= 0) {

                    int idx = y * gWidth + x;
                    gPixels[idx] = Vec3(1, 1, 1);
                }
            }
        }
    }

    glDrawPixels(gWidth, gHeight, GL_RGB, GL_FLOAT, gPixels);
}

void display() {

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glEnable(GL_DEPTH_TEST);
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    render();

    glutSwapBuffers();
}

int main(int argc, char** argv) {

    glutInit(&argc, argv);
    glutInitDisplayMode(GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);
    glutInitWindowSize(gWidth, gHeight);
    glutInitWindowPosition(100, 100);

    glutCreateWindow("HW4 TNR");

    createScene();

    glutDisplayFunc(display);

    glutMainLoop();

    return 0;
}