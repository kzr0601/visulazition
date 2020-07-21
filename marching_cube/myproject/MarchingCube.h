//refer: https://github.com/CHINA-JIGE/MarchingCube

#pragma once
#ifndef __MARCHING_CUBE
#define __MARCHING_CUBE

#include <vector>
#include <map>
using namespace std;

struct Float3 {
	Float3() : x(0), y(0), z(0) {}
	Float3(float x, float y, float z) : x(x), y(y), z(z) {}

	Float3 operator+(const Float3& vec)
	{
		return Float3(x + vec.x, y + vec.y, z + vec.z);
	}

	Float3 operator-(const Float3& vec)
	{
		return Float3(x - vec.x, y - vec.y, z - vec.z);
	}

	Float3 operator*(const float& scale)
	{
		return Float3(scale*x, y*scale, z*scale);
	}

	Float3 operator*(const double& scale)
	{
		return Float3(scale*x, y*scale, z*scale);
	}

	float length()
	{
		return sqrt(x*x + y*y + z*z);
	}

	Float3 normalize()
	{
		float len = length();
		if (len == 0)
			return Float3(0, 0, 0);
		return Float3(x / len, y / len, z / len);
	}

	float x;
	float y;
	float z;
};

struct Int3 {
	Int3(): x(0), y(0), z(0){}
	Int3(int x, int y, int z) : x(x), y(y), z(z) {}

	void assign(int v, int index) {
		if (index == 0)
			x = v;
		else if (index == 1)
			y = v;
		else
			z = v;
	}

	bool operator < (const Int3& right) const  {
		if (x < right.x)
			return true;
		else if (x > right.x)
			return false;

		if (y < right.y)
			return true;
		else if (y > right.y)
			return false;

		if (z < right.z)
			return true;
		else if (z > right.z)
			return false;

		return false;
	}

	int x;
	int y;
	int z;
};

// at most four triangles in a marching cube
struct TriangleCase {
	int index[12];
};

//CT切层
struct CTSlice
{
	CTSlice() = delete;
	CTSlice(int pxWidth, int pxHeight) { pixelWidth = pxWidth; pixelHeight = pxHeight; };

	char GetPixel(int pixelX, int pixelY)
	{
		if (pixelX >= pixelWidth)
			pixelX--;
		if (pixelY >= pixelHeight)
			pixelY--;
		return bitArray.at(pixelY *pixelWidth + pixelX);
	}

	int pixelWidth;
	int pixelHeight;
	std::vector<char> bitArray;
};

// marching cube reconstructor
class MCReconstructor {
public:
	MCReconstructor();
	MCReconstructor(float xlen, float ylen, float zlen, int xcount, int ycount, int zcount, char value);

	bool LoadCTSlicer(std::string fileName, int pixelWidth=256, int pixelHeight=256, int depth=253); 
	void reconstruct(char inputv);
	void writeObjFile(std::string fileName);

	vector<Float3> meshpoint;
private:
	void calNormal();

	vector<CTSlice> mCTSlices;   // raw data
	vector<CTSlice> mCTSlices_binary;  // binaried data

	char value; //given value, 等值面的值
	vector< vector< vector<Float3> > > vertexNormal;

	// cube information
	float cubeXlen, cubeYlen, cubeZlen;
	int cubeCountX, cubeCountY, cubeCountZ;  
	
	// lookup table
	static const TriangleCase marchingCube_LookUp_Table[128];

	// save in obj
	int vertexCount;  // current vertex id
	map<Int3, int> vertex_id_table; // map [cubeCountX, cubeCountY, cubeCountZ] to vertex id
	map<int, Float3> vertex_table;  // map vertex id to vertex coordinate
	map<int, Float3> normal_table;  // map vertex id to normal
	vector<Int3> face_list;  // face list: [ [v0_index, v1_index, v2_index], ...]

	vector<Float3> vertex_list;
	vector<Float3> normal_list;
	//vector<Int3> face_list;
};

#endif // !__MARCHING_CUBE

