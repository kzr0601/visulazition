#include "MarchingCube.h"
#include <fstream>
#include <iostream>
using namespace std;

MCReconstructor::MCReconstructor()
{
	vertexCount = 1;
	value = 130;
}

MCReconstructor::MCReconstructor(float xlen, float ylen, float zlen, int xcount, int ycount, int zcount, char value)
{
	cubeXlen = xlen;
	cubeYlen = ylen;
	cubeZlen = zlen;

	cubeCountX = xcount;
	cubeCountY = ycount;
	cubeCountZ = zcount;

	value = value;
	vertexCount = 1;
}

void MCReconstructor::calNormal()
{
	vector<Float3> vvv;
	Float3 f(1.0, 2, 3);
	vvv.resize(cubeCountX+1, f);
	vector<vector<Float3> > vv;
	vv.resize(cubeCountZ+1, vvv);
	vertexNormal.resize(cubeCountY+1, vv);

	CTSlice down = mCTSlices.at(0);
	CTSlice up = mCTSlices.at(1);
	for (int r = 0; r <= cubeCountZ; r++) {
		for (int c = 0; c <= cubeCountX; c++) {
			float gx, gy, gz;

			int px = 1.0 * c / cubeCountX * down.pixelWidth;
			int py = 1.0 * r / cubeCountZ * down.pixelHeight;
			int density = down.GetPixel(px, py);

			int d_y_add_1 = up.GetPixel(px, py);
			gy = 1.0 *(d_y_add_1 - density) / cubeYlen;

			if (c == 0) {
				int d_xjia1 = down.GetPixel(1.0 * (c + 1) / cubeCountX * down.pixelWidth, 1.0 * r / cubeCountZ * down.pixelHeight);
				gx = 1.0 * (d_xjia1 - density) / cubeXlen;
			}
			else if (c == cubeCountX) {
				int d_xjian1 = down.GetPixel(1.0 * (c - 1) / cubeCountX * down.pixelWidth, 1.0 * r / cubeCountZ * down.pixelHeight);
				gx = 1.0 * (density - d_xjian1) / cubeXlen;
			}
			else {
				int d_x_add_1 = down.GetPixel(1.0 * (c + 1) / cubeCountX * down.pixelWidth, 1.0 * r / cubeCountZ * down.pixelHeight);
				int d_x_minus_1 = down.GetPixel(1.0 * (c - 1) / cubeCountX * down.pixelWidth, 1.0 * r / cubeCountZ * down.pixelHeight);
				gx = 1.0 * (d_x_add_1 - d_x_minus_1) / (2*cubeXlen);
			}

			if (r == 0) {
				int d_zjia1 = down.GetPixel(1.0 * c / cubeCountX * down.pixelWidth, 1.0 * (r+1) / cubeCountZ * down.pixelHeight);
				gz = 1.0 * (d_zjia1 - density) / cubeZlen;
			}
			else if (r == cubeCountZ) {
				int d_zjian1 = down.GetPixel(1.0 * c / cubeCountX * down.pixelWidth, 1.0 * (r - 1) / cubeCountZ * down.pixelHeight);
				gz = 1.0 * (density - d_zjian1) / cubeZlen;
			}
			else {
				int d_z_add_1 = down.GetPixel(1.0 * c / cubeCountX * down.pixelWidth, 1.0 * (r+1) / cubeCountZ * down.pixelHeight);
				int d_z_minus_1 = down.GetPixel(1.0 * c / cubeCountX * down.pixelWidth, 1.0 * (r-1) / cubeCountZ * down.pixelHeight);
				gz = 1.0 * (d_z_add_1 - d_z_minus_1) / (2 * cubeZlen);
			}

			vertexNormal[0][r][c] = Float3(gx, gy, gz).normalize();
		}
	}

	down = mCTSlices.at(cubeCountY);
	CTSlice down_down = mCTSlices.at(cubeCountY - 1);
	for (int r = 0; r <= cubeCountZ; r++) {
		for (int c = 0; c <= cubeCountX; c++) {
			float gx, gy, gz;

			int px = 1.0 * c / cubeCountX * down.pixelWidth;
			int py = 1.0 * r / cubeCountZ * down.pixelHeight;
			int density = down.GetPixel(px, py);

			int d_y_minus_1 = down_down.GetPixel(px, py);
			gy = 1.0 *(density - d_y_minus_1) / cubeYlen;

			if (c == 0) {
				int d_xjia1 = down.GetPixel(1.0 * (c + 1) / cubeCountX * down.pixelWidth, 1.0 * r / cubeCountZ * down.pixelHeight);
				gx = 1.0 * (d_xjia1 - density) / cubeXlen;
			}
			else if (c == cubeCountX ) {
				int d_xjian1 = down.GetPixel(1.0 * (c - 1) / cubeCountX * down.pixelWidth, 1.0 * r / cubeCountZ * down.pixelHeight);
				gx = 1.0 * (density - d_xjian1) / cubeXlen;
			}
			else {
				int d_x_add_1 = down.GetPixel(1.0 * (c + 1) / cubeCountX * down.pixelWidth, 1.0 * r / cubeCountZ * down.pixelHeight);
				int d_x_minus_1 = down.GetPixel(1.0 * (c - 1) / cubeCountX * down.pixelWidth, 1.0 * r / cubeCountZ * down.pixelHeight);
				gx = 1.0 * (d_x_add_1 - d_x_minus_1) / (2 * cubeXlen);
			}

			if (r == 0) {
				int d_zjia1 = down.GetPixel(1.0 * c / cubeCountX * down.pixelWidth, 1.0 * (r + 1) / cubeCountZ * down.pixelHeight);
				gz = 1.0 * (d_zjia1 - density) / cubeZlen;
			}
			else if (r == cubeCountZ ) {
				int d_zjian1 = down.GetPixel(1.0 * c / cubeCountX * down.pixelWidth, 1.0 * (r - 1) / cubeCountZ * down.pixelHeight);
				gz = 1.0 * (density - d_zjian1) / cubeZlen;
			}
			else {
				int d_z_add_1 = down.GetPixel(1.0 * c / cubeCountX * down.pixelWidth, 1.0 * (r + 1) / cubeCountZ * down.pixelHeight);
				int d_z_minus_1 = down.GetPixel(1.0 * c / cubeCountX * down.pixelWidth, 1.0 * (r - 1) / cubeCountZ * down.pixelHeight);
				gz = 1.0 * (d_z_add_1 - d_z_minus_1) / (2 * cubeZlen);
			}

			vertexNormal[cubeCountY][r][c] = Float3(gx, gy, gz).normalize();
		}
	}

	for (int i = 1; i < cubeCountY; i++) {
		CTSlice down = mCTSlices.at(i);
		CTSlice up = mCTSlices.at(i+1);
		CTSlice down_down = mCTSlices.at(i - 1);
		for (int r = 0; r <= cubeCountZ; r++) {
			for (int c = 0; c <= cubeCountX; c++) {
				float gx, gy, gz;

				int px = 1.0 * c / cubeCountX * down.pixelWidth;
				int py = 1.0 * r / cubeCountZ * down.pixelHeight;
				int density = down.GetPixel(px, py);

				int d_y_add_1 = up.GetPixel(px, py);
				int d_y_minus_1 = down_down.GetPixel(px, py);
				gy = 1.0 *(d_y_add_1 - d_y_minus_1) / (2*cubeYlen);

				if (c == 0) {
					int d_xjia1 = down.GetPixel(1.0 * (c + 1) / cubeCountX * down.pixelWidth, 1.0 * r / cubeCountZ * down.pixelHeight);
					gx = 1.0 * (d_xjia1 - density) / cubeXlen;
				}
				else if (c == cubeCountX ) {
					int d_xjian1 = down.GetPixel(1.0 * (c - 1) / cubeCountX * down.pixelWidth, 1.0 * r / cubeCountZ * down.pixelHeight);
					gx = 1.0 * (density - d_xjian1) / cubeXlen;
				}
				else {
					int d_x_add_1 = down.GetPixel(1.0 * (c + 1) / cubeCountX * down.pixelWidth, 1.0 * r / cubeCountZ * down.pixelHeight);
					int d_x_minus_1 = down.GetPixel(1.0 * (c - 1) / cubeCountX * down.pixelWidth, 1.0 * r / cubeCountZ * down.pixelHeight);
					gx = 1.0 * (d_x_add_1 - d_x_minus_1) / (2 * cubeXlen);
				}

				if (r == 0) {
					int d_zjia1 = down.GetPixel(1.0 * c / cubeCountX * down.pixelWidth, 1.0 * (r + 1) / cubeCountZ * down.pixelHeight);
					gz = 1.0 * (d_zjia1 - density) / cubeZlen;
				}
				else if (r == cubeCountZ ) {
					int d_zjian1 = down.GetPixel(1.0 * c / cubeCountX * down.pixelWidth, 1.0 * (r - 1) / cubeCountZ * down.pixelHeight);
					gz = 1.0 * (density - d_zjian1) / cubeZlen;
				}
				else {
					int d_z_add_1 = down.GetPixel(1.0 * c / cubeCountX * down.pixelWidth, 1.0 * (r + 1) / cubeCountZ * down.pixelHeight);
					int d_z_minus_1 = down.GetPixel(1.0 * c / cubeCountX * down.pixelWidth, 1.0 * (r - 1) / cubeCountZ * down.pixelHeight);
					gz = 1.0 * (d_z_add_1 - d_z_minus_1) / (2 * cubeZlen);
				}

				vertexNormal[i][r][c] = Float3(gx, gy, gz).normalize();
			}
		}
	}

	cout << "calculate normal successful" << endl;
}

bool MCReconstructor::LoadCTSlicer(std::string fileName, int pixelWidth, int pixelHeight, int depth)
{
	cubeCountY = depth-1;

	CTSlice slice(pixelWidth, pixelHeight);
	mCTSlices.resize(depth, slice);
	mCTSlices_binary.resize(depth, slice);

	std::ifstream inFile(fileName, std::ios::binary);
	if (!inFile.is_open())
	{
		cout << "File Load Error!! " << endl;
		return false;
	}
	for (int i = 0; i < depth; i++) {
		char* fileBuff = new char[pixelWidth*pixelHeight];
		inFile.read(fileBuff, pixelWidth*pixelHeight);

		mCTSlices.at(i).bitArray.assign(fileBuff, fileBuff + pixelWidth*pixelHeight);

		for (int r = 0; r < pixelHeight; r++) {
			for (int c = 0; c < pixelWidth; c++) {
				if (fileBuff[r*pixelWidth + c] >= value)
					fileBuff[r*pixelWidth + c] = 1;
				else
					fileBuff[r*pixelWidth + c] = -1;
				//fileBuff[r*pixelWidth + c] = fileBuff[r*pixelWidth + c] - value;
			}
		}

		mCTSlices_binary.at(i).bitArray.assign(fileBuff, fileBuff + pixelWidth*pixelHeight);
	}

	inFile.close();

	calNormal();

	cout << "load CTSlice successful" << endl;
	return true;
}

void MCReconstructor::reconstruct(char inputv)
{
	value = inputv;

	// binary data with inputv
	for (int i = 0; i < mCTSlices_binary.size(); i++) {
		CTSlice ctslice = mCTSlices.at(i);
		int pixelHeight = ctslice.pixelHeight;
		int pixelWidth  = ctslice.pixelWidth;
		vector<char> fileBuff = ctslice.bitArray;

		for (int r = 0; r < pixelHeight; r++) {
			for (int c = 0; c < pixelWidth; c++) {
				if (fileBuff[r*pixelWidth + c] >= value)
					fileBuff[r*pixelWidth + c] = 1;
				else
					fileBuff[r*pixelWidth + c] = -1;
			}
		}

		mCTSlices_binary.at(i).bitArray = fileBuff;
	}

	// reconstruct face and vertex normal
	for (int cubeIndexY = 0; cubeIndexY < cubeCountY; cubeIndexY++) {
		for (int cubeIndexX = 0; cubeIndexX < cubeCountX; cubeIndexX++) {
			for (int cubeIndexZ = 0; cubeIndexZ < cubeCountZ; cubeIndexZ++) {
				Float3 v1 = Float3(cubeIndexX * cubeXlen, cubeIndexY * cubeYlen, cubeIndexZ * cubeZlen);

				// 8 vertexs of cube
				Float3 v[8] = {
					v1 + Float3(0, 0, 0),
					v1 + Float3(cubeXlen, 0, 0),
					v1 + Float3(cubeXlen, cubeYlen, 0),
					v1 + Float3(0, cubeYlen, 0),
					v1 + Float3(0, 0, cubeZlen),
					v1 + Float3(cubeXlen, 0, cubeZlen),
					v1 + Float3(cubeXlen, cubeYlen, cubeZlen),
					v1 + Float3(0, cubeYlen, cubeZlen)
				};

				// vertex value
				char vertex_value[8];

				// calculate case index
				int caseIndex = 0;
				for (int i = 0; i<8; ++i)
				{
					int currentCubeIndexY = 0;
					if (i == 0 || i == 1 || i == 4 || i == 5)
					{
						currentCubeIndexY = cubeIndexY; //立方体下表面
					}
					else
					{
						currentCubeIndexY = cubeIndexY + 1; //立方体上表面
					}
					
					CTSlice&  currentSlice = mCTSlices_binary.at(currentCubeIndexY);

					int pixelCoordX = v[i].x / cubeXlen / cubeCountX * currentSlice.pixelWidth;
					int pixelCoordY = v[i].z / cubeZlen / cubeCountZ * currentSlice.pixelHeight;
					if (pixelCoordX == currentSlice.pixelWidth)pixelCoordX--;
					if (pixelCoordY == currentSlice.pixelHeight)pixelCoordY--;

					if (currentSlice.GetPixel(pixelCoordX, pixelCoordY) == -1)
					{
						//设置第i个二进制位
						caseIndex |= (1 << i);
					}

					vertex_value[i] = mCTSlices.at(currentCubeIndexY).GetPixel(pixelCoordX, pixelCoordY);
				}

				// complementary case
				if (caseIndex >= 128)
					caseIndex = 255 - caseIndex;
				// rotational-symmetry?

				// median point on 12 edges
				Float3 medianPointOnEdge[12] =
				{
					(v[0] + v[1]) * 0.5f,
					(v[1] + v[2]) * 0.5f,
					(v[3] + v[2]) * 0.5f,
					(v[0] + v[3]) * 0.5f,
					(v[4] + v[5]) * 0.5f,
					(v[5] + v[6]) * 0.5f,
					(v[7] + v[6]) * 0.5f,
					(v[4] + v[7]) * 0.5f,
					(v[0] + v[4]) * 0.5f,
					(v[1] + v[5]) * 0.5f,
					(v[3] + v[7]) * 0.5f,
					(v[2] + v[6]) * 0.5f
				};

				// "value = given value" point on 12 edges
				// interpolation ratio
				float ratio[12] = 
				{
					(vertex_value[0] - value)*(vertex_value[1] - value)<0 ? 1.0*(inputv - vertex_value[0]) / (vertex_value[1] - vertex_value[0]) : 0,
					(vertex_value[1] - value)*(vertex_value[2] - value)<0 ? 1.0*(inputv - vertex_value[1]) / (vertex_value[2] - vertex_value[1]) : 0,
					(vertex_value[2] - value)*(vertex_value[3] - value)<0 ? 1.0*(inputv - vertex_value[2]) / (vertex_value[3] - vertex_value[2]) : 0,
					(vertex_value[0] - value)*(vertex_value[3] - value)<0 ? 1.0*(inputv - vertex_value[0]) / (vertex_value[3] - vertex_value[0]) : 0,
					(vertex_value[4] - value)*(vertex_value[5] - value)<0 ? 1.0*(inputv - vertex_value[4]) / (vertex_value[5] - vertex_value[4]) : 0,
					(vertex_value[6] - value)*(vertex_value[5] - value)<0 ? 1.0*(inputv - vertex_value[5]) / (vertex_value[6] - vertex_value[5]) : 0,
					(vertex_value[6] - value)*(vertex_value[7] - value)<0 ? 1.0*(inputv - vertex_value[6]) / (vertex_value[7] - vertex_value[6]) : 0,
					(vertex_value[4] - value)*(vertex_value[7] - value)<0 ? 1.0*(inputv - vertex_value[4]) / (vertex_value[7] - vertex_value[4]) : 0,
					(vertex_value[0] - value)*(vertex_value[4] - value)<0 ? 1.0*(inputv - vertex_value[0]) / (vertex_value[4] - vertex_value[0]) : 0,
					(vertex_value[1] - value)*(vertex_value[5] - value)<0 ? 1.0*(inputv - vertex_value[1]) / (vertex_value[5] - vertex_value[1]) : 0,
					(vertex_value[3] - value)*(vertex_value[7] - value)<0 ? 1.0*(inputv - vertex_value[3]) / (vertex_value[7] - vertex_value[3]) : 0,
					(vertex_value[2] - value)*(vertex_value[6] - value)<0 ? 1.0*(inputv - vertex_value[2]) / (vertex_value[6] - vertex_value[2]) : 0
				};

				//{
				//	(vertex_value[0] - value)*(vertex_value[1] - value)<0 ? v[0] + 1.0*(inputv-vertex_value[0])/(vertex_value[1]-vertex_value[0])*(v[1]-v[0]) : v[0],
				//	(vertex_value[1] - value)*(vertex_value[2] - value)<0 ? v[1] + 1.0*(inputv - vertex_value[1]) / (vertex_value[2] - vertex_value[1])*(v[2] - v[1]) : v[1],					
				//	(vertex_value[2] - value)*(vertex_value[3] - value)<0 ? v[2] + 1.0*(inputv - vertex_value[2]) / (vertex_value[3] - vertex_value[2])*(v[3] - v[2]) : v[2],
				//	(vertex_value[0] - value)*(vertex_value[3] - value)<0 ? v[0] + 1.0*(inputv - vertex_value[0]) / (vertex_value[3] - vertex_value[0])*(v[3] - v[0]) : v[0],
				//	(vertex_value[4] - value)*(vertex_value[5] - value)<0 ? v[4] + 1.0*(inputv - vertex_value[4]) / (vertex_value[5] - vertex_value[4])*(v[5] - v[4]) : v[4],
				//	(vertex_value[6] - value)*(vertex_value[5] - value)<0 ? v[6] + 1.0*(inputv - vertex_value[6]) / (vertex_value[5] - vertex_value[6])*(v[5] - v[6]) : v[6],
				//	(vertex_value[6] - value)*(vertex_value[7] - value)<0 ? v[6] + 1.0*(inputv - vertex_value[6]) / (vertex_value[7] - vertex_value[6])*(v[7] - v[6]) : v[6],
				//	(vertex_value[4] - value)*(vertex_value[7] - value)<0 ? v[4] + 1.0*(inputv - vertex_value[4]) / (vertex_value[7] - vertex_value[4])*(v[7] - v[4]) : v[4],
				//	(vertex_value[0] - value)*(vertex_value[4] - value)<0 ? v[0] + 1.0*(inputv - vertex_value[0]) / (vertex_value[4] - vertex_value[0])*(v[4] - v[0]) : v[0],
				//	(vertex_value[1] - value)*(vertex_value[5] - value)<0 ? v[1] + 1.0*(inputv - vertex_value[1]) / (vertex_value[5] - vertex_value[1])*(v[5] - v[1]) : v[1],
				//	(vertex_value[3] - value)*(vertex_value[7] - value)<0 ? v[3] + 1.0*(inputv - vertex_value[3]) / (vertex_value[7] - vertex_value[3])*(v[7] - v[3]) : v[3],
				//	(vertex_value[2] - value)*(vertex_value[6] - value)<0 ? v[2] + 1.0*(inputv - vertex_value[2]) / (vertex_value[6] - vertex_value[2])*(v[6] - v[2]) : v[2]
				//};
				
				// calculate normal of 8 vertex
				Float3 normal_vertex[8] = {
					vertexNormal[cubeIndexY][cubeIndexZ][cubeIndexX],
					vertexNormal[cubeIndexY][cubeIndexZ][cubeIndexX+1],
					vertexNormal[cubeIndexY+1][cubeIndexZ][cubeIndexX + 1],
					vertexNormal[cubeIndexY][cubeIndexZ][cubeIndexX],
					vertexNormal[cubeIndexY][cubeIndexZ+1][cubeIndexX],
					vertexNormal[cubeIndexY][cubeIndexZ+1][cubeIndexX + 1],
					vertexNormal[cubeIndexY+1][cubeIndexZ+1][cubeIndexX + 1],
					vertexNormal[cubeIndexY+1][cubeIndexZ+1][cubeIndexX]
				};

				Float3 normalmedianPointOnEdge[12] =
				{
					(normal_vertex[0] + normal_vertex[1]) * 0.5f,
					(normal_vertex[1] + normal_vertex[2]) * 0.5f,
					(normal_vertex[3] + normal_vertex[2]) * 0.5f,
					(normal_vertex[0] + normal_vertex[3]) * 0.5f,
					(normal_vertex[4] + normal_vertex[5]) * 0.5f,
					(normal_vertex[5] + normal_vertex[6]) * 0.5f,
					(normal_vertex[7] + normal_vertex[6]) * 0.5f,
					(normal_vertex[4] + normal_vertex[7]) * 0.5f,
					(normal_vertex[0] + normal_vertex[4]) * 0.5f,
					(normal_vertex[1] + normal_vertex[5]) * 0.5f,
					(normal_vertex[3] + normal_vertex[7]) * 0.5f,
					(normal_vertex[2] + normal_vertex[6]) * 0.5f
				};

				// look up table
				if (caseIndex == 0)
					continue;
				TriangleCase tc = marchingCube_LookUp_Table[caseIndex];

				//record obj vertex/face/normal
				for (int i = 0; i < 12; i += 3) {
					if (tc.index[i] == -1)
						break;
					for (int k = 0; k < 3; k++) {
						vertex_list.push_back(medianPointOnEdge[tc.index[i + k]]);
						normal_list.push_back(normalmedianPointOnEdge[tc.index[i + k]]);
					}
				}


				for (int i = 0; i < 12; i += 3)
				{
					if (tc.index[i] == -1)break;
					//减的那个东西是为了让模型靠近原点
					meshpoint.push_back(medianPointOnEdge[tc.index[i]]);
					meshpoint.push_back(medianPointOnEdge[tc.index[i+1]]);
					meshpoint.push_back(medianPointOnEdge[tc.index[i+2]]);
				}
			}
		}
	}

	cout << "reconstruct successful" << endl;
}

void MCReconstructor::writeObjFile(std::string filename)
{
	// create obj file
	ofstream outfile;
	outfile.open(filename, ios::out);

	// write vertex
	for (int i = 0; i < vertex_list.size(); i++) {
		outfile << "v " << vertex_list[i].x << " " << vertex_list[i].y << " " << vertex_list[i].z << "\n";
	}
	outfile << "\n";

	// write normal
	for (int i = 0; i < vertex_list.size(); i++) {
		outfile << "vn " << normal_list[i].x << " " << normal_list[i].y << " " << normal_list[i].z << "\n";
	}
	outfile << "\n";

	//write face
	int face_len = vertex_list.size() / 3;

	cout << meshpoint.size() << " " << vertex_list.size() << " " << face_len << endl;

	for (int i = 0; i < face_len; i++) {
		outfile << "f " << 3*i+1 << "//" << 3 * i + 1 << " " << 3 * i + 2 << "//" << 3 * i + 2 << " " << 3 * i + 3 << "//" << 3 * i + 3 << "\n";
	}
	outfile << "\n";

	//cout << vertexCount << " = 1 + " << vertex_table.size() << endl;
	//// write vertex
	//for (int i = 1; i < vertexCount; i++) {
	//	outfile << "v " << vertex_table[i].x << " " << vertex_table[i].y << " " << vertex_table[i].z << "\n";
	//}
	//outfile << "\n";

	//// write normal
	//for (int i = 1; i < vertexCount; i++) {
	//	outfile << "vn " << normal_table[i].x << " " << normal_table[i].y << " " << normal_table[i].z << "\n";
	//}
	//outfile << "\n";

	////write face
	//for (int i = 0; i < face_list.size(); i++) {
	//	Int3 int3 = face_list[i];
	//	outfile << "f " << int3.x << "//" << int3.x << " " << int3.y << "//" << int3.y << " " << int3.z << "//" << int3.z << "\n";
	//}
	//outfile << "\n";

	outfile.close();

	cout << "write obj file successful" << endl;
}

const TriangleCase MCReconstructor::marchingCube_LookUp_Table[128] =
{
	{ { -1, -1, -1, -1, -1, -1, -1, -1} }, /* 0 0 */

	{ { 0, 3, 8, -1, -1, -1, -1, -1, -1, -1, -1, -1} }, /* 1 1 */
	{ { 0, 9, 1, -1, -1, -1, -1, -1, -1, -1, -1, -1} }, /* 2 1 */
	{ { 1, 3, 8, 9, 1, 8, -1, -1} }, /* 3 2 */
	{ { 1, 11, 2, -1, -1, -1, -1, -1, -1, -1, -1, -1} }, /* 4 1 */
	{ { 0, 3, 8, 1, 11, 2, -1, -1, -1, -1, -1, -1} }, /* 5 3 */
	{ { 9, 11, 2, 0, 9, 2, -1, -1, -1, -1, -1, -1} }, /* 6 2 */
	{ { 2, 3, 8, 2, 8, 11, 11, 8, 9, -1, -1, -1} }, /* 7 5 */
	{ { 3, 2, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1} }, /* 8 1 */
	{ { 0, 2, 10, 8, 0, 10, -1, -1, -1, -1, -1, -1} }, /* 9 2 */
	{ { 1, 0, 9, 2, 10, 3, -1, -1, -1, -1, -1, -1} }, /* 10 3 */
	{ { 1, 2, 10, 1, 10, 9, 9, 10, 8, -1, -1, -1} }, /* 11 5 */
	{ { 3, 1, 11, 10, 3, 11, -1, -1, -1, -1, -1, -1} }, /* 12 2 */
	{ { 0, 1, 11, 0, 11, 8, 8, 11, 10, -1, -1, -1} }, /* 13 5 */
	{ { 3, 0, 9, 3, 9, 10, 10, 9, 11, -1, -1, -1} }, /* 14 5 */
	{ { 9, 11, 8, 11, 10, 8, -1, -1, -1, -1, -1, -1} }, /* 15 8 */
	{ { 4, 8, 7, -1, -1, -1, -1, -1} }, /* 16 1 */
	{ { 4, 0, 3, 7, 4, 3, -1, -1} }, /* 17 2 */
	{ { 0, 9, 1, 8, 7, 4, -1, -1, -1, -1, -1, -1} }, /* 18 3 */
	{ { 4, 9, 1, 4, 1, 7, 7, 1, 3, -1, -1, -1 } }, /* 19 5 */
	{ { 1, 11, 2, 8, 7, 4, -1, -1, -1, -1, -1, -1} }, /* 20 4 */
	{ { 3, 7, 4, 3, 4, 0, 1, 11, 2, -1, -1, -1} }, /* 21 7 */
	{ { 9, 11, 2, 9, 2, 0, 8, 7, 4, -1, -1, -1 } }, /* 22 7 */
	{ { 2, 9, 11, 2, 7, 9, 2, 3, 7, 7, 4, 9 } }, /* 23 14 */
	{ { 8, 7, 4, 3, 2, 10, -1, -1} }, /* 24 3 */
	{ { 10, 7, 4, 10, 4, 2, 2, 4, 0, -1, -1, -1 } }, /* 25 5 */
	{ { 9, 1, 0, 8, 7, 4, 2, 10, 3, -1, -1, -1 } }, /* 26 6 */
	{ { 4, 10, 7, 9, 10, 4, 9, 2, 10, 9, 1, 2 } }, /* 27 9 */
	{ { 3, 1, 11, 3, 11, 10, 7, 4, 8, -1, -1, -1 } }, /* 28 7 */
	{ { 1, 11, 10, 1, 10, 4, 1, 4, 0, 7, 4, 10 } }, /* 29 11 */
	{ { 4, 8, 7, 9, 10, 0, 9, 11, 10, 10, 3, 0 } }, /* 30 12 */
	{ { 4, 10, 7, 4, 9, 10, 9, 11, 10, -1, -1, -1 } }, /* 31 5 */
	{ { 9, 4, 5, -1, -1, -1, -1, -1, -1, -1, -1, -1} }, /* 32 1 */
	{ { 9, 4, 5, 0, 3, 8, -1, -1, -1, -1, -1, -1} }, /* 33 3 */
	{ { 0, 4, 5, 1, 0, 5, -1, -1, -1, -1, -1, -1} }, /* 34 2 */
	{ { 8, 4, 5, 8, 5, 3, 3, 5, 1, -1, -1, -1} }, /* 35 5 */
	{ { 1, 11, 2, 9, 4, 5, -1, -1, -1, -1, -1, -1} }, /* 36 3 */
	{ { 3, 8, 0, 1, 11, 2, 4, 5, 9, -1, -1, -1} }, /* 37 6 */
	{ { 5, 11, 2, 5, 2, 4, 4, 2, 0, -1, -1, -1} }, /* 38 5 */
	{ { 2, 5, 11, 3, 5, 2, 3, 4, 5, 3, 8, 4} }, /* 39 9 */
	{ { 9, 4, 5, 2, 10, 3, -1, -1, -1, -1, -1, -1} }, /* 40 4 */
	{ { 0, 2, 10, 0, 10, 8, 4, 5, 9, -1, -1, -1} }, /* 41 7 */
	{ { 0, 4, 5, 0, 5, 1, 2, 10, 3, -1, -1, -1} }, /* 42 7 */
	{ { 2, 5, 1, 2, 8, 5, 2, 10, 8, 4, 5, 8} }, /* 43 11 */
	{ { 11, 10, 3, 11, 3, 1, 9, 4, 5, -1, -1, -1} }, /* 44 7 */
	{ { 4, 5, 9, 0, 1, 8, 8, 1, 11, 8, 11, 10} }, /* 45 12 */
	{ { 5, 0, 4, 5, 10, 0, 5, 11, 10, 10, 3, 0} }, /* 46 14 */
	{ { 5, 8, 4, 5, 11, 8, 11, 10, 8, -1, -1, -1} }, /* 47 5 */
	{ { 9, 8, 7, 5, 9, 7, -1, -1, -1, -1, -1, -1} }, /* 48 2 */
	{ { 9, 0, 3, 9, 3, 5, 5, 3, 7, -1, -1, -1} }, /* 49 5 */
	{ { 0, 8, 7, 0, 7, 1, 1, 7, 5, -1, -1, -1} }, /* 50 5 */
	{ { 1, 3, 5, 3, 7, 5, -1, -1, -1, -1, -1, -1} }, /* 51 8 */
	{ { 9, 8, 7, 9, 7, 5, 11, 2, 1, -1, -1, -1} }, /* 52 7 */
	{ { 11, 2, 1, 9, 0, 5, 5, 0, 3, 5, 3, 7} }, /* 53 12 */
	{ { 8, 2, 0, 8, 5, 2, 8, 7, 5, 11, 2, 5} }, /* 54 11 */
	{ { 2, 5, 11, 2, 3, 5, 3, 7, 5, -1, -1, -1} }, /* 55 5 */
	{ { 7, 5, 9, 7, 9, 8, 3, 2, 10, -1, -1, -1} }, /* 56 7 */
	{ { 9, 7, 5, 9, 2, 7, 9, 0, 2, 2, 10, 7} }, /* 57 14 */
	{ { 2, 10, 3, 0, 8, 1, 1, 8, 7, 1, 7, 5} }, /* 58 12 */
	{ { 10, 1, 2, 10, 7, 1, 7, 5, 1, -1, -1, -1} }, /* 59 5 */
	{ { 9, 8, 5, 8, 7, 5, 11, 3, 1, 11, 10, 3} }, /* 60 10 */
	{ { 0, 1, 9, 7, 10, 11, 5, 7, 11, -1, -1, -1} }, /* 61 7 */
	{ { 0, 3, 9, 7, 10, 11, 5, 7, 11, -1, -1, -1 } }, /* 62 7 */
	{ { 10, 5, 11, 7, 5, 10, -1, -1, -1, -1, -1, -1} }, /* 63 2 */
	{ { 11, 5, 6, -1, -1, -1, -1, -1, -1, -1, -1, -1} }, /* 64 1 */
	{ { 0, 3, 8, 5, 6, 11, -1, -1, -1, -1, -1, -1} }, /* 65 4 */
	{ { 9, 1, 0, 5, 6, 11, -1, -1, -1, -1, -1, -1} }, /* 66 3 */
	{ { 1, 3, 8, 1, 8, 9, 5, 6, 11, -1, -1, -1} }, /* 67 7 */
	{ { 1, 5, 6, 2, 1, 6, -1, -1, -1, -1, -1, -1} }, /* 68 2 */
	{ { 1, 5, 6, 1, 6, 2, 3, 8, 0, -1, -1, -1} }, /* 69 7 */
	{ { 9, 5, 6, 9, 6, 0, 0, 6, 2, -1, -1, -1} }, /* 70 5 */
	{ { 5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2} }, /* 71 11 */
	{ { 2, 10, 3, 11, 5, 6, -1, -1, -1, -1, -1, -1} }, /* 72 3 */
	{ { 10, 8, 0, 10, 0, 2, 11, 5, 6, -1, -1, -1} }, /* 73 7 */
	{ { 0, 9, 1, 2, 10, 3, 5, 6, 11, -1, -1, -1} }, /* 74 6 */
	{ { 5, 6, 11, 1, 2, 9, 9, 2, 10, 9, 10, 8} }, /* 75 12 */
	{ { 6, 10, 3, 6, 3, 5, 5, 3, 1, -1, -1, -1} }, /* 76 5 */
	{ { 0, 10, 8, 0, 5, 10, 0, 1, 5, 5, 6, 10} }, /* 77 14 */
	{ { 3, 6, 10, 0, 6, 3, 0, 5, 6, 0, 9, 5} }, /* 78 9 */
	{ { 6, 9, 5, 6, 10, 9, 10, 8, 9, -1, -1, -1} }, /* 79 5 */
	{ { 5, 6, 11, 4, 8, 7, -1, -1, -1, -1, -1, -1} }, /* 80 3 */
	{ { 4, 0, 3, 4, 3, 7, 6, 11, 5, -1, -1, -1} }, /* 81 7 */
	{ { 1, 0, 9, 5, 6, 11, 8, 7, 4, -1, -1, -1} }, /* 82 6 */
	{ { 11, 5, 6, 1, 7, 9, 1, 3, 7, 7, 4, 9} }, /* 83 12 */
	{ { 6, 2, 1, 6, 1, 5, 4, 8, 7, -1, -1, -1} }, /* 84 7 */
	{ { 1, 5, 2, 5, 6, 2, 3, 4, 0, 3, 7, 4} }, /* 85 10 */
	{ { 8, 7, 4, 9, 5, 0, 0, 5, 6, 0, 6, 2} }, /* 86 12 */
	{ { 4, 5, 9, 2, 6, 7, 2, 3, 7, -1, -1, -1} }, /* 87 7 */
	{ { 3, 2, 10, 7, 4, 8, 11, 5, 6, -1, -1, -1} }, /* 88 6 */
	{ { 5, 6, 11, 4, 2, 7, 4, 0, 2, 2, 10, 7} }, /* 89 12 */
	{ { 0, 9, 1, 4, 8, 7, 2, 10, 3, 5, 6, 11} }, /* 90 13 */
	{ { 1, 2, 11, 4, 5, 9, 6, 7, 10, -1, -1, -1} }, /* 91 6 */  //91 8?
	{ { 8, 7, 4, 3, 5, 10, 3, 1, 5, 5, 6, 10} }, /* 92 12 */
	{ { 6, 7, 10, 0, 1, 4, 1, 4, 5, -1, -1, -1} }, /* 93 7 */
	{ { 0, 3, 8, 4, 5, 9, 6, 7, 10, -1, -1, -1 } }, /* 94 6 */  //94 8?
	{ { 6, 9, 5, 6, 10, 9, 4, 9, 7, 7, 9, 10} }, /* 95 3 */
	{ { 11, 9, 4, 6, 11, 4, -1, -1, -1, -1, -1, -1} }, /* 96 2 */
	{ { 4, 6, 11, 4, 11, 9, 0, 3, 8, -1, -1, -1} }, /* 97 7 */
	{ { 11, 1, 0, 11, 0, 6, 6, 0, 4, -1, -1, -1} }, /* 98 5 */
	{ { 8, 1, 3, 8, 6, 1, 8, 4, 6, 6, 11, 1} }, /* 99 14 */
	{ { 1, 9, 4, 1, 4, 2, 2, 4, 6, -1, -1, -1} }, /* 100 5 */
	{ { 3, 8, 0, 1, 9, 2, 2, 9, 4, 2, 4, 6} }, /* 101 12 */
	{ { 0, 4, 2, 4, 6, 2, -1, -1, -1, -1, -1, -1} }, /* 102 8 */
	{ { 8, 2, 3, 8, 4, 2, 4, 6, 2, -1, -1, -1} }, /* 103 5 */
	{ { 11, 9, 4, 11, 4, 6, 10, 3, 2, -1, -1, -1} }, /* 104 7 */
	{ { 0, 2, 8, 2, 10, 8, 4, 11, 9, 4, 6, 11} }, /* 105 10 */
	{ { 3, 2, 10, 0, 6, 1, 0, 4, 6, 6, 11, 1} }, /* 106 12 */
	{ { 1, 2, 11, 4, 8, 10, 4, 6, 10, -1, -1, -1 } }, /* 107 7 */
	{ { 9, 4, 6, 9, 6, 3, 9, 3, 1, 10, 3, 6} }, /* 108 11 */
	{ { 0, 1, 9, 4, 8, 10, 4, 6, 10, -1, -1, -1 } }, /* 109 7 */
	{ { 3, 6, 10, 3, 0, 6, 0, 4, 6, -1, -1, -1} }, /* 110 5 */
	{ { 6, 8, 4, 10, 8, 6, -1, -1, -1, -1, -1, -1} }, /* 111 2 */
	{ { 7, 6, 11, 7, 11, 8, 8, 11, 9, -1, -1, -1} }, /* 112 5 */
	{ { 0, 3, 7, 0, 7, 11, 0, 11, 9, 6, 11, 7} }, /* 113 11 */
	{ { 11, 7, 6, 1, 7, 11, 1, 8, 7, 1, 0, 8} }, /* 114 9 */
	{ { 11, 7, 6, 11, 1, 7, 1, 3, 7, -1, -1, -1} }, /* 115 5 */
	{ { 1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6} }, /* 116 14 */
	{ { 0, 1, 9, 2, 3, 7, 2, 6, 7, -1, -1, -1 } }, /* 117 7 */
	{ { 7, 0, 8, 7, 6, 0, 6, 2, 0, -1, -1, -1} }, /* 118 5 */
	{ { 7, 2, 3, 6, 2, 7, -1, -1, -1, -1, -1, -1} }, /* 119 2 */
	{ { 2, 10, 3, 11, 8, 6, 11, 9, 8, 8, 7, 6} }, /* 120 12 */
	{ { 6, 7, 10, 0, 2, 11, 0, 9, 11, -1, -1, -1 } }, /* 121 7 */
	{ { 0, 3, 8, 1, 2, 11, 6, 7, 10, -1, -1, -1 } }, /* 122 6 */
	{ { 10, 1, 2, 10, 7, 1, 11, 1, 6, 6, 1, 7} }, /* 123 3 */
	{ { 6, 7, 10, 1, 3, 8, 1, 8, 9, -1, -1, -1 } }, /* 124 7 */
	{ { 0, 1, 9, 10, 7, 6, -1, -1, -1, -1, -1, -1} }, /* 125 4 */
	{ { 7, 0, 8, 7, 6, 0, 3, 0, 10, 10, 0, 6} }, /* 126 3 */
	{ { 7, 6, 10, -1, -1, -1, -1, -1, -1, -1, -1, -1} }, /* 127 1 */
};