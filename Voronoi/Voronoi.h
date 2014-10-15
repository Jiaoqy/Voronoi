#ifndef _VORONOI_H_
#define _VORONOI_H_

// INCLUDES ///////////////////////////////////////////////
#include "stdio.h"
#include "math.h"
#include <stdlib.h>
#include "glut.h"
#include <ctime>
#include <iostream>
#include "windows.h" // for time statistics

using namespace std;

// DEFINES ////////////////////////////////////////////////
#define MAX_VERTEX_NUM 4092

#ifdef SINGLE
#define REAL float
#else
#define REAL double
#endif


// TYPES //////////////////////////////////////////////////
typedef struct POINT2D_TYP
{
	REAL x;
	REAL y;

} POINT2D, *POINT2D_PTR;

typedef struct EDGE_TYP
{
	POINT2D center1;
	POINT2D center2;
	EDGE_TYP* pNext;

} EDGE, *EDGE_PTR;

typedef struct TRIANGLE_TYP
{
	int i1; // vertex index
	int i2;
	int i3;

	POINT2D center;
	TRIANGLE_TYP* pNext;
	TRIANGLE_TYP* pPrev;

} TRIANGLE, *TRIANGLE_PTR;

typedef struct MESH_TYP
{
	int vertex_num;
	int triangle_num;

	POINT2D_PTR pVerArr; // POINT2D to outer vertices arrary
	TRIANGLE_PTR pTriArr; // POINT2D to outer triangles arrary
	EDGE_PTR pEdgeArr; // POINT2D to voronoi edges arrary

} MESH, *MESH_PTR;


typedef struct VERTEX_TYP
{

	int value;
	int ID;
	POINT2D_PTR vernode;
	TRIANGLE_PTR pTriArr; // POINT2D to outer triangles arrary
	VERTEX_TYP* next;

} VERTEX, *VERTEX_PTR;
typedef struct ordertri_TYP
{

	int vertex_num;

	VERTEX_PTR VerArr; // POINT2D to outer vertices arrary

} orderVer, *ORDERVER_PTR;



// PROTOTYPES ///////////////////////////////////////////
// Delaunay triangulation functions
void InitMesh(MESH_PTR pMesh, int ver_num);
void UnInitMesh(MESH_PTR pMesh);

void AddBoundingBox(MESH_PTR pMesh);
void RemoveBoundingBox(MESH_PTR pMesh);;
void IncrementalDelaunay(MESH_PTR pMesh);

void Insert(MESH_PTR pMesh, int ver_index);
bool FlipTest(MESH_PTR pMesh, TRIANGLE_PTR pTestTri);

REAL InCircle(POINT2D_PTR pa, POINT2D_PTR pb, POINT2D_PTR pp, POINT2D_PTR  pd);
REAL InTriangle(MESH_PTR pMesh, POINT2D_PTR pVer, TRIANGLE_PTR pTri);

void InsertInTriangle(MESH_PTR pMesh, TRIANGLE_PTR pTargetTri, int ver_index);
void InsertOnEdge(MESH_PTR pMesh, TRIANGLE_PTR pTargetTri, int ver_index);

// Helper functions
void RemoveTriangleNode(MESH_PTR pMesh, TRIANGLE_PTR pTri);
TRIANGLE_PTR AddTriangleNode(MESH_PTR pMesh, TRIANGLE_PTR pPrevTri, int i1, int i2, int i3);

// I/O functions
void Input(char* pFile, MESH_PTR pMesh);
void CreatePionts(MESH_PTR pMesh);
void Output(char* pFile, MESH_PTR pMesh);
void printDelaunay(MESH_PTR pMesh);
void printVoronoi(MESH_PTR pMesh);

REAL CounterClockWise(POINT2D_PTR pa, POINT2D_PTR pb, POINT2D_PTR pc);

/////////////////////////////////////////////////////////////////////////////////////
void CalCenter(MESH_PTR pMesh);
void FindEdge(MESH_PTR pMesh);

#endif