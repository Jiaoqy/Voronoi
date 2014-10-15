#include "Voronoi.h"

void Input(char* pFile, MESH_PTR pMesh)
{
	FILE* fp = fopen(pFile, "r");

	if (!fp)
	{
		fprintf(stderr,"Error:%s open failed\n", pFile);
		exit(1);
	}

	//int face;
	int amount;

	//fscanf( fp, "%d", &face);
	fscanf( fp, "%d", &amount);
	if (amount < 3)
	{
		fprintf(stderr,"Error:vertex amount should be greater than 2, but it is %d \n",amount);
		exit(1);
	}

	InitMesh(pMesh, amount);

	REAL x,y,z;
	for ( int j=3; j<amount+3; ++j)
	{
		fscanf( fp, "%lg %lg %lg", &x, &y, &z);
		((POINT2D_PTR)(pMesh->pVerArr+j))->x = x;
		((POINT2D_PTR)(pMesh->pVerArr+j))->y = y;
		glBegin(GL_POINTS);
		glVertex2f(x, y);
		glEnd();
	}


	fclose(fp);
}

void CreatePionts(MESH_PTR pMesh){
	cout<<"请输入要插入点的数目： "<<endl;
	int num;
	cin>>num;
	InitMesh(pMesh, num);
	srand(time(0));
	for (int j=3; j<num+3; j++)
	{
		((POINT2D_PTR)(pMesh->pVerArr+j))->x = (double)(rand()%400) + (double)(rand()%10)/10.0+200;
		((POINT2D_PTR)(pMesh->pVerArr+j))->y = (double)(rand()%200) + (double)(rand()%10)/10.0+200;
		glBegin(GL_POINTS);
		glVertex2f(((POINT2D_PTR)(pMesh->pVerArr+j))->x, ((POINT2D_PTR)(pMesh->pVerArr+j))->y);
		glEnd();
	}
}

// Algorithm IncrementalDelaunay(V)
// Input: 由n个点组成的二维点集V
// Output: Delaunay三角剖分DT
//	1.add a appropriate triangle boudingbox to contain V ( such as: we can use triangle abc, a=(0, 3M), b=(-3M,-3M), c=(3M, 0), M=Max({|x1|,|x2|,|x3|,...} U {|y1|,|y2|,|y3|,...}))
//	2.initialize DT(a,b,c) as triangle abc
//	3.for i <- 1 to n 
//		do (Insert(DT(a,b,c,v1,v2,...,vi-1), vi))   
//	4.remove the boundingbox and relative triangle which cotains any vertex of triangle abc from DT(a,b,c,v1,v2,...,vn) and return DT(v1,v2,...,vn).
void IncrementalDelaunay(MESH_PTR pMesh)
{
	// Add a appropriate triangle boudingbox to contain V
	AddBoundingBox(pMesh);

	// Get a vertex/POINT2D vi from V and Insert(vi)
	for (int i=3; i<pMesh->vertex_num+3; i++)
	{
		Insert(pMesh, i);
	}

	// Remove the bounding box
	RemoveBoundingBox(pMesh);
}

//求三角网每个三角形的的外心
void CalCenter(MESH_PTR pMesh)
{
	int i1,i2,i3;
	REAL centerx,centery;
	POINT2D_PTR p1 = (POINT2D_PTR)malloc(sizeof(POINT2D));
	POINT2D_PTR p2 =(POINT2D_PTR)malloc(sizeof(POINT2D));
	POINT2D_PTR p3 = (POINT2D_PTR)malloc(sizeof(POINT2D));
	TRIANGLE_PTR pTri = (TRIANGLE_PTR)malloc(sizeof(TRIANGLE));
	pTri =  pMesh->pTriArr;
	while(pTri!=NULL)
	{
		i1=pTri->i1;
		i2=pTri->i2;
		i3=pTri->i3;
		POINT2D_PTR pVer1 = (POINT2D_PTR)(pMesh->pVerArr+i1);
		POINT2D_PTR pVer2 = (POINT2D_PTR)(pMesh->pVerArr+i2);
		POINT2D_PTR pVer3 = (POINT2D_PTR)(pMesh->pVerArr+i3);

		REAL   Xmove=pVer1->x;
		REAL   Ymove=pVer1->y;
		p2->x=pVer2->x-pVer1->x;
		p2->y=pVer2->y-pVer1->y;
		p3->x=pVer3->x-pVer1->x;
		p3->y=pVer3->y-pVer1->y;
		p1->x=0;
		p1->y=0;


		REAL   x1=p2->x,y1=p2->y,x2=p3->x,y2=p3->y;
		REAL   m=2.0*(x1*y2-y1*x2);
		//pTri->center.x=(x1*x1*y2-x2*x2*y1+y1*y2*(y1-y2))/m;
		//pTri->center.y=(x1*x2*(x2-x1)-y1*y1*x2+x1*y2*y2)/m;
		centerx=(x1*x1*y2-x2*x2*y1+y1*y2*(y1-y2))/m;
		centery=(x1*x2*(x2-x1)-y1*y1*x2+x1*y2*y2)/m;
		//radius=distance(center,p[0]);
		centerx+=Xmove;
		centery+=Ymove;
		pTri->center.x=centerx;
		pTri->center.y=centery;

		pTri=pTri->pNext;


		//        printf("%lg   %lg",centerx,centery);
	}

}

//寻找每条维诺边，把每条维诺边存入pMesh结构体的pEdgeArr成员变量中
void FindEdge(MESH_PTR pMesh)
{
	TRIANGLE_PTR pTri,pFor;
	POINT2D_PTR pVer;
	EDGE_PTR pEdge = pMesh->pEdgeArr;
	EDGE_PTR pNewEdge;
	int Va,Vb,Vp,*p,value =0;
	REAL Xm,Ym,k,k2,X0,Y0,t =0;
	pTri = pMesh->pTriArr;
	//遍历三角链表的每个三角形，分别找出每个三角形的相邻三角形，并把相邻三角形的外心连接形成维诺边
	while(pTri!=NULL)
	{
		p = &(pTri->i1);
		for(int j=0;j<3;j++)
		{
			if(j==0)
			{
				Va = *p++;Vb = *p;Vp = pTri->i3;
			}else if(j==1)
			{
				Va = *p++;Vb = *p;Vp = pTri->i1;

			}else
			{
				Va = *p;Vb = pTri->i1;Vp = pTri->i2;
			}
			pFor=pMesh->pTriArr;
			//在三角链表中寻找与其有共边VaVb的相邻三角形
			while(pFor!=NULL)
			{
				if((pFor->i1 == Va||pFor->i2 == Va||pFor->i3 == Va)&&(pFor->i1 == Vb||pFor->i2== Vb||pFor->i3 == Vb))//找到
				{
					pNewEdge = (EDGE_PTR)malloc(sizeof(EDGE));
					pNewEdge->center1.x = pTri->center.x;
					pNewEdge->center1.y = pTri->center.y;
					pNewEdge->center2.x = pFor->center.x;
					pNewEdge->center2.y = pFor->center.y;
					pEdge->pNext = pNewEdge;
					pEdge = pNewEdge;
					pEdge->pNext =NULL;
					if(pFor!=pTri)
						value++;

				}

				pFor = pFor->pNext;

			}
			if(value==0)//如果在三角链表中找不到共边VaVb的相邻三角形，则该边是三角网的外边，应向外作中垂线射线
			{
				//Xm和Ym是边VaVb的中点坐标

				Xm = ((pMesh->pVerArr+Va)->x+(pMesh->pVerArr+Vb)->x)/2;
				Ym =((pMesh->pVerArr+Va)->y+(pMesh->pVerArr+Vb)->y)/2;

				k = (pTri->center.y -Ym)/(pTri->center.x -Xm);

				pNewEdge = (EDGE_PTR)malloc(sizeof(EDGE));
				pNewEdge->center1.x = pTri->center.x;
				pNewEdge->center1.y = pTri->center.y;
				REAL yy =(pMesh->pVerArr+Va)->y - (pMesh->pVerArr+Vb)->y;
				REAL xx = (pMesh->pVerArr+Va)->x - (pMesh->pVerArr+Vb)->x;
				//如果边VaVb斜率不为0则k2是它的斜率
				if(yy&&xx)
				{k2 = (yy)/(xx);}
				else
				{k2 = 0;}
				X0 = (pMesh->pVerArr+Vp)->x;
				Y0 = X0*k2-k2*((pMesh->pVerArr+Va)->x)+(pMesh->pVerArr+Va)->y;
				pVer = &pTri->center;
				if(k2!=0)
				{
					if(k>0) //
					{
						if((pMesh->pVerArr+Vp)->y>Y0)
							t =10;
						else
							t=1000;
					}
					else
					{
						if((pMesh->pVerArr+Vp)->y>Y0)
							t =1000;
						else
							t=10;

					}
				}
				else
				{
					if(xx == 0)
					{
						if((pMesh->pVerArr+Vp)->x>Xm)
							t=10;
						else
							t=1000;
					}
					else
					{
						if((pMesh->pVerArr+Vp)->y>Ym)
							t=10;
						else
							t=1000;
					}

				}
				pNewEdge->center2.x = t;
				pNewEdge->center2.y = t*k-k*Xm+Ym;

				pEdge->pNext = pNewEdge;
				pEdge = pNewEdge;
				pEdge->pNext =NULL;

			}

			value = 0;

		}

		pTri = pTri->pNext;
	}

}

// The format of output file should be as follows:
// triangle index
// x1 y1 (the coordinate of first vertex of triangle)
// x2 y2 (the coordinate of second vertex of triangle)
// x3 y3 (the coordinate of third vertex of triangle)
void Output(char* pFile, MESH_PTR pMesh)
{
	FILE* fp = fopen(pFile, "w");
	if (!fp)
	{
		fprintf(stderr,"Error:%s open failed\n", pFile);

		UnInitMesh(pMesh);
		exit(1);
	}

	TRIANGLE_PTR pTri = pMesh->pTriArr;
	int* pi;
	int vertex_index;
	int tri_index = 0;
	while(pTri != NULL)	
	{
		fprintf(fp, "Triangle: %d\n", ++tri_index);

		pi = &(pTri->i1);
		for (int j=0; j<3; j++)	
		{	
			vertex_index = *pi++;		
			fprintf(fp, "%lg %lg\n", ((POINT2D_PTR)(pMesh->pVerArr+vertex_index))->x, ((POINT2D_PTR)(pMesh->pVerArr+vertex_index))->y);
		}

		pTri = pTri->pNext;
	}

	fclose(fp);

	UnInitMesh(pMesh);
}

// Allocate memory to store vertices and triangles
void InitMesh(MESH_PTR pMesh, int ver_num )
{
	// Allocate memory for vertex array
	pMesh->pVerArr = (POINT2D_PTR)malloc((ver_num+3)*sizeof(POINT2D));
	if (pMesh->pVerArr == NULL)
	{
		fprintf(stderr,"Error:Allocate memory for mesh failed\n");
		exit(1);
	}

	pMesh->vertex_num = ver_num;
	pMesh->pEdgeArr = (EDGE_PTR)malloc(sizeof(EDGE));
	pMesh->pEdgeArr->pNext =NULL;

}

// Deallocate memory
void UnInitMesh(MESH_PTR pMesh)
{
	// free vertices
	if(pMesh->pVerArr != NULL)
		free(pMesh->pVerArr);

	// free triangles
	TRIANGLE_PTR pTri = pMesh->pTriArr;
	TRIANGLE_PTR pTemp = NULL;
	while (pTri != NULL)
	{
		pTemp = pTri->pNext;
		free(pTri);
		pTri = pTemp;
	}
}

void AddBoundingBox(MESH_PTR pMesh)
{
	REAL max = 0;
	REAL max_x = 0;
	REAL max_y = 0;
	REAL t;

	for (int i=3; i<pMesh->vertex_num+3; i++)
	{
		t = abs(((POINT2D_PTR)(pMesh->pVerArr+i))->x);
		if (max_x < t)
		{
			max_x = t;
		}

		t = abs(((POINT2D_PTR)(pMesh->pVerArr+i))->y);
		if (max_y < t)
		{
			max_y = t;
		}
	}

	max = max_x > max_y ? max_x:max_y;


	POINT2D v1 = {0, 4*max};
	POINT2D v2 = {-4*max, -4*max};
	POINT2D v3 = {4*max, 0};

	// Assign to Vertex array
	*(pMesh->pVerArr) = v1;
	*(pMesh->pVerArr + 1) = v2;
	*(pMesh->pVerArr + 2) = v3;

	// add the Triangle boundingbox
	AddTriangleNode(pMesh, NULL, 0, 1, 2);
}

void RemoveBoundingBox(MESH_PTR pMesh)
{
	int statify[3]={0,0,0};
	int vertex_index;
	int* pi;
	//    int k = 1;

	// Remove the first triangle-boundingbox
	//pMesh->pTriArr = pMesh->pTriArr->pNext;
	//pMesh->pTriArr->pPrev = NULL; // as head

	TRIANGLE_PTR pTri = pMesh->pTriArr;
	TRIANGLE_PTR pNext = NULL;
	while (pTri != NULL)
	{
		pNext = pTri->pNext;

		statify[0] = 0;
		statify[1] = 0;
		statify[2] = 0;

		pi = &(pTri->i1);
		for (int j=0, k = 1; j<3; j++, k*= 2)
		{
			vertex_index = *pi++;

			if(vertex_index == 0 || vertex_index == 1 || vertex_index == 2) // bounding box vertex
			{
				statify[j] = k;
			}
		}

		switch(statify[0] | statify[1] | statify[2] )
		{
		case 0: // no statify
			break;
		case 1:
		case 2:
		case 4: // 1 statify, remove 1 triangle, 1 vertex
			RemoveTriangleNode(pMesh, pTri);
			break;
		case 3:
		case 5:
		case 6: // 2 statify, remove 1 triangle, 2 vertices
			RemoveTriangleNode(pMesh, pTri);
			break;
		case 7: // 3 statify, remove 1 triangle, 3 vertices
			RemoveTriangleNode(pMesh, pTri);
			break;
		default:
			break;
		}

		// go to next item
		pTri = pNext;
	}
}


// Return a positive value if the POINT2Ds pa, pb, and
// pc occur in counterclockwise order; a negative
// value if they occur in clockwise order; and zero
// if they are collinear. The result is also a rough
// approximation of twice the signed area of the
// triangle defined by the three POINT2Ds.
REAL CounterClockWise(POINT2D_PTR pa, POINT2D_PTR pb, POINT2D_PTR pc)
{
	return ((pb->x - pa->x)*(pc->y - pb->y) - (pc->x - pb->x)*(pb->y - pa->y));
}

// Adjust if the POINT2D lies in the triangle abc
REAL InTriangle(MESH_PTR pMesh, POINT2D_PTR pVer, TRIANGLE_PTR pTri)
{
	int vertex_index;
	POINT2D_PTR pV1, pV2, pV3;

	vertex_index =pTri->i1;
	pV1 = (POINT2D_PTR)(pMesh->pVerArr+vertex_index);
	vertex_index =pTri->i2;
	pV2 = (POINT2D_PTR)(pMesh->pVerArr+vertex_index);
	vertex_index =pTri->i3;
	pV3 = (POINT2D_PTR)(pMesh->pVerArr+vertex_index);

	REAL ccw1 = CounterClockWise(pV1, pV2, pVer);
	REAL ccw2 = CounterClockWise(pV2, pV3, pVer);
	REAL ccw3 = CounterClockWise(pV3, pV1, pVer);

	REAL r = -1;
	if (ccw1>0 && ccw2>0 && ccw3>0)
	{
		r = 1;
	}
	else if(ccw1*ccw2*ccw3 == 0 && (ccw1*ccw2 > 0 || ccw1*ccw3 > 0 || ccw2*ccw3 > 0) )
	{
		r = 0;
	}

	return r;
}

// Algorithm Insert(DT(a,b,c,v1,v2,...,vi-1), vi)
// 1.find the triangle vavbvc which contains vi // FindTriangle()
// 2.if (vi located at the interior of vavbvc)  
// 3.    then add triangle vavbvi, vbvcvi and vcvavi into DT // UpdateDT()
// FlipTest(DT, va, vb, vi)
// FlipTest(DT, vb, vc, vi)
// FlipTest(DT, vc, va, vi)
// 4.else if (vi located at one edge (E.g. edge vavb) of vavbvc) 
// 5.    then add triangle vavivc, vivbvc, vavdvi and vivdvb into DT (here, d is the third vertex of triangle which contains edge vavb) // UpdateDT()
// FlipTest(DT, va, vd, vi)
// FlipTest(DT, vc, va, vi)
// FlipTest(DT, vd, vb, vi)
// FlipTest(DT, vb, vc, vi)
// 6.return DT(a,b,c,v1,v2,...,vi)
void Insert(MESH_PTR pMesh, int ver_index)
{
	POINT2D_PTR pVer = (POINT2D_PTR)(pMesh->pVerArr+ver_index);
	TRIANGLE_PTR pTargetTri = NULL;
	TRIANGLE_PTR pEqualTri1 = NULL;
	TRIANGLE_PTR pEqualTri2 = NULL;

	int j = 0;
	TRIANGLE_PTR pTri = pMesh->pTriArr;
	while (pTri != NULL)
	{
		REAL r = InTriangle(pMesh, pVer, pTri);
		if(r > 0) // should be in triangle
		{
			pTargetTri = pTri;
		}
		else if (r == 0) // should be on edge
		{
			if(j == 0)
			{
				pEqualTri1 = pTri;
				j++;
			}
			else
			{
				pEqualTri2 = pTri;
			}

		}

		pTri = pTri->pNext;
	}

	if (pEqualTri1 != NULL && pEqualTri2 != NULL)
	{
		InsertOnEdge(pMesh, pEqualTri1, ver_index);
		InsertOnEdge(pMesh, pEqualTri2, ver_index);
	}
	else
	{
		InsertInTriangle(pMesh, pTargetTri, ver_index);
	}
}

void InsertInTriangle(MESH_PTR pMesh, TRIANGLE_PTR pTargetTri, int ver_index)
{
	int index_a, index_b, index_c;
	TRIANGLE_PTR pTri = NULL;
	TRIANGLE_PTR pNewTri = NULL;

	pTri = pTargetTri;
	if(pTri == NULL)
	{
		return;
	}

	// Inset p into target triangle
	index_a = pTri->i1;
	index_b = pTri->i2;
	index_c = pTri->i3;

	// Insert edge pa, pb, pc
	for(int i=0; i<3; i++)
	{
		// allocate memory
		if(i == 0)
		{
			pNewTri = AddTriangleNode(pMesh, pTri, index_a, index_b, ver_index);
		}
		else if(i == 1)
		{
			pNewTri = AddTriangleNode(pMesh, pTri, index_b, index_c, ver_index);
		}
		else
		{
			pNewTri = AddTriangleNode(pMesh, pTri, index_c, index_a, ver_index);
		}

		// go to next item
		if (pNewTri != NULL)
		{
			pTri = pNewTri;
		}
		else
		{
			pTri = pTri;
		}
	}

	// Get the three sub-triangles
	pTri = pTargetTri;
	TRIANGLE_PTR pTestTri[3];
	for (int i=0; i< 3; i++)
	{
		pTestTri[i] = pTri->pNext;

		pTri = pTri->pNext;
	}

	// remove the Target Triangle
	RemoveTriangleNode(pMesh, pTargetTri);

	for (int i=0; i< 3; i++)
	{
		// Flip test
		FlipTest(pMesh, pTestTri[i]);
	}
}

void InsertOnEdge(MESH_PTR pMesh, TRIANGLE_PTR pTargetTri, int ver_index)
{
	int index_a, index_b, index_c;
	TRIANGLE_PTR pTri = NULL;
	TRIANGLE_PTR pNewTri = NULL;

	pTri = pTargetTri;
	if(pTri == NULL)
	{
		return;
	}

	// Inset p into target triangle
	index_a = pTri->i1;
	index_b = pTri->i2;
	index_c = pTri->i3;

	// Insert edge pa, pb, pc
	for(int i=0; i<3; i++)
	{
		// allocate memory
		if(i == 0)
		{
			pNewTri = AddTriangleNode(pMesh, pTri, index_a, index_b, ver_index);
		}
		else if(i == 1)
		{
			pNewTri = AddTriangleNode(pMesh, pTri, index_b, index_c, ver_index);
		}
		else
		{
			pNewTri = AddTriangleNode(pMesh, pTri, index_c, index_a, ver_index);
		}

		// go to next item
		if (pNewTri != NULL)
		{
			pTri = pNewTri;
		}
		else
		{
			pTri = pTri;
		}
	}

	// Get the two sub-triangles
	pTri = pTargetTri;
	TRIANGLE_PTR pTestTri[2];
	for (int i=0; i< 2; i++)
	{
		pTestTri[i] = pTri->pNext;
		pTri = pTri->pNext;
	}

	// remove the Target Triangle
	RemoveTriangleNode(pMesh, pTargetTri);

	for (int i=0; i< 2; i++)
	{
		// Flip test
		FlipTest(pMesh, pTestTri[i]);
	}
}

// Precondition: the triangle satisfies CCW order
// Algorithm FlipTest(DT(a,b,c,v1,v2,...,vi), va, vb, vi)
// 1.find the third vertex (vd) of triangle which contains edge vavb // FindThirdVertex()
// 2.if(vi is in circumcircle of abd)  // InCircle()
// 3.    then remove edge vavb, add new edge vivd into DT // UpdateDT()
//		  FlipTest(DT, va, vd, vi)
//		  FlipTest(DT, vd, vb, vi)

bool FlipTest(MESH_PTR pMesh, TRIANGLE_PTR pTestTri)
{
	bool flipped = false;

	int index_a = pTestTri->i1;
	int index_b = pTestTri->i2;
	int index_p = pTestTri->i3;

	int statify[3]={0,0,0};
	int vertex_index;
	int* pi;
	//    int k = 1;

	// find the triangle which has edge consists of start and end
	TRIANGLE_PTR pTri = pMesh->pTriArr;

	int index_d = -1;
	while (pTri != NULL)
	{
		statify[0] = 0;
		statify[1] = 0;
		statify[2] = 0;

		pi = &(pTri->i1);
		for (int j=0, k = 1; j<3; j++, k*= 2)
		{
			vertex_index = *pi++;
			if(vertex_index == index_a || vertex_index == index_b)
			{
				statify[j] = k;
			}
		}

		switch(statify[0] | statify[1] | statify[2] )
		{
		case 3:
			if(CounterClockWise((POINT2D_PTR)(pMesh->pVerArr+index_a), (POINT2D_PTR)(pMesh->pVerArr+index_b), (POINT2D_PTR)(pMesh->pVerArr+pTri->i3)) < 0)
			{
				index_d = pTri->i3;
			}

			break;
		case 5:
			if(CounterClockWise((POINT2D_PTR)(pMesh->pVerArr+index_a), (POINT2D_PTR)(pMesh->pVerArr+index_b), (POINT2D_PTR)(pMesh->pVerArr+pTri->i2)) < 0)
			{
				index_d = pTri->i2;
			}

			break;
		case 6:
			if(CounterClockWise((POINT2D_PTR)(pMesh->pVerArr+index_a), (POINT2D_PTR)(pMesh->pVerArr+index_b), (POINT2D_PTR)(pMesh->pVerArr+pTri->i1)) < 0)
			{
				index_d = pTri->i1;
			}

			break;

		default:
			break;
		}

		if (index_d != -1)
		{
			POINT2D_PTR pa = (POINT2D_PTR)(pMesh->pVerArr+index_a);
			POINT2D_PTR pb = (POINT2D_PTR)(pMesh->pVerArr+index_b);
			POINT2D_PTR pd = (POINT2D_PTR)(pMesh->pVerArr+index_d);
			POINT2D_PTR pp = (POINT2D_PTR)(pMesh->pVerArr+index_p);

			if(InCircle( pa, pb, pp, pd) < 0) // not local Delaunay
			{
				flipped = true;

				// add new triangle adp,  dbp, remove abp, abd.
				// allocate memory for adp
				TRIANGLE_PTR pT1 = AddTriangleNode(pMesh, pTestTri, pTestTri->i1, index_d, pTestTri->i3);
				// allocate memory for dbp
				TRIANGLE_PTR pT2 = AddTriangleNode(pMesh, pT1, index_d, pTestTri->i2, index_p);
				// remove abp
				RemoveTriangleNode(pMesh, pTestTri);
				// remove abd
				RemoveTriangleNode(pMesh, pTri);

				FlipTest(pMesh, pT1); // pNewTestTri satisfies CCW order
				FlipTest(pMesh, pT2); // pNewTestTri2  satisfies CCW order

				break;
			}
		}

		// go to next item
		pTri = pTri->pNext;
	}

	return flipped;
}

// In circle test, use vector cross product
REAL InCircle(POINT2D_PTR pa, POINT2D_PTR pb, POINT2D_PTR pp, POINT2D_PTR  pd)
{
	REAL det;
	REAL alift, blift, plift, bdxpdy, pdxbdy, pdxady, adxpdy, adxbdy, bdxady;

	REAL adx = pa->x - pd->x;
	REAL ady = pa->y - pd->y;

	REAL bdx = pb->x - pd->x;
	REAL bdy = pb->y - pd->y;

	REAL pdx = pp->x - pd->x;
	REAL pdy = pp->y - pd->y;

	bdxpdy = bdx * pdy;
	pdxbdy = pdx * bdy;
	alift = adx * adx + ady * ady;

	pdxady = pdx * ady;
	adxpdy = adx * pdy;
	blift = bdx * bdx + bdy * bdy;

	adxbdy = adx * bdy;
	bdxady = bdx * ady;
	plift = pdx * pdx + pdy * pdy;

	det = alift * (bdxpdy - pdxbdy)
		+ blift * (pdxady - adxpdy)
		+ plift * (adxbdy - bdxady);

	return -det;
}

// Remove a node from the triangle list and deallocate the memory
void RemoveTriangleNode(MESH_PTR pMesh, TRIANGLE_PTR pTri)
{
	if (pTri == NULL)
	{
		return;
	}

	// remove from the triangle list
	if (pTri->pPrev != NULL)
	{
		pTri->pPrev->pNext = pTri->pNext;
	}
	else // remove the head, need to reset the root node
	{
		pMesh->pTriArr = pTri->pNext;
	}

	if (pTri->pNext != NULL)
	{
		pTri->pNext->pPrev = pTri->pPrev;
	}

	// deallocate memory
	free(pTri);
}

// Create a new node and add it into triangle list
TRIANGLE_PTR AddTriangleNode(MESH_PTR pMesh, TRIANGLE_PTR pPrevTri, int i1, int i2, int i3)
{
	// test if 3 vertices are co-linear
	if(CounterClockWise((POINT2D_PTR)(pMesh->pVerArr+i1), (POINT2D_PTR)(pMesh->pVerArr+i2), (POINT2D_PTR)(pMesh->pVerArr+i3)) == 0)
	{
		return NULL;
	}

	// allocate memory
	TRIANGLE_PTR pNewTestTri = (TRIANGLE_PTR)malloc(sizeof(TRIANGLE));

	pNewTestTri->i1 = i1;
	pNewTestTri->i2 = i2;
	pNewTestTri->i3 = i3;

	// insert after prev triangle
	if (pPrevTri == NULL) // add root
	{
		pMesh->pTriArr = pNewTestTri;
		pNewTestTri->pNext = NULL;
		pNewTestTri->pPrev = NULL;
	}
	else
	{
		pNewTestTri->pNext = pPrevTri->pNext;
		pNewTestTri->pPrev = pPrevTri;

		if(pPrevTri->pNext != NULL)
		{
			pPrevTri->pNext->pPrev = pNewTestTri;
		}

		pPrevTri->pNext = pNewTestTri;
	}

	return pNewTestTri;
}

void printDelaunay(MESH_PTR pMesh){
	//MESH mesh;
	//Input("D://input_POINT2Ds.txt", &mesh);
	// 	CreatePionts(pMesh);
	// 	IncrementalDelaunay(pMesh);

	//MESH_PTR pMesh = &mesh;
	TRIANGLE_PTR pTri = pMesh->pTriArr;
	int* pi;
	int vertex_index;
	int tri_index = 0;
	REAL vex[3][2];
	while (pTri != NULL)
	{
		pi = &(pTri->i1);
		for (int j = 0; j<3; j++)
		{
			vertex_index = *pi++;
			vex[j][0] = ((POINT2D_PTR)(pMesh->pVerArr+vertex_index))->x;
			vex[j][1] = ((POINT2D_PTR)(pMesh->pVerArr+vertex_index))->y;
		}
		// 		REAL circumcenter[2];
		// 		REAL temp = 2*(vex[0][0]*(vex[1][1]-vex[2][1])+vex[1][0]*(vex[2][1]-vex[0][1])+vex[2][0]*(vex[0][1]-vex[1][1]));
		// 		circumcenter[0] = ((vex[0][0]*vex[0][0]*(vex[1][1]-vex[2][1])+vex[1][0]*vex[1][0]*(vex[2][1]-vex[0][1])+vex[2][0]*vex[2][0]*(vex[0][1]-vex[1][1])) - ((vex[0][1]-vex[1][1])*(vex[1][1]-vex[2][1])*(vex[2][1]-vex[0][1]))) / temp;
		// 		circumcenter[1] = (-(vex[0][1]*vex[0][1]*(vex[1][0]-vex[2][0])+vex[1][1]*vex[1][1]*(vex[2][0]-vex[0][0])+vex[2][1]*vex[2][1]*(vex[0][0]-vex[1][0])) + ((vex[0][0]-vex[1][0])*(vex[1][0]-vex[2][0])*(vex[2][0]-vex[0][0]))) / temp;
		// 
		// 		glBegin(GL_POINT2DS);
		// 		  glVertex2dv(circumcenter);
		// 		glEnd();

		glBegin(GL_LINES);
		glVertex2dv(vex[0]);
		glVertex2dv(vex[1]);
		glVertex2dv(vex[0]);
		glVertex2dv(vex[2]);
		glVertex2dv(vex[2]);
		glVertex2dv(vex[1]);
		glEnd();
		pTri = pTri->pNext;
	}
	UnInitMesh(pMesh);
}

void printVoronoi(MESH_PTR pMesh){

	TRIANGLE_PTR pTri = pMesh->pTriArr;
	int* pi;
	int vertex_index;
	int tri_index = 0;
	REAL vex[3][2];
	while (pTri != NULL)
	{
		pi = &(pTri->i1);
		for (int j = 0; j<3; j++)
		{
			vertex_index = *pi++;
			vex[j][0] = ((POINT2D_PTR)(pMesh->pVerArr+vertex_index))->x;
			vex[j][1] = ((POINT2D_PTR)(pMesh->pVerArr+vertex_index))->y;
		}

		glLineWidth(1.0);
		glColor3f(0.0, 1.0, 0.0);
		glBegin(GL_LINES);
		glVertex2dv(vex[0]);
		glVertex2dv(vex[1]);
		glVertex2dv(vex[0]);
		glVertex2dv(vex[2]);
		glVertex2dv(vex[2]);
		glVertex2dv(vex[1]);
		glEnd();
		pTri = pTri->pNext;
	}

	CalCenter(pMesh);
	FindEdge(pMesh);

	EDGE_PTR pEdge = pMesh->pEdgeArr->pNext;

	glLineWidth(2.0);
	glColor3f(1.0, 0.0, 0.0);

	//画出维诺边
	while(pEdge!=NULL)
	{
		glBegin(GL_LINES);
		glVertex2d(pEdge->center1.x,pEdge->center1.y);
		glVertex2d(pEdge->center2.x,pEdge->center2.y);
		glEnd();
		pEdge = pEdge->pNext;

	}

	UnInitMesh(pMesh);
}