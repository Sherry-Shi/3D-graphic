/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"

/* NOT part of API - just for general assistance */

short	ctoi(float color)		/* convert float color to GzIntensity short */
{
  return(short)((int)(color * ((1 << 12) - 1)));
}

int GzRotXMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along x axis
// Pass back the matrix using mat value
	float pi=3.14159265359;
	float sin1=sin(degree/180*pi);
	float cosin1=cos(degree/180*pi);
	mat[1][1]=cosin1;
	mat[1][2]=-sin1;
	mat[2][1]=sin1;
	mat[2][2]=cosin1;
	return GZ_SUCCESS;
}


int GzRotYMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along y axis
// Pass back the matrix using mat value
	float pi=3.14159265359;
	float sin1=sin(degree/180*pi);
	float cosin1=cos(degree/180*pi);
	mat[0][0]=cosin1;
	mat[0][2]=sin1;
	mat[2][0]=-sin1;
	mat[2][2]=cosin1;
	return GZ_SUCCESS;
}


int GzRotZMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along z axis
// Pass back the matrix using mat value
	float pi=3.14159265359;
	float sin1=sin(degree/180*pi);
	float cosin1=cos(degree/180*pi);
	mat[0][0]=cosin1;
	mat[0][1]=-sin1;
	mat[1][0]=sin1;
	mat[1][1]=cosin1;

	return GZ_SUCCESS;
}


int GzTrxMat(GzCoord translate, GzMatrix mat)
{
// Create translation matrix
// Pass back the matrix using mat value
	mat[0][3]=translate[0];
	mat[1][3]=translate[1];
	mat[2][3]=translate[2];

	return GZ_SUCCESS;
}


int GzScaleMat(GzCoord scale, GzMatrix mat)
{
// Create scaling matrix
// Pass back the matrix using mat value
	mat[0][0]=scale[0];
	mat[1][1]=scale[1];
	mat[2][2]=scale[2];

	return GZ_SUCCESS;
}


//----------------------------------------------------------
// Begin main functions

int GzNewRender(GzRender **render, GzDisplay *display)
{
/*  
- malloc a renderer struct 
- setup Xsp and anything only done once 
- save the pointer to display 
- init default camera 
*/ 
	GzRender *p = new GzRender;
	memset(p, 0, sizeof(GzRender));
	*render = p;
	p->display = display;

	p->Xsp[0][0] = (display->xres)/2;
	p->Xsp[1][1] = -(display->yres)/2;
	p->Xsp[2][2] = 1;
	p->Xsp[3][3] = 1;
	p->Xsp[0][3] = (display->xres)/2;
	p->Xsp[1][3] = (display->yres)/2;

	p->matlevel = -1;

	p->camera.position[0] = DEFAULT_IM_X;
	p->camera.position[1] = DEFAULT_IM_Y;
	p->camera.position[2] = DEFAULT_IM_Z;
	p->camera.lookat[0] = 0;
	p->camera.lookat[1] = 0;
	p->camera.lookat[2] = 0;
	p->camera.worldup[0] = 0;
	p->camera.worldup[1] = 1;
	p->camera.worldup[2] = 0;
	p->camera.FOV = DEFAULT_FOV;

	return GZ_SUCCESS;
}


int GzFreeRender(GzRender *render)
{
/* 
-free all renderer resources
*/
	if(render)
		delete render;

	return GZ_SUCCESS;
}

int GzBeginRender(GzRender *render)
{
/*  
- setup for start of each frame - init frame buffer color,alpha,z
- compute Xiw and projection xform Xpi from camera definition 
- init Ximage - put Xsp at base of stack, push on Xpi and Xiw 
- now stack contains Xsw and app can push model Xforms when needed
*/ 
	if (render->display == NULL) return GZ_FAILURE;
	for (int i = 0; i < render->display->yres*render->display->xres; i++)
	{
		render->display->fbuf[i].red = 2500;
		render->display->fbuf[i].green = 970;
		render->display->fbuf[i].blue = 0;
		render->display->fbuf[i].alpha = 1;
		render->display->fbuf[i].z = 0;
	}
	int bytes = sizeof(int);
	int int_max = int(pow(2.0,(int)(bytes*8))/2-1);
	int i,j;
	for(j=0;j<render->display->yres;j++)
	{
		for(i=0;i<render->display->xres;i++)
		{
			render->display->fbuf[j*render->display->xres+i].z=int_max;
		}
	}
	float pi=3.14159265359;
	float d = 1/tan(render->camera.FOV/180*pi/2);
	render->Xsp[2][2] = int_max ;
	GzPushMatrix(render, render->Xsp);
	GzMatrix Xpi;
	memset(Xpi, 0, sizeof(GzMatrix));
	Xpi[0][0] = 1;
	Xpi[1][1] = 1;
	Xpi[2][2] = 1/d;
	Xpi[3][3] = 1;
	Xpi[3][2] = 1/d;
	GzPushMatrix(render, Xpi);
	GzCoord x, y, z;
	GzCoord C, l, up;
	memcpy(up, render->camera.worldup, sizeof(GzCoord));
	memcpy(l, render->camera.lookat, sizeof(GzCoord));
	memcpy(C, render->camera.position, sizeof(GzCoord));
	float cl = sqrt((l[0] - C[0])*(l[0] - C[0])+(l[1] - C[1])*(l[1] - C[1])+(l[2] - C[2])*(l[2] - C[2]));
	
	z[0] = (l[0] - C[0]) / cl;
	z[1] = (l[1] - C[1]) / cl;
	z[2] = (l[2] - C[2]) / cl;

	float y0 = up[0] - (up[0]*z[0]+up[1]*z[1]+up[2]*z[2])*z[0];
	float y1 = up[1] - (up[0]*z[0]+up[1]*z[1]+up[2]*z[2])*z[1];
	float y2 = up[2] - (up[0]*z[0]+up[1]*z[1]+up[2]*z[2])*z[2];
	float up1 = sqrt(y0*y0+y1*y1+y2*y2);
	y[0] =y0/up1;
	y[1] =y1/up1;
	y[2] =y2/up1;

	float x0 = y[1]*z[2] - y[2]*z[1];
	float x1 = y[2]*z[0] - y[0]*z[2];
	float x2 = y[0]*z[1] - y[1]*z[0];
	float xx = sqrt(x0*x0+x1*x1+x2*x2);
	x[0] =x0 /xx;
	x[1] =x1/ xx;
	x[2] =x2/ xx;
	GzMatrix Xiw;
	memset(Xiw, 0, sizeof(GzMatrix));
	memcpy(Xiw[0], x, sizeof(GzCoord));
	memcpy(Xiw[1], y, sizeof(GzCoord));
	memcpy(Xiw[2], z, sizeof(GzCoord));
	Xiw[0][3] = -(x[0]*C[0]+x[1]*C[1]+x[2]*C[2]);
	Xiw[1][3] = -(y[0]*C[0]+y[1]*C[1]+y[2]*C[2]);
	Xiw[2][3] = -(z[0]*C[0]+z[1]*C[1]+z[2]*C[2]);
	Xiw[3][3] = 1;
	GzPushMatrix(render, Xiw);
	
	return GZ_SUCCESS;
}

int GzPutCamera(GzRender *render, GzCamera *camera)
{
/*
- overwrite renderer camera structure with new camera definition
*/
	memcpy(render->camera.position, camera->position, sizeof(GzCoord));
	memcpy(render->camera.lookat, camera->lookat, sizeof(GzCoord));
	memcpy(render->camera.worldup, camera->worldup, sizeof(GzCoord));
	render->camera.FOV = camera->FOV;

	return GZ_SUCCESS;	
}

int GzPushMatrix(GzRender *render, GzMatrix	matrix)
{
/*
- push a matrix onto the Ximage stack
- check for stack overflow
*/
	if(render->matlevel >= MATLEVELS-1) return GZ_FAILURE;
	else{
		if(render->matlevel != -1)
		{
			int top = render->matlevel + 1;
			for(int k1=0; k1<4; k1++)
			{
				for(int k2=0; k2<4; k2++)
				{
					float i1=render->Ximage[top-1][k1][0]*matrix[0][k2];
					float i2=render->Ximage[top-1][k1][1]*matrix[1][k2];
					float i3=render->Ximage[top-1][k1][2]*matrix[2][k2];
					float i4=render->Ximage[top-1][k1][3]*matrix[3][k2];
					render->Ximage[top][k1][k2] =i1 +i2+i3+i4;
				}
			}
			render->matlevel = top;
		}
		else{
			memcpy(render->Ximage[0], matrix, sizeof(GzMatrix));
			render->matlevel = render->matlevel+1;
		}	
		return GZ_SUCCESS;
	}
}

int GzPopMatrix(GzRender *render)
{
/*
- pop a matrix off the Ximage stack
- check for stack underflow
*/
	if(render->matlevel <= -1) return GZ_FAILURE;
	else{
		int top =render->matlevel-1;
		render->matlevel=top;
		return GZ_SUCCESS;
	}
		
}


int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer	*valueList) /* void** valuelist */
{
/*
- set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
- later set shaders, interpolaters, texture maps, and lights
*/
	if(nameList[0]==GZ_RGB_COLOR)
	{
		render->flatcolor[0]=((float*)valueList[0])[0];
		render->flatcolor[1]=((float*)valueList[0])[1];
		render->flatcolor[2]=((float*)valueList[0])[2];
	}

	return GZ_SUCCESS;
}

int GzPutTriangle(GzRender	*render, int numParts, GzToken *nameList, GzPointer	*valueList)
/* numParts : how many names and values */
{
/*  
- pass in a triangle description with tokens and values corresponding to 
      GZ_POSITION:3 vert positions in model space 
- Xform positions of verts using matrix on top of stack 
- Clip - just discard any triangle with any vert(s) behind view plane 
       - optional: test for triangles with all three verts off-screen (trivial frustum cull)
- invoke triangle rasterizer  
*/ 
	float x1,x2,x3,y1,y2,y3,z1,z2,z3,temp;

	x1=((GzCoord*)valueList[0])[0][0];
	y1=((GzCoord*)valueList[0])[0][1];
	z1=((GzCoord*)valueList[0])[0][2];
	x2=((GzCoord*)valueList[0])[1][0];
	y2=((GzCoord*)valueList[0])[1][1];
	z2=((GzCoord*)valueList[0])[1][2];
	x3=((GzCoord*)valueList[0])[2][0];
	y3=((GzCoord*)valueList[0])[2][1];
	z3=((GzCoord*)valueList[0])[2][2];
	GzMatrix mat;
	memcpy(mat, render->Ximage[render->matlevel], sizeof(GzMatrix));
	float z11, z12, z13;

	z11 = mat[2][0]*x1+mat[2][1]*y1+mat[2][2]*z1+mat[2][3];
	z12 = mat[2][0]*x2+mat[2][1]*y2+mat[2][2]*z2+mat[2][3];
	z13 = mat[2][0]*x3+mat[2][1]*y3+mat[2][2]*z3+mat[2][3];	

	if(z11>0&&z12>0&&z13>0)//cull z<0
	{
		float x11, x12, x13, y11, y12, y13;
		float w1, w2, w3;
		w1 = mat[3][0]*x1+mat[3][1]*y1+mat[3][2]*z1+mat[3][3];
		x11 = (mat[0][0]*x1+mat[0][1]*y1+mat[0][2]*z1+mat[0][3])/w1;
		y11 = (mat[1][0]*x1+mat[1][1]*y1+mat[1][2]*z1+mat[1][3])/w1;
		z11 =z11/ w1;
		x1=x11;
		y1=y11;
		z1=z11;
		w2 = mat[3][0]*x2+mat[3][1]*y2+mat[3][2]*z2+mat[3][3];
		x12 = (mat[0][0]*x2+mat[0][1]*y2+mat[0][2]*z2+mat[0][3])/w2;
		y12 = (mat[1][0]*x2+mat[1][1]*y2+mat[1][2]*z2+mat[1][3])/w2;
		z12 =z12/w2;
		x2=x12;
		y2=y12;
		z2=z12;
		w3 = mat[3][0]*x3+mat[3][1]*y3+mat[3][2]*z3+mat[3][3];
		x13 = (mat[0][0]*x3+mat[0][1]*y3+mat[0][2]*z3+mat[0][3])/w3;
		y13 = (mat[1][0]*x3+mat[1][1]*y3+mat[1][2]*z3+mat[1][3])/w3;
		z13 =z13/w3;
		x3=x13;
		y3=y13;
		z3=z13;
	}

	if(y1<y2)
	{
		temp=y1;
		y1=y2;
		y2=temp;
		temp=x1;
		x1=x2;
		x2=temp;
		temp=z1;
		z1=z2;
		z2=temp;
	}
	if(y1<y3)
	{
		temp=y1;
		y1=y3;
		y3=temp;
		temp=x1;
		x1=x3;
		x3=temp;
		temp=z1;
		z1=z3;
		z3=temp;
	}
	if(y2<y3)
	{
		temp=y2;
		y2=y3;
		y3=temp;
		temp=x2;
		x2=x3;
		x3=temp;
		temp=z2;
		z2=z3;
		z3=temp;
	}
	else if(y2==y3)
	{
		if(x2>x3)
		{
			temp=x2;
			x2=x3;
			x3=temp;
			temp=y2;
			y2=y3;
			y3=temp;
			temp=z2;
			z2=z3;
			z3=temp;
		}
	}
	if(y1==y2)
	{
		if(x1>x2)
		{
			temp=y1;
			y1=y2;
			y2=temp;
			temp=x1;
			x1=x2;
			x2=temp;
			temp=z1;
			z1=z2;
			z2=temp;
		}
	}
	// Switch points over


	float xK12, xK23, xK13, xB12=0, xB13=0, xB23=0;
	int hori_edge = 0;
	int xScan, yScan, xL, xR;
	GzIntensity r, g, b, a;
	GzDepth zScan_new, zScan_old;
	float xM;

	if(y2==y3)
	{
		hori_edge=23;
	}
	if(y1==y3)
	{
		hori_edge=13;
	}
	if(y1==y2)
	{
		hori_edge=12;
	}
	
	if(hori_edge!=12)
	{
		xK12 = (x2-x1)/(y2-y1);
		xB12 = x1-y1*xK12;
	}
	if(hori_edge!=13)
	{
		xK13 = (x3-x1)/(y3-y1);
		xB13 = x1-y1*xK13;
	}
	if(hori_edge!=23)
	{
		xK23 = (x3-x2)/(y3-y2);
		xB23 = x2-y2*xK23;
	}

	// Calculate the parameters of plane's equation Ax+By+Cz+D=0 according to point 1 2 3
	float A, B, C, D;
	A = y1*z2 - y1*z3 - y2*z1 + y2*z3 + y3*z1 - y3*z2;
	B = -x1*z2 + x1*z3 + x2*z1 - x2*z3 - x3*z1 + x3*z2;
	C = x1*y2 - x1*y3 - x2*y1 + x2*y3 + x3*y1 - x3*y2;
	D = -x1*y2*z3 + x1*y3*z2 + x2*y1*z3 - x2*y3*z1 - x3*y1*z2 + x3*y2*z1;

	if(hori_edge!=0)		// Horizontal edge exists
	{
		if(hori_edge==12)		// Edge 12 is horizontal. Then the scan line should from edge 13 to edge 23.
		{
			yScan=ceil(y3);
			while(yScan<=y1)
			{
				xL=ceil(xK13*yScan+xB13);
				xR=floor(xK23*yScan+xB23);
				xScan=xL;
				while(xScan<=xR)
				{
					if(xScan>=0&&xScan<render->display->xres&&yScan>=0&&yScan<render->display->yres)
					{
						zScan_new=(GzDepth)((-A*xScan-B*yScan-D)/C);
						GzGetDisplay(render->display, xScan, yScan, &r, &g, &b, &a, &zScan_old);
						if(zScan_new<zScan_old)
							GzPutDisplay(render->display, xScan, yScan, ctoi(render->flatcolor[0]), ctoi(render->flatcolor[1]), ctoi(render->flatcolor[2]), 1, zScan_new);
					}
					xScan++;
				}
				yScan++;
			}
		}
		else if(hori_edge==23)	// Edge 23 is horizontal. Then the scan line should from edge 12 to edge 13
		{
			yScan = ceil(y3);
			while(yScan<=y1)
			{
				xL=ceil(xK12*yScan+xB12);
				xR=floor(xK13*yScan+xB13);
				xScan=xL;
				while(xScan<=xR)
				{
					if(xScan>=0&&xScan<render->display->xres&&yScan>=0&&yScan<render->display->yres)
					{
						zScan_new=(GzDepth)((-A*xScan-B*yScan-D)/C);
						GzGetDisplay(render->display, xScan, yScan, &r, &g, &b, &a, &zScan_old);
						if(zScan_new<zScan_old)
							GzPutDisplay(render->display, xScan, yScan, ctoi(render->flatcolor[0]), ctoi(render->flatcolor[1]), ctoi(render->flatcolor[2]), 1, zScan_new);
					}
					xScan++;
				}
				yScan++;
			}
		}
		// Edge 13 can not be horizontal, because 1 is the top-left point and 3 is the bottom-right point. If edge 13 is horizontal, the three points will be in a horizontal line
	}

	else if(hori_edge==0)	// No horizontal edges
	{
		xM = y2*xK13+xB13;
	 if (xM<x2)		// Middle point 2 is on the right side of the long edge. Then the left point of scan line should be on edge 13, the right point of the scan line should be on edge 12 then edge 23
	{
		yScan = ceil(y3);
		while (yScan <= y2)	// Scan line form edge 13 to edge 23
		{
			xL = ceil(xK13*yScan + xB13);
			xR = floor(xK23*yScan + xB23);
			xScan = xL;
			while (xScan <= xR)
			{
				if (xScan >= 0 && xScan<render->display->xres&&yScan >= 0 && yScan<render->display->yres)
				{
					zScan_new = (GzDepth)((-A*xScan - B*yScan - D) / C);
					GzGetDisplay(render->display, xScan, yScan, &r, &g, &b, &a, &zScan_old);
					if (zScan_new<zScan_old)
						GzPutDisplay(render->display, xScan, yScan, ctoi(render->flatcolor[0]), ctoi(render->flatcolor[1]), ctoi(render->flatcolor[2]), 1, zScan_new);
				}
				xScan++;
			}
			yScan++;
		}
		yScan = ceil(y2);
		while (yScan <= y1)	// Scan line from edge 13 to edge 12
		{
			xL = ceil(xK13*yScan + xB13);
			xR = floor(xK12*yScan + xB12);
			xScan = xL;
			while (xScan <= xR)
			{
				if (xScan >= 0 && xScan<render->display->xres&&yScan >= 0 && yScan<render->display->yres)
				{
					zScan_new = (GzDepth)((-A*xScan - B*yScan - D) / C);
					GzGetDisplay(render->display, xScan, yScan, &r, &g, &b, &a, &zScan_old);
					if (zScan_new<zScan_old)
						GzPutDisplay(render->display, xScan, yScan, ctoi(render->flatcolor[0]), ctoi(render->flatcolor[1]), ctoi(render->flatcolor[2]), 1, zScan_new);
				}
				xScan++;
			}
			yScan++;
		}

	}
		else if(xM>x2)		// Middle point 2 is on the left side of the long edge. Then the left point of scan line should be on edge 12 then edge 23, the right point of the scan line should be on edge 13
		{
			yScan=ceil(y3);
			while (yScan <=y2)	// Scan line from edge 23 to edge 13
			{
				xL = ceil(xK23*yScan + xB23);
				xR = floor(xK13*yScan + xB13);
				xScan = xL;
				while (xScan <= xR)
				{
					if (xScan >= 0 && xScan<render->display->xres&&yScan >= 0 && yScan<render->display->yres)
					{
						zScan_new = (GzDepth)((-A*xScan - B*yScan - D) / C);
						GzGetDisplay(render->display, xScan, yScan, &r, &g, &b, &a, &zScan_old);
						if (zScan_new<zScan_old)
							GzPutDisplay(render->display, xScan, yScan, ctoi(render->flatcolor[0]), ctoi(render->flatcolor[1]), ctoi(render->flatcolor[2]), 1, zScan_new);
					}
					xScan++;
				}
				yScan++;
			}
			yScan =ceil( y2);
			while(yScan<=y1)	// Scan line from edge 12 to edge 13
			{
				xL=ceil(xK12*yScan+xB12);
				xR = floor(xK13*yScan + xB13);
				xScan=xL;
				while(xScan<=xR)
				{
					if(xScan>=0&&xScan<render->display->xres&&yScan>=0&&yScan<render->display->yres)
					{
						zScan_new=(GzDepth)((-A*xScan-B*yScan-D)/C);
						GzGetDisplay(render->display, xScan, yScan, &r, &g, &b, &a, &zScan_old);
						if(zScan_new<zScan_old)		//  overwrite the pixel
							GzPutDisplay(render->display, xScan, yScan, ctoi(render->flatcolor[0]), ctoi(render->flatcolor[1]), ctoi(render->flatcolor[2]), 1, zScan_new);
					}
					xScan++;
				}
				yScan++;
			}
		
		}
		
	}
	

	return GZ_SUCCESS;
}

