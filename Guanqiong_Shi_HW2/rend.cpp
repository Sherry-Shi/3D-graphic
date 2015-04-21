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

int GzNewRender(GzRender **render, GzDisplay *display)
{
/* 
- malloc a renderer struct
- span interpolator needs pointer to display for pixel writes
*/
	GzRender *p = new GzRender;
	*render = p;
	p->display = display;

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


int GzBeginRender(GzRender	*render)
{
/* 
- set up for start of each frame - init frame buffer
*/
	//memset(render->display->fbuf, 255, sizeof(GzPixel)*render->display->xres*render->display->yres);
	if (render->display == NULL) return GZ_FAILURE;
	for (int i = 0; i < render->display->yres*render->display->xres; i++)
	{
		render->display->fbuf[i].red = 2500;
		render->display->fbuf[i].green = 970;
		render->display->fbuf[i].blue = 0;
		render->display->fbuf[i].alpha = 1;
		render->display->fbuf[i].z = 0;
	}
	// Initialize the fbuf and Z value
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

	return GZ_SUCCESS;
}


int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer *valueList) /* void** valuelist */
{
/*
- set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
- later set shaders, interpolaters, texture maps, and lights
*/
	render->flatcolor[0]=((float*)valueList[0])[0];
	render->flatcolor[1]=((float*)valueList[0])[1];
	render->flatcolor[2]=((float*)valueList[0])[2];
	
	return GZ_SUCCESS;
}


int GzPutTriangle(GzRender *render, int	numParts, GzToken *nameList,
	GzPointer *valueList) 
/* numParts - how many names and values */
{
/* 
- pass in a triangle description with tokens and values corresponding to
      GZ_NULL_TOKEN:		do nothing - no values
      GZ_POSITION:		3 vert positions in model space
- Invoke the scan converter and return an error code
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

	// Switch points. Ensuring the top-left point is point1; the bottom-right point is point3.
	// Switch points. Ensuring the top-left point is point1; the bottom-right point is point3.
	if(y1<y2)
	{
		temp=y1;y1=y2;y2=temp;
		temp=x1;x1=x2;x2=temp;
		temp=z1;z1=z2;z2=temp;
	}
	if(y1<y3)
	{
		temp=y1;y1=y3;y3=temp;
		temp=x1;x1=x3;x3=temp;
		temp=z1;z1=z3;z3=temp;
	}
	if(y2<y3)
	{
		temp=y2;y2=y3;y3=temp;
		temp=x2;x2=x3;x3=temp;
		temp=z2;z2=z3;z3=temp;
	}
	else if(y2==y3)
	{
		if(x2>x3)
		{
			temp=x2;x2=x3;x3=temp;
			temp=y2;y2=y3;y3=temp;
			temp=z2;z2=z3;z3=temp;
		}
	}
	if(y1==y2)
	{
		if(x1>x2)
		{
			temp=y1;y1=y2;y2=temp;
			temp=x1;x1=x2;x2=temp;
			temp=z1;z1=z2;z2=temp;
		}
	}
	// Switch points over


	float xK12, xK23, xK13, xB12=0, xB13=0, xB23=0;
	int hori_edge = 0;
	int xScan, yScan, xL, xR;
	GzIntensity r, g, b, a;
	GzDepth zScan_new, zScan_old;
	float xM;

	if(y1==y2)
	{
		hori_edge=12;
	}
	if(y2==y3)
	{
		hori_edge=23;
	}
	if(y1==y3)
	{
		hori_edge=13;
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


	if(hori_edge==0)	// No horizontal edges
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
	else if(hori_edge!=0)		// Horizontal edge exists
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

	return GZ_SUCCESS;
}



