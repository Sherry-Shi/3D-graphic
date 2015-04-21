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
	float sine, cosine;
	sine=sin(degree/180*3.14159265358979323846);
	cosine=cos(degree/180*3.14159265358979323846);
	mat[1][1]=cosine;
	mat[1][2]=-sine;
	mat[2][1]=sine;
	mat[2][2]=cosine;
	return GZ_SUCCESS;
}

int GzRotYMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along y axis
// Pass back the matrix using mat value
	float sine, cosine;
	sine=sin(degree/180*3.14159265358979323846);
	cosine=cos(degree/180*3.14159265358979323846);
	mat[0][0]=cosine;
	mat[0][2]=sine;
	mat[2][0]=-sine;
	mat[2][2]=cosine;
	return GZ_SUCCESS;
}


int GzRotZMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along z axis
// Pass back the matrix using mat value
	float sine, cosine;
	sine=sin(degree/180*3.14159265358979323846);
	cosine=cos(degree/180*3.14159265358979323846);
	mat[0][0]=cosine;
	mat[0][1]=-sine;
	mat[1][0]=sine;
	mat[1][1]=cosine;

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

int validateFBColor(GzColor color, GzIntensity *r, GzIntensity *g, GzIntensity *b)
{
	GzColor tempColor;
	memcpy(tempColor, color, sizeof(GzColor));
	if(tempColor[0]>1.0)
		tempColor[0] = 1.0;
	if(tempColor[1]>1.0)
		tempColor[1] = 1.0;
	if(tempColor[2]>1.0)
		tempColor[2] = 1.0;
	*r = ctoi(tempColor[0]);
	*g = ctoi(tempColor[1]);
	*b = ctoi(tempColor[2]);
	return 0;
}


int putColor(float left, float right, float current, GzColor colorLeft, GzColor colorRight, GzColor colorCurrent)
{
	colorCurrent[0] = colorRight[0] - (colorRight[0]-colorLeft[0])*(right-current)/(right-left);
	colorCurrent[1] = colorRight[1] - (colorRight[1]-colorLeft[1])*(right-current)/(right-left);
	colorCurrent[2] = colorRight[2] - (colorRight[2]-colorLeft[2])*(right-current)/(right-left);
	return 0;
}
int normalize(GzMatrix mat)
{
	float magnitude = sqrt(mat[0][0]*mat[0][0] + mat[0][1]*mat[0][1] + mat[0][2]*mat[0][2]);
	mat[0][0] /= magnitude;
	mat[0][1] /= magnitude;
	mat[0][2] /= magnitude;
	magnitude = sqrt(mat[1][0]*mat[1][0] + mat[1][1]*mat[1][1] + mat[1][2]*mat[1][2]);
	mat[1][0] /= magnitude;
	mat[1][1] /= magnitude;
	mat[1][2] /= magnitude;
	magnitude = sqrt(mat[2][0]*mat[2][0] + mat[2][1]*mat[2][1] + mat[2][2]*mat[2][2]);
	mat[2][0] /= magnitude;
	mat[2][1] /= magnitude;
	mat[2][2] /= magnitude;
	return 0;
}


int combinedLPA(GzRender *render, GzCoord normal, GzColor color)
{
	GzColor spe, dif, amb;
	GzCoord R, tempNormal;
	float r, g, b;
	GzLight light;
	memset(spe, 0, sizeof(GzColor));
	memset(dif, 0, sizeof(GzColor));
	memset(amb, 0, sizeof(GzColor));

	for(int i=0;i<render->numlights;i++)
	{
		memcpy(&light, &render->lights[i], sizeof(GzLight));
		
		// spe
		// check N.L N.E
		float nl = normal[0]*light.direction[0]+normal[1]*light.direction[1]+normal[2]*light.direction[2];
		float ne = normal[0]*0 + normal[1]*0 + normal[2]*-1;
		if(nl*ne<0)
			continue;
		else if(nl<0 && ne<0)
		{
			tempNormal[0] = -normal[0];
			tempNormal[1] = -normal[1];
			tempNormal[2] = -normal[2];

			nl = tempNormal[0]*light.direction[0]+tempNormal[1]*light.direction[1]+tempNormal[2]*light.direction[2];
			ne = tempNormal[0]*0 + tempNormal[1]*0 + tempNormal[2]*(-1);
		}
		// R = 2(N.L)N-L
		R[0] = 2*nl*normal[0] - light.direction[0];
		R[1] = 2*nl*normal[1] - light.direction[1];
		R[2] = 2*nl*normal[2] - light.direction[2];
		// R, E unit??? re[0, 1]???
		float re = R[0]*0 + R[1]*0 + R[2]*(-1);
		if(re<0)
			re = 0;
		if(re>1)
			re=1;
		spe[0] += light.color[0]*pow(re, render->spec);
		spe[1] += light.color[1]*pow(re, render->spec);
		spe[2] += light.color[2]*pow(re, render->spec);

		// dif
		// N.L
		dif[0] += light.color[0]*nl;
		dif[1] += light.color[1]*nl;
		dif[2] += light.color[2]*nl;
	}

	amb[0] = render->Ka[0]*render->ambientlight.color[0];
	amb[1] = render->Ka[1]*render->ambientlight.color[1];
	amb[2] = render->Ka[2]*render->ambientlight.color[2];

	color[0] = render->Ks[0]*spe[0] + render->Kd[0]*dif[0] + amb[0];
	color[1] = render->Ks[1]*spe[1] + render->Kd[1]*dif[1] + amb[1];
	color[2] = render->Ks[2]*spe[2] + render->Kd[2]*dif[2] + amb[2];

	return 0;
}
int addNorm(float *normX, float *normY, float *normZ, int x, int y, GzCoord normCurrent)
{
	float XX,YY,ZZ;
	XX=(-normX[0]*x-normX[1]*y-normX[3])/normX[2];
	YY=(-normY[0]*x-normY[1]*y-normY[3])/normY[2];
	ZZ=(-normZ[0]*x-normZ[1]*y-normZ[3])/normZ[2];
	float magg=sqrt(XX*XX+YY*YY+ZZ*ZZ);
	XX/=magg;
	YY/=magg;
	ZZ/=magg;
	normCurrent[0]=XX;normCurrent[1]=YY;normCurrent[2]=ZZ;

	return 0;
}
void shade(GzCoord norm, GzCoord color)
{
  GzCoord	light;
  float		coef;

  light[0] = 0.707f;
  light[1] = 0.5f;
  light[2] = 0.5f;

  coef = light[0]*norm[0] + light[1]*norm[1] + light[2]*norm[2];
  if (coef < 0) 	coef *= -1;

  if (coef > 1.0)	coef = 1.0;
  color[0] = coef*0.95f;
  color[1] = coef*0.65f;
  color[2] = coef*0.88f;
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
	GzMatrix tempMatrix;
	memcpy(tempMatrix, matrix, sizeof(GzMatrix));
	 if(render->matlevel == -1)
		{
			memcpy(render->Ximage[0], matrix, sizeof(GzMatrix));
			render->matlevel += 1;
		}
	else if(render->matlevel < MATLEVELS-1)
	{
		if(render->matlevel != -1)
		{
			short tos = render->matlevel + 1;
			for(int i=0; i<4; i++)
			{
				for(int j=0; j<4; j++)
				{
					render->Ximage[tos][i][j] = render->Ximage[tos-1][i][0]*matrix[0][j]+render->Ximage[tos-1][i][1]*matrix[1][j]+render->Ximage[tos-1][i][2]*matrix[2][j]+render->Ximage[tos-1][i][3]*matrix[3][j];
				}
			}
			normalize(tempMatrix);
			if(render->matlevel > 1)
			{
				tempMatrix[0][3] = 0;
				tempMatrix[1][3] = 0;
				tempMatrix[2][3] = 0;
				for(int i=0; i<4; i++)
				{
					for(int j=0; j<4; j++)
					{
						float tmmmp= render->Xnorm[render->matlevel-2][i][0]*tempMatrix[0][j]+
									 render->Xnorm[render->matlevel-2][i][1]*tempMatrix[1][j]+
									 render->Xnorm[render->matlevel-2][i][2]*tempMatrix[2][j]+
									 render->Xnorm[render->matlevel-2][i][3]*tempMatrix[3][j];
						render->Xnorm[render->matlevel-1][i][j] = tmmmp;
					}
				}
			}
			if(render->matlevel == 1)
			{
				tempMatrix[0][3] = 0;
				tempMatrix[1][3] = 0;
				tempMatrix[2][3] = 0;
				memcpy(render->Xnorm[0], tempMatrix, sizeof(GzMatrix));
			}

			render->matlevel = tos;
		}
	
		
		return GZ_SUCCESS;
	}
		return GZ_FAILURE;	
}

int GzPopMatrix(GzRender *render)
{
/*
- pop a matrix off the Ximage stack
- check for stack underflow
*/
	if(render->matlevel > -1)
	{
		render->matlevel -= 1;
		return GZ_SUCCESS;
	}
	else
		return GZ_FAILURE;
}


int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer	*valueList) /* void** valuelist */
{
/*
- set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
- later set shaders, interpolaters, texture maps, and lights
*/
	for(int i=0;i<numAttributes;i++)
	{
		switch(nameList[i])
		{	
		case GZ_DIFFUSE_COEFFICIENT:
			render->Kd[0] = (*(GzColor*)valueList[i])[0];
			render->Kd[1] = (*(GzColor*)valueList[i])[1];
			render->Kd[2] = (*(GzColor*)valueList[i])[2];
			break;
		case GZ_INTERPOLATE:
			render->interp_mode = *(int*)valueList[i];
			break;
		case GZ_AMBIENT_COEFFICIENT:
			render->Ka[0] = (*(GzColor*)valueList[i])[0];
			render->Ka[1] = (*(GzColor*)valueList[i])[1];
			render->Ka[2] = (*(GzColor*)valueList[i])[2];
			break;
		
		case GZ_DIRECTIONAL_LIGHT:
			render->lights[render->numlights] = *(GzLight*)valueList[i];
			render->numlights++;
			break;
		case GZ_AMBIENT_LIGHT:
			render->ambientlight = *(GzLight*)valueList[i];
			break;
		case GZ_SPECULAR_COEFFICIENT:
			render->Ks[0] = (*(GzColor*)valueList[i])[0];
			render->Ks[1] = (*(GzColor*)valueList[i])[1];
			render->Ks[2] = (*(GzColor*)valueList[i])[2];
			break;
		case GZ_DISTRIBUTION_COEFFICIENT:
			render->spec = *(float*)valueList[i];
			break;
		case GZ_RGB_COLOR:
			memcpy(render->flatcolor, valueList[i], sizeof(GzColor));
			break;
		}
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
	float x1=((GzCoord*)valueList[0])[0][0];
	float y1=((GzCoord*)valueList[0])[0][1];
	float z1=((GzCoord*)valueList[0])[0][2];
	float nx1 = ((GzCoord*)valueList[1])[0][0];
	float ny1 = ((GzCoord*)valueList[1])[0][1];
	float nz1 = ((GzCoord*)valueList[1])[0][2];
	float x2=((GzCoord*)valueList[0])[1][0];
	float y2=((GzCoord*)valueList[0])[1][1];
	float z2=((GzCoord*)valueList[0])[1][2];
	float nx2 = ((GzCoord*)valueList[1])[1][0];
	float ny2 = ((GzCoord*)valueList[1])[1][1];
	float nz2 = ((GzCoord*)valueList[1])[1][2];
	float x3=((GzCoord*)valueList[0])[2][0];
	float y3=((GzCoord*)valueList[0])[2][1];
	float z3=((GzCoord*)valueList[0])[2][2];
	float nx3 = ((GzCoord*)valueList[1])[2][0];
	float ny3 = ((GzCoord*)valueList[1])[2][1];
	float nz3 = ((GzCoord*)valueList[1])[2][2];
	float temp;
	if(render->interp_mode == GZ_FLAT)
	{
		GzCoord normalFlat;
		normalFlat[0]=nx1;
		normalFlat[1]=ny1;
		normalFlat[2]=nz1;
		shade(normalFlat, render->flatcolor);
	}
	else
		memset(render->flatcolor, 0, sizeof(GzColor));

	// mat stands for the Matrix on top of the Stack
	GzMatrix mat;
	float d = 1/tan(render->camera.FOV/180*3.1416/2);
	memcpy(mat, render->Ximage[render->matlevel], sizeof(GzMatrix));
	GzMatrix matNorm;
	memcpy(matNorm, render->Xnorm[render->matlevel-2], sizeof(GzMatrix));

	float x01, x02, x03, y01, y02, y03, z01, z02, z03;
	float nx01, ny01, nz01, nx02, ny02, nz02, nx03, ny03, nz03;

	z01 = mat[2][0]*x1+mat[2][1]*y1+mat[2][2]*z1+mat[2][3];
	z02 = mat[2][0]*x2+mat[2][1]*y2+mat[2][2]*z2+mat[2][3];
	z03 = mat[2][0]*x3+mat[2][1]*y3+mat[2][2]*z3+mat[2][3];	

	if(z01>0&&z02>0&&z03>0)
	{
		// If all z values are greater than 0, calculate the x and y
		float w1, w2, w3;
		w1 = mat[3][0]*x1+mat[3][1]*y1+mat[3][2]*z1+mat[3][3];
		x01 = (mat[0][0]*x1+mat[0][1]*y1+mat[0][2]*z1+mat[0][3])/w1;
		y01 = (mat[1][0]*x1+mat[1][1]*y1+mat[1][2]*z1+mat[1][3])/w1;
		z01 /= w1;
		w2 = mat[3][0]*x2+mat[3][1]*y2+mat[3][2]*z2+mat[3][3];
		x02 = (mat[0][0]*x2+mat[0][1]*y2+mat[0][2]*z2+mat[0][3])/w2;
		y02 = (mat[1][0]*x2+mat[1][1]*y2+mat[1][2]*z2+mat[1][3])/w2;
		z02 /= w2;
		w3 = mat[3][0]*x3+mat[3][1]*y3+mat[3][2]*z3+mat[3][3];
		x03 = (mat[0][0]*x3+mat[0][1]*y3+mat[0][2]*z3+mat[0][3])/w3;
		y03 = (mat[1][0]*x3+mat[1][1]*y3+mat[1][2]*z3+mat[1][3])/w3;
		z03 /= w3;
		x1=x01;y1=y01;z1=z01;
		nx01 = matNorm[0][0]*nx1+matNorm[0][1]*ny1+matNorm[0][2]*nz1;
		ny01 = matNorm[1][0]*nx1+matNorm[1][1]*ny1+matNorm[1][2]*nz1;
		nz01 = matNorm[2][0]*nx1+matNorm[2][1]*ny1+matNorm[2][2]*nz1;
		x2=x02;y2=y02;z2=z02;
		nx02 = matNorm[0][0]*nx2+matNorm[0][1]*ny2+matNorm[0][2]*nz2;
		ny02 = matNorm[1][0]*nx2+matNorm[1][1]*ny2+matNorm[1][2]*nz2;
		nz02 = matNorm[2][0]*nx2+matNorm[2][1]*ny2+matNorm[2][2]*nz2;
		x3=x03;y3=y03;z3=z03;
		nx03 = matNorm[0][0]*nx3+matNorm[0][1]*ny3+matNorm[0][2]*nz3;
		ny03 = matNorm[1][0]*nx3+matNorm[1][1]*ny3+matNorm[1][2]*nz3;
		nz03 = matNorm[2][0]*nx3+matNorm[2][1]*ny3+matNorm[2][2]*nz3;
		nx1=nx01;ny1=ny01;nz1=nz01;
		nx2=nx02;ny2=ny02;nz2=nz02;
		nx3=nx03;ny3=ny03;nz3=nz03;

		// Switch points. Ensuring the top-left point is point1; the bottom-right point is point3.
		if(y1<y2)
		{
			temp=y1;y1=y2;y2=temp;
			temp=x1;x1=x2;x2=temp;
			temp=z1;z1=z2;z2=temp;

			temp=ny1;ny1=ny2;ny2=temp;
			temp=nx1;nx1=nx2;nx2=temp;
			temp=nz1;nz1=nz2;nz2=temp;
		}
		if(y1<y3)
		{
			temp=y1;y1=y3;y3=temp;
			temp=x1;x1=x3;x3=temp;
			temp=z1;z1=z3;z3=temp;

			temp=ny1;ny1=ny3;ny3=temp;
			temp=nx1;nx1=nx3;nx3=temp;
			temp=nz1;nz1=nz3;nz3=temp;
		}
		if(y2<y3)
		{
			temp=y2;y2=y3;y3=temp;
			temp=x2;x2=x3;x3=temp;
			temp=z2;z2=z3;z3=temp;

			temp=ny2;ny2=ny3;ny3=temp;
			temp=nx2;nx2=nx3;nx3=temp;
			temp=nz2;nz2=nz3;nz3=temp;
		}
		
		GzColor color1, color2, color3, flatColor;
		GzCoord point1, normal1, point2, normal2, point3, normal3;
		point1[0] = x1; point1[1] = y1; point1[2] = z1;
		normal1[0] = nx1; normal1[1] = ny1; normal1[2] = nz1;
		point2[0] = x2; point2[1] = y2; point2[2] = z2;
		normal2[0] = nx2; normal2[1] = ny2; normal2[2] = nz2;
		point3[0] = x3; point3[1] = y3; point3[2] = z3;
		normal3[0] = nx3; normal3[1] = ny3; normal3[2] = nz3;
		float A, B, C, D;
		A=y1*z2-y1*z3-y2*z1+y2*z3+y3*z1-y3*z2;
		B=-x1*z2+x1*z3+x2*z1-x2*z3-x3*z1+x3*z2;
		C=x1*y2-x1*y3-x2*y1+x2*y3+x3*y1-x3*y2;
		D=-x1*y2*z3+x1*y3*z2+x2*y1*z3-x2*y3*z1-x3*y1*z2+x3*y2*z1;
		float xK12, xK23, xK13, xB12=0, xB13=0, xB23=0;
		int hori_edge = 0;
		int xScan, yScan, xL, xR;
		GzIntensity r, g, b, a;
		GzDepth zScan_new, zScan_old;
		float xM;
		GzColor colorLeft, colorRight, colorCurrent;
		GzCoord normLeft, normRight, normCurrent;
		float magnitude;
		if(y1==y3)
		{
			hori_edge=13;
		}
		if(y1==y2)
		{
			hori_edge=12;
		}
		if(y2==y3)
		{
			hori_edge=23;
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
		combinedLPA(render, normal1, color1);
		combinedLPA(render, normal2, color2);
		combinedLPA(render, normal3, color3);
		float ANormZ, BNormZ, CNormZ, DNormZ;
		ANormZ=y1*nz2-y1*nz3-y2*nz1+y2*nz3+y3*nz1-y3*nz2;
		BNormZ=-x1*nz2+x1*nz3+x2*nz1-x2*nz3-x3*nz1+x3*nz2;
		CNormZ=x1*y2-x1*y3-x2*y1+x2*y3+x3*y1-x3*y2;
		DNormZ=-x1*y2*nz3+x1*y3*nz2+x2*y1*nz3-x2*y3*nz1-x3*y1*nz2+x3*y2*nz1;
		float normZ[4] = {ANormZ, BNormZ, CNormZ, DNormZ};
		float ANormX, BNormX, CNormX, DNormX;
		ANormX=y1*nx2-y1*nx3-y2*nx1+y2*nx3+y3*nx1-y3*nx2;
		BNormX=-x1*nx2+x1*nx3+x2*nx1-x2*nx3-x3*nx1+x3*nx2;
		CNormX=x1*y2-x1*y3-x2*y1+x2*y3+x3*y1-x3*y2;
		DNormX=-x1*y2*nx3+x1*y3*nx2+x2*y1*nx3-x2*y3*nx1-x3*y1*nx2+x3*y2*nx1;
		float normX[4] = {ANormX, BNormX, CNormX, DNormX};
		float ANormY, BNormY, CNormY, DNormY;
		ANormY=y1*ny2-y1*ny3-y2*ny1+y2*ny3+y3*ny1-y3*ny2;
		BNormY=-x1*ny2+x1*ny3+x2*ny1-x2*ny3-x3*ny1+x3*ny2;
		CNormY=x1*y2-x1*y3-x2*y1+x2*y3+x3*y1-x3*y2;
		DNormY=-x1*y2*ny3+x1*y3*ny2+x2*y1*ny3-x2*y3*ny1-x3*y1*ny2+x3*y2*ny1;
		float normY[4] = {ANormY, BNormY, CNormY, DNormY};
//	    if(hori_edge==0)	// No horizontal edges
//		{
//			xM = y2*xK13+xB13;
//			if(xM<x2)		// Middle point 2 is on the right side of the long edge. Then the left point of scan line should be on edge 13, the right point of the scan line should be on edge 12 then edge 23
//			{
//				yScan=int(y1);
//				while(yScan>=y2)	// Scan line from edge 13 to edge 12
//				{
//					xL=ceil(xK13*yScan+xB13);
//					xR=floor(xK12*yScan+xB12);
//
//					if(render->interp_mode == GZ_COLOR)
//					{
//						colorLeft[0] = color1[0]-(color1[0]-color2[0])*(y1-yScan)/(y1-y2);
//						colorLeft[1] = color1[1]-(color1[1]-color2[1])*(y1-yScan)/(y1-y2);
//						colorLeft[2] = color1[2]-(color1[2]-color2[2])*(y1-yScan)/(y1-y2);
//
//						colorRight[0] = color1[0]-(color1[0]-color3[0])*(y1-yScan)/(y1-y3);
//						colorRight[1] = color1[1]-(color1[1]-color3[1])*(y1-yScan)/(y1-y3);
//						colorRight[2] = color1[2]-(color1[2]-color3[2])*(y1-yScan)/(y1-y3);
//					}
//
//					xScan=xL;
//					while(xScan<=xR)
//					{
//						if(yScan==120&&xScan==148)
//							yScan=yScan;
//						if(xScan>=0&&xScan<=render->display->xres&&yScan>=0&&yScan<=render->display->yres)
//						{
//							zScan_new=(GzDepth)((-A*xScan-B*yScan-D)/C);
//							GzGetDisplay(render->display, xScan, yScan, &r, &g, &b, &a, &zScan_old);
//							if(zScan_new<=zScan_old)
//							{
//								if(render->interp_mode == GZ_COLOR)
//								{
//									putColor(xK12*yScan+xB12, xK13*yScan+xB13, (float)xScan, colorLeft, colorRight, colorCurrent);
//								}
//								else if(render->interp_mode == GZ_NORMALS)
//								{
//									// addNorm(); combinedLPA();
//									addNorm(normX, normY, normZ, xScan, yScan, normCurrent);
//									combinedLPA(render, normCurrent, colorCurrent);
//								}
//								else
//								{
//									// colorCurrent = render->flatcolor
//									memcpy(colorCurrent, render->flatcolor, sizeof(GzColor));
//								}
//								validateFBColor(colorCurrent, &r, &g, &b);
//								GzPutDisplay(render->display, xScan, yScan, r, g, b, 1, zScan_new);
//							}
//						}
//						xScan++;
//					}
//					yScan--;
//				}
//				while(yScan>=y3)	// Scan line form edge 13 to edge 23
//				{
//					xL=ceil(xK13*yScan+xB13);
//					xR=floor(xK23*yScan+xB23);
//
//					if(render->interp_mode == GZ_COLOR)
//					{
//						colorLeft[0] = color1[0]-(color1[0]-color2[0])*(y1-yScan)/(y1-y2);
//						colorLeft[1] = color1[1]-(color1[1]-color2[1])*(y1-yScan)/(y1-y2);
//						colorLeft[2] = color1[2]-(color1[2]-color2[2])*(y1-yScan)/(y1-y2);
//
//						colorRight[0] = color1[0]-(color1[0]-color3[0])*(y1-yScan)/(y1-y3);
//						colorRight[1] = color1[1]-(color1[1]-color3[1])*(y1-yScan)/(y1-y3);
//						colorRight[2] = color1[2]-(color1[2]-color3[2])*(y1-yScan)/(y1-y3);
//					}
//
//					xScan=xL;
//					while(xScan<=xR)
//					{
//						if(yScan==120&&xScan==148)
//							yScan=yScan;
//						if(xScan>=0&&xScan<=render->display->xres&&yScan>=0&&yScan<=render->display->yres)
//						{
//							zScan_new=(GzDepth)((-A*xScan-B*yScan-D)/C);
//							GzGetDisplay(render->display, xScan, yScan, &r, &g, &b, &a, &zScan_old);
//							if(zScan_new<=zScan_old)
//							{
//								if(render->interp_mode == GZ_COLOR)
//								{
//									putColor(xK12*yScan+xB12, xK13*yScan+xB13, (float)xScan, colorLeft, colorRight, colorCurrent);
//								}
//								else if(render->interp_mode == GZ_NORMALS)
//								{
//									// addNorm(); combinedLPA();
//									addNorm(normX, normY, normZ, xScan, yScan, normCurrent);
//									combinedLPA(render, normCurrent, colorCurrent);
//								}
//								else
//								{
//									// colorCurrent = render->flatcolor
//									memcpy(colorCurrent, render->flatcolor, sizeof(GzColor));
//								}
//								validateFBColor(colorCurrent, &r, &g, &b);
//								GzPutDisplay(render->display, xScan, yScan, r, g, b, 1, zScan_new);
//							}
//						}
//						xScan++;
//					}
//					yScan--;
//				}
//			}
//			else if(xM>x2)		// Middle point 2 is on the left side of the long edge. Then the left point of scan line should be on edge 12 then edge 23, the right point of the scan line should be on edge 13
//			{
//				yScan=int(y1);
//				while(yScan>=y2)	// Scan line from edge 12 to edge 13
//				{
//					xL=ceil(xK12*yScan+xB12);
//					xR=floor(xK13*yScan+xB13);
//
//					if(render->interp_mode == GZ_COLOR)
//					{
//						colorLeft[0] = color1[0]-(color1[0]-color2[0])*(y1-yScan)/(y1-y2);
//						colorLeft[1] = color1[1]-(color1[1]-color2[1])*(y1-yScan)/(y1-y2);
//						colorLeft[2] = color1[2]-(color1[2]-color2[2])*(y1-yScan)/(y1-y2);
//
//						colorRight[0] = color1[0]-(color1[0]-color3[0])*(y1-yScan)/(y1-y3);
//						colorRight[1] = color1[1]-(color1[1]-color3[1])*(y1-yScan)/(y1-y3);
//						colorRight[2] = color1[2]-(color1[2]-color3[2])*(y1-yScan)/(y1-y3);
//					}
//
//					xScan=xL;
//					while(xScan<=xR)
//					{
//						if(yScan==120&&xScan==148)
//							yScan=yScan;
//						if(xScan>=0&&xScan<=render->display->xres&&yScan>=0&&yScan<=render->display->yres)
//						{
//							zScan_new=(GzDepth)((-A*xScan-B*yScan-D)/C);
//							GzGetDisplay(render->display, xScan, yScan, &r, &g, &b, &a, &zScan_old);
//							if(zScan_new<=zScan_old)		// If the new Z value is not bigger than the old one, overwrite the pixel
//							{
//								if(render->interp_mode == GZ_COLOR)
//								{
//									putColor(xK12*yScan+xB12, xK13*yScan+xB13, (float)xScan, colorLeft, colorRight, colorCurrent);
//								}
//								else if(render->interp_mode == GZ_NORMALS)
//								{
//									addNorm(normX, normY, normZ, xScan, yScan, normCurrent);
//									combinedLPA(render, normCurrent, colorCurrent);
//								}
//								else
//								{
//									// colorCurrent = render->flatcolor
//									memcpy(colorCurrent, render->flatcolor, sizeof(GzColor));
//								}
//								validateFBColor(colorCurrent, &r, &g, &b);
//								GzPutDisplay(render->display, xScan, yScan, r, g, b, 1, zScan_new);
//							}
//							
//						}
//						xScan++;
//					}
//					yScan--;
//				}
//				while(yScan>=y3)	// Scan line from edge 23 to edge 13
//				{
//					xL=ceil(xK23*yScan+xB23);
//					xR=floor(xK13*yScan+xB13);
//
//					if(render->interp_mode == GZ_COLOR)
//					{
//						colorLeft[0] = color1[0]-(color1[0]-color2[0])*(y1-yScan)/(y1-y2);
//						colorLeft[1] = color1[1]-(color1[1]-color2[1])*(y1-yScan)/(y1-y2);
//						colorLeft[2] = color1[2]-(color1[2]-color2[2])*(y1-yScan)/(y1-y2);
//
//						colorRight[0] = color1[0]-(color1[0]-color3[0])*(y1-yScan)/(y1-y3);
//						colorRight[1] = color1[1]-(color1[1]-color3[1])*(y1-yScan)/(y1-y3);
//						colorRight[2] = color1[2]-(color1[2]-color3[2])*(y1-yScan)/(y1-y3);
//					}
//
//					xScan=xL;
//					while(xScan<=xR)
//					{
//						if(yScan==120&&xScan==148)
//							yScan=yScan;
//						if(xScan>=0&&xScan<=render->display->xres&&yScan>=0&&yScan<=render->display->yres)
//						{
//							zScan_new=(GzDepth)((-A*xScan-B*yScan-D)/C);
//							GzGetDisplay(render->display, xScan, yScan, &r, &g, &b, &a, &zScan_old);
//							if(zScan_new<=zScan_old)
//							{
//								if(render->interp_mode == GZ_COLOR)
//								{
//									putColor(xK12*yScan+xB12, xK13*yScan+xB13, (float)xScan, colorLeft, colorRight, colorCurrent);
//								}
//								else if(render->interp_mode == GZ_NORMALS)
//								{
//									// addNorm(); combinedLPA();
//									addNorm(normX, normY, normZ, xScan, yScan, normCurrent);
//									combinedLPA(render, normCurrent, colorCurrent);
//								}
//								else
//								{
//									// colorCurrent = render->flatcolor
//									memcpy(colorCurrent, render->flatcolor, sizeof(GzColor));
//								}
//								validateFBColor(colorCurrent, &r, &g, &b);
//								GzPutDisplay(render->display, xScan, yScan, r, g, b, 1, zScan_new);
//							}
//						}
//						xScan++;
//					}
//					yScan--;
//				}
//			}
//		else if(hori_edge!=0)		// Horizontal edge exists
//		{
//		
//			if(hori_edge==23)	// Edge 23 is horizontal. Then the scan line should from edge 12 to edge 13
//			{
//				yScan=int(y1);
//				while(yScan>=y3)
//				{
//					xL=ceil(xK12*yScan+xB12);
//					xR=floor(xK13*yScan+xB13);
//
//					if(render->interp_mode == GZ_COLOR)
//					{
//						colorLeft[0] = color1[0]-(color1[0]-color2[0])*(y1-yScan)/(y1-y2);
//						colorLeft[1] = color1[1]-(color1[1]-color2[1])*(y1-yScan)/(y1-y2);
//						colorLeft[2] = color1[2]-(color1[2]-color2[2])*(y1-yScan)/(y1-y2);
//
//						colorRight[0] = color1[0]-(color1[0]-color3[0])*(y1-yScan)/(y1-y3);
//						colorRight[1] = color1[1]-(color1[1]-color3[1])*(y1-yScan)/(y1-y3);
//						colorRight[2] = color1[2]-(color1[2]-color3[2])*(y1-yScan)/(y1-y3);
//					}
//
//					xScan=xL;
//					while(xScan<=xR)
//					{
//						if(xScan>=0&&xScan<=render->display->xres&&yScan>=0&&yScan<=render->display->yres)
//						{
//							zScan_new=(GzDepth)((-A*xScan-B*yScan-D)/C);
//							GzGetDisplay(render->display, xScan, yScan, &r, &g, &b, &a, &zScan_old);
//							if(zScan_new<=zScan_old)
//							{
//								if(render->interp_mode == GZ_COLOR)
//								{
//									putColor(xK12*yScan+xB12, xK13*yScan+xB13, (float)xScan, colorLeft, colorRight, colorCurrent);
//								}
//								else if(render->interp_mode == GZ_NORMALS)
//								{
//									// addNorm(); combinedLPA();
//									addNorm(normX, normY, normZ, xScan, yScan, normCurrent);
//									combinedLPA(render, normCurrent, colorCurrent);
//								}
//								else
//								{
//									// colorCurrent = render->flatcolor
//									memcpy(colorCurrent, render->flatcolor, sizeof(GzColor));
//								}
//								validateFBColor(colorCurrent, &r, &g, &b);
//								GzPutDisplay(render->display, xScan, yScan, r, g, b, 1, zScan_new);
//							}
//						}
//						xScan++;
//					}
//					yScan--;
//				}
//			}
//				else if(hori_edge==12)		// Edge 12 is horizontal. Then the scan line should from edge 13 to edge 23.
//			{
//				yScan=int(y1);
//				while(yScan>=y3)
//				{
//					xL=ceil(xK13*yScan+xB13);
//					xR=floor(xK23*yScan+xB23);
//
//					if(render->interp_mode == GZ_COLOR)
//					{
//						colorLeft[0] = color1[0]-(color1[0]-color2[0])*(y1-yScan)/(y1-y2);
//						colorLeft[1] = color1[1]-(color1[1]-color2[1])*(y1-yScan)/(y1-y2);
//						colorLeft[2] = color1[2]-(color1[2]-color2[2])*(y1-yScan)/(y1-y2);
//
//						colorRight[0] = color1[0]-(color1[0]-color3[0])*(y1-yScan)/(y1-y3);
//						colorRight[1] = color1[1]-(color1[1]-color3[1])*(y1-yScan)/(y1-y3);
//						colorRight[2] = color1[2]-(color1[2]-color3[2])*(y1-yScan)/(y1-y3);
//					}
//
//					xScan=xL;
//					while(xScan<=xR)
//					{
//						if(xScan>=0&&xScan<=render->display->xres&&yScan>=0&&yScan<=render->display->yres)
//						{
//							zScan_new=(GzDepth)((-A*xScan-B*yScan-D)/C);
//							GzGetDisplay(render->display, xScan, yScan, &r, &g, &b, &a, &zScan_old);
//							if(zScan_new<=zScan_old)
//							{
//								if(render->interp_mode == GZ_COLOR)
//								{
//									putColor(xK12*yScan+xB12, xK13*yScan+xB13, (float)xScan, colorLeft, colorRight, colorCurrent);
//								}
//								else if(render->interp_mode == GZ_NORMALS)
//								{
//									// addNorm(); combinedLPA();
//									addNorm(normX, normY, normZ, xScan, yScan, normCurrent);
//									combinedLPA(render, normCurrent, colorCurrent);
//								}
//								else
//								{
//									// colorCurrent = render->flatcolor
//									memcpy(colorCurrent, render->flatcolor, sizeof(GzColor));
//								}
//								validateFBColor(colorCurrent, &r, &g, &b);
//								GzPutDisplay(render->display, xScan, yScan, r, g, b, 1, zScan_new);
//							}
//						}
//						xScan++;
//					}
//					yScan--;
//				}
//			}
//			// Edge 13 can not be horizontal, because 1 is the top-left point and 3 is the bottom-right point. If edge 13 is horizontal, the three points will be in a horizontal line
//		}
//			
//		}
//		
//	}
//
//	return GZ_SUCCESS;
//}
//
if(hori_edge!=0)		// Horizontal edge exists
	{
		if(hori_edge==12)		// Edge 12 is horizontal. Then the scan line should from edge 13 to edge 23.
		{
			yScan=ceil(y3);
			while(yScan<=y1)
			{
			
				if(render->interp_mode == GZ_COLOR)
				{
						colorLeft[0] = color1[0]-(color1[0]-color2[0])*(y1-yScan)/(y1-y2);
						colorLeft[1] = color1[1]-(color1[1]-color2[1])*(y1-yScan)/(y1-y2);
						colorLeft[2] = color1[2]-(color1[2]-color2[2])*(y1-yScan)/(y1-y2);

						colorRight[0] = color1[0]-(color1[0]-color3[0])*(y1-yScan)/(y1-y3);
						colorRight[1] = color1[1]-(color1[1]-color3[1])*(y1-yScan)/(y1-y3);
						colorRight[2] = color1[2]-(color1[2]-color3[2])*(y1-yScan)/(y1-y3);
				}
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
						{
							//GzPutDisplay(render->display, xScan, yScan, ctoi(render->flatcolor[0]), ctoi(render->flatcolor[1]), ctoi(render->flatcolor[2]), 1, zScan_new);
							if(render->interp_mode == GZ_COLOR)
							{
									putColor(xK12*yScan+xB12, xK13*yScan+xB13, (float)xScan, colorLeft, colorRight, colorCurrent);
							}
							else if(render->interp_mode == GZ_NORMALS)
							{
									// addNorm(); combinedLPA();
									addNorm(normX, normY, normZ, xScan, yScan, normCurrent);
									combinedLPA(render, normCurrent, colorCurrent);
							}
							else
							{
									// colorCurrent = render->flatcolor
									memcpy(colorCurrent, render->flatcolor, sizeof(GzColor));
							}
								validateFBColor(colorCurrent, &r, &g, &b);
								GzPutDisplay(render->display, xScan, yScan, r, g, b, 1, zScan_new);
						}
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
				if(render->interp_mode == GZ_COLOR)
					{
						colorLeft[0] = color1[0]-(color1[0]-color2[0])*(y1-yScan)/(y1-y2);
						colorLeft[1] = color1[1]-(color1[1]-color2[1])*(y1-yScan)/(y1-y2);
						colorLeft[2] = color1[2]-(color1[2]-color2[2])*(y1-yScan)/(y1-y2);

						colorRight[0] = color1[0]-(color1[0]-color3[0])*(y1-yScan)/(y1-y3);
						colorRight[1] = color1[1]-(color1[1]-color3[1])*(y1-yScan)/(y1-y3);
						colorRight[2] = color1[2]-(color1[2]-color3[2])*(y1-yScan)/(y1-y3);
					}
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
						{
							//GzPutDisplay(render->display, xScan, yScan, ctoi(render->flatcolor[0]), ctoi(render->flatcolor[1]), ctoi(render->flatcolor[2]), 1, zScan_new);
							if(render->interp_mode == GZ_COLOR)
								{
									putColor(xK12*yScan+xB12, xK13*yScan+xB13, (float)xScan, colorLeft, colorRight, colorCurrent);
								}
								else if(render->interp_mode == GZ_NORMALS)
								{
									// addNorm(); combinedLPA();
									addNorm(normX, normY, normZ, xScan, yScan, normCurrent);
									combinedLPA(render, normCurrent, colorCurrent);
								}
								else
								{
									// colorCurrent = render->flatcolor
									memcpy(colorCurrent, render->flatcolor, sizeof(GzColor));
								}
								validateFBColor(colorCurrent, &r, &g, &b);
								GzPutDisplay(render->display, xScan, yScan, r, g, b, 1, zScan_new);
						}
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
			if(render->interp_mode == GZ_COLOR)
					{
						colorLeft[0] = color1[0]-(color1[0]-color2[0])*(y1-yScan)/(y1-y2);
						colorLeft[1] = color1[1]-(color1[1]-color2[1])*(y1-yScan)/(y1-y2);
						colorLeft[2] = color1[2]-(color1[2]-color2[2])*(y1-yScan)/(y1-y2);

						colorRight[0] = color1[0]-(color1[0]-color3[0])*(y1-yScan)/(y1-y3);
						colorRight[1] = color1[1]-(color1[1]-color3[1])*(y1-yScan)/(y1-y3);
						colorRight[2] = color1[2]-(color1[2]-color3[2])*(y1-yScan)/(y1-y3);
					}
			xL = ceil(xK13*yScan + xB13);
			xR = floor(xK23*yScan + xB23);
			xScan = xL;
			while (xScan <= xR)
			{
				if(yScan==120&&xScan==148)
							yScan=yScan;
				if (xScan >= 0 && xScan<render->display->xres&&yScan >= 0 && yScan<render->display->yres)
				{
					zScan_new = (GzDepth)((-A*xScan - B*yScan - D) / C);
					GzGetDisplay(render->display, xScan, yScan, &r, &g, &b, &a, &zScan_old);
					if (zScan_new<zScan_old)
					{//
						if(render->interp_mode == GZ_COLOR)
								{
									putColor(xK12*yScan+xB12, xK13*yScan+xB13, (float)xScan, colorLeft, colorRight, colorCurrent);
								}
								else if(render->interp_mode == GZ_NORMALS)
								{
									// addNorm(); combinedLPA();
									addNorm(normX, normY, normZ, xScan, yScan, normCurrent);
									combinedLPA(render, normCurrent, colorCurrent);
								}
								else
								{
									// colorCurrent = render->flatcolor
									memcpy(colorCurrent, render->flatcolor, sizeof(GzColor));
								}
								validateFBColor(colorCurrent, &r, &g, &b);
								GzPutDisplay(render->display, xScan, yScan, r, g, b, 1, zScan_new);
					}			//GzPutDisplay(render->display, xScan, yScan, ctoi(render->flatcolor[0]), ctoi(render->flatcolor[1]), ctoi(render->flatcolor[2]), 1, zScan_new);
				}
				xScan++;
			}
			yScan++;
		}
		yScan = ceil(y2);
		while (yScan <= y1)	// Scan line from edge 13 to edge 12
		{
			if(render->interp_mode == GZ_COLOR)
					{
						colorLeft[0] = color1[0]-(color1[0]-color2[0])*(y1-yScan)/(y1-y2);
						colorLeft[1] = color1[1]-(color1[1]-color2[1])*(y1-yScan)/(y1-y2);
						colorLeft[2] = color1[2]-(color1[2]-color2[2])*(y1-yScan)/(y1-y2);

						colorRight[0] = color1[0]-(color1[0]-color3[0])*(y1-yScan)/(y1-y3);
						colorRight[1] = color1[1]-(color1[1]-color3[1])*(y1-yScan)/(y1-y3);
						colorRight[2] = color1[2]-(color1[2]-color3[2])*(y1-yScan)/(y1-y3);
					}
			xL = ceil(xK13*yScan + xB13);
			xR = floor(xK12*yScan + xB12);
			xScan = xL;
			while (xScan <= xR)
			{
			
				if(yScan==120&&xScan==148)
					yScan=yScan;
				if (xScan >= 0 && xScan<render->display->xres&&yScan >= 0 && yScan<render->display->yres)
				{
					zScan_new = (GzDepth)((-A*xScan - B*yScan - D) / C);
					GzGetDisplay(render->display, xScan, yScan, &r, &g, &b, &a, &zScan_old);
					if (zScan_new<zScan_old)
					{
						//GzPutDisplay(render->display, xScan, yScan, ctoi(render->flatcolor[0]), ctoi(render->flatcolor[1]), ctoi(render->flatcolor[2]), 1, zScan_new);
						if(zScan_new<=zScan_old)
							{
								if(render->interp_mode == GZ_COLOR)
								{
									putColor(xK12*yScan+xB12, xK13*yScan+xB13, (float)xScan, colorLeft, colorRight, colorCurrent);
								}
								else if(render->interp_mode == GZ_NORMALS)
								{
									// addNorm(); combinedLPA();
									addNorm(normX, normY, normZ, xScan, yScan, normCurrent);
									combinedLPA(render, normCurrent, colorCurrent);
								}
								else
								{
									// colorCurrent = render->flatcolor
									memcpy(colorCurrent, render->flatcolor, sizeof(GzColor));
								}
								validateFBColor(colorCurrent, &r, &g, &b);
								GzPutDisplay(render->display, xScan, yScan, r, g, b, 1, zScan_new);
							}
						}
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
				if(render->interp_mode == GZ_COLOR)
					{
						colorLeft[0] = color1[0]-(color1[0]-color2[0])*(y1-yScan)/(y1-y2);
						colorLeft[1] = color1[1]-(color1[1]-color2[1])*(y1-yScan)/(y1-y2);
						colorLeft[2] = color1[2]-(color1[2]-color2[2])*(y1-yScan)/(y1-y2);

						colorRight[0] = color1[0]-(color1[0]-color3[0])*(y1-yScan)/(y1-y3);
						colorRight[1] = color1[1]-(color1[1]-color3[1])*(y1-yScan)/(y1-y3);
						colorRight[2] = color1[2]-(color1[2]-color3[2])*(y1-yScan)/(y1-y3);
					}
				xL = ceil(xK23*yScan + xB23);
				xR = floor(xK13*yScan + xB13);
				xScan = xL;
				while (xScan <= xR)
				{
					if(yScan==120&&xScan==148)
							yScan=yScan;
					if (xScan >= 0 && xScan<render->display->xres&&yScan >= 0 && yScan<render->display->yres)
					{
						zScan_new = (GzDepth)((-A*xScan - B*yScan - D) / C);
						GzGetDisplay(render->display, xScan, yScan, &r, &g, &b, &a, &zScan_old);
						if(zScan_new<=zScan_old)		// If the new Z value is not bigger than the old one, overwrite the pixel
							{
								if(render->interp_mode == GZ_COLOR)
								{
									putColor(xK12*yScan+xB12, xK13*yScan+xB13, (float)xScan, colorLeft, colorRight, colorCurrent);
								}
								else if(render->interp_mode == GZ_NORMALS)
								{
									addNorm(normX, normY, normZ, xScan, yScan, normCurrent);
									combinedLPA(render, normCurrent, colorCurrent);
								}
								else
								{
									// colorCurrent = render->flatcolor
									memcpy(colorCurrent, render->flatcolor, sizeof(GzColor));
								}
								validateFBColor(colorCurrent, &r, &g, &b);
								GzPutDisplay(render->display, xScan, yScan, r, g, b, 1, zScan_new);
							}
					}
					xScan++;
				}
				yScan++;
			}
			yScan =ceil( y2);
			while(yScan<=y1)	// Scan line from edge 12 to edge 13
			{
			if(render->interp_mode == GZ_COLOR)
					{
						colorLeft[0] = color1[0]-(color1[0]-color2[0])*(y1-yScan)/(y1-y2);
						colorLeft[1] = color1[1]-(color1[1]-color2[1])*(y1-yScan)/(y1-y2);
						colorLeft[2] = color1[2]-(color1[2]-color2[2])*(y1-yScan)/(y1-y2);

						colorRight[0] = color1[0]-(color1[0]-color3[0])*(y1-yScan)/(y1-y3);
						colorRight[1] = color1[1]-(color1[1]-color3[1])*(y1-yScan)/(y1-y3);
						colorRight[2] = color1[2]-(color1[2]-color3[2])*(y1-yScan)/(y1-y3);
					}
				xL=ceil(xK12*yScan+xB12);
				xR = floor(xK13*yScan + xB13);
				xScan=xL;
				while(xScan<=xR)
				{
					if(yScan==120&&xScan==148)
							yScan=yScan;
						if(xScan>=0&&xScan<=render->display->xres&&yScan>=0&&yScan<=render->display->yres)
						{
							zScan_new=(GzDepth)((-A*xScan-B*yScan-D)/C);
							GzGetDisplay(render->display, xScan, yScan, &r, &g, &b, &a, &zScan_old);
							if(zScan_new<=zScan_old)
							{
								if(render->interp_mode == GZ_COLOR)
								{
									putColor(xK12*yScan+xB12, xK13*yScan+xB13, (float)xScan, colorLeft, colorRight, colorCurrent);
								}
								else if(render->interp_mode == GZ_NORMALS)
								{
									// addNorm(); combinedLPA();
									addNorm(normX, normY, normZ, xScan, yScan, normCurrent);
									combinedLPA(render, normCurrent, colorCurrent);
								}
								else
								{
									// colorCurrent = render->flatcolor
									memcpy(colorCurrent, render->flatcolor, sizeof(GzColor));
								}
								validateFBColor(colorCurrent, &r, &g, &b);
								GzPutDisplay(render->display, xScan, yScan, r, g, b, 1, zScan_new);
							}
						}
					xScan++;
				}
				yScan++;
			}
		
		}
		
	}
	

	return GZ_SUCCESS;
}