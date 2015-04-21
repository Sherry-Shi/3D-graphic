/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"

GzColor	*image=NULL;
int xs, ys;
int reset = 1;

/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
	unsigned char		pixel[3];
	unsigned char     dummy;
	char  		foo[8];
	int   		i, j;
	FILE			*fd;

	if (reset) {          /* open and load texture file */
	fd = fopen ("texture", "rb");
	if (fd == NULL) {
		fprintf (stderr, "texture file not found\n");
		exit(-1);
	}
	fscanf (fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
	image = (GzColor*)malloc(sizeof(GzColor)*(xs+1)*(ys+1));
	if (image == NULL) {
		fprintf (stderr, "malloc for texture image failed\n");
		exit(-1);
	}

	for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */	// ??? (xs+1)*(ys+1)
		fread(pixel, sizeof(pixel), 1, fd);
		image[i][RED] = (float)((int)pixel[RED]) * (1.0 / 255.0);
		image[i][GREEN] = (float)((int)pixel[GREEN]) * (1.0 / 255.0);
		image[i][BLUE] = (float)((int)pixel[BLUE]) * (1.0 / 255.0);
		}

	reset = 0;          /* init is done */
	fclose(fd);
	}

/* bounds-test u,v to make sure nothing will overflow image array bounds */
/* determine texture cell corner values and perform bilinear interpolation */
/* set color to interpolated GzColor value and return */
	u *= (xs-1);
	v *= (ys-1);	// ??? xs, ys
	if(u>=0 && u<=xs && v>=0 && v<=ys)
	{
		int x, y;
		x = (int)u;
		y = (int)v;
		float s = u-x, t = v-y;
		GzColor up_left, up_right, down_left, down_right;
		up_left[0]=image[y*xs+x][0];up_left[1]=image[y*xs+x][1];up_left[2]=image[y*xs+x][2];
		up_right[0]=image[y*xs+x+1][0];up_right[1]=image[y*xs+x+1][1];up_right[2]=image[y*xs+x+1][2];
		down_left[0]=image[(y+1)*xs+x][0];down_left[1]=image[(y+1)*xs+x][1];down_left[2]=image[(y+1)*xs+x][2];
		down_right[0]=image[(y+1)*xs+x+1][0];down_right[1]=image[(y+1)*xs+x+1][1];down_right[2]=image[(y+1)*xs+x+1][2];

		color[0] = s*t*down_right[0] + (1-s)*t*down_left[0] + s*(1-t)*up_right[0] + (1-s)*(1-t)*up_left[0];
		color[1] = s*t*down_right[1] + (1-s)*t*down_left[1] + s*(1-t)*up_right[1] + (1-s)*(1-t)*up_left[1];
		color[2] = s*t*down_right[2] + (1-s)*t*down_left[2] + s*(1-t)*up_right[2] + (1-s)*(1-t)*up_left[2];

		return GZ_SUCCESS;
	}
	else
		return GZ_FAILURE;
}

/* Procedural texture function */
int ptex_fun(float u, float v, GzColor color)
{
	float u0 = u * 4;
	int x = (int)u0;
	switch(x%4)
	{
	case 0:
		color[0]=1;color[1]=0;color[2]=0;break;
	case 1:
		color[0]=0;color[1]=1;color[2]=0;break;
	case 2:
		color[0]=0;color[1]=0;color[2]=1;break;
	case 3:
		color[0]=1;color[1]=1;color[2]=1;break;
	}
	
	return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeTexture()
{
	if(image!=NULL)
		free(image);
	return GZ_SUCCESS;
}

