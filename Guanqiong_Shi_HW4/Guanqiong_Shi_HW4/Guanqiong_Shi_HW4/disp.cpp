/*   CS580 HW   */
#include    "stdafx.h"  
#include	"Gz.h"
#include	"disp.h"
#include	"math.h"

int GzNewFrameBuffer(char** framebuffer, int width, int height)
{
/* create a framebuffer:
 -- allocate memory for framebuffer : (sizeof)GzPixel x width x height
 -- pass back pointer 
 -- NOTE: this function is optional and not part of the API, but you may want to use it within the display function.
*/
	//framebuffer = new char[3*width*height];
	char *p=new char[3*sizeof(char)*width*height];
	*framebuffer = p;
	return GZ_SUCCESS;
}

int GzNewDisplay(GzDisplay	**display, int xRes, int yRes)
{
/* create a display:
  -- allocate memory for indicated resolution
  -- pass back pointer to GzDisplay object in display
*/
	if(xRes>MAXXRES || yRes>MAXYRES)
		return GZ_FAILURE;
	GzDisplay *p=new GzDisplay;
	*display=p;
	p->xres=xRes;
	p->yres=yRes;
	p->fbuf=new GzPixel[sizeof(GzPixel)*xRes*yRes];
		
	return GZ_SUCCESS;
}


int GzFreeDisplay(GzDisplay	*display)
{
/* clean up, free memory */
	if(display)
		delete display;
	return GZ_SUCCESS;
}


int GzGetDisplayParams(GzDisplay *display, int *xRes, int *yRes)
{
/* pass back values for a display */
	*xRes=(int)(display->xres);
	*yRes=(int)(display->yres);
	return GZ_SUCCESS;
}


int GzInitDisplay(GzDisplay	*display)
{
/* set everything to some default values - start a new frame */
	//memset(display->fbuf, 0, sizeof(display->xres*display->yres*sizeof(GzPixel)));
int bytes = sizeof(int);
	int int_max = int(pow(2.0,(int)(bytes*8))/2-1);
	int i,j;
	for(j=0;j<display->yres;j++)
	{
		for(i=0;i<display->xres;i++)
		{
			display->fbuf[i].red=2500;
			display->fbuf[i].green=970;
			display->fbuf[i].blue=0;
			display->fbuf[i].alpha=1;
			display->fbuf[j*display->xres+i].z=int_max;
		}
	}
	return GZ_SUCCESS;
}


int GzPutDisplay(GzDisplay *display, int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
/* write pixel values into the display */
	if(i>=0&&i<=display->xres&&j>=0&&j<=display->yres)
	{
		if(r>4095)
			r=4095;
		if(g>4095)
			g=4095;
		if(b>4095)
			b=4095;
		if(a>4095)
			a=4095;
		int offset = j*display->xres+i;
		display->fbuf[offset].red = r;
		display->fbuf[offset].green = g;
		display->fbuf[offset].blue = b;
		display->fbuf[offset].alpha = a;
		display->fbuf[offset].z = z;
		return GZ_SUCCESS;
	}
	else
		return GZ_FAILURE;
}


int GzGetDisplay(GzDisplay *display, int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
	/* pass back pixel value in the display */
	if(i>=0&&i<=display->xres&&j>=0&&j<=display->yres)
	{
		int offset = j*display->xres+i;
		*r = display->fbuf[offset].red;
		*g = display->fbuf[offset].green;
		*b = display->fbuf[offset].blue;
		*a = display->fbuf[offset].alpha;
		*z = display->fbuf[offset].z;
		return GZ_SUCCESS;
	}
	else
		return GZ_FAILURE;
}


int GzFlushDisplay2File(FILE* outfile, GzDisplay *display)
{
	/* write pixels to ppm file -- "P6 %d %d 255\r" */
	fprintf(outfile, "P6 %d %d 255\n", display->xres, display->yres);
	int	i, j;
	GzIntensity r, g, b, a;
	GzDepth z;
	unsigned char r_byte, g_byte, b_byte;
		for (j = 0; j < display->yres; j++) {
			for (i = 0; i < display->xres; i++) {
				GzGetDisplay(display, i, j, &r, &g, &b, &a, &z);
				r_byte = (unsigned char)(r>>4);
				fwrite(&r_byte, sizeof(unsigned char), 1, outfile);
				g_byte = (unsigned char)(g>>4);
				fwrite(&g_byte, sizeof(unsigned char), 1, outfile);
				b_byte = (unsigned char)(b>>4);
				fwrite(&b_byte, sizeof(unsigned char), 1, outfile);
			}
		}
	return GZ_SUCCESS;
}

int GzFlushDisplay2FrameBuffer(char* framebuffer, GzDisplay *display)
{

	/* write pixels to framebuffer: 
		- Put the pixels into the frame buffer
		- Caution: store the pixel to the frame buffer as the order of blue, green, and red 
		- Not red, green, and blue !!!
	*/
	//memset(framebuffer, 0, display->xres*display->yres*3);
	int i=0, j=0;
	char *offset = framebuffer;
	GzIntensity r, g, b, a;
	GzDepth z;
	for (j = 0; j < display->yres; j++) {
		for (i = 0; i < display->xres; i++) {
			GzGetDisplay(display, i, j, &r, &g, &b, &a, &z);
			*offset = (char)(b>>4);
			offset++;
			*offset = (char)(g>>4);
			offset++;
			*offset = (char)(r>>4);
			offset++;
		}
	}

	return GZ_SUCCESS;
}