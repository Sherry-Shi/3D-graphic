// Application5.cpp: implementation of the Application5 class.
//
//////////////////////////////////////////////////////////////////////

/*
 * application test code for homework assignment #5
*/

#include "stdafx.h"
#include "CS580HW.h"
#include "Application5.h"
#include "Gz.h"
#include "disp.h"
#include "rend.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[]=__FILE__;
#define new DEBUG_NEW
#endif

#define PHONG
#define TEX
//#define PROCEDUAL_TEX

#ifdef TEX
#define INFILE  "ppot.asc"
#define OUTFILE "output-texture.ppm"
#define APPCAMERA
#else
#define INFILE "pot4.asc"
#define OUTFILE "output-notexture.ppm"
#endif


extern int tex_fun(float u, float v, GzColor color); /* image texture function */
extern int ptex_fun(float u, float v, GzColor color); /* procedural texture function */

void shade(GzCoord norm, GzCoord color);

float AAFilter[AAKERNEL_SIZE][3] =  	/* each sample is defined by Xshift, Yshift, weight*/
	{  -0.52, 0.38, 0.128,                  0.41, 0.56, 0.119,                     0.27, 0.08, 0.294,
	-0.17, -0.29, 0.249,                    0.58, -0.55, 0.104,                   -0.31, -0.71, 0.106    };

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Application5::Application5()
{

}

Application5::~Application5()
{
	Clean();
}

int Application5::Initialize()
{
	GzCamera	camera;  
	int		    xRes, yRes;	/* display parameters */ 

	GzToken		nameListShader[9]; 	    /* shader attribute names */
	GzPointer   valueListShader[9];		/* shader attribute pointers */
	GzToken     nameListLights[10];		/* light info */
	GzPointer   valueListLights[10];
	int			shaderType, interpStyle;
	float		specpower;
	int		status; 
 
	status = 0; 

	/* 
	 * Allocate memory for user input
	 */
	m_pUserInput = new GzInput;

	/* 
	 * initialize the display and the renderer 
	 */ 
 	m_nWidth = 256;		// frame buffer and display width
	m_nHeight = 256;    // frame buffer and display height
	m_renderNum = AAKERNEL_SIZE;

	GzDisplay *pd;
	GzRender *pr;
	for(int i=0;i<m_renderNum;i++)
	{
		pd = new GzDisplay;
		m_pDisplay[i] = pd;
		pr = new GzRender;
		m_pRender[i] = pr;
	}
	//*m_pDisplay = new GzDisplay[m_renderNum];
	//*m_pRender = new GzRender[m_renderNum];

	status |= GzNewFrameBuffer(&m_pFrameBuffer, m_nWidth, m_nHeight);

	for(int i=0; i<m_renderNum; i++)
	{
		status |= GzNewDisplay(&m_pDisplay[i], m_nWidth, m_nHeight);

		status |= GzGetDisplayParams(m_pDisplay[i], &xRes, &yRes); 
 
		status |= GzNewRender(&m_pRender[i], m_pDisplay[i]); 
	}

/* Translation matrix */
GzMatrix	scale = 
{ 
	3.25,	0.0,	0.0,	0.0, 
	0.0,	3.25,	0.0,	-3.25, 
	0.0,	0.0,	3.25,	3.5, 
	0.0,	0.0,	0.0,	1.0 
}; 
 
GzMatrix	rotateX = 
{ 
	1.0,	0.0,	0.0,	0.0, 
	0.0,	.7071,	.7071,	0.0, 
	0.0,	-.7071,	.7071,	0.0, 
	0.0,	0.0,	0.0,	1.0 
}; 
 
GzMatrix	rotateY = 
{ 
	.866,	0.0,	-0.5,	0.0, 
	0.0,	1.0,	0.0,	0.0, 
	0.5,	0.0,	.866,	0.0, 
	0.0,	0.0,	0.0,	1.0 
}; 

#ifdef APPCAMERA 	/* set up app-defined camera if desired, else use camera defaults */
    camera.position[X] = -3;
    camera.position[Y] = -25;
    camera.position[Z] = -4;

    camera.lookat[X] = 7.8;
    camera.lookat[Y] = 0.7;
    camera.lookat[Z] = 6.5;

    camera.worldup[X] = -0.2;
    camera.worldup[Y] = 1.0;
    camera.worldup[Z] = 0.0;

    camera.FOV = 63.7;              /* degrees *              /* degrees */

	for(int i=0; i<m_renderNum;i++)
	{
		status |= GzPutCamera(m_pRender[i], &camera); 
	}
#endif 

	/* Start Renderer */
	for(int i=0;i<m_renderNum;i++)
	{
		status |= GzBeginRender(m_pRender[i]);
	}

	/* Light */
	GzLight	light1 = { {-0.7071, 0.7071, 0}, {0.5, 0.5, 0.9} };
	GzLight	light2 = { {0, -0.7071, -0.7071}, {0.9, 0.2, 0.3} };
	GzLight	light3 = { {0.7071, 0.0, -0.7071}, {0.2, 0.7, 0.3} };
	GzLight	ambientlight = { {0, 0, 0}, {0.3, 0.3, 0.3} };

	/* Material property */
	GzColor specularCoefficient = { 0.3, 0.3, 0.3 };
	GzColor ambientCoefficient = { 0.1, 0.1, 0.1 };
	GzColor diffuseCoefficient = {0.7, 0.7, 0.7};

/* 
  renderer is ready for frame --- define lights and shader at start of frame 
*/

        /*
         * Tokens associated with light parameters
         */
        nameListLights[0] = GZ_DIRECTIONAL_LIGHT;
        valueListLights[0] = (GzPointer)&light1;
        nameListLights[1] = GZ_DIRECTIONAL_LIGHT;
        valueListLights[1] = (GzPointer)&light2;
        nameListLights[2] = GZ_DIRECTIONAL_LIGHT;
        valueListLights[2] = (GzPointer)&light3;
		for(int i=0; i<m_renderNum; i++)
			status |= GzPutAttribute(m_pRender[i], 3, nameListLights, valueListLights);

        nameListLights[0] = GZ_AMBIENT_LIGHT;
        valueListLights[0] = (GzPointer)&ambientlight;
        for(int i=0; i<m_renderNum; i++)
			status |= GzPutAttribute(m_pRender[i], 1, nameListLights, valueListLights);

        /*
         * Tokens associated with shading 
         */
		int tokenNum = 0;
        nameListShader[0]  = GZ_DIFFUSE_COEFFICIENT;
        valueListShader[0] = (GzPointer)diffuseCoefficient;
		tokenNum++;
	/* 
	* Select either GZ_COLOR or GZ_NORMALS as interpolation mode  
	*/
        nameListShader[1]  = GZ_INTERPOLATE;
#ifdef PHONG
        interpStyle = GZ_NORMALS;         /* Phong shading */
#else
		interpStyle = GZ_COLOR;
#endif
        valueListShader[1] = (GzPointer)&interpStyle;
		tokenNum++;
        nameListShader[2]  = GZ_AMBIENT_COEFFICIENT;
        valueListShader[2] = (GzPointer)ambientCoefficient;
		tokenNum++;
        nameListShader[3]  = GZ_SPECULAR_COEFFICIENT;
        valueListShader[3] = (GzPointer)specularCoefficient;
		tokenNum++;
        nameListShader[4]  = GZ_DISTRIBUTION_COEFFICIENT;
        specpower = 32;
        valueListShader[4] = (GzPointer)&specpower;
		tokenNum++;

#ifdef TEX
        nameListShader[5]  = GZ_TEXTURE_MAP;
   /* set up null texture function or valid pointer */
#ifdef PROCEDUAL_TEX
        valueListShader[5] = (GzPointer)(ptex_fun);
#else
        valueListShader[5] = (GzPointer)(tex_fun);	/* or use ptex_fun */
#endif
		tokenNum++;
#endif
        for(int i=0; i<m_renderNum; i++)
			status |= GzPutAttribute(m_pRender[i], tokenNum, nameListShader, valueListShader);

		for(int i=0; i<m_renderNum; i++)
		{
			nameListShader[0]  = GZ_AASHIFTX;
			valueListShader[0] = (GzPointer)&AAFilter[i][0];
			nameListShader[1]  = GZ_AASHIFTY;
			valueListShader[1] = (GzPointer)&AAFilter[i][1];
			status |= GzPutAttribute(m_pRender[i], 2, nameListShader, valueListShader);
		}

	for(int i=0; i<m_renderNum; i++)
	{
		status |= GzPushMatrix(m_pRender[i], scale);  
		status |= GzPushMatrix(m_pRender[i], rotateY); 
		status |= GzPushMatrix(m_pRender[i], rotateX); 
	}

	if (status) exit(GZ_FAILURE); 

	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS); 
}

int Application5::Render() 
{
	GzToken		nameListTriangle[3]; 	/* vertex attribute names */
	GzPointer	valueListTriangle[3]; 	/* vertex attribute pointers */
	GzCoord		vertexList[3];	/* vertex position coordinates */ 
	GzCoord		normalList[3];	/* vertex normals */ 
	GzTextureIndex  	uvList[3];		/* vertex texture map indices */ 
	char		dummy[256]; 
	int			status; 


	/* Initialize Display */
	for(int i=0; i<m_renderNum; i++)
		status |= GzInitDisplay(m_pDisplay[i]); 
	
	/* 
	* Tokens associated with triangle vertex values 
	*/ 
	nameListTriangle[0] = GZ_POSITION; 
	nameListTriangle[1] = GZ_NORMAL; 
	nameListTriangle[2] = GZ_TEXTURE_INDEX;  

	// I/O File open
	FILE *infile;
	if( (infile  = fopen( INFILE , "r" )) == NULL )
	{
         AfxMessageBox( "The input file was not opened\n" );
		 return GZ_FAILURE;
	}

	FILE *outfile;
	if( (outfile  = fopen( OUTFILE , "wb" )) == NULL )
	{
         AfxMessageBox( "The output file was not opened\n" );
		 return GZ_FAILURE;
	}

	/* 
	* Walk through the list of triangles, set color 
	* and render each triangle 
	*/ 
	while( fscanf(infile, "%s", dummy) == 1) { 	/* read in tri word */
	    fscanf(infile, "%f %f %f %f %f %f %f %f", 
		&(vertexList[0][0]), &(vertexList[0][1]),  
		&(vertexList[0][2]), 
		&(normalList[0][0]), &(normalList[0][1]), 	
		&(normalList[0][2]), 
		&(uvList[0][0]), &(uvList[0][1]) ); 
	    fscanf(infile, "%f %f %f %f %f %f %f %f", 
		&(vertexList[1][0]), &(vertexList[1][1]), 	
		&(vertexList[1][2]), 
		&(normalList[1][0]), &(normalList[1][1]), 	
		&(normalList[1][2]), 
		&(uvList[1][0]), &(uvList[1][1]) ); 
	    fscanf(infile, "%f %f %f %f %f %f %f %f", 
		&(vertexList[2][0]), &(vertexList[2][1]), 	
		&(vertexList[2][2]), 
		&(normalList[2][0]), &(normalList[2][1]), 	
		&(normalList[2][2]), 
		&(uvList[2][0]), &(uvList[2][1]) ); 

	    /* 
	     * Set the value pointers to the first vertex of the 	
	     * triangle, then feed it to the renderer 
	     * NOTE: this sequence matches the nameList token sequence
	     */ 
	     valueListTriangle[0] = (GzPointer)vertexList; 
		 valueListTriangle[1] = (GzPointer)normalList; 
		 valueListTriangle[2] = (GzPointer)uvList; 
		 for(int i=0; i<m_renderNum; i++)
			GzPutTriangle(m_pRender[i], 3, nameListTriangle, valueListTriangle);
	}

	GzIntensity rFilter=0, gFilter=0, bFilter=0;
	GzDisplay *filterDisplay;
	GzNewDisplay(&filterDisplay, m_nWidth, m_nHeight);
	int i=0, j=0, k=0;
	for(i=0;i<m_nHeight;i++)
	{
		for(j=0; j<m_nWidth; j++)
		{
			for(k=0; k<m_renderNum; k++)
			{
				rFilter += AAFilter[k][2] * m_pDisplay[k]->fbuf[i*m_nWidth+j].red;
				gFilter += AAFilter[k][2] * m_pDisplay[k]->fbuf[i*m_nWidth+j].green;
				bFilter += AAFilter[k][2] * m_pDisplay[k]->fbuf[i*m_nWidth+j].blue;
			}
			filterDisplay->fbuf[i*m_nWidth+j].red = rFilter;
			filterDisplay->fbuf[i*m_nWidth+j].green = gFilter;
			filterDisplay->fbuf[i*m_nWidth+j].blue = bFilter;

			rFilter = 0; gFilter = 0; bFilter = 0;
		}
	}

	GzFlushDisplay2File(outfile, filterDisplay); 	/* write out or update display to file*/
	GzFlushDisplay2FrameBuffer(m_pFrameBuffer, filterDisplay);	// write out or update display to frame buffer

	/* 
	 * Close file
	 */ 

	if( fclose( infile ) )
      AfxMessageBox( "The input file was not closed\n" );

	if( fclose( outfile ) )
      AfxMessageBox( "The output file was not closed\n" );
 
	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS); 
}

int Application5::Clean()
{
	/* 
	 * Clean up and exit 
	 */ 
	int	status = 0; 

	for(int i=0;i<m_renderNum;i++)
	{
		status |= GzFreeRender(m_pRender[i]); 
		status |= GzFreeDisplay(m_pDisplay[i]);
	}
	status |= GzFreeTexture();
	
	if (status) 
		return(GZ_FAILURE); 
	else 
		return(GZ_SUCCESS);
}



