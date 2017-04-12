/* CS580 Homework 3 */

#include	"stdafx.h"
#include	"stdio.h"
#include	"math.h"
#include	"Gz.h"
#include	"rend.h"
#include    <fstream>
#include    <algorithm>
#include	<cmath>
using namespace std;
//ofstream cmp("cmp.txt");
//ofstream Matrix("Matrix_stack.txt");
/* NOT part of API - just for general assistance */
extern int tex_fun(float u, float v, GzColor color); /* image texture function */
extern int ptex_fun(float u, float v, GzColor color); /* procedural texture function */
ofstream uv_debug("uv_debug.txt");
int num = 0;
short	ctoi(float color)		/* convert float color to GzIntensity short */
{
	return(short)((int)(color * ((1 << 12) - 1)));
}

int MyCompare(const void * a, const void * b)
{

	float *a1;
	a1 = ((float*)a + 1);
	float *a2;
	a2 = ((float*)b + 1);

	return (*a2 - *a1);

}



int GzRotXMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along x axis
// Pass back the matrix using mat value
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (i == 0 && j == 0)mat[i][j] = 1.0;
			else if (i == 1 && j == 1)mat[i][j] = cos(degree/180*3.1415926);
			else if (i == 1 && j == 2)mat[i][j] = -sin(degree/180*3.1415926);
			else if (i == 2 && j == 1)mat[i][j] = sin(degree/180*3.1415926);
			else if (i == 2 && j == 2)mat[i][j] = cos(degree/180*3.1415926);
			else if (i == 3 && j == 3)mat[i][j] = 1.0;
			else mat[i][j] = 0;
		}
	}


	return GZ_SUCCESS;
}


int GzRotYMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along y axis
// Pass back the matrix using mat value

	for(int i=0;i<4;i++)
		for (int j = 0; j < 4; j++)
		{
			if (i == 0 && j == 0)mat[i][j] = cos(degree / 180 * 3.1415926);
			else if (i == 0 && j == 2)mat[i][j] = sin(degree / 180 * 3.1415926);
			else if (i == 2 && j == 0)mat[i][j] = -sin(degree / 180 * 3.1415926);
			else if (i == 2 && j == 2)mat[i][j] = cos(degree / 180 * 3.1415926);
			else if (i == 3 && j == 3)mat[i][j] = 1.0;
			else if (i == 1 && j == 1)mat[i][j] = 1.0;
			else mat[i][j] = 0;
		}
	return GZ_SUCCESS;
}


int GzRotZMat(float degree, GzMatrix mat)
{
// Create rotate matrix : rotate along z axis
// Pass back the matrix using mat value
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (i == 0 && j == 0)mat[i][j] = cos(degree / 180 * 3.1415926);
			else if (i == 0 && j == 1)mat[i][j] = -sin(degree / 180 * 3.1415926);
			else if (i == 1 && j == 0)mat[i][j] = sin(degree / 180 * 3.1415926);
			else if (i == 1 && j == 1)mat[i][j] = cos(degree / 180 * 3.1415926);
			else if (i == 2 && j == 2)mat[i][j] = 1.0;
			else if (i == 3 && j == 3)mat[i][j] = 1.0;
			else mat[i][j] = 0;
		}
	}
	return GZ_SUCCESS;
}


int GzTrxMat(GzCoord translate, GzMatrix mat)
{
// Create translation matrix
// Pass back the matrix using mat value
	for(int i=0;i<4;i++)
		for (int j = 0; j < 4; j++)
		{
			if (i == j)mat[i][j] = 1.0;
			else if (j ==3 )
			{
				if (i == 0)mat[i][j] = translate[0];
				if (i == 1)mat[i][j] = translate[1];
				if (i == 2)mat[i][j] = translate[2];
				if (i == 3)mat[i][j] = 1;
			}
			else mat[i][j] = 0;

		}
	return GZ_SUCCESS;
}


int GzScaleMat(GzCoord scale, GzMatrix mat)
{
// Create scaling matrix
// Pass back the matrix using mat value
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (i == 0 && j == 0)mat[i][j] = scale[0];
			else if (i == 1 && j == 1)mat[i][j] = scale[1];
			else if (i == 2 && j == 2)mat[i][j] = scale[2];
			else if (i == 3 && j == 3)mat[i][j] = 1;
			else mat[i][j] = 0;
		}
	}
	return GZ_SUCCESS;
}


//----------------------------------------------------------
// Begin main functions

int GzNewRender(GzRender **render, GzDisplay	*display)
{
/*  
- malloc a renderer struct 
- setup Xsp and anything only done once 
- save the pointer to display 
- init default camera 
*/ 
	int xRes, yRes;
	GzGetDisplayParams(display, &xRes, &yRes);
	*render = (GzRender *)malloc(sizeof(GzRender));
	(*render)->display = display;
	GzMatrix Xsp;
	//ofstream Xform1("Xsp.txt");

	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (i == 0 && j == 0)Xsp[i][j] = float(xRes / 2);
			else if (i == 0 && j == 3)Xsp[i][j] =float(xRes / 2);
			else if (i == 1 && j == 1)Xsp[i][j] =float( -yRes / 2);
			else if (i == 1 && j == 3)Xsp[i][j] = float(yRes / 2);
			else if (i == 2 && j == 2)Xsp[i][j] = float(MAXINT);
			else if (i == 3 && j == 3)Xsp[i][j] = 1.0;
			else Xsp[i][j] = 0;
			//Xform1 << Xsp[i][j] << " ";
		}
		//Xform1 << endl;
	}
	for(int i=0;i<4;i++)
		for (int j = 0; j < 4; j++)
		{
			(*render)->Xsp[i][j] = Xsp[i][j];
		}
	(*render)->matlevel = -1;
	(*render)->numlights = 0;
	//defualt camera
	(*render)->camera.lookat[0] = 0;
	(*render)->camera.lookat[1] = 0;
	(*render)->camera.lookat[2] = 0;
	(*render)->camera.position[0] = DEFAULT_IM_X;
	(*render)->camera.position[1] = DEFAULT_IM_Y;
	(*render)->camera.position[2] = DEFAULT_IM_Z;
	(*render)->camera.worldup[0] = 0;
	(*render)->camera.worldup[1] = 1;
	(*render)->camera.worldup[2] = 0;
	(*render)->camera.FOV=DEFAULT_FOV;

	if (*render != NULL)
	{
		return GZ_SUCCESS;
	}
	else
	{
		return GZ_FAILURE;
	}

}


int GzFreeRender(GzRender *render)
{
/* 
-free all renderer resources
*/
	free(render);
	if (render == NULL)
	{
		return GZ_SUCCESS;
	}
	else
	{
		return GZ_FAILURE;
	}

}


int GzBeginRender(GzRender *render)
{
/*  
- setup for start of each frame - init frame buffer color,alpha,z
- compute Xiw and projection xform Xpi from camera definition 
- init Ximage - put Xsp at base of stack, push on Xpi and Xiw 
- now stack contains Xsw and app can push model Xforms when needed 
*/ 
	GzInitDisplay(render->display);
	GzMatrix Xiw, Xpi;
	double CL[3],UP[3],temp[3],Xxes[3];
	double scale,Pmulti;
	//ofstream Xform2("Xpi.txt");
	//ofstream Xform3("Xiw.txt");
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (i == 0 && j == 0)Xpi[i][j] = 1.0;
			else if (i == 1 && j == 1)Xpi[i][j] = 1.0;
			else if (i == 2 && j == 2)Xpi[i][j] = tan(((render->camera.FOV) / 2)/180*3.1415926);
			else if (i == 3 && j == 2)Xpi[i][j] = tan(((render->camera.FOV / 2))/180*3.1415926);
			else if (i == 3 && j == 3)Xpi[i][j] = 1.0;
			else Xpi[i][j] = 0;
			//Xform2 << Xpi[i][j] << " ";
		}
		//Xform2 << endl;
	}
	//Calculate Xiw,CL for Z, UP for Y, Xxes for X
	CL[X] = (render->camera.lookat[X] - render->camera.position[X]);
	CL[Y] = (render->camera.lookat[Y] - render->camera.position[Y]);
	CL[Z] = (render->camera.lookat[Z] - render->camera.position[Z]);
	//Xform3 <<"c[X]="<< render->camera.position[X] << " " <<"c[Y]="<< (render->camera.position[Y]) << " " <<"c[Z]="<< render->camera.position[Z] << endl;
	scale = (double)sqrt(CL[X] * CL[X] + CL[Y] * CL[Y] + CL[Z] * CL[Z]);
	//Xform3 << "scale=" << scale << endl;
	CL[X] /= scale;
	CL[Y] /= scale;
	CL[Z] /= scale;
	//Xform3 << "CL[X]=" << double(CL[X]) << " " << "CL[Y]=" << CL[Y] << " " << "CL[Z]=" << CL[Z] << endl;
	Pmulti = double(((render->camera.worldup[X]) * CL[X]) + ((render->camera.worldup[Y]) * CL[Y]) +( (render->camera.worldup[Z]) * CL[Z]));
	//Xform3 << "Pmulti=" <<Pmulti << endl;
	temp[X] = Pmulti*CL[X];
	temp[Y] = Pmulti*CL[Y];
	temp[Z] = Pmulti*CL[Z];
	UP[X] = render->camera.worldup[X] - temp[X];
	UP[Y] = render->camera.worldup[Y] - temp[Y];
	UP[Z] = render->camera.worldup[Z] - temp[Z];
	scale = (double)sqrt(UP[X] * UP[X] + UP[Y] * UP[Y] + UP[Z] * UP[Z]);
	UP[X] /= scale;
	UP[Y] /= scale;
	UP[Z] /= scale;

	Xxes[X] = CL[Z] * UP[Y] - UP[Z] * CL[Y];
	Xxes[Y] = UP[Z] * CL[X]	- CL[Z] * UP[X];
	Xxes[Z] = CL[Y] * UP[X] - UP[Y] * CL[X];
	scale = (double)sqrt(Xxes[X] * Xxes[X] + Xxes[Y] * Xxes[Y] + Xxes[Z] * Xxes[Z]);
	Xxes[X] /= scale;
	Xxes[Y] /= scale;
	Xxes[Z] /= scale;
	//Xform3 << endl;
	//Xform3 << Xxes[X] << " " << Xxes[Y] << " " << Xxes[Z] << endl;
	//Xform3 << endl;
	for (int i = 0; i < 4; i++)
	{
		for (int j = 0; j < 4; j++)
		{
			if (i == 0)
			{
				if (j == 0)Xiw[i][j] = Xxes[X];
				if (j == 1)Xiw[i][j] = Xxes[Y];
				if (j == 2)Xiw[i][j] = Xxes[Z];
				if (j == 3)Xiw[i][j] = ((-Xxes[X]) * render->camera.position[X]) + ((-Xxes[Y] )* render->camera.position[Y]) + ((-Xxes[Z]) * render->camera.position[Z]);

		}
			if (i == 1)
			{
				if (j == 0)Xiw[i][j] = UP[X];
				if (j == 1)Xiw[i][j] = UP[Y];
				if (j == 2)Xiw[i][j] = UP[Z];
				if (j == 3)Xiw[i][j] = (-UP[X] * render->camera.position[X]) + (-UP[Y] * render->camera.position[Y]) + (-UP[Z] * render->camera.position[Z]);

			}
			if (i == 2)
			{
				if (j == 0)Xiw[i][j] = CL[X];
				if (j == 1)Xiw[i][j] = CL[Y];
				if (j == 2)Xiw[i][j] = CL[Z];
				if (j == 3)Xiw[i][j] = (-CL[X] * render->camera.position[X]) + (-CL[Y] * render->camera.position[Y]) + (-CL[Z] * render->camera.position[Z]);
			}
			if (i == 3)
			{
				if (j == 3)Xiw[i][j] = 1.0;
				else Xiw[i][j] = 0;
			}
			//Xform3 << Xiw[i][j] << " ";
		}
		//Xform3 << endl;
	}
	for(int i=0;i<4;i++)
		for (int j = 0; j < 4; j++)
		{
			(render->camera.Xiw[i][j]) = Xiw[i][j];
			(render->camera.Xpi[i][j]) = Xpi[i][j];
		}
	GzPushMatrix(render, render->Xsp);
	GzPushMatrix(render, render->camera.Xpi);
	GzPushMatrix(render, render->camera.Xiw);

	return GZ_SUCCESS;
}

int GzPutCamera(GzRender *render, GzCamera *camera)
{
/*
- overwrite renderer camera structure with new camera definition
*/
	render->camera.FOV = camera->FOV;
	render->camera.lookat[X] = camera->lookat[X];
	render->camera.lookat[Y] = camera->lookat[Y];
	render->camera.lookat[Z] = camera->lookat[Z];
	render->camera.position[X] = camera->position[X];
	render->camera.position[Y] = camera->position[Y];
	render->camera.position[Z] = camera->position[Z];
	render->camera.worldup[X] = camera->worldup[X];
	render->camera.worldup[Y] = camera->worldup[Y];
	render->camera.worldup[Z] = camera->worldup[Z];
	

	return GZ_SUCCESS;	
}

int GzPushMatrix(GzRender *render, GzMatrix	matrix)
{
/*
- push a matrix onto the Ximage stack
- check for stack overflow
*/
	ofstream Matrix("Xnorm.txt");
	if (render->matlevel == MATLEVELS-1)return GZ_FAILURE;
	else if (render->matlevel == -1)
	{
		render->matlevel++;
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				(render->Ximage[render->matlevel])[j][k] = matrix[j][k];
				//Matrix << matrix[j][k] << " ";
			}
			//Matrix << endl;
		}
		
	}
	else
	{
		GzMatrix temp;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				float t=0.0;
				for (int k = 0; k < 4; k++)
				{
					t += (render->Ximage[render->matlevel])[i][k] * matrix[k][j];
				}
				temp[i][j] = t;
			}
		}

		render->matlevel++;
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				(render->Ximage[render->matlevel])[j][k] = temp[j][k];
				Matrix << temp[j][k] << " ";
			}
			Matrix << endl;
		}
		
	
	}

	Matrix << "Xnorm" << endl;
	if ((render->matlevel) == 2)
	{
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				(render->Xnorm[0])[j][k] = matrix[j][k];
				if (k == 3)(render->Xnorm[0])[j][k] = 0;
				
			}
		}
		(render->Xnorm[0])[3][3] = 1;
	}
	else if ((render->matlevel) > 2)
	{
		float scale;
		scale = sqrt((matrix[0][0])*(matrix[0][0]) + (matrix[0][1])*(matrix[0][1]) + (matrix[0][2])*(matrix[0][2]));
		Matrix << "scale=" << scale << endl;
		for(int i=0;i<4;i++)
			for (int j = 0; j < 4; j++)
			{
				matrix[i][j] /= scale;
				if (j == 3)matrix[i][j] = 0;
			}
		matrix[3][3] = 1;
		GzMatrix temp;
		for (int i = 0; i < 4; i++)
		{
			for (int j = 0; j < 4; j++)
			{
				float t = 0.0;
				for (int k = 0; k < 4; k++)
				{
					t += (render->Xnorm[(render->matlevel)-3])[i][k] * matrix[k][j];
				}
				temp[i][j] = t;
			}			
		}
		for (int j = 0; j < 4; j++)
		{
			for (int k = 0; k < 4; k++)
			{
				(render->Xnorm[render->matlevel-2])[j][k] = temp[j][k];
				Matrix << temp[j][k] << " ";
			}
			Matrix << endl;
		}

	}

	return GZ_SUCCESS;
}

int GzPopMatrix(GzRender *render)
{
/*
- pop a matrix off the Ximage stack
- check for stack underflow
*/
	if (render->matlevel == -1)return GZ_FAILURE;
	else render->matlevel--;
	
	return GZ_SUCCESS;
}


int GzPutAttribute(GzRender	*render, int numAttributes, GzToken	*nameList, 
	GzPointer	*valueList) /* void** valuelist */
{
/*
- set renderer attribute states (e.g.: GZ_RGB_COLOR default color)
- later set shaders, interpolaters, texture maps, and lights
*/
	for (int i = 0; i < numAttributes; i++)
	{

		if (nameList[i] == GZ_RGB_COLOR)
		{
			//GzColor a = (GzColor *)valueList[0];

			(render->flatcolor)[0] = (*((GzColor *)(valueList[i])))[0];
			(render->flatcolor)[1] = (*((GzColor *)(valueList[i])))[1];
			(render->flatcolor)[2] = (*((GzColor *)(valueList[i])))[2];

		}
		if (nameList[i] == GZ_DIRECTIONAL_LIGHT)
		{
			(render->lights)[(render->numlights)++] = *(GzLight *)(valueList[i]);
		}
		if (nameList[i] == GZ_AMBIENT_LIGHT)
		{
			render->ambientlight = *(GzLight *)(valueList[i]);
		}
		if (nameList[i] == GZ_DIFFUSE_COEFFICIENT)
		{

			(render->Kd)[0] = (*((GzColor *)(valueList[i])))[0];
			(render->Kd)[1] = (*((GzColor *)(valueList[i])))[1];
			(render->Kd)[2] = (*((GzColor *)(valueList[i])))[2];
		}
		if (nameList[i] == GZ_AMBIENT_COEFFICIENT)
		{
			(render->Ka)[0] = (*((GzColor *)(valueList[i])))[0];
			(render->Ka)[1] = (*((GzColor *)(valueList[i])))[1];
			(render->Ka)[2] = (*((GzColor *)(valueList[i])))[2];
		}
		if (nameList[i] == GZ_SPECULAR_COEFFICIENT)
		{
			(render->Ks)[0] = (*((GzColor *)(valueList[i])))[0];
			(render->Ks)[1] = (*((GzColor *)(valueList[i])))[1];
			(render->Ks)[2] = (*((GzColor *)(valueList[i])))[2];
		}
		if (nameList[i] == GZ_DISTRIBUTION_COEFFICIENT)
		{
			(render->spec) = *((float *)(valueList[i]));
		}
		if (nameList[i] == GZ_INTERPOLATE)
		{
			(render->interp_mode) = *((int *)(valueList[i]));
		}
		if (nameList[i] == GZ_TEXTURE_MAP)
		{
			(render->tex_fun) = (GzTexture)(valueList[i]);
		}
	}
	
	return GZ_SUCCESS;
}
float* GzGetColor(GzRender *render, GzCoord NormalList)
{
	
	GzColor KS;
	GzColor KD;
	GzColor RGB;
	GzCoord EYE = { 0,0,-1 };
	memset(RGB, 0, sizeof(RGB));
	memset(KS, 0, sizeof(KS));
	memset(KD, 0, sizeof(KD));
	for (int i = 0; i < render->numlights; i++)
	{
		float NL, RE, RN;
		GzCoord R;
		NL = NormalList[0] * ((render->lights)[i].direction[0])
			+ NormalList[1] * ((render->lights)[i].direction[1])
			+ NormalList[2] * ((render->lights)[i].direction[2]);
		RN = NormalList[0] * EYE[0] + NormalList[1] * EYE[1] + NormalList[2] * EYE[2];
		if (NL < 0 && RN < 0)
		{
			NormalList[0] *= -1;
			NormalList[1] *= -1;
			NormalList[2] *= -1;
			NL = NormalList[0] * ((render->lights)[i].direction[0])
				+ NormalList[1] * ((render->lights)[i].direction[1])
				+ NormalList[2] * ((render->lights)[i].direction[2]);
			RN = NormalList[0] * EYE[0] + NormalList[1] * EYE[1] + NormalList[2] * EYE[2];
		}
		if (NL*RN < 0)continue;
		R[0] = 2 * NL*NormalList[0] - ((render->lights)[i].direction)[0];
		R[1] = 2 * NL*NormalList[1] - ((render->lights)[i].direction)[1];
		R[2] = 2 * NL*NormalList[2] - ((render->lights)[i].direction)[2];
		RE = R[0] * EYE[0] + R[1] * EYE[1] + R[2] * EYE[2];
		if (RE < 0)	  RE *= -1.0;
		if (RE > 1.0)	RE = 1.0;
		if (NL > 1.0)	NL = 1.0;
		KS[0] += (((render->lights)[i].color)[0] * pow(RE, render->spec));
		KS[1] += (((render->lights)[i].color)[1] * pow(RE, render->spec));
		KS[2] += (((render->lights)[i].color)[2] * pow(RE, render->spec));
		KD[0] += (((render->lights)[i].color)[0] * NL);
		KD[1] += (((render->lights)[i].color)[1] * NL);
		KD[2] += (((render->lights)[i].color)[2] * NL);
		for (int k = 0; k < 3; k++)
		{
			if (KD[k] > 1)KD[k] = 1;
			if (KD[k] > 1)KS[k] = 1;
		}

	}
	RGB[0] = (render->Ks)[0] * KS[0] + (render->Kd)[0] * KD[0] + (render->Ka)[0] * (render->ambientlight).color[0];
	RGB[1] = (render->Ks)[1] * KS[1] + (render->Kd)[1] * KD[1] + (render->Ka)[1] * (render->ambientlight).color[1];
	RGB[2] = (render->Ks)[2] * KS[2] + (render->Kd)[2] * KD[2] + (render->Ka)[2] * (render->ambientlight).color[2];

	return RGB;
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
	//0.X(model->screen)

	float vertex_homo[3][4];
	float normal_homo[3][4];
	float uv_list[3][2];
	float temp[3][4];
	for (int k = 0; k < numParts; k++)
	{
		if (nameList[k] == GZ_TEXTURE_INDEX)
		{
	
			for (int i = 0; i < 3; i++)
			{				
				for (int j = 0; j < 2; j++)
				{
					uv_list[i][j] = *((float *)valueList[k] + i * 2 + j);
				}
				
			}
		}
		if (nameList[k] == GZ_POSITION)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					vertex_homo[i][j] = *((float *)valueList[k] + i * 3 + j);

				}
				vertex_homo[i][3] = 1.0;
			}
		}
		if (nameList[k] == GZ_NORMAL)
		{
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 3; j++)
				{
					normal_homo[i][j] = *((float *)valueList[k] + i * 3 + j);

				}
				normal_homo[i][3] = 1.0;			
			}
		}
	}
	for (int k = 0; k < 3; k++)
	{
		for (int i = 0; i < 4; i++)
		{
			temp[k][i] = 0.0;
			for (int j = 0; j < 4; j++)
			{
				temp[k][i] += (((render->Ximage[render->matlevel])[i][j])* (vertex_homo[k][j]));
				
			}

		}
		if (temp[k][2] <=0 )return GZ_SUCCESS;
	}

	float vertex[3][3];
	float ABC[3][3];
	float ABCD[4],R4[4],G4[4],B4[4],X4[4],Y4[4],Z4[4],U4[4],V4[4];
	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			vertex[i][j] = (temp[i][j] / temp[i][3]);
		}

	}
	//Transform uv_list to screen level for each vertex
	for (int i = 0; i < 3; i++)
	{
		double v=0.0;
		float x = vertex[i][2];
		v = (x) / (MAXINT - (x));
		//uv_debug << v << endl;
		uv_list[i][0] = uv_list[i][0] / (v + 1.0);
		uv_list[i][1] = uv_list[i][1] / (v + 1.0);
	}
	uv_debug << endl;
	//Transform normal_list to screen level
	for (int k = 0; k < 3; k++)
	{
		for (int i = 0; i < 4; i++)
		{
			temp[k][i] = 0.0;
			for (int j = 0; j < 4; j++)
			{
				temp[k][i] += (((render->Xnorm[(render->matlevel)-2])[i][j])* (normal_homo[k][j]));

			}

		}
	}

	for (int i = 0; i < 3; i++)
	{
		for (int j = 0; j < 3; j++)
		{
			normal_homo[i][j] = temp[i][j];
		}

	}
/*---------------------------------------Color------------------------------------------------*/
	GzColor ks_color[3];
	GzColor kd_color[3];
	GzColor vertex_color[3];
	if (render->interp_mode == GZ_COLOR)
	{
		//1.Calculate the color for each vertex
	

		memset(vertex_color, 0, sizeof(vertex_color));
		GzCoord E = { 0.0,0.0,-1.0 };
		ofstream color_debug("color_debug.txt");
		//color_debug << "numlight:" << render->numlights << endl;
		for (int j = 0; j < 3; j++)
		{
			memset(ks_color, 0, sizeof(vertex_color));
			memset(kd_color, 0, sizeof(vertex_color));
			for (int i = 0; i < render->numlights; i++)
			{
				GzCoord R;
				float coef;//N point L
				float coef2;//N point E
				float coef3;//R point E
				coef = normal_homo[j][0] * ((render->lights)[i].direction)[0] +
					normal_homo[j][1] * ((render->lights)[i].direction)[1] +
					normal_homo[j][2] * ((render->lights)[i].direction)[2];
				//color_debug << "coef=" << coef << endl;
				coef2 = normal_homo[j][0] * E[0] + normal_homo[j][1] * E[1] + normal_homo[j][2] * E[2];
				if (coef2 < 0 && coef < 0)
				{
					normal_homo[j][0] *= -1;
					normal_homo[j][1] *= -1;
					normal_homo[j][2] *= -1;
					coef = normal_homo[j][0] * ((render->lights)[i].direction)[0] +
						normal_homo[j][1] * ((render->lights)[i].direction)[1] +
						normal_homo[j][2] * ((render->lights)[i].direction)[2];
					coef2 = normal_homo[j][0] * E[0] + normal_homo[j][1] * E[1] + normal_homo[j][2] * E[2];
				}
				else if (coef*coef2 < 0)continue;
				//color_debug << "N*L=" << coef << endl;
				R[0] = 2 * coef*normal_homo[j][0] - ((render->lights)[i].direction)[0];
				R[1] = 2 * coef*normal_homo[j][1] - ((render->lights)[i].direction)[1];
				R[2] = 2 * coef*normal_homo[j][2] - ((render->lights)[i].direction)[2];
				//color_debug << "R={" << R[0] << "," << R[1] << "," << R[2] << "}" << endl;
				coef3 = R[0] * E[0] + R[1] * E[1] + R[2] * E[2];
				if (coef3 < 0) 	coef3 *= -1;
				if (coef3 > 1.0)	coef3 = 1.0;
				if (coef > 1.0)coef = 1.0;

				ks_color[j][0] += (((render->lights)[i].color)[0] * pow(coef3, render->spec));
				ks_color[j][1] += (((render->lights)[i].color)[1] * pow(coef3, render->spec));
				ks_color[j][2] += (((render->lights)[i].color)[2] * pow(coef3, render->spec));

				kd_color[j][0] += (((render->lights)[i].color)[0] * coef);
				kd_color[j][1] += (((render->lights)[i].color)[1] * coef);
				kd_color[j][2] += (((render->lights)[i].color)[2] * coef);

			}
			for (int k = 0; k < 3; k++)
			{
				if (ks_color[j][k] > 1)ks_color[j][k] = 0;
				if (kd_color[j][k] > 1)ks_color[j][k] = 0;
			}
	
			//with texture
			vertex_color[j][0] = ks_color[j][0] + kd_color[j][0] + (render->ambientlight).color[0];
			vertex_color[j][1] = ks_color[j][1] + kd_color[j][1] + (render->ambientlight).color[1];
			vertex_color[j][2] = ks_color[j][2] + kd_color[j][2] + (render->ambientlight).color[2];

		}

	}


/*---------------------------------Gouraud shade&& Flat shade-----------------------------------------------*/
	//1.Sorts verts by Y

	//qsort(vertex, 3, sizeof(vertex[0]), MyCompare);
	float change=0.0;
	float change_color = 0.0;
	float change_normal = 0.0;
	float change_uv = 0.0;
	if (vertex[0][1] < vertex[1][1])
	{
		for (int i = 0; i < 3; i++)
		{
			change = vertex[0][i];		
			change_color = vertex_color[0][i];
			change_normal = normal_homo[0][i];

			vertex[0][i] = vertex[1][i];
			vertex_color[0][i] = vertex_color[1][i];
			normal_homo[0][i] = normal_homo[1][i];

			vertex[1][i] = change;
			vertex_color[1][i] = change_color;
			normal_homo[1][i] = change_normal;

			if (i < 2)
			{
				change_uv = uv_list[0][i];
				uv_list[0][i] = uv_list[1][i];
				uv_list[1][i] = change_uv;
			}

		}
	}
	if (vertex[0][1] < vertex[2][1])
	{
		for (int i = 0; i < 3; i++)
		{
			change = vertex[0][i];
			change_color = vertex_color[0][i];
			change_normal = normal_homo[0][i];

			vertex[0][i] = vertex[2][i];
			vertex_color[0][i] = vertex_color[2][i];
			normal_homo[0][i] = normal_homo[2][i];

			vertex[2][i] = change;
			vertex_color[2][i] = change_color;
			normal_homo[2][i] = change_normal;

			if (i < 2)
			{
				change_uv = uv_list[0][i];
				uv_list[0][i] = uv_list[2][i];
				uv_list[2][i] = change_uv;
			}
		}
	}

	if (vertex[1][1] < vertex[2][1])
	{
		for (int i = 0; i < 3; i++)
		{
			change = vertex[1][i];
			change_color = vertex_color[1][i];
			change_normal = normal_homo[1][i];

			vertex[1][i] = vertex[2][i];
			vertex_color[1][i] = vertex_color[2][i];
			normal_homo[1][i] = normal_homo[2][i];

			vertex[2][i] = change;
			vertex_color[2][i] = change_color;
			normal_homo[2][i] = change_normal;

			if (i < 2)
			{
				change_uv = uv_list[1][i];
				uv_list[1][i] = uv_list[2][i];
				uv_list[2][i] = change_uv;
			}
		}
	}

	//2.calculate A,B,C
	ABC[0][0] = vertex[1][1] - vertex[0][1];
	ABC[0][1] = -vertex[1][0] + vertex[0][0];
	ABC[2][0] = vertex[2][1] - vertex[0][1];
	ABC[2][1] = -vertex[2][0] + vertex[0][0];
	ABC[1][0] = vertex[2][1] - vertex[1][1];
	ABC[1][1] = -vertex[2][0] + vertex[1][0];
	ABC[0][2] = -ABC[0][0] * vertex[0][0] - ABC[0][1] * vertex[0][1];
	ABC[1][2] = -ABC[1][0] * vertex[1][0] - ABC[1][1] * vertex[1][1];
	ABC[2][2] = -ABC[2][0] * vertex[2][0] - ABC[2][1] * vertex[2][1];

	//3.sort line
	int right = -1, left1 = -1, left2 = -1, left = -1, right1 = -1, right2 = -1;
	float x = (float)((-ABC[2][2] - ABC[2][1] * vertex[1][1]) / ABC[2][0]);

	if (x >  vertex[1][0])
	{
		right = 2;
		left1 = 0;
		left2 = 1;
	}
	else
	{
		left = 2;
		right1 = 0;
		right2 = 1;
	}

	//4.get Z plane

	ABCD[0] = (vertex[1][1] - vertex[0][1])*(vertex[2][2] - vertex[0][2]) - (vertex[2][1] - vertex[0][1])*(vertex[1][2] - vertex[0][2]);
	ABCD[1] = (vertex[1][2] - vertex[0][2])*(vertex[2][0] - vertex[0][0]) - (vertex[2][2] - vertex[0][2])*(vertex[1][0] - vertex[0][0]);
	ABCD[2] = (vertex[1][0] - vertex[0][0])*(vertex[2][1] - vertex[0][1]) - (vertex[2][0] - vertex[0][0])*(vertex[1][1] - vertex[0][1]);
	ABCD[3] = 0 - (ABCD[0] * vertex[0][0] + ABCD[1] * vertex[0][1] + ABCD[2] * vertex[0][2]);
	//4.1 get color
	//4.1.1 R
	R4[0] = (vertex[1][1] - vertex[0][1])*(vertex_color[2][0] - vertex_color[0][0]) - (vertex[2][1] - vertex[0][1])*(vertex_color[1][0] - vertex_color[0][0]);
	R4[1] = (vertex_color[1][0] - vertex_color[0][0])*(vertex[2][0] - vertex[0][0]) - (vertex_color[2][0] - vertex_color[0][0])*(vertex[1][0] - vertex[0][0]);
	R4[2] = (vertex[1][0] - vertex[0][0])*(vertex[2][1] - vertex[0][1]) - (vertex[2][0] - vertex[0][0])*(vertex[1][1] - vertex[0][1]);
	R4[3] = 0 - (R4[0] * vertex[0][0] + R4[1] * vertex[0][1] + R4[2] * vertex_color[0][0]);
	//4.1.2 G
	G4[0] = (vertex[1][1] - vertex[0][1])*(vertex_color[2][1] - vertex_color[0][1]) - (vertex[2][1] - vertex[0][1])*(vertex_color[1][1] - vertex_color[0][1]);
	G4[1] = (vertex_color[1][1] - vertex_color[0][1])*(vertex[2][0] - vertex[0][0]) - (vertex_color[2][1] - vertex_color[0][1])*(vertex[1][0] - vertex[0][0]);
	G4[2] = (vertex[1][0] - vertex[0][0])*(vertex[2][1] - vertex[0][1]) - (vertex[2][0] - vertex[0][0])*(vertex[1][1] - vertex[0][1]);
	G4[3] = 0 - (G4[0] * vertex[0][0] + G4[1] * vertex[0][1] + G4[2] * vertex_color[0][1]);
	//4.1.3 B	
	B4[0] = (vertex[1][1] - vertex[0][1])*(vertex_color[2][2] - vertex_color[0][2]) - (vertex[2][1] - vertex[0][1])*(vertex_color[1][2] - vertex_color[0][2]);
	B4[1] = (vertex_color[1][2] - vertex_color[0][2])*(vertex[2][0] - vertex[0][0]) - (vertex_color[2][2] - vertex_color[0][2])*(vertex[1][0] - vertex[0][0]);
	B4[2] = (vertex[1][0] - vertex[0][0])*(vertex[2][1] - vertex[0][1]) - (vertex[2][0] - vertex[0][0])*(vertex[1][1] - vertex[0][1]);
	B4[3] = 0 - (B4[0] * vertex[0][0] + B4[1] * vertex[0][1] + B4[2] * vertex_color[0][2]);
	//4.2 get normal
	//4.2.1 X
	X4[0] = (vertex[1][1] - vertex[0][1])*(normal_homo[2][0] - normal_homo[0][0]) - (vertex[2][1] - vertex[0][1])*(normal_homo[1][0] - normal_homo[0][0]);
	X4[1] = (normal_homo[1][0] - normal_homo[0][0])*(vertex[2][0] - vertex[0][0]) - (normal_homo[2][0] - normal_homo[0][0])*(vertex[1][0] - vertex[0][0]);
	X4[2] = (vertex[1][0] - vertex[0][0])*(vertex[2][1] - vertex[0][1]) - (vertex[2][0] - vertex[0][0])*(vertex[1][1] - vertex[0][1]);
	X4[3] = 0 - (X4[0] * vertex[0][0] + X4[1] * vertex[0][1] + X4[2] * normal_homo[0][0]);
	//4.2.2 Y
	Y4[0] = (vertex[1][1] - vertex[0][1])*(normal_homo[2][1] - normal_homo[0][1]) - (vertex[2][1] - vertex[0][1])*(normal_homo[1][1] - normal_homo[0][1]);
	Y4[1] = (normal_homo[1][1] - normal_homo[0][1])*(vertex[2][0] - vertex[0][0]) - (normal_homo[2][1] - normal_homo[0][1])*(vertex[1][0] - vertex[0][0]);
	Y4[2] = (vertex[1][0] - vertex[0][0])*(vertex[2][1] - vertex[0][1]) - (vertex[2][0] - vertex[0][0])*(vertex[1][1] - vertex[0][1]);
	Y4[3] = 0 - (Y4[0] * vertex[0][0] + Y4[1] * vertex[0][1] + Y4[2] * normal_homo[0][1]);
	//4.2.3 Z
	Z4[0] = (vertex[1][1] - vertex[0][1])*(normal_homo[2][2] - normal_homo[0][2]) - (vertex[2][1] - vertex[0][1])*(normal_homo[1][2] - normal_homo[0][2]);
	Z4[1] = (normal_homo[1][2] - normal_homo[0][2])*(vertex[2][0] - vertex[0][0]) - (normal_homo[2][2] - normal_homo[0][2])*(vertex[1][0] - vertex[0][0]);
	Z4[2] = (vertex[1][0] - vertex[0][0])*(vertex[2][1] - vertex[0][1]) - (vertex[2][0] - vertex[0][0])*(vertex[1][1] - vertex[0][1]);
	Z4[3] = 0 - (Z4[0] * vertex[0][0] + Z4[1] * vertex[0][1] + Z4[2] * normal_homo[0][2]);

	//4.3 get UV

	//4.3.1 get U
	U4[0] = (vertex[1][1] - vertex[0][1])*(uv_list[2][0] - uv_list[0][0]) - (vertex[2][1] - vertex[0][1])*(uv_list[1][0] - uv_list[0][0]);
	U4[1] = (uv_list[1][0] - uv_list[0][0])*(vertex[2][0] - vertex[0][0]) - (uv_list[2][0] - uv_list[0][0])*(vertex[1][0] - vertex[0][0]);
	U4[2] = (vertex[1][0] - vertex[0][0])*(vertex[2][1] - vertex[0][1]) - (vertex[2][0] - vertex[0][0])*(vertex[1][1] - vertex[0][1]);
	U4[3] = 0 - (U4[0] * vertex[0][0] + U4[1] * vertex[0][1] + U4[2] * uv_list[0][0]);	
	//4.3.2 get V
	V4[0] = (vertex[1][1] - vertex[0][1])*(uv_list[2][1] - uv_list[0][1]) - (vertex[2][1] - vertex[0][1])*(uv_list[1][1] - uv_list[0][1]);
	V4[1] = (uv_list[1][1] - uv_list[0][1])*(vertex[2][0] - vertex[0][0]) - (uv_list[2][1] - uv_list[0][1])*(vertex[1][0] - vertex[0][0]);
	V4[2] = (vertex[1][0] - vertex[0][0])*(vertex[2][1] - vertex[0][1]) - (vertex[2][0] - vertex[0][0])*(vertex[1][1] - vertex[0][1]);
	V4[3] = 0 - (V4[0] * vertex[0][0] + V4[1] * vertex[0][1] + V4[2] * uv_list[0][1]);

	for (int j = 0; j< 4; j++)
	{
		//uv_debug << U4[0] << " " << U4[1] << " " << U4[2] << " " << U4[3] << endl;
	}
	//5.start scan
	int  y = floor(vertex[0][1]);
	GzIntensity R = 0, G = 0, B = 0;
	float x1 = 0.0, x2 = 0.0;
	int xRes, yRes;
	GzGetDisplayParams(render->display, &xRes, &yRes);
	R = ctoi(((render->flatcolor)[0])) >> 4;
	G = ctoi(((render->flatcolor)[1])) >> 4;
	B = ctoi(((render->flatcolor)[2])) >> 4;

	if (left == 2)
	{
		x1 = float((0 - (ABC[left][2] + ABC[left][1] * y)) / ABC[left][0]);
		x2 = float((0 - (ABC[right1][2] + ABC[right1][1] * y)) / ABC[right1][0]);
	}
	else if (right == 2)
	{
		x1 = float((0 - (ABC[left1][2] + ABC[left1][1] * y)) / ABC[left1][0]);
		x2 = float((0 - (ABC[right][2] + ABC[right][1] * y)) / ABC[right][0]);
	}

	{
		while (y >= vertex[1][1])
		{

			for (int i = 0; i < yRes; i++)
				for (int j = 0; j < xRes; j++)
				{
					int x0 = 0, y0 = 0;
					float u0 = 0.0, v0 = 0.0;
					float v_p = 0.0;
					float z0 = MAXINT, z00;
					GzColor color0;

					x0 = ((render->display) + i*yRes + j)->xres;
					y0 = ((render->display) + i*yRes + j)->yres;
					z0 = (0 - (ABCD[0] * x0 + ABCD[1] * y + ABCD[3])) / ABCD[2];
					z00 = z0 / MAXINT;
					u0 = (0 - (U4[0] * x0 + U4[1] * y + U4[3])) / U4[2];
					v0 = (0 - (V4[0] * x0 + V4[1] * y + V4[3])) / V4[2];

					//retransform
					v_p = (z0) / (MAXINT - (z0 ));
					u0 *= (v_p + 1.0);
					v0 *= (v_p + 1.0);
					
					if (ceil(x1) <= floor(x2))
					{
						if ((x0 <= floor(x2)) && (x0 >= ceil(x1)) && (y0 == y))
						{
							if (z0 < (((render->display) + i*yRes + j)->fbuf->z))
							{
								(render->tex_fun)(u0, v0, color0);
								if (render->interp_mode == GZ_COLOR)
								{
									float RGB[3];
									RGB[0] = color0[0]*((0 - (R4[0] * x0 + R4[1] * y + R4[3])) / R4[2]);
									RGB[1] = color0[1]*((0 - (G4[0] * x0 + G4[1] * y + G4[3])) / G4[2]);
									RGB[2] = color0[2]*((0 - (B4[0] * x0 + B4[1] * y + B4[3])) / B4[2]);
									if (RGB[0] > 1)RGB[0] = 1.0;
									if (RGB[1] > 1)RGB[1] = 1.0;
									if (RGB[2] > 1)RGB[2] = 1.0;
									
									R = ctoi(RGB[0]) >> 4;
									G = ctoi(RGB[1]) >> 4;
									B = ctoi(RGB[2]) >> 4;


								}
								if (render->interp_mode == GZ_NORMALS)
								{
									GzCoord Normal;
									float *RGB;
									Normal[0]= ((0 - (X4[0] * x0 + X4[1] * y + X4[3])) / X4[2])*(v_p + 1.0);
									Normal[1]= ((0 - (Y4[0] * x0 + Y4[1] * y + Y4[3])) / Y4[2])*(v_p + 1.0);
									Normal[2] = ((0 - (Z4[0] * x0 + Z4[1] * y + Z4[3])) / Z4[2])*(v_p + 1.0);
									float scale;
									scale = sqrt((Normal[0] * Normal[0]) + (Normal[1] * Normal[1]) + (Normal[2] * Normal[2]));
									Normal[0] /= scale;
									Normal[1] /= scale;
									Normal[2] /= scale;

									render->Ka[0] = color0[0];
									render->Ka[1] = color0[1];
									render->Ka[2] = color0[2];

									render->Kd[0] = color0[0];
									render->Kd[1] = color0[1];
									render->Kd[2] = color0[2];

									RGB=GzGetColor(render, Normal);
									if (RGB[0] > 1)RGB[0] = 1.0;
									if (RGB[1] > 1)RGB[1] = 1.0;
									if (RGB[2] > 1)RGB[2] = 1.0;
									R = ctoi(RGB[0]) >> 4;
									G = ctoi(RGB[1]) >> 4;
									B = ctoi(RGB[2]) >> 4;

								}
								((render->display) + i*yRes + j)->fbuf->z = z0;
								(((render->display) + i*yRes + j)->fbuf)->red = R ;
								(((render->display) + i*yRes + j)->fbuf)->green = G;
								(((render->display) + i*yRes + j)->fbuf)->blue = B ;
							}

						}
						
					}
					
				}
			if (right == 2)
			{
				x2 += (float)(ABC[right][1]) / ABC[right][0];
				x1 += (float)(ABC[left1][1]) / ABC[left1][0];
			}
			else if (left == 2)
			{
				x2 += (float)(ABC[right1][1]) / ABC[right1][0];
				x1 += (float)(ABC[left][1]) / ABC[left][0];
			}
			y--;

		}
	}
	if (left == 2)
	{
		x1 = float((0 - (ABC[left][2] + ABC[left][1] * y)) / ABC[left][0]);
		x2 = float((0 - (ABC[right2][2] + ABC[right2][1] * y)) / ABC[right2][0]);
	}
	else if (right == 2)
	{
		x1 = float((0 - (ABC[left2][2] + ABC[left2][1] * y)) / ABC[left2][0]);
		x2 = float((0 - (ABC[right][2] + ABC[right][1] * y)) / ABC[right][0]);
	}

	{
		while (y >= vertex[2][1])
		{


			for (int i = 0; i < xRes; i++)
				for (int j = 0; j < yRes; j++)
				{
					int x0 = 0, y0 = 0;
					float u0 = 0.0, v0 = 0.0;
					float v_p = 0.0;
					float z0 = MAXINT;
					GzColor color0;

					x0 = ((render->display) + i*yRes + j)->xres;
					y0 = ((render->display) + i*yRes + j)->yres;
					z0 = (0 - (ABCD[0] * x0 + ABCD[1] * y + ABCD[3])) / ABCD[2];
					u0 = (0 - (U4[0] * x0 + U4[1] * y + U4[3])) / U4[2];
					v0 = (0 - (V4[0] * x0 + V4[1] * y + V4[3])) / V4[2];

					//retransform
					v_p = (z0) / (MAXINT - (z0));
					u0 *= (v_p + 1.0);
					v0 *= (v_p + 1.0);


					if (ceil(x1) <= floor(x2))
					{
						if ((x0 <= floor(x2)) && (x0 >= ceil(x1)) && (y0 == y))
						{
							if (z0 < ((render->display + i*yRes + j)->fbuf->z))
							{
								(render->tex_fun)(u0, v0, color0);
								uv_debug << v_p + 1.0 << endl;
								if (render->interp_mode == GZ_COLOR)
								{
									float RGB[3];
									RGB[0] = color0[0] * ((0 - (R4[0] * x0 + R4[1] * y + R4[3])) / R4[2]);
									RGB[1] = color0[1] * ((0 - (G4[0] * x0 + G4[1] * y + G4[3])) / G4[2]);
									RGB[2] = color0[2] * ((0 - (B4[0] * x0 + B4[1] * y + B4[3])) / B4[2]);
									if (RGB[0] > 1)RGB[0] = 1.0;
									if (RGB[1] > 1)RGB[1] = 1.0;
									if (RGB[2] > 1)RGB[2] = 1.0;

									R = ctoi(RGB[0]) >> 4;
									G = ctoi(RGB[1]) >> 4;
									B = ctoi(RGB[2]) >> 4;


								}
								if (render->interp_mode == GZ_NORMALS)
								{
									GzCoord Normal;
									float *RGB;
									Normal[0] = ((0 - (X4[0] * x0 + X4[1] * y + X4[3])) / X4[2])*(v_p + 1.0);
									Normal[1] = ((0 - (Y4[0] * x0 + Y4[1] * y + Y4[3])) / Y4[2])*(v_p + 1.0);
									Normal[2] = ((0 - (Z4[0] * x0 + Z4[1] * y + Z4[3])) / Z4[2])*(v_p + 1.0);
									float scale;
									scale = sqrt((Normal[0] * Normal[0]) + (Normal[1] * Normal[1]) + (Normal[2] * Normal[2]));
									Normal[0] /= scale;
									Normal[1] /= scale;
									Normal[2] /= scale;
									render->Kd[0] = color0[0];
									render->Kd[1] = color0[1];
									render->Kd[2] = color0[2];

									render->Ka[0] = color0[0];
									render->Ka[1] = color0[1];
									render->Ka[2] = color0[2];

									RGB=GzGetColor(render, Normal);
									if (RGB[0] > 1)RGB[0] = 1.0;
									if (RGB[1] > 1)RGB[1] = 1.0;
									if (RGB[2] > 1)RGB[2] = 1.0;
									R = ctoi(RGB[0]) >> 4;
									G = ctoi(RGB[1]) >> 4;
									B = ctoi(RGB[2]) >> 4;

								}
								(render->display + i*yRes + j)->fbuf->z = z0;
								(render->display + i*yRes + j)->fbuf->red = R ;
								(render->display + i*yRes + j)->fbuf->green = G ;
								(render->display + i*yRes + j)->fbuf->blue = B ;
							}

						}
					}
				}
			if (right == 2)
			{
				x2 += (float)((ABC[right][1]) / ABC[right][0]);
				x1 += (float)((ABC[left2][1]) / ABC[left2][0]);
			}
			else if (left == 2)
			{
				x2 += (float)((ABC[right2][1]) / ABC[right2][0]);
				x1 += (float)((ABC[left][1]) / ABC[left][0]);
			}
			y--;


		}
	}
	
	return GZ_SUCCESS;
}


