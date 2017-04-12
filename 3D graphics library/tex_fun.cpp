/* Texture functions for cs580 GzLib	*/
#include    "stdafx.h" 
#include	"stdio.h"
#include	"Gz.h"
#include    "math.h"
#include    <fstream>
#include   <complex>
using namespace std;
GzColor	*image=NULL;
int xs, ys;
int flag = 0;
int reset = 1;
ofstream real_texture("real_texture.txt");
/* Image texture function */
int tex_fun(float u, float v, GzColor color)
{
	unsigned char		pixel[3];
	unsigned char     dummy;
	char  		foo[8];
	int   		i, j;
	FILE			*fd;

	if (reset) {          /* open and load texture file */
		fd = fopen("texture.ppm", "rb");
		if (fd == NULL) {
			fprintf(stderr, "texture file not found\n");
			exit(-1);
		}
		fscanf(fd, "%s %d %d %c", foo, &xs, &ys, &dummy);
		image = (GzColor*)malloc(sizeof(GzColor)*(xs + 1)*(ys + 1));
		if (image == NULL) {
			fprintf(stderr, "malloc for texture image failed\n");
			exit(-1);
		}

		for (i = 0; i < xs*ys; i++) {	/* create array of GzColor values */
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
	if (u <= 1&&u>=0 && v <=1&&v>=0)
	{
		u *= (xs-1);
		v *= (ys-1);
		int A[2], B[2], C[2], D[2];
		GzColor Ac, Bc, Cc, Dc;
		double intpart;
		float s, t;
		s = u-floor(u);
		A[0] = floor(u);
		D[0] = A[0];
		B[0] = A[0] + 1;
		C[0] = B[0];
		t = v-floor(v);
		A[1] = floor(v);
		B[1] = A[1];
		C[1] = A[1] + 1;
		D[1] = C[1];
		real_texture << s << " " << t << " "<< endl;
		Ac[0] = image[A[1] *xs+ A[0]][0];
		Ac[1] = image[A[1] *xs+ A[0]][1];
		Ac[2] = image[A[1] *xs+ A[0]][2];
		Bc[0] = image[B[1] *xs+ B[0]][0];
		Bc[1] = image[B[1] *xs+ B[0]][1];
		Bc[2] = image[B[1] *xs+ B[0]][2];
		Cc[0] = image[C[1] *xs+ C[0]][0];
		Cc[1] = image[C[1] *xs+ C[0]][1];
		Cc[2] = image[C[1] *xs+ C[0]][2];
		Dc[0] = image[D[1] *xs+ D[0]][0];
		Dc[1] = image[D[1] *xs+ D[0]][1];
		Dc[2] = image[D[1] *xs+ D[0]][2];
		color[0] = s*t*Cc[0] + (1 - s)*t*Dc[0] + s*(1 - t)*Bc[0] + (1 - s)*(1 - t)*Ac[0];
		color[1] = s*t*Cc[1] + (1 - s)*t*Dc[1] + s*(1 - t)*Bc[1] + (1 - s)*(1 - t)*Ac[1];
		color[2] = s*t*Cc[2] + (1 - s)*t*Dc[2] + s*(1 - t)*Bc[2] + (1 - s)*(1 - t)*Ac[2];


		
		for (int k = 0; k < 3; k++)
		{
			if (color[k] > 1)color[k] = 0;
			else if (color[k] < 0)color[k] = 1;
		}
	}
	return GZ_SUCCESS;
}
/* Procedural texture function */
float reverse(short color)
{
	float x;
	x = float(color << 4) / ((1 << 12) - 1);
	return x;
}

int ptex_fun(float u, float v, GzColor color)
{	
	int n = 100;
	int W = 16, H = 16;
	complex<float> C(-0.8, 0.156);

	complex<float> x((u*W-W/2)/(W/2), (v*H-H/2)/(H/2));
	float len;
	float z=0.0;
	float key;
	for (int i = 0; i < n; i++)
	{
		//complex<float> xx(x.real()*x.real() - x.imag()*x.imag(), x.real()*x.real() + x.imag()*x.imag());
		x = x*x + C;
		len = sqrt(norm(x));
		if (len > 2)
		{
			z = i*1.0;
			//real_texture << z << endl;
			break;
		}
	}	
	key =float( z / n);
	
	color[2] = key;
	color[0] = key;
	color[1] = key;
	//real_texture << color[0] << " " << color[1] << " " << color[2] << endl;;

	return GZ_SUCCESS;
}

/* Free texture memory */
int GzFreeTexture()
{
	if(image!=NULL)
		free(image);
	return GZ_SUCCESS;
}

