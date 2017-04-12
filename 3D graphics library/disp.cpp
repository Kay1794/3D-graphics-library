/*   CS580 HW1 display functions to be completed   */

#include   "StdAfx.h"  
#include	"Gz.h"
#include	"disp.h"
#include	"stdio.h"

int GzNewFrameBuffer(char** framebuffer, int width, int height)
{
/* HW1.1 create a for MS Windows display:
 -- allocate memory for framebuffer : 3 bytes(b, g, r) x width x height
 -- pass back pointer 
 */	
	*framebuffer =(char *)malloc(3*sizeof(char)*width*height+100);//To devoid heap corruption;
	if (framebuffer != NULL)
	{
		return GZ_SUCCESS;
	}
	else return GZ_FAILURE;
}

int GzNewDisplay(GzDisplay	**display, int xRes, int yRes)
{
/* HW1.2 create a display:
  -- allocate memory for indicated resolution
  -- pass back pointer to GzDisplay object in display
*/
	*display = (GzDisplay *)malloc(sizeof(GzDisplay)*xRes*yRes);
	for (int a = 0; a < yRes; a++)
	{
		for (int b = 0; b < xRes; b++)
		{
			(*display + a*xRes + b)->xres = b;
			(*display + a*xRes + b)->yres = a;
			(*display + a*xRes + b)->fbuf = (GzPixel *)malloc(sizeof(GzPixel));
		}
	}
	if (display != NULL)
	{
		return GZ_SUCCESS;
	}
	else return GZ_FAILURE;
}


int GzFreeDisplay(GzDisplay	*display)
{
/* HW1.3 clean up, free memory */
	free(display);
	if (display == NULL)
	{
		return GZ_SUCCESS;
	}
	else return GZ_FAILURE;
}


int GzGetDisplayParams(GzDisplay *display, int *xRes, int *yRes)
{
/* HW1.4 pass back values for a display */
	int scale;
	scale = _msize(display) / sizeof(GzDisplay);
	*xRes = (display + scale-1)->xres+1;
	*yRes = (display + scale-1)->yres+1;
	return GZ_SUCCESS;
}


int GzInitDisplay(GzDisplay	*display)
{
/* HW1.5 set everything to some default values - start a new frame */

	int xRes, yRes;
	GzGetDisplayParams(display, &xRes, &yRes);
	for(int count=0;count<(xRes*yRes);count++)
	{
		display->fbuf = (GzPixel *)malloc(sizeof(GzPixel));
		display->fbuf->blue = 255;
		display->fbuf->red = 255;
		display->fbuf->green = 255;
		display->fbuf->alpha = 0;
		display->fbuf->z = MAXINT;
		display++;

	}

	return GZ_SUCCESS;
}


int GzPutDisplay(GzDisplay *display, int i, int j, GzIntensity r, GzIntensity g, GzIntensity b, GzIntensity a, GzDepth z)
{
	/* HW1.6 write pixel values into the display */
	

	int xRes, yRes;
	GzGetDisplayParams(display, &xRes, &yRes);
	
	// coordinates bounds check;
	if (i >= 0 && j >= 0 && i < xRes && j < yRes)
	{
			display = display +( j*xRes + i);
			// RGB bounds check

			//red
			if (r > 4095)display->fbuf->red = 255;
			else if (r < 0)display->fbuf->red = 0;
			else
			{
				display->fbuf->red =(GzIntensity)r >> 4;
			}
			//green
			if (g > 4095)display->fbuf->green = 255;
			else if (g < 0)display->fbuf->green = 0;
			else
			{
				display->fbuf->green = (GzIntensity)g >> 4;
			}
			//blue
			if (b > 4095)display->fbuf->blue = 255;
			else if (b < 0)display->fbuf->blue = 0;
			else
			{
				display->fbuf->blue = (GzIntensity)b >> 4;
			}	
			//aplpha & z
			display->fbuf->alpha = a;
			display->fbuf->z = z;
			
	}
	

	return GZ_SUCCESS;
}


int GzGetDisplay(GzDisplay *display, int i, int j, GzIntensity *r, GzIntensity *g, GzIntensity *b, GzIntensity *a, GzDepth *z)
{
/* HW1.7 pass back a pixel value to the display */

	int xRes, yRes;
	GzGetDisplayParams(display, &xRes, &yRes);
	if (i >= 0 && j >= 0 && i < xRes && j < yRes)
	{
		display = display + j*xRes + i;
		*r = display->fbuf->red;
		*g = display->fbuf->green;
		*b = display->fbuf->blue;
		*a = display->fbuf->alpha;
		*z = display->fbuf->z;
	}

	return GZ_SUCCESS;
}


int GzFlushDisplay2File(FILE* outfile, GzDisplay *display)
{

/* HW1.8 write pixels to ppm file -- "P6 %d %d 255\r" */
	int xRes, yRes;
	GzGetDisplayParams(display, &xRes, &yRes);
	fprintf(outfile, "P6 %d %d 255\r", xRes, yRes);
	for(int count=0;count<xRes*yRes;count++)
	{
		fprintf(outfile, "%c%c%c",display->fbuf->red,display->fbuf->green,display->fbuf->blue);
		display++;
	}
	return GZ_SUCCESS;
}

int GzFlushDisplay2FrameBuffer(char* framebuffer, GzDisplay *display)
{

	/* HW1.9 write pixels to framebuffer:
		- put the pixels into the frame buffer
		- CAUTION: when storing the pixels into the frame buffer, the order is blue, green, and red
		- NOT red, green, and blue !!!
	*/
	int xRes, yRes;
	GzGetDisplayParams(display, &xRes, &yRes);
	int pp = 0;
	for (int count = 0; count < xRes*yRes; count++)
	{
		framebuffer += pp;
		pp = sprintf(framebuffer, "%c%c%c", display->fbuf->blue, display->fbuf->green, display->fbuf->red);
		display++;
	}
	return GZ_SUCCESS;
}