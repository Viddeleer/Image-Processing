#include <math.h>

#define MAXTHREADS 64						// maximum 64 threads limit in windows

static inline double magic_kernel_sharp_2021(double x)			// original
{
	if      (x < 0.0)  x = -x;	
	if      (x <= 0.5) return (577.0 / 576.0) - (239.0 / 144.0) * x * x;
	else if (x <= 1.5) return (1.0 / 144.0) * (140.0 * x * x - 379.0 * x + 239.0);
	else if (x <= 2.5) return -(1.0 / 144.0) * (24.0 * x * x - 113.0 * x + 130.0);
	else if (x <= 3.5) return (1.0 / 144.0) * (4.0 * x * x - 27.0 * x + 45.0);
	else if (x <= 4.5) return -(1.0 / 1152.0) * (2.0 * x - 9.0) * (2.0 * x - 9.0);
	else return 0.0;
}

struct ThreadParams
{
	float *src;
	float *dst;
	int src_w;
	int src_h;
	int src_c;
	int dst_w;
	int dst_h;
	int dst_c;
	int from;
	int to;	
} ThreadParameters[MAXTHREADS];			

void WINAPI ResizeImageMagicKernelSharp2021Thread(LPVOID lpParameters)
{
	float* dst;
	float* src;
	int src_w, src_h, src_c, dst_w, dst_h, dst_c, from, to;

	ThreadParams inputparams;
	memcpy(&inputparams, lpParameters, sizeof(ThreadParams));
	src = (float*)inputparams.src;
	dst = (float*)inputparams.dst;
	src_w = inputparams.src_w;
	src_h = inputparams.src_h;
	src_c = inputparams.src_c;
	dst_w = inputparams.dst_w;
	dst_h = inputparams.dst_h;
	dst_c = inputparams.dst_c;
	from = inputparams.from;
	to = inputparams.to;

	if (!src || !dst) return;
	if (src_w <= 0 || src_h <= 0 || dst_w <= 0 || dst_h <= 0) return;

	// continuous ratios
	const double fx_ratio = (double)src_w / (double)dst_w;
	const double fy_ratio = (double)src_h / (double)dst_h;

	// scale used to compute integer footprint radius (>= 1)
	const double scale_x = (fx_ratio > 1.0) ? fx_ratio : 1.0;
	const double scale_y = (fy_ratio > 1.0) ? fy_ratio : 1.0;

	int radius_x = 5;
	int radius_y = 5;

	if (dst_w < src_w) radius_x = (int)ceil(3.5 * scale_x);
	if (dst_h < src_h) radius_y = (int)ceil(3.5 * scale_y);
	
	for (int dy = from; dy < to; ++dy)
	{
		// map dst pixel center to source continuous coordinate:		
		double src_y_f = ((dy + 0.5) * (double)src_h / (double)dst_h) - 0.5;		
		int iy = (int)floor(src_y_f);
		double frac_y = src_y_f - (double)iy;   // in [0,1)

		for (int dx = 0; dx < dst_w; ++dx)
		{
			double src_x_f = ((dx + 0.5) * (double)src_w / (double)dst_w) - 0.5;			
			int ix = (int)floor(src_x_f);
			double frac_x = src_x_f - (double)ix;

			double sum_r = 0.0;
			double sum_g = 0.0;
			double sum_b = 0.0;
			double sum_a = 0.0;
			double wsum = 0.0;

			//int isedge = 0;

			// iterate over contributing source pixels
			for (int ky = -radius_y; ky <= radius_y; ++ky)
			{
				int sy = iy + ky;
				if ((unsigned)sy >= (unsigned)src_h) { /* isedge = 1; */ continue; }

				// distance in Y between sample pos and that row
				double dy_dist = fabs(frac_y - (double)ky);				
				if (dy_dist >= 4.5) continue;										
				double wy = magic_kernel_sharp_2021(dy_dist);

				for (int kx = -radius_x; kx <= radius_x; ++kx)
				{
					int sx = ix + kx;
					if ((unsigned)sx >= (unsigned)src_w) { /* isedge = 1; */ continue; }

					double dx_dist = fabs(frac_x - (double)kx);
					if (dx_dist >= 4.5) continue;									
					double wx = magic_kernel_sharp_2021(dx_dist);

					double w = wx * wy;
					if (w == 0.0) continue;

					int pixel_index = (sy * src_w + sx) * src_c;
					
					sum_r += src[pixel_index++] * w;
					sum_g += src[pixel_index++] * w;
					sum_b += src[pixel_index++] * w;
					if (src_c >= 4) sum_a += src[pixel_index] * w;
					wsum += w;
				}
			}

			int idx = (dy * dst_w + dx) * dst_c;

			//if (!isedge) wsum = 1.0;			// for upscaling no normalization is needed except for edge pixels, ideally put upscale and downscale in separate paths so division by wsum is not needed for upscaling

			if (wsum > 0.0)
			{								
					// normalize				
					dst[idx++] = (float)(sum_r / wsum);
					dst[idx++] = (float)(sum_g / wsum);
					dst[idx++] = (float)(sum_b / wsum);
					if (dst_c >= 4) dst[idx] = (float)(sum_a / wsum);			
			}
			else
			{
				// no contributors
				dst[idx++] = dst[idx++] = dst[idx++] = 0;
				if (dst_c >= 4) dst[idx] = 0;
			}
		}
	}
}

static inline float srgb_to_linear(float s)  // s in [0,1]
{
	if (s <= 0.04045f) return s / 12.92f;
	return powf((s + 0.055f) / 1.055f, 2.4f);
}

static inline float linear_to_srgb(float v)  // v in [0,1]
{
	if (v <= 0.0031308f) return v * 12.92f;
	return 1.055f * powf(v, 1.0f / 2.4f) - 0.055f;
}

unsigned char *ResizeImageMagicKernelSharp2021MultiThreaded(unsigned char *src, int src_w, int src_h, int src_c, int dst_w, int dst_h, int dst_c, int nrthreads)
{
int curthread = 0;
unsigned char *dst;
float* tmpin, * tmpout;
HANDLE threadhandles[MAXTHREADS];
unsigned i, insize, outsize;

	if (!src) return NULL;
	if (!src_w || !src_h || !src_c || !dst_w || !dst_h || !dst_c) return NULL;

	insize = src_w * src_h * src_c;
	outsize = dst_w * dst_h * dst_c;
	dst		= (unsigned char *)malloc(outsize);
	tmpin	= (float*)malloc(insize * sizeof(float));
	tmpout	= (float*)malloc(outsize * sizeof(float));

	if (!dst || !tmpin || !tmpout)
	{
		if (dst) free(dst);
		if (tmpin) free(tmpin);
		if (tmpout) free(tmpout);
		return NULL;
	}

	// populate the temporary float input buffer	
	for (i = 0; i<insize; i++) tmpin[i] = src[i];

	// if necessary, perform colorspace conversion from sRGB to linear RGB
	for (i = 0; i<insize; i++) tmpin[i] = 255.0f * srgb_to_linear(tmpin[i] / 255.0f);
	
	if (nrthreads < 2) nrthreads = 1;
	else if (nrthreads > MAXTHREADS) nrthreads = MAXTHREADS;	

	if (nrthreads == 1)	// Single thread
	{
		ThreadParameters[0].src = tmpin;
		ThreadParameters[0].dst = tmpout;
		ThreadParameters[0].src_w = src_w;
		ThreadParameters[0].src_h = src_h;
		ThreadParameters[0].src_c = src_c;
		ThreadParameters[0].dst_w = dst_w;
		ThreadParameters[0].dst_h = dst_h;
		ThreadParameters[0].dst_c = dst_c;
		ThreadParameters[0].from = (curthread * dst_h) ;
		ThreadParameters[0].to = ((curthread + 1) * dst_h);
		ThreadParameters[0].curthread = -1;
		ResizeImageMagicKernelSharp2021Thread((LPVOID)&ThreadParameters);
	}
	else				// Multithreaded version
	{
		for (curthread = 0; curthread < nrthreads; curthread++)
		{			
			ThreadParameters[curthread].src = tmpin;
			ThreadParameters[curthread].dst = tmpout;
			ThreadParameters[curthread].src_w = src_w;
			ThreadParameters[curthread].src_h = src_h;
			ThreadParameters[curthread].src_c = src_c;
			ThreadParameters[curthread].dst_w = dst_w;
			ThreadParameters[curthread].dst_h = dst_h;
			ThreadParameters[curthread].dst_c = dst_c;
			ThreadParameters[curthread].from = (curthread * dst_h) / nrthreads;
			ThreadParameters[curthread].to = ((curthread + 1) * dst_h) / nrthreads;
			ThreadParameters[curthread].curthread = curthread;
			do
			{
				threadhandles[curthread] = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)ResizeImageMagicKernelSharp2021Thread, (LPVOID)&ThreadParameters[curthread], 0, NULL);
				if (threadhandles[curthread] == NULL) Sleep(50);
			} while (threadhandles[curthread] == NULL);						// may hang infinitely if CreateThread keeps failing, consider using a limited number of attempts
		}
		WaitForMultipleObjects(nrthreads, threadhandles, TRUE, INFINITE);

		for (curthread = 0; curthread < nrthreads; curthread++) CloseHandle(threadhandles[curthread]);		
	}	

	// if necessary, perform colorspace conversion from linear to sRGB
	for (i = 0; i<outsize; i++) tmpout[i] = 255.0f * linear_to_srgb(tmpout[i] / 255.0f);

	// populate the final output buffer	
	float fval;
	for (i = 0; i < outsize; i++)
	{
		//fval = round(tmpout[i]);
		fval = tmpout[i];
		dst[i] = (fval < 0) ? 0 : ((fval > 255) ? 255 : fval);
	}
	return dst;

}





