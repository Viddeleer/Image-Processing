#include <math.h>

#define MAXTHREADS 64						// maximum 64 threads limit in windows

#define MAGIC_KERNEL_RADIUS_2013 			2
#define MAGIC_KERNEL_FAST_PRECISION_BITS	18
#define MAGIC_KERNEL_FAST_PRECISION_FACTOR	(unsigned)(1<<MAGIC_KERNEL_FAST_PRECISION_BITS)
#define MAGIC_KERNEL_FAST_PRECISION_MASK	MAGIC_KERNEL_FAST_PRECISION_FACTOR-1

#define FIXED_MUL_FAST(a, b)				(((int64_t)(a) * (b)) >> MAGIC_KERNEL_FAST_PRECISION_BITS)
#define CLAMP_255_FAST(x)					(unsigned char)((x) <= 0 ? 0 : ((x) >= (255 << MAGIC_KERNEL_FAST_PRECISION_BITS) ? 255 : ((x) >> MAGIC_KERNEL_FAST_PRECISION_BITS)))
#define DIV_INT64_FAST(a, b)				(((int64_t)(a) << MAGIC_KERNEL_FAST_PRECISION_BITS) / (b))

#define MAGIC_KERNEL_LUT_RES_BITS_FAST		9
#define MAGIC_KERNEL_LUT_RES_FAST			(unsigned)(1<<MAGIC_KERNEL_LUT_RES_BITS_FAST)							// 512 for value 9
#define MAGIC_KERNEL_LUT_SIZE_2013_FAST		((int)(MAGIC_KERNEL_RADIUS_2013 * MAGIC_KERNEL_LUT_RES_FAST + 1))		// lookup table size

static int   magic_kernel_lut_2013_INT_FAST[MAGIC_KERNEL_LUT_SIZE_2013_FAST];

int MKS2013_INITIALIZED = 0;

// Calculate the kernel
static inline float magic_kernel_sharp_2013(float x) // magic_kernel_sharp_2013_official_implementation, , uses MAGIC_KERNEL_RADIUS = 2
{
	if (x < 0.0f) x = -x;
	if (x <= 1.0f) return 1.0f - (2.5f * x * x) + (1.5f * x * x * x);
	else if (x < 2.0f) return 0.5f * (2.0f - x) * (2.0f - x) * (1.0f - x);
	else return 0.0f;
}

// Initialize the lookup table
static void init_magic_kernel_lut_2013_INT_FAST()
{
	if (MKS2013_INITIALIZED) return;
	for (int i = 0; i < MAGIC_KERNEL_LUT_SIZE_2013_FAST; ++i) 
	{
		float x = (float)i / MAGIC_KERNEL_LUT_RES_FAST;
		magic_kernel_lut_2013_INT_FAST[i] = (int)round(MAGIC_KERNEL_FAST_PRECISION_FACTOR * magic_kernel_sharp_2013(x));						
	}
	MKS2013_INITIALIZED = 1;
}

struct ThreadParams
{
	unsigned char *src;
	unsigned char *dst;
	int src_w;
	int src_h;
	int src_c;
	int dst_w;
	int dst_h;
	int dst_c;
	int from;
	int to;	
} ThreadParameters[MAXTHREADS];			

void WINAPI ResizeImageMagicKernelSharp2013Thread(LPVOID lpParameters)
{
	int x_ratio, y_ratio;
	int dy, ix, iy, sx, sy, src_x, src_y, frac_x, frac_y, frac_xi, frac_yi;
	int dxf, dyf, lut_x, lut_y;
	int sum_r = 0, sum_g = 0, sum_b = 0, sum_a = 0, wsum = 0;
	int w, wx0, wx1, wx, wy0, wy1, wy;
	unsigned char* dst;	
	unsigned char* src;
	int src_w, src_h, src_c, dst_w, dst_h, dst_c, from, to;
		
	ThreadParams inputparams;
	memcpy(&inputparams, lpParameters, sizeof(ThreadParams));	
	src = inputparams.src;
	dst = inputparams.dst;
	src_w = inputparams.src_w;
	src_h = inputparams.src_h;
	src_c = inputparams.src_c;
	dst_w = inputparams.dst_w;
	dst_h = inputparams.dst_h;
	dst_c = inputparams.dst_c;
	from = inputparams.from;
	to = inputparams.to;
		
	x_ratio = (src_w << 16) / dst_w;
	y_ratio = (src_h << 16) / dst_h;

	int scale_x = (x_ratio > 65536) ? x_ratio : 65536;
	int scale_y = (y_ratio > 65536) ? y_ratio : 65536;

	int radius_x = 1 + ((MAGIC_KERNEL_RADIUS_2013 * scale_x) >> 16); // ceil
	int radius_y = 1 + ((MAGIC_KERNEL_RADIUS_2013 * scale_y) >> 16); // ceil

	for (dy = from; dy < to; ++dy)
	{
		src_y = FIXED_MUL((dy << 16) + 32768, y_ratio) - 32768;
		iy = src_y >> 16;
		frac_y = src_y & 0xFFFF;

		for (int dx = 0; dx < dst_w; ++dx)
		{
			src_x = FIXED_MUL((dx << 16) + 32768, x_ratio) - 32768;
			ix = src_x >> 16;
			frac_x = src_x & 0xFFFF;

			sum_r = 0, sum_g = 0, sum_b = 0; wsum = 0;

			for (int ky = -radius_y; ky <= radius_y; ++ky)
			{
				sy = iy + ky;
				if ((unsigned)sy >= (unsigned)src_h) continue;				

				dyf = frac_y - (ky << 16);
				if (dyf < 0) dyf = -dyf;
				lut_y = dyf >> (16 - MAGIC_KERNEL_LUT_RES_BITS_FAST);																		
				if (lut_y >= MAGIC_KERNEL_LUT_SIZE_2013_FAST) continue;

				// get interpolated LUT value
				frac_yi = dyf & MAGIC_KERNEL_FAST_PRECISION_MASK;							
				wy0 = magic_kernel_lut_2013_INT_FAST[lut_y];
				wy1 = magic_kernel_lut_2013_INT_FAST[lut_y + 1];
				wy = wy0 + FIXED_MUL_FAST(frac_yi, (wy1 - wy0));

				// slightly faster alternative without int64, only works with MAGIC_KERNEL_FAST_PRECISION_BITS <= 14
				//wy = wy0 + ((frac_yi * (wy1 - wy0)) >> MAGIC_KERNEL_FAST_PRECISION_BITS);				

				// fastest LUT alternative without interpolation, but with small rounding errors
				//wy = magic_kernel_lut_2013_INT_FAST[lut_y];

				// now compute weight w
				for (int kx = -radius_x; kx <= radius_x; ++kx)
				{
					sx = ix + kx;
					if ((unsigned)sx >= (unsigned)src_w) continue;

					dxf = frac_x - (kx << 16);
					if (dxf < 0) dxf = -dxf;
					lut_x = dxf >> (16 - MAGIC_KERNEL_LUT_RES_BITS_FAST);																	
					if (lut_x >= MAGIC_KERNEL_LUT_SIZE_2013_FAST) continue;

					// get interpolated LUT value
					frac_xi = dxf & MAGIC_KERNEL_FAST_PRECISION_MASK;
					wx0 = magic_kernel_lut_2013_INT_FAST[lut_x];
					wx1 = magic_kernel_lut_2013_INT_FAST[lut_x + 1];
					wx = wx0 + FIXED_MUL_FAST(frac_xi, (wx1 - wx0));

					// slightly faster alternative without int64, only works with MAGIC_KERNEL_FAST_PRECISION_BITS <= 14
					//wx = wx0 + ((frac_xi * (wx1 - wx0)) >> MAGIC_KERNEL_FAST_PRECISION_BITS);			

					// fastest LUT alternative without interpolation, but with small rounding errors
					//wx = magic_kernel_lut_2013_INT_FAST[lut_x];

					w = FIXED_MUL_FAST(wx, wy);
					//w = (wx * wy) >> MAGIC_KERNEL_FAST_PRECISION_BITS;

					int pixel_index = (sy * src_w + sx) * src_c;
					// currently only supports RGB, enable sum_a for RGBA input, and use src_c = 4
					sum_r += src[pixel_index + 0] * w;
					sum_g += src[pixel_index + 1] * w;
					sum_b += src[pixel_index + 2] * w;
					//sum_a += src[pixel_index + 3] * w;
					wsum += w;
				}
			}

			int idx = (dy * dst_w + dx) * dst_c;

			if (wsum)
			{
				// currently only supports RGB, enable sum_a for RGBA output, and use dst_c = 4
				dst[idx++] = CLAMP_255_FAST(DIV_INT64_FAST(sum_r, wsum));
				dst[idx++] = CLAMP_255_FAST(DIV_INT64_FAST(sum_g, wsum));
				dst[idx++] = CLAMP_255_FAST(DIV_INT64_FAST(sum_b, wsum));
				//dst[idx++] = CLAMP_255_FAST(DIV_INT64_FAST(sum_a, wsum));
			}
			else dst[idx++] = dst[idx++] = dst[idx++] = 0;
		}
	}	
}


unsigned char *ResizeImageMagicKernelSharp2013MultiThreaded(unsigned char *src, int src_w, int src_h, int src_c, int dst_w, int dst_h, int dst_c, int nrthreads)
{
int curthread = 0;
unsigned char *dst;
HANDLE threadhandles[MAXTHREADS];

	// input sanity checks
	if (!src) return NULL;
	if (!src_w || !src_h || !src_c || !dst_w || !dst_h || !dst_c) return NULL;
	dst = (unsigned char *)malloc(dst_w * dst_h * dst_c);
	if (!dst) return NULL;

	init_magic_kernel_lut_2013_INT_FAST();

	if (nrthreads < 2) nrthreads = 1;
	else if (nrthreads > MAXTHREADS) nrthreads = MAXTHREADS;
	
	if (nrthreads == 1)	// Single thread
	{
		ThreadParameters[0].src = src;
		ThreadParameters[0].dst = dst;
		ThreadParameters[0].src_w = src_w;
		ThreadParameters[0].src_h = src_h;
		ThreadParameters[0].src_c = src_c;
		ThreadParameters[0].dst_w = dst_w;
		ThreadParameters[0].dst_h = dst_h;
		ThreadParameters[0].dst_c = dst_c;
		ThreadParameters[0].from = (curthread * dst_h) ;
		ThreadParameters[0].to = ((curthread + 1) * dst_h);		
		ResizeImageMagicKernelSharp2013Thread((LPVOID)&ThreadParameters);
	}
	else				// Multithreaded version
	{
		for (curthread = 0; curthread < nrthreads; curthread++)
		{
			ThreadParameters[curthread].src = src;
			ThreadParameters[curthread].dst = dst;
			ThreadParameters[curthread].src_w = src_w;
			ThreadParameters[curthread].src_h = src_h;
			ThreadParameters[curthread].src_c = src_c;
			ThreadParameters[curthread].dst_w = dst_w;
			ThreadParameters[curthread].dst_h = dst_h;
			ThreadParameters[curthread].dst_c = dst_c;
			ThreadParameters[curthread].from = (curthread * dst_h) / nrthreads;
			ThreadParameters[curthread].to = ((curthread + 1) * dst_h) / nrthreads;			
			do
			{
				threadhandles[curthread] = CreateThread(NULL, 0, (LPTHREAD_START_ROUTINE)ResizeImageMagicKernelSharp2013Thread, (LPVOID)&ThreadParameters[curthread], 0, NULL);
				if (threadhandles[curthread] == NULL) Sleep(50);
			} while (threadhandles[curthread] == NULL);						// may hang infinitely if CreateThread keeps failing, consider using a limited number of attempts
		}
		WaitForMultipleObjects(nrthreads, threadhandles, TRUE, INFINITE);

		for (curthread = 0; curthread < nrthreads; curthread++) CloseHandle(threadhandles[curthread]);		
	}	
	return dst;
}
