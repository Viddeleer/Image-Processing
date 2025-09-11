# Simple Magic Kernel Sharp 2021 Implementation
Simple C/C++ Implementation of the Magic Kernel Sharp 2021 resizing algorithm by John Costella ( https://johncostella.com/magic/ ). The algorithm provides sharper results than Bicubic / Lanczos 3 interpolation and less aliasing, but at the cost of computation time. This version is slower than the original code, but I needed a simple plug and play version. It doesn't have separate paths for upscaling and downscaling, but its results are almost identical to the imagemagick implementation (apart from small rounding errors). This algorithm works for both upscaling and downscaling, also when changing aspect ratio.

Before resizing, images are converted from sRGB to linear color space, and converted back after resizing.

The algorithm is sped up by multithreading. Other possible optimizations like using a lookup table or using integer arithmetic are not implemented (yet?)

```C++
Usage:

unsigned char *outputbuffer = ResizeMagicKernelSharp2021MultiThreaded(inputbuffer, inputwidth, inputheight, inputchannels, outputwidth, outputheight, outputchannels, nrofthreads);
```

where input/outputchannels=3 for RGB images and input/outputchannels=4 for RGBA images (not fully implemented yet).

