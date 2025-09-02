# Fast Magic Kernel Sharp 2021
Fast C/C++ Implementation of the Magic Kernel Sharp 2021 resizing algorithm by John Costella ( https://johncostella.com/magic/ ). The algorithm provides sharper results than Bicubic / Lanczos 3 interpolation and less aliasing, but at the cost of computation time. This version is optimized for speed to make it more usable.

The original algorithm is sped up by applying 3 optimizations:
1. Using a lookup table
2. Converting most float operations to integer arithmetic
3. Using multithreading

The algorithm works for both upscaling and downscaling.

```C++
Usage:

unsigned char *outputbuffer = ResizeMagicKernelSharp2021MultiThreaded(inputbuffer, inputwidth, inputheight, inputchannels, outputwidth, outputheight, outputchannels, nrofthreads);
```

where input/outputchannels=3 for RGB images and input/outputchannels=4 for RGBA images (not fully implemented yet, see source code comments).

*CAUTION* Work in progress, the calculations do not match the original exactly 
