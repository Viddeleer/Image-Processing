# Fast Magic Kernel Classic
Fast C/C++ Implementation of the Magic Kernel resizing algorithm by John Costella ( https://johncostella.com/magic/ ). The algorithm provides better results than Bicubic / Lanczos 3 interpolation with less aliasing.

This version use 3 optimizations:
1. Using a lookup table
2. Converting most float operations to integer arithmetic
3. Using multithreading

The algorithm works for both upscaling and downscaling.

```C++
Usage:

unsigned char *outputbuffer = ResizeMagicKernelClassicMultiThreaded(inputbuffer, inputwidth, inputheight, inputchannels, outputwidth, outputheight, outputchannels, nrofthreads);
```

where input/outputchannels=3 for RGB images and input/outputchannels=4 for RGBA images (not fully implemented yet, see source code comments).


