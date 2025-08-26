# Fast Magic Kernel Sharp 2021
Fast C/C++ Implementation of the Magic Kernel Sharp 2021 resizing algorithm by John Costella ( https://johncostella.com/magic/ ), used by Facebook and Instagram.

The original algorithm is sped up by applying 3 optimizations:
1. Using a lookup table
2. Converting most float operations to integer arithmetic
3. Using multithreading

The algorithm works for both upscaling and downscaling.

Usage:

unsigned char *outputbuffer = ResizeMagicKernelSharp2021MultiThreaded(inputbuffer, inputwidth, inputheight, inputchannels, outputwidth, outputheight, outputchannels, nrofthreads);


