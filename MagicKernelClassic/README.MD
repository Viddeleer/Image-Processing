# Simple C/C++ Magic Kernel Classic Implementation

Implementation of the original Magic Kernel resizing algorithm as proposed by John Costella ( https://johncostella.com/magic/ ). This version is less efficient than the original reference code, but I needed a simple plug and play "drop-in" version.This algorithm works for both upscaling and downscaling, also when changing aspect ratio.

Before resizing, images are converted from sRGB to linear color space, and converted back to sRGB after resizing.

The algorithm is sped up by multithreading. Other possible optimizations like using separate filter steps, a lookup table or using integer arithmetic are not implemented (yet?)

```C++
Usage:

unsigned char *outputbuffer = ResizeMagicKernelClassicMultiThreaded(inputbuffer, inputwidth, inputheight, inputchannels, outputwidth, outputheight, outputchannels, nrofthreads);

Example:

outputbuffer = ResizeMagicKernelClassicMultiThreaded(inputbuffer, 5, 5, 3, 500, 500, 3, 6);
```

where input/outputchannels=3 for RGB images and input/outputchannels=4 for RGBA images (not fully implemented yet).



