# Simple C/C++ Magic Kernel Sharp 2013 Implementation
Simple C/C++ Implementation of the Magic Kernel Sharp 2013 resizing algorithm by John Costella ( https://johncostella.com/magic/ ). The algorithm provides sharper results than Bicubic / Lanczos 3 interpolation and less aliasing, but at the cost of computation time. This version is less efficient than the original reference code, but I needed a simple plug and play "drop-in" version. It doesn't have separate paths for upscaling and downscaling, but its results are (almost) identical to the ImageMagick implementation (except for some small rounding errors). This algorithm works for both upscaling and downscaling, also when changing aspect ratio.

Before resizing, images are converted from sRGB to linear color space, and converted back to sRGB after resizing.

The algorithm is sped up by multithreading. Other possible optimizations like using separate filter steps, a lookup table or using integer arithmetic are not implemented (yet?)

```C++
Usage:

unsigned char *outputbuffer = ResizeMagicKernelSharp2013MultiThreaded(inputbuffer, inputwidth, inputheight, inputchannels, outputwidth, outputheight, outputchannels, nrofthreads);

Example:

outputbuffer = ResizeMagicKernelSharp2021MultiThreaded(inputbuffer, 5, 5, 3, 500, 500, 3, 6);
```

where input/outputchannels=3 for RGB images and input/outputchannels=4 for RGBA images (not fully implemented yet).

**Input image (5x5 pixels):**<BR>
![Input](https://github.com/Viddeleer/Image-Processing/blob/main/MagicKernelSharp2013/images/mks-2013-input.png)

**Output image (resized to 500x500 pixels):**<BR>
![Output](https://github.com/Viddeleer/Image-Processing/blob/main/MagicKernelSharp2013/images/mks-2013-output.png)

**For comparison, the 500x500 output image by ImageMagick:**<BR>
```
magick mks-2013-input.png -colorspace RGB -filter MagicKernelSharp2013 -resize 500x500 -colorspace sRGB -dither none mks-2013-magick_output.png
```

![Output](https://github.com/Viddeleer/Image-Processing/blob/main/MagicKernelSharp2013/images/mks-2013-magick_output.png)
