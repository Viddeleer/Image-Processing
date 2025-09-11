# Simple Magic Kernel Sharp 2021 Implementation
Simple C/C++ Implementation of the Magic Kernel Sharp 2021 resizing algorithm by John Costella ( https://johncostella.com/magic/ ). The algorithm provides sharper results than Bicubic / Lanczos 3 interpolation and less aliasing, but at the cost of computation time. This version is slower than the original code, but I needed a simple plug and play version. It doesn't have separate paths for upscaling and downscaling, but its results are almost identical to the imagemagick implementation (apart from small rounding errors). This algorithm works for both upscaling and downscaling, also when changing aspect ratio.

Before resizing, images are converted from sRGB to linear color space, and converted back after resizing.

The algorithm is sped up by multithreading. Other possible optimizations like using a lookup table or using integer arithmetic are not implemented (yet?)

```C++
Usage:

unsigned char *outputbuffer = ResizeMagicKernelSharp2021MultiThreaded(inputbuffer, inputwidth, inputheight, inputchannels, outputwidth, outputheight, outputchannels, nrofthreads);
```

where input/outputchannels=3 for RGB images and input/outputchannels=4 for RGBA images (not fully implemented yet).

Input (9x9 pixels):<BR>
![Input](https://github.com/Viddeleer/Image-Processing/blob/main/MagicKernelSharp2021/images/mks-2021-input.png)

Output (resized to 900x900 pixels):<BR>
![Output](https://github.com/Viddeleer/Image-Processing/blob/main/MagicKernelSharp2021/images/mks-2021-output.PNG)

For comparison, the 900x900 output by ImageMagick:<BR>

magick mks-2021-input.png -colorspace RGB -filter MagicKernelSharp2021 -resize 900x900 -colorspace sRGB -dither none mks-2021-magick_output.png

![Output](https://github.com/Viddeleer/Image-Processing/blob/main/MagicKernelSharp2021/images/mks-2021-magick_output.png)
