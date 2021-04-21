#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#define TWOPI 6.2831853

/******************************** Resizing *****************************
  To resize we'll need some interpolation methods and a function to create
  a new image and fill it in with our interpolation methods.
************************************************************************/

float nn_interpolate(image im, float x, float y, int c)
{
    /***********************************************************************
      This function performs nearest-neighbor interpolation on image "im"
      given a floating column value "x", row value "y" and integer channel "c",
      and returns the interpolated value.
    ************************************************************************/

    return get_pixel(im, (int) roundf(x), (int) roundf(y), c);
}

image nn_resize(image im, int w, int h)
{
    // TODO Fill in (also fix the return line)
    /***********************************************************************
      This function uses nearest-neighbor interpolation on image "im" to a new
      image of size "w x h"
    ************************************************************************/
    image resized_im = make_image(w,h,im.c);
    float scale_factor_w = ((float)im.w) / ((float)w);
    float scale_factor_h = ((float)im.h) / ((float)h);
    printf("scale_factor_w: %f, scale_factor_h: %f\n", scale_factor_w, scale_factor_h);
    printf("h %d, w: %d, new_h: %d, new_w: %d\n", im.h, im.w, h, w);

    for (int i = 0; i < w; ++i) {
        for (int j = 0; j < h; ++j) {
            for (int ci = 0; ci < im.c; ++ci) {

                int index = i + j * h + ci * w *h;
//                if (index %5000 == 0) {
//                    printf("resized: h %d, w: %d, c: %d, index: %d\n", i, j, ci, index);
//                    printf("original: h %f, w: %f, c: %d\n", (float) i * scale_factor_w, (float) j * scale_factor_h,
//                           ci);
//                }
                resized_im.data[index] = nn_interpolate(im, (float)i * scale_factor_w, (float)j * scale_factor_h, ci);
            }
        }
    }
    return resized_im;
}

int index_to_coords(int i, int j, int k, image im) {
    int index = 0;
    index += im.w * im.h * k;
    index += i * im.w;
    index += j;
    return index;
}

int index_to_c_coord(int index, image im) {
    return index / (im.h * im.w);
}

int index_to_i_coord(int index, image im) {
    return (index - index_to_c_coord(index, im) * im.w * im.h) / im.w;
}

int index_to_j_coord(int index, image im) {
    index % im.w;
}

float bilinear_interpolate(image im, float x, float y, int c)
{
    // TODO
    /***********************************************************************
      This function performs bilinear interpolation on image "im" given
      a floating column value "x", row value "y" and integer channel "c".
      It interpolates and returns the interpolated value.
    ************************************************************************/

    float left = roundf(x - 0.5f);
    float right = roundf(x + 0.5f);
    float down = roundf(y + 0.5f);
    float up = roundf(y - 0.5f);

    float upperLeft = get_pixel(im, left, up, c);
    float upperRight = get_pixel(im, right, up, c);
    float lowerLeft = get_pixel(im, left, down, c);
    float lowerRight = get_pixel(im, right, down, c);

    float d1 = x - left;
    float d2 = right - x;
    float d3 = y - up;
    float d4 = down - y;

    float A1 = d2*d4;
    float A2 = d1*d4;
    float A3 = d2*d3;
    float A4 = d1*d3;

    return upperLeft*A1 + upperRight*A2 + lowerLeft*A3 + lowerRight*A4;
}

image bilinear_resize(image im, int w, int h)
{
    // TODO
    /***********************************************************************
      This function uses bilinear interpolation on image "im" to a new image
      of size "w x h". Algorithm is same as nearest-neighbor interpolation.
    ************************************************************************/
    image resized_im = make_image(w,h,im.c);
    float scale_factor_w = ((float)im.w) / ((float)w);
    float scale_factor_h = ((float)im.h) / ((float)h);
    printf("scale_factor_w: %f, scale_factor_h: %f\n", scale_factor_w, scale_factor_h);
    printf("h %d, w: %d, new_h: %d, new_w: %d\n", im.h, im.w, h, w);

    for (int i = 0; i < w; ++i) {
        for (int j = 0; j < h; ++j) {
            for (int ci = 0; ci < im.c; ++ci) {

                int index = i + j * h + ci * w *h;
//                if (index %5000 == 0) {
//                    printf("resized: h %d, w: %d, c: %d, index: %d\n", i, j, ci, index);
//                    printf("original: h %f, w: %f, c: %d\n", (float) i * scale_factor_w, (float) j * scale_factor_h,
//                           ci);
//                }
                resized_im.data[index] = bilinear_interpolate(im, (float)i * scale_factor_w, (float)j * scale_factor_h, ci);
            }
        }
    }
    return resized_im;
}


/********************** Filtering: Box filter ***************************
  We want to create a box filter. We will only use square box filters.
************************************************************************/

void l1_normalize(image im)
{
    // TODO
    /***********************************************************************
      This function divides each value in image "im" by the sum of all the
      values in the image and modifies the image in place.
    ************************************************************************/
    int pixel_count = im.h * im.w * im.c;
    float sum_values = 0;
    for(int i = 0; i < pixel_count; ++i) {
        sum_values += im.data[i];
    }
    for(int i = 0; i < pixel_count; ++i) {
        im.data[i] /= sum_values;
    }
}

image make_box_filter(int w)
{
    // TODO
    /***********************************************************************
      This function makes a square filter of size "w x w". Make an image of
      width = height = w and number of channels = 1, with all entries equal
      to 1. Then use "l1_normalize" to normalize your filter.
    ************************************************************************/
    image im = make_image(w,w,1);
    int pixel_count = im.h * im.w * im.c;
    for(int i = 0; i < pixel_count; ++i) {
        im.data[i] = 1;
    }
    l1_normalize(im);
    return im;
}

image convolve_image(image im, image filter, int preserve)
{
    // TODO
    /***********************************************************************
      This function convolves the image "im" with the "filter". The value
      of preserve is 1 if the number of input image channels need to be 
      preserved. Check the detailed algorithm given in the README.  
    ************************************************************************/
    return make_image(1,1,1);
}

image make_highpass_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with highpass filter values using image.data[]
    ************************************************************************/
    return make_image(1,1,1);
}

image make_sharpen_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with sharpen filter values using image.data[]
    ************************************************************************/
    return make_image(1,1,1);
}

image make_emboss_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with emboss filter values using image.data[]
    ************************************************************************/
    return make_image(1,1,1);
}

// Question 2.3.1: Which of these filters should we use preserve when we run our convolution and which ones should we not? Why?
// Answer: TODO

// Question 2.3.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer: TODO

image make_gaussian_filter(float sigma)
{
    // TODO
    /***********************************************************************
      sigma: a float number for the Gaussian.
      Create a Gaussian filter with the given sigma. Note that the kernel size 
      is the next highest odd integer from 6 x sigma. Return the Gaussian filter.
    ************************************************************************/
    return make_image(1,1,1);
}

void assert_images_have_same_dimensions(image a, image b) {
    assert(a.h == b.h && a.w == b.w && a.c == b.c);
}

image add_image(image a, image b)
{
    // TODO
    /***********************************************************************
      The input images a and image b have the same height, width, and channels.
      Sum the given two images and return the result, which should also have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
//    assert_images_have_same_dimensions(a, b);
    return make_image(1,1,1);

//    image im = make_image(a.w,a.h,a.c);
//    int pixel_count = im.h * im.w * im.c;
//    for(int i = 0; i < pixel_count; ++i) {
//        im.data[i] = a.data[i] + b.data[i];
//    }
//    return im;
}

image sub_image(image a, image b)
{
    /***********************************************************************
      The input image a and image b have the same height, width, and channels.
      Subtract the given two images and return the result, which should have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
//    assert_images_have_same_dimensions(a, b);
    return make_image(1,1,1);

//    image im = make_image(a.w,a.h,a.c);
//    int pixel_count = im.h * im.w * im.c;
//    for(int i = 0; i < pixel_count; ++i) {
//        im.data[i] = a.data[i] - b.data[i];
//    }
//    return im;
}

image make_gx_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gx filter and return it
    ************************************************************************/
    return make_image(1,1,1);
}

image make_gy_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gy filter and return it
    ************************************************************************/
    return make_image(1,1,1);
}

void feature_normalize(image im)
{
    // TODO
    /***********************************************************************
      Calculate minimum and maximum pixel values. Normalize the image by
      subtracting the minimum and dividing by the max-min difference.
    ************************************************************************/
}

image *sobel_image(image im)
{
    // TODO
    /***********************************************************************
      Apply Sobel filter to the input image "im", get the magnitude as sobelimg[0]
      and gradient as sobelimg[1], and return the result.
    ************************************************************************/
    image *sobelimg = calloc(2, sizeof(image));
    return sobelimg;
}

image colorize_sobel(image im)
{
  // TODO
  /***********************************************************************
    Create a colorized version of the edges in image "im" using the 
    algorithm described in the README.
  ************************************************************************/
  return make_image(1,1,1);
}

// EXTRA CREDIT: Median filter

/*
image apply_median_filter(image im, int kernel_size)
{
  return make_image(1,1,1);
}
*/

// SUPER EXTRA CREDIT: Bilateral filter

/*
image apply_bilateral_filter(image im, float sigma1, float sigma2)
{
  return make_image(1,1,1);
}
*/
