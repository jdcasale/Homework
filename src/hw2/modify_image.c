#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"

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
//    printf("scale_factor_w: %f, scale_factor_h: %f\n", scale_factor_w, scale_factor_h);
//    printf("h %d, w: %d, new_h: %d, new_w: %d\n", im.h, im.w, h, w);

    float a_w = ((float)im.w - 0.5f) / ((float)w - 0.5f);
    float b_w = (a_w - 1) / 2;

    float a_h = ((float)im.h - 0.5f) / ((float)h - 0.5f);
    float b_h = (a_h - 1) / 2;

    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            for (int ci = 0; ci < im.c; ++ci) {
                int index = i * w + j + ci * w *h;
                float i_old = a_h * (float)i + b_h;
                float j_old = a_w * (float)j + b_w;
//                resized_im.data[index] = nn_interpolate(im, ((float)j + 0.5f) * scale_factor_w ,((float)i + 0.5f) * scale_factor_h, ci);
                resized_im.data[index] = nn_interpolate(im, j_old ,i_old, ci);
            }
        }
    }
    return resized_im;
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
      This function uses bilinear interpolation on image "im" to a_w new image
      of size "w x h". Algorithm is same as nearest-neighbor interpolation.
    ************************************************************************/
    image resized_im = make_image(w,h,im.c);
//    float a_w, b_w;
//    float x, y;


    float a_w = ((float)im.w - 0.5f) / ((float)w - 0.5f);
    float b_w = (a_w - 1) / 2;

    float a_h = ((float)im.h - 0.5f) / ((float)h - 0.5f);
    float b_h = (a_h - 1) / 2;

//    a_w*(-0.5f) + b_w = -0.5f;
//    a_w*(w - 0.5f) + b_w =  im.w - 0.5f;
//    a_w*(h - 0.5f) + b_w =  im.h - 0.5f;

    float scale_factor_w = ((float)im.w) / ((float)w);
    float scale_factor_h = ((float)im.h) / ((float)h);
//    printf("scale_factor_w: %f, scale_factor_h: %f\n", scale_factor_w, scale_factor_h);
//    printf("h %d, w: %d, new_h: %d, new_w: %d\n", im.h, im.w, h, w);

    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            for (int ci = 0; ci < im.c; ++ci) {

                float i_old = a_h * (float)i + b_h;
                float j_old = a_w * (float)j + b_w;

                int index = j + i * w + ci * w *h;
//                resized_im.data[index] = bilinear_interpolate(im, (float)i * scale_factor_w, (float)j * scale_factor_h, ci);
                resized_im.data[index] = bilinear_interpolate(im, j_old, i_old, ci);
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

float convolve_internal(
        image im,
        image filter,
        int kernel_width,
        int i,
        int j,
        int c_img,
        int c_filter) {
    int left = i - kernel_width;
    int up = j - kernel_width;
    float sum = 0;
    for (int offset_i = 0; offset_i < filter.h; ++offset_i) {
        for (int offset_j = 0; offset_j < filter.w; ++offset_j) {
            sum += get_pixel(im, up + offset_j, left + offset_i, c_img) * get_pixel(filter, offset_j, offset_i, c_filter);
        }
    }
    return sum;
}

image sum_channels(image im)
{
    int channel_size = im.w * im.h;
    image summed = make_image(im.w, im.h, 1);
    for(int i = 0; i < channel_size; ++i) {
        float r_val = im.data[i];
        float g_val = im.data[i+ channel_size];
        float b_val = im.data[i+ channel_size*2];
        float y_val = r_val + g_val + b_val;
        summed.data[i] = y_val;
    }
    return summed;
}

image convolve_image(image im, image filter, int preserve)
{
    assert(im.c == filter.c || filter.c == 1);
    /***********************************************************************
      This function convolves the image "im" with the "filter". The value
      of preserve is 1 if the number of input image channels need to be 
      preserved. Check the detailed algorithm given in the README.  
    ************************************************************************/
    image convolved = make_image(im.w,im.h,im.c);
    int kernel_width = filter.w / 2;
    if (filter.c == im.c) {
        for (int cx = 0; cx < im.c; ++cx) {
            for (int i = 0; i < im.h; ++i) {
                for (int j = 0; j < im.w; ++j) {
                    set_pixel(convolved, j, i, cx, convolve_internal(im, filter, kernel_width, i, j, cx, cx));
                }
            }
        }
    } else if (im.c != 1 && filter.c == 1) {
        for (int cx = 0; cx < im.c; ++cx) {
            for (int i = 0; i < im.h; ++i) {
                for (int j = 0; j < im.w; ++j) {
                    set_pixel(convolved, j, i, cx, convolve_internal(im, filter, kernel_width, i, j, cx, 0));
                }
            }
        }
    } else {
        assert(0);
    }
    if (preserve) {
        return convolved;
    } else {
        image summed = sum_channels(convolved);
        free_image(convolved);
        return summed;
    }
}

image make_highpass_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with highpass filter values using image.data[]
    ************************************************************************/
    image filter = make_image(3,3,1);
    set_pixel(filter, 0, 0, 0, 0);
    set_pixel(filter, 2, 0, 0, 0);
    set_pixel(filter, 0, 2, 0, 0);
    set_pixel(filter, 2, 2, 0, 0);

    set_pixel(filter, 1, 0, 0, -1);
    set_pixel(filter, 0, 1, 0, -1);
    set_pixel(filter, 2, 1, 0, -1);
    set_pixel(filter, 1, 2, 0, -1);

    set_pixel(filter, 1, 1, 0, 4);

    return filter;
}

image make_sharpen_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with sharpen filter values using image.data[]
    ************************************************************************/
    image filter = make_highpass_filter();
    set_pixel(filter, 1,1,0,5);
    return filter;
}

image make_emboss_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 filter with emboss filter values using image.data[]
    ************************************************************************/
    image filter = make_image(3,3,1);
    set_pixel(filter, 0, 0, 0, -2);
    set_pixel(filter, 2, 0, 0, 0);
    set_pixel(filter, 0, 2, 0, 0);
    set_pixel(filter, 2, 2, 0, 2);

    set_pixel(filter, 1, 0, 0, -1);
    set_pixel(filter, 0, 1, 0, -1);
    set_pixel(filter, 2, 1, 0, 1);
    set_pixel(filter, 1, 2, 0, 1);

    set_pixel(filter, 1, 1, 0, 1);

    return filter;
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
    int kernal_size_f = (6.0f * sigma);
    int kernal_size = (int)(6.0f * sigma);
//    if ()
    if (kernal_size % 2 == 0) {
        kernal_size++;
    }
    int start = - kernal_size/2;
    int end = kernal_size/2;
    image filter = make_image(kernal_size,kernal_size,1);

    float denominator = sigma * sigma * 2.0f;
    for (int x = start; x <= end; x++) {
        for (int y = start; y <= end; y++) {
//            float x2 = (((float)kernal_size) / 2.0f) - (x + 0.5f);
//            float y2 = (((float)kernal_size) / 2.0f) - (y + 0.5f);
//            printf("kernel_size %d, x %d, y %d, x2 %f , y2 %f\n\n", kernal_size, x, y, x2, y2);
            float x2 = (float)x;
            float y2 = (float)y;
            float val = (expf(-(x2*x2 + y2*y2) / denominator)) / (denominator * M_PI);
            set_pixel(filter, x + end, y + end, 0, val);
        }
    }

    l1_normalize(filter);
    return filter;
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
    assert_images_have_same_dimensions(a, b);

    image im = make_image(a.w,a.h,a.c);
    int pixel_count = im.h * im.w * im.c;
    for(int i = 0; i < pixel_count; ++i) {
        im.data[i] = a.data[i] + b.data[i];
    }
    return im;
}

image sub_image(image a, image b)
{
    /***********************************************************************
      The input image a and image b have the same height, width, and channels.
      Subtract the given two images and return the result, which should have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
    assert_images_have_same_dimensions(a, b);

    image im = make_image(a.w,a.h,a.c);
    int pixel_count = im.h * im.w * im.c;
    for(int i = 0; i < pixel_count; ++i) {
        im.data[i] = a.data[i] - b.data[i];
    }
    return im;
}

image make_gx_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gx filter and return it
    ************************************************************************/
    image filter = make_image(3,3,1);
    set_pixel(filter, 0, 0, 0, -1);

    set_pixel(filter, 1, 0, 0, 0);
    set_pixel(filter, 2, 0, 0, 1);

    set_pixel(filter, 0, 1, 0, -2);
    set_pixel(filter, 0, 2, 0, -1);

    set_pixel(filter, 1, 1, 0, 0);
    set_pixel(filter, 1, 2, 0, 0);

    set_pixel(filter, 2, 1, 0, 2);
    set_pixel(filter, 2, 2, 0, 1);

    return filter;
}

image make_gy_filter()
{
    // TODO
    /***********************************************************************
      Create a 3x3 Sobel Gy filter and return it
    ************************************************************************/
    image filter = make_image(3,3,1);
    set_pixel(filter, 0, 0, 0, -1);

    set_pixel(filter, 1, 0, 0, -2);
    set_pixel(filter, 2, 0, 0, -1);

    set_pixel(filter, 0, 1, 0, 0);
    set_pixel(filter, 0, 2, 0, 1);

    set_pixel(filter, 1, 1, 0, 0);
    set_pixel(filter, 2, 2, 0, 1);

    set_pixel(filter, 1, 2, 0, 2);
    set_pixel(filter, 2, 1, 0, 0);

    return filter;
}

void feature_normalize(image im)
{
    // TODO
    /***********************************************************************
      Calculate minimum and maximum pixel values. Normalize the image by
      subtracting the minimum and dividing by the max-min difference.
    ************************************************************************/
    float max = -5000.0f;
    float min = 5000.0f;

    for (int i = 0; i < im.h * im.w * im.c; ++i) {
        if (im.data[i] < min) {
            min = im.data[i];
        }
        if (im.data[i] > max) {
            max = im.data[i];
        }
    }
    float range = max - min;
    for (int i = 0; i < im.h * im.w * im.c; ++i) {
        im.data[i] -= min;
        if (range) {
            im.data[i] /= range;
        }
    }
    printf("feature_normalize -- im.c: %d ", im.c);
}

image *sobel_image(image im)
{
    // TODO
    /***********************************************************************
      Apply Sobel filter to the input image "im", get the magnitude as sobelimg[0]
      and gradient as sobelimg[1], and return the result.
    ************************************************************************/
    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();
    image gx = convolve_image(im, gx_filter, 0);
    image gy = convolve_image(im, gy_filter, 0);
    image magnitude = make_image(im.w, im.h, 1);
    image direction = make_image(im.w, im.h, 1);
//    feature_normalize(gx);
//    feature_normalize(gy);
    for (int i = 0; i < im.h * im.w; ++i) {
        magnitude.data[i] = sqrtf(gy.data[i]*gy.data[i] + gx.data[i]*gx.data[i]);
        direction.data[i] = atanf(gy.data[i] / gx.data[i]);
    }

    image *sobelimg = calloc(2, sizeof(image));
    sobelimg[1] = magnitude;
//    feature_normalize(direction);
    sobelimg[0] = direction;
    free_image(gx);
    free_image(gy);

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
