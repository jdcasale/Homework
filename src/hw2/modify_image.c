#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"

// EXTRA CREDIT: apply_median_filter implemented at line 524
// DOUBLE EXTRA CREDIT: not attempted

/******************************** Resizing *****************************
  To resize we'll need some interpolation methods and a function to create
  a new image and fill it in with our interpolation methods.
************************************************************************/

float nn_interpolate(image im, float x, float y, int c) {
    /***********************************************************************
      This function performs nearest-neighbor interpolation on image "im"
      given a floating column value "x", row value "y" and integer channel "c",
      and returns the interpolated value.
    ************************************************************************/

    return get_pixel(im, (int) roundf(x), (int) roundf(y), c);
}

image nn_resize(image im, int w, int h) {
    // TODO Fill in (also fix the return line)
    /***********************************************************************
      This function uses nearest-neighbor interpolation on image "im" to a new
      image of size "w x h"
    ************************************************************************/
    image resized_im = make_image(w, h, im.c);

    // did some math on paper to get here
    float a_w = ((float) im.w) / ((float) w);
    float b_w = (a_w - 1) / 2;

    float a_h = ((float) im.h) / ((float) h);
    float b_h = (a_h - 1) / 2;

    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            for (int ci = 0; ci < im.c; ++ci) {
                int index = i * w + j + ci * w * h;
                float i_old = a_h * (float) i + b_h;
                float j_old = a_w * (float) j + b_w;
                resized_im.data[index] = nn_interpolate(im, j_old, i_old, ci);
            }
        }
    }
    return resized_im;
}

float bilinear_interpolate(image im, float x, float y, int c) {
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

    float A1 = d2 * d4;
    float A2 = d1 * d4;
    float A3 = d2 * d3;
    float A4 = d1 * d3;

    return upperLeft * A1 + upperRight * A2 + lowerLeft * A3 + lowerRight * A4;
}

image bilinear_resize(image im, int w, int h) {
    // TODO
    /***********************************************************************
      This function uses bilinear interpolation on image "im" to a_w new image
      of size "w x h". Algorithm is same as nearest-neighbor interpolation.
    ************************************************************************/
    image resized_im = make_image(w, h, im.c);

    // did some math on paper to get here
    float a_w = ((float) im.w) / ((float) w);
    float b_w = (a_w - 1) / 2;

    float a_h = ((float) im.h) / ((float) h);
    float b_h = (a_h - 1) / 2;

    for (int i = 0; i < h; ++i) {
        for (int j = 0; j < w; ++j) {
            for (int ci = 0; ci < im.c; ++ci) {
                int index = i * w + j + ci * w * h;
                float i_old = a_h * (float) i + b_h;
                float j_old = a_w * (float) j + b_w;
                resized_im.data[index] = bilinear_interpolate(im, j_old, i_old, ci);
            }
        }
    }
    return resized_im;
}


/********************** Filtering: Box filter ***************************
  We want to create a box filter. We will only use square box filters.
************************************************************************/

void l1_normalize(image im) {
    // TODO
    /***********************************************************************
      This function divides each value in image "im" by the sum of all the
      values in the image and modifies the image in place.
    ************************************************************************/
    int pixel_count = im.h * im.w * im.c;
    float sum_values = 0;
    for (int i = 0; i < pixel_count; ++i) {
        sum_values += im.data[i];
    }
    for (int i = 0; i < pixel_count; ++i) {
        im.data[i] /= sum_values;
    }
}

image make_box_filter(int w) {
    /***********************************************************************
      This function makes a square filter of size "w x w". Make an image of
      width = height = w and number of channels = 1, with all entries equal
      to 1. Then use "l1_normalize" to normalize your filter.
    ************************************************************************/
    image im = make_image(w, w, 1);
    int pixel_count = im.h * im.w * im.c;
    for (int i = 0; i < pixel_count; ++i) {
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
            sum += get_pixel(im, up + offset_j, left + offset_i, c_img) *
                   get_pixel(filter, offset_j, offset_i, c_filter);
        }
    }
    return sum;
}

image sum_channels(image im) {
    int channel_size = im.w * im.h;
    image summed = make_image(im.w, im.h, 1);
    for (int i = 0; i < channel_size; ++i) {
        float r_val = im.data[i];
        float g_val = im.data[i + channel_size];
        float b_val = im.data[i + channel_size * 2];
        float y_val = r_val + g_val + b_val;
        summed.data[i] = y_val;
    }
    return summed;
}

image convolve_image(image im, image filter, int preserve) {
    assert(im.c == filter.c || filter.c == 1);
    /***********************************************************************
      This function convolves the image "im" with the "filter". The value
      of preserve is 1 if the number of input image channels need to be 
      preserved. Check the detailed algorithm given in the README.  
    ************************************************************************/
    image convolved = make_image(im.w, im.h, im.c);
    int kernel_width = filter.w / 2;

    if (filter.c == im.c) {
        // From assignment: "If filter and im have the same number of channels then it's just a normal convolution.
        // We sum over spatial and channel dimensions and produce a 1 channel image. UNLESS
        // If preserve is set to 1 we should produce an image with the same number of channels as the input.
        //
        // Either way, we need to sum over all three channels first:
        for (int cx = 0; cx < im.c; ++cx) {
            for (int i = 0; i < im.h; ++i) {
                for (int j = 0; j < im.w; ++j) {
                    set_pixel(convolved, j, i, cx, convolve_internal(im, filter, kernel_width, i, j, cx, cx));
                }
            }
        }
    } else if (im.c != 1 && filter.c == 1) {
        // From assignment: If the filter only has one channel but im has multiple channels we want to apply the
        // filter to each of those channels. Then we either sum between channels or not depending on if preserve is set.
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
        // if preserve is set, we return the multichannel image.
        return convolved;
    } else {
        // otherwise, we sum over the channels and return a single-channel image
        image summed = sum_channels(convolved);
        free_image(convolved);
        return summed;
    }
}

image make_highpass_filter() {
    /***********************************************************************
      Create a 3x3 filter with highpass filter values using image.data[]
    ************************************************************************/
    image filter = make_image(3, 3, 1);
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

image make_sharpen_filter() {
    /***********************************************************************
      Create a 3x3 filter with sharpen filter values using image.data[]
    ************************************************************************/
    // sharpen filter is the same as the highpass filter except the center pixel is 5 instead of 4
    image filter = make_highpass_filter();
    set_pixel(filter, 1, 1, 0, 5);
    return filter;
}

image make_emboss_filter() {
    /***********************************************************************
      Create a 3x3 filter with emboss filter values using image.data[]
    ************************************************************************/
    image filter = make_image(3, 3, 1);
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
// Answer:
// I think the answer here generally depends on what we're trying to accomplish with the filter, but there are some
// common-sense rules we can apply:

// We DON'T need preserve if the information we're trying to capture can be represented by a single, scalar value for each pixel.
// Things like the gradient and magnitude calculations produced by the Sobel operator fall into this image.
// Examples:
// - Highpass, if we only care about flagging pixels where _any_ channel passes the filter
// - GX and GY filters within Sobel operator -- these calculate a scalar value based on all three channels, so it probably
// does not make sense to calculate this for each channel separately.

// We DO need preserve if the information we're trying to capture information about each channel individually, or if
// we're trying to make stylistic modifications to an image and produce a viewable image as a result.
// Examples:
// - Box filter
// - Sharpen
// - Emboss
// - Gaussian blur

// Question 2.3.2: Do we have to do any post-processing for the above filters? Which ones and why?
// Answer:
// Yes, we need to do some postprocessing for Highpass, Sharpen and Emboss filters. When convolved with an image,
// all three of the filters can result in pixel multipliers greater than 1, which means they can create results outside
// of the range of valid pixel values. In the assignment, we handle this by applying our clamp_pixel() function, which
// clamps all of these above-or-below-range pixel values to either MAX or MIN, depending on which way they overflow.
// Because the box filter is simply averaging together pixel values (each of which is within the valid range) it cannot
// produce an invalid overflow in the same way -- thus we need no postprocessing for the Box filter.
//
// We can verify this -- each of the tests for these filters applies clamp_pixel() to the result before making comparisons.
// If we comment out that line on any of the Highpass, Sharpen, or Emboss tests, the test fails, showing a pixel value
// that overflows out of the valid range.
// If we comment out the clamp_pixel() call in the test for the box filter (test_convolution()), the test still passes.

image make_gaussian_filter(float sigma) {
    /***********************************************************************
      sigma: a float number for the Gaussian.
      Create a Gaussian filter with the given sigma. Note that the kernel size 
      is the next highest odd integer from 6 x sigma. Return the Gaussian filter.
    ************************************************************************/
    float kernel_size_f = (6.0f * sigma);
    int kernel_size = (int) (kernel_size_f);
    // paranoid about rounding errors
    if (fabsf(roundf(kernel_size_f) - kernel_size_f) > 0.00001) {
        ++kernel_size;
    }
    if (kernel_size % 2 == 0) {
        kernel_size++;
    }

    if (sigma == 4) {
        assert(kernel_size == 25);
    }
    int start = -kernel_size / 2;
    int end = kernel_size / 2;
    image filter = make_image(kernel_size, kernel_size, 1);

    float denominator = sigma * sigma * 2.0f;

    // instead of adding/subtracting 0.5 and calculating distange to the center axis of the kernel, in effect adding
    // a mu term to the gaussian, shift the window to the left so that mu = 0 and we don't need to do extra math
    for (int x = start; x <= end; x++) {
        for (int y = start; y <= end; y++) {
            float x2 = (float) x;
            float y2 = (float) y;
            float val = (expf(-(x2 * x2 + y2 * y2) / denominator)) / (denominator * M_PI);
            set_pixel(filter, x + end, y + end, 0, val);
        }
    }

    l1_normalize(filter);
    return filter;
}

void assert_images_have_same_dimensions(image a, image b) {
    assert(a.h == b.h && a.w == b.w && a.c == b.c);
}

image add_image(image a, image b) {
    /***********************************************************************
      The input images a and image b have the same height, width, and channels.
      Sum the given two images and return the result, which should also have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
    assert_images_have_same_dimensions(a, b);

    image im = make_image(a.w, a.h, a.c);
    int pixel_count = im.h * im.w * im.c;
    for (int i = 0; i < pixel_count; ++i) {
        im.data[i] = a.data[i] + b.data[i];
    }
    return im;
}

image sub_image(image a, image b) {
    /***********************************************************************
      The input image a and image b have the same height, width, and channels.
      Subtract the given two images and return the result, which should have
      the same height, width, and channels as the inputs. Do necessary checks.
    ************************************************************************/
    assert_images_have_same_dimensions(a, b);

    image im = make_image(a.w, a.h, a.c);
    int pixel_count = im.h * im.w * im.c;
    for (int i = 0; i < pixel_count; ++i) {
        im.data[i] = a.data[i] - b.data[i];
    }
    return im;
}

image make_gx_filter() {
    /***********************************************************************
      Create a 3x3 Sobel Gx filter and return it
    ************************************************************************/
    image filter = make_image(3, 3, 1);
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

image make_gy_filter() {
    /***********************************************************************
      Create a 3x3 Sobel Gy filter and return it
    ************************************************************************/
    image filter = make_image(3, 3, 1);
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

void feature_normalize(image im) {
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
}

image *sobel_image(image im) {
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
    for (int i = 0; i < im.h * im.w; ++i) {
        magnitude.data[i] = sqrtf(gy.data[i] * gy.data[i] + gx.data[i] * gx.data[i]);
        direction.data[i] = atan2f(gy.data[i], gx.data[i]);
    }

    image *sobelimg = calloc(2, sizeof(image));
    sobelimg[0] = magnitude;
    sobelimg[1] = direction;
    free_image(gx);
    free_image(gy);

    return sobelimg;
}

image colorize_sobel(image im) {
//  /***********************************************************************
//    Create a colorized version of the edges in image "im" using the
//    algorithm described in the README.
//  ************************************************************************/

    image f = make_gaussian_filter(4);
    image blur = convolve_image(im, f, 1);

    clamp_image(blur);
    feature_normalize(blur);

    image *sobels = sobel_image(blur);

    // need 3-channel output
    image output = make_image(im.w, im.h, 3);
    feature_normalize(sobels[0]);
    feature_normalize(sobels[1]);
    for (int x = 0; x < im.w; ++x) {
        for (int y = 0; y < im.h; ++y) {
            set_pixel(output, x, y, 1, get_pixel(sobels[0], x, y, 0));
            set_pixel(output, x, y, 2, get_pixel(sobels[0], x, y, 0));
            set_pixel(output, x, y, 0, get_pixel(sobels[1], x, y, 0));
        }
    }

    hsv_to_rgb(output);
    return output;
}

// EXTRA CREDIT: Median filter


int cmpfunc(const void *a, const void *b) {
    return (*(int *) a - *(int *) b);
}

/**
 * This is an implementation of the algorithm described on Wikipedia
 */

image apply_median_filter(image im, int kernel_size) {
    // kernel size must be an odd, positive number
    assert(kernel_size % 2 && kernel_size > 0);

    image kernel = make_image(kernel_size, kernel_size, 1);
    image output = make_image(im.w, im.h, im.c);

    int edge = kernel_size / 2;

    // this window will grab pixels from outside the image, but we rely on clamping
    // (returning the value of the pixel at the edge) to give a full kernel anyway
    for (int ci = 0; ci < im.c; ++ci) {
        for (int x = 0; x < im.w; ++x) {
            for (int y = 0; y < im.h; ++y) {
                int counter = 0;
                for (int fx = 0; fx < kernel_size; ++fx) {
                    for (int fy = 0; fy < kernel_size; ++fy) {
                        kernel.data[counter] = get_pixel(im, x + fx - 1*edge,y + fy - 1*edge, ci);
                        ++counter;
                    }
                }
                qsort(kernel.data, kernel_size * kernel_size, sizeof(float), cmpfunc);
                set_pixel(output, x, y, ci, kernel.data[kernel_size * kernel_size / 2]);
            }
        }
    }
    free_image(kernel);
    return output;
}

// SUPER EXTRA CREDIT: Bilateral filter

//image apply_bilateral_filter(image im, float sigma1, float sigma2)
//{
//    return make_image(1,1,1);
//}
