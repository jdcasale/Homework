#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "image.h"
#include "matrix.h"
#include <time.h>

// Frees an array of descriptors.
// descriptor *d: the array.
// int n: number of elements in array.
void free_descriptors(descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        free(d[i].data);
    }
    free(d);
}

// Create a feature descriptor for an index in an image.
// image im: source image.
// int i: index in image for the pixel we want to describe.
// returns: descriptor for that index.
descriptor describe_index(image im, int i)
{
    int w = 5;
    descriptor d;
    d.p.x = i%im.w;
    d.p.y = i/im.w;
    d.data = calloc(w*w*im.c, sizeof(float));
    d.n = w*w*im.c;
    int c, dx, dy;
    int count = 0;
    // If you want you can experiment with other descriptors
    // This subtracts the central value from neighbors
    // to compensate some for exposure/lighting changes.
    for(c = 0; c < im.c; ++c){
        float cval = im.data[c*im.w*im.h + i];
        for(dx = -w/2; dx < (w+1)/2; ++dx){
            for(dy = -w/2; dy < (w+1)/2; ++dy){
                float val = get_pixel(im, i%im.w+dx, i/im.w+dy, c);
                d.data[count++] = cval - val;
            }
        }
    }
    return d;
}

// Marks the spot of a point in an image.
// image im: image to mark.
// ponit p: spot to mark in the image.
void mark_spot(image im, point p)
{
    int x = p.x;
    int y = p.y;
    int i;
    for(i = -9; i < 10; ++i){
        set_pixel(im, x+i, y, 0, 1);
        set_pixel(im, x, y+i, 0, 1);
        set_pixel(im, x+i, y, 1, 0);
        set_pixel(im, x, y+i, 1, 0);
        set_pixel(im, x+i, y, 2, 1);
        set_pixel(im, x, y+i, 2, 1);
    }
}

// Marks corners denoted by an array of descriptors.
// image im: image to mark.
// descriptor *d: corners in the image.
// int n: number of descriptors to mark.
void mark_corners(image im, descriptor *d, int n)
{
    int i;
    for(i = 0; i < n; ++i){
        mark_spot(im, d[i].p);
    }
}

// Creates a 1d Gaussian filter.
// float sigma: standard deviation of Gaussian.
// returns: single row image of the filter.
image make_1d_gaussian(float sigma)
{
    float edge_len = (6.0f * sigma);
    int edge_len_int = (int) (edge_len);
    // paranoid about rounding errors
    if (fabsf(roundf(edge_len) - edge_len) > 0.00001) {
        ++edge_len_int;
    }
    if (edge_len_int % 2 == 0) {
        ++edge_len_int;
    }

    image filter = make_image(edge_len_int,1,1);

    int half_edge = edge_len_int/2;
    float two_sigma_squared = sigma * sigma * 2.0f;

    for (int i = -1 * half_edge; i < half_edge; ++i) {
        float val = 1.0f/sqrtf(two_sigma_squared * M_PI)* expf(-powf((float)i, 2)/two_sigma_squared);
        set_pixel(filter, i, 0, 0, val);
    }
    return filter;
}

image reflect_img(image im) {
    image reflected = make_image(im.h, im.w, im.c);
    for (int x = 0; x < im.w; ++x) {
        for (int y = 0; y < im.h; ++y) {
            for (int c = 0; c < im.c; ++c) {
                set_pixel(reflected, y, x, c, get_pixel(im, x, y, c));
            }
        }
    }
    return reflected;
}

// Smooths an image using separable Gaussian filter.
// image im: image to smooth.
// float sigma: std dev. for Gaussian.
// returns: smoothed image.
image smooth_image(image im, float sigma)
{
    image filter = make_1d_gaussian(sigma);
    image convolved_once = convolve_image(im, filter, 1);
    image reflected_filter = reflect_img(filter);
    image convolved_twice = convolve_image(convolved_once, reflected_filter, 1);
    free_image(filter);
    free_image(reflected_filter);
    free_image(convolved_once);
    return convolved_twice;
}


void square_image_values(image im) {
    for(int x = 0; x < im.w; ++x) {
        for(int y = 0; y < im.h; ++y) {
            for(int c = 0; c < im.c; ++c) {
                set_pixel(im, x, y, c, powf(get_pixel(im, x, y, c), 2));
            }
        }
    }
}

image multiply_images_elementwise(image im, image im2) {
    assert(im.h == im2.h);
    assert(im.w == im2.w);
    assert(im.c == im2.c);
    image elementwise_product = make_image(im.w, im.h, im.c);
    for(int x = 0; x < im.w; ++x) {
        for(int y = 0; y < im.h; ++y) {
            for(int c = 0; c < im.c; ++c) {
                set_pixel(elementwise_product, x, y, c, get_pixel(im, x, y, c) * get_pixel(im2, x, y, c));
            }
        }
    }
    return elementwise_product;
}


// Calculate the structure matrix of an image.
// image im: the input image.
// float sigma: std dev. to use for weighted sum.
// returns: structure matrix. 1st channel is Ix^2, 2nd channel is Iy^2,
//          third channel is IxIy.
image structure_matrix(image im, float sigma)
{
    image S = make_image(im.w, im.h, 3);

    image gx_filter = make_gx_filter();
    image gy_filter = make_gy_filter();
    image gauss = make_gaussian_filter(sigma);
    image gauss2 = make_gaussian_filter(sigma);
    image gauss3 = make_gaussian_filter(sigma);

    image ix = convolve_image(im, gx_filter, 0);

    image iy = convolve_image(im, gy_filter, 0);
    image ix_iy = multiply_images_elementwise(ix, iy);

    square_image_values(ix);
    image ix_weighted = convolve_image(ix, gauss, 0);

    square_image_values(iy);
    image iy_weighted = convolve_image(iy, gauss2, 0);
    image ix_iy_weighted = convolve_image(ix_iy, gauss3, 0);

    int c = 0;
    for(int x = 0; x < im.w; ++x) {
        for(int y = 0; y < im.h; ++y) {
                set_pixel(S, x, y, c, get_pixel(ix_weighted, x, y, 0));
        }
    }
    ++c;
    for(int x = 0; x < im.w; ++x) {
        for(int y = 0; y < im.h; ++y) {
            set_pixel(S, x, y, c, get_pixel(iy_weighted, x, y, 0));
        }
    }
    ++c;
    for(int x = 0; x < im.w; ++x) {
        for(int y = 0; y < im.h; ++y) {
            set_pixel(S, x, y, c, get_pixel(ix_iy_weighted, x, y, 0));
        }
    }
    free_image(gx_filter);
    free_image(gy_filter);
    free_image(ix);
    free_image(iy);
    free_image(ix_weighted);
    free_image(iy_weighted);
    free_image(ix_iy);
    free_image(ix_iy_weighted);
    free_image(gauss);
    free_image(gauss2);
    free_image(gauss3);
    return S;
}

// Estimate the cornerness of each pixel given a structure matrix S.
// image S: structure matrix for an image.
// returns: a response map of cornerness calculations.
image cornerness_response(image S)
{
    image R = make_image(S.w, S.h, 1);
    // TODO: fill in R, "cornerness" for each pixel using the structure matrix.
    // We'll use formulation det(S) - alpha * trace(S)^2, alpha = .06.
    float alpha = 0.06f;
    for (int y = 0; y < S.h; ++y) {
        for (int x = 0; x < S.w; ++x) {
            float ad = get_pixel(S, x, y, 0) * get_pixel(S, x, y, 1);
            float bc = get_pixel(S, x, y, 3) * get_pixel(S, x, y, 3);
            float trace = get_pixel(S, x, y, 0) + get_pixel(S, x, y, 1);;
            float val = (ad - bc) - alpha * trace * trace;
            set_pixel(R, x, y, 0, val);
        }
    }
    return R;
}

// Perform non-max supression on an image of feature responses.
// image im: 1-channel image of feature responses.
// int w: distance to look for larger responses.
// returns: image with only local-maxima responses within w pixels.
image nms_image(image im, int w)
{
    image r = copy_image(im);
    // TODO: perform NMS on the response map.
    // for every pixel in the image:
    //     for neighbors within w:
    //         if neighbor response greater than pixel response:
    //             set response to be very low (I use -999999 [why not 0??])
    return r;
}

// Perform harris corner detection and extract features from the corners.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
// int *n: pointer to number of corners detected, should fill in.
// returns: array of descriptors of the corners in the image.
descriptor *harris_corner_detector(image im, float sigma, float thresh, int nms, int *n)
{
    // Calculate structure matrix
    image S = structure_matrix(im, sigma);

    // Estimate cornerness
    image R = cornerness_response(S);

    // Run NMS on the responses
    image Rnms = nms_image(R, nms);


    //TODO: count number of responses over threshold
    int count = 1; // change this

    
    *n = count; // <- set *n equal to number of corners in image.
    descriptor *d = calloc(count, sizeof(descriptor));
    //TODO: fill in array *d with descriptors of corners, use describe_index.


    free_image(S);
    free_image(R);
    free_image(Rnms);
    return d;
}

// Find and draw corners on an image.
// image im: input image.
// float sigma: std. dev for harris.
// float thresh: threshold for cornerness.
// int nms: distance to look for local-maxes in response map.
void detect_and_draw_corners(image im, float sigma, float thresh, int nms)
{
    int n = 0;
    descriptor *d = harris_corner_detector(im, sigma, thresh, nms, &n);
    mark_corners(im, d, n);
}
