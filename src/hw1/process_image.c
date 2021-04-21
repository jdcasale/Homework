#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "image.h"

// EXTRA CREDIT: scale_image IS IMPLEMENTED AT LINE 79
// DOUBLE EXTRA CREDIT: not attempted

float min(float a, float b) {
    if (a <= b) {
        return a;
    }
    return b;
}

float max(float a, float b) {
    if (a >= b) {
        return a;
    }
    return b;
}

float get_pixel(image im, int x, int y, int c)
{
//    assert (x <= im.w);
//    assert (y <= im.h);
//    assert (c <= im.c);
    int clamped_x = max(min(x, im.w-1), 0);
    int clamped_y = max(min(y, im.h-1), 0);
    int clamped_c = max(min(c, 2), 0);
    int channel_offset = clamped_c * im.w * im.h;
    int in_channel_offset = clamped_y * im.w + clamped_x;
    float val = im.data[channel_offset + in_channel_offset];
    return val;
}

void set_pixel(image im, int x, int y, int c, float v)
{
    int channel_offset = c * im.w * im.h;
    int in_channel_offset = y * im.w + x;
    im.data[channel_offset + in_channel_offset] = v;
}

image copy_image(image im)
{
    image copy = make_image(im.w, im.h, im.c);
    for(int i = 0; i < im.w * im.h * im.c; ++i) {
        copy.data[i] = im.data[i];
    }
    return copy;
}

image rgb_to_grayscale(image im)
{
    int channel_size = im.w * im.h;
    image grayscale = make_image(im.w, im.h, 1);
    for(int i = 0; i < channel_size; ++i) {
        float r_val = im.data[i];
        float g_val = im.data[i+ channel_size];
        float b_val = im.data[i+ channel_size*2];
        float y_val = 0.299 * r_val + 0.587 * g_val + .114 * b_val;
        grayscale.data[i] = y_val;
    }
    return grayscale;
}

void shift_image(image im, int c, float v)
{
    for(int i = im.w * im.h * c; i < im.w * im.h * (c+1); ++i) {
        im.data[i] += v;
    }
}

void clamp_image(image im)
{
    for(int i = 0; i < im.w * im.h * im.c; ++i) {
        float f = min(max(im.data[i], 0), 1);
        im.data[i] = f;
    }
}

void scale_image(image im, int c, float v) {
    for(int i = im.w * im.h * c; i < im.w * im.h * (c+1); ++i) {
        im.data[i] *= v;
    }
}


// These might be handy
float three_way_max(float a, float b, float c)
{
    return (a > b) ? ( (a > c) ? a : c) : ( (b > c) ? b : c) ;
}

float three_way_min(float a, float b, float c)
{
    return (a < b) ? ( (a < c) ? a : c) : ( (b < c) ? b : c) ;
}

float get_h_prime(float v_val, float c_val, float r_val, float g_val, float b_val) {
    if (v_val == r_val) { return (g_val - b_val) / c_val; }
    else if (v_val == g_val) { return (b_val - r_val) / c_val + 2.0f; }
    else if (v_val == b_val) { return (r_val - g_val) / c_val + 4.0f; }
    else {
        assert(0);
    }
}

void rgb_to_hsv(image im)
{
    int channel_size = im.w * im.h;
    for(int i = 0; i < channel_size; ++i) {
        int idx1 = i;
        int idx2 = i + channel_size;
        int idx3 = i + 2 * channel_size;

        float r_val = im.data[idx1];
        float g_val = im.data[idx2];
        float b_val = im.data[idx3];
        float val = three_way_max(r_val, g_val, b_val);
        float min_val = three_way_min(r_val, g_val, b_val);
        float c_val = val - min_val;

        if (c_val == 0) {
            im.data[idx1] = 0.0f;
            im.data[idx2] = 0.0f;
            im.data[idx3] = val;
        } else {
            float h_prime = get_h_prime(val, c_val, r_val, g_val, b_val);
            float h_val = h_prime < 0.0f ? (h_prime / 6.0f) + 1.0f : h_prime / 6.0f;
            im.data[idx1] = h_val;
            float sat = c_val / val;
            im.data[idx2] = sat;
            im.data[idx3] = val;
        }

    }
}

void hsv_to_rgb(image im)
{
    int channel_size = im.w * im.h;
    for(int i = 0; i < channel_size; ++i) {
        int idx1 = i;
        int idx2 = i + channel_size;
        int idx3 = i + 2 * channel_size;

        float h_val = im.data[idx1];
        float s_val = im.data[idx2];
        float v_val = im.data[idx3];

        float h_prime = h_val*6.0f;
        float hi = fmod(floor(h_prime), 6.0f);
        float f_val = h_prime - hi;
        int hi_int = (int) hi;

        float p_val = v_val * (1 - s_val);
        float q_val = v_val * (1- f_val * s_val);
        float t_val = v_val * (1 - (1-f_val) * s_val);
        switch(hi_int) {
            case 0  :
                im.data[idx1] = v_val;
                im.data[idx2] = t_val;
                im.data[idx3] = p_val;
                break;
            case 1  :
                im.data[idx1] = q_val;
                im.data[idx2] = v_val;
                im.data[idx3] = p_val;
                break;
            case 2  :
                im.data[idx1] = p_val;
                im.data[idx2] = v_val;
                im.data[idx3] = t_val;
                break;
            case 3  :
                im.data[idx1] = p_val;
                im.data[idx2] = q_val;
                im.data[idx3] = v_val;
                break;
            case 4  :
                im.data[idx1] = t_val;
                im.data[idx2] = p_val;
                im.data[idx3] = v_val;
                break;
            case 5  :
                im.data[idx1] = v_val;
                im.data[idx2] = p_val;
                im.data[idx3] = q_val;
                break;
            default:
                printf("invalid hi value: %d, %f",hi_int, hi);
                assert(0);
        }
    }
}