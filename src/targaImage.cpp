///////////////////////////////////////////////////////////////////////////////
//
//      TargaImage.cpp                          Author:     Stephen Chenney
//                                              Modified:   Eric McDaniel
//                                              Date:       Fall 2004
//
//      Implementation of TargaImage methods.  You must implement the image
//  modification functions.
//
///////////////////////////////////////////////////////////////////////////////

#include "globals.h"
#include "targaImage.h"
#include "libtarga.h"
#include <stdlib.h>
#include <assert.h>
#include <memory.h>
#include <math.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <algorithm>

using namespace std;

// constants
const int           RED             = 0;                // red channel
const int           GREEN           = 1;                // green channel
const int           BLUE            = 2;                // blue channel
const unsigned char BACKGROUND[3]   = { 0, 0, 0 };      // background color


// Computes n choose s, efficiently
double Binomial(int n, int s)
{
    double        res;

    res = 1;
    for (int i = 1 ; i <= s ; i++)
        res = (n - i + 1) * res / i ;

    return res;
}// Binomial


///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage() : width(0), height(0), data(NULL)
{}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h) : width(w), height(h)
{
   data = new unsigned char[width * height * 4];
   ClearToBlack();
}// TargaImage



///////////////////////////////////////////////////////////////////////////////
//
//      Constructor.  Initialize member variables to values given.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(int w, int h, unsigned char *d)
{
    int i;

    width = w;
    height = h;
    data = new unsigned char[width * height * 4];

    for (i = 0; i < width * height * 4; i++)
	    data[i] = d[i];
}// TargaImage

///////////////////////////////////////////////////////////////////////////////
//
//      Copy Constructor.  Initialize member to that of input
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::TargaImage(const TargaImage& image) 
{
   width = image.width;
   height = image.height;
   data = NULL; 
   if (image.data != NULL) {
      data = new unsigned char[width * height * 4];
      memcpy(data, image.data, sizeof(unsigned char) * width * height * 4);
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Destructor.  Free image memory.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage::~TargaImage()
{
    if (data)
        delete[] data;
}// ~TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Converts an image to RGB form, and returns the rgb pixel data - 24 
//  bits per pixel. The returned space should be deleted when no longer 
//  required.
//
///////////////////////////////////////////////////////////////////////////////
unsigned char* TargaImage::To_RGB(void)
{
    unsigned char   *rgb = new unsigned char[width * height * 3];
    int		    i, j;

    if (! data)
	    return NULL;

    // Divide out the alpha
    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = i * width * 4;
	    int out_offset = i * width * 3;

	    for (j = 0 ; j < width ; j++)
        {
	        RGBA_To_RGB(data + (in_offset + j*4), rgb + (out_offset + j*3));
	    }
    }

    return rgb;
}// TargaImage


///////////////////////////////////////////////////////////////////////////////
//
//      Save the image to a targa file. Returns 1 on success, 0 on failure.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Save_Image(const char *filename)
{
    TargaImage	*out_image = Reverse_Rows();

    if (! out_image)
	    return false;

    if (!tga_write_raw(filename, width, height, out_image->data, TGA_TRUECOLOR_32))
    {
	    cout << "TGA Save Error: %s\n", tga_error_string(tga_get_last_error());
	    return false;
    }

    delete out_image;

    return true;
}// Save_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Load a targa image from a file.  Return a new TargaImage object which 
//  must be deleted by caller.  Return NULL on failure.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Load_Image(char *filename)
{
    unsigned char   *temp_data;
    TargaImage	    *temp_image;
    TargaImage	    *result;
    int		        width, height;

    if (!filename)
    {
        cout << "No filename given." << endl;
        return NULL;
    }// if

    temp_data = (unsigned char*)tga_load(filename, &width, &height, TGA_TRUECOLOR_32);
    if (!temp_data)
    {
        cout << "TGA Error: %s\n", tga_error_string(tga_get_last_error());
	    width = height = 0;
	    return NULL;
    }
    temp_image = new TargaImage(width, height, temp_data);
    free(temp_data);

    result = temp_image->Reverse_Rows();

    delete temp_image;

    return result;
}// Load_Image


///////////////////////////////////////////////////////////////////////////////
//
//      Convert image to grayscale.  Red, green, and blue channels should all 
//  contain grayscale value.  Alpha channel shoould be left unchanged.  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::To_Grayscale()
{
    if (data == NULL) {
        return false;
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;
            uchar gray_value = *(data + pos) * 0.299f + *(data + pos + 1) * 0.587f + *(data + pos + 2) * 0.114f;
            for (int k = 0; k < 3; k++) {
                *(data + pos + k) = gray_value;
            }
        }
    }

	return true;
}// To_Grayscale


///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using uniform quantization.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Uniform()
{
    if (data == NULL) {
        return false;
    }

    const int shades[] = { 5, 5, 6 };
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;
            for (int k = 0; k < 3; k++) {
                *(data + pos + k) = (*(data + pos + k) >> shades[k]) << shades[k];
            }
        }
    }

    return true;
}// Quant_Uniform


///////////////////////////////////////////////////////////////////////////////
//
//      Convert the image to an 8 bit image using populosity quantization.  
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Quant_Populosity()
{
    if (data == NULL) {
        return false;
    }

    // uniform quantization down to 32 shades each color
    // then build histogram
    int shades[] = { 3, 3, 3 };
    rgb_color histogram[32768];
    for (int i = 0; i < 32768; i++) {
        histogram[i].r = i % 32;
        histogram[i].g = (i % 1024) / 32;
        histogram[i].b = i / 1024;
        histogram[i].count = 0;
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;
            for (int k = 0; k < 3; k++) {
                *(data + pos + k) = *(data + pos + k) >> shades[k];
            }

            int index = *(data + pos) + *(data + pos + 1) * 32 + *(data + pos + 2) * 1024;
            histogram[index].count++;
        }
    }

    // sort histogram by count desc
    sort(histogram, histogram+32768, compare);

    // map to a color by index (32 * 32 * 32)
    rgb_color indexMap[32768];
    // most popular 256 colors
    for (int i = 0; i < 256; i++) {
        int index = histogram[i].r + histogram[i].g * 32 + histogram[i].b * 1024;
        indexMap[index].r = histogram[i].r << 3;
        indexMap[index].g = histogram[i].g << 3;
        indexMap[index].b = histogram[i].b << 3;
    }

    // find closest color for rest of colors
    for (int i = 256; i < 32768; i++) {
        int minIndex = -1;
        double minValue = 1e9f;
        for (int k = 0; k < 256; k++) {
            double distance = sqrt(pow(histogram[i].r - histogram[k].r, 2) + pow(histogram[i].g - histogram[k].g, 2) + pow(histogram[i].b - histogram[k].b, 2));
            if (distance < minValue) {
                minIndex = k;
                minValue = distance;
            }
        }
    
        int index = histogram[i].r + histogram[i].g * 32 + histogram[i].b * 1024;
        indexMap[index].r = histogram[minIndex].r << 3;
        indexMap[index].g = histogram[minIndex].g << 3;
        indexMap[index].b = histogram[minIndex].b << 3;
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;

            // replace color with the closest color
            int index = *(data + pos) + *(data + pos + 1) * 32 + *(data + pos + 2) * 1024;
            *(data + pos) = indexMap[index].r;
            *(data + pos + 1) = indexMap[index].g;
            *(data + pos + 2) = indexMap[index].b;
        }
    }

    return true;
}// Quant_Populosity


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image using a threshold of 1/2.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Threshold()
{
    if (data == NULL) {
        return false;
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;
            double alpha = (double)*(data + pos + 3);
            double value = (*(data + pos) / alpha) * 0.299 + (*(data + pos + 1) / alpha) * 0.587 + (*(data + pos + 2) / alpha) * 0.114;
            
            for (int k = 0; k < 3; k++) {
                *(data + pos + k) = (value >= 0.5) ? 255 : 0;
            }
        }
    }

    return true;
}// Dither_Threshold


///////////////////////////////////////////////////////////////////////////////
//
//      Dither image using random dithering.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Random()
{
    srand(time(NULL));
    if (data == NULL) {
        return false;
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;
            double alpha = (double)*(data + pos + 3);
            double value = (0.4 * rand() / RAND_MAX - 0.2) + (*(data + pos) / alpha) * 0.299 + (*(data + pos + 1) / alpha) * 0.587 + (*(data + pos + 2) / alpha) * 0.114;

            if (value > 1) {
                value = 1;
            } else if (value < 0) {
                value = 0;
            }

            // thresh using 0.5 threshold
            for (int k = 0; k < 3; k++) {
                *(data + pos + k) = (value >= 0.5) ? 255 : 0;
            }
        }
    }

    return true;
}// Dither_Random

///////////////////////////////////////////////////////////////////////////////
//
//      Chech if given point is inbound. Return success of assesstion
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::InBound(int h, int w) const
{
    if (h >= 0 && h < height && w >= 0 && w < width) {
        return true;
    }

    return false;
}

///////////////////////////////////////////////////////////////////////////////
//
//      Perform Floyd-Steinberg dithering on the image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_FS()
{
    if (data == NULL) {
        return false;
    }

    double image[height][width];
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;
            image[i][j] = (*(data + pos) * 0.299 + *(data + pos + 1) * 0.587 + *(data + pos + 2) * 0.114) / 255.0;
        }
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i & 1) ? (width - j - 1) : j;

            double old_value = image[i][pos];
            if (image[i][pos] > 0.5) {
                image[i][pos] = 1.0;
            } else {
                image[i][pos] = 0.0;
            }
            double error = old_value - image[i][pos];

            if (i & 1) {
                if (InBound(i, pos-1)) {
                    image[i][pos-1] += 0.4375 * error;
                }
                if (InBound(i+1, pos-1)) {
                    image[i+1][pos-1] += 0.0625 * error;
                }
                if (InBound(i+1, pos+1)) {
                    image[i+1][pos+1] += 0.1875 * error; 
                }
            } else {
                if (InBound(i, pos+1)) {
                    image[i][pos+1] += 0.4375 * error;
                }
                if (InBound(i+1, pos-1)) {
                    image[i+1][pos-1] += 0.1875 * error;
                }
                if (InBound(i+1, pos+1)) {
                    image[i+1][pos+1] += 0.0625 * error; 
                }
            }
            if (InBound(i+1, pos)) {
                image[i+1][pos] += 0.3125 * error;
            }
        }
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;
            
            for (int k = 0; k < 3; k++) {
                *(data + pos + k) = (uchar)(image[i][j] * 255.0);
            }
        }
    }

    return true;
}// Dither_FS


///////////////////////////////////////////////////////////////////////////////
//
//      Dither the image while conserving the average brightness.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Bright()
{
    if (data == NULL) {
        return false;
    }

    double gray_scale[] = { 0.299, 0.587, 0.114 };
    uint32_t sum = 0;
    int histogram[256] = { 0 };
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;
            double gray_value = 0;
            for (int k = 0; k < 3; k++) {
                gray_value += *(data + pos + k) * gray_scale[k];
            }

            sum += gray_value;
            histogram[(int)gray_value]++;
        }
    }

    double avg = 1 - ((double)sum / (255.0) / (double)(height * width));

    int count = 0, thresh = 0;
    do {
        count += histogram[thresh];
        thresh++;
    } while (count <= avg * (height * width));

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;
            uchar gray_value = 0;
            for (int k = 0; k < 3; k++) {
                gray_value += *(data + pos + k) * gray_scale[k];
            }

            for (int k = 0; k < 3; k++) {
                *(data + pos + k) = (gray_value >= thresh) ? 255 : 0;
            }
        }
    }

    return true;
}// Dither_Bright


///////////////////////////////////////////////////////////////////////////////
//
//      Perform clustered differing of the image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Cluster()
{
    if (data == NULL) {
        return false;
    }

    double mask[4][4] = {
        { 0.7059, 0.3529, 0.5882, 0.2353 },
        { 0.0588, 0.9412, 0.8235, 0.4118 },
        { 0.4706, 0.7647, 0.8824, 0.1176 },
        { 0.1765, 0.5294, 0.2941, 0.6471 }
    };

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;
            double alpha = (double)*(data + pos + 3);
            double value = (*(data + pos) / alpha) * 0.299 + (*(data + pos + 1) / alpha) * 0.587 + (*(data + pos + 2) / alpha) * 0.114;

            uchar result = 0;
            if (value >= mask[i%4][j%4]) {
                result = 255;
            }

            for (int k = 0; k < 3; k++) {
                *(data + pos + k) = result;
            }
        }
    }

    return true;
}// Dither_Cluster

uchar TargaImage::Trim(double v) const {
    if (v < 0) {
        return (uchar)0;
    } else if (v >= 256) {
        return (uchar)255;
    } else {
        return v;
    }
}

///////////////////////////////////////////////////////////////////////////////
//
//  Convert the image to an 8 bit image using Floyd-Steinberg dithering over
//  a uniform quantization - the same quantization as in Quant_Uniform.
//  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Dither_Color()
{
    if (data == NULL) {
        return false;
    }

    const int shades[] = { 8-1, 8-1, 4-1 };
    double temp_image[height][width][3];
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;
            
            for (int k = 0; k < 3; k++) {
                temp_image[i][j][k] = *(data + pos + k);
            }
        }
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t curr = (i & 1) ? (width - j - 1) : j;

            double old_values[3] = { 0.0 };
            double errors[3] = { 0.0 };
            for (int k = 0; k < 3; k++) {
                old_values[k] = temp_image[i][curr][k];
                temp_image[i][curr][k] = round(shades[k] * temp_image[i][curr][k] / 255.0) * (255.0 / shades[k]);
                errors[k] = old_values[k] - temp_image[i][curr][k];
            }

            if (i & 1) {
                if (InBound(i, curr-1)) {
                    for (int k = 0; k < 3; k++) {
                        temp_image[i][curr-1][k] += errors[k] * 0.4375;
                    }
                }
                if (InBound(i+1, curr-1)) {
                    for (int k = 0; k < 3; k++) {
                        temp_image[i+1][curr-1][k] += errors[k] * 0.0625;
                    }
                }
                if (InBound(i+1, curr+1)) {
                    for (int k = 0; k < 3; k++) {
                        temp_image[i+1][curr+1][k] += errors[k] * 0.1875;
                    }
                }
            } else {
                if (InBound(i, curr+1)) {
                    for (int k = 0; k < 3; k++) {
                        temp_image[i][curr+1][k] += errors[k] * 0.4375;
                    }
                }
                if (InBound(i+1, curr-1)) {
                    for (int k = 0; k < 3; k++) {
                        temp_image[i][curr-1][k] += errors[k] * 0.1875;
                    }
                }
                if (InBound(i+1, curr+1)) {
                    for (int k = 0; k < 3; k++) {
                        temp_image[i+1][curr+1][k] += errors[k] * 0.0625;
                    }
                }
            }
            if (InBound(i+1, curr)) {
                for (int k = 0; k < 3; k++) {
                    temp_image[i+1][curr][k] += errors[k] * 0.3125;
                }
            }
        }
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;
            
            for (int k = 0; k < 3; k++) {
                *(data + pos + k) = Trim(temp_image[i][j][k]);
            }
        }
    }

    return true;
}// Dither_Color


///////////////////////////////////////////////////////////////////////////////
//
//      Composite the current image over the given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Over(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout <<  "Comp_Over: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Over


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "in" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_In(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_In: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_In


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image "out" the given image.  See lecture notes for 
//  details.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Out(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Out: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Out


///////////////////////////////////////////////////////////////////////////////
//
//      Composite current image "atop" given image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Atop(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Atop: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Atop


///////////////////////////////////////////////////////////////////////////////
//
//      Composite this image with given image using exclusive or (XOR).  Return
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Comp_Xor(TargaImage* pImage)
{
    if (width != pImage->width || height != pImage->height)
    {
        cout << "Comp_Xor: Images not the same size\n";
        return false;
    }

    ClearToBlack();
    return false;
}// Comp_Xor


///////////////////////////////////////////////////////////////////////////////
//
//      Calculate the difference bewteen this imag and the given one.  Image 
//  dimensions must be equal.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Difference(TargaImage* pImage)
{
    if (!pImage)
        return false;

    if (width != pImage->width || height != pImage->height)
    {
        cout << "Difference: Images not the same size\n";
        return false;
    }// if

    for (int i = 0 ; i < width * height * 4 ; i += 4)
    {
        unsigned char        rgb1[3];
        unsigned char        rgb2[3];

        RGBA_To_RGB(data + i, rgb1);
        RGBA_To_RGB(pImage->data + i, rgb2);

        data[i] = abs(rgb1[0] - rgb2[0]);
        data[i+1] = abs(rgb1[1] - rgb2[1]);
        data[i+2] = abs(rgb1[2] - rgb2[2]);
        data[i+3] = 255;
    }

    return true;
}// Difference


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 box filter on this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Box()
{
    if (data == NULL) {
        return false;
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;

            double values[3] = { 0.0 };
            for (int m = -2; m <= 2; m++) {
                for (int n = -2; n <= 2; n++) {
                    int r = i + m, c = j + n;
                    if (r < 0) {
                        r = -r;
                    } else if (r >= height) {
                        r = 2 * (height - 1) - r;
                    }
                    if (c < 0) {
                        c = -c;
                    } else if (c >= width) {
                        c = 2 * (width - 1) - c;
                    }

                    for (int k = 0; k < 3; k++) {
                        values[k] += *(data + (r * width + c) * 4 + k) * 1/25.0;
                    }
                }
            }

            for (int k = 0; k < 3; k++) {
                *(data + pos + k) = Trim(values[k]);
            }
        }
    }

    return true;
}// Filter_Box


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Bartlett filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Bartlett()
{
    if (data == NULL) {
        return false;
    }

    const double mask[5][5] = {
        { 1/81.0, 2/81.0, 3/81.0, 2/81.0, 1/81.0 },
        { 2/81.0, 4/81.0, 6/81.0, 4/81.0, 2/81.0 },
        { 3/81.0, 6/81.0, 9/81.0, 6/81.0, 3/81.0 },
        { 2/81.0, 4/81.0, 6/81.0, 4/81.0, 2/81.0 },
        { 1/81.0, 2/81.0, 3/81.0, 2/81.0, 1/81.0 },
    };

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;

            double values[3] = { 0.0 };
            for (int m = -2; m <= 2; m++) {
                for (int n = -2; n <= 2; n++) {
                    int r = i + m, c = j + n;
                    if (r < 0) {
                        r = -r;
                    } else if (r >= height) {
                        r = 2 * (height - 1) - r;
                    }
                    if (c < 0) {
                        c = -c;
                    } else if (c >= width) {
                        c = 2 * (width - 1) - c;
                    }

                    for (int k = 0; k < 3; k++) {
                        values[k] += *(data + (r * width + c) * 4 + k) * mask[m+2][n+2];
                    }
                }
            }

            for (int k = 0; k < 3; k++) {
                *(data + pos + k) = Trim(values[k]);
            }      
        }
    }

    return true;
}// Filter_Bartlett


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Gaussian()
{
    if (data == NULL) {
        return false;
    }

    const double mask[5][5] = {
        { 1/256.0, 4/256.0, 6/256.0, 4/256.0, 1/256.0 },
        { 4/256.0, 16/256.0, 24/256.0, 16/256.0, 4/256.0 },
        { 6/256.0, 24/256.0, 36/256.0, 24/256.0, 6/256.0 },
        { 4/256.0, 16/256.0, 24/256.0, 16/256.0, 4/256.0 },
        { 1/256.0, 4/256.0, 6/256.0, 4/256.0, 1/256.0 },
    };

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;

            double values[3] = { 0.0 };
            for (int m = -2; m <= 2; m++) {
                for (int n = -2; n <= 2; n++) {
                    int r = i + m, c = j + n;
                    if (r < 0) {
                        r = -r;
                    } else if (r >= height) {
                        r = 2 * (height - 1) - r;
                    }
                    if (c < 0) {
                        c = -c;
                    } else if (c >= width) {
                        c = 2 * (width - 1) - c;
                    }

                    for (int k = 0; k < 3; k++) {
                        values[k] += *(data + (r * width + c) * 4 + k) * mask[m+2][n+2];
                    }
                }
            }

            for (int k = 0; k < 3; k++) {
                *(data + pos + k) = Trim(values[k]);
            }      
        }
    }

    return true;
}// Filter_Gaussian

///////////////////////////////////////////////////////////////////////////////
//
//      Perform NxN Gaussian filter on this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////

bool TargaImage::Filter_Gaussian_N(unsigned int N)
{
    if (data == NULL) {
        return false;
    }

    double mask[N][N];
    double sum = 0;
    for (int i = 0; i < N; i++) {
        double v = Binomial(N-1, i);
        mask[0][i] = v;
        mask[i][0] = v;
        sum += v;
    }
    sum = pow(sum, 2);

    for (int i = 1; i < N; i++) {
        for (int j = 1; j < N; j++) {
            mask[i][j] = (mask[i][0] * mask[0][j]) / (double)sum;
        }
    }
    
    mask[0][0] /= sum;
    for (int i = 1; i < N; i++) {
        mask[i][0] /= (double)sum;
        mask[0][i] /= (double)sum;
    }

    int v = (int)N/2;
    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;

            double values[3] = { 0.0 };
            for (int m = -v; m <= v; m++) {
                for (int n = -v; n <= v; n++) {
                    int r = i + m, c = j + n;
                    if (r < 0) {
                        r = -r;
                    } else if (r >= height) {
                        r = 2 * (height - 1) - r;
                    }
                    if (c < 0) {
                        c = -c;
                    } else if (c >= width) {
                        c = 2 * (width - 1) - c;
                    }

                    for (int k = 0; k < 3; k++) {
                        values[k] += *(data + (r * width + c) * 4 + k) * (double)mask[m+v][n+v];
                    }
                }
            }

            for (int k = 0; k < 3; k++) {
                *(data + pos + k) = Trim(values[k]);
            }
        }
    }

    // print masks
    // for (int i = 0; i < N; i++) {
    //     for (int j = 0; j < N; j++) {
    //         cout << mask[i][j] << ",";
    //     }
    //     cout << "\n";
    // }
    // cout << sum << endl;

   return true;
}// Filter_Gaussian_N


///////////////////////////////////////////////////////////////////////////////
//
//      Perform 5x5 edge detect (high pass) filter on this image.  Return 
//  success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Edge()
{
    if (data == NULL) {
        return false;
    }

    const double mask[5][5] = {
        { -1/256.0, -4/256.0, -6/256.0, -4/256.0, -1/256.0 },
        { -4/256.0, -16/256.0, -24/256.0, -16/256.0, -4/256.0 },
        { -6/256.0, -24/256.0, 220/256.0, -24/256.0, -6/256.0 },
        { -4/256.0, -16/256.0, -24/256.0, -16/256.0, -4/256.0 },
        { -1/256.0, -4/256.0, -6/256.0, -4/256.0, -1/256.0 },
    };

    double rgb_max[3] = { 0.0 };
    double rgb_min[3] = { 0.0 };
    double range[3] = { 0.0 };
    double temp_image[height][width][3];

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;

            double values[3] = { 0.0 };
            for (int m = -2; m <= 2; m++) {
                for (int n = -2; n <= 2; n++) {
                    int r = i + m, c = j + n;
                    if (r < 0) {
                        r = -r;
                    } else if (r >= height) {
                        r = 2 * (height - 1) - r;
                    }
                    if (c < 0) {
                        c = -c;
                    } else if (c >= width) {
                        c = 2 * (width - 1) - c;
                    }

                    for (int k = 0; k < 3; k++) {
                        values[k] += *(data + (r * width + c) * 4 + k) * mask[m+2][n+2];
                    }
                }
            }

            for (int k = 0; k < 3; k++) {
                temp_image[i][j][k] = Trim(values[k]);
                if (temp_image[i][j][k] > rgb_max[k]) {
                    rgb_max[k] = temp_image[i][j][k];
                }
                if (temp_image[i][j][k] < rgb_min[k]) {
                    rgb_min[k] = temp_image[i][j][k];
                }
            }
        }
    }

    for (int i = 0; i < 3; i++) {
        range[i] = (double)rgb_max[i] - (double)rgb_min[i];
        cout << rgb_max[i] << " " << rgb_min[i] << "\n";
    }
    cout << "---" << endl;

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;

            for (int k = 0; k < 3; k++) {
                *(data + pos + k) = (temp_image[i][j][k] - rgb_min[k]) / range[k] * 255.0;
            }
        }
    }

    return true;
}// Filter_Edge


///////////////////////////////////////////////////////////////////////////////
//
//      Perform a 5x5 enhancement filter to this image.  Return success of 
//  operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Filter_Enhance()
{
    if (data == NULL) {
        return false;
    }

        const double mask[5][5] = {
        { -1/256.0, -4/256.0, -6/256.0, -4/256.0, -1/256.0 },
        { -4/256.0, -16/256.0, -24/256.0, -16/256.0, -4/256.0 },
        { -6/256.0, -24/256.0, 476/256.0, -24/256.0, -6/256.0 },
        { -4/256.0, -16/256.0, -24/256.0, -16/256.0, -4/256.0 },
        { -1/256.0, -4/256.0, -6/256.0, -4/256.0, -1/256.0 },
    };

    double rgb_max[3] = { 0.0 };
    double rgb_min[3] = { 0.0 };
    double range[3] = { 0.0 };
    double temp_image[height][width][3];

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;

            double values[3] = { 0.0 };
            for (int m = -2; m <= 2; m++) {
                for (int n = -2; n <= 2; n++) {
                    int r = i + m, c = j + n;
                    if (r < 0) {
                        r = -r;
                    } else if (r >= height) {
                        r = 2 * (height - 1) - r;
                    }
                    if (c < 0) {
                        c = -c;
                    } else if (c >= width) {
                        c = 2 * (width - 1) - c;
                    }

                    for (int k = 0; k < 3; k++) {
                        values[k] += *(data + (r * width + c) * 4 + k) * mask[m+2][n+2];
                    }
                }
            }

            for (int k = 0; k < 3; k++) {
                temp_image[i][j][k] = Trim(values[k]);
                if (temp_image[i][j][k] > rgb_max[k]) {
                    rgb_max[k] = temp_image[i][j][k];
                }
                if (temp_image[i][j][k] < rgb_min[k]) {
                    rgb_min[k] = temp_image[i][j][k];
                }
            }
        }
    }

    for (int i = 0; i < 3; i++) {
        range[i] = (double)rgb_max[i] - (double)rgb_min[i];
    }

    for (int i = 0; i < height; i++) {
        for (int j = 0; j < width; j++) {
            size_t pos = (i * width + j) * 4;

            for (int k = 0; k < 3; k++) {
                *(data + pos + k) = Trim((temp_image[i][j][k] - rgb_min[k]) / range[k] * 255.0);
            }
        }
    }

    return true;
}// Filter_Enhance


///////////////////////////////////////////////////////////////////////////////
//
//      Run simplified version of Hertzmann's painterly image filter.
//      You probably will want to use the Draw_Stroke funciton and the
//      Stroke class to help.
// Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::NPR_Paint()
{
    ClearToBlack();
    return false;
}



///////////////////////////////////////////////////////////////////////////////
//
//      Halve the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Half_Size()
{
    if (data == NULL) {
        return false;
    }

    const double mask[3][3] = {
        { 1/16.0, 1/8.0, 1/16.0 },
        { 1/8.0, 1/4.0, 1/8.0 },
        { 1/16.0, 1/8.0, 1/16.0 },
    };

    int new_height = height / 2;
    int new_width = width / 2;
    uchar *new_data = new uchar[new_height * new_width * 4];

    for (int i = 0; i < new_height; i++) {
        for (int j = 0; j < new_width; j++) {
            size_t pos = (i * new_width + j) * 4;

            int origin_i = 2 * i;
            int origin_j = 2 * j;
            // size_t origin_pos = (origin_i * width + origin_j) * 4;
            double temp_rgb[3] = { 0.0 };
            for (int m = -1; m <= 1; m++) {
                for (int n = -1; n <= 1; n++) {
                    int r = origin_i + m, c = origin_j + n;

                    if (r < 0) {
                        r = -r;
                    } else if (r >= height) {
                        r = 2 * (height - 1) - r;
                    }

                    if (c < 0) {
                        c = -c;
                    } else if (c >= width) {
                        c = 2 * (width - 1) - c;
                    }

                    for (int k = 0; k < 3; k++) {
                        temp_rgb[k] += *(data + (r * width + c) * 4 + k) * mask[m+1][n+1];
                    }
                }
            }

            for (int k = 0; k < 3; k++) {
                *(new_data + pos + k) = temp_rgb[k];
            }
            *(new_data + pos + 3) = 255;
        }

    }
    delete[] data;
    data = new_data;
    height = new_height;
    width = new_width;

    return true;
}// Half_Size

/**
 * 
 * 
 * 
*/
void TargaImage::handleOverbound(int* r, int* c, int& h, int w) {
    if (*r < 0) {
        *r = -*r;
    } else if (*r >= h) {
        *r = 2 * (h - 1) - *r;
    }

    if (*c < 0) {
        *c = -*c;
    } else if (*c >= w) {
        *c = 2 * (w - 1) - *c;
    }
}

///////////////////////////////////////////////////////////////////////////////
//
//      Double the dimensions of this image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Double_Size()
{
    if (data == NULL) {
        return false;
    }

    int new_height = height * 2;
    int new_width = width * 2;
    uchar *new_data = new uchar[new_height * new_width * 4];

    const double both_odd[4][4] = {
        { 1/64.0, 3/64.0, 3/64.0, 1/64.0 },
        { 3/64.0, 9/64.0, 9/64.0, 3/64.0 },
        { 3/64.0, 9/64.0, 9/64.0, 3/64.0 },
        { 1/64.0, 3/64.0, 3/64.0, 1/64.0 }
    };
    const double both_even[3][3] = {
        { 1/16.0, 1/8.0, 1/16.0 },
        { 1/8.0, 1/4.0, 1/8.0 },
        { 1/16.0, 1/8.0, 1/16.0 }
    };
    const double i_even_j_odd[4][3] = {
        { 1/32.0, 2/32.0, 1/32.0 },
        { 3/32.0, 6/32.0, 3/32.0 },
        { 3/32.0, 6/32.0, 3/32.0 },
        { 1/32.0, 2/32.0, 1/32.0 }
    };
    const double i_odd_j_even[3][4] = {
        { 1/32.0, 3/32.0, 3/32.0, 1/32.0 },
        { 2/32.0, 6/32.0, 6/32.0, 2/32.0 },
        { 1/32.0, 3/32.0, 3/32.0, 1/32.0 }
    };

    for (int i = 0; i < new_height; i++) {
        for (int j = 0; j < new_width; j++) {
            size_t pos = (i * new_width + j) * 4;

            double temp_rgb[3] = { 0.0 };
            // i, j odd
            if (i & 1 && j & 1) {
                int origin_i = i / 2 - 1;
                int origin_j = j / 2 - 1;

                for (int m = 0; m < 4; m++) {
                    for (int n = 0; n < 4; n++) {
                        int r = origin_i + m;
                        int c = origin_j + n;
                        handleOverbound(&r, &c, height, width);

                        for (int k = 0; k < 3; k++) {
                            temp_rgb[k] += *(data + (r * width + c) * 4 + k) * both_odd[m][n];
                        }
                    }
                }
            }
            // i, j even
            else if (!(i & 1) && !(j & 1)) {
                int origin_i = i / 2;
                int origin_j = j / 2;
                
                for (int m = 0; m < 3; m++) {
                    for (int n = 0; n < 3; n++) {
                        int r = origin_i + m;
                        int c = origin_j + n;
                        handleOverbound(&r, &c, height, width);

                        for (int k = 0; k < 3; k++) {
                            temp_rgb[k] += *(data + (r * width + c) * 4 + k) * both_even[m][n];
                        }
                    }
                }
            }
            // i odd, j even
            else if (i & 1 && !(j & 1)) {
                int origin_i = i / 2 - 1;
                int origin_j = j / 2 - 1;

                for (int m = 0; m < 3; m++) {
                    for (int n = 0; n < 4; n++) {
                        int r = origin_i + m;
                        int c = origin_j + n;
                        handleOverbound(&r, &c, height, width);

                        for (int k = 0; k < 3; k++) {
                            temp_rgb[k] += *(data + (r * width + c) * 4 + k) * i_odd_j_even[m][n];
                        }
                    }
                }
            }
            // i even, j odd
            else if (!(i & 1) && j & 1) {
                int origin_i = i / 2 - 1;
                int origin_j = j / 2 - 1;

                for (int m = 0; m < 4; m++) {
                    for (int n = 0; n < 3; n++) {
                        int r = origin_i + m;
                        int c = origin_j + n;
                        handleOverbound(&r, &c, height, width);

                        for (int k = 0; k < 3; k++) {
                            temp_rgb[k] += *(data + (r * width + c) * 4 + k) * i_even_j_odd[m][n];
                        }
                    }
                }
            }

            for (int k = 0; k < 3; k++) {
                *(new_data + pos + k) = *(temp_rgb + k);
            }
            *(new_data + pos + 3) = 255;
        }
    }

    delete[] data;
    data = new_data;
    height = new_height;
    width = new_width;

    return true;
}// Double_Size

///////////////////////////////////////////////////////////////////////////////
//
//      Scale the image dimensions by the given factor.  The given factor is 
//  assumed to be greater than one.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Resize(float scale)
{
    ClearToBlack();
    return false;
}// Resize


//////////////////////////////////////////////////////////////////////////////
//
//      Rotate the image clockwise by the given angle.  Do not resize the 
//  image.  Return success of operation.
//
///////////////////////////////////////////////////////////////////////////////
bool TargaImage::Rotate(float angleDegrees)
{
    ClearToBlack();
    return false;
}// Rotate


//////////////////////////////////////////////////////////////////////////////
//
//      Given a single RGBA pixel return, via the second argument, the RGB
//      equivalent composited with a black background.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::RGBA_To_RGB(unsigned char *rgba, unsigned char *rgb)
{
    const unsigned char	BACKGROUND[3] = { 0, 0, 0 };

    unsigned char  alpha = rgba[3];

    if (alpha == 0)
    {
        rgb[0] = BACKGROUND[0];
        rgb[1] = BACKGROUND[1];
        rgb[2] = BACKGROUND[2];
    }
    else
    {
	    float	alpha_scale = (float)255 / (float)alpha;
	    int	val;
	    int	i;

	    for (i = 0 ; i < 3 ; i++)
	    {
	        val = (int)floor(rgba[i] * alpha_scale);
	        if (val < 0)
		    rgb[i] = 0;
	        else if (val > 255)
		    rgb[i] = 255;
	        else
		    rgb[i] = val;
	    }
    }
}// RGA_To_RGB


///////////////////////////////////////////////////////////////////////////////
//
//      Copy this into a new image, reversing the rows as it goes. A pointer
//  to the new image is returned.
//
///////////////////////////////////////////////////////////////////////////////
TargaImage* TargaImage::Reverse_Rows(void)
{
    unsigned char   *dest = new unsigned char[width * height * 4];
    TargaImage	    *result;
    int 	        i, j;

    if (! data)
    	return NULL;

    for (i = 0 ; i < height ; i++)
    {
	    int in_offset = (height - i - 1) * width * 4;
	    int out_offset = i * width * 4;

	    for (j = 0 ; j < width ; j++)
        {
	        dest[out_offset + j * 4] = data[in_offset + j * 4];
	        dest[out_offset + j * 4 + 1] = data[in_offset + j * 4 + 1];
	        dest[out_offset + j * 4 + 2] = data[in_offset + j * 4 + 2];
	        dest[out_offset + j * 4 + 3] = data[in_offset + j * 4 + 3];
        }
    }

    result = new TargaImage(width, height, dest);
    delete[] dest;
    return result;
}// Reverse_Rows


///////////////////////////////////////////////////////////////////////////////
//
//      Clear the image to all black.
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::ClearToBlack()
{
    memset(data, 0, width * height * 4);
}// ClearToBlack


///////////////////////////////////////////////////////////////////////////////
//
//      Helper function for the painterly filter; paint a stroke at
// the given location
//
///////////////////////////////////////////////////////////////////////////////
void TargaImage::Paint_Stroke(const Stroke& s) {
   int radius_squared = (int)s.radius * (int)s.radius;
   for (int x_off = -((int)s.radius); x_off <= (int)s.radius; x_off++) {
      for (int y_off = -((int)s.radius); y_off <= (int)s.radius; y_off++) {
         int x_loc = (int)s.x + x_off;
         int y_loc = (int)s.y + y_off;
         // are we inside the circle, and inside the image?
         if ((x_loc >= 0 && x_loc < width && y_loc >= 0 && y_loc < height)) {
            int dist_squared = x_off * x_off + y_off * y_off;
            if (dist_squared <= radius_squared) {
               data[(y_loc * width + x_loc) * 4 + 0] = s.r;
               data[(y_loc * width + x_loc) * 4 + 1] = s.g;
               data[(y_loc * width + x_loc) * 4 + 2] = s.b;
               data[(y_loc * width + x_loc) * 4 + 3] = s.a;
            } else if (dist_squared == radius_squared + 1) {
               data[(y_loc * width + x_loc) * 4 + 0] = 
                  (data[(y_loc * width + x_loc) * 4 + 0] + s.r) / 2;
               data[(y_loc * width + x_loc) * 4 + 1] = 
                  (data[(y_loc * width + x_loc) * 4 + 1] + s.g) / 2;
               data[(y_loc * width + x_loc) * 4 + 2] = 
                  (data[(y_loc * width + x_loc) * 4 + 2] + s.b) / 2;
               data[(y_loc * width + x_loc) * 4 + 3] = 
                  (data[(y_loc * width + x_loc) * 4 + 3] + s.a) / 2;
            }
         }
      }
   }
}


///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke() {}

///////////////////////////////////////////////////////////////////////////////
//
//      Build a Stroke
//
///////////////////////////////////////////////////////////////////////////////
Stroke::Stroke(unsigned int iradius, unsigned int ix, unsigned int iy,
               unsigned char ir, unsigned char ig, unsigned char ib, unsigned char ia) :
   radius(iradius),x(ix),y(iy),r(ir),g(ig),b(ib),a(ia)
{
}

