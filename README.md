# Fundamental-Matrix

### About

#### Corner Detection

Load images and detect corners using the code from my [Corner Detection](https://github.com/wxue/Corner-Detection), which I improved here by adding local maximum.

#### Extract Patches

To extract patches (by 5*5), and form a descriptor simply by vectorizing the image pixel value in raster scan order. Store data in “Patches_left” and “Patches_right”.

#### Get Fundamental Matrix

* Make pairs and imply the corre2 function to evaluate the correlation, and select pairs with lager correlation.

* Filter the pairs to make sure each point in one image only be paired with one point in another image.

* Imply 8 point RANSAC 10 times as an example(30~50 would be better. I choose 10 only for running time sake).

* Let the obtain fundamental matrix test all the paired points, evaluate its weights.

* Use the obtained fundamental matrix to select the correlation pairs once again. 

* Run 8 point RANSAC again to get a better fundamental matrix.

> Credits to [Professor Hua](http://www.cs.stevens.edu/~ghua/)

### License

The MIT License (MIT)

Copyright (c) 2014 Weiyu Xue

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
