﻿// PCL-test-3.h : Include file for standard system include files,
// or project specific include files.

#pragma once

//C++/C libs
#include <iostream>
#include <thread>
#include <vector>
#include <utility>
#include <cmath>
#include <cstring>


//PCL components
#include <pcl/point_types.h>
#include <pcl/ModelCoefficients.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/vtk_io.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/features/normal_3d_omp.h>
#include <pcl/filters/extract_indices.h>
#include <pcl/sample_consensus/ransac.h>
#include <pcl/sample_consensus/sac_model_plane.h>
#include <pcl/sample_consensus/method_types.h>
#include <pcl/sample_consensus/model_types.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/segmentation/extract_clusters.h>
#include <pcl/surface/gp3.h>

//Other stuff
#include "simple_fft/fft_settings.h"
#include "simple_fft/fft.h"

unsigned int pow2_ceil(size_t number) {
    unsigned int nearest_power = 1;//Find nearest power of 2 for fft
    while (nearest_power < number) {
        nearest_power *= 2;
    }
    return nearest_power;
}

double noise() {
    return (rand() / (RAND_MAX + 1.0f) - 0.5);
}