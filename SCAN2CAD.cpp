#include "SCAN2CAD.h"

using namespace std;
using namespace pcl;

int main(int argc, char** argv) {
    //Create neat tidy example
    PointCloud<PointXYZ> testbeam;

    //Test beam parameters
    vector<double> testbeam_position{ 30.0, 30.0, 30.0 }; // Must be away from origin
    vector<int> testbeam_size{ 240, 120, 120 }; //must be divisible by 6
    double testbeam_noiselevel = 0.5;

    //Create beam (todo: factor code)
    int size_x = testbeam_size[0];
    int size_y = testbeam_size[1];
    int size_z = testbeam_size[2];
    cout << (size_x % 6 != 0 || size_y % 6 != 0 || size_z % 6 != 0 ? "Warning: Sides of test beam must be divisible by 6" : "Creating Test Beam") << endl;
    testbeam.width = size_x * size_y * size_z;
    testbeam.height = 1;
    testbeam.is_dense = false;
    testbeam.points.resize(testbeam.width);

    unsigned int beam_size = 0;
    for (int i = 0; i < size_x; i++) {//Loop through all grid points in the (x, y, z) box
        for (int j = 0; j < size_y; j++) {
            for (int k = 0; k < size_z; k++) {
                if (//Include only these points in the beam:
                    k == 0 || k == (size_z - 1) || //top and bottom
                    (j == 0 || j == (size_y - 1)) && (k <= size_z / 6 || k >= size_z - (size_z / 6)) || //side edges of bottom/top
                    (abs((j - size_y / 2)) == size_y / 6) && (!(k <= size_z / 6 || k >= size_z - (size_z / 6))) || //Inside edge of botto/top
                    (abs((k - size_z / 2)) == size_z / 3 && (abs((j - size_y / 2)) >= size_y / 6)) || //Stem 
                    (abs((j - size_y / 2)) <= size_y / 6 || abs((k - size_z / 2)) >= size_z / 3) && (i == 0 || i == size_x - 1)//End caps
                    ) {
                    testbeam[beam_size].x = i + testbeam_position[0] + testbeam_noiselevel * noise();
                    testbeam[beam_size].y = j + testbeam_position[1] + testbeam_noiselevel * noise();
                    testbeam[beam_size].z = k + testbeam_position[2] + testbeam_noiselevel * noise();
                    beam_size++;
                }
            }
        }
    }

    testbeam.width = beam_size;
    testbeam.points.resize(beam_size);

    cout << "Saving pcd..." << endl;
    io::savePCDFileBinary("test_beam.pcd", testbeam);
    cout << "Saved " << testbeam.size() << " data points.\n" << endl;


    //Main part 1: RANSAC Plane fitting
    cout << "RANSAC Plane Fitting Started!" << endl << endl;
    //Load file
    pcl::PointCloud<pcl::PointXYZ>::Ptr beam(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PCDReader reader;
    reader.read("test_beam.pcd", *beam);


    //Objects for storing algorithm data
    string filename = "";
    ModelCoefficients::Ptr coefficients(new ModelCoefficients);
    PointIndices::Ptr inliers(new PointIndices);
    PointCloud<PointXYZ>::Ptr plane_points(new PointCloud<PointXYZ>);
    PointCloud<PointXYZ>::Ptr clustered(new PointCloud<PointXYZ>);
    PointCloud<PointXYZ>::Ptr output(new PointCloud<PointXYZ>);

    //Object for removing points from cloud
    ExtractIndices<PointXYZ> extract;
    extract.setInputCloud(beam);

    //Load segmentation algorithm
    SACSegmentation<PointXYZ> seg;
    seg.setInputCloud(beam);
    seg.setModelType(SACMODEL_PLANE); //Fitted model to dataset
    seg.setMethodType(SAC_PROSAC); //Method used (RANSAC is the simplest)
    seg.setDistanceThreshold(0.6); //Maximum distance to model, should be greater than noiselevel
    seg.setMaxIterations(500); //Iterations tried. More means better plane

    int cloud_index = 0;
    while (true) {
        cout << "Plane " << cloud_index+1 << ":" << endl;
        cout << "\tFitting plane" << endl;
        seg.segment(*inliers, *coefficients); 
        
        if (inliers->indices.size() > 0) { //If plane found
            cout << "\tPlane found, copying points" << endl;
            copyPointCloud(*beam, inliers->indices, *plane_points); 
  
            //Keep only points belonging to plane
            extract.setIndices(inliers); 
            extract.setNegative(false);
            extract.filter(*plane_points);

            //Write plane to file
            filename = "test_ransac_";
            filename += string(to_string(cloud_index));
            filename += ".pcd";
            cout << "\tSaving " << filename << " with " << inliers->indices.size() << " points" << endl;
            io::savePCDFileBinary(filename, *plane_points);

            cout << "\tDeleting found points from dataset" << endl;
            extract.setNegative(true);
            extract.filter(*beam);

            cout << "\tCalculating mesh normals" << endl;
            search::KdTree<PointXYZ>::Ptr normal_tree(new search::KdTree<PointXYZ>);
            NormalEstimationOMP<PointXYZ, Normal> neomp; //OMP version multithreaded
            PointCloud<Normal>::Ptr normals(new PointCloud<Normal>);
            normal_tree->setInputCloud(plane_points);
            neomp.setInputCloud(plane_points);
            neomp.setViewPoint(0, 0, 0);
            neomp.setSearchMethod(normal_tree);
            neomp.setKSearch(32);//Search depth - higher means smoother output with less detail (shouldn't matter for planes)

            //Compute normals
            neomp.compute(*normals);
            PointCloud<PointNormal>::Ptr plane_with_normals(new PointCloud<PointNormal>);
            concatenateFields(*plane_points, *normals, *plane_with_normals);

            cout << "\tTriangulating mesh" << endl;
            search::KdTree<PointNormal>::Ptr mesh_tree(new search::KdTree<PointNormal>);
            mesh_tree->setInputCloud(plane_with_normals);
            GreedyProjectionTriangulation<PointNormal> GPT;
            GPT.setInputCloud(plane_with_normals);
            GPT.setSearchMethod(mesh_tree);
            
            //Set values for the parameters
            GPT.setSearchRadius(2.0);
            GPT.setMu(2.5);
            GPT.setMaximumNearestNeighbors(100);
            GPT.setMaximumSurfaceAngle(M_PI / 2);
            GPT.setMinimumAngle(0);
            GPT.setMaximumAngle(M_PI);
            GPT.setNormalConsistency(false);

            //Compute mesh
            PolygonMesh triangles;
            GPT.reconstruct(triangles);

            //Write mesh to file
            filename = "mesh_";
            filename += string(to_string(cloud_index));
            filename += ".vtk";
            cout << "\tSaving " << filename << endl << endl;
            io::saveVTKFile(filename, triangles);
        }


        if (beam->size() == 0) {//If no more points are left, object is segmented, break loop
            cout << "RANSAC Plane Fitting Ended!" << endl << endl;
            break;
        }
        cloud_index++;
    }

    /*TODO: Mesh projection, plane segmentation mesh stitching */

    return (0);
}
