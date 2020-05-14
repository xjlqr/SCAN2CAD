#include "SCAN2CAD.h"

using namespace std;
using namespace pcl;

class BoundingBox {
public:
    Eigen::Vector3d min_point;
    Eigen::Vector3d max_point;
    bool contains(Eigen::Vector3d);
};

bool BoundingBox::contains(Eigen::Vector3d point) {
    return point[0] > min_point[0] && point[0] < max_point[0] &&
        point[1] > min_point[1] && point[1] < max_point[1] &&
        point[2] > min_point[2] && point[2] < max_point[2];
}

Eigen::Vector3d find_point_intersection(ModelCoefficients plane_1, ModelCoefficients plane_2, ModelCoefficients plane_3) {
    Eigen::Vector3d point(0, 0, 0);
    Eigen::Vector3d n1(plane_1.values[0], plane_1.values[1], plane_1.values[2]);
    Eigen::Vector3d n2(plane_2.values[0], plane_2.values[1], plane_2.values[2]);
    Eigen::Vector3d n3(plane_3.values[0], plane_3.values[1], plane_3.values[2]);

    Eigen::Matrix3d plane_matrix;
    plane_matrix << n1[0], n2[0], n3[0], n1[1], n2[1], n3[1], n1[2], n2[2], n3[2];
    double det = plane_matrix.determinant();
    if (abs(det) > 0) {
        point = (n2.cross(n3) * -plane_1.values[3] + n3.cross(n1) * -plane_2.values[3] + n1.cross(n2) * -plane_3.values[3]) / det;
        return point;
    }
    else {
        point = n1 * (NAN, NAN, NAN);
        return point;
    }
}

line find_line_intersection(ModelCoefficients plane_1, ModelCoefficients plane_2) {
    line intersection;
    Eigen::Vector3d n1(plane_1.values[0], plane_1.values[1], plane_1.values[2]);
    Eigen::Vector3d n2(plane_2.values[0], plane_2.values[1], plane_2.values[2]);

    intersection.direction = n1.cross(n2);

    ModelCoefficients plane_3;
    plane_3.values.push_back(intersection.direction[0]);
    plane_3.values.push_back(intersection.direction[1]);
    plane_3.values.push_back(intersection.direction[2]);
    plane_3.values.push_back(0.0);


    intersection.point = find_point_intersection(plane_1, plane_2, plane_3);

    return intersection;
}

const double EPSILON = 0.0001;

int main(int argc, char** argv) {


    //Create test dataset (lines 60-142)
    PointCloud<PointXYZ> testbeam;

    //Test beam parameters
    vector<double> testbeam_position{ 30.0, 30.0, 30.0 }; // Must be away from origin
    vector<int> testbeam_size{ 120, 30, 30 }; //must be divisible by 6
    double testbeam_noiselevel = 0.4;

    //Create beam (todo: factor code)
    int size_x = testbeam_size[0];
    int size_y = testbeam_size[1];
    int size_z = testbeam_size[2];
    cout << (size_x % 6 != 0 || size_y % 6 != 0 || size_z % 6 != 0 ? "Warning: Sides of test beam must be divisible by 6" : "Creating Test Beam") << endl;
    testbeam.width = size_x * size_y * size_z;
    testbeam.height = 1;
    testbeam.is_dense = false;
    testbeam.points.resize(testbeam.width);

    Eigen::Vector3d minpoint(0.0, 0.0, 0.0);
    Eigen::Vector3d maxpoint(0.0, 0.0, 0.0);



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

                    //God-ugly but PointXYZs can't be looped 
                    if (testbeam[beam_size].x > maxpoint[0]) {
                        maxpoint[0] = testbeam[beam_size].x;
                    }
                    if (testbeam[beam_size].y > maxpoint[1]) {
                        maxpoint[1] = testbeam[beam_size].y;
                    }
                    if (testbeam[beam_size].z > maxpoint[2]) {
                        maxpoint[2] = testbeam[beam_size].z;
                    }
                    if (testbeam[beam_size].x < minpoint[0]) {
                        minpoint[0] = testbeam[beam_size].x;
                    }
                    if (testbeam[beam_size].y < minpoint[1]) {
                        minpoint[1] = testbeam[beam_size].y;
                    }
                    if (testbeam[beam_size].z < minpoint[2]) {
                        minpoint[2] = testbeam[beam_size].z;
                    }

                    beam_size++;
                }
            }
        }
    }

    maxpoint[0] += 5.0;
    maxpoint[1] += 5.0;
    maxpoint[2] += 5.0;
    minpoint[0] -= 5.0;
    minpoint[1] -= 5.0;
    minpoint[2] -= 5.0;

    BoundingBox bbox;
    bbox.min_point = minpoint;
    bbox.max_point = maxpoint;

    cout << "Boundingbox: (" << bbox.min_point[0] << ", " << bbox.min_point[1] << ", " << bbox.min_point[2] << "); (" << bbox.max_point[0] << ", " << bbox.max_point[1] << ", " << bbox.max_point[2] << ")" << endl; 

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

    double res = 0.0;
    int n_points = 0;
    int nres;
    std::vector<int> indices(2);
    std::vector<float> sqr_distances(2);
    pcl::search::KdTree<PointXYZ> res_tree;
    res_tree.setInputCloud(beam);

    for (std::size_t i = 0; i < beam->size(); ++i)
    {
        //Considering the second neighbor since the first is the point itself.
        nres = res_tree.nearestKSearch(i, 2, indices, sqr_distances);
        if (nres == 2)
        {
            res += sqrt(sqr_distances[1]);
            ++n_points;
        }
    }
    if (n_points != 0)
    {
        res /= n_points;
    }


    cout << "cloud resolution: " << res << endl;
    //Objects for storing algorithm data
    string filename = "";
    ModelCoefficients::Ptr coefficients(new ModelCoefficients);
    vector<ModelCoefficients> planes(0);
    PointIndices::Ptr inliers(new PointIndices);
    PointCloud<PointXYZ>::Ptr plane_points(new PointCloud<PointXYZ>);
    PointCloud<PointXYZ>::Ptr clustered(new PointCloud<PointXYZ>);
    PointCloud<PointXYZ>::Ptr output(new PointCloud<PointXYZ>);
    PointCloud<PointXYZ>::Ptr filtered(new PointCloud<PointXYZ>);
    //Object for removing points from cloud
    ExtractIndices<PointXYZ> extract;
    extract.setInputCloud(beam);

    //Load RANSAC algorithm
    SACSegmentation<PointXYZ> seg;
    seg.setInputCloud(beam);
    seg.setModelType(SACMODEL_PLANE); //Fitted model to dataset
    seg.setMethodType(SAC_RRANSAC); //Method used (RANSAC is the simplest)
    seg.setDistanceThreshold(0.4); //Maximum distance to model, should be greater than noiselevel
    seg.setMaxIterations(500); //Iterations tried. More means better plane
    seg.setOptimizeCoefficients(true);


    //Clustering
    EuclideanClusterExtraction<PointXYZ> ece;
    search::KdTree<PointXYZ>::Ptr ece_tree(new search::KdTree<PointXYZ>);
    ece.setSearchMethod(ece_tree);
    ece.setClusterTolerance(2.0*sqrt(2.0)*res);
    ece.setMinClusterSize(20);
    ece.setMaxClusterSize(1000000000);

    RadiusOutlierRemoval<PointXYZ> ror1(true);
    ror1.setNegative(false);
    ror1.setRadiusSearch((1.0+testbeam_noiselevel)*sqrt(2.0));
    ror1.setMinNeighborsInRadius(3);

    RadiusOutlierRemoval<PointXYZ> ror2(true);
    ror2.setNegative(false);
    ror2.setRadiusSearch(1.0+2.0*testbeam_noiselevel+0.01);
    ror2.setMinNeighborsInRadius(2);

    // The indices_rem array indexes all points of cloud_in that have 5 or more neighbors within the 0.1 search radius




    int cloud_index = 0;
    while (true) {
        cout << "Plane " << cloud_index+1 << ":" << endl;
        cout << "\tFitting plane" << endl;
        seg.segment(*inliers, *coefficients);
        vector<PointIndices> clusters;
        
        
        if (inliers->indices.size() > 0) { //If plane found
            planes.push_back(*coefficients);
            cout << "\tPlane found with coefficients " << planes[cloud_index].values[0] << ", " << planes[cloud_index].values[1] << ", " << planes[cloud_index].values[2] << ", " << planes[cloud_index].values[3] << " copying points" << endl;
            copyPointCloud(*beam, inliers->indices, *plane_points);

            filename = "test_ransac_";
            filename += string(to_string(cloud_index));
            filename += ".pcd";
            cout << "\tSaving " << filename << " with " << inliers->indices.size() << " points" << endl;
            io::savePCDFileBinary(filename, *plane_points);

            //Keep only points belonging to plane
            extract.setIndices(inliers);
            extract.setNegative(false);
            extract.filter(*plane_points);
            
            // Large radius search isolates corners
            PointCloud<PointXYZ>::Ptr complete_plane(new PointCloud<PointXYZ>);
            copyPointCloud(*plane_points, *complete_plane);
            ror1.setInputCloud(plane_points);
            ror1.filter(*plane_points);
            PointCloud<PointXYZ>::Ptr removedPoints1(new PointCloud<PointXYZ>);
            auto removedIndices1 = ror1.getRemovedIndices();
            copyPointCloud(*complete_plane, *removedIndices1, *removedPoints1);
            cout << "\tLarge radius KNN search found " << removedPoints1->size() << " points" << endl;

            //Small radius snips corners
            copyPointCloud(*plane_points, *complete_plane);
            ror2.setInputCloud(plane_points);
            ror2.filter(*plane_points);
            PointCloud<PointXYZ>::Ptr removedPoints2(new PointCloud<PointXYZ>);
            auto removedIndices2 = ror2.getRemovedIndices();
            copyPointCloud(*complete_plane, *removedIndices2, *removedPoints2);
            cout << "\tSmall radius KNN search found " << removedPoints1->size() << " points" << endl;
                                                                
            //Cluster remaining points into surfaces            
            ece.setInputCloud(plane_points);
            ece.extract(clusters);
            cout << "\tClustering done, " << clusters.size() << " cluster(s) found" << endl;
            for (int i = 0; i < clusters.size(); i++) {
                cout << "\t\tCluster " << i << endl;

                PointCloud<PointXYZ>::Ptr cluster(new PointCloud<PointXYZ>);

                //Make_shared doesn't work because
                PointIndices::Ptr indices(new PointIndices);
                indices->header = clusters[i].header;
                indices->indices = clusters[i].indices;

                copyPointCloud(*plane_points, clusters[i].indices, *cluster);
                filename = "test_cluster_";
                filename += string(to_string(cloud_index));
                filename += string(to_string(i));
                filename += ".pcd";
                cout << "\t\tSaving " << filename << endl;
                io::savePCDFileBinary(filename, *cluster);

                cout << "\t\tCalculating mesh normals" << endl;
                search::KdTree<PointXYZ>::Ptr normal_tree(new search::KdTree<PointXYZ>);
                NormalEstimationOMP<PointXYZ, Normal> neomp; //OMP version multithreaded
                PointCloud<Normal>::Ptr normals(new PointCloud<Normal>);
                normal_tree->setInputCloud(cluster);
                neomp.setInputCloud(cluster);
                neomp.setViewPoint(0, 0, 0);
                neomp.setSearchMethod(normal_tree);
                neomp.setKSearch(32);//Search depth - higher means smoother output with less detail (shouldn't matter for planes)

                //Compute normals
                neomp.compute(*normals);
                PointCloud<PointNormal>::Ptr cluster_with_normals(new PointCloud<PointNormal>);
                concatenateFields(*cluster, *normals, *cluster_with_normals);

                cout << "\t\tTriangulating mesh" << endl;
                search::KdTree<PointNormal>::Ptr mesh_tree(new search::KdTree<PointNormal>);
                mesh_tree->setInputCloud(cluster_with_normals);

                //Greedy Triangulation Method
                GreedyProjectionTriangulation<PointNormal> GPT;
                GPT.setInputCloud(cluster_with_normals);
                GPT.setSearchMethod(mesh_tree);
                GPT.setSearchRadius((1.0+testbeam_noiselevel)*sqrt(2.0)-EPSILON);
                GPT.setMaximumAngle(2.0 * M_PI / 3.0);
                GPT.setMinimumAngle(2 * 0 * M_PI / 36.0);
                GPT.setMaximumSurfaceAngle(2 * M_PI / 8);
                GPT.setMaximumNearestNeighbors(100);
                GPT.setMu(3.0);
                GPT.setNormalConsistency(false);
                PolygonMesh triangles;
                GPT.reconstruct(triangles);
                
                //ConcaveHull Method
                //ConcaveHull<PointNormal> cchul;
                //cchul.setInputCloud(cluster_with_normals);
                //cchul.setSearchMethod(mesh_tree);
                //cchul.setAlpha(res*2);
                //cchul.setDimension(2);
                //PolygonMesh triangles;
                //cchul.reconstruct(triangles);

                //Write mesh to file
                filename = "mesh_";
                filename += string(to_string(cloud_index));
                filename += string(to_string(i));
                filename += ".vtk";
                cout << "\t\tSaving " << filename << endl << endl;
                io::saveVTKFile(filename, triangles);

            }

            if (clusters.size() == 0) {
                break;
            }

            //Remove spent points from dataset
            cout << "\tDeleting RANSAC points from dataset" << endl;
            extract.setIndices(inliers);
            extract.setNegative(true);
            extract.filter(*beam);

            cout << "\tAdding KNN-removed points back" << endl;
            bool success = concatenate(*beam, *removedPoints1, *beam);
            success = concatenate(*beam, *removedPoints2, *beam);
        }


        if (beam->size() == 0 || clusters.size() == 0) {//If no more points are left, or no clusters have been found, object is segmented, break loop
            cout << "RANSAC Plane Fitting Ended!" << endl << endl;
            break;
        }
        else {
        }
        cloud_index++;
    }

    //Find intersections between three planes (corners of the model)
    vector<vector<vector<bool>>> corner_connection_map(planes.size(), vector<vector<bool>>(planes.size(), vector<bool>(planes.size(), false)));
    vector<PointXYZ> corners;
    if (planes.size() > 2) {
        cout << "Checking plane-point intersections" << endl;
        for (int i = 0; i < planes.size(); i++) {
            for (int j = i + 1; j < planes.size(); j++) {
                for (int k = j + 1; k < planes.size(); k++){
                    Eigen::Vector3d point = find_point_intersection(planes[i], planes[j], planes[k]);
                    if (!isnan(point[0]) && bbox.contains(point)) {
                        corner_connection_map[i][j][k] = true;
                        cout << "(" << i << ", " << j << ", " << k << ") have corner (" << point[0] << ", " << point[1] << ", " << point[2] << ")" << endl;
                        corners.push_back(PointXYZ(point[0], point[1], point[2]));
                    }
                }
            }
        }
    }

    //find intersections between two planes (edges of the model)
    vector<vector<bool>> edge_connection_map(planes.size(), vector<bool>(planes.size(), false));
    vector<line> edges;
    if (planes.size() > 1) {
        cout << "Checking plane-line intersections" << endl;
        for (int i = 0; i < planes.size(); i++) {
            for (int j = 0; j < planes.size(); j++) {
                if (j != i) {
                    line edge = find_line_intersection(planes[i], planes[j]);
                    if (!isnan(edge.point[0])) {
                        edge_connection_map[i][j] = true;
                        edges.push_back(edge);
                        cout << "(" << edge.point[0] << ", " << edge.point[1] << ", " << edge.point[2] << "), (" << edge.direction[0] << ", " << edge.direction[1] << ", " << edge.direction[2] << ")" << endl;
                    }
                }
            }
        }
    }
    /*TODO: Mesh projection, plane segmentation mesh stitching */

    return (0);
}
