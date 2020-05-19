#include "SCAN2CAD.h"

//Print during RANSAC, save PCD files before output
#define DEBUG
//#define DEBUG_VERBOSE

const double EPSILON = 0.0001; //Boundingbox check error margin
const int CLUSTER_MIN_SIZE = 64; //Minimum amount of points in a cluster
const double CLUSTER_MESH_DISTANCE_THRESHOLD = 5.0; //Maximum distance to cluster

using namespace std;
using namespace pcl;

class BoundingBox {
public:
    Eigen::Vector3d min_point;
    Eigen::Vector3d max_point;
    bool contains(Eigen::Vector3d);
    void findDimensions(PointCloud<PointXYZ>);
};

void BoundingBox::findDimensions(PointCloud<PointXYZ> cloud){
    Eigen::Vector3d minpoint(cloud[0].x, cloud[0].y, cloud[0].z);
    Eigen::Vector3d maxpoint(cloud[0].x, cloud[0].y, cloud[0].z);

    for (int i = 0; i < cloud.size(); i++) {
        if (cloud[i].x > maxpoint[0]) {
            maxpoint[0] = cloud[i].x;
        }
        if (cloud[i].y > maxpoint[1]) {
            maxpoint[1] = cloud[i].y;
        }
        if (cloud[i].z > maxpoint[2]) {
            maxpoint[2] = cloud[i].z;
        }
        if (cloud[i].x < minpoint[0]) {
            minpoint[0] = cloud[i].x;
        }
        if (cloud[i].y < minpoint[1]) {
            minpoint[1] = cloud[i].y;
        }
        if (cloud[i].z < minpoint[2]) {
            minpoint[2] = cloud[i].z;
        }
    }

    maxpoint[0] += EPSILON;
    maxpoint[1] += EPSILON;
    maxpoint[2] += EPSILON;
    minpoint[0] -= EPSILON;
    minpoint[1] -= EPSILON;
    minpoint[2] -= EPSILON;

    min_point = minpoint;
    max_point = maxpoint;

    #ifdef DEBUG_VERBOSE 
        cout << "BoundingBox.findDimensions: (" << min_point[0] << ", " << min_point[1] << ", " << min_point[2] << "); (" << max_point[0] << ", " << max_point[1] << ", " << max_point[2] << ")" << endl;
    #endif
}

bool BoundingBox::contains(Eigen::Vector3d point) {
    bool result = point[0] > min_point[0] && point[0] < max_point[0] && point[1] > min_point[1] && point[1] < max_point[1] && point[2] > min_point[2] && point[2] < max_point[2];
    #ifdef DEBUG_VERBOSE 
        cout << "BoundingBox.contains: (" << point[0] << ", " << point[1] << ", " << point[2] << ") = " << (result?"true":"false") << endl;
    #endif
    return result;
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
        #ifdef DEBUG_VERBOSE
            cout << "find_point_intersection: (" << point[0] << "," << point[1] << "," << point[2] << ")" << endl;
        #endif // DEBUG

        return point;
    }
    else {
        point = n1 * (NAN, NAN, NAN);
        #ifdef DEBUG_VERBOSE
            cout << "find_point_intersection: NAN" << endl;
        #endif // DEBUG
        return point;
    }
}

vector<PointXYZ> find_corner_points(vector<ModelCoefficients> planes, BoundingBox bbox) {
    vector<PointXYZ> corners;

    for (int i = 0; i < planes.size(); i++) {
        for (int j = i + 1; j < planes.size(); j++) {
            for (int k = j + 1; k < planes.size(); k++) {
                Eigen::Vector3d point = find_point_intersection(planes[i], planes[j], planes[k]);
                if (!isnan(point[0]) && bbox.contains(point)) {
                    PointXYZ corner_point = PointXYZ(point[0], point[1], point[2]);
                    corners.push_back(corner_point);
                }
            }
        }

    }

    #ifdef DEBUG_VERBOSE
        cout << "find_corner_points: " << corners.size() << " points found" << endl;
    #endif // DEBUG

    return corners;
}

double closest_point_dist(PointXYZ point, PointCloud<PointXYZ> cloud ) {
    double dist = ((double) point.x - (double) cloud[0].x) * ((double) point.x - (double) cloud[0].x) + 
                  ((double) point.y - (double) cloud[0].y) * ((double) point.y - (double) cloud[0].y) + 
                  ((double) point.z - (double) cloud[0].z) * ((double) point.z - (double) cloud[0].z);

    for (int i = 0; i < cloud.size(); i++) {
        double tmp = ((double) point.x - (double) cloud[i].x) * ((double) point.x - (double) cloud[i].x) +
                     ((double) point.y - (double) cloud[i].y) * ((double) point.y - (double) cloud[i].y) +
                     ((double) point.z - (double) cloud[i].z) * ((double) point.z - (double) cloud[i].z);
        if (tmp < dist) {
            dist = tmp;
        }
    }

    dist = sqrt(dist);

    #ifdef DEBUG_VERBOSE
    cout << "closest_point_dist: " << dist << " from (" << point.x << ", " << point.y << ", " << point.y << ") to cloud" << endl;
    #endif // DEBUG

    return dist;
}

int main(int argc, char** argv) {
    
    cout << "Loading file... ";

    //Load file
    pcl::PointCloud<pcl::PointXYZ>::Ptr beam(new pcl::PointCloud<pcl::PointXYZ>);
    pcl::PCDReader reader;
    reader.read("test_beam.pcd", *beam);
    BoundingBox bbox;
    bbox.findDimensions(*beam);

    cout << "done!" << endl;


    //Main part 1: RANSAC Plane fitting
    cout << "Initiating RANSAC deconstruction..." << endl;

    //Objects for storing algorithm data
    string filename = "";
    ModelCoefficients::Ptr coefficients(new ModelCoefficients);
    vector<ModelCoefficients> planes(0);
    PointIndices::Ptr inlier_indices(new PointIndices);
    PointCloud<PointXYZ>::Ptr inliers(new PointCloud<PointXYZ>);

    //Object for removing points from cloud
    ExtractIndices<PointXYZ> extract;
    extract.setInputCloud(beam);

    //Load RANSAC algorithm
    SACSegmentation<PointXYZ> seg;
    seg.setInputCloud(beam);
    seg.setModelType(SACMODEL_PLANE);       //Fitted model to dataset (plane)
    seg.setMethodType(SAC_RANSAC);          //Method used (RANSAC is the simplest)
    seg.setDistanceThreshold(0.4);          //Maximum distance to model for inliers (should be greater than noiselevel of data)
    seg.setMaxIterations(500);              //Iterations tried. More means better plane; ~500 seems sufficient for test model.
    seg.setOptimizeCoefficients(true);      //Fit new plane to inliers after inliers found?

    //Clustering
    EuclideanClusterExtraction<PointXYZ> ece;
    search::KdTree<PointXYZ>::Ptr ece_tree(new search::KdTree<PointXYZ>);
    ece.setSearchMethod(ece_tree);
    ece.setClusterTolerance(2.0*sqrt(2.0)); //Max distance between any two points to be considered cluster. Should be around ~2*
    ece.setMinClusterSize(CLUSTER_MIN_SIZE);          //Min number of points to be considered cluster
    ece.setMaxClusterSize(1000000000);  //Max number of points to be considered cluster

    RadiusOutlierRemoval<PointXYZ> ror(true);
    ror.setNegative(false);             //false: Exclude points found; true: Keeps only points found (should be false)
    ror.setRadiusSearch(1.2*sqrt(2.0)); //Max search radius. Stops NN search when distance to neighbor exceeds this radius. Should be > max KNN distance between point (test data is about ~1.6)
    ror.setMinNeighborsInRadius(8);     //Minimum neighbors found for keeping. In planes: ~1 excludes outliers only, ~3-4 excludes corners, ~8-9 excludes edge points, more than ~9-10 deletes entire planes



    vector<vector<PointCloud<PointXYZ>>> point_clusters;
    int cloud_index = 0;
    int cluster_count = 0;
    while (true) {


        #ifdef DEBUG
            cout << "Iteration " << cloud_index + 1 << ": Fitting plane..." << endl;
        #endif

        seg.segment(*inlier_indices, *coefficients);
        vector<PointIndices> clusters;       
        
        if (inlier_indices->indices.size() > 0) { //If plane found...
            planes.push_back(*coefficients); //Save coefficient to list of planes

            #ifdef DEBUG
                cout << "\tPlane found with coefficients " << planes[cloud_index].values[0] << ", " << planes[cloud_index].values[1] << ", " << planes[cloud_index].values[2] << ", " << planes[cloud_index].values[3] << ", copying points" << endl;
            #endif // DEBUG

            copyPointCloud(*beam, inlier_indices->indices, *inliers); //Copy inliers to new cloud

            #ifdef DEBUG
                filename = "test_ransac_";
                filename += string(to_string(cloud_index));
                filename += ".pcd";
                cout << "\tSaving " << filename << " with " << inlier_indices->indices.size() << " points" << endl;
                io::savePCDFileBinary(filename, *inliers);

                //copy inliers to new cloud for saving removed outliers
                PointCloud<PointXYZ>::Ptr complete_plane(new PointCloud<PointXYZ>);
                copyPointCloud(*inliers, *complete_plane); 
            #endif

            //KNN search keeps only central points in cluster
            ror.setInputCloud(inliers);
            ror.filter(*inliers);

            #ifdef DEBUG
                PointCloud<PointXYZ>::Ptr removed_points(new PointCloud<PointXYZ>);
                auto removed_indices = ror.getRemovedIndices();
                copyPointCloud(*complete_plane, *removed_indices, *removed_points);
                cout << "\tKNN search found " << removed_indices->size() << " points" << endl;
            #endif        
                                                                
            //Cluster remaining points         
            ece.setInputCloud(inliers);
            ece.extract(clusters);

            #ifdef DEBUG
                cout << "\tClustering done, " << clusters.size() << " cluster(s) found" << endl;
            #endif      

            //Point clusters stored in 2D vector indexed with (plane, cluster)
            vector<PointCloud<PointXYZ>> current_clusters(0);
            point_clusters.push_back(current_clusters);
            
            #ifdef DEBUG
                cout << "\t\tSaving:" << endl;
            #endif // DEBUG


            for (int i = 0; i < clusters.size(); i++) {
                //todo: Excessive copying, replace with make_shared/etc?
                PointCloud<PointXYZ>::Ptr cluster(new PointCloud<PointXYZ>);
                PointIndices::Ptr indices(new PointIndices);
                indices->header = clusters[i].header;
                indices->indices = clusters[i].indices;

                copyPointCloud(*inliers, clusters[i].indices, *cluster);
                point_clusters[cloud_index].push_back(*cluster);
                cluster_count++;

                #ifdef DEBUG
                    filename = "test_cluster_";
                    filename += string(to_string(cloud_index));
                    filename += string(to_string(i));
                    filename += ".pcd";
                    cout << "\t\t" << filename << endl;
                    io::savePCDFileBinary(filename, *cluster);
                #endif
            }

            #ifdef DEBUG
                cout << "\tDeleting RANSAC points from dataset" << endl;
            #endif

            //Remove spent points from dataset
            extract.setIndices(inlier_indices);
            extract.setNegative(true);
            extract.filter(*beam);
        }

        #ifndef DEBUG   
                cout << planes.size() << " plane(s) found with " << cluster_count << " cluster(s)" << endl;
        #endif

        if (beam->size() < CLUSTER_MIN_SIZE ) {//If there aren't enough points to make a cluster, we're done
            cout << "RANSAC deconstruction done!" << endl;
            break;
        }
        cloud_index++;
    }



    cout << "Building list of corner points... ";
    vector<vector<PointCloud<PointXYZ>>> face_points(point_clusters.size());
    for (int i = 0; i < point_clusters.size(); i++) {
        face_points[i] = vector<PointCloud<PointXYZ>>(point_clusters[i].size());
    }

    vector<PointXYZ> corner_points = find_corner_points(planes, bbox);
    cout << corner_points.size() << " possible mesh edge points found." << endl;

    cout << "Grouping mesh edges by cluster... ";
    for (int i = 0; i < point_clusters.size(); i++ ) {
        for (int j = 0; j < point_clusters[i].size(); j++) {
            for (int k = 0; k < corner_points.size(); k++) {
                if (closest_point_dist(corner_points[k], point_clusters[i][j]) < CLUSTER_MESH_DISTANCE_THRESHOLD) {
                    face_points[i][j].push_back(corner_points[k]);
                }
            }
        }
    }
    cout << "done!" << endl;



    cout << "Building mesh faces... ";
    vector < PolygonMesh > mesh;
    for (int i = 0; i < face_points.size(); i++) {
        for (int j = 0; j < face_points[i].size(); j++) {
            #ifdef DEBUG
                cout << "Face " << i << "." << j << ": " << face_points[i][j].size() << " points" << endl;
                filename = "test_face_";
                filename += string(to_string(i));
                filename += string(to_string(j));
                filename += ".pcd";
                cout << "\t\tSaving " << filename << endl;
                io::savePCDFileBinary(filename, face_points[i][j]);
            #endif

            PointCloud<PointXYZ>::Ptr face(new PointCloud<PointXYZ>);
            copyPointCloud(face_points[i][j], *face);

            #ifdef DEBUG
                cout << "\t\tTriangulating mesh face " << i << "." << j << endl;
            #endif

            search::KdTree<PointXYZ>::Ptr mesh_tree(new search::KdTree<PointXYZ>);
            mesh_tree->setInputCloud(face);
                
            ConvexHull<PointXYZ> chul;
            chul.setInputCloud(face);
            chul.setSearchMethod(mesh_tree);
            chul.setDimension(2);
            PolygonMesh triangles;
            chul.reconstruct(triangles);
            mesh.push_back(triangles);


            #ifdef DEBUG
                //Write mesh to file
                filename = "mesh_";
                filename += string(to_string(i));
                filename += string(to_string(j));
                filename += ".vtk";
                cout << "\t\tSaving " << filename << endl << endl;
                io::saveVTKFile(filename, triangles);
            #endif                     
        }
    }
    cout << "done!" << endl;
    
    cout << "Merging mesh faces... ";
    while (mesh.size() > 1) {
        //Take two meshes off
        PolygonMesh mesh_1 = mesh.back();
        mesh.pop_back();

        PolygonMesh mesh_2 = mesh.back();
        mesh.pop_back();

        //Concatenate them
        PolygonMesh mesh_merged = mesh_1 + mesh_2; //Concat fields, moves and duplicates some data
        
        //Move the concat'd mesh vertices to their correct location
        int poly_index = mesh_merged.polygons.size() - 1;
        for (int j = 0; j < mesh_merged.polygons[poly_index].vertices.size(); j++) {
            mesh_merged.polygons[poly_index].vertices[j] -= mesh_2.cloud.width;
        }

        //Copy mesh vertices to cloud
        PointCloud<PointXYZ> merged_points;
        fromPCLPointCloud2(mesh_merged.cloud, merged_points);
        for (int i = merged_points.size()-1; i >= 0; i--) {
            for (int j = i-1; j >= 0; j--) {
                if ( (abs(merged_points[i].x - merged_points[j].x) < EPSILON) && //If two points match...
                     (abs(merged_points[i].y - merged_points[j].y) < EPSILON) &&
                     (abs(merged_points[i].z - merged_points[j].z) < EPSILON) ){

                    poly_index = 0; //Iterate through mesh until one of the matching pair is found
                    for (int k = 0; k < mesh_merged.polygons.size(); k++) {
                        for (int w = 0; w < mesh_merged.polygons[k].vertices.size(); w++) {
                            if (poly_index == i) {
                                if (mesh_merged.polygons[k].vertices[w] > j) { //Always pick the lowest index vertex to connect
                                    mesh_merged.polygons[k].vertices[w] = j;
                                }
                            }
                            poly_index++;
                        }
                    }

                    //todo: figure out how to delete the now-unused verteces

                }
            }
        }

        mesh.push_back(mesh_merged);   
    }
    cout << "done!" << endl;

    #ifdef DEBUG
        cout << "Saving various formats for debugging" << endl;
        io::savePolygonFile("mesh.vtk", mesh[0], false);
        io::savePolygonFile("mesh.stl", mesh[0], true);
        io::savePolygonFile("mesh.ply", mesh[0], true);
    #endif

    cout << "Saving output: mesh.stl" << endl;
    io::savePolygonFileSTL("mesh.stl", mesh[0], true);

    return (0);
}
