#include <iostream>
#include "parser.h"
#include "ppm.h"
#include <math.h>


using namespace std;
using namespace parser;

typedef unsigned char RGB[3];

typedef struct Ray{
    Vec3f origin;
    Vec3f direction;
} ray;

typedef struct Hit{
    bool isHit;
    float t;
    Vec3f normal_vector;
    Vec3f intersection_point;
    int material_id;
} hit;

Vec3f normalization(const Vec3f vector){
    float length = sqrt(vector.x*vector.x + vector.y*vector.y + vector.z*vector.z);
    Vec3f normalized_vector = {vector.x/length, vector.y/length, vector.z/length};
    return normalized_vector;
}

Vec3f crossProduct(const Vec3f vector1, const Vec3f vector2){
    Vec3f cross_product = {vector1.y*vector2.z - vector1.z*vector2.y, vector1.z*vector2.x - vector1.x*vector2.z, vector1.x*vector2.y - vector1.y*vector2.x};
    return cross_product;
}

Vec3f subtractVectors(const Vec3f vector1, const Vec3f vector2){
    Vec3f subtracted_vector = {vector1.x-vector2.x, vector1.y-vector2.y, vector1.z-vector2.z};
    return subtracted_vector;
}

Ray createRay(int i, int j, const Camera &camera, Vec3f image_top_left, float px, float py, Vec3f normalized_camera_u, Vec3f normalized_camera_v){
    Ray ray;
    Vec3f pixel_position;
    
    float u_offset = px*(j+0.5);
    float v_offset = py*(i+0.5);

    pixel_position.x = image_top_left.x + (normalized_camera_u.x*u_offset) - (normalized_camera_v.x*v_offset);
    pixel_position.y = image_top_left.y + (normalized_camera_u.y*u_offset) - (normalized_camera_v.y*v_offset);
    pixel_position.z = image_top_left.z + (normalized_camera_u.z*u_offset) - (normalized_camera_v.z*v_offset);

    ray.origin = camera.position;
    ray.direction = normalization(subtractVectors(pixel_position, camera.position));
    
    return ray;
}

vector<Vec3f> calculateTrianglesNormalVectors(const Scene &scene){
    vector<Vec3f> triange_normal_vectors;
    int number_of_triangles = scene.triangles.size();
    Triangle current_triangle;
    Vec3f v0, v1, v2;
    for(int triangle_index=0; triangle_index<number_of_triangles; triangle_index++){
        current_triangle = scene.triangles[triangle_index];
        v0 = scene.vertex_data[current_triangle.indices.v0_id - 1];
        v1 = scene.vertex_data[current_triangle.indices.v1_id - 1];
        v2 = scene.vertex_data[current_triangle.indices.v2_id - 1];
        Vec3f normal_vector = normalization(crossProduct(subtractVectors(v1, v0), subtractVectors(v2, v0)));
        triange_normal_vectors.push_back(normal_vector);
    }
    return triange_normal_vectors;
}

vector<vector<Vec3f>> calculateMeshesNormalVectors(const Scene &scene){
    vector<vector<Vec3f>> mesh_normal_vectors;
    vector<Vec3f> face_normal_vector;
    int number_of_meshes = scene.meshes.size();
    Mesh current_mesh;
    Vec3f v0, v1, v2;
    // FOR EACH MESH
    for(int mesh_index=0; mesh_index<number_of_meshes; mesh_index++){
        current_mesh = scene.meshes[mesh_index];
        int number_of_faces = current_mesh.faces.size();
        // FOR EACH FACE IN Current Mesh
        for(int face_index=0; face_index<number_of_faces; face_index++){
            Face current_face = current_mesh.faces[face_index];
            v0 = scene.vertex_data[current_face.v0_id - 1];
            v1 = scene.vertex_data[current_face.v1_id - 1];
            v2 = scene.vertex_data[current_face.v2_id - 1];
            Vec3f normal_vector = normalization(crossProduct(subtractVectors(v1, v0), subtractVectors(v2, v0)));
            face_normal_vector.push_back(normal_vector);
        }
        mesh_normal_vectors.push_back(face_normal_vector);
    }
    return mesh_normal_vectors;
}


int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);
    //std::cout<<scene.cameras[0].position.x<<std::endl;

    int camera_number = scene.cameras.size();
    vector<Vec3f> triangle_normal_vectors = calculateTrianglesNormalVectors(scene);
    vector<vector<Vec3f>> mesh_normal_vectors = calculateMeshesNormalVectors(scene);

    // FOR EACH CAMERA
    for(int camera_index=0; camera_index<camera_number; camera_index++){
        
        Camera camera = scene.cameras[camera_index];

        float left = camera.near_plane.x;
        float right = camera.near_plane.y;
        float top = camera.near_plane.z;
        float bottom = camera.near_plane.w;
        Vec3f normalized_gaze = normalization(camera.gaze);
        Vec3f normalized_camera_v = normalization(camera.up);
        Vec3f normalized_camera_u = crossProduct(normalized_gaze, normalized_camera_v);
        Vec3f image_top_left;
        Vec3f image_plane_center;

        image_plane_center.x=camera.position.x + (normalized_camera_u.x*left) + (normalized_camera_v.x*top);
        image_plane_center.y=camera.position.y + (normalized_camera_u.y*left) + (normalized_camera_v.y*top);
        image_plane_center.z=camera.position.z + (normalized_camera_u.z*left) + (normalized_camera_v.z*top);

        image_top_left.x=image_plane_center.x + (normalized_camera_u.x*left) + (normalized_camera_v.x*top);
        image_top_left.y=image_plane_center.y + (normalized_camera_u.y*left) + (normalized_camera_v.y*top);
        image_top_left.z=image_plane_center.z + (normalized_camera_u.z*left) + (normalized_camera_v.z*top);

        float px = (right-left)/camera.image_width;
        float py = (top-bottom)/camera.image_height;
        

        //create image for current camera
        unsigned char* image = new unsigned char [camera.image_width * camera.image_height * 3];

        int pixel_index = 0;

        // FOR EACH PIXEL
        for(int i=0; i<camera.image_height; i++){
            for(int j=0; j<camera.image_width; j++){
                Ray ray = createRay(i, j, camera, image_top_left, px, py, normalized_camera_u, normalized_camera_v);
                
            }
        }



    }
    

}
