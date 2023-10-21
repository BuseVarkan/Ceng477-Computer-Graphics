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

Vec3f subtract_vectors(const Vec3f vector1, const Vec3f vector2){
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
    ray.direction = normalization(subtract_vectors(pixel_position, camera.position));
    
    return ray;
}

int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);
    //std::cout<<scene.cameras[0].position.x<<std::endl;

    int camera_number = scene.cameras.size();

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
