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
    bool isShadow;
} ray;

typedef struct Hit{
    bool isHit;
    float t;
    Vec3f normal_vector;
    Vec3f intersection_point;
    int material_id;
    int object_id;
    // 0-> Triangle
    // 1-> Sphere
    // 2-> Mesh
    int object_type;
} hit;

Vec3f normalization(const Vec3f &vector){
    float length = sqrt(vector.x*vector.x + vector.y*vector.y + vector.z*vector.z);
    Vec3f normalized_vector = {vector.x/length, vector.y/length, vector.z/length};
    
    return normalized_vector;
}

Vec3f crossProduct(const Vec3f &vector1, const Vec3f &vector2){
    Vec3f cross_product = {vector1.y*vector2.z - vector1.z*vector2.y, vector1.z*vector2.x - vector1.x*vector2.z, vector1.x*vector2.y - vector1.y*vector2.x};
    return cross_product;
}

Vec3f subtractVectors(const Vec3f &vector1, const Vec3f &vector2){
    Vec3f subtracted_vector;
    subtracted_vector.x = vector1.x-vector2.x;
    subtracted_vector.y = vector1.y-vector2.y;
    subtracted_vector.z = vector1.z-vector2.z;
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
    ray.isShadow = false;
    return ray;
}

vector<Vec3f> calculateTrianglesNormalVectors(const Scene &scene){
    vector<Vec3f> triangle_normal_vectors;
    int number_of_triangles = scene.triangles.size();
    Triangle current_triangle;
    Vec3f v0, v1, v2;
    for(int triangle_index=0; triangle_index<number_of_triangles; triangle_index++){
        current_triangle = scene.triangles[triangle_index];
        v0 = scene.vertex_data[current_triangle.indices.v0_id - 1];
        v1 = scene.vertex_data[current_triangle.indices.v1_id - 1];
        v2 = scene.vertex_data[current_triangle.indices.v2_id - 1];
        Vec3f normal_vector = normalization(crossProduct(subtractVectors(v1, v0), subtractVectors(v2, v0)));
        triangle_normal_vectors.push_back(normal_vector);
    }
    return triangle_normal_vectors;
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

float dotProduct(const Vec3f vector1, const Vec3f vector2){
    float dot_product = vector1.x*vector2.x + vector1.y*vector2.y + vector1.z*vector2.z;
    return dot_product;
}

vector<float> findRoots(float A, float B, float C){
    vector<float> roots;
    float discriminant = B*B - 4*A*C;
    roots.push_back(discriminant);
    if(discriminant >= 0){
        float t1 = (-B + sqrtf(discriminant))/(2*A);
        float t2 = (-B - sqrtf(discriminant))/(2*A);
        roots.push_back(t1);
        roots.push_back(t2);
    }
    return roots;
}

float calculateDeterminant(const Vec3f &vector1, const Vec3f &vector2, const Vec3f &vector3){
    
    float determinant = (vector1.x * (vector2.y*vector3.z - vector3.y*vector2.z)) + (vector1.y * (vector3.x*vector2.z - vector2.x*vector3.z)) + (vector1.z * (vector2.x*vector3.y - vector2.y*vector3.x));    
    return determinant;
}

Hit findClosestHit(vector<Hit> &all_hits){
    Hit closest_hit;

    int number_of_hits = all_hits.size();
    if(number_of_hits > 0){
        closest_hit.isHit = true;
        closest_hit = all_hits[0];
        for(int hit_index=1; hit_index<number_of_hits; hit_index++){
            if(all_hits[hit_index].t < closest_hit.t){
                closest_hit = all_hits[hit_index];
            }
        }
    }
    else{
        closest_hit.isHit = false;
    }
    
    return closest_hit;
}

Hit operateHit(const Scene &scene, const Ray &ray, vector<Vec3f> &triangle_normal_vectors, vector<vector<Vec3f>> &mesh_normal_vectors){
    Triangle current_triangle;
    Mesh current_mesh;
    Sphere current_sphere;
    Vec3f normal;
    Vec3f intersection_point;

    int number_of_spheres = scene.spheres.size();
    int number_of_triangles = scene.triangles.size();
    int number_of_meshes = scene.meshes.size();

    vector<Hit> all_hits;

    // Check hit for each triangle
    for(int triangle_index = 0; triangle_index < number_of_triangles; triangle_index++){
        current_triangle = scene.triangles[triangle_index];
        Vec3f v0 = scene.vertex_data[current_triangle.indices.v0_id - 1];
        Vec3f v1 = scene.vertex_data[current_triangle.indices.v1_id - 1];
        Vec3f v2 = scene.vertex_data[current_triangle.indices.v2_id - 1];

        Hit tri_hit;
        tri_hit.isHit = true;

        Vec3f o = ray.origin;
        Vec3f d = ray.direction;

        Vec3f a_minus_b = subtractVectors(v0, v1);
        Vec3f a_minus_c = subtractVectors(v0, v2);
        Vec3f a_minus_o = subtractVectors(v0, o);

        float t;
        float detA = calculateDeterminant(a_minus_b, a_minus_c, d);

        
        if(detA != 0.0){
            t = calculateDeterminant(a_minus_b, a_minus_c, a_minus_o) / detA;
            if(t>0){
                float beta = calculateDeterminant(a_minus_o, a_minus_c, d) / detA;
                if(beta >= 0 && beta <= 1){
                    float gamma = calculateDeterminant(a_minus_b, a_minus_o, d) / detA;
                    if(gamma >= 0 && gamma <= 1-beta){
                        tri_hit.t = t;
                        tri_hit.intersection_point.x = o.x + t*d.x;
                        tri_hit.intersection_point.y = o.y + t*d.y;
                        tri_hit.intersection_point.z = o.z + t*d.z;
                        //tri_hit.normal_vector = triangle_normal_vectors[triangle_index];
                        tri_hit.normal_vector = normalization(crossProduct(subtractVectors(v1, v0), subtractVectors(v2, v0)));
                        tri_hit.material_id = current_triangle.material_id;
                        tri_hit.object_id = triangle_index;
                        tri_hit.object_type = 0;
                        if((ray.isShadow == true)&& (tri_hit.t > 0) && (tri_hit.t < 1)){
                            return tri_hit;
                        }
                    }
                    else{
                        tri_hit.isHit=false;
                    }
                }
                else{
                    tri_hit.isHit=false;
                }
            }
            else{
                tri_hit.isHit=false;
            }
        }

        else{
            tri_hit.isHit=false;
        }

        if(tri_hit.isHit && tri_hit.t >= 0){
            all_hits.push_back(tri_hit);
        }
        
    }


    // Check hit for each sphere
    for(int sphere_index = 0; sphere_index < number_of_spheres; sphere_index++){
        current_sphere = scene.spheres[sphere_index];
        Vec3f center = scene.vertex_data[current_sphere.center_vertex_id - 1];
        float radius = current_sphere.radius;

        Hit sphere_hit;

        const float A = dotProduct(ray.direction, ray.direction);
        Vec3f e_minus_c = subtractVectors(ray.origin, center);
        const float B = 2 * dotProduct(ray.direction, e_minus_c);
        const float C = dotProduct(e_minus_c, e_minus_c) - radius * radius;

        vector<float> roots = findRoots(A, B, C);
        if(roots.size() > 1){	
            sphere_hit.isHit = true;

            const float t1 = roots[1];
            const float t2 = roots[2];

            sphere_hit.material_id = current_sphere.material_id;
            
            //hit.objType = SPHERE;
            sphere_hit.object_id = sphere_index;

            const float t = (t1 < t2) ? t1 : t2;
            sphere_hit.intersection_point.x = ray.origin.x + t*ray.direction.x;
            sphere_hit.intersection_point.y = ray.origin.y + t*ray.direction.y;
            sphere_hit.intersection_point.z = ray.origin.z + t*ray.direction.z;
            sphere_hit.normal_vector = normalization(subtractVectors(sphere_hit.intersection_point, center)); //can be optimized with radius
            sphere_hit.object_type = 1;
            sphere_hit.t = t;
            if((ray.isShadow == true) && (sphere_hit.t > 0) && (sphere_hit.t < 1)){
                return sphere_hit;
            }
            
        }
        else{
            sphere_hit.isHit = false;
        }
        if(sphere_hit.isHit && sphere_hit.t >= 0){
            all_hits.push_back(sphere_hit);
        }
        
    }



    // Check hit for each mesh
    for(int mesh_index = 0; mesh_index < number_of_meshes; mesh_index++){
        current_mesh = scene.meshes[mesh_index];

        Hit mesh_hit;
	    mesh_hit.isHit = false;
        vector<Hit> all_mesh_hits;

        // FOR EACH FACE IN Current Mesh
        for(int face_index = 0; face_index < current_mesh.faces.size(); face_index++){

            Vec3f v0 = scene.vertex_data[current_mesh.faces[face_index].v0_id - 1];	
            Vec3f v1 = scene.vertex_data[current_mesh.faces[face_index].v1_id - 1];
            Vec3f v2 = scene.vertex_data[current_mesh.faces[face_index].v2_id - 1];

            Hit face_hit;
            face_hit.isHit = true;

            Vec3f o = ray.origin;
            Vec3f d = ray.direction;

            Vec3f a_minus_b = subtractVectors(v0, v1);
            Vec3f a_minus_c = subtractVectors(v0, v2);
            Vec3f a_minus_o = subtractVectors(v0, o);

            float t;

            float detA = calculateDeterminant(a_minus_b, a_minus_c, d);
            

            if(detA != 0.0){
                t = calculateDeterminant(a_minus_b, a_minus_c, a_minus_o) / detA;
                if(t>0){
                    float beta = calculateDeterminant(a_minus_o, a_minus_c, d) / detA;
                    if(beta >= 0 && beta <= 1){
                        float gamma = calculateDeterminant(a_minus_b, a_minus_o, d) / detA;
                        if(gamma >= 0 && gamma <= 1-beta){
                            face_hit.t = t;
                            face_hit.intersection_point.x = o.x + t*d.x;
                            face_hit.intersection_point.y = o.y + t*d.y;
                            face_hit.intersection_point.z = o.z + t*d.z;
                            //face_hit.normal_vector = mesh_normal_vectors[mesh_index][face_index];
                            face_hit.normal_vector = normalization(crossProduct(subtractVectors(v1, v0), subtractVectors(v2, v0)));
                            face_hit.material_id = current_triangle.material_id;
                            face_hit.object_id = mesh_index;
                            face_hit.object_type = 2;
                            if((ray.isShadow == true) && (face_hit.t > 0) && (face_hit.t < 1)){
                                return face_hit;
                            }
                            
                        }
                        else{
                            face_hit.isHit=false;
                        }
                    }
                    else{
                        face_hit.isHit=false;
                    }
                }
                else{
                    face_hit.isHit=false;
                }
            }

            else{
                face_hit.isHit=false;
            }
            if(face_hit.isHit && face_hit.t >= 0){
                face_hit.material_id = current_mesh.material_id;
                //face_hit.objType = MESH;
                face_hit.object_id = mesh_index;
                face_hit.intersection_point.x = o.x + t*d.x;
                face_hit.intersection_point.y = o.y + t*d.y;
                face_hit.intersection_point.z = o.z + t*d.z;
                //face_hit.normal_vector = mesh_normal_vectors[mesh_index][face_index];
                face_hit.normal_vector = normalization(crossProduct(subtractVectors(v1, v0), subtractVectors(v2, v0)));
                all_mesh_hits.push_back(face_hit);
            }
            
        }
        mesh_hit = findClosestHit(all_mesh_hits);

        if(mesh_hit.isHit && mesh_hit.t >= 0){
            all_hits.push_back(mesh_hit);
        }
    }

    Hit hitResult = findClosestHit(all_hits);

    return hitResult;
}

Vec3f computeColor(const Scene &scene, const Ray &ray, const Hit &hit, int max_recursion_depth, vector<Vec3f>&triangle_normal_vectors, vector<vector<Vec3f>>&mesh_normal_vectors){
    Vec3f pixel_value;
    // add ambient
    pixel_value.x = 0;
    pixel_value.y = 0;
    pixel_value.z = 0;

    if(hit.isHit){
        // add ambient 
        pixel_value.x = scene.materials[hit.material_id - 1].ambient.x * scene.ambient_light.x;
        pixel_value.y = scene.materials[hit.material_id - 1].ambient.y * scene.ambient_light.y;
        pixel_value.z = scene.materials[hit.material_id - 1].ambient.z * scene.ambient_light.z;
        //create ray from object hit point+epsilon to the all lights, determine if shadow occurs
            // if shadow true, pixel_value = ambient
            // if shadow false, pixel_value = ambient + diffuse + specular 
        Vec3f hit_point = hit.intersection_point;
        Vec3f normal_vector = hit.normal_vector;
        // add epsilon to hit point in normal vector direction
        hit_point.x += scene.shadow_ray_epsilon*normal_vector.x;
        hit_point.y += scene.shadow_ray_epsilon*normal_vector.y;
        hit_point.z += scene.shadow_ray_epsilon*normal_vector.z;

        Vec3f reflection;
        reflection.x = 0;
        reflection.y = 0;
        reflection.z = 0;
        if(max_recursion_depth>0 && scene.materials[hit.material_id - 1].is_mirror){
            float dot_product_wi = -2 * dotProduct(normal_vector, ray.direction);
            Vec3f normal_wi;
            normal_wi.x = dot_product_wi * normal_vector.x + ray.direction.x;
            normal_wi.y = dot_product_wi * normal_vector.y + ray.direction.y;
            normal_wi.z = dot_product_wi * normal_vector.z + ray.direction.z;
            normal_wi = normalization(normal_wi);

            Vec3f epsilon;
            epsilon.x = scene.shadow_ray_epsilon * normal_wi.x;
            epsilon.y = scene.shadow_ray_epsilon * normal_wi.y;
            epsilon.z = scene.shadow_ray_epsilon * normal_wi.z;

            Ray reflection_ray;

            hit_point.x += epsilon.x;
            hit_point.y += epsilon.y;
            hit_point.z += epsilon.z;

            reflection_ray.origin = hit_point;
            reflection_ray.direction = normal_wi;
            reflection_ray.isShadow = false;

            Hit reflection_hit = operateHit(scene, reflection_ray, triangle_normal_vectors, mesh_normal_vectors);

            if(reflection_hit.isHit){
                reflection = computeColor(scene, reflection_ray, reflection_hit, max_recursion_depth-1, triangle_normal_vectors, mesh_normal_vectors);
            }

            pixel_value.x += scene.materials[hit.material_id - 1].mirror.x * reflection.x;
            pixel_value.y += scene.materials[hit.material_id - 1].mirror.y * reflection.y;
            pixel_value.z += scene.materials[hit.material_id - 1].mirror.z * reflection.z;

        }
        // for each point lights
        for(int light_index=0; light_index<scene.point_lights.size(); light_index++){
            PointLight current_light = scene.point_lights[light_index];
            Vec3f light_position = current_light.position;
            Vec3f light_intensity = current_light.intensity;
            Vec3f light_direction = subtractVectors(light_position, hit_point);
            Ray shadow_ray;
            shadow_ray.origin = hit_point;
            shadow_ray.direction = light_direction;
            shadow_ray.isShadow = true;
            Hit shadow_hit = operateHit(scene, shadow_ray, triangle_normal_vectors, mesh_normal_vectors);
            if( (shadow_hit.isHit && (shadow_hit.t <= 0 && shadow_hit.t >= 1)) || !shadow_hit.isHit){
                //diffuse
                float distance_square = dotProduct(light_direction, light_direction);
                Vec3f irradiance = {light_intensity.x/distance_square, light_intensity.y/distance_square, light_intensity.z/distance_square};
                light_direction = normalization(light_direction);
                float dot_product = dotProduct(normal_vector, light_direction);
                if(dot_product > 0){
                    pixel_value.x += scene.materials[hit.material_id - 1].diffuse.x * irradiance.x * dot_product;
                    pixel_value.y += scene.materials[hit.material_id - 1].diffuse.y * irradiance.y * dot_product;
                    pixel_value.z += scene.materials[hit.material_id - 1].diffuse.z * irradiance.z * dot_product;
                }

                //specular
                Vec3f vector_h = normalization(subtractVectors(light_direction, ray.direction));
                float dot_product_hn = dotProduct(normal_vector, vector_h);
                if(dot_product_hn > 0){
                    pixel_value.x += scene.materials[hit.material_id - 1].specular.x * irradiance.x * powf(dot_product_hn, scene.materials[hit.material_id - 1].phong_exponent);
                    pixel_value.y += scene.materials[hit.material_id - 1].specular.y * irradiance.y * powf(dot_product_hn, scene.materials[hit.material_id - 1].phong_exponent);
                    pixel_value.z += scene.materials[hit.material_id - 1].specular.z * irradiance.z * powf(dot_product_hn, scene.materials[hit.material_id - 1].phong_exponent);
                }
                
                
            }
            
        }
        
        
    }
    // first recursion
    else if(max_recursion_depth == scene.max_recursion_depth){
        // primary ray, no hit
        pixel_value.x = scene.background_color.x;
        pixel_value.y = scene.background_color.y;
        pixel_value.z = scene.background_color.z;
    }
    else{
        return {0,0,0};
    }
    return pixel_value;
}
int main(int argc, char* argv[])
{
    // Sample usage for reading an XML scene file
    parser::Scene scene;

    scene.loadFromXml(argv[1]);

    int camera_number = scene.cameras.size();
    vector<Vec3f> triangle_normal_vectors = calculateTrianglesNormalVectors(scene);
    vector<vector<Vec3f>> mesh_normal_vectors = calculateMeshesNormalVectors(scene);
 
    // FOR EACH CAMERA
    for(int camera_index=0; camera_index<camera_number; camera_index++){
        
        Camera camera = scene.cameras[camera_index];

        float left = camera.near_plane.x;
        float right = camera.near_plane.y;
        float bottom = camera.near_plane.z;
        float top = camera.near_plane.w;
        Vec3f normalized_gaze = normalization(camera.gaze);
        Vec3f normalized_camera_v = normalization(camera.up);
        Vec3f normalized_camera_u = crossProduct(normalized_gaze, normalized_camera_v);
        Vec3f image_top_left;
        Vec3f image_plane_center;

        image_plane_center.x=camera.position.x + (normalized_gaze.x*camera.near_distance);
        image_plane_center.y=camera.position.y + (normalized_gaze.y*camera.near_distance);
        image_plane_center.z=camera.position.z + (normalized_gaze.z*camera.near_distance);

        image_top_left.x=image_plane_center.x + (normalized_camera_u.x*left) + (normalized_camera_v.x*top);
        image_top_left.y=image_plane_center.y + (normalized_camera_u.y*left) + (normalized_camera_v.y*top);
        image_top_left.z=image_plane_center.z + (normalized_camera_u.z*left) + (normalized_camera_v.z*top);

        float px = (right-left)/camera.image_width;
        float py = (top-bottom)/camera.image_height;
        

        //create image for current camera
        unsigned char* image = new unsigned char [camera.image_width * camera.image_height * 3];

        int pixel_index = 0;
        Vec3f pixel_value;
        // FOR EACH PIXEL
        for(int i=0; i<camera.image_height; i++){
            for(int j=0; j<camera.image_width; j++){
                Ray ray = createRay(i, j, camera, image_top_left, px, py, normalized_camera_u, normalized_camera_v);
                Hit hit = operateHit(scene, ray, triangle_normal_vectors, mesh_normal_vectors);


                pixel_value = computeColor(scene, ray, hit, scene.max_recursion_depth,triangle_normal_vectors, mesh_normal_vectors);
                
                image[pixel_index] = (pixel_value.x > 255) ? 255 : round(pixel_value.x);
                image[pixel_index+1] = (pixel_value.y > 255) ? 255 : round(pixel_value.y);
                image[pixel_index+2] = (pixel_value.z > 255) ? 255 : round(pixel_value.z);
                pixel_index += 3;
            }
        }


        write_ppm(camera.image_name.c_str(),image, camera.image_width, camera.image_height);
    }
    

}
