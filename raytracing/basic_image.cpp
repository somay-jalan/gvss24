#define STB_IMAGE_IMPLEMENTATION
#include "../src/color.h"
#include "../src/vec3.h"
#include "../src/ray.h"
#include "../src/sphere.hpp"
#include "../src/plane.hpp"
#include "../src/pointLight.hpp"
#include "../src/stb_image.h"

#include<bits/stdc++.h>
#include <iostream>

using namespace std;
//THIS SECTION OF CODE HAS BEEN TAKEN WITH THE HELP OF GPT
//.--------------------------------------------------------------------------
int texture_width, texture_height, channels;
unsigned char *texture_image = nullptr;

void load_texture() {
    // Load texture image from file
    texture_image = stbi_load("../texture.png", &texture_width, &texture_height, &channels, 0);
    if (!texture_image) {
        cerr << "Failed to load texture image." << endl;
        exit(EXIT_FAILURE);
    }
}
color get_texture_color(double u, double v) {
    if (!texture_image) {
        cerr << "Texture image not loaded." << endl;
        return color(1, 1, 1); // Default color if texture not loaded
    }

    int i = min(max(int(u * texture_width), 0), texture_width - 1);
    int j = min(max(int((1 - v) * texture_height), 0), texture_height - 1); // Flip v to match image coordinates

    // Assuming each pixel is 3 channels (RGB)
    int pixel_index = (j * texture_width + i) * channels;
    double r = texture_image[pixel_index] / 255.0;
    double g = texture_image[pixel_index + 1] / 255.0;
    double b = texture_image[pixel_index + 2] / 255.0;

    return color(r, g, b);
}

//.--------------------------------------------------------------------------


bool shadow_check(const ray& r,
    vector<PointLight> &point_light,
    vector<Sphere> &sphere_light,
    vector<Plane> &plane_light,
    vector<Plane> &plane_object,
    vector<Sphere> &sphere_object,
    double t
    ){
    for(int i=0;i<point_light.size();i++){
        double temp_t = point_light[i].hit_PointLight(r);
        // if(temp_t!=INT_MAX) cout<<"i:"<<i<<" temp_t"<<temp_t<<"||";
        if(temp_t<t and temp_t>0){
            return false;
        }
    }
    for(int i=0;i<sphere_light.size();i++){
        PointLight tempLight = PointLight(sphere_light[i].getCenter(),sphere_light[i].getObjectColor());
        double temp_t = tempLight.hit_PointLight(r);
        // if(temp_t!=INT_MAX) cout<<"i:"<<i<<" temp_t"<<temp_t<<"||";
        if(temp_t<t and temp_t>0){
            return false;
        }
    }

    for(int i=0;i<plane_light.size();i++){
        PointLight tempLight = PointLight(plane_light[i].getCenter(),plane_light[i].getObjectColor());
        double temp_t = tempLight.hit_PointLight(r);
        // if(temp_t!=INT_MAX) cout<<"i:"<<i<<" temp_t"<<temp_t<<"||";
        if(temp_t<t and temp_t>0){
            return false;
        }
    }
    for(int i=0;i<plane_object.size();i++){
        double temp_t = plane_object[i].hit_plane(r);
        // if(temp_t!=INT_MAX) cout<<"i:"<<i<<" temp_t"<<temp_t<<"||";
        if(temp_t<t and temp_t>0){
            return false;
        }
    }
    for(int i=0;i<sphere_object.size();i++){
        double temp_t = sphere_object[i].hit_sphere(r);
        // if(temp_t!=INT_MAX) cout<<"i:"<<i<<" temp_t"<<temp_t<<"||";
        if(temp_t<t and temp_t>0){
            return false;
        }
    }
    return true;
}

color ray_color(const ray& r,
    vector<PointLight> &point_light,
    vector<Sphere> &sphere_light,
    vector<Plane> &plane_light,
    vector<Plane> &plane_object,
    vector<Sphere> &sphere_object,
    color &ambient_Color,
    int recDepth,
    int maxDepth) {
    
    Sphere hitobjectSphere;
    Plane hitobjectPlane;
    PointLight hitLightPoint;
    string check="None";
    double t = 1e8;;
    for(int i=0;i<sphere_object.size();i++){
        double temp_t = sphere_object[i].hit_sphere(r);
        // if(temp_t!=INT_MAX) cout<<"i:"<<i<<" temp_t"<<temp_t<<"||";
        if(temp_t<t and temp_t>1e-6){
            t =temp_t;
            hitobjectSphere = sphere_object[i];
            check = "Sphere";
        }
    }
    for(int i=0;i<plane_object.size();i++){
        double temp_t = plane_object[i].hit_plane(r);
        if(temp_t<t and temp_t>1e-6){
            t =temp_t;
            hitobjectPlane = plane_object[i];
            check = "Plane";
        }
    }
    for(int i=0;i<plane_light.size();i++){
        double temp_t = plane_light[i].hit_plane(r);
        if(temp_t<t and temp_t>1e-6){
            t =temp_t;
            hitobjectPlane = plane_light[i];
            check = "PlaneLight";
        }
    }
    for(int i=0;i<sphere_light.size();i++){
        double temp_t = sphere_light[i].hit_sphere(r);
        // if(temp_t!=INT_MAX) cout<<"i:"<<i<<" temp_t"<<temp_t<<"||";
        if(temp_t<t and temp_t>1e-6){
            t =temp_t;
            hitobjectSphere = sphere_light[i];
            check = "SphereLight";
        }
    }
    for(int i=0;i<point_light.size();i++){
        double temp_t = point_light[i].hit_PointLight(r);
        // if(temp_t!=INT_MAX) cout<<"i:"<<i<<" temp_t"<<temp_t<<"||";
        if(temp_t<t and temp_t>1e-6){
            t =temp_t;
            hitLightPoint = point_light[i];
            check = "PointLight";
        }
    }
    
    


    color rColor = color(0,0,0);
    if(t>1e-6 and check!="None") {
        if(check=="Sphere"){
            point3 hitPt = r.at(t);
            // Calculate UV coordinates for the hit point on the sphere
            vec3 normal = unit_vector(hitPt - hitobjectSphere.getCenter());
            double u = 0.5 + atan2(normal.z(), normal.x()) / (2 * M_PI);
            double v = 0.5 - asin(normal.y()) / M_PI;

            color texture_color = get_texture_color(u, v);
            rColor = texture_color;  // Apply texture color

            if(recDepth<maxDepth){
                if(hitobjectSphere.getKd()==0 and hitobjectSphere.getKa()==0){
                    vec3 Normal = unit_vector(hitPt - hitobjectSphere.getCenter());
                    ray reflected = ray(hitPt+1e-6*(r.direction()-2*dot(Normal,r.direction())*Normal),r.direction()-2*dot(Normal,r.direction())*Normal);
                    recDepth+=1;
                    rColor+=hitobjectSphere.getKr()*ray_color(reflected,point_light,sphere_light,plane_light,plane_object,sphere_object,ambient_Color,recDepth,maxDepth);
                }
                if(hitobjectSphere.getRefIndex()>0){
                    vec3 Normal = unit_vector(hitPt - hitobjectSphere.getCenter());
                    double relativeRef = hitobjectSphere.getRefIndex();
                    double Ro = pow((relativeRef-1)/(relativeRef+1),2);
                    double refSchlick = Ro + (1-Ro)*(1-dot(Normal,(-1)*r.direction()));
                    if(1-(pow(refSchlick,2)*(1-pow((dot(Normal,(-1)*r.direction())),2)))>0){
                        vec3 Direction = (refSchlick*(dot(Normal,(-1)*r.direction()))-sqrt(1-(pow(refSchlick,2)*(1-pow((dot(Normal,(-1)*r.direction())),2)))))*Normal - refSchlick*(-1)*r.direction();
                        hitobjectSphere.setRefIndex(1/relativeRef);
                        ray reflected = ray(hitPt+1e-6*Direction,Direction);
                        recDepth+=1;
                        rColor+=hitobjectSphere.getKt()*ray_color(reflected,point_light,sphere_light,plane_light,plane_object,sphere_object,ambient_Color,recDepth,maxDepth);
                    }
                }
            }
            for(int i=0;i<point_light.size();i++){
                vec3 lightDir = unit_vector(point_light[i].getCenter() - hitPt);
                double t = ((point_light[i].getCenter() - hitPt).x()/lightDir.x())+ ((point_light[i].getCenter() - hitPt).y()/lightDir.y())+  ((point_light[i].getCenter() - hitPt).z()/lightDir.z());
                // cout << "t:"<<t<<endl;
                ray shadow_ray = ray(hitPt + 1e-6*lightDir,lightDir);
                if(shadow_check(shadow_ray,point_light,sphere_light,plane_light,plane_object,sphere_object,t)){
                    vec3 Normal = unit_vector(hitPt - hitobjectSphere.getCenter());
                    vec3 half = unit_vector(r.direction()*(-1)+lightDir);
                    rColor+= point_light[i].getLightColor() * hitobjectSphere.getKd() * max(0.0,dot(Normal, lightDir));
                    rColor+= point_light[i].getLightColor() * hitobjectSphere.getKs() * pow(max(0.0,dot(Normal, half)),hitobjectSphere.getphongConst());
                }
            }
            for(int i=0;i<sphere_light.size();i++){
                vec3 lightDir = unit_vector(sphere_light[i].getCenter() - hitPt);
                double t = ((sphere_light[i].getCenter() - hitPt).x()/lightDir.x())+ ((sphere_light[i].getCenter() - hitPt).y()/lightDir.y())+  ((sphere_light[i].getCenter() - hitPt).z()/lightDir.z());
                // cout << "t:"<<t<<endl;
                ray shadow_ray = ray(hitPt + 1e-6*lightDir,lightDir);
                if(shadow_check(shadow_ray,point_light,sphere_light,plane_light,plane_object,sphere_object,t)){
                    vec3 Normal = unit_vector(hitPt - hitobjectSphere.getCenter());
                    vec3 half = unit_vector(r.direction()+lightDir);
                    rColor+= sphere_light[i].getObjectColor() * hitobjectSphere.getKd() * max(0.0,dot(Normal, lightDir));
                    rColor+= sphere_light[i].getObjectColor() * hitobjectSphere.getKs() * pow(max(0.0,dot(Normal, half)),hitobjectSphere.getphongConst());
                }
            }
            for(int i=0;i<plane_light.size();i++){
                vec3 lightDir = unit_vector(plane_light[i].getCenter() - hitPt);
                double t = ((plane_light[i].getCenter() - hitPt).x()/lightDir.x())+ ((plane_light[i].getCenter() - hitPt).y()/lightDir.y())+  ((plane_light[i].getCenter() - hitPt).z()/lightDir.z());
                // cout << "t:"<<t<<endl;
                ray shadow_ray = ray(hitPt + 1e-6*lightDir,lightDir);
                if(shadow_check(shadow_ray ,point_light,sphere_light,plane_light,plane_object,sphere_object,t)){
                    vec3 Normal = unit_vector(hitPt - hitobjectSphere.getCenter());
                    vec3 half = unit_vector(r.direction()+lightDir);
                    rColor+= plane_light[i].getObjectColor() * hitobjectSphere.getKd() * max(0.0,dot(Normal, lightDir));
                    rColor+= plane_light[i].getObjectColor() * hitobjectSphere.getKs() * pow(max(0.0,dot(Normal, half)),hitobjectSphere.getphongConst());
                }
            }
            rColor+=hitobjectSphere.getKa()*ambient_Color;
            rColor*=hitobjectSphere.getObjectColor();

        }else if(check=="Plane"){
            // cout <<"t:"<< t<<endl;
            point3 hitPt = r.at(t);
            if(recDepth<maxDepth){
                if(hitobjectPlane.getKd()==0 and hitobjectPlane.getKa()==0){
                    vec3 Normal = unit_vector(hitobjectPlane.getNormal());
                    ray reflected = ray(hitPt+1e-6*Normal,r.direction()-2*dot(Normal,r.direction())*Normal);
                    recDepth+=1;
                    rColor+=hitobjectPlane.getKr()*ray_color(reflected,point_light,sphere_light,plane_light,plane_object,sphere_object,ambient_Color,recDepth,maxDepth);
                }
                if(hitobjectPlane.getRefIndex()>0){
                    vec3 Normal = unit_vector(hitobjectPlane.getNormal());
                    double relativeRef = hitobjectPlane.getRefIndex();
                    double Ro = pow((relativeRef-1)/(relativeRef+1),2);
                    double refSchlick = Ro + (1-Ro)*(1-dot(Normal,(-1)*r.direction()));
                    if(1-(pow(refSchlick,2)*(1-pow((dot(Normal,(-1)*r.direction())),2)))>0){
                        vec3 Direction = (refSchlick*(dot(Normal,(-1)*r.direction()))-sqrt(1-(pow(refSchlick,2)*(1-pow((dot(Normal,(-1)*r.direction())),2)))))*Normal - refSchlick*(-1)*r.direction();
                        hitobjectPlane.setRefIndex(1/relativeRef);
                        ray reflected = ray(hitPt+1e-6*Direction,Direction);
                        recDepth+=1;
                        rColor+=hitobjectPlane.getKt()*ray_color(reflected,point_light,sphere_light,plane_light,plane_object,sphere_object,ambient_Color,recDepth,maxDepth);
                    }
                }
            }
            for(int i=0;i<point_light.size();i++){
                vec3 lightDir = unit_vector(point_light[i].getCenter() - hitPt);
                double t = ((point_light[i].getCenter() - hitPt).x()/lightDir.x())+ ((point_light[i].getCenter() - hitPt).y()/lightDir.y())+  ((point_light[i].getCenter() - hitPt).z()/lightDir.z());
                // cout << "t:"<<t<<endl;
                ray shadow_ray = ray(hitPt+ 1e-6*lightDir,lightDir);
                if(shadow_check(shadow_ray,point_light,sphere_light,plane_light,plane_object,sphere_object,t)){
                    vec3 Normal = unit_vector(hitobjectPlane.getNormal());
                    vec3 half = unit_vector(r.direction()*(-1)+lightDir);
                    rColor+= point_light[i].getLightColor() * hitobjectPlane.getKd() * max(0.0,dot(Normal, lightDir));
                    rColor+= point_light[i].getLightColor() * hitobjectPlane.getKs() * pow(max(0.0,dot(Normal, half)),hitobjectPlane.getphongConst());
                }
            }
            for(int i=0;i<sphere_light.size();i++){
                vec3 lightDir = unit_vector(sphere_light[i].getCenter() - hitPt);
                double t = ((sphere_light[i].getCenter() - hitPt).x()/lightDir.x())+ ((sphere_light[i].getCenter() - hitPt).y()/lightDir.y())+  ((sphere_light[i].getCenter() - hitPt).z()/lightDir.z());
                // cout << "t:"<<t<<endl;
                ray shadow_ray = ray(hitPt+ 1e-6*lightDir,lightDir);
                if(shadow_check(shadow_ray,point_light,sphere_light,plane_light,plane_object,sphere_object,t)){
                    vec3 Normal = unit_vector(hitobjectPlane.getNormal());
                    vec3 half = unit_vector(r.direction()+lightDir);
                    rColor+= sphere_light[i].getObjectColor() * hitobjectPlane.getKd() * max(0.0,dot(Normal, lightDir));
                    rColor+= sphere_light[i].getObjectColor() * hitobjectPlane.getKs() * pow(max(0.0,dot(Normal, half)),hitobjectPlane.getphongConst());
                }
            }
            for(int i=0;i<plane_light.size();i++){
                vec3 lightDir = unit_vector(plane_light[i].getCenter() - hitPt);
                double t = ((plane_light[i].getCenter() - hitPt).x()/lightDir.x())+ ((plane_light[i].getCenter() - hitPt).y()/lightDir.y())+  ((plane_light[i].getCenter() - hitPt).z()/lightDir.z());
                // cout << "t:"<<t<<endl;
                ray shadow_ray = ray(hitPt+ 1e-6*lightDir,lightDir);
                if(shadow_check(shadow_ray,point_light,sphere_light,plane_light,plane_object,sphere_object,t)){
                    vec3 Normal = unit_vector(hitobjectPlane.getNormal());
                    vec3 half = unit_vector(r.direction()+lightDir);
                    rColor+= plane_light[i].getObjectColor() * hitobjectPlane.getKd() * max(0.0,dot(Normal, lightDir));
                    rColor+= plane_light[i].getObjectColor() * hitobjectPlane.getKs() * pow(max(0.0,dot(Normal, half)),hitobjectSphere.getphongConst());
                }
            }
            rColor+=hitobjectPlane.getKa()*ambient_Color;
            rColor*=hitobjectPlane.getObjectColor();
        }else if(check=="PlaneLight"){
            return hitobjectPlane.getObjectColor();
        }else if(check=="SphereLight"){
            return hitobjectSphere.getObjectColor();
        }else if(check=="PointLight"){
            return hitLightPoint.getLightColor();
        }

    }else{
        return color(0.882,0.962,1);
    }


    return rColor;
}

int main() {
    vector<PointLight> point_light;
    vector<Sphere> sphere_light;
    vector<Plane> plane_light;
    vector<Plane> plane_object;
    vector<Sphere> sphere_object;
    // load_texture();
    vector<vector<color>> texture_image; 
    // Image

    auto aspect_ratio = 16.0 / 9.0;
    int image_width = 4000;

    // Calculate the image height, and ensure that it's at least 1.
    int image_height = int(image_width / aspect_ratio);
    image_height = (image_height < 1) ? 1 : image_height;
    
    // Camera
    float focal_length = 5;
    float viewport_height = 2.0;
    float viewport_width = viewport_height * (double(image_width) / image_height);
    point3 camera_center = point3(0,0,10);

    // Compute Vu and Vv
    vec3 viewport_u = vec3(viewport_width, 0, 0);
    vec3 viewport_v = vec3(0, -viewport_height, 0);
    
    // Compute du and dv
    vec3 pixel_delta_u = viewport_u / image_width;
    vec3 pixel_delta_v = viewport_v / image_height;

    // Calculate the location of upper left pixel
    point3 viewport_upper_left = camera_center - vec3(0,0,focal_length) - viewport_u / 2 - viewport_v /2;
    point3 pixel00_loc = viewport_upper_left + 0.5 * (pixel_delta_u + pixel_delta_v);

    //Adding point light
    point3 clight = point3(5,10,10);;
    color colorLight = color(1,1,1);
    PointLight light = PointLight(clight,colorLight);
    // point_light.push_back(light);

    //Adding sphere light
    point3 center_sphereLight = point3(0,0,-1);
    double radiusLight = 1;
    color colorSphereLight = color(1,1,0);
    double kd_sphereLight = 0.0;
    double ks_sphereLight = 0;
    double ka_sphereLight = 0.0;
    double kr_sphereLight = 0.0;
    double kt_sphereLight = 0.0;
    double refIndex_sphereLight = 0.0;
    double phongConst_sphereLight = 0;
    Sphere sphereLight = Sphere(center_sphereLight,radiusLight,kd_sphereLight,ks_sphereLight,ka_sphereLight,kr_sphereLight,kt_sphereLight,refIndex_sphereLight,phongConst_sphereLight,colorSphereLight);
    sphere_light.push_back(sphereLight);

    //Adding plane light
    point3 center_planeLight = point3(0,2,0);
    vec3 Normal_planeLight = vec3(0,-1,0);
    vec3 XminLight = vec3(-1,1,-1);
    vec3 XmaxLight = vec3(1,3,1);
    double kd_planeLight = 0.0;
    double ks_planeLight = 0.0;
    double ka_planeLight = 0.0;
    double kr_planeLight = 0.0;
    double kt_planeLight = 0.0;
    double refIndex_planeLight = 0.0;
    double phongConst_planeLight = 0.0;
    color color_planeLight = color(1,1,1);
    Plane planeLight = Plane(center_planeLight,Normal_planeLight,XminLight,XmaxLight,kd_planeLight,ks_planeLight,ka_planeLight,kr_planeLight,kt_planeLight,refIndex_planeLight,phongConst_planeLight,color_planeLight);
    // plane_light.push_back(planeLight);

    point3 clight2 = point3(10,5,0);;
    color colorLight2 = color(1,1,1);
    PointLight light2 = PointLight(clight2,colorLight2);
    point_light.push_back(light2);

    //Adding object Sphere
    point3 center_sphere = point3(-2,1,-2);
    double radius = 1;
    color colorSphere = color(0.9254,0.64,0.729);
    double kd_sphere = 0.0;
    double ks_sphere = 1;
    double ka_sphere = 0.0;
    double kr_sphere = 1.0;
    double kt_sphere = 1.0;
    double refIndex_sphere = 0.0;
    double phongConst_sphere = 500;
    Sphere sphere = Sphere(center_sphere,radius,kd_sphere,ks_sphere,ka_sphere,kr_sphere,kt_sphere,refIndex_sphere,phongConst_sphere,colorSphere);
    sphere_object.push_back(sphere);
    point3 center_sphere2 = point3(2,0,-2.5);
    double radius2 = 1;
    color colorSphere2 = color(0.596,0.960,0.709);
    double kd_sphere2 = 0.0;
    double ks_sphere2 = 1.0;
    double ka_sphere2 = 0.0;
    double kr_sphere2 = 1.0;
    double kt_sphere2 = 1.0;
    double refIndex_sphere2 = 0.0;
    double phongConst_sphere2 = 500;
    Sphere sphere2 = Sphere(center_sphere2,radius2,kd_sphere2,ks_sphere2,ka_sphere2,kr_sphere2,kt_sphere2,refIndex_sphere2,phongConst_sphere2,colorSphere2);
    sphere_object.push_back(sphere2);
    point3 center_sphere3 = point3(1,0,0);
    double radius3 = 0.5;
    color colorSphere3 = color(0.517,0.549,0.968);
    double kd_sphere3 = 0.0;
    double ks_sphere3 = 1;
    double ka_sphere3 = 0.0;
    double kr_sphere3 = 0.5;
    double kt_sphere3 = 1;
    double refIndex_sphere3 = 0.25;
    double phongConst_sphere3 = 500;
    Sphere sphere3 = Sphere(center_sphere3,radius3,kd_sphere3,ks_sphere3,ka_sphere3,kr_sphere3,kt_sphere3,refIndex_sphere3,phongConst_sphere3,colorSphere3);
    sphere_object.push_back(sphere3);

    //Adding object Plane;
    point3 center_plane = point3(-10,-1,-1);
    vec3 Normal_plane = vec3(0,1,0);
    vec3 Xmin = vec3(-100,-5,-200);
    vec3 Xmax = vec3(100,-0.5,100);
    double kd_plane = 0.8;
    double ks_plane = 0.5;
    double ka_plane = 1;
    double kr_plane = 0.0;
    double kt_plane = 0.0;
    double refIndex_plane = 0.0;
    double phongConst_plane = 0.9;
    color color_plane = color(0.4,0.4,0.4);
    Plane plane = Plane(center_plane,Normal_plane,Xmin,Xmax,kd_plane,ks_plane,ka_plane,kr_plane,kt_plane,refIndex_plane,phongConst_plane,color_plane);
    plane_object.push_back(plane);

    point3 center_plane2 = point3(0,0,0);
    vec3 Normal_plane2 = vec3(1,0,1);
    vec3 Xmin2 = vec3(-1,-1,-2);
    vec3 Xmax2 = vec3(1,1,2);
    double kd_plane2 = 0.0;
    double ks_plane2 = 0.5;
    double ka_plane2 = 0.0;
    double kr_plane2 = 1.0;
    double kt_plane2 = 0.5;
    double refIndex_plane2 = 0.9;
    double phongConst_plane2 = 0.9;
    color color_plane2 = color(0.5,0.5,0.5);
    Plane plane2 = Plane(center_plane2,Normal_plane2,Xmin2,Xmax2,kd_plane2,ks_plane2,ka_plane2,kr_plane2,kt_plane2,refIndex_plane2,phongConst_plane2,color_plane2);
    // plane_object.push_back(plane2);

    //Ambient Color or Light;
    color ambient_Color = color(0.5,0.5,0.5);
    // Render
    int maxDepth = 39;

    std::cout << "P3\n" << image_width << ' ' << image_height << "\n255\n";
    // i controls the columns and j controls the rows
    for(int j=0; j<image_height; j++) {
        std::clog << "\rScanlines remaining: " << (image_height - j) << ' ' << std::flush;
        for(int i=0; i< image_width; i++) {
                 vec3 pixel_sample = pixel00_loc
                                + ((i) * pixel_delta_u)
                                + ((j) * pixel_delta_v);
                vec3 ray_direction = unit_vector(pixel_sample - camera_center); 
                ray r(camera_center, ray_direction);

                auto pixel_color = ray_color(r,point_light,sphere_light,plane_light,plane_object,sphere_object,ambient_Color,1,maxDepth);
            write_color(std::cout, pixel_color);
        }
    }
    std::clog << "\rDone.                                \n";

}