#include "../src/color.h"
#include "../src/vec3.h"
#include "../src/ray.h"
#include "../src/sphere.hpp"
#include "../src/plane.hpp"
#include "../src/pointLight.hpp"

#include<bits/stdc++.h>
#include <iostream>

using namespace std;

color ray_color(const ray& r,
    vector<PointLight> &point_light,
    vector<Sphere> &sphere_light,
    vector<Plane> &plane_light,
    vector<Plane> &plane_object,
    vector<Sphere> &sphere_object,
    color &ambient_Color,
    double ka) {
    
    Sphere hitobjectSphere;
    Plane hitobjectPlane;
    string check="None";
    double t = 1e8;;
    for(int i=0;i<sphere_object.size();i++){
        double temp_t = sphere_object[i].hit_sphere(r);
        // if(temp_t!=INT_MAX) cout<<"i:"<<i<<" temp_t"<<temp_t<<"||";
        if(temp_t<t){
            t =temp_t;
            hitobjectSphere = sphere_object[i];
            check = "Sphere";
        }
    }
    for(int i=0;i<plane_object.size();i++){
        double temp_t = plane_object[i].hit_plane(r);
        if(temp_t<t){
            t =temp_t;
            hitobjectPlane = plane_object[i];
            check = "Plane";
        }
    }
    // cout <<"t:"<< t<<endl;
    // cout <<"check:"<< check<<endl;


    color rColor = color(0,0,0);
    if(t >0 and t<1e8) {
        if(check=="Sphere"){
            point3 hitPt = r.at(t);
            for(int i=0;i<point_light.size();i++){
                vec3 lightDir = unit_vector(point_light[i].getCenter() - hitPt);
                vec3 Normal = unit_vector(hitPt - hitobjectSphere.getCenter());
                vec3 half = unit_vector(r.direction()+lightDir);
                rColor+= point_light[i].getLightColor() * hitobjectSphere.getKd() * max(0.0,dot(Normal, lightDir));
                rColor+= point_light[i].getLightColor() * hitobjectSphere.getKs() * pow(max(0.0,dot(Normal, half)),hitobjectSphere.getphongConst());
            }
            for(int i=0;i<sphere_light.size();i++){
                vec3 lightDir = unit_vector(sphere_light[i].getCenter() - hitPt);
                vec3 Normal = unit_vector(hitPt - hitobjectSphere.getCenter());
                vec3 half = unit_vector(r.direction()+lightDir);
                rColor+= sphere_light[i].getObjectColor() * hitobjectSphere.getKd() * max(0.0,dot(Normal, lightDir));
                rColor+= sphere_light[i].getObjectColor() * hitobjectSphere.getKs() * pow(max(0.0,dot(Normal, half)),hitobjectSphere.getphongConst());
            }
            for(int i=0;i<plane_light.size();i++){
                vec3 lightDir = unit_vector(plane_light[i].getNormal() - hitPt);
                vec3 Normal = unit_vector(hitPt - hitobjectSphere.getCenter());
                vec3 half = unit_vector(r.direction()+lightDir);
                rColor+= plane_light[i].getObjectColor() * hitobjectSphere.getKd() * max(0.0,dot(Normal, lightDir));
                rColor+= plane_light[i].getObjectColor() * hitobjectSphere.getKs() * pow(max(0.0,dot(Normal, half)),hitobjectSphere.getphongConst());
            }

        }else if(check=="Plane"){
            // cout <<"t:"<< t<<endl;
            point3 hitPt = r.at(t);
            for(int i=0;i<point_light.size();i++){
                vec3 lightDir = unit_vector(point_light[i].getCenter() - hitPt);
                vec3 Normal = unit_vector(hitobjectPlane.getNormal());
                vec3 half = unit_vector(r.direction()+lightDir);
                rColor+= point_light[i].getLightColor() * hitobjectPlane.getKd() * max(0.0,dot(Normal, lightDir));
                rColor+= point_light[i].getLightColor() * hitobjectPlane.getKs() * pow(max(0.0,dot(Normal, half)),hitobjectPlane.getphongConst());
            }
            for(int i=0;i<sphere_light.size();i++){
                vec3 lightDir = unit_vector(sphere_light[i].getCenter() - hitPt);
                vec3 Normal = unit_vector(hitobjectPlane.getNormal());
                vec3 half = unit_vector(r.direction()+lightDir);
                rColor+= sphere_light[i].getObjectColor() * hitobjectPlane.getKd() * max(0.0,dot(Normal, lightDir));
                rColor+= sphere_light[i].getObjectColor() * hitobjectPlane.getKs() * pow(max(0.0,dot(Normal, half)),hitobjectPlane.getphongConst());
            }
            for(int i=0;i<plane_light.size();i++){
                vec3 lightDir = unit_vector(plane_light[i].getNormal() - hitPt);
                vec3 Normal = unit_vector(hitobjectPlane.getNormal());
                vec3 half = unit_vector(r.direction()+lightDir);
                rColor+= plane_light[i].getObjectColor() * hitobjectPlane.getKd() * max(0.0,dot(Normal, lightDir));
                rColor+= plane_light[i].getObjectColor() * hitobjectPlane.getKs() * pow(max(0.0,dot(Normal, half)),hitobjectSphere.getphongConst());
            }
        }
        rColor+=ka*ambient_Color;
        rColor*=hitobjectSphere.getObjectColor();

    }else{
        return color(0.99,0.99,0.99);

    }


    return rColor;
}

int main() {
    vector<PointLight> point_light;
    vector<Sphere> sphere_light;
    vector<Plane> plane_light;
    vector<Plane> plane_object;
    vector<Sphere> sphere_object;

    // Image

    auto aspect_ratio = 16.0 / 9.0;
    int image_width = 800;

    // Calculate the image height, and ensure that it's at least 1.
    int image_height = int(image_width / aspect_ratio);
    image_height = (image_height < 1) ? 1 : image_height;
    
    // Camera
    float focal_length = 1.0;
    float viewport_height = 2.0;
    float viewport_width = viewport_height * (double(image_width) / image_height);
    point3 camera_center = point3(0,0,1);

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
    point3 clight = point3(10,10,10);
    color colorLight = color(1,1,1);
    PointLight light = PointLight(clight,colorLight);
    point_light.push_back(light);

    //Adding object Sphere
    point3 center_sphere = point3(0,0,-5);
    double radius = 1;
    color colorSphere = color(0.8,0.4,0.4);
    double kd_sphere = 0.8;
    double ks_sphere = 0.5;
    double phongConst_sphere = 0.9;
    Sphere sphere = Sphere(center_sphere,radius,kd_sphere,ks_sphere,phongConst_sphere,colorSphere);
    sphere_object.push_back(sphere);
    point3 center_sphere2 = point3(2,2,-3);
    double radius2 = 1;
    color colorSphere2 = color(0.4,0.4,0.4);
    double kd_sphere2 = 0.8;
    double ks_sphere2 = 0.5;
    double phongConst_sphere2 = 0.9;
    Sphere sphere2 = Sphere(center_sphere2,radius2,kd_sphere2,ks_sphere2,phongConst_sphere2,colorSphere2);
    sphere_object.push_back(sphere2);

    //Adding object Plane;
    point3 center_plane = point3(0,0,-1.5);
    vec3 Normal_plane = vec3(0,0,1);
    vec3 Xmin = vec3(-1,-1,-2);
    vec3 Xmax = vec3(1,1,-1);
    double kd_plane = 0.8;
    double ks_plane = 0.5;
    double phongConst_plane = 0.9;
    color color_plane = color(0.5,01,0.1);
    Plane plane = Plane(center_plane,Normal_plane,Xmin,Xmax,kd_plane,ks_plane,phongConst_plane,color_plane);
    plane_object.push_back(plane);

    //Ambient Color or Light;
    color ambient_Color = color(1,1,1);
    double ka = 0.1; 
    // Render

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

                auto pixel_color = ray_color(r,point_light,sphere_light,plane_light,plane_object,sphere_object,ambient_Color,ka);
            write_color(std::cout, pixel_color);
        }
    }
    std::clog << "\rDone.                                \n";

}