#ifndef POINTLIGHT_H
#define POINTLIGHT_H

#include "vec3.h"
#include "color.h"
#include <climits>



class PointLight {
  public:
    PointLight() {}

    PointLight(const point3& center, const color& lightColor) : center(center), lightColor(lightColor) {}

    const point3& getCenter() const  { return center; }
    const color& getLightColor() const  { return lightColor; }
    
    double hit_PointLight(const ray& r) {
        vec3 oc = center - r.origin();
        double a = dot(oc, r.direction());
        if(a==-1 or a==1){
            return ((oc.x()/r.direction().x())+(oc.y()+r.direction().y())+(oc.z()+r.direction().z()));
        }
        else
            return INT_MAX;
    }
  private:
    point3 center;
    double radius;
    color lightColor;
    double Kd;
    double Ks;
    double Ka;

};

#endif