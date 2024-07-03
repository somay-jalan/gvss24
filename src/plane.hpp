#ifndef PLANE_H
#define PLANE_H

#include "vec3.h"
#include "color.h"
#include <climits>


class Plane {
  public:
    Plane() {}

    Plane(const point3& center, const vec3& Normal, const vec3& Xmin, const vec3& Xmax, const double Kd, const double Ks, const double phongConst,  const color& objColor) 
    : center(center), Normal(Normal), Xmin(Xmin), Xmax(Xmax), Kd(Kd), Ks(Ks), phongConst(phongConst), objectColor(objColor) {}

    const point3& getCenter() const  { return center; }
    const vec3& getNormal() const { return Normal; }
    const double getKd() const { return Kd; }
    const double getKs() const { return Ks; }
    const double getphongConst() const { return phongConst; }
    const color& getObjectColor() const  { return objectColor; }
    
    double hit_plane(const ray& r) {
        double t;
        t = dot(Normal, center - r.origin()) / (dot(Normal, r.direction()));
        if(t>0){
            point3 hitPoint = r.at(t);
            if(hitPoint[0] >= Xmin[0] && hitPoint[1] >= Xmin[1] && hitPoint[2] >= Xmin[2] && hitPoint[0] <= Xmax[0] && hitPoint[1] <= Xmax[1] && hitPoint[2] <= Xmax[2])
              return t;
            else return INT_MAX;
        }
        else
            return INT_MAX;
    }
  private:
    point3 center;
    vec3 Normal;
    color objectColor;
    double Kd;
    double Ks;
    double phongConst;
    vec3 Xmin;
    vec3 Xmax;

};

#endif