//Sean M. Law

#include "Vector.hpp"

#include <cmath>

Vector::Vector(){
  xcoor=0.0;
  ycoor=0.0;
  zcoor=0.0;
}

Vector::Vector(double xcoorin, double ycoorin, double zcoorin){
  xcoor=xcoorin;
  ycoor=ycoorin;
  zcoor=zcoorin;
}

Vector::Vector(const Vector& vec){
  xcoor=vec.xcoor;
  ycoor=vec.ycoor;
  zcoor=vec.zcoor;
}

Vector& Vector::operator= (const Vector& vec){
  xcoor=vec.xcoor;
  ycoor=vec.ycoor;
  zcoor=vec.zcoor;
  return(*this);
}

Vector& Vector::operator= (const double val){
  xcoor=val;
  ycoor=val;
  zcoor=val;
  return(*this);
}

//Addition
Vector Vector::operator+ (const Vector& vec) const{
	return Vector(xcoor+vec.xcoor,ycoor+vec.ycoor,zcoor+vec.zcoor);
}

Vector& Vector::operator+= (const Vector& vec){
  xcoor+=vec.xcoor;
  ycoor+=vec.ycoor;
  zcoor+=vec.zcoor;
  return(*this);
}

Vector Vector::operator+ (const double val) const{
	return Vector(xcoor+val,ycoor+val,zcoor+val);
}

Vector& Vector::operator+= (const double val){
  xcoor+=val;
  ycoor+=val;
  zcoor+=val;
  return(*this);
}

//Subtraction
Vector Vector::operator- (const Vector& vec) const{
	return Vector(xcoor-vec.xcoor,ycoor-vec.ycoor,zcoor-vec.zcoor);
}

Vector& Vector::operator-= (const Vector& vec){
  xcoor-=vec.xcoor;
  ycoor-=vec.ycoor;
  zcoor-=vec.zcoor;
  return(*this);
}

Vector Vector::operator- (const double val) const{
	return Vector(xcoor-val,ycoor-val,zcoor-val);
}

Vector& Vector::operator-= (const double val){
  xcoor-=val;
  ycoor-=val;
  zcoor-=val;
  return(*this);
}

//Multiplication
Vector Vector::operator* (const Vector& vec) const{
	return Vector(xcoor*vec.xcoor,ycoor*vec.ycoor,zcoor*vec.zcoor);
}

Vector& Vector::operator*= (const Vector& vec){
  xcoor*=vec.xcoor;
  ycoor*=vec.ycoor;
  zcoor*=vec.zcoor;
  return(*this);
}

Vector Vector::operator* (const double val) const{
	return Vector(xcoor*val,ycoor*val,zcoor*val);
}

Vector& Vector::operator*= (const double val){
  xcoor*=val;
  ycoor*=val;
  zcoor*=val;
  return(*this);
}

//Division
Vector Vector::operator/ (const Vector& vec) const{
	return Vector(xcoor/vec.xcoor,ycoor/vec.ycoor,zcoor/vec.zcoor);
}

Vector& Vector::operator/= (const Vector& vec){
  xcoor/=vec.xcoor;
  ycoor/=vec.ycoor;
  zcoor/=vec.zcoor;
  return(*this);
}

Vector Vector::operator/ (const double val) const{
	return Vector(xcoor/val,ycoor/val,zcoor/val);
}

Vector& Vector::operator/= (const double val){
  xcoor/=val;
  ycoor/=val;
  zcoor/=val;
  return(*this);
}


Vector Vector::operator- () const {
  return Vector(-xcoor,-ycoor,-zcoor);
}

double Vector::dot (const Vector& vec) const { //Dot Product
	return xcoor*vec.xcoor+ycoor*vec.ycoor+zcoor*vec.zcoor;
}

Vector Vector::cross (const Vector& vec) const { //Cross Product
  return Vector(ycoor*vec.zcoor - zcoor*vec.ycoor,
                zcoor*vec.xcoor - xcoor*vec.zcoor,
                xcoor*vec.ycoor - ycoor*vec.xcoor);
}

double Vector::norm () const { //Normal
  return sqrt(xcoor*xcoor+ycoor*ycoor+zcoor*zcoor);
}

