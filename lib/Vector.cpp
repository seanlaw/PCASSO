//Sean M. Law
//Aaron T. Frank

/*
This file is part of PCASSO.

PCASSO is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

PCASSO is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with PCASSO.  If not, see <http://www.gnu.org/licenses/>.

This file is part of MoleTools.

MoleTools is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

MoleTools is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with MoleTools.  If not, see <http://www.gnu.org/licenses/>.
*/

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

