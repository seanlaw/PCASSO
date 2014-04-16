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

#ifndef VECTOR_H
#define VECTOR_H

//#include <cmath>
//#include <iostream>

class Vector {
  private:
    double xcoor;
    double ycoor;
    double zcoor;

  public:
    Vector();
    Vector(double xcoorin, double ycoorin, double zcoorin); //Constructor
    Vector(const Vector& vec); //Overload Constructor 

    Vector& operator= (const Vector& vec);
    Vector& operator= (const double val);
    //Addition
    Vector operator+ (const Vector& vec) const;
    Vector& operator+= (const Vector& vec);
    Vector operator+ (const double val) const;
    Vector& operator+= (const double val);
    //Subtraction
    Vector operator- (const Vector& vec) const;
    Vector& operator-= (const Vector& vec);
    Vector operator- (const double val) const;
    Vector& operator-= (const double val);
    //Multiplication
    Vector operator* (const Vector& vec) const;
    Vector& operator*= (const Vector& vec);
    Vector operator* (const double val) const;
    Vector& operator*= (const double val);
    //Division
    Vector operator/ (const Vector& vec) const;
    Vector& operator/= (const Vector& vec);
    Vector operator/ (const double val) const;
    Vector& operator/= (const double val);

    double& x(){return xcoor;};
    double& y(){return ycoor;};
    double& z(){return zcoor;};

    Vector operator- () const;
    double dot (const Vector& vec) const; //Dot Product
    Vector cross (const Vector& vec) const; //Cross Product
    double norm () const;
};

#endif
