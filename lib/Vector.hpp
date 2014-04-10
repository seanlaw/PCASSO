//Sean M. Law

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
