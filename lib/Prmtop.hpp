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

#ifndef PRMTOP_H
#define PRMTOP_H

#include <string>
#include <map>

class Prmtop {
  private:
    //Maps are referenced by pair(resname, atmname)
    std::map<std::pair<std::string, std::string>, double> mass;
    std::map<std::pair<std::string, std::string>, double> charge;

  public:
    Prmtop(); //Constructor
    void readTopology(const std::string& topin);
    void readParameter(const std::string& prmin);
    double getMass(const std::string& resnamein, const std::string& atmnamein);
    double getCharge(const std::string& resnamein, const std::string& atmnamein);
};

#endif
