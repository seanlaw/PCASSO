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

#include "Prmtop.hpp"
#include "Misc.hpp"

#include <fstream>

Prmtop::Prmtop (){
  mass.clear();
  charge.clear();
}

void Prmtop::readTopology(const std::string& topin){
  std::ifstream topFile;
  std::istream* topinp;
  std::string line;
  std::vector<std::string> s;
  double m; //mass
  double c; //charge
  std::map<std::string, double> massRef;
  std::string resname;
  std::string word;

  resname.clear();

  if (topin.length() > 0){
    topFile.open(topin.c_str(), std::ios::in);
    topinp=&topFile;

    while (topinp->good() && !(topinp->eof())){
      getline(*topinp, line);
      word=Misc::trim(line).substr(0,4);
      Misc::toupper(word);
      if (word.compare(0,4,"MASS") == 0){
        Misc::splitStr(line, " \t", s, false);
        if (s.size() >= 4){
          std::stringstream(s.at(3)) >> m;
          massRef.insert(std::make_pair(s.at(2), m));
        }
      }
      else if (word.compare(0,4,"RESI") == 0 || Misc::trim(line).compare(0,4,"PRES") == 0){
        Misc::splitStr(line, " \t", s, false);
        if (s.size() >= 3){
          resname=s.at(1);
        }
      }
      else if (word.compare(0,4,"ATOM") == 0){
        Misc::splitStr(line, " \t", s, false);
        if (s.size() >= 3){
          mass.insert(std::make_pair(std::make_pair(resname, s.at(1)), massRef[s.at(2)]));
          std::stringstream(s.at(3)) >> c;
          charge.insert(std::make_pair(std::make_pair(resname, s.at(1)), c));
        }
      }
      else{
        //Do nothing
      }
    }
  }

}


void Prmtop::readParameter(const std::string& prmin){

}


double Prmtop::getMass(const std::string& resnamein, const std::string& atmnamein){
  if (mass.find(std::make_pair(resnamein, atmnamein)) != mass.end()){
    return  mass[std::make_pair(resnamein, atmnamein)];
  }
  else{
    std::cerr << "Warning: Could not a find mass for atom " << resnamein;
    std::cerr << " " << atmnamein << " and was set to 1.0" << std::endl;
    return 1.0;
  }
}

double Prmtop::getCharge(const std::string& resnamein, const std::string& atmnamein){
  if (charge.find(std::make_pair(resnamein, atmnamein)) != charge.end()){
    return charge[std::make_pair(resnamein, atmnamein)];
  }
  else{
    std::cerr << "Warning: Could not find a charge for atom " << resnamein;
    std::cerr << " " << atmnamein << " and was set to 0.0." << std::endl;
    return 0.0;
  }
}

