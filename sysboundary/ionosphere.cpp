/*
 * This file is part of Vlasiator.
 * 
 * Copyright 2010, 2011, 2012 Finnish Meteorological Institute
 * 
 * Vlasiator is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License version 3
 * as published by the Free Software Foundation.
 * 
 * Vlasiator is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program. If not, see <http://www.gnu.org/licenses/>.
 */

#include <cstdlib>
#include <mpi.h>
#include <iostream>
#include <limits>

#include "sysboundarycondition.h"
#include "ionosphere.h"
#include "../parameters.h"

using namespace std;

namespace SBC {
   Ionosphere::Ionosphere(): SysBoundaryCondition() { }
   Ionosphere::~Ionosphere() { }
   
   bool Ionosphere::initSysBoundary() {
      bool result = setCenter();
      result = result & setRadius();      
      return result;
   }
   
   int Ionosphere::assignSysBoundary(creal* cellParams) {
      creal dx = cellParams[CellParams::DX];
      creal dy = cellParams[CellParams::DY];
      creal dz = cellParams[CellParams::DZ];
      creal x = cellParams[CellParams::XCRD] + 0.5*dx;
      creal y = cellParams[CellParams::YCRD] + 0.5*dy;
      creal z = cellParams[CellParams::ZCRD] + 0.5*dz;
      
      creal r = sqrt((x-center[0])*(x-center[0]) + (y-center[1])*(y-center[1]) + (z-center[2])*(z-center[2]));
      Real rx, ry, rz;
      
      int typeToAssign = sysboundarytype::NOT_SYSBOUNDARY;
      
      if(r < radius) {
         // Determine the sign of the quadrant in each direction and calculate the radius for a cell located one cell length further that direction.
         // If at least one of these further cells is outside the ionospheric radius then the current one is an ionosphere cell.
         // If not, then the cell is inside the ionosphere radius and should not need to be computed.
         creal quadrant[3] = {x/abs(x), y/abs(y), z/abs(z)};
         rx = sqrt((x+quadrant[0]*dx-center[0])*(x+quadrant[0]*dx-center[0]) + (y-center[1])*(y-center[1]) + (z-center[2])*(z-center[2]));
         ry = sqrt((x-center[0])*(x-center[0]) + (y+quadrant[1]*dy-center[1])*(y+quadrant[1]*dy-center[1]) + (z-center[2])*(z-center[2]));
         rz = sqrt((x-center[0])*(x-center[0]) + (y-center[1])*(y-center[1]) + (z+quadrant[2]*dz-center[2])*(z+quadrant[2]*dz-center[2]));
         if(rx>radius || ry>radius || rz>radius) {
            typeToAssign = getIndex();
         } else {
            typeToAssign = sysboundarytype::DO_NOT_COMPUTE;
         }
      }
      
      return typeToAssign;
   }
   
   bool Ionosphere::setCenter() {
      for(uint i=0; i<3; i++) center[i] = Parameters::ionoCenter[1];
      return true;
   }
   
   bool Ionosphere::setRadius() {
      radius = Parameters::ionoRadius;
      return true;
   }
   
   std::string Ionosphere::getName() const {return "Ionosphere";}
   
   int Ionosphere::getIndex() const {
      return sysboundarytype::IONOSPHERE;
   }
}
