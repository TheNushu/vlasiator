/*
 * This file is part of Vlasiator.
 * Copyright 2010-2016 Finnish Meteorological Institute
 *
 * For details of usage, see the COPYING file and read the "Rules of the Road"
 * at http://www.physics.helsinki.fi/vlasiator/
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with this program; if not, write to the Free Software Foundation, Inc.,
 * 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */

#include <cmath>
#include <cstdlib>
#include <iostream>

#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"
#include "../../common.h"
#include "../../object_wrapper.h"
#include "../../readparameters.h"

#include "MultiCircularAlfven.h"

using namespace spatial_cell;

namespace projects {
MultiCircularAlfven::MultiCircularAlfven() : TriAxisSearch() {}
MultiCircularAlfven::~MultiCircularAlfven() {}

bool MultiCircularAlfven::initialize(void) {
   bool success = Project::initialize();

   creal m = physicalconstants::MASS_PROTON;
   creal e = physicalconstants::CHARGE;
   creal kB = physicalconstants::K_B;
   creal gamma = 5.0 / 3.0;
   creal mu0 = physicalconstants::MU_0;

   n = rho0 / m; // number density
   p0 = n * kB * T;
   // n should be 1e6 with rho0 = 1.6726219e-21
   
   cosalpha = cos(alpha);
   sinalpha = sin(alpha);
   // wave vector

   kwave = 2 * M_PI / lambda;
   // Alfven speed
   creal VA = B / sqrt(mu0 * rho0); 
   B1 = v1 * B / VA;

   if (verbose) {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      if (myRank == MASTER_RANK) {
         if (P::ohmHallTerm == 0) {
            std::cerr << "The Hall term is needed for accurate Alfven waves!";
         }
         std::cout << "wavelength = " << lambda << " m\n";
         std::cout << "Alfven speed = " << VA << " m/s\n";
         std::cout << "velocity perturbation amplitude = " << v1 << " m/s\n";
         std::cout << "magnetic perturbation amplitude = " << B1 << " T\n";
         std::cout << "angle = " << alpha * 180 / M_PI << " [degree]\n";
         std::cout << "-------------------------------------------------\n";
         std::cout << "rho0: " << rho0 << " kg/m^3\n";
         std::cout << "pressure (p0): " << p0 << " Pa\n";
         std::cout << "magnetic field (B): " << B << " T\n";
         std::cout << "Temperature (T): " << T << " K\n";
         std::cout << "-------------------------------------------------\n";
         creal Lx = P::xmax - P::xmin;
         //How to use this with SI:std::cout << "wavelength / di = " << Lx / lRef << "\n";
         // Useful for setting the Vspace ~ 10*Vth width
         Real Vth = sqrt(gamma * p0 / rho0);
         std::cout << "Vth: " << Vth << " [m/s]\n";
         const vmesh::MeshParameters& mp = getObjectWrapper().velocityMeshes[0];
         Real dvx = (mp.meshLimits[1] - mp.meshLimits[0]) / (mp.gridLength[0] * mp.blockLength[0]);
         std::cout << "Vmax:" << mp.meshLimits[1] << "\n";
         std::cout << "Vmin:" << mp.meshLimits[0] << "\n";
         std::cout << "nvcells:" << mp.gridLength[0] * mp.blockLength[0] << "\n";
         std::cout << "dvx = " << dvx << " [m/s]\n";
         std::cout << "Vth/dvx = " << Vth / dvx << "\n";
         std::cout << "Vmax/Vth = " << mp.meshLimits[1] / Vth << "\n";
         if (10 * dvx > Vth) {
            std::cout << "vmesh maybe too coarse!\n";
         }
         if (5 * Vth > mp.meshLimits[1]) {
            std::cout << "vmesh extent maybe too small!\n";
         }
         std::cout << "------------------------------------------------\n";
      }
   }

   return success;
}


void MultiCircularAlfven::addParameters() {
   typedef Readparameters RP;

   // Define all parameters with SI units directly
   RP::add("MultiCircularAlfven.v1", "Perturbation velocity amplitude (m/s)", 0.1);
   RP::add("MultiCircularAlfven.lambda", "Wavelength (m)", 32.0);
   RP::add("MultiCircularAlfven.rho0", "Density (kg/m^3)", 1.6726219e-21);
   RP::add("MultiCircularAlfven.p0", "Pressure (Pa)", 0.1);
   RP::add("MultiCircularAlfven.B", "Magnetic field strength (T)", 1e-8);
   RP::add("MultiCircularAlfven.alpha", "Wave vector angle (radians)", 0.4636476090008061);
   RP::add("MultiCircularAlfven.T", "Temperature (K)", 1e6);  // Example temperature
   RP::add("MultiCircularAlfven.verbose", "Verbose output", 1);
}

void MultiCircularAlfven::getParameters() {
   Project::getParameters();

   typedef Readparameters RP;
   RP::get("MultiCircularAlfven.v1", v1);
   RP::get("MultiCircularAlfven.lambda", lambda);
   RP::get("MultiCircularAlfven.alpha", alpha);
   RP::get("MultiCircularAlfven.rho0", rho0);
   RP::get("MultiCircularAlfven.p0", p0);
   RP::get("MultiCircularAlfven.B", B);
   RP::get("MultiCircularAlfven.T", T);
   RP::get("MultiCircularAlfven.verbose", verbose);
}

Real MultiCircularAlfven::getMaxwellian(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx,
                                   creal& dvy, creal& dvz, const uint popID) const {
   creal m = getObjectWrapper().particleSpecies[popID].mass;
   creal kB = physicalconstants::K_B;

   Real xpar = x * cosalpha + y * sinalpha;
   //phase1 initially 0, but can be changed if wanted to
   Real uperp1 = v1 * sin(kwave * xpar + phase1); 
   Real uz1 = v1 * cos(kwave * xpar + phase1);

   Real uperp2 = -v1 * sin(kwave * xpar + phase2);
   Real uz2 = -v1 * cos(kwave * xpar + phase2);

   Real ux = (-uperp1 * sinalpha) + (-uperp2 * sinalpha);
   Real uy = (uperp1 * cosalpha) + (uperp2 * cosalpha);
   Real uz = uz1 + uz2;
   
   creal coef = m / (2 * M_PI * kB * T);
   creal f = n * sqrt(coef) * coef * exp(-coef * M_PI * (sqr(vx - ux) + sqr(vy - uy) + sqr(vz - uz)));

   return f;
}


Real MultiCircularAlfven::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx,
                                           creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,
                                           const uint popID) const {
   Real f = 0.0;
   f = getMaxwellian(x + 0.5 * dx, y + 0.5 * dy, z + 0.5 * dz, vx + 0.5 * dvx, vy + 0.5 * dvy, vz + 0.5 * dvz, dvx, dvy,
                     dvz, popID);

   return f;
}

std::vector<std::array<Real, 3>> MultiCircularAlfven::getV0(creal x, creal y, creal z, const uint popID) const {
   std::vector<std::array<Real, 3>> V0;
   std::array<Real, 3> v = {{0.0, 0.0, 0.0}};
   V0.push_back(v);
   return V0;
}

void MultiCircularAlfven::calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) {}

void MultiCircularAlfven::setProjectBField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                                      FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                                      FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
   // Background field
   ConstantField bgField;
   bgField.initialize(B * cosalpha, B * sinalpha, 0.0);
   setBackgroundField(bgField, BgBGrid);

   if (!P::isRestart) {
      auto localSize = perBGrid.getLocalSize().data();

#pragma omp parallel for collapse(3)
      for (int i = 0; i < localSize[0]; ++i) {
         for (int j = 0; j < localSize[1]; ++j) {
            for (int k = 0; k < localSize[2]; ++k) {
               const std::array<Real, 3> x = perBGrid.getPhysicalCoords(i, j, k);
               std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(i, j, k);
               Real xpar = x[0] * cosalpha + x[1] * sinalpha;
               
               Real Bperp1 = B1 * sin(kwave * xpar + phase1);
               Real Bperp2 = B1 * sin(kwave * xpar + phase2);

               Real Bx1 = (-Bperp1 * sinalpha);
               Real Bx2 = (-Bperp2 * sinalpha);

               Real By1 = Bperp1 * cosalpha;
               Real By2 = Bperp2 * cosalpha;

               Real Bz1 = B1 * cos(kwave * xpar + phase1);
               Real Bz2 = B1 * cos(kwave * xpar + phase2);

               cell->at(fsgrids::bfield::PERBX) = Bx1 + Bx2;
               cell->at(fsgrids::bfield::PERBY) = By1 + By2;
               cell->at(fsgrids::bfield::PERBZ) = Bz1 + Bz2;
            }
         }
      }
   }
}


} // namespace projects
