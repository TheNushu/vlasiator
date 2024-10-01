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

#include "CircularAlfven.h"

using namespace spatial_cell;

namespace projects {
CircularAlfven::CircularAlfven() : TriAxisSearch() {}
CircularAlfven::~CircularAlfven() {}

bool CircularAlfven::initialize(void) {
   bool success = Project::initialize();

   creal m = physicalconstants::MASS_PROTON;
   creal e = physicalconstants::CHARGE;
   creal kB = physicalconstants::K_B;
   creal gamma = 5.0 / 3.0;
   creal mu0 = physicalconstants::MU_0;

   // Three reference values: BRef, nRef, mRef=m
   // All other reference values can be derived.
   lRef = sqrt(m / (mu0 * nRef)) / e;
   uRef = BRef / sqrt(mu0 * m * nRef);
   tRef = lRef / uRef;
   pRef = m * nRef * sqr(uRef);
   TRef = pRef / (kB * nRef);

   cosalpha = cos(alpha);
   sinalpha = sin(alpha);
   // Convert to SI units
   rho0 *= nRef * m;
   p0 *= pRef;
   Bkpar0 *= BRef;
   v1 *= uRef;
   // wave vector
   kwave = 2 * M_PI / (lambda * lRef);
   // Alfven speed
   creal VA = Bkpar0 / sqrt(mu0 * rho0);
   // Magnetic field perturbation amplitude
   B1 = v1 * Bkpar0 / VA;

   if (verbose) {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      if (myRank == MASTER_RANK) {
         if (P::ohmHallTerm == 0){
            std::cerr << "The Hall term is needed for accurate Alfven waves!";
         }
         std::cout << "wavelength = " << lambda << "\n";
         std::cout << "Lx / lRef = " << (P::xmax - P::xmin) / lRef << "\n";
         std::cout << "VA = " << VA << "\n";
         std::cout << "velocity perturbation amplitude = " << v1 << "\n";
         std::cout << "magnetic perturbation amplitude = " << B1 << "\n";
         std::cout << "angle = " << alpha / 3.1415926 * 180 << " [degree]\n";
         std::cout << "-------------------------------------------------\n";
         std::cout << "Unit conversion factors from dimensionless to SI:\n";
         std::cout << "lRef: " << lRef << " [m]\n";
         std::cout << "tRef: " << tRef << " [s]\n";
         std::cout << "vRef: " << uRef << " [m/s]\n";
         std::cout << "nRef: " << nRef << " [#]\n";
         std::cout << "TRef: " << TRef << " [K]\n";
         std::cout << "pRef: " << pRef << " [Pa]\n";
         std::cout << "BRef: " << BRef << " [T]\n";
         std::cout << "-------------------------------------------------\n";
         creal Lx = P::xmax - P::xmin;
         std::cout << "wavelength / di = " << Lx / lRef << "\n";
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

void CircularAlfven::addParameters() {
   typedef Readparameters RP;

   RP::add("CircularAlfven.v1", "Perturbation velocity amplitude", 0.1);
   RP::add("CircularAlfven.lambda", "Dimensionless wavelength", 32.0);
   RP::add("CircularAlfven.rho0", "Dimensionless density", 1.0);
   RP::add("CircularAlfven.p0", "Dimensionless pressure", 0.1);
   RP::add("CircularAlfven.Bkpar0", "Dimensionless parallel magnetic field to wave vector", 1.0);
   RP::add("CircularAlfven.alpha", "In-plane wave vector tilted angle w.r.t. +x in radian", 0.4636476090008061);
   // SI units inputs
   RP::add("CircularAlfven.BRef", "Reference magnetic field strength", 1e-8);
   RP::add("CircularAlfven.nRef", "Reference number density", 1e6);
   RP::add("CircularAlfven.verbose", "Turn on/off detailed information", 0);
}

void CircularAlfven::getParameters() {
   Project::getParameters();

   typedef Readparameters RP;
   RP::get("CircularAlfven.v1", v1);
   RP::get("CircularAlfven.lambda", lambda);
   RP::get("CircularAlfven.alpha", alpha);
   RP::get("CircularAlfven.rho0", rho0);
   RP::get("CircularAlfven.p0", p0);
   RP::get("CircularAlfven.Bkpar0", Bkpar0);
   RP::get("CircularAlfven.BRef", BRef);
   RP::get("CircularAlfven.nRef", nRef);
   RP::get("CircularAlfven.verbose", verbose);
}

Real CircularAlfven::getMaxwellian(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, creal& dvx,
                                   creal& dvy, creal& dvz, const uint popID) const {
   creal m = getObjectWrapper().particleSpecies[popID].mass;
   creal kB = physicalconstants::K_B;

   Real xpar = x * cosalpha + y * sinalpha;
   Real uperp = v1 * sin(kwave * xpar);
   Real n, ux, uy, uz, T;

   n = rho0 / m;
   ux = -uperp * sinalpha;
   uy = uperp * cosalpha;
   uz = v1 * cos(kwave * xpar);
   T = p0 / (n * kB);

   creal coef = m / (2 * M_PI * kB * T);
   // Maxwellian f(v) = f(v;n,u,T)
   creal f = n * sqrt(coef) * coef * exp(-coef * M_PI * (sqr(vx - ux) + sqr(vy - uy) + sqr(vz - uz)));

   return f;
}

Real CircularAlfven::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx,
                                           creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,
                                           const uint popID) const {
   Real f = 0.0;
   f = getMaxwellian(x + 0.5 * dx, y + 0.5 * dy, z + 0.5 * dz, vx + 0.5 * dvx, vy + 0.5 * dvy, vz + 0.5 * dvz, dvx, dvy,
                     dvz, popID);

   return f;
}

std::vector<std::array<Real, 3>> CircularAlfven::getV0(creal x, creal y, creal z, const uint popID) const {
   std::vector<std::array<Real, 3>> V0;
   std::array<Real, 3> v = {{0.0, 0.0, 0.0}};
   V0.push_back(v);
   return V0;
}

void CircularAlfven::calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) {}

void CircularAlfven::setProjectBField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                                      FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                                      FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
   // Background field
   ConstantField bgField;
   bgField.initialize(Bkpar0 * cosalpha, Bkpar0 * sinalpha, 0.0);
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
               Real Bperp = B1 * sin(kwave * xpar);
               cell->at(fsgrids::bfield::PERBX) = -Bperp * sinalpha;
               cell->at(fsgrids::bfield::PERBY) = Bperp * cosalpha;
               cell->at(fsgrids::bfield::PERBZ) = B1 * cos(kwave * xpar);
            }
         }
      }
   }
}

} // namespace projects
