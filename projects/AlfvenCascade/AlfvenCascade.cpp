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

#include <vector>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <random>

#include "../../backgroundfield/backgroundfield.h"
#include "../../backgroundfield/constantfield.hpp"
#include "../../common.h"
#include "../../object_wrapper.h"
#include "../../readparameters.h"

#include "AlfvenCascade.h"

using namespace spatial_cell;

// Structure to hold wave parameters
struct WaveParameters {
    Real wavelength;
    Real amplitude;
    Real phase;
    Real angle;
};

namespace projects {
AlfvenCascade::AlfvenCascade() : TriAxisSearch() {}
AlfvenCascade::~AlfvenCascade() {}

// Vector to store multiple waves
std::vector<WaveParameters> waves;

bool AlfvenCascade::initialize(void) {
   bool success = Project::initialize();

   creal m = physicalconstants::MASS_PROTON;
   creal e = physicalconstants::CHARGE;
   creal kB = physicalconstants::K_B;
   creal gamma = 5.0 / 3.0;
   creal mu0 = physicalconstants::MU_0;

   n = rho0 / m; // number density
   p0 = n * kB * T; // pressure
   
   // Initialize waves based on parameters
   waves.clear();
   for (size_t i = 0; i < wavelengths.size(); i++) {
       WaveParameters wave;
       wave.wavelength = wavelengths[i];
       wave.amplitude = amplitudes[i];
       wave.phase = phases[i];
       wave.angle = angles[i];
       waves.push_back(wave);
   }

   // Calculate Alfvén speed
   VA = B / sqrt(mu0 * rho0);

   if (verbose) {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      if (myRank == MASTER_RANK) {
         std::cout << "Initialized multi-wave turbulence simulation\n";
         std::cout << "Number of waves: " << waves.size() << "\n";
         std::cout << "Background field strength: " << B << " T\n";
         std::cout << "Alfvén speed: " << VA << " m/s\n";
         
         for (size_t i = 0; i < waves.size(); i++) {
             std::cout << "\nWave " << i + 1 << ":\n";
             std::cout << "Wavelength: " << waves[i].wavelength << " m\n";
             std::cout << "Amplitude: " << waves[i].amplitude << " m/s\n";
             std::cout << "Phase: " << waves[i].phase << " rad\n";
             std::cout << "Angle: " << waves[i].angle * 180/M_PI << " degrees\n";
         }
      }
   }

   return success;
}

void AlfvenCascade::addParameters() {
   typedef Readparameters RP;
   
   // For repeated parameters, just add the base parameter once
   RP::add("AlfvenCascade.numberOfWaves", "Number of waves in the simulation", 1);
   RP::add("AlfvenCascade.wavelength", "Wavelength of wave (m)", 32.0);
   RP::add("AlfvenCascade.amplitude", "Velocity amplitude (m/s)", 0.1);
   RP::add("AlfvenCascade.phase", "Initial phase (rad)", 0.0);
   RP::add("AlfvenCascade.angle", "Wave angle (rad)", 0.4636476090008061);
   
   RP::add("AlfvenCascade.rho0", "Background density (kg/m^3)", 1.6726219e-21);
   RP::add("AlfvenCascade.B", "Background magnetic field strength (T)", 1e-8);
   RP::add("AlfvenCascade.T", "Temperature (K)", 1e6);
   RP::add("AlfvenCascade.spectralIndex", "Power law index for initial spectrum", -5.0/3.0);
   RP::add("AlfvenCascade.randomSeed", "Seed for random phase generation", 12345);
   RP::add("AlfvenCascade.verbose", "Verbose output", 1);
}

void AlfvenCascade::getParameters() {
   typedef Readparameters RP;
   Project::getParameters();

   int nWaves;
   RP::get("AlfvenCascade.numberOfWaves", nWaves);

   // Clear vectors in case they have any default values
   wavelengths.clear();
   amplitudes.clear();
   phases.clear();
   angles.clear();

   // Read each parameter multiple times
   Real value;
   for (int i = 0; i < nWaves; i++) {
      RP::get("AlfvenCascade.wavelength", value);
      wavelengths.push_back(value);
      
      RP::get("AlfvenCascade.amplitude", value);
      amplitudes.push_back(value);
      
      RP::get("AlfvenCascade.phase", value);
      phases.push_back(value);
      
      RP::get("AlfvenCascade.angle", value);
      angles.push_back(value);
   }

   // Get scalar parameters
   RP::get("AlfvenCascade.rho0", rho0);
   RP::get("AlfvenCascade.B", B);
   RP::get("AlfvenCascade.T", T);
   RP::get("AlfvenCascade.spectralIndex", spectralIndex);
   RP::get("AlfvenCascade.randomSeed", randomSeed);
   RP::get("AlfvenCascade.verbose", verbose);
}

Real AlfvenCascade::getMaxwellian(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz, 
                                 creal& dvx, creal& dvy, creal& dvz, const uint popID) const {
   creal m = getObjectWrapper().particleSpecies[popID].mass;
   creal kB = physicalconstants::K_B;

   // Calculate total perturbation velocity from all waves
   Real ux = 0.0, uy = 0.0, uz = 0.0;
   
   for (const auto& wave : waves) {
       Real cosalpha = cos(wave.angle);
       Real sinalpha = sin(wave.angle);
       Real kwave = 2 * M_PI / wave.wavelength;
       Real xpar = x * cosalpha + y * sinalpha;
       
       Real uperp = wave.amplitude * sin(kwave * xpar + wave.phase);
       Real upara = wave.amplitude * cos(kwave * xpar + wave.phase);
       
       ux += -uperp * sinalpha;
       uy += uperp * cosalpha;
       uz += upara;
   }
   
   creal coef = m / (2 * M_PI * kB * T);
   creal f = n * sqrt(coef) * coef * exp(-coef * M_PI * (sqr(vx - ux) + sqr(vy - uy) + sqr(vz - uz)));

   return f;
}

Real AlfvenCascade::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx,
                                           creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,
                                           const uint popID) const {
   Real f = 0.0;
   f = getMaxwellian(x + 0.5 * dx, y + 0.5 * dy, z + 0.5 * dz, vx + 0.5 * dvx, vy + 0.5 * dvy, vz + 0.5 * dvz, dvx, dvy,
                     dvz, popID);

   return f;
}

std::vector<std::array<Real, 3>> AlfvenCascade::getV0(creal x, creal y, creal z, const uint popID) const {
   std::vector<std::array<Real, 3>> V0;
   std::array<Real, 3> v = {{0.0, 0.0, 0.0}};
   V0.push_back(v);
   return V0;
}

void AlfvenCascade::calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) {}

void AlfvenCascade::setProjectBField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                                    FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                                    FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
   // Set background field
   ConstantField bgField;
   bgField.initialize(B, 0.0, 0.0); // Background field in x-direction
   setBackgroundField(bgField, BgBGrid);

   if (!P::isRestart) {
      auto localSize = perBGrid.getLocalSize().data();

#pragma omp parallel for collapse(3)
      for (int i = 0; i < localSize[0]; ++i) {
         for (int j = 0; j < localSize[1]; ++j) {
            for (int k = 0; k < localSize[2]; ++k) {
               const std::array<Real, 3> x = perBGrid.getPhysicalCoords(i, j, k);
               std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(i, j, k);
               
               Real Bx = 0.0, By = 0.0, Bz = 0.0;
               
               // Sum contributions from all waves
               for (const auto& wave : waves) {
                   Real cosalpha = cos(wave.angle);
                   Real sinalpha = sin(wave.angle);
                   Real kwave = 2 * M_PI / wave.wavelength;
                   Real xpar = x[0] * cosalpha + x[1] * sinalpha;
                   creal mu0 = physicalconstants::MU_0;

                   // Calculate B1 from v1 using Alfvén wave relation
                   Real B1 = wave.amplitude * sqrt(mu0 * rho0);
                   
                   Real Bperp = B1 * sin(kwave * xpar + wave.phase);
                   Real Bpara = B1 * cos(kwave * xpar + wave.phase);
                   
                   Bx += -Bperp * sinalpha;
                   By += Bperp * cosalpha;
                   Bz += Bpara;
               }
               
               cell->at(fsgrids::bfield::PERBX) = Bx;
               cell->at(fsgrids::bfield::PERBY) = By;
               cell->at(fsgrids::bfield::PERBZ) = Bz;
            }
         }
      }
   }
}

} // namespace projects