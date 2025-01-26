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

#include "Alfvencascadev2.h"

using namespace spatial_cell;

namespace projects {
Alfvencascadev2::Alfvencascadev2() : TriAxisSearch() {}
Alfvencascadev2::~Alfvencascadev2() {}

bool Alfvencascadev2::initialize(void) {
   bool success = Project::initialize();

   // ======== PERRONE 2013 PARAMETERS ========
   creal m_proton = physicalconstants::MASS_PROTON;
   creal mu0 = physicalconstants::MU_0;
   
   // Plasma parameters
   rho0 = 1.6726219e-21;              // Proton density ~1e6 m⁻³
   T_proton = 1e6;                    // Proton temperature (K)
   n_proton = rho0 / m_proton;        // Number density
   n_alpha = 0.05 * n_proton;         // 5% alpha density
   T_alpha = T_proton;                // Initial equal temperature
   B0 = 1e-8;                         // Background magnetic field (T)
   eta = 2e-4;                        // Resistivity (Ω·m)
   
   // Derived parameters
   creal VA = B0 / sqrt(mu0 * rho0);  // Alfvén speed
   creal d_p = sqrt(m_proton/(mu0 * n_proton * pow(physicalconstants::CHARGE,2))); // Proton skin depth
   lambda = 20*M_PI*d_p;              // Domain size ~20πd_p

   // Turbulent spectrum parameters
   k_min = 0.1/d_p;                   // Normalized k_min = 0.1
   k_max = 0.3/d_p;                   // Normalized k_max = 0.3
   delta_B = 0.3 * B0;                // Perturbation amplitude

   // Initialize random phases
   std::mt19937 gen(randomSeed);
   std::uniform_real_distribution<Real> phase_dist(0, 2*M_PI);
   
   // Generate wave spectrum
   waves.clear();
   for (Real k = k_min; k <= k_max; k += (k_max-k_min)/10) {
       WaveParameters wave;
       wave.wavelength = 2*M_PI/k;    // Physical wavelength
       wave.amplitude = delta_B / (B0 * sqrt(mu0 * rho0)) * VA; // δB/B0=0.3
       wave.phase = phase_dist(gen);  // Random phase
       wave.angle = atan2(1.0, 0.0);  // Perpendicular propagation
       waves.push_back(wave);
   }

   if (verbose) {
      int myRank;
      MPI_Comm_rank(MPI_COMM_WORLD, &myRank);
      if (myRank == MASTER_RANK) {
         std::cout << "===== PERRONE 2013 SETUP =====\n";
         std::cout << "Domain size: " << lambda << " m\n";
         std::cout << "Proton skin depth: " << d_p << " m\n";
         std::cout << "Alfvén speed: " << VA << " m/s\n";
         std::cout << "Alpha density: " << n_alpha << " m⁻³\n";
         std::cout << "Resistivity: " << eta << " Ω·m\n";
         std::cout << "Number of waves: " << waves.size() << "\n";
      }
   }

   return success;
}

void Alfvencascadev2::addParameters() {
   typedef Readparameters RP;
   
   // ======== CORE PARAMETERS ========
   RP::add("Alfvencascadev2.rho0", "Proton mass density (kg/m³)", 1.6726219e-21);
   RP::add("Alfvencascadev2.B", "Background B field (T)", 1e-8);
   RP::add("Alfvencascadev2.T", "Temperature (K)", 1e6);
   RP::add("Alfvencascadev2.eta", "Resistivity (Ω·m)", 2e-4);
   RP::add("Alfvencascadev2.randomSeed", "Random seed", 12345);
   RP::add("Alfvencascadev2.verbose", "Verbose output", 1);
}

void Alfvencascadev2::getParameters() {
   typedef Readparameters RP;
   Project::getParameters();

   RP::get("Alfvencascadev2.rho0", rho0);
   RP::get("Alfvencascadev2.B", B0);
   RP::get("Alfvencascadev2.T", T_proton);
   RP::get("Alfvencascadev2.eta", eta);
   RP::get("Alfvencascadev2.randomSeed", randomSeed);
   RP::get("Alfvencascadev2.verbose", verbose);
}

Real Alfvencascadev2::getMaxwellian(creal& x, creal& y, creal& z, creal& vx, creal& vy, creal& vz,
                                 creal& dvx, creal& dvy, creal& dvz, const uint popID) const {
   // ======== MULTI-ION DISTRIBUTIONS ========
   Real mass, density, temperature;
   if (popID == 0) {  // Protons
      mass = physicalconstants::MASS_PROTON;
      density = rho0 / mass;
      temperature = T_proton;
   } else {           // Alpha particles
      mass = 4*physicalconstants::MASS_PROTON;
      density = n_alpha;
      temperature = T_alpha;
   }

   // Thermal speed
   creal v_th = sqrt(2 * physicalconstants::K_B * temperature / mass);
   
   // Turbulent velocity perturbations
   Real ux = 0.0, uy = 0.0, uz = 0.0;
   for (const auto& wave : waves) {
       Real cosalpha = cos(wave.angle);
       Real sinalpha = sin(wave.angle);
       Real kwave = 2*M_PI/wave.wavelength;
       Real xpar = x*cosalpha + y*sinalpha;
       
       ux += -wave.amplitude * sin(kwave*xpar + wave.phase) * sinalpha;
       uy += wave.amplitude * sin(kwave*xpar + wave.phase) * cosalpha;
       uz += wave.amplitude * cos(kwave*xpar + wave.phase);
   }

   // Bi-Maxwellian distribution
   creal coef = mass / (2*M_PI*physicalconstants::K_B*temperature);
   return density * pow(coef, 1.5) * exp(-coef*(pow(vx-ux,2) + pow(vy-uy,2) + pow(vz-uz,2)));
}

Real Alfvencascadev2::calcPhaseSpaceDensity(creal& x, creal& y, creal& z, creal& dx, creal& dy, creal& dz, creal& vx,
                                           creal& vy, creal& vz, creal& dvx, creal& dvy, creal& dvz,
                                           const uint popID) const {
   return getMaxwellian(x+0.5*dx, y+0.5*dy, z+0.5*dz, 
                        vx+0.5*dvx, vy+0.5*dvy, vz+0.5*dvz,
                        dvx, dvy, dvz, popID);
}

void Alfvencascadev2::setProjectBField(FsGrid<std::array<Real, fsgrids::bfield::N_BFIELD>, FS_STENCIL_WIDTH>& perBGrid,
                                    FsGrid<std::array<Real, fsgrids::bgbfield::N_BGB>, FS_STENCIL_WIDTH>& BgBGrid,
                                    FsGrid<fsgrids::technical, FS_STENCIL_WIDTH>& technicalGrid) {
   // Background field
   ConstantField bgField;
   bgField.initialize(B0, 0.0, 0.0);  // B0 along x-axis
   setBackgroundField(bgField, BgBGrid);

   if (!P::isRestart) {
      auto localSize = perBGrid.getLocalSize().data();
      creal mu0 = physicalconstants::MU_0;

#pragma omp parallel for collapse(3)
      for (int i = 0; i < localSize[0]; ++i) {
         for (int j = 0; j < localSize[1]; ++j) {
            for (int k = 0; k < localSize[2]; ++k) {
               const std::array<Real, 3> xyz = perBGrid.getPhysicalCoords(i, j, k);
               std::array<Real, fsgrids::bfield::N_BFIELD>* cell = perBGrid.get(i, j, k);

               // Initialize turbulent fluctuations
               Real Bx = 0.0, By = 0.0, Bz = 0.0;
               for (const auto& wave : waves) {
                   Real cosalpha = cos(wave.angle);
                   Real sinalpha = sin(wave.angle);
                   Real kwave = 2*M_PI/wave.wavelength;
                   Real xpar = xyz[0]*cosalpha + xyz[1]*sinalpha;

                   Real B1 = wave.amplitude * sqrt(mu0 * rho0);
                   Bx += -B1 * sin(kwave*xpar + wave.phase) * sinalpha;
                   By += B1 * sin(kwave*xpar + wave.phase) * cosalpha;
                   Bz += B1 * cos(kwave*xpar + wave.phase);
               }

               // Apply resistive term
               cell->at(fsgrids::bfield::PERBX) = Bx + eta*currentX;
               cell->at(fsgrids::bfield::PERBY) = By + eta*currentY;
               cell->at(fsgrids::bfield::PERBZ) = Bz;
            }
         }
      }
   }
}

// Remaining functions unchanged from original implementation
void Alfvencascadev2::calcCellParameters(spatial_cell::SpatialCell* cell, creal& t) {}
std::vector<std::array<Real, 3>> Alfvencascadev2::getV0(creal x, creal y, creal z, const uint popID) const {
   return {std::array<Real,3>{0,0,0}};
}

} // namespace projects