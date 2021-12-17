/*****************************************************************************/
/***************** Copyright (C) 2020-2021, Richard Spinney. *****************/
/*****************************************************************************/
//                                                                           //
//    This program is free software: you can redistribute it and/or modify   //
//    it under the terms of the GNU General Public License as published by   //
//    the Free Software Foundation, either version 3 of the License, or      //
//    (at your option) any later version.                                    //
//                                                                           //
//    This program is distributed in the hope that it will be useful,        //
//    but WITHOUT ANY WARRANTY; without even the implied warranty of         //
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          //
//    GNU General Public License for more details.                           //
//                                                                           //
//    You should have received a copy of the GNU General Public License      //
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.  //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////

#ifndef DEFINITIONS_H
#define DEFINITIONS_H

//Optimisaton constants
#define MAXIMUM_ITERATIONS 1000
#define ERR_TOL 1E-6

//Box constraints for the optimisation
#define MIN_ALPHA -10.0
#define MAX_ALPHA 10.0
#define MIN_GAMMA -500.0
#define MAX_GAMMA 500.0
#define MIN_KAPPA 0.0
#define MAX_KAPPA 1e30

//Default number of threads in case std::thread::hardware_concurrency() fails
#define DEFAULT_THREADS 4

//Initial Seed value
#define SEED_BASE 1234

#endif /*DEFINITIONS_H*/
