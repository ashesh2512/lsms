#ifndef REPLICA_EXCHANGE_WL_H
#define REPLICA_EXCHANGE_WL_H

#include <mpi.h>

#include <random>

#include "Communication/REWLCommunication.hpp"

class REWL {
 public:
  REWL(MPI_Comm comm, int numberOfWindows, int numberOfWalkersPerWindow,
       unsigned rngSeed);
  ~REWL();

  void assignSwapPartner();

  void swapEnergy(Real &energyForSwap);

  bool determineAcceptance(Real myDOSRatio);

  void swapConfig(Real evecsForSwap[], int numElements);

  void swapPotentialShifts(Real potentialShiftForSwap[], int numElements);

 private:
  // Random number generator
  std::mt19937 rng;
  std::uniform_real_distribution<Real>
      rnd;  // rnd11(-1.0,1.0),rnd0pi(0.0,2.0*M_PI);

  REWLCommunication REWLcomm;
  int numWalkers;
  int numWindows;
  int numWalkersPerWindow;

  int myID;
  int myWindow;
  int partnerID;
  int partnerWindow;

  int swapDirection;

  int upExchanges;
  int downExchanges;
};

#endif
