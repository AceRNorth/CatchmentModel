#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <sstream>
#include <math.h>
#include <time.h>
#include <vector>
#include <assert.h>
#include <stdio.h>
#include <ctime>
#include <cstdlib> // for exit function

#ifndef RANDOMC_H
#define RANDOMC_H

using namespace std;
//using std::ifstream;

const double PI=3.14159265;

	struct totals
	{
		int W;
		int H;
		int E;
		int NumPat;
		double GseS;
		double GsiI;
		double GieI;
		double totalruntime;
	};	

	struct initials
	{
		int W,H;
		double heg_time;
		string inputfile;
	};	
			
	struct Patch
	{
		double x;
		double y;
		double GseS;
		double GsiI;
		double GieI;
		char type;
		vector<int> connecIND;
		vector<double> connecW;
		double TotW;
		int group;
		int groupindex;
	};

	struct Vil
	{
		double x;
		double y;
	};

//	struct 	PatchArray {vector<Patch> Pat[nx][ny];};
//	struct 	VilArray {vector<Vil> V[nx][ny];};

	struct PatchGroup
	{
		vector<int> members;
		int NumPat;
		int W;
		int H;
		int E;
		double GseS;
		double GsiI;
		double GieI;
	};

	struct Times
	{
	double interval;
	double maxT;
	int MoveRec;
	double rec;
	int N;
	};

	struct Pars
	{
	double m,delta,b,xw,xh,LD,U;
	int set;
	int ngroups;
	};

void RunOnceInt(double);

int random_poisson(double);
double random_exp(double);
double chop(double);
		void record(int);
		double dist(double,double,double,double);
		double distB(double,double,double,double,double);
		double abso(double);
		void RunMaxT(int);
	void EstRStar();
		void RunNReps(int);
		void initiate();
		void ColPressure(int,int,int);
		
		void UpdateConnec();
		void VilCreate();
		void WPatCreate();
		void HPatCreate();
		void PatchDestroy();

		void WtoH();
		void WtoE();
		void HtoE();
		void EtoW();
		void EtoH();
		void PutW();
		void PutH();

		double OneStep();
		
		int pick(double*,int,double);
		int pick(int*,int,int);
		


// Define 32 bit signed and unsigned integers.
// GieIange these definitions, if necessary, to match a particular platform
#if defined(_WIN16) || defined(__MSDOS__) || defined(_MSDOS)
   // 16 bit systems use long int for 32 bit integer
   typedef long int           int32;   // 32 bit signed integer
   typedef unsigned long int  uint32;  // 32 bit unsigned integer
#else
   // Most other systems use int for 32 bit integer
   typedef int                int32;   // 32 bit signed integer
   typedef unsigned int       uint32;  // 32 bit unsigned integer
#endif

// Define 64 bit signed and unsigned integers, if possible
#if (defined(__WINDOWS__) || defined(_WIN32)) && (defined(_MSC_VER) || defined(__INTEL_COMPILER))
   // Microsoft and other compilers under Windows use __int64
   typedef __int64            int64;   // 64 bit signed integer
   typedef unsigned __int64   uint64;  // 64 bit unsigned integer
   #define INT64_DEFINED               // Remember that int64 is defined
#elif defined(__unix__) && (defined(_M_IX86) || defined(_M_X64))
   // Gnu and other compilers under Linux etc. use long long
   typedef long long          int64;   // 64 bit signed integer
   typedef unsigned long long uint64;  // 64 bit unsigned integer
   #define INT64_DEFINED               // Remember that int64 is defined
#else
   // 64 bit integers not defined
   // You may include definitions for other platforms here
#endif


/***********************************************************************
System-specific user interface functions
***********************************************************************/

void EndOfProgram(void);               // System-specific exit code (userintf.cpp)

void FatalError(char * ErrorText);     // System-specific error reporting (userintf.cpp)


/***********************************************************************
Define random number generator classes
***********************************************************************/

class CRandomMersenne {                // Encapsulate random number generator
#if 0
   // Define constants for type MT11213A:
#define MERS_N   351
#define MERS_M   175
#define MERS_R   19
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   17
#define MERS_A   0xE4BD75F5
#define MERS_B   0x655E5280
#define MERS_C   0xFFD58000
#else
   // or constants for type MT19937:
#define MERS_N   624
#define MERS_M   397
#define MERS_R   31
#define MERS_U   11
#define MERS_S   7
#define MERS_T   15
#define MERS_L   18
#define MERS_A   0x9908B0DF
#define MERS_B   0x9D2C5680
#define MERS_C   0xEFC60000
#endif
public:
   CRandomMersenne(uint32 seed) {      // Constructor
      RandomInit(seed); LastInterval = 0;}
   void RandomInit(uint32 seed);       // Re-seed
   void RandomInitByArray(uint32 seeds[], int length); // Seed by more than 32 bits
   int IRandom (int min, int max);     // Output random integer
   int IRandomX(int min, int max);     // Output random integer, exact
   double Random();                    // Output random float
   uint32 BRandom();                   // Output random bits
private:
   void Init0(uint32 seed);            // Basic initialization procedure
   uint32 mt[MERS_N];                  // State vector
   int mti;                            // Index into mt
   uint32 LastInterval;                // Last interval length for IRandomX
   uint32 RLimit;                      // Rejection limit used by IRandomX
   enum TArch {LITTLE_ENDIAN1, BIG_ENDIAN1, NONIEEE}; // Definition of architecture
   TArch Architecture;                 // Conversion to float depends on architecture
};


class CRandomMother {             // Encapsulate random number generator
public:
   void RandomInit(uint32 seed);       // Initialization
   int IRandom(int min, int max);      // Get integer random number in desired interval
   double Random();                    // Get floating point random number
   uint32 BRandom();                   // Output random bits
   CRandomMother(uint32 seed) {   // Constructor
      RandomInit(seed);}
protected:
   uint32 x[5];                        // History buffer
};

#endif

// body file: mersenne.cpp
