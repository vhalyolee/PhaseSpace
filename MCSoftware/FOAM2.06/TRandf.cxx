using namespace std;
#include"TRandf.h"
#ifdef ROOT_DEF
ClassImp(TRanluxEngine)
ClassImp(TRMersenneTwister)
ClassImp(TRanmarEngine)
#endif

///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          TRandf                                           //
//         Universal package of the 3 random number generators               //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          TRanluxEngine                                    //
//                                                                           //
// The best available  random number generator of Luscher,
// with the user controlled "luxury level"
// See "A review of pseudorandom number generators", F. James
// Computer Physics Communications 60 (1990), pages 329-344.
//
///////////////////////////////////////////////////////////////////////////////
//                                                                           //
//                          TRMersenneTwister                                //
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
//
// A "fast, compact, huge-period generator" based on M. Matsumoto and 
// T. Nishimura, "Mersenne Twister: A 623-dimensionally equidistributed 
// uniform pseudorandom number generator", to appear in ACM Trans. on
// Modeling and Computer Simulation.  It is a twisted GFSR generator
// with a Mersenne-prime period of 2^19937-1, uniform on open interval (0,1)
//
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//                          TRanmarEngine                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//       PseudoRandom number generator used by Foam 2.x                    //
// Fortran version rewritten by S. Jadach, Nov 1997.                       //
// C++ translation by M. Slusarczyk, Feb 2000, and S. Jadach, Oct 2000	   // 
//                                                                         //
// Universal random number generator proposed by MARSAGLIA and ZAMAN       //
// in report FSU-SCRI-87-50                                                //
//        modified by F. James, 1988 and 1989, to generate a vector        //
//        of pseudorandom numbers rvec of length lenv, and to put in       //
//        the COMMON block everything needed to specify currrent state,    //
//        and to add input and output entry points rmarin, rmarut.         //
//
// ????????????????????? info below to be actualized ????????????????????????
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   //
// ++  CALLing sequences for TPSEMAR:                                 ++   //
// ++      TPSEMAR()          INITIALIZES the generator               ++   //
// ++      TPSEMAR(i1,n1,n2)  INITIALIZES the generator               ++   //
// ++                  32-bit integer i1, and number counts n1,n2     ++   //
// ++                  (for initializing, set n1=n2=0, but to restart ++   //
// ++                    a previously generated sequence, use values  ++   //
// ++                    output by rmarut)                            ++   //
// ++      Initialize(i1,n1,n2) used by the two above constructors    ++   //
// ++      FlatArray(lenv,RVec) GENERATES random vector RVec          ++   //
// ++      Out(i1,n1,n2)   outputs the value of the original          ++   //
// ++                  seed and the two number counts, to be used     ++   //
// ++                  for restarting by initializing to i1 and       ++   //
// ++                  skipping n2*100000000+n1 numbers.              ++   //
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++   //
//???????????????????????????????????????????????????????????????????????????
//                                                                         //
//         Initializing routine for ranmar, may be called before           //
//         generating pseudorandom numbers with ranmar. the input          //
//         values should be in the ranges:  0<=ijklin<=900 000 000         //
//                                          0<=ntotin<=999 999 999         //
//                                          0<=ntot2n<<999 999 999!        //
// to get the standard values in MARSAGLIA's paper, ijklin=54217137        //
//                                            ntotin,ntot2n=0              //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

//_____________________________________________________________________________
TRanluxEngine::TRanluxEngine(void)
: int_modulus(0x1000000),
  mantissa_bit_24( pow(0.5,24.) ),
  mantissa_bit_12( pow(0.5,12.) )
{
  SetSeedLux(314159265,3);
}

//_____________________________________________________________________________
TRanluxEngine::TRanluxEngine(long seed, int lux)
: int_modulus(0x1000000),
  mantissa_bit_24( pow(0.5,24.) ),
  mantissa_bit_12( pow(0.5,12.) )
 {
   luxury = lux;
   SetSeedLux(seed, luxury);
 }

//_____________________________________________________________________________
void TRanluxEngine::SetSeed(long seed) {
  // implements abstrac method
  int lux=3;
  SetSeedLux(seed, lux);
}

//_____________________________________________________________________________
void TRanluxEngine::SetSeedLux(long seed, int lux) {
// The initialisation is carried out using a Multiplicative
// Congruential generator using formula constants of L'Ecuyer 
// as described in "A review of pseudorandom number generators"
// (Fred James) published in Computer Physics Communications 60 (1990)
// pages 329-344
  const int ecuyer_a = 53668;
  const int ecuyer_b = 40014;
  const int ecuyer_c = 12211;
  const int ecuyer_d = 2147483563;
  const int lux_levels[5] = {0,24,73,199,365};  
  long int_seed_table[24];
  long next_seed = seed;
  long k_multiple;
  int i;
// number of additional random numbers that need to be 'thrown away'
// every 24 numbers is set using luxury level variable.
  theSeed = seed;
  if( (lux > 4)||(lux < 0) ){
     if(lux >= 24){
        nskip = lux - 24;
     }else{
        nskip = lux_levels[3]; // corresponds to default luxury level
     }
  }else{
     luxury = lux;
     nskip = lux_levels[luxury];
  }
  for(i = 0;i != 24;i++){
     k_multiple = next_seed / ecuyer_a;
     next_seed = ecuyer_b * (next_seed - k_multiple * ecuyer_a) 
                                       - k_multiple * ecuyer_c ;
     if(next_seed < 0) next_seed += ecuyer_d;
     int_seed_table[i] = next_seed % int_modulus;
  }     
  for(i = 0;i != 24;i++)
    float_seed_table[i] = int_seed_table[i] * mantissa_bit_24;
  i_lag = 23;
  j_lag = 9;
  carry = 0. ;
  if( float_seed_table[23] == 0. ) carry = mantissa_bit_24;
  count24 = 0;
}

//_____________________________________________________________________________
void TRanluxEngine::showStatus() 
{
    cout << endl;
    cout << "--------- Ranlux engine status ---------" << endl;
    cout << " Initial seed = " << theSeed << endl;
    cout << " float_seed_table[] = ";
    for (int i=0; i<24; ++i)
      cout << float_seed_table[i] << " ";
    cout << endl;
    cout << " i_lag = " << i_lag << ", j_lag = " << j_lag << endl;
    cout << " carry = " << carry << ", count24 = " << count24 << endl;
    cout << " luxury = " << luxury << " nskip = " << nskip << endl;
    cout << "----------------------------------------" << endl;
}

//_____________________________________________________________________________
double TRanluxEngine::Flat() {
   float next_random;
   float uni;
   int i;
   uni = float_seed_table[j_lag] - float_seed_table[i_lag] - carry;
   if(uni < 0. ){
      uni += 1.0;
      carry = mantissa_bit_24;
   }else{
      carry = 0.;
   }
   float_seed_table[i_lag] = uni;
   i_lag --;
   j_lag --;
   if(i_lag < 0) i_lag = 23;
   if(j_lag < 0) j_lag = 23;

   if( uni < mantissa_bit_12 ){
      uni += mantissa_bit_24 * float_seed_table[j_lag];
      if( uni == 0) uni = mantissa_bit_24 * mantissa_bit_24;
   }
   next_random = uni;
   count24 ++;

 // every 24th number generation, several random numbers are generated
 // and wasted depending upon the luxury level.

   if(count24 == 24 ){
      count24 = 0;
      for( i = 0; i != nskip ; i++){
          uni = float_seed_table[j_lag] - float_seed_table[i_lag] - carry;
          if(uni < 0. ){
             uni += 1.0;
             carry = mantissa_bit_24;
          }else{
             carry = 0.;
          }
          float_seed_table[i_lag] = uni;
          i_lag --;
          j_lag --;
          if(i_lag < 0)i_lag = 23;
          if(j_lag < 0) j_lag = 23;
       }
   }
   return (double) next_random;
 }

//_____________________________________________________________________________
void TRanluxEngine::FlatArray(const int size, double* vect)
 {
   float next_random;
   float uni;
   int i;
   int index;

   for (index=0; index<size; ++index) {
     uni = float_seed_table[j_lag] - float_seed_table[i_lag] - carry;
     if(uni < 0. ){
        uni += 1.0;
        carry = mantissa_bit_24;
     }else{
        carry = 0.;
     }

     float_seed_table[i_lag] = uni;
     i_lag --;
     j_lag --;
     if(i_lag < 0) i_lag = 23;
     if(j_lag < 0) j_lag = 23;

     if( uni < mantissa_bit_12 ){
        uni += mantissa_bit_24 * float_seed_table[j_lag];
        if( uni == 0) uni = mantissa_bit_24 * mantissa_bit_24;
     }
     next_random = uni;
     vect[index] = (double)next_random;
     count24 ++;

 // every 24th number generation, several random numbers are generated
 // and wasted depending upon the luxury level.
     if(count24 == 24 ){
        count24 = 0;
        for( i = 0; i != nskip ; i++){
            uni = float_seed_table[j_lag] - float_seed_table[i_lag] - carry;
            if(uni < 0. ){
               uni += 1.0;
               carry = mantissa_bit_24;
            }else{
               carry = 0.;
            }
            float_seed_table[i_lag] = uni;
            i_lag --;
            j_lag --;
            if(i_lag < 0)i_lag = 23;
            if(j_lag < 0) j_lag = 23;
         }
     }
   }
 }
/////////////////////////////////////////////////////////////////////////////
//                      End of TRanluxEngine                               //
/////////////////////////////////////////////////////////////////////////////




//_____________________________________________________________________________
TRMersenneTwister::TRMersenneTwister() {
  powersOfTwo();
  SetSeed(0);
}

//_____________________________________________________________________________
TRMersenneTwister::TRMersenneTwister(long seed) {
  powersOfTwo();
  SetSeed(seed);
}

//_____________________________________________________________________________
void TRMersenneTwister::powersOfTwo() {
  twoToMinus_32 = ldexp (1.0, -32);
  twoToMinus_53 = ldexp (1.0, -53);
  nearlyTwoToMinus_54 = ldexp (1.0, -54) - ldexp (1.0, -100);
}

//_____________________________________________________________________________
void TRMersenneTwister::SetSeed(long seed) {
  theSeed = seed ? seed : 4357;
  mt[0] = (unsigned int)theSeed;
  int i;
  for( i=1; i < 624; ++i ) {
    mt[i] = (69069 * mt[i-1]) & 0xffffffff;
  }
  count624=N;
}

//_____________________________________________________________________________
double TRMersenneTwister::Flat() {
  unsigned int y;
  if( count624 >= N ) {
    register int i;
    for( i=0; i < NminusM; ++i ) {
      y = (mt[i] & 0x80000000) | (mt[i+1] & 0x7fffffff);
      mt[i] = mt[i+M]       ^ (y >> 1) ^ ((y & 0x1) ? 0x9908b0df : 0x0 );
    }
    for(    ; i < N-1    ; ++i ) {
      y = (mt[i] & 0x80000000) | (mt[i+1] & 0x7fffffff);
      mt[i] = mt[i-NminusM] ^ (y >> 1) ^ ((y & 0x1) ? 0x9908b0df : 0x0 );
    }
    y = (mt[i] & 0x80000000) | (mt[0] & 0x7fffffff);
    mt[i] = mt[M-1] ^ (y >> 1) ^ ((y & 0x1) ? 0x9908b0df : 0x0 );
    count624 = 0;
  }
  y = mt[count624];
  y ^= ( y >> 11);
  y ^= ((y << 7 ) & 0x9d2c5680);
  y ^= ((y << 15) & 0xefc60000);
  y ^= ( y >> 18);

  return                      y * twoToMinus_32  +    // Scale to range 
         (mt[count624++] >> 11) * twoToMinus_53  +    // fill remaining bits
                	    nearlyTwoToMinus_54;      // make sure non-zero
}

//_____________________________________________________________________________
void TRMersenneTwister::FlatArray( const int size, double *vect ) {
  for( int i=0; i < size; ++i) vect[i] = Flat();
}

//_____________________________________________________________________________
void TRMersenneTwister::ShowStatus() 
{
   cout << endl;
   cout << "--------- MTwist engine status ---------" << endl;
   cout <<  setprecision(20);
   cout << " Initial seed      = " << theSeed <<  endl;
   cout << " Current index     = " << count624 <<  endl;
   cout << " Array status mt[] = " <<  endl;
   for (int i=0; i<624; i+=5) {
      cout << mt[i]   << " " << mt[i+1] << " " << mt[i+2] << " " 
	       << mt[i+3] << " " << mt[i+4] <<  endl;
   }
    cout << "----------------------------------------" <<  endl;
}
/////////////////////////////////////////////////////////////////////////////
//                        End of TRMersenneTwister                         //
/////////////////////////////////////////////////////////////////////////////




//_____________________________________________________________________________
TRanmarEngine::TRanmarEngine(void){
  //   Constructor without internal seed
  m_modcns = 1000000000;
  long ijkl_new  = 54217137;
  long ntot_new  = 0;
  long ntot2_new = 0;
  Initialize(ijkl_new, ntot_new,ntot2_new);
}// TRanmarEngine

//_____________________________________________________________________________
TRanmarEngine::TRanmarEngine(long ijkl_new){
  // Constructor with one external seed
  long ntot_new  = 0;
  long ntot2_new = 0;
  Initialize(ijkl_new, ntot_new,ntot2_new);
}

//_____________________________________________________________________________
TRanmarEngine::TRanmarEngine(long ijkl_new, long ntot_new, long ntot2_new){
  // Constructor with original external seed and series length
  Initialize(ijkl_new, ntot_new,ntot2_new);
}

//_____________________________________________________________________________
TRanmarEngine::~TRanmarEngine(){
  // Destructor
}

//_____________________________________________________________________________
void TRanmarEngine::SetSeed(long ijkl_new){
  // (re-)start from one external seed
  long ntot_new  = 0;
  long ntot2_new = 0;
  Initialize(ijkl_new, ntot_new,ntot2_new);
}

//_____________________________________________________________________________
void TRanmarEngine::Initialize(long ijkl_new, long ntot_new, long ntot2_new){
  // Initialize or RE-initialize random number generator
  double t,uni,s;
  long m,i24,jj,idum,loop2,now;
  long i,j,ij,k,l,ii,kl;
  
  m_ijkl = ijkl_new;
  m_ntot = max(ntot_new,0);
  m_ntot2= max(ntot2_new,0);
  
  ij = m_ijkl/30082;
  kl = m_ijkl - 30082*ij;
  i = MOD(ij/177, 177) + 2;
  j = MOD(ij, 177)     + 2;
  k = MOD(kl/169, 178) + 1;
  l = MOD(kl, 169);
  
  for( ii= 1; ii<=97 ; ii++)
    {
      s = 0.0;
      t = 0.5;
      for(jj= 1; jj<= 24 ;jj++)
	{
	  m = MOD(MOD(i*j,179)*k, 179);
	  i = j;
	  j = k;
	  k = m;
	  l = MOD(53*l+1, 169);
	  if (MOD(l*m,64) >=  32)  s = s+t;
	  t = 0.5*t;
	}
      m_U[ii] = s;
    }
  
  m_twom24 = 1.0;
  
  for(i24 = 1 ; i24 <= 24 ; i24++)
    m_twom24 = 0.5*m_twom24;
  
  m_C  =   362436.*m_twom24;
  m_CD =  7654321.*m_twom24;
  m_CM = 16777213.*m_twom24;
  m_i97 = 97;
  m_j97 = 33;
  
  // complete initialization by skipping
  // (ntot2*m_modcns + ntot) random numbers  
  for(loop2= 1; loop2<= m_ntot2+1 ; loop2++)
    {
      now = m_modcns;
      if (loop2  ==  m_ntot2+1)  now=m_ntot;
      if (now  >  0)
	{
	  for(idum = 1; idum <= m_ntot ; idum++)
	    {
	      uni = m_U[m_i97]-m_U[m_j97];
	      if (uni  <  0.)  uni=uni+1.0;
	      m_U[m_i97] = uni;
	      m_i97 = m_i97-1;
	      if (m_i97  ==  0)  m_i97=97;
	      m_j97 = m_j97-1;
	      if (m_j97  ==  0)  m_j97=97;
	      m_C = m_C - m_CD;
	      if (m_C  <  0.)   m_C=m_C+m_CM;
	    }
	}
    }
}// Initialize

//_____________________________________________________________________________
double TRanmarEngine::Flat(){
  //  Generate single random numbers
  double RVec[1];
  int lenv=1;
  FlatArray(lenv, RVec);
  return RVec[0];
}

//_____________________________________________________________________________
void TRanmarEngine::FlatArray(const int lenv, double *RVec){
  //  Generate lenv random numbers
  double   zuni,uni;
  int ivec;
  // Normal entry to generate lenv random numbers
  for( ivec= 0; ivec < lenv; ivec++)
    {
      uni = m_U[m_i97]-m_U[m_j97];
      if (uni  <  0.)  uni=uni+1.0;
      m_U[m_i97] = uni;
      m_i97 = m_i97-1;
      if (m_i97  ==  0)  m_i97=97;
      m_j97 = m_j97-1;
      if (m_j97  ==  0)  m_j97=97;
      m_C = m_C - m_CD;
      if (m_C  <  0.)   m_C=m_C+m_CM;
      uni = uni-m_C;
      if (uni  <  0.) uni=uni+1.0;
      *(RVec +ivec) = uni;
      // Replace exact zeros by uniform distr. *2**-24
      if (uni  ==  0.)
	{
	  zuni = m_twom24*m_U[2];
	  // an exact zero here is very unlikely, but let's be safe.
	  if (zuni  ==  0.) zuni= m_twom24*m_twom24;
	  *(RVec +ivec) = zuni;
	}
    }
  m_ntot  = m_ntot + lenv;
  if (m_ntot  >=  m_modcns)
    {
      m_ntot2  =  m_ntot2 + 1;
      m_ntot   =  m_ntot - m_modcns;
    }
}//FlatArray

//_____________________________________________________________________________
void TRanmarEngine::ShowStatus()
{
    cout <<  endl;
    cout << "----- HepJamesRandom engine status -----" <<  endl;
    cout << " Initial seed = " << m_ijkl <<  endl;
    cout << " u[] = ";
   for (int i=0; i<97; ++i)
      cout << m_U[i] << " ";
    cout <<  endl;
    cout << " c = " << m_C << ", cd = " << m_CD << ", cm = " << m_CM
             <<  endl;
    cout << " i97 = " << m_i97 << ", u[i97] = " << m_U[m_i97] <<  endl;
    cout << " j97 = " << m_j97 << ", u[j97] = " << m_U[m_j97] <<  endl;
    cout << "----------------------------------------" <<  endl;
}
/////////////////////////////////////////////////////////////////////////////
//          End of class  TRanmarEngine                                    //
/////////////////////////////////////////////////////////////////////////////
