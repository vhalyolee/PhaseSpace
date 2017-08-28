#include <stdlib.h>
#include <stdio.h>

#include "sobol.h"

#ifdef __cplusplus
extern "C" {
#endif

typedef int Cbool;

#define true 1
#define false 0

//******************************************************************************

double d_uniform_01 ( int *seed )

//******************************************************************************
//
//  Purpose:
//
//    D_UNIFORM_01 returns a unit double precision pseudorandom number.
//
//  Discussion:
//
//    This routine implements the recursion
//
//      seed = 16807 * seed mod ( 2**31 - 1 )
//      d_uniform_01 = seed / ( 2**31 - 1 )
//
//    The integer arithmetic never requires more than 32 bits,
//    including a sign bit.
//
//    If the initial seed is 12345, then the first three computations are
//
//      Input     Output      D_UNIFORM_01
//      SEED      SEED
//
//         12345   207482415  0.096616
//     207482415  1790989824  0.833995
//    1790989824  2035175616  0.947702
//
//  Modified:
//
//    11 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Paul Bratley, Bennett Fox, Linus Schrage,
//    A Guide to Simulation,
//    Springer Verlag, pages 201-202, 1983.
//
//    Pierre L'Ecuyer,
//    Random Number Generation,
//    in Handbook of Simulation
//    edited by Jerry Banks,
//    Wiley Interscience, page 95, 1998.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    P A Lewis, A S Goodman, J M Miller,
//    A Pseudo-Random Number Generator for the System/360,
//    IBM Systems Journal,
//    Volume 8, pages 136-143, 1969.
//
//  Parameters:
//
//    Input/output, int *SEED, the "seed" value.  Normally, this
//    value should not be 0.  On output, SEED has been updated.
//
//    Output, double D_UNIFORM_01, a new pseudorandom variate, strictly between
//    0 and 1.
//
{
  int k;
  double r;

  k = *seed / 127773;

  *seed = 16807 * ( *seed - k * 127773 ) - k * 2836;

  if ( *seed < 0 )
  {
    *seed = *seed + 2147483647;
  }
//
//  Although SEED can be represented exactly as a 32 bit integer,
//  it generally cannot be represented exactly as a 32 bit real number!
//
  r = ( double ) ( *seed ) * 4.656612875E-10;

  return r;
}
//******************************************************************************

int i4_bit_hi1 ( int n )

//******************************************************************************
//
//  Purpose:
//
//    I4_BIT_HI1 returns the position of the high 1 bit base 2 in an integer.
//
//  Example:
//
//       N    Binary    Hi 1
//    ----    --------  ----
//       0           0     0
//       1           1     1
//       2          10     2
//       3          11     2 
//       4         100     3
//       5         101     3
//       6         110     3
//       7         111     3
//       8        1000     4
//       9        1001     4
//      10        1010     4
//      11        1011     4
//      12        1100     4
//      13        1101     4
//      14        1110     4
//      15        1111     4
//      16       10000     5
//      17       10001     5
//    1023  1111111111    10
//    1024 10000000000    11
//    1025 10000000001    11
//
//  Modified:
//
//    13 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the integer to be measured.
//    N should be nonnegative.  If N is nonpositive, I4_BIT_HI1
//    will always be 0.
//
//    Output, int I4_BIT_HI1, the location of the high order bit.
//
{
  int bit;

  bit = 0;

  while ( 0 < n )
  {
    bit = bit + 1;
    n = n / 2;
  }

  return bit;
}
//******************************************************************************

int i4_bit_lo0 ( int n )

//******************************************************************************
//
//  Purpose:
//
//    I4_BIT_LO0 returns the position of the low 0 bit base 2 in an integer.
//
//  Example:
//
//       N    Binary    Lo 0
//    ----    --------  ----
//       0           0     1
//       1           1     2
//       2          10     1
//       3          11     3 
//       4         100     1
//       5         101     2
//       6         110     1
//       7         111     4
//       8        1000     1
//       9        1001     2
//      10        1010     1
//      11        1011     3
//      12        1100     1
//      13        1101     2
//      14        1110     1
//      15        1111     5
//      16       10000     1
//      17       10001     2
//    1023  1111111111     1
//    1024 10000000000     1
//    1025 10000000001     1
//
//  Modified:
//
//    13 March 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int N, the integer to be measured.
//    N should be nonnegative.
//
//    Output, int I4_BIT_LO0, the position of the low 1 bit.
//
{
  int bit;
  int n2;

  bit = 0;

  while ( true )
  {
    bit = bit + 1;
    n2 = n / 2;

    if ( n == 2 * n2 )
    {
      break;
    }

    n = n2;

  }

  return bit;
}
//******************************************************************************

void i4_sobol ( int dim_num, int *seed, float quasi[ ] )

//******************************************************************************
//
//  Purpose:
//
//    I4_SOBOL generates a new quasirandom Sobol vector with each call.
//
//  Discussion:
//
//    The routine adapts the ideas of Antonov and Saleev.
//
//    This routine uses INT for integers and FLOAT for real values.
//
//  Modified:
//
//    05 July 2006
//
//  Author:
//
//    Bennett Fox
//
//    C++ translation by John Burkardt
//
//  Reference:
//
//    IA Antonov, VM Saleev,
//    USSR Computational Mathematics and Mathematical Physics,
//    Volume 19, 1980, pages 252 - 256.
//
//    Paul Bratley and Bennett Fox,
//    Algorithm 659:
//    Implementing Sobol's Quasirandom Sequence Generator,
//    ACM Transactions on Mathematical Software,
//    Volume 14, Number 1, pages 88-100, 1988.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom 
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Stephen Joe, Frances Kuo
//    Remark on Algorithm 659:
//    Implementing Sobol's Quasirandom Sequence Generator,
//    ACM Transactions on Mathematical Software,
//    Volume 29, Number 1, pages 49-57, March 2003.
//
//    Ilya Sobol,
//    USSR Computational Mathematics and Mathematical Physics,
//    Volume 16, pages 236-242, 1977.
//
//    Ilya Sobol, YL Levitan, 
//    The Production of Points Uniformly Distributed in a Multidimensional 
//    Cube (in Russian),
//    Preprint IPM Akad. Nauk SSSR, 
//    Number 40, Moscow 1976.
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//    DIM_NUM must satisfy 2 <= DIM_NUM <= 1111.
//
//    Input/output, int *SEED, the "seed" for the sequence.
//    This is essentially the index in the sequence of the quasirandom
//    value to be generated.  On output, SEED has been set to the
//    appropriate next value, usually simply SEED+1.
//    If SEED is less than 0 on input, it is treated as though it were 0.
//    An input value of 0 requests the first (0-th) element of the sequence.
//
//    Output, float QUASI[DIM_NUM], the next quasirandom vector.
//
{
# define DIM_MAX 40
# define DIM_MAX2 1111
# define LOG_MAX 30

  static int atmost;
  static int dim_num_save = 0;
  int i;
  Cbool includ[LOG_MAX];
  static Cbool initialized = false;
  int j;
  int j2;
  int k;
  int l = 0;
  static int lastq[DIM_MAX2];
  int m;
  static int maxcol;
  int newv;
  static int poly[DIM_MAX2] =
  {
        1,    3,    7,   11,   13,   19,   25,   37,   59,   47,
       61,   55,   41,   67,   97,   91,  109,  103,  115,  131,
      193,  137,  145,  143,  241,  157,  185,  167,  229,  171,
      213,  191,  253,  203,  211,  239,  247,  285,  369,  299,
      301,  333,  351,  355,  357,  361,  391,  397,  425,  451,
      463,  487,  501,  529,  539,  545,  557,  563,  601,  607,
      617,  623,  631,  637,  647,  661,  675,  677,  687,  695, 
      701,  719,  721,  731,  757,  761,  787,  789,  799,  803,
      817,  827,  847,  859,  865,  875,  877,  883,  895,  901,
      911,  949,  953,  967,  971,  973,  981,  985,  995, 1001,
     1019, 1033, 1051, 1063, 1069, 1125, 1135, 1153, 1163, 1221,
     1239, 1255, 1267, 1279, 1293, 1305, 1315, 1329, 1341, 1347,
     1367, 1387, 1413, 1423, 1431, 1441, 1479, 1509, 1527, 1531,
     1555, 1557, 1573, 1591, 1603, 1615, 1627, 1657, 1663, 1673, 
     1717, 1729, 1747, 1759, 1789, 1815, 1821, 1825, 1849, 1863,
     1869, 1877, 1881, 1891, 1917, 1933, 1939, 1969, 2011, 2035,
     2041, 2053, 2071, 2091, 2093, 2119, 2147, 2149, 2161, 2171,
     2189, 2197, 2207, 2217, 2225, 2255, 2257, 2273, 2279, 2283,
     2293, 2317, 2323, 2341, 2345, 2363, 2365, 2373, 2377, 2385,
     2395, 2419, 2421, 2431, 2435, 2447, 2475, 2477, 2489, 2503, 
     2521, 2533, 2551, 2561, 2567, 2579, 2581, 2601, 2633, 2657,
     2669, 2681, 2687, 2693, 2705, 2717, 2727, 2731, 2739, 2741,
     2773, 2783, 2793, 2799, 2801, 2811, 2819, 2825, 2833, 2867,
     2879, 2881, 2891, 2905, 2911, 2917, 2927, 2941, 2951, 2955,
     2963, 2965, 2991, 2999, 3005, 3017, 3035, 3037, 3047, 3053,
     3083, 3085, 3097, 3103, 3159, 3169, 3179, 3187, 3205, 3209,
     3223, 3227, 3229, 3251, 3263, 3271, 3277, 3283, 3285, 3299,
     3305, 3319, 3331, 3343, 3357, 3367, 3373, 3393, 3399, 3413,
     3417, 3427, 3439, 3441, 3475, 3487, 3497, 3515, 3517, 3529,
     3543, 3547, 3553, 3559, 3573, 3589, 3613, 3617, 3623, 3627,
     3635, 3641, 3655, 3659, 3669, 3679, 3697, 3707, 3709, 3713,
     3731, 3743, 3747, 3771, 3791, 3805, 3827, 3833, 3851, 3865,
     3889, 3895, 3933, 3947, 3949, 3957, 3971, 3985, 3991, 3995,
     4007, 4013, 4021, 4045, 4051, 4069, 4073, 4179, 4201, 4219,
     4221, 4249, 4305, 4331, 4359, 4383, 4387, 4411, 4431, 4439,
     4449, 4459, 4485, 4531, 4569, 4575, 4621, 4663, 4669, 4711,
     4723, 4735, 4793, 4801, 4811, 4879, 4893, 4897, 4921, 4927,
     4941, 4977, 5017, 5027, 5033, 5127, 5169, 5175, 5199, 5213,
     5223, 5237, 5287, 5293, 5331, 5391, 5405, 5453, 5523, 5573,
     5591, 5597, 5611, 5641, 5703, 5717, 5721, 5797, 5821, 5909,
     5913, 5955, 5957, 6005, 6025, 6061, 6067, 6079, 6081, 6231,
     6237, 6289, 6295, 6329, 6383, 6427, 6453, 6465, 6501, 6523,
     6539, 6577, 6589, 6601, 6607, 6631, 6683, 6699, 6707, 6761,
     6795, 6865, 6881, 6901, 6923, 6931, 6943, 6999, 7057, 7079,
     7103, 7105, 7123, 7173, 7185, 7191, 7207, 7245, 7303, 7327, 
     7333, 7355, 7365, 7369, 7375, 7411, 7431, 7459, 7491, 7505, 
     7515, 7541, 7557, 7561, 7701, 7705, 7727, 7749, 7761, 7783,
     7795, 7823, 7907, 7953, 7963, 7975, 8049, 8089, 8123, 8125,
     8137, 8219, 8231, 8245, 8275, 8293, 8303, 8331, 8333, 8351,
     8357, 8367, 8379, 8381, 8387, 8393, 8417, 8435, 8461, 8469,
     8489, 8495, 8507, 8515, 8551, 8555, 8569, 8585, 8599, 8605,
     8639, 8641, 8647, 8653, 8671, 8675, 8689, 8699, 8729, 8741,
     8759, 8765, 8771, 8795, 8797, 8825, 8831, 8841, 8855, 8859,
     8883, 8895, 8909, 8943, 8951, 8955, 8965, 8999, 9003, 9031,
     9045, 9049, 9071, 9073, 9085, 9095, 9101, 9109, 9123, 9129,
     9137, 9143, 9147, 9185, 9197, 9209, 9227, 9235, 9247, 9253,
     9257, 9277, 9297, 9303, 9313, 9325, 9343, 9347, 9371, 9373,
     9397, 9407, 9409, 9415, 9419, 9443, 9481, 9495, 9501, 9505,
     9517, 9529, 9555, 9557, 9571, 9585, 9591, 9607, 9611, 9621,
     9625, 9631, 9647, 9661, 9669, 9679, 9687, 9707, 9731, 9733,
     9745, 9773, 9791, 9803, 9811, 9817, 9833, 9847, 9851, 9863,
     9875, 9881, 9905, 9911, 9917, 9923, 9963, 9973,10003,10025,
    10043,10063,10071,10077,10091,10099,10105,10115,10129,10145,
    10169,10183,10187,10207,10223,10225,10247,10265,10271,10275,
    10289,10299,10301,10309,10343,10357,10373,10411,10413,10431,
    10445,10453,10463,10467,10473,10491,10505,10511,10513,10523,
    10539,10549,10559,10561,10571,10581,10615,10621,10625,10643,
    10655,10671,10679,10685,10691,10711,10739,10741,10755,10767,
    10781,10785,10803,10805,10829,10857,10863,10865,10875,10877,
    10917,10921,10929,10949,10967,10971,10987,10995,11009,11029,
    11043,11045,11055,11063,11075,11081,11117,11135,11141,11159,
    11163,11181,11187,11225,11237,11261,11279,11297,11307,11309,
    11327,11329,11341,11377,11403,11405,11413,11427,11439,11453,
    11461,11473,11479,11489,11495,11499,11533,11545,11561,11567,
    11575,11579,11589,11611,11623,11637,11657,11663,11687,11691,
    11701,11747,11761,11773,11783,11795,11797,11817,11849,11855,
    11867,11869,11873,11883,11919,11921,11927,11933,11947,11955,
    11961,11999,12027,12029,12037,12041,12049,12055,12095,12097,
    12107,12109,12121,12127,12133,12137,12181,12197,12207,12209,
    12239,12253,12263,12269,12277,12287,12295,12309,12313,12335,
    12361,12367,12391,12409,12415,12433,12449,12469,12479,12481,
    12499,12505,12517,12527,12549,12559,12597,12615,12621,12639,
    12643,12657,12667,12707,12713,12727,12741,12745,12763,12769,
    12779,12781,12787,12799,12809,12815,12829,12839,12857,12875,
    12883,12889,12901,12929,12947,12953,12959,12969,12983,12987,
    12995,13015,13019,13031,13063,13077,13103,13137,13149,13173,
    13207,13211,13227,13241,13249,13255,13269,13283,13285,13303,
    13307,13321,13339,13351,13377,13389,13407,13417,13431,13435,
    13447,13459,13465,13477,13501,13513,13531,13543,13561,13581,
    13599,13605,13617,13623,13637,13647,13661,13677,13683,13695,
    13725,13729,13753,13773,13781,13785,13795,13801,13807,13825,
    13835,13855,13861,13871,13883,13897,13905,13915,13939,13941,
    13969,13979,13981,13997,14027,14035,14037,14051,14063,14085,
    14095,14107,14113,14125,14137,14145,14151,14163,14193,14199,
    14219,14229,14233,14243,14277,14287,14289,14295,14301,14305,
    14323,14339,14341,14359,14365,14375,14387,14411,14425,14441,
    14449,14499,14513,14523,14537,14543,14561,14579,14585,14593,
    14599,14603,14611,14641,14671,14695,14701,14723,14725,14743,
    14753,14759,14765,14795,14797,14803,14831,14839,14845,14855,
    14889,14895,14909,14929,14941,14945,14951,14963,14965,14985,
    15033,15039,15053,15059,15061,15071,15077,15081,15099,15121,
    15147,15149,15157,15167,15187,15193,15203,15205,15215,15217,
    15223,15243,15257,15269,15273,15287,15291,15313,15335,15347,
    15359,15373,15379,15381,15391,15395,15397,15419,15439,15453,
    15469,15491,15503,15517,15527,15531,15545,15559,15593,15611,
    15613,15619,15639,15643,15649,15661,15667,15669,15681,15693,
    15717,15721,15741,15745,15765,15793,15799,15811,15825,15835,
    15847,15851,15865,15877,15881,15887,15899,15915,15935,15937,
    15955,15973,15977,16011,16035,16061,16069,16087,16093,16097,
    16121,16141,16153,16159,16165,16183,16189,16195,16197,16201,
    16209,16215,16225,16259,16265,16273,16299,16309,16355,16375,
    16381 };
  static float recipd;
  static int seed_save = 0;
  int seed_temp;
  static int v[DIM_MAX][LOG_MAX];
//
  if ( !initialized || dim_num != dim_num_save )
  {
    initialized = true;
    for ( i = 0; i < DIM_MAX; i++ )
    {
      for ( j = 0; j < LOG_MAX; j++ )
      {
        v[i][j] = 0;
      }
    }
//
//  Initialize (part of) V.
//
    v[ 0][0] = 1;
    v[ 1][0] = 1;
    v[ 2][0] = 1;
    v[ 3][0] = 1;
    v[ 4][0] = 1;
    v[ 5][0] = 1;
    v[ 6][0] = 1;
    v[ 7][0] = 1;
    v[ 8][0] = 1;
    v[ 9][0] = 1;
    v[10][0] = 1;
    v[11][0] = 1;
    v[12][0] = 1;
    v[13][0] = 1;
    v[14][0] = 1;
    v[15][0] = 1;
    v[16][0] = 1;
    v[17][0] = 1;
    v[18][0] = 1;
    v[19][0] = 1;
    v[20][0] = 1;
    v[21][0] = 1;
    v[22][0] = 1;
    v[23][0] = 1;
    v[24][0] = 1;
    v[25][0] = 1;
    v[26][0] = 1;
    v[27][0] = 1;
    v[28][0] = 1;
    v[29][0] = 1;
    v[30][0] = 1;
    v[31][0] = 1;
    v[32][0] = 1;
    v[33][0] = 1;
    v[34][0] = 1;
    v[35][0] = 1;
    v[36][0] = 1;
    v[37][0] = 1;
    v[38][0] = 1;
    v[39][0] = 1;

    v[ 2][1] = 1;
    v[ 3][1] = 3;
    v[ 4][1] = 1;
    v[ 5][1] = 3;
    v[ 6][1] = 1;
    v[ 7][1] = 3;
    v[ 8][1] = 3;
    v[ 9][1] = 1;
    v[10][1] = 3;
    v[11][1] = 1;
    v[12][1] = 3;
    v[13][1] = 1;
    v[14][1] = 3;
    v[15][1] = 1;
    v[16][1] = 1;
    v[17][1] = 3;
    v[18][1] = 1;
    v[19][1] = 3;
    v[20][1] = 1;
    v[21][1] = 3;
    v[22][1] = 1;
    v[23][1] = 3;
    v[24][1] = 3;
    v[25][1] = 1;
    v[26][1] = 3;
    v[27][1] = 1;
    v[28][1] = 3;
    v[29][1] = 1;
    v[30][1] = 3;
    v[31][1] = 1;
    v[32][1] = 1;
    v[33][1] = 3;
    v[34][1] = 1;
    v[35][1] = 3;
    v[36][1] = 1;
    v[37][1] = 3;
    v[38][1] = 1;
    v[39][1] = 3;

    v[ 3][2] = 7;
    v[ 4][2] = 5;
    v[ 5][2] = 1;
    v[ 6][2] = 3;
    v[ 7][2] = 3;
    v[ 8][2] = 7;
    v[ 9][2] = 5;
    v[10][2] = 5;
    v[11][2] = 7;
    v[12][2] = 7;
    v[13][2] = 1;
    v[14][2] = 3;
    v[15][2] = 3;
    v[16][2] = 7;
    v[17][2] = 5;
    v[18][2] = 1;
    v[19][2] = 1;
    v[20][2] = 5;
    v[21][2] = 3;
    v[22][2] = 3;
    v[23][2] = 1;
    v[24][2] = 7;
    v[25][2] = 5;
    v[26][2] = 1;
    v[27][2] = 3;
    v[28][2] = 3;
    v[29][2] = 7;
    v[30][2] = 5;
    v[31][2] = 1;
    v[32][2] = 1;
    v[33][2] = 5;
    v[34][2] = 7;
    v[35][2] = 7;
    v[36][2] = 5;
    v[37][2] = 1;
    v[38][2] = 3;
    v[39][2] = 3;

    v[ 5][3] =  1;
    v[ 6][3] =  7;
    v[ 7][3] =  9;
    v[ 8][3] = 13;
    v[ 9][3] = 11;
    v[10][3] =  1;
    v[11][3] =  3;
    v[12][3] =  7;
    v[13][3] =  9;
    v[14][3] =  5;
    v[15][3] = 13;
    v[16][3] = 13;
    v[17][3] = 11;
    v[18][3] =  3;
    v[19][3] = 15;
    v[20][3] =  5;
    v[21][3] =  3;
    v[22][3] = 15;
    v[23][3] =  7;
    v[24][3] =  9;
    v[25][3] = 13;
    v[26][3] =  9;
    v[27][3] =  1;
    v[28][3] = 11;
    v[29][3] =  7;
    v[30][3] =  5;
    v[31][3] = 15;
    v[32][3] =  1;
    v[33][3] = 15;
    v[34][3] = 11;
    v[35][3] =  5;
    v[36][3] =  3;
    v[37][3] =  1;
    v[38][3] =  7;
    v[39][3] =  9;
  
    v[ 7][4] =  9;
    v[ 8][4] =  3;
    v[ 9][4] = 27;
    v[10][4] = 15;
    v[11][4] = 29;
    v[12][4] = 21;
    v[13][4] = 23;
    v[14][4] = 19;
    v[15][4] = 11;
    v[16][4] = 25;
    v[17][4] =  7;
    v[18][4] = 13;
    v[19][4] = 17;
    v[20][4] =  1;
    v[21][4] = 25;
    v[22][4] = 29;
    v[23][4] =  3;
    v[24][4] = 31;
    v[25][4] = 11;
    v[26][4] =  5;
    v[27][4] = 23;
    v[28][4] = 27;
    v[29][4] = 19;
    v[30][4] = 21;
    v[31][4] =  5;
    v[32][4] =  1;
    v[33][4] = 17;
    v[34][4] = 13;
    v[35][4] =  7;
    v[36][4] = 15;
    v[37][4] =  9;
    v[38][4] = 31;
    v[39][4] =  9;

    v[13][5] = 37;
    v[14][5] = 33;
    v[15][5] =  7;
    v[16][5] =  5;
    v[17][5] = 11;
    v[18][5] = 39;
    v[19][5] = 63;
    v[20][5] = 27;
    v[21][5] = 17;
    v[22][5] = 15;
    v[23][5] = 23;
    v[24][5] = 29;
    v[25][5] =  3;
    v[26][5] = 21;
    v[27][5] = 13;
    v[28][5] = 31;
    v[29][5] = 25;
    v[30][5] =  9;
    v[31][5] = 49;
    v[32][5] = 33;
    v[33][5] = 19;
    v[34][5] = 29;
    v[35][5] = 11;
    v[36][5] = 19;
    v[37][5] = 27;
    v[38][5] = 15;
    v[39][5] = 25;

    v[19][6] =  13;
    v[20][6] =  35;
    v[21][6] = 115;
    v[22][6] =  41;
    v[23][6] =  79;
    v[24][6] =  17;
    v[25][6] =  29;
    v[26][6] = 119;
    v[27][6] =  75;
    v[28][6] =  73;
    v[29][6] = 105;
    v[30][6] =   7;
    v[31][6] =  59;
    v[32][6] =  65;
    v[33][6] =  21;
    v[34][6] =   3;
    v[35][6] = 113;
    v[36][6] =  61;
    v[37][6] =  89;
    v[38][6] =  45;
    v[39][6] = 107;

    v[37][7] =  7;
    v[38][7] = 23;
    v[39][7] = 39;
//
//  Check parameters.
//
    if ( dim_num < 2 || DIM_MAX < dim_num )
    {
      printf("\n");
      printf("I4_SOBOL - Fatal error!\n");
      printf("  The spatial dimension DIM_NUM should satisfy:\n");
      printf("    2 <= DIM_NUM <= %d\n", DIM_MAX);
      printf("  But this input value is DIM_NUM = %d\n", dim_num);
      exit(1);
    }

    dim_num_save = dim_num;
//
//  Set ATMOST = 2^LOG_MAX - 1.
//
    atmost = 0;
    for ( i = 1; i <= LOG_MAX; i++ )
    {
      atmost = 2 * atmost + 1;
    }
//
//  Find the highest 1 bit in ATMOST (should be LOG_MAX).
//
    maxcol = i4_bit_hi1 ( atmost );
//
//  Initialize row 1 of V.
//
    for ( j = 0; j < maxcol; j++ )
    {
      v[0][j] = 1;
    }
//
//  Initialize the remaining rows of V.
//
    for ( i = 1; i < dim_num; i++ )
    {
//
//  The bit pattern of the integer POLY(I) gives the form
//  of polynomial I.
//
//  Find the degree of polynomial I from binary encoding.
//
      j = poly[i];
      m = 0;

      while ( true )
      {
        j = j / 2;
        if ( j <= 0 )
        {
          break;
        }
        m = m + 1;
      }
//
//  We expand this bit pattern to separate components
//  of the logical array INCLUD.
//
      j = poly[i];
      for ( k = m-1; 0 <= k; k-- )
      {
        j2 = j / 2;
        includ[k] = ( j != ( 2 * j2 ) );
        j = j2;
      }
//
//  Calculate the remaining elements of row I as explained
//  in Bratley and Fox, section 2.
//
//  Some tricky indexing here.  Did I change it correctly?
//
      for ( j = m; j < maxcol; j++ )
      {
        newv = v[i][j-m];
        l = 1;

        for ( k = 0; k < m; k++ )
        {
          l = 2 * l;

          if ( includ[k] )
          {
            newv = ( newv ^ ( l * v[i][j-k-1] ) );
          }
        }
        v[i][j] = newv;
      }
    }
//
//  Multiply columns of V by appropriate power of 2.
//
    l = 1;
    for ( j = maxcol-2; 0 <= j; j-- )
    {
      l = 2 * l;
      for ( i = 0; i < dim_num; i++ )
      {
        v[i][j] = v[i][j] * l;
      }
    }
//
//  RECIPD is 1/(common denominator of the elements in V).
//
    recipd = 1.0E+00 / ( ( float ) ( 2 * l ) );
  }

  if ( *seed < 0 )
  {
    *seed = 0;
  }

  if ( *seed == 0 )
  {
    l = 1;
    for ( i = 0; i < dim_num; i++ )
    {
      lastq[i] = 0;
    }
  }
  else if ( *seed == seed_save + 1 )
  {
    l = i4_bit_lo0 ( *seed );
  }
  else if ( *seed <= seed_save )
  {
    seed_save = 0;
    l = 1;
    for ( i = 0; i < dim_num; i++ )
    {
      lastq[i] = 0;
    }

    for ( seed_temp = seed_save; seed_temp <= (*seed)-1; seed_temp++ )
    {

      l = i4_bit_lo0 ( seed_temp );

      for ( i = 0; i < dim_num; i++ )
      {
        lastq[i] = ( lastq[i] ^ v[i][l-1] );
      }
    }
    l = i4_bit_lo0 ( *seed );
  }
  else if ( seed_save+1 < *seed )
  {
    for ( seed_temp = seed_save+1; seed_temp <= (*seed)-1; seed_temp++ )
    {

      l = i4_bit_lo0 ( seed_temp );

      for ( i = 0; i < dim_num; i++ )
      {
        lastq[i] = ( lastq[i] ^ v[i][l-1] );
      }
    }
    l = i4_bit_lo0 ( *seed );
  }
//
//  Check that the user is not calling too many times!
//
  if ( maxcol < l )
  {
    printf("\n");
    printf("I4_SOBOL - Fatal error!\n");
    printf("  The value of SEED seems to be too large!\n");
    printf("  SEED =   %d\n", *seed);
    printf("  MAXCOL = %d\n", maxcol);
    printf("  L =      %d\n", l);
    exit ( 2 );
  }
//
//  Calculate the new components of QUASI.
//  The caret indicates the bitwise exclusive OR.
//
  for ( i = 0; i < dim_num; i++ )
  {
    quasi[i] = ( ( float ) lastq[i] ) * recipd;

    lastq[i] = ( lastq[i] ^ v[i][l-1] );
  }

  seed_save = *seed;
  *seed = *seed + 1;

  return;
# undef DIM_MAX
# undef DIM_MAX2
# undef LOG_MAX
}
//******************************************************************************

int i4_uniform ( int b, int c, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    I4_UNIFORM returns an integer pseudorandom number.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Modified:
//
//    27 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int B, C, the limits of the interval.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, int I4_UNIFORM, a number between A and B.
//
{
  float d;
  int value;

  if ( b <= c )
  {
    d = ( float ) ( b ) + ( float ) ( 1 + c - b ) * d_uniform_01 ( seed );

    value = ( int ) ( d );

    if ( value < b )
    {
      value = b;
    }
    if ( c < value )
    {
      value = c;
    }
  }
  else
  {
    d = ( float ) ( c ) + ( float ) ( 1 + b - c ) * d_uniform_01 ( seed );

    value = ( int ) ( d );

    if ( value < c )
    {
      value = c;
    }
    if ( b < value )
    {
      value = b;
    }
  }

  return value;
}
//******************************************************************************

unsigned int i4_xor ( unsigned int i, unsigned int j )

//******************************************************************************
//
//  Purpose:
//
//    I4_XOR calculates the exclusive OR of two integers.
//
//  Modified:
//
//    16 February 2005
//
//  Author:
//
//   John Burkardt
//
//  Parameters:
//
//    Input, unsigned int I, J, two values whose exclusive OR is needed.
//
//    Output, unsigned int I4_XOR, the exclusive OR of I and J.
//
{
  unsigned int i2;
  unsigned int j2;
  unsigned int k;
  unsigned int l;

  k = 0;
  l = 1;

  while ( i != 0 || j != 0 )
  {
    i2 = i / 2;
    j2 = j / 2;

    if ( 
      ( ( i == 2 * i2 ) && ( j != 2 * j2 ) ) ||
      ( ( i != 2 * i2 ) && ( j == 2 * j2 ) ) )
    {
      k = k + l;
    }

    i = i2;
    j = j2;
    l = 2 * l;
  }

  return k;
}
//******************************************************************************

int i8_bit_hi1 ( long int n )

//******************************************************************************
//
//  Purpose:
//
//    I8_BIT_HI1 returns the position of the high 1 bit base 2 in an integer.
//
//  Example:
//
//       N    Binary    Hi 1
//    ----    --------  ----
//       0           0     0
//       1           1     1
//       2          10     2
//       3          11     2 
//       4         100     3
//       5         101     3
//       6         110     3
//       7         111     3
//       8        1000     4
//       9        1001     4
//      10        1010     4
//      11        1011     4
//      12        1100     4
//      13        1101     4
//      14        1110     4
//      15        1111     4
//      16       10000     5
//      17       10001     5
//    1023  1111111111    10
//    1024 10000000000    11
//    1025 10000000001    11
//
//  Modified:
//
//    03 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, long int N, the integer to be measured.
//    N should be nonnegative.  If N is nonpositive, I8_BIT_HI1
//    will always be 0.
//
//    Output, int I8_BIT_HI1, the number of bits base 2.
//
{
  int bit;

  bit = 0;

  while ( 0 < n )
  {
    bit = bit + 1;
    n = n / 2;
  }

  return bit;
}
//******************************************************************************

int i8_bit_lo0 ( long int n )

//******************************************************************************
//
//  Purpose:
//
//    I8_BIT_LO0 returns the position of the low 0 bit base 2 in an integer.
//
//  Example:
//
//       N    Binary    Lo 0
//    ----    --------  ----
//       0           0     1
//       1           1     2
//       2          10     1
//       3          11     3 
//       4         100     1
//       5         101     2
//       6         110     1
//       7         111     4
//       8        1000     1
//       9        1001     2
//      10        1010     1
//      11        1011     3
//      12        1100     1
//      13        1101     2
//      14        1110     1
//      15        1111     5
//      16       10000     1
//      17       10001     2
//    1023  1111111111     1
//    1024 10000000000     1
//    1025 10000000001     1
//
//  Modified:
//
//    03 August 2004
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, long int N, the integer to be measured.
//    N should be nonnegative.
//
//    Output, int I8_BIT_LO0, the position of the low 1 bit.
//
{
  int bit;
  long int n2;

  bit = 0;

  while ( true )
  {
    bit = bit + 1;
    n2 = n / 2;

    if ( n == 2 * n2 )
    {
      break;
    }

    n = n2;

  }

  return bit;
}
//******************************************************************************

void i8_sobol ( int dim_num, long int *seed, double quasi[ ] )

//******************************************************************************
//
//  Purpose:
//
//    I8_SOBOL generates a new quasirandom Sobol vector with each call.
//
//  Discussion:
//
//    The routine adapts the ideas of Antonov and Saleev.
//
//  Modified:
//
//    27 November 2005
//
//  Author:
//
//    Bennett Fox
//
//    C++ translation by John Burkardt
//
//  Reference:
//
//    Antonov and Saleev,
//    USSR Computational Mathematics and Mathematical Physics,
//    Volume 19, 1980, pages 252 - 256.
//
//    Paul Bratley and Bennett Fox,
//    Algorithm 659:
//    Implementing Sobol's Quasirandom Sequence Generator,
//    ACM Transactions on Mathematical Software,
//    Volume 14, Number 1, pages 88-100, 1988.
//
//    Bennett Fox,
//    Algorithm 647:
//    Implementation and Relative Efficiency of Quasirandom 
//    Sequence Generators,
//    ACM Transactions on Mathematical Software,
//    Volume 12, Number 4, pages 362-376, 1986.
//
//    Ilya Sobol,
//    USSR Computational Mathematics and Mathematical Physics,
//    Volume 16, pages 236-242, 1977.
//
//    Ilya Sobol and Levitan, 
//    The Production of Points Uniformly Distributed in a Multidimensional 
//    Cube (in Russian),
//    Preprint IPM Akad. Nauk SSSR, 
//    Number 40, Moscow 1976.
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//    DIM_NUM must satisfy 2 <= DIM_NUM <= 40.
//
//    Input/output, long int *SEED, the "seed" for the sequence.
//    This is essentially the index in the sequence of the quasirandom
//    value to be generated.  On output, SEED has been set to the
//    appropriate next value, usually simply SEED+1.
//    If SEED is less than 0 on input, it is treated as though it were 0.
//    An input value of 0 requests the first (0-th) element of the sequence.
//
//    Output, double QUASI[DIM_NUM], the next quasirandom vector.
//
{
# define DIM_MAX 40
# define LOG_MAX 62
//
//  Here, we have commented out the definition of ATMOST, because
//  in some cases, a compiler was complaining that the value of ATMOST could not
//  seem to be properly stored.  We only need ATMOST in order to specify MAXCOL,
//  so as long as we set MAXCOL (below) to what we expect it should be, we
//  may be able to get around this difficulty.
//  JVB, 24 January 2006.
//
//static long int atmost = 4611686018427387903;
//
  static int dim_num_save = 0;
  long int i;
  Cbool includ[8];
  static Cbool initialized = false;
  long int j;
  long int j2;
  long int k;
  long int l = 0;
  static long int lastq[DIM_MAX];
  long int m;
  static long int maxcol;
  long int newv;
  static long int poly[DIM_MAX] =
  {
        1,   3,   7,  11,  13,  19,  25,  37,  59,  47,
       61,  55,  41,  67,  97,  91, 109, 103, 115, 131,
      193, 137, 145, 143, 241, 157, 185, 167, 229, 171,
      213, 191, 253, 203, 211, 239, 247, 285, 369, 299 
  };
  static double recipd;
  static long int seed_save = 0;
  long int seed_temp;
  static long int v[DIM_MAX][LOG_MAX];

  if ( !initialized || dim_num != dim_num_save )
  {
    initialized = true;
//
//  Initialize (part of) V.
//
    v[ 0][0] = 1;
    v[ 1][0] = 1;
    v[ 2][0] = 1;
    v[ 3][0] = 1;
    v[ 4][0] = 1;
    v[ 5][0] = 1;
    v[ 6][0] = 1;
    v[ 7][0] = 1;
    v[ 8][0] = 1;
    v[ 9][0] = 1;
    v[10][0] = 1;
    v[11][0] = 1;
    v[12][0] = 1;
    v[13][0] = 1;
    v[14][0] = 1;
    v[15][0] = 1;
    v[16][0] = 1;
    v[17][0] = 1;
    v[18][0] = 1;
    v[19][0] = 1;
    v[20][0] = 1;
    v[21][0] = 1;
    v[22][0] = 1;
    v[23][0] = 1;
    v[24][0] = 1;
    v[25][0] = 1;
    v[26][0] = 1;
    v[27][0] = 1;
    v[28][0] = 1;
    v[29][0] = 1;
    v[30][0] = 1;
    v[31][0] = 1;
    v[32][0] = 1;
    v[33][0] = 1;
    v[34][0] = 1;
    v[35][0] = 1;
    v[36][0] = 1;
    v[37][0] = 1;
    v[38][0] = 1;
    v[39][0] = 1;

    v[ 2][1] = 1;
    v[ 3][1] = 3;
    v[ 4][1] = 1;
    v[ 5][1] = 3;
    v[ 6][1] = 1;
    v[ 7][1] = 3;
    v[ 8][1] = 3;
    v[ 9][1] = 1;
    v[10][1] = 3;
    v[11][1] = 1;
    v[12][1] = 3;
    v[13][1] = 1;
    v[14][1] = 3;
    v[15][1] = 1;
    v[16][1] = 1;
    v[17][1] = 3;
    v[18][1] = 1;
    v[19][1] = 3;
    v[20][1] = 1;
    v[21][1] = 3;
    v[22][1] = 1;
    v[23][1] = 3;
    v[24][1] = 3;
    v[25][1] = 1;
    v[26][1] = 3;
    v[27][1] = 1;
    v[28][1] = 3;
    v[29][1] = 1;
    v[30][1] = 3;
    v[31][1] = 1;
    v[32][1] = 1;
    v[33][1] = 3;
    v[34][1] = 1;
    v[35][1] = 3;
    v[36][1] = 1;
    v[37][1] = 3;
    v[38][1] = 1;
    v[39][1] = 3;

    v[ 3][2] = 7;
    v[ 4][2] = 5;
    v[ 5][2] = 1;
    v[ 6][2] = 3;
    v[ 7][2] = 3;
    v[ 8][2] = 7;
    v[ 9][2] = 5;
    v[10][2] = 5;
    v[11][2] = 7;
    v[12][2] = 7;
    v[13][2] = 1;
    v[14][2] = 3;
    v[15][2] = 3;
    v[16][2] = 7;
    v[17][2] = 5;
    v[18][2] = 1;
    v[19][2] = 1;
    v[20][2] = 5;
    v[21][2] = 3;
    v[22][2] = 3;
    v[23][2] = 1;
    v[24][2] = 7;
    v[25][2] = 5;
    v[26][2] = 1;
    v[27][2] = 3;
    v[28][2] = 3;
    v[29][2] = 7;
    v[30][2] = 5;
    v[31][2] = 1;
    v[32][2] = 1;
    v[33][2] = 5;
    v[34][2] = 7;
    v[35][2] = 7;
    v[36][2] = 5;
    v[37][2] = 1;
    v[38][2] = 3;
    v[39][2] = 3;

    v[ 5][3] =  1;
    v[ 6][3] =  7;
    v[ 7][3] =  9;
    v[ 8][3] = 13;
    v[ 9][3] = 11;
    v[10][3] =  1;
    v[11][3] =  3;
    v[12][3] =  7;
    v[13][3] =  9;
    v[14][3] =  5;
    v[15][3] = 13;
    v[16][3] = 13;
    v[17][3] = 11;
    v[18][3] =  3;
    v[19][3] = 15;
    v[20][3] =  5;
    v[21][3] =  3;
    v[22][3] = 15;
    v[23][3] =  7;
    v[24][3] =  9;
    v[25][3] = 13;
    v[26][3] =  9;
    v[27][3] =  1;
    v[28][3] = 11;
    v[29][3] =  7;
    v[30][3] =  5;
    v[31][3] = 15;
    v[32][3] =  1;
    v[33][3] = 15;
    v[34][3] = 11;
    v[35][3] =  5;
    v[36][3] =  3;
    v[37][3] =  1;
    v[38][3] =  7;
    v[39][3] =  9;
  
    v[ 7][4] =  9;
    v[ 8][4] =  3;
    v[ 9][4] = 27;
    v[10][4] = 15;
    v[11][4] = 29;
    v[12][4] = 21;
    v[13][4] = 23;
    v[14][4] = 19;
    v[15][4] = 11;
    v[16][4] = 25;
    v[17][4] =  7;
    v[18][4] = 13;
    v[19][4] = 17;
    v[20][4] =  1;
    v[21][4] = 25;
    v[22][4] = 29;
    v[23][4] =  3;
    v[24][4] = 31;
    v[25][4] = 11;
    v[26][4] =  5;
    v[27][4] = 23;
    v[28][4] = 27;
    v[29][4] = 19;
    v[30][4] = 21;
    v[31][4] =  5;
    v[32][4] =  1;
    v[33][4] = 17;
    v[34][4] = 13;
    v[35][4] =  7;
    v[36][4] = 15;
    v[37][4] =  9;
    v[38][4] = 31;
    v[39][4] =  9;

    v[13][5] = 37;
    v[14][5] = 33;
    v[15][5] =  7;
    v[16][5] =  5;
    v[17][5] = 11;
    v[18][5] = 39;
    v[19][5] = 63;
    v[20][5] = 27;
    v[21][5] = 17;
    v[22][5] = 15;
    v[23][5] = 23;
    v[24][5] = 29;
    v[25][5] =  3;
    v[26][5] = 21;
    v[27][5] = 13;
    v[28][5] = 31;
    v[29][5] = 25;
    v[30][5] =  9;
    v[31][5] = 49;
    v[32][5] = 33;
    v[33][5] = 19;
    v[34][5] = 29;
    v[35][5] = 11;
    v[36][5] = 19;
    v[37][5] = 27;
    v[38][5] = 15;
    v[39][5] = 25;

    v[19][6] =  13;
    v[20][6] =  35;
    v[21][6] = 115;
    v[22][6] =  41;
    v[23][6] =  79;
    v[24][6] =  17;
    v[25][6] =  29;
    v[26][6] = 119;
    v[27][6] =  75;
    v[28][6] =  73;
    v[29][6] = 105;
    v[30][6] =   7;
    v[31][6] =  59;
    v[32][6] =  65;
    v[33][6] =  21;
    v[34][6] =   3;
    v[35][6] = 113;
    v[36][6] =  61;
    v[37][6] =  89;
    v[38][6] =  45;
    v[39][6] = 107;

    v[37][7] =  7;
    v[38][7] = 23;
    v[39][7] = 39;
//
//  Check parameters.
//
    if ( dim_num < 2 || DIM_MAX < dim_num )
    {
      printf("\n");
      printf("I8_SOBOL - Fatal error!\n");
      printf("  The spatial dimension DIM_NUM should satisfy:\n");
      printf("    2 <= DIM_NUM <= %d\n", DIM_MAX);
      printf("  But this input value is DIM_NUM = %d\n", dim_num);
      exit ( 1 );
    }

    dim_num_save = dim_num;
//
//  Find the number of bits in ATMOST.
//
//  Here, we have short-circuited the computation of MAXCOL from ATMOST, because
//  in some cases, a compiler was complaining that the value of ATMOST could not
//  seem to be properly stored.  We only need ATMOST in order to specify MAXCOL,
//  so if we know what the answer should be we can try to get it this way!
//  JVB, 24 January 2006.
//
//  maxcol = i8_bit_hi1 ( atmost );
//
    maxcol = 62;
//
//  Initialize row 1 of V.
//
    for ( j = 0; j < maxcol; j++ )
    {
      v[0][j] = 1;
    }
//
//  Initialize the remaining rows of V.
//
    for ( i = 1; i < dim_num; i++ )
    {
//
//  The bit pattern of the integer POLY(I) gives the form
//  of polynomial I.
//
//  Find the degree of polynomial I from binary encoding.
//
      j = poly[i];
      m = 0;

      while ( true )
      {
        j = j / 2;
        if ( j <= 0 )
        {
          break;
        }
        m = m + 1;
      }
//
//  We expand this bit pattern to separate components
//  of the logical array INCLUD.
//
      j = poly[i];
      for ( k = m-1; 0 <= k; k-- )
      {
        j2 = j / 2;
        includ[k] = ( j != ( 2 * j2 ) );
        j = j2;
      }
//
//  Calculate the remaining elements of row I as explained
//  in Bratley and Fox, section 2.
//
//  Some tricky indexing here.  Did I change it correctly?
//
      for ( j = m; j < maxcol; j++ )
      {
        newv = v[i][j-m];
        l = 1;

        for ( k = 0; k < m; k++ )
        {
          l = 2 * l;

          if ( includ[k] )
          {
            newv = ( newv ^ ( l * v[i][j-k-1] ) );
          }

        }

        v[i][j] = newv;

      }

    }
//
//  Multiply columns of V by appropriate power of 2.
//
    l = 1;
    for ( j = maxcol - 2; 0 <= j; j-- )
    {
      l = 2 * l;
      for ( i = 0; i < dim_num; i++ )
      {
        v[i][j] = v[i][j] * l;
      }
    }
//
//  RECIPD is 1/(common denominator of the elements in V).
//
    recipd = 1.0E+00 / ( ( double ) ( 2 * l ) );
  }

  if ( *seed < 0 )
  {
    *seed = 0;
  }

  if ( *seed == 0 )
  {
    l = 1;
    for ( i = 0; i < dim_num; i++ )
    {
      lastq[i] = 0;
    }
  }
  else if ( *seed == seed_save + 1 )
  {
    l = i8_bit_lo0 ( *seed );
  }
  else if ( *seed <= seed_save )
  {
    seed_save = 0;
    l = 1;
    for ( i = 0; i < dim_num; i++ )
    {
      lastq[i] = 0;
    }

    for ( seed_temp = seed_save; seed_temp <= (*seed)-1; seed_temp++ )
    {

      l = i8_bit_lo0 ( seed_temp );

      for ( i = 0; i < dim_num; i++ )
      {
        lastq[i] = ( lastq[i] ^ v[i][l-1] );
      }

    }

    l = i8_bit_lo0 ( *seed );
  }
  else if ( seed_save+1 < *seed )
  {
    for ( seed_temp = seed_save+1; seed_temp <= (*seed)-1; seed_temp++ )
    {

      l = i8_bit_lo0 ( seed_temp );

      for ( i = 0; i < dim_num; i++ )
      {
        lastq[i] = ( lastq[i] ^ v[i][l-1] );
      }
    }
    l = i8_bit_lo0 ( *seed );
  }
//
//  Check that the user is not calling too many times!
//
  if ( maxcol < l )
  {
    printf("\n");
    printf("I8_SOBOL - Fatal error!\n");
    printf("  The value of SEED seems to be too large!\n");
    printf("  SEED =   %ld\n", *seed);
    printf("  MAXCOL = %ld\n", maxcol);
    printf("  L =      %ld\n", l);
    exit ( 2 );
  }
//
//  Calculate the new components of QUASI.
//  The caret indicates the bitwise exclusive OR.
//
  for ( i = 0; i < dim_num; i++ )
  {
    quasi[i] = ( ( double ) lastq[i] ) * recipd;

    lastq[i] = ( lastq[i] ^ v[i][l-1] );
  }

  seed_save = *seed;
  *seed = *seed + 1;

  return;
# undef DIM_MAX
# undef LOG_MAX
}
//******************************************************************************

long int i8_uniform ( long int b, long int c, int *seed )

//******************************************************************************
//
//  Purpose:
//
//    I8_UNIFORM returns an integer pseudorandom number.
//
//  Discussion:
//
//    The pseudorandom number should be uniformly distributed
//    between A and B.
//
//  Modified:
//
//    27 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, long int B, C, the limits of the interval.
//
//    Input/output, int *SEED, the "seed" value, which should NOT be 0.
//    On output, SEED has been updated.
//
//    Output, long int I8_UNIFORM, a number between A and B.
//
{
  double d;
  long int value;

  if ( b <= c )
  {
    d = ( double ) ( b ) + ( double ) ( 1 + c - b ) * d_uniform_01 ( seed );

    value = ( long int ) ( d );

    if ( value < b )
    {
      value = b;
    }
    if ( c < value )
    {
      value = c;
    }
  }
  else
  {
    d = ( double ) ( c ) + ( double ) ( 1 + b - c ) * d_uniform_01 ( seed );

    value = ( long int ) ( d );

    if ( value < c )
    {
      value = c;
    }
    if ( b < value )
    {
      value = b;
    }
  }

  return value;
}
//******************************************************************************

unsigned long int i8_xor ( unsigned long int i, unsigned long int j )

//******************************************************************************
//
//  Purpose:
//
//    I8_XOR calculates the exclusive OR of two integers.
//
//  Modified:
//
//    16 February 2005
//
//  Author:
//
//   John Burkardt
//
//  Parameters:
//
//    Input, unsigned long int I, J, two values whose exclusive OR is needed.
//
//    Output, unsigned long int I8_XOR, the exclusive OR of I and J.
//
{
  unsigned long int i2;
  unsigned long int j2;
  unsigned long int k;
  unsigned long int l;

  k = 0;
  l = 1;

  while ( i != 0 || j != 0 )
  {
    i2 = i / 2;
    j2 = j / 2;

    if ( 
      ( ( i == 2 * i2 ) && ( j != 2 * j2 ) ) ||
      ( ( i != 2 * i2 ) && ( j == 2 * j2 ) ) )
    {
      k = k + l;
    }

    i = i2;
    j = j2;
    l = 2 * l;
  }

  return k;
}

#ifdef __cplusplus
}
#endif
