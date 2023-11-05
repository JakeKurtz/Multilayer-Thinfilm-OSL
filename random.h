#include "stdosl.h"

void rng_seed(output int rng, int seed)
{
  int chash = seed;
  if (chash == 0) chash = 1;
  rng = chash * 891694213;
}

float rng_uniform(output int rng)
{
  float res = rng / float(2137483647) * 0.5 + 0.5;
  rng *= 891694213;
  return res;
}

float rand_range(float start, float end, vector p)
{
    return (end-start)*hashnoise(p) + start;
}