#define D 200
#define M 100
#define SCALE 100000000

typedef struct protocolIO
{
  long long beta1[M][D], beta2[M][D], beta_avg[D];
}protocolIO;

void aggregate(void* args);
