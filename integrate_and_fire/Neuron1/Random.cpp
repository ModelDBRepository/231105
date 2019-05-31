#include "Random.h"

Random::Random()
{
	count = 0;
	//unsigned seed = time(NULL);
	//generator.seed(seed);

   //rand=Ran(time(NULL));
	rand = Ran(175737);
   normaldist = Normaldev_BM(0, 1, 2546247);
}
Random::Random(int seed)
{
	//unsigned seed = time(NULL);
	//generator.seed(seed);
   rand=Ran(seed);
   normaldist = Normaldev_BM(0, 1, seed+431112);
}

/*void Random::resetRandom(int seed) {
	//rand.resetRan(seed);
	rand = Ran(seed);
	cout << "resetRandom" << endl;
	//rand = Ran(seed);
}*/


Doub Random::RandomUniform0_to_1() {
	//Doub x_rand = distribution(generator);
   //Doub x_rand = (rand()%1000+0.5)/1000;
	count++;
   Doub x_rand = rand.doub();
	return x_rand;
}

Random::~Random()
{
}
