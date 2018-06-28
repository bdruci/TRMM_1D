#include <iostream>
#include <vector>
#include <stdlib.h>
#include <cmath>

using namespace std;

class Particle
{
   public:
      int slab, numEvents;
      double x, mu;
      bool alive;
};

double fRand(double fMin, double fMax)
{
   double f = (double)rand() / RAND_MAX;
   return fMin + f*(fMax-fMin);
}

int main()
{
   //parameters 
   int numHis = 1000000;
   int numSlabs = 2;
   double zmax = 1.0; //cm
   double slabTh = zmax / numSlabs;
   double Siga = 0.5; //cm^-1
   double Sigs = 0.5; 
   double Sigt = Siga + Sigs;

   //tallies
   vector<double> flux;
   vector<double> time;
   vector<double> rightLeak;
   vector<double> leftLeak;
   vector<double> collisionProb;
   for(int i = 0; i < numSlabs; i++)
   {
      flux.push_back(0);
      time.push_back(0);
      rightLeak.push_back(0);
      leftLeak.push_back(0);
      collisionProb.push_back(0);
   }

   //seed random num generator
   srand(1.0);


   //Run Transport
   for(int j = 0 ; j < numHis ; j++)
   {  
      Particle p;
      p.x = fRand(0.0, zmax);
      p.mu = fRand(-1.0,1.0);
      p.numEvents = 0;
      p.alive = 1;
      p.slab = (int) floor(p.x / slabTh);
     // cout << "HISTORY #" << j << endl;
     // cout << "particle has position " << p.x << " and is in slab " << p.slab << " and has mu = " << p.mu <<  endl;
      vector<double> CHFlux; //current history
      vector<double> CHTime;
      vector<double> CHRL; //right leak
      vector<double> CHLL; //left leak
      vector<double> CHCP; //collision prob
      for(int i = 0; i < numSlabs; i++)
      {
         CHFlux.push_back(0);
         CHTime.push_back(0);
         CHRL.push_back(0);
         CHLL.push_back(0);
         CHCP.push_back(0);
      }
      while(p.alive)
      {
         double d2c = -log(fRand(0.0,1.0))/Sigt;
         double d2s;
         if(p.mu > 0)
            d2s = ((1+p.slab)*slabTh - p.x)/p.mu;
         else
            d2s = (p.slab*slabTh - p.x)/p.mu;
       //  cout << "d2c = " << d2c << " d2s = " << d2s << endl;
 
        if(d2c < 0 || d2s < 0)
         {
            cout << "ERROR: d2c = " << d2c << " d2s = " << d2s << endl;
	    cout << " in history " << j << " where the particle was at position x = " << p.x << endl;
            cout << " in slab " << p.slab << " and mu = " << p.mu << endl;
            exit(1);
         }

         if(d2c > d2s)
         {
         // cout << "hit surface!" << endl;
            //score tally
            CHFlux.at(p.slab) += d2s;
            CHTime.at(p.slab) += d2s;

            p.x = p.x + p.mu*d2s;
            if(p.mu < 0)
            {
               CHLL.at(p.slab) += 1;
               p.numEvents += 1;
               p.x -=  0.00001;
               p.slab -= 1;
            }
            if(p.mu > 0)
            {
               CHRL.at(p.slab) += 1;
               p.numEvents += 1;
               p.x += 0.00001;
               p.slab += 1;
            }
            //did it escape?
            if(p.slab == -1) 
            {
               p.alive = 0;
            }
            if(p.slab == numSlabs)
            {
               p.alive = 0;
            }
         }
         else //collision
         {
         // cout << "Collision!" << endl;
            CHFlux.at(p.slab) += d2c;
            CHTime.at(p.slab) += d2c;
            CHCP.at(p.slab) += 1;
            p.numEvents += 1;
            double r = fRand(0.0,1.0);
            if(r < Siga/Sigt) //absorbed
            {
               p.x += p.mu*d2c;
               p.alive = 0;
            }
            else //scatters
            {
               p.x += p.mu*d2c;
               p.mu = fRand(-1.0,1.0);
            }
         }

      }
      //update global tallies
      for(int i = 0; i < numSlabs; i++)
      {
         flux.at(i)      += CHFlux.at(i);
         time.at(i)      += CHTime.at(i) / p.numEvents;
         rightLeak.at(i) += CHRL.at(i) / p.numEvents;
         leftLeak.at(i)  += CHLL.at(i) / p.numEvents;
         collisionProb.at(i) += CHCP.at(i) / p.numEvents;
      }
      CHFlux.clear();
      CHTime.clear();
      CHRL.clear();
      CHLL.clear();
      CHCP.clear();
   }
   
   //scale by numHis   
   for(int i = 0; i < numSlabs; i++)
   {
      flux.at(i) = flux.at(i) / numHis / slabTh;
      time.at(i) = time.at(i) / numHis;
      rightLeak.at(i) = rightLeak.at(i) / numHis;
      leftLeak.at(i)  = leftLeak.at(i) / numHis;
      collisionProb.at(i) = collisionProb.at(i) / numHis;
   }
   
   //Calculate Parameters
   vector<double> sigtdag, sigLdagL, sigLdagR;
   for(int i = 0 ; i < numSlabs; i++)
   {
      sigtdag.push_back(collisionProb.at(i) / time.at(i));
      sigLdagL.push_back(leftLeak.at(i) / time.at(i));
      sigLdagR.push_back(rightLeak.at(i) / time.at(i));
   }

   //OUTPUT
   cout << "FLUX" << endl;
   cout << "=======" << endl;
   for(int i = 0; i < numSlabs; i++)
   {
      cout << flux.at(i) << endl;
   }

   cout << endl; 
   cout << "sigtdag" << endl;
   cout << "==========" << endl;
   for(int i = 0; i < numSlabs; i++)
   {
      cout << sigtdag.at(i) << endl;
   }

   cout << endl;

   cout << "sigLdagL" << endl;
   cout << "============" << endl;
   for(int i = 0 ; i < numSlabs; i++)
   {
      cout << sigLdagL.at(i) << endl;
   }

   cout << endl;
 
   cout << "sigLdagR" << endl;
   cout << "============" << endl;
   for(int i = 0 ; i < numSlabs; i++)
   {
      cout << sigLdagR.at(i) << endl;
   }
   return 0;
}
