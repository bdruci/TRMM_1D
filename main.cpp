#include <iostream>
#include <vector>
#include <stdlib.h>
#include <cmath>

using namespace std;

class Particle
{
   public:
      int slab;
      double x, mu;
      bool alive;
};

double fRand(double fMin, double fMax)
//PURPOSE: generate a random double between fMin and fMax
//INPUT: fMin: the minumum value you wish to set the range of random numbers to
//       fmax: the maximum ''
//RETURNS: a random number between these two value
//MODIFIES: N/A
{
   double f = (double)rand() / RAND_MAX;
   return fMin + f*(fMax-fMin);
}

int main()
{
   //parameters 
   int numHis = 100000; //number of histories per slab
   int numSlabs = 2; 
   double slabTh = 1.0; //how thick each slab is in cm
   double zmax = slabTh * numSlabs; //total thickness of entire object
   double Siga = 0.5; //cm^-1 cross section of each slab TODO: make these into a vector 
   double Sigs = 0.5; 
   double Sigt = Siga + Sigs;

   //tallies TODO: could these be packaged into a "slab" class? Then I could just have a vector of slabs
   vector<double> flux; //vector containing the flux of each slab -> can this actually be done slab by slab?
   vector<double> time; //the total time each particle takes to meet its end (1 x numSlabs) 
   vector<double> rightLeakProb; //total leakage probability from slab to the right (1 x numSlabs)
   vector<double> leftLeakProb;  //''  to the left ''
   vector<double> collisionProb; //total collision probability in slab. A collision here means a collision that would end the particle track
   vector<double> inscatterProb; //if a collision occurs, what is the likelihood it was an inscatter
   for(int i = 0; i < numSlabs; i++) //fill each vector with zeros TODO: find a more effective way to do this
   {
      flux.push_back(0);
      time.push_back(0);
      rightLeakProb.push_back(0);
      leftLeakProb.push_back(0);
      collisionProb.push_back(0);
      inscatterProb.push_back(0);
   }

   //seed random num generator
   srand(1.0);

   
   //Run Transport in a slab
   for(int k = 0; k < numSlabs ; k++)
   {
     for(int j = 0 ; j < numHis ; j++)
     {  
         Particle p;
         p.x = fRand(slabTh*k, slabTh*(k+1));
         p.mu = fRand(-1.0,1.0);
         p.alive = 1;
         p.slab = k;
         // cout << "HISTORY #" << j << endl;
         // cout << "particle has position " << p.x << " and is in slab " << p.slab << " and has mu = " << p.mu <<  endl;
        // vector<double> CHFlux; //current history
        // vector<double> CHTime;
        // vector<double> CHRL; //right leak
        // vector<double> CHLL; //left leak
        // vector<double> CHCP; //collision prob
        // for(int i = 0; i < numSlabs; i++)
        // {
         //   CHFlux.push_back(0);
           // CHTime.push_back(0);
           // CHRL.push_back(0);
           // CHLL.push_back(0);
           // CHCP.push_back(0);
        // }
         while(p.alive)
         {
            double dist2Coll = -log(fRand(0.0,1.0))/Sigt;
            double dist2Surf;
            if(p.mu > 0)
               dist2Surf = ((1+p.slab)*slabTh - p.x)/p.mu; 
            else
               dist2Surf = (p.slab*slabTh - p.x)/p.mu;
          //  cout << "d2c = " << d2c << " d2s = " << d2s << endl;
 
           if(dist2Coll < 0 || dist2Surf < 0) //TODO: add throw error so we know way more information. 
            {
               cout << "ERROR: d2c = " << dist2Coll << " d2s = " << dist2Surf << endl;
               cout << " in history " << j << " where the particle was at position x = " << p.x << endl;
               cout << " in slab " << p.slab << " and mu = " << p.mu << endl;
               exit(1);
            }

            if(dist2Coll > dist2Surf) //particle leaked from the current slab! 
            {
            // cout << "hit surface!" << endl;
               //score tally 
               flux.at(k) += dist2Surf;
               time.at(k) += dist2Surf;

               p.x = p.x + p.mu*dist2Surf;
               if(p.mu < 0) //leaked to the left
               {
                  leftLeakProb.at(k) += 1;
                  p.x -=  0.00001; //if I wish to adapt this code to do entire slab transport
                  p.slab -= 1; //""
                  p.alive = 0; //TODO: put this as a function in particle
               }
               if(p.mu > 0) //leaked to the right
               {
                  rightLeakProb.at(k) += 1;
                  p.x += 0.00001;
                  p.slab += 1;
                  p.alive = 0; 
               }
            }
            else //collision!
            {
            // cout << "Collision!" << endl;
               flux.at(k) += dist2Coll;
               time.at(k) += dist2Coll; 
               
               //which event occurred? 
               double r = fRand(0.0,1.0);
               if(r < Siga/Sigt) //absorbed
               {
                  p.x += p.mu*dist2Coll;
                  collisionProb.at(k) += 1; //since the collision ended the particle track.
                  p.alive = 0;
               }
               else //scatters
               {
                  p.x += p.mu*dist2Coll;
                  collisionProb.at(k) += 1;
                  inscatterProb.at(k) += 1;
                  p.mu = fRand(-1.0,1.0);
                 // p.alive = 0; //normally we wouldn't kill it, but for this calculation
               }
            }

         }
        /*
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
        */ //remnants of a bygone era
    }
   }
   //scale by numHis TODO: find a better way to do this. vector of slabs and auto would work better  
   for(int i = 0; i < numSlabs; i++)
   {
      flux.at(i) = flux.at(i) / numHis / slabTh;
      time.at(i) = time.at(i) / numHis;
      rightLeakProb.at(i) = rightLeakProb.at(i) / numHis;
      leftLeakProb.at(i)  = leftLeakProb.at(i) / numHis;
      collisionProb.at(i) = collisionProb.at(i) / numHis;
      inscatterProb.at(i) = inscatterProb.at(i) / numHis;
   }
   
   //Calculate Parameters
   vector<double> sigtdag, sigLdagL, sigLdagR, sigSS; //sigma_t, sig leakage to the right, " to the left, sigma self scatter
   for(int i = 0 ; i < numSlabs; i++)
   {
      sigtdag.push_back(collisionProb.at(i) / time.at(i));
      sigLdagL.push_back(leftLeakProb.at(i) / time.at(i));
      sigLdagR.push_back(rightLeakProb.at(i) / time.at(i));
      sigSS.push_back(inscatterProb.at(i) / time.at(i));
   }

   //OUTPUT TODO: again, make a class called "slab"
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

   cout << "sigSS" << endl;
   cout << "=============" << endl;
   for(int i = 0 ; i < numSlabs; i++)
   {
      cout << sigSS.at(i) << endl;
   }
  
   cout << endl;
   return 0;
}
