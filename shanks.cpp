// Shank's algorithm for recovering secret key or nonce for El Gamal 
// HW1
// Jessica Sorrell, Adv Crypto with Prof Radziszowski
//
// Command line args should be [ modulus, alpha, beta, infile, outfile ]
// where infile is the file containing the ciphertext pairs and outfile
// is the file to which the plaintext will be written

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <cstring>
#include <errno.h>
#include <vector>
#include <list>
#include <algorithm>

typedef struct {
  int first;
  int second;
} tTuple;

bool compare (tTuple t1, tTuple t2){
  return t1.second < t2.second;
}

int modularExp ( int base, int exponent, int modulus){
  
  int result = 1;
  while ( exponent > 0){

    if (exponent % 2 == 1){
      result = (result * base) % modulus;
    }
    exponent = exponent >> 1;
    base = (base * base) % modulus;
  }
  return result;
}

int inverse ( int a, int modulus){
  
  int b0 = modulus, t, q;
  int x0 = 0, x1 = 1;
  if (modulus == 1) return 1;
  while (a > 1) {
    q = a / modulus;
    t = modulus, modulus = a % modulus, a = t;
    t = x0, x0 = x1 - q * x0, x1 = t;
  }
  if (x1 < 0) x1 += b0;
  return x1;
}

int shanks ( int modulus, int alpha, int beta){
  
  std::list<tTuple> alphaMJ;
  std::list<tTuple> betaAlphaI;
  std::list<tTuple>::iterator it;
  
  int m = (int) ceil( sqrt( (float)modulus ) );   
      
      for (int j = 0; j < m; j++){
	tTuple aMJ;
	aMJ.first = j;
	aMJ.second = modularExp(alpha, m*j, modulus);
	alphaMJ.push_back(aMJ);
      }

      for (int i = 0; i < m; i++){
	tTuple baI;
	baI.first = i;
	baI.second = (beta*modularExp(alpha, modulus - (i+1), modulus)) % modulus;
	betaAlphaI.push_back(baI);

      }

      // sort the lists by the second element
      alphaMJ.sort(compare); 
      betaAlphaI.sort(compare);
      
      // list to hold the intersection
      std::list<tTuple> intersection(2*m);
      
      // find the intersection of the two lists
      it = std::set_intersection( alphaMJ.begin(), alphaMJ.end(), 
				  betaAlphaI.begin(), betaAlphaI.end(), 
				  intersection.begin(), compare); 
      
      int j = (*intersection.begin()).first;
      int shared = (*intersection.begin()).second;
      int i = 0;
      
      for (it = betaAlphaI.begin(); it != betaAlphaI.end(); ++it){
	if ((*it).second == shared){
	  i = (*it).first;
	}
      }
      
      int a = (m*j + i) % (modulus);

      return a;
}

int main(int argc, char* argv[]){
  int modulus = atoi(argv[1]);
  int alpha = atoi(argv[2]);
  int beta = atoi(argv[3]);

  std::ofstream of;
  of.open( argv[5] );
  
  std::ifstream infile( argv[4] );

  std::string line;

  if (infile.is_open()){

    // lists of ciphertext pairs
    std::list<tTuple> ciphertext;
    std::list<tTuple> alphaMJ;
    std::list<tTuple> betaAlphaI;
    std::list<tTuple>::iterator it;
  
    while ( getline (infile, line) ){
      tTuple ct;
      char* cstr = new char [line.length() + 1];
      std::strcpy (cstr, line.c_str());
      ct.first =  atoi (strtok( cstr, " "));
      ct.second = atoi (strtok (NULL, " "));
      ciphertext.push_back( ct );
    }
    
    int k = 0;
    int plaintext = 0;
    std::string alphabet = "abcdefghijklmnopqrstuvwxyz";

    int char1, char2, char3;
    for (it = ciphertext.begin(); it!= ciphertext.end(); ++it){
     
      k = shanks(modulus, alpha, (*it).first);
      // of << k << " \\\\";
      plaintext = ((*it).second*
		   inverse(modularExp(beta, k, modulus), modulus)) % modulus ;
    
      char1 = (int) floor(plaintext/(pow(26, 2)));
      char2 = (int) floor((plaintext - char1*pow(26, 2)) / 26);
      char3 = plaintext - (char1*pow(26, 2)) - (char2*26);
      std::cout << char1 << ", " << char2 << ", " << char3 << std::endl;
      of << alphabet[char1] << alphabet[char2] << alphabet[char3] << "\t";
      }

  }
  else {
    std::cout << "Unable to open file";
  }
  
  return 1;
}
