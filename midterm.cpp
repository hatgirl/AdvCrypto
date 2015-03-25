// Code for the take-home Advanced Cryptography midterm with Professor
// Radziszowski at R.I.T.
// Most of these questions pertain to elliptic curves
//
// author: Jessica Sorrell
// date: 24/3/15

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
#include <limits>



typedef struct{
  int x;
  int y;
} pair;



int invSub (int a, int modulus){
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



// modular inverse
int inverse ( int a, int modulus){
    
 
  while ( a >= modulus ){
    a -= modulus;
  }
  while (a < 0){
    a += modulus;
  }
  
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

// Addition as defined for elliptic curve groups
// Takes the two points to be added, as well as the coefficient of the 
// first degree term. This coefficient is used to calculate lambda in the 
// case of point doubling.
pair ecAdd( pair a, pair b, int coeff, int modulus ){

  pair sum;
  int lam = 0;
  if (a.x == std::numeric_limits<int>::max()){
    return b;
  }
  else if (b.x == std::numeric_limits<int>::max()){
    return a;
  }
  else if ((a.y+modulus)%modulus == (-b.y + modulus)%modulus){
    sum.x = std::numeric_limits<int>::max();
    sum.y = std::numeric_limits<int>::max();
    return sum;
  }
  else if ( ((a.x+modulus)%modulus == (b.x+modulus)%modulus) 
	    && ( (a.y+modulus)%modulus == (b.y+modulus)%modulus)){
    lam = ( int(pow(a.x, 2)) *3 + coeff) * inverse(2*a.y, modulus) % modulus;
    }
  else {
    lam = ((b.y - a.y) * inverse(b.x - a.x, modulus)) % modulus;
    }

  sum.x = (lam*lam - a.x - b.x)  % modulus;
  sum.y = (lam*(a.x - sum.x) - a.y) % modulus;
  return sum;
}


// elliptic curve multiplication
// double and add
pair ecMult( pair a, int k, int coeff, int modulus){

  pair product;

  if (k == 0){
    product.x = std::numeric_limits<int>::max();
    product.y = std::numeric_limits<int>::max();
    return product;
  }
  else if (k == 1){
    return a;
  }
  else if (k % 2 == 1){
    return ecAdd(a, ecMult(a, k-1, coeff, modulus), coeff, modulus);
  }
  return ecMult( ecAdd(a, a, coeff, modulus), k/2, coeff, modulus );
}


// Check to ensure that the modulus is a prime
// code from Wikipedia's article on primality testing
bool isPrime (int modulus){
    if (modulus <= 3) {
        return modulus > 1;
    } else if (modulus % 2 == 0 || modulus % 3 == 0) {
        return false;
    } else {
        for (unsigned short i = 5; i * i <= modulus; i += 6) {
            if (modulus % i == 0 || modulus % (i + 2) == 0) {
                return false;
            }
        }
        return true;
    }
}

// Find the order of elements in the group.
void elmtOrder ( int curve[4], int modulus ){

  // First thing we need to do is find a generator of the group
  // Z_modulus. 
  // We'll need this to determine quadratic reduosity of the y^2
  // values we calculate, thus determining whether or not a 
  // point is a member of our elliptic curve group.
  if ( !isPrime(modulus) ){
    printf("This given modulus is not prime. Exiting. \n");
    exit(1);
  }  
  int g = 1;
  int order = 1; 
  while (order != modulus-1){
    g += 1;
    int temp = g;
    order = 1;
    while ( temp != 1 ){
      temp = temp*g % modulus;
      order ++;
    }
  }
  
  int intElmts[modulus-1];
  intElmts[0] = 1;
  intElmts[1] = g;

  int temp = g;
  for (int i = 2; i < modulus -1; i++){
    temp = temp*g % modulus;
    intElmts[i] = temp;
    printf("g^%d is %d \n", i, intElmts[i]);
  }

  // list of elements in the group and the list's iterator  
  std::list<pair> elmts;
  std::list<pair>::iterator it;

  // the coefficient of the x term in the curve
  int coeff = curve[1];

  // for temp storing y-coordinates
  int y2 = 0;
  int y = 0;
  for (int i = 0; i < modulus; i++){
    
    pair elmt;
    elmt.x = i;
    // calculate the y^2 value
    y2 = int(pow(i, 3)*curve[3] + pow(i, 2)*curve[2] + pow(i, 1)*curve[1] +
	     curve[0])  % modulus;
      
    // if y^2 is 0, don't bother looking it up, just add it
    if (y2 == 0 ){
      elmt.y = 0;
      elmts.push_back(elmt);
    }
    else{
      // determine whether y2 is a valid square in Z_modulus
      int index = std::distance(intElmts, 
				std::find(intElmts, intElmts + modulus, y2));
      // we've stored elements in the intElmts array such that even indices 
      // correspond to quadratic residues.
      if ( index % 2 == 0){
	pair elmt;
	elmt.x = i;
	//the y coordinate will be the square root of y2
	elmt.y = intElmts[index/2];
	elmts.push_back(elmt);
	
	elmt.y = -elmt.y;
	elmts.push_back(elmt);
      }
    }
  }
  printf("Actual number of elements is %d. \n",
  	 (int)elmts.size());
  
  // print out the elements of the group
  for (it = elmts.begin(); it!= elmts.end(); ++it){
   printf("(%d, %d) \n", (*it).x, (*it).y);
  }

  // find the order of all elements and print them out
  it = elmts.begin();
  for (it = elmts.begin(); it!= elmts.end(); ++it){
 
    pair elmt = *it;
    pair temp = elmt;
    int order = 1;
    int counter = 0;
    while (temp.x != std::numeric_limits<int>::max() ){
      temp = ecAdd(temp, elmt, coeff, modulus);
      order ++;
      //   printf("elmt.x is %d \n", elmt.x);
      //printf("elmt.y is %d \n", elmt.y);
      //printf("temp.x is %d \n", temp.x);
      //printf("temp.y is %d \n", temp.y);
   
      counter ++;
    }
    printf("The order of (%d, %d) is %d. \n", elmt.x, elmt.y, order);
  }

}




int main (int argc, char* argv[]){

  int curve[] = {28, 1, 0, 1};
  int modulus = 71;
 
  elmtOrder(curve, modulus);
  return 1;
}
