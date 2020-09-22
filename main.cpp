/*
Taylor Ropiequet
2349522
ropiequet@chapman.edu
CPSC 350-02
Assignment 1

This is my main method file that implements my methods from my .cpp file.
*/


#include <iostream>
#include <string>
#include <fstream>
#include "assign1.h"

using namespace std;

int main(int argc, char** argv){

  dnaProcessor dna;

  remove("taylor.out"); //makes sure output file is clean each time program is executed

  string fileName = argv[1];

  string answer;

  ofstream outfile;
  outfile.open("taylor.out", ios::out | ios::app);

  outfile << "Taylor Ropiequet" << endl;
  outfile << "ID: 2349522" << endl;
  outfile << "CPSC 350-02" << endl;
  outfile << "Assignment 1" << endl;
  outfile << endl;


  while (true){
    if (dna.checkFile(fileName)){
      dna.writeResults(fileName);
      dna.generateStrings(fileName);
      cout << "Would you like to process another file? (Y/N)" << endl;
      cin >> answer;
      while (true){
        if (answer == "y" || answer == "Y"){
          cout << "Please enter file name" << endl;
          cin >> fileName;
          break;
        }else if (answer == "n" || answer == "N"){
          exit(0);
        }else{
          cout << "Invalid input. Please try again! (Y/N)" << endl;
          cin >> answer;
        }
      }
    }else{
      cout << "Error! File does not exist! Please enter a new file name" << endl;
      cin >> fileName;
    }
  }

  return 0;
}
