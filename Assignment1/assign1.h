/*
Taylor Ropiequet
2349522
ropiequet@chapman.edu
CPSC 350-02
Assignment 1

This is my .h header file with a class "dnaProcessor" listing my variables and methods.
*/

#include <iostream>

using namespace std;

class dnaProcessor{

public:
  //list variables/constructors
  dnaProcessor(); //default constructor
  //void openFile(string fileName);
  bool checkLine(string line);
  double calculateSum(string fileName);
  double countStrings(string fileName);
  double calculateMean(string fileName);
  double calculateVariance(string fileName);
  double calculateDeviation(string fileName);
  double calculateNucleotideProbability(string fileName, char c);
  double calculateBigramProbability(string fileName, string bigram);
  void writeResults(string fileName);
  void generateStrings(string fileName);
  string getNucleotides(string fileName);
  bool checkFile(string fileName);

private:
  //list variables/constructors
};
