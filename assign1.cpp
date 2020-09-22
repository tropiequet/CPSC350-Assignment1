/*
Taylor Ropiequet
2349522
ropiequet@chapman.edu
CPSC 350-02
Assignment 1

This is my .cpp implementation file with all my methods to calculate statistics
regarding the DNA strings in the given file. It also generates 1000 new DNA strings
following the mean, variance, and frequency of what was calculated from the given file.
*/


#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h>
#include <ctime>
#include<cstdlib>
#include <string>
#include <random>
#include "assign1.h"
using namespace std;

dnaProcessor::dnaProcessor(){
}

bool dnaProcessor::checkLine(string line){ //makes sure the line is a proper dna sequence
  bool check = true;
  for (int i = 0; i < line.size() - 1; ++i){ //iterates through each char of the line
    if (tolower(line[i]) == 'a'){
      check = true;
    }else if (tolower(line[i]) == 'c'){
      check = true;
    }else if (tolower(line[i]) == 't'){
      check = true;
    }else if (tolower(line[i]) == 'g'){
      check = true;
    }else{
      check = false;
      break;
    }
  }
  return check;
}


double dnaProcessor::calculateSum(string fileName){ //counts the number of nucleotides in the valid DNA strings
  ifstream infile;
  infile.open(fileName);
  string line;
  double count = 0;
  while (infile.good()){
     getline(infile, line);
     if (checkLine(line)){
       count = count + (line.size() - 1); //the -1 gets rid of the new line
     }
  }
  infile.close();
  return count;
}


double dnaProcessor::countStrings(string fileName){ //counts the number of valid DNA strings
  ifstream infile;
  infile.open(fileName);
  string line;
  double count = 0;
  while (infile.good()){
    getline(infile, line);
    if (checkLine(line)){
      count++;
    }
  }
  infile.close();
  return count;
}


double dnaProcessor::calculateMean(string fileName){
  ifstream infile;
  infile.open(fileName);
  double mean = calculateSum(fileName) / countStrings(fileName);
  return mean;
}


double dnaProcessor::calculateVariance(string fileName){
  ifstream infile;
  infile.open(fileName);
  string line;
  double variance = 0;
  double length = 0;
  while (infile.good()){
    getline(infile, line);
    if (checkLine(line)){
      length = length + pow((line.size() - 1) - calculateMean(fileName), 2);
    }
  }

  variance = length / countStrings(fileName);
  infile.close();
  return variance;
}


double dnaProcessor::calculateDeviation(string fileName){
  ifstream infile;
  infile.open(fileName);
  double variance = calculateVariance(fileName);
  double deviation = sqrt(variance);
  infile.close();
  return deviation;
}


double dnaProcessor::calculateNucleotideProbability(string fileName, char c){
  ifstream infile;
  infile.open(fileName);
  double count = 0;
  double prob = 0;
  string line;
  switch (c){
    case 'a':
      while (infile.good()){
        getline(infile, line);
        if (checkLine(line)){
          for (int i = 0; i < line.size() - 1; ++i){
            if (tolower(line[i]) == 'a'){
              count++;
            }
          }

        }
      }
      break;

    case 'c':
      while (infile.good()){
        getline(infile, line);
        if (checkLine(line)){
          for (int i = 0; i < line.size() - 1; ++i){
            if (tolower(line[i]) == 'c'){
              count++;
            }
          }

        }
      }
      break;

    case 't':
      while (infile.good()){
        getline(infile, line);
        if (checkLine(line)){
          for (int i = 0; i < line.size() - 1; ++i){
            if (tolower(line[i]) == 't'){
              count++;
            }
          }

        }
      }
      break;

    case 'g':
      while (infile.good()){
        getline(infile, line);
        if (checkLine(line)){
          for (int i = 0; i < line.size() - 1; ++i){
            if (tolower(line[i]) == 'g'){
              count++;
            }
          }

        }
      }
      break;
  }
  prob = count / calculateSum(fileName);

  infile.close();

  return prob;

}


double dnaProcessor::calculateBigramProbability(string fileName, string bigram){
  ifstream infile;
  infile.open(fileName);
  double count = 0;
  double prob = 0;
  string line;
  if (bigram == "aa"){
    while (infile.good()){
      getline(infile, line);
      if (checkLine(line)){
        for (int i = 0; i < line.size() - 1; i+=2){
          if (tolower(line[i]) =='a' && tolower(line[i+1]) == 'a'){
            count++;
          }
        }
      }
    }
  }else if (bigram == "ac"){
    while (infile.good()){
      getline(infile, line);
      if (checkLine(line)){
        for (int i = 0; i < line.size() - 1; i+=2){
          if (tolower(line[i]) =='a' && tolower(line[i+1]) == 'c'){
            count++;
          }
        }
      }
    }
  }else if (bigram == "at"){
    while (infile.good()){
      getline(infile, line);
      if (checkLine(line)){
        for (int i = 0; i < line.size() - 1; i+=2){
          if (tolower(line[i]) =='a' && tolower(line[i+1]) == 't'){
            count++;
          }
        }
      }
    }
  }else if (bigram == "ag"){
    while (infile.good()){
      getline(infile, line);
      if (checkLine(line)){
        for (int i = 0; i < line.size() - 1; i+=2){
          if (tolower(line[i]) =='a' && tolower(line[i+1]) == 'g'){
            count++;
          }
        }
      }
    }
  }else if (bigram == "ca"){
    while (infile.good()){
      getline(infile, line);
      if (checkLine(line)){
        for (int i = 0; i < line.size() - 1; i+=2){
          if (tolower(line[i]) =='c' && tolower(line[i+1]) == 'a'){
            count++;
          }
        }
      }
    }
  }else if (bigram == "cc"){
    while (infile.good()){
      getline(infile, line);
      if (checkLine(line)){
        for (int i = 0; i < line.size() - 1; i+=2){
          if (tolower(line[i]) =='c' && tolower(line[i+1]) == 'c'){
            count++;
          }
        }
      }
    }
  }else if (bigram == "ct"){
    while (infile.good()){
      getline(infile, line);
      if (checkLine(line)){
        for (int i = 0; i < line.size() - 1; i+=2){
          if (tolower(line[i]) =='c' && tolower(line[i+1]) == 't'){
            count++;
          }
        }
      }
    }
  }else if (bigram == "cg"){
    while (infile.good()){
      getline(infile, line);
      if (checkLine(line)){
        for (int i = 0; i < line.size() - 1; i+=2){
          if (tolower(line[i]) =='c' && tolower(line[i+1]) == 'g'){
            count++;
          }
        }
      }
    }
  }else if (bigram == "ta"){
    while (infile.good()){
      getline(infile, line);
      if (checkLine(line)){
        for (int i = 0; i < line.size() - 1; i+=2){
          if (tolower(line[i]) =='t' && tolower(line[i+1]) == 'a'){
            count++;
          }
        }
      }
    }
  }else if (bigram == "tc"){
    while (infile.good()){
      getline(infile, line);
      if (checkLine(line)){
        for (int i = 0; i < line.size() - 1; i+=2){
          if (tolower(line[i]) =='t' && tolower(line[i+1]) == 'c'){
            count++;
          }
        }
      }
    }
  }else if (bigram == "tt"){
    while (infile.good()){
      getline(infile, line);
      if (checkLine(line)){
        for (int i = 0; i < line.size() - 1; i+=2){
          if (tolower(line[i]) =='t' && tolower(line[i+1]) == 't'){
            count++;
          }
        }
      }
    }
  }else if (bigram == "tg"){
    while (infile.good()){
      getline(infile, line);
      if (checkLine(line)){
        for (int i = 0; i < line.size() - 1; i+=2){
          if (tolower(line[i]) =='t' && tolower(line[i+1]) == 'g'){
            count++;
          }
        }
      }
    }
  }else if (bigram == "ga"){
    while (infile.good()){
      getline(infile, line);
      if (checkLine(line)){
        for (int i = 0; i < line.size() - 1; i+=2){
          if (tolower(line[i]) =='g' && tolower(line[i+1]) == 'a'){
            count++;
          }
        }
      }
    }
  }else if (bigram == "gc"){
    while (infile.good()){
      getline(infile, line);
      if (checkLine(line)){
        for (int i = 0; i < line.size() - 1; i+=2){
          if (tolower(line[i]) =='g' && tolower(line[i+1]) == 'c'){
            count++;
          }
        }
      }
    }
  }else if (bigram == "gt"){
    while (infile.good()){
      getline(infile, line);
      if (checkLine(line)){
        for (int i = 0; i < line.size() - 1; i+=2){
          if (tolower(line[i]) =='g' && tolower(line[i+1]) == 't'){
            count++;
          }
        }
      }
    }
  }else if (bigram == "gg"){
    while (infile.good()){
      getline(infile, line);
      if (checkLine(line)){
        for (int i = 0; i < line.size() - 1; i+=2){
          if (tolower(line[i]) =='g' && tolower(line[i+1]) == 'g'){
            count++;
          }
        }
      }
    }
  }
  infile.close();

 return (count / (calculateSum(fileName) / 2.00000));
}


void dnaProcessor::writeResults(string fileName){
  ifstream infile;
  infile.open(fileName);

  ofstream outfile;
  outfile.open("taylor.out", ios::out | ios::app);


  outfile << "The Sum of the length of the DNA strings is: " << calculateSum(fileName) << endl;
  outfile << "The Mean of the length of the DNA strings is: " << calculateMean(fileName) << endl;
  outfile << "The Variance of the length of the DNA strings is: " << calculateVariance(fileName) << endl;
  outfile << "The Standard Deviation of the length of the DNA strings is: " << calculateDeviation(fileName) << endl;

  outfile << endl;

  outfile << "Here is the relative probability of each nucleotide:" << endl;
  outfile << "A: " << calculateNucleotideProbability(fileName,'a') << endl;
  outfile << "C: " << calculateNucleotideProbability(fileName,'c') << endl;
  outfile << "T: " << calculateNucleotideProbability(fileName,'t') << endl;
  outfile << "G: " << calculateNucleotideProbability(fileName,'g') << endl;

  outfile << endl;
  outfile << "Here is the relative probability of each bigram:" << endl;
  outfile << "AA: " << calculateBigramProbability(fileName,"aa") << endl;
  outfile << "AC: " << calculateBigramProbability(fileName,"ac") << endl;
  outfile << "AT: " << calculateBigramProbability(fileName,"at") << endl;
  outfile << "AG: " << calculateBigramProbability(fileName,"ag") << endl;
  outfile << "CA: " << calculateBigramProbability(fileName,"ca") << endl;
  outfile << "CC: " << calculateBigramProbability(fileName,"cc") << endl;
  outfile << "CT: " << calculateBigramProbability(fileName,"ct") << endl;
  outfile << "CG: " << calculateBigramProbability(fileName,"cg") << endl;
  outfile << "TA: " << calculateBigramProbability(fileName,"ta") << endl;
  outfile << "TC: " << calculateBigramProbability(fileName,"tc") << endl;
  outfile << "TT: " << calculateBigramProbability(fileName,"tt") << endl;
  outfile << "TG: " << calculateBigramProbability(fileName,"tg") << endl;
  outfile << "GA: " << calculateBigramProbability(fileName,"ga") << endl;
  outfile << "GC: " << calculateBigramProbability(fileName,"gc") << endl;
  outfile << "GT: " << calculateBigramProbability(fileName,"gt") << endl;
  outfile << "GG: " << calculateBigramProbability(fileName,"gg") << endl;

  outfile << endl;

  infile.close();
  outfile.close();
}



string dnaProcessor::getNucleotides(string fileName){ //adds all nucleotides to a string
  ifstream infile;
  infile.open(fileName);
  string dna;
  string line;
  while (infile.good()){
    getline(infile, line);
    if (checkLine(line)){
      line.erase(line.size() - 1); //gets rid off newline character
      dna.append(line);
    }
  }
  infile.close();
  return dna;
}



void dnaProcessor::generateStrings(string fileName){
  ifstream infile;
  infile.open(fileName);
  ofstream outfile;
  outfile.open("taylor.out", ios::out | ios::app);

  srand (time(NULL));


  for (int i = 0; i < 1000; ++i){
    double a = (rand() + 1.0) / (RAND_MAX + 1.0);
    double b = (rand() + 1.0) / (RAND_MAX + 1.0);
    double c = (sqrt(-2.0 * log(a))) * (cos(2.0 * M_PI * b));
    double d = round((calculateDeviation(fileName) * c) + calculateMean(fileName));


    string dna;
    string nucleotides = getNucleotides(fileName);
    int max = calculateSum(fileName);


    for (int j = 0; j < d; ++j){
      int random = rand() % max;
      if (tolower(nucleotides.at(random)) == 'a'){
        dna.append("A");
      }else if (tolower(nucleotides.at(random)) == 'c'){
        dna.append("C");
      }else if (tolower(nucleotides.at(random)) == 't'){
        dna.append("T");
      }else if (tolower(nucleotides.at(random)) == 'g'){
        dna.append("G");
      }
       // nucleotides.erase(random, 1);
       // max--;

    }

    outfile << dna << endl;
  }

  outfile << "--------------------------------------------------------------------" << endl;
  infile.close();
  outfile.close();
}



bool dnaProcessor::checkFile(string fileName){
  ifstream infile(fileName);
  string line;
  getline(infile, line);
  if (line.size() == 0){
    return false;
  }else{
    return true;
  }
}
