#pragma once

// system headers
#include <vector>

// Chromosom
#include "Chromosom.h"
using namespace std;

class GeneticAlgorithm {
public:
    GeneticAlgorithm(const vector<int> &tasksTime, const int &processorMaxTime, const int &maxGenerationsNumber,
                     const int &populationSize, const int &mutationProbability, const int& crossOverProbability,
                     const int& topKChromosoes);
    void run();

private:
    vector<int> tasksTime;
    int processorMaxTime;
    int maxGenerationsNumber;
    int populationSize;
    int mutationProbability;
    int crossOverProbability;
    int topKChromosoes;
    int tasksTotalTime;


    pair<Chromosome, Chromosome> crossOver(const Chromosome &firstChromosome, const Chromosome &secondChromosome);
    int calcFitness(const Chromosome &chromosome);
    void mutation(Chromosome &chromosome);
    vector<Chromosome> generateRandomPopulation(const int &popSize, const int& chromosomeSize);
    pair<unsigned int, unsigned int> makeSelection(const vector<unsigned int> &cumFitnessTable);
    vector<unsigned int> buildCumulativeFitnessTable(const vector<Chromosome>& chromosomes);
    pair<vector<Chromosome>,vector<Chromosome>> extractKBest(const int& k, vector<Chromosome>& originalVector);
    vector<Chromosome> filterInfeasableOffSprings(const vector<Chromosome>& offSprings);
    bool isFeasableFitness(const int& fitness);
};