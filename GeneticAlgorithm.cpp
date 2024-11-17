// system imports
#include <functional>
#include <algorithm>

//
#include "GeneticAlgorithm.h"
#include "RandomNumGeneratorHelper.h"
using namespace std;

#define OFFSPRING_PERCENTAGE 0.7

//! 1. generate random population
//! 2. elitism (divide into best , population - best)
//!     2.i get K best chromosomes
//!     2.ii divide into best , remaining
//! 3. makeSelection (roullette Wheel)
//! 4. crossOverMutation
//! 5. mutation
//! 6. adjust best array in comparison to the new offSprings
//! 7. SSGA
//! 8. create new population
//! 9. REPEAT from step 2


bool GeneticAlgorithm::isFeasableFitness(const int &fitness)
{
    return tasksTotalTime - fitness <= processorMaxTime;
}

vector<Chromosome> GeneticAlgorithm::filterInfeasableOffSprings(const vector<Chromosome> &offSprings)
{
    vector<Chromosome> filteredChromosomes;
    for (const auto &o : offSprings)
    {
        int fitness = calcFitness(o);
        if (isFeasableFitness(fitness))
        {
            filteredChromosomes.push_back(o);
        }
    }
    return filteredChromosomes;
}

pair<vector<Chromosome>, vector<Chromosome>> GeneticAlgorithm::extractKBest(
    const int &k,
    vector<Chromosome> &originalVector)
{
    sort(originalVector.begin(), originalVector.end(),
         [&](const Chromosome &first, const Chromosome &second)
         {
             return this->calcFitness(first) > this->calcFitness(second);
         });
    vector<Chromosome> best;
    vector<Chromosome> remaining;
    for (int i = 0; i < originalVector.size(); i++)
    {
        if (i < k)
        {
            best.push_back(originalVector[i]);
        }
        else
        {
            remaining.push_back(originalVector[i]);
        }
    }
    return {best, remaining};
}

template <typename T>
std::vector<T> concatenateVectors(const std::vector<T> &first)
{
    return first;
}

template <typename T, typename... Vectors>
std::vector<T> concatenateVectors(const std::vector<T> &first, const Vectors &...vectors)
{
    std::vector<T> result = first;
    (result.insert(result.end(), vectors.begin(), vectors.end()), ...);
    return result;
}

void GeneticAlgorithm::run()
{
    // 1.generate random population
    vector<Chromosome> population = generateRandomPopulation(this->populationSize, this->tasksTime.size());
    for (int currGeneration = 0; currGeneration < this->maxGenerationsNumber; currGeneration++)
    {
        //! 2.elitism

        //! 2.i sort original population according to fitness value in descending order to get K best choromosomes
        //! 2.ii separate into KBest , remaining
        pair<vector<Chromosome>, vector<Chromosome>> extracted = extractKBest(this->topKChromosoes, population);
        vector<Chromosome> bestChromosomes(extracted.first);
        vector<Chromosome> parentChromosomes(extracted.second);

        //! 3. makeSelection using roullette wheel
        vector<unsigned int> cumulativeFreqTable = buildCumulativeFitnessTable(parentChromosomes);
        vector<Chromosome> offSprings;
        for (int i = 0; i < parentChromosomes.size(); i++)
        {
            pair<unsigned int, unsigned int> indecies = makeSelection(cumulativeFreqTable);

            //! 4. crossover mutation
            pair<Chromosome, Chromosome> offspring = crossOver(parentChromosomes[indecies.first], parentChromosomes[indecies.second]);

            //! 5. mutation
            mutation(offspring.first);
            mutation(offspring.second);

            offSprings.push_back(offspring.first);
            offSprings.push_back(offspring.second);
        }
        offSprings = filterInfeasableOffSprings(offSprings);

        //! 6. picking most promosing among offspring , best

        vector<Chromosome> bestWithOffSprings = concatenateVectors(bestChromosomes, offSprings);
        extracted = extractKBest(this->topKChromosoes, bestWithOffSprings);
        bestChromosomes = extracted.first;
        offSprings = extracted.second;

        //! 7. SSGA picking subset of offspring (70%) , subset of parents(30%)
        const int remainingSize = parentChromosomes.size();
        const int offSpringsPercentage = min(int(remainingSize * OFFSPRING_PERCENTAGE), int(offSprings.size()));
        const int parentPercentage = remainingSize - offSpringsPercentage;
        offSprings = extractKBest(offSpringsPercentage, offSprings).first;
        parentChromosomes = extractKBest(parentPercentage, parentChromosomes).first;

        //! 8. updating population
        population = concatenateVectors(bestChromosomes, offSprings, parentChromosomes);
    }
    Chromosome optimalChromosome = population[0];
    cout << "optimal chromosome is: ";
    optimalChromosome.print();

    cout << "core 1 tasks: ";
    for(int i = 0; i < optimalChromosome.genes.size(); i++){
        if(optimalChromosome.genes[i] == 1){
            cout << tasksTime[i] << ' ';
        }
    }
    cout << '\n';
    cout << "core 2 tasks: ";
    for(int i = 0; i < optimalChromosome.genes.size(); i++){
        if(optimalChromosome.genes[i] == 0){
            cout << tasksTime[i] << ' ';
        }
    }
    cout << '\n';
    cout <<"Total Time "<< tasksTotalTime - calcFitness(optimalChromosome) << "s \n";
    cout <<"FitnessScore "<< calcFitness(optimalChromosome) << '\n';
}

// ################################## private Functions ##################################
GeneticAlgorithm::GeneticAlgorithm(const vector<int> &tasksTime, const int &processorMaxTime,
                                   const int &maxGenerationsNumber, const int &populationSize,
                                   const int &mutationProbability, const int &crossOverProbability,
                                   const int &topKChromosoes)
{
    this->tasksTime = tasksTime;
    this->processorMaxTime = processorMaxTime;
    this->maxGenerationsNumber = maxGenerationsNumber;
    this->populationSize = populationSize;
    this->mutationProbability = mutationProbability;
    this->crossOverProbability = crossOverProbability;
    this->topKChromosoes = topKChromosoes;

    this->tasksTotalTime = 0;
    for (auto i : tasksTime)
        this->tasksTotalTime += i;
}

int GeneticAlgorithm::calcFitness(const Chromosome &chromosome)
{
    int processor1TotalTime = 0;
    int processor2TotalTime = 0;
    //! sum time taken by first processor, time taken by second processor
    for (int i = 0; i < chromosome.genes.size(); i++)
    {
        chromosome.genes[i] ? processor1TotalTime += this->tasksTime[i] : processor2TotalTime += this->tasksTime[i];
    }
    //! time is the max taken by the two processors
    return min(this->tasksTotalTime - processor1TotalTime, this->tasksTotalTime - processor2TotalTime);
}

void GeneticAlgorithm::mutation(Chromosome &chromosome)
{
    for (int i = 0; i < chromosome.genes.size(); i++)
    {
        int probability = RandomGenerator::generateRandomNumber();
        //! Flip bits(gene)
        if (probability <= this->mutationProbability)
        {
            chromosome.genes[i] = !chromosome.genes[i];
        }
    }
}

pair<Chromosome, Chromosome> GeneticAlgorithm::crossOver(const Chromosome &firstChromosome, const Chromosome &secondChromosome)
{
    unsigned int size = firstChromosome.genes.size();
    int probability = RandomGenerator::generateRandomNumber();

    //! no cross over happens
    if (probability > this->crossOverProbability)
    {
        return {firstChromosome, secondChromosome};
    }

    // point at which crossOver between 2 chromosome happens
    unsigned int crossOverPoint = RandomGenerator::generateRandomNumber(1, size - 1);
    Chromosome firstOffSpring{size};
    Chromosome secondOffSpring{size};
    for (int i = 0; i < size; i++)
    {
        //! if crossOverPoint is not yet passed , the put bits inf firstParent
        if (i < crossOverPoint)
        {
            firstOffSpring.genes[i] = firstChromosome.genes[i];
            secondOffSpring.genes[i] = secondChromosome.genes[i];
        }
        else
        {
            //! crossOverPoint is passed , then exchange bits,  perform actual crossOver
            firstOffSpring.genes[i] = secondChromosome.genes[i];
            secondOffSpring.genes[i] = firstChromosome.genes[i];
        } // 10101,
    }
    return make_pair(firstOffSpring, secondOffSpring);
}

pair<unsigned int, unsigned int> GeneticAlgorithm::makeSelection(const vector<unsigned int> &cummulativeFitnessTable)
{
    //! total sum is at the last index of the vector , since it is a prefix sum of all previous elements
    unsigned int totalSum = cummulativeFitnessTable.back();

    //! indecies of picked choromosomes
    int firstParentIndex = -1;
    int secondParentIndex = -2;

    //! find first prefixSome that is >= the required prefix sum
    function<unsigned int(const unsigned int)> findIndex = [&](const unsigned int requiredPrefixSum) -> unsigned int
    {
        for (int i = 0; i < cummulativeFitnessTable.size(); i++)
        {
            const unsigned int prefix = cummulativeFitnessTable[i];
            if (prefix >= requiredPrefixSum)
            {
                return i;
            }
        }
        return -1;
    };

    //! keep picking index for the selected chromosomes until 2 are picked, and not equal each other
    while (firstParentIndex == secondParentIndex || firstParentIndex == -1 || secondParentIndex == -2)
    {
        //! random numbers for required prefix some for both choromosomes to be picked
        unsigned int firstParentPrefixSum = RandomGenerator::generateRandomNumber(0, totalSum - 1);
        unsigned int secondParentPrefixSum = RandomGenerator::generateRandomNumber(0, totalSum - 1);
        if (firstParentIndex == -1)
        {
            firstParentIndex = findIndex(firstParentPrefixSum);
        }
        else
        {
            secondParentIndex = findIndex(secondParentPrefixSum);
        }
    }
    return make_pair(firstParentIndex, secondParentIndex);
}

vector<Chromosome> GeneticAlgorithm::generateRandomPopulation(const int &popSize, const int &chromosomeSize)
{
    vector<Chromosome> population;
    while (population.size() != popSize)
    {
        Chromosome chromosome(chromosomeSize);
        for (int j = 0; j < chromosomeSize; j++)
        {
            chromosome.genes[j] = RandomGenerator::generateRandomNumber(0, 1);
        }
        if (isFeasableFitness(calcFitness(chromosome)))
        {
            population.push_back(chromosome);
        }
    }
    return population;
}

vector<unsigned int> GeneticAlgorithm::buildCumulativeFitnessTable(const vector<Chromosome> &chromosomes)
{
    vector<unsigned int> fitnesses;
    for (const auto &chromosome : chromosomes)
    {
        unsigned int fitness = calcFitness(chromosome);
        fitnesses.push_back(fitness);
    }

    vector<unsigned int> prefixSums(fitnesses.size());
    prefixSums[0] = fitnesses[0];
    for (int i = 1; i < fitnesses.size(); i++)
    {
        prefixSums[i] = fitnesses[i] + prefixSums[i - 1];
    }
    return prefixSums;
}