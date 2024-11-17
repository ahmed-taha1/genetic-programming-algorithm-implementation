#include <iostream>
#include "GeneticAlgorithm.h"
#include <fstream>
#include <vector>

using namespace std;

int main() {
    ifstream inputFile("test.txt");  // Open the input file

    if (!inputFile) {
        cerr << "Error: Unable to open file test.txt" << endl;
        return 1;
    }

    int testsNumber;
    inputFile >> testsNumber;  // Read number of test cases from the file

    for(int i = 1; i <= testsNumber; i++) {
        int processorMaxTime;
        inputFile >> processorMaxTime;

        int numOfTasks;
        inputFile >> numOfTasks;

        vector<int> tasksTime(numOfTasks);
        for (auto& t : tasksTime)
            inputFile >> t;

        int maxGenerationsNumber = 20;
        int populationSize = 50;
        int mutationProbability = 50;
        int crossOverProbability = 50;
        int topKChromosomes = 3;

        GeneticAlgorithm geneticAlgorithm(
            tasksTime, processorMaxTime, maxGenerationsNumber, populationSize,
            mutationProbability, crossOverProbability, topKChromosomes
        );
        cout << "TEST #" << i << " \n";
        geneticAlgorithm.run();
        cout << "\n";
    }

    inputFile.close();  // Close the file when done
    return 0;
}
