#pragma once
#include <iostream>
#include <vector>

class Chromosome{
public :
    Chromosome(unsigned int chromosomeSize): genes(chromosomeSize){}
    std::vector<bool> genes;

    void print(){
        for(int i = 0 ;i<genes.size();i++){
            std::cout<<genes[i];
        }
        std::cout<<"\n";
    }
};