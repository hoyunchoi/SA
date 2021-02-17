#pragma once

#include <vector>
#include <string>
#include <random>
#include <cmath>
#include <algorithm>
#include <iostream>
#include <fstream>

#include "../library/pcg_random.hpp"
#include "../library/CSV.hpp"

struct EuclideanNode{
    int m_index{0};
    double m_xPos{0.0};
    double m_yPos{0.0};
    int m_nextIndex{0};
    int m_previousIndex{0};

    //* Constructor
    EuclideanNode(){}
    EuclideanNode(const int& t_index, const double& t_xPos, const double& t_yPos) : m_index(t_index), m_xPos(t_xPos), m_yPos(t_yPos){}
};


/*
    Solve travelling salesman problem with simulated annealing method
*/
namespace SA{
    //* Network parameters
    int networkSize;
    std::vector<EuclideanNode> nodes;

    //* Random parameters
    pcg32 randomEngine;
    std::uniform_int_distribution<int> nodeDistribution;
    std::uniform_real_distribution<double> probDistribution(0,1);

    //* Simulation parameters
    double energy;
    double temperature;
    double cooling_power;
    unsigned maxIteration;
    unsigned maxIteration_constT ;
    std::string directory;

    //* Calculate energy : total length of the path
    double getEnergy(){
        double energy = 0;
        for (int index=0; index<networkSize; ++index){
            const int nextIndex = nodes[index].m_nextIndex;
            energy += std::sqrt(std::pow(nodes[index].m_xPos - nodes[nextIndex].m_xPos, 2.0) + std::pow(nodes[index].m_yPos - nodes[nextIndex].m_yPos, 2.0));
        }
        return energy;
    }

    double getEnergy(const int& t_index1, const int& t_index2){
        double energy = 0;
        int nextIndex;
        for (int index=0; index<networkSize; ++index){
            if (index == t_index1){
                nextIndex = nodes[t_index2].m_nextIndex;
            }
            else if (index == t_index2){
                nextIndex = nodes[t_index1].m_nextIndex;
            }
            else if (index == nodes[t_index1].m_previousIndex){
                nextIndex = t_index2;
            }
            else if (index == nodes[t_index2].m_previousIndex){
                nextIndex = t_index1;
            }
            else{
                nextIndex = nodes[index].m_nextIndex;
            }
            energy += std::sqrt(std::pow(nodes[index].m_xPos - nodes[nextIndex].m_xPos, 2.0) + std::pow(nodes[index].m_yPos - nodes[nextIndex].m_yPos, 2.0));
        }
        return energy;
    }

    //* save functions
    void savePos(const std::string& t_fileName = ""){
        //* File name is given: print into the file
        if (t_fileName.size()){
            std::ofstream file(t_fileName);
            for (const EuclideanNode& node : nodes){
                file << std::fixed << std::setprecision(10) << node.m_xPos << "," << node.m_yPos << "\n";
            }
            file.close();
        }
        //* File name is not given: print to terminal
        else{
            for (const EuclideanNode& node : nodes){
                std::cout << std::fixed << std::setprecision(10) << node.m_xPos << "," << node.m_yPos << "\n";
            }
        }
    }

    void saveData(const std::string& t_fileName = ""){
        //* File name is given: print into the file
        if (t_fileName.size()){
            std::ofstream file(t_fileName);
            for (const EuclideanNode& node : nodes){
                file << node.m_nextIndex;
                node.m_index != networkSize-1 ? file << "," : file << "," << temperature << "," << std::setprecision(10) << energy << "\n";
            }
            file.close();
        }
        //* File name is not given: print to terminal
        else{
            for (const EuclideanNode& node : nodes){
                std::cout << node.m_nextIndex;
                node.m_index != networkSize-1 ? std::cout << "," : std::cout << "," << temperature << "," << std::setprecision(10) << energy << "\n";
            }
        }
    }

    void saveData(std::ofstream& t_file){
        for (const EuclideanNode& node : nodes){
            t_file << node.m_nextIndex;
            node.m_index != networkSize-1 ? t_file << "," : t_file << "," << temperature << "," << std::setprecision(10) << energy << "\n";
        }
    }

    //* Get data
    const std::vector<EuclideanNode> setInitial(){
        std::vector<EuclideanNode> initial;

        //* Read position file
        std::vector<std::vector<double>> posVec;
        CSV::read(directory + "position.txt", posVec);
        for (int index=0; index<networkSize; ++index){
            EuclideanNode node(index, posVec[index][0], posVec[index][1]);
            initial.emplace_back(node);
        }

        //* Read last part of order file
        std::ifstream orderFile(directory + "data.txt");
        std::string line;
        orderFile.seekg(-2, std::ios_base::end);
        bool findingNewLine = true;
        while (findingNewLine){
            char temp;
            orderFile.get(temp);
            if (temp == '\n'){
                findingNewLine = false;
            }
            else{
                orderFile.seekg(-2, std::ios_base::cur);
            }
        }
        getline(orderFile, line);
        for (int index=0; index<networkSize; ++index){
            const int nextIndex = std::stoi(line.substr(0, line.find(",")));
            initial[index].m_nextIndex = nextIndex;
            initial[nextIndex].m_previousIndex = index;
            line = line.substr(line.find(",")+1);
        }
        orderFile.close();

        return initial;
    }



    //* Initialize
    void initialize(const int& t_randomEngineSeed, const std::vector<EuclideanNode>& t_initial={}){
        //* Generate random engine
        randomEngine.seed(t_randomEngineSeed);
        nodeDistribution.param(std::uniform_int_distribution<int>::param_type(0, networkSize-1));

        //* Set nodes
        nodes.clear(); nodes.reserve(networkSize);
        std::vector<int> randomIndex(networkSize);
        for (int index=0; index<networkSize; ++index){
            randomIndex[index] = index;
        }
        std::shuffle(randomIndex.begin(), randomIndex.end(), randomEngine);

        //* No initial condition is given. Set the ramdom position and order to default
        if (t_initial.size() != networkSize){
            for (int index=0; index<networkSize; ++index){
                const double theta = 2.0*M_PI/networkSize*randomIndex[index];
                EuclideanNode node(index, std::cos(theta), std::sin(theta));
                node.m_nextIndex = index+1;
                node.m_previousIndex = index-1;
                nodes.emplace_back(node);
            }
            nodes[networkSize-1].m_nextIndex = 0;
            nodes[0].m_previousIndex = networkSize-1;
            energy = getEnergy();

            savePos(directory + "position.txt");
            saveData(directory + "data.txt");
        }
        //* Initial condition is given. Set the position and order
        else{
            nodes = t_initial;
            energy = getEnergy();
        }

        //* Calculate initial energy
    }//* End of function SA::initialize


    //* Run Simulated Annealing with constant temperature
    void run_constT(std::ofstream& t_file){
        for (unsigned iter_constT=0; iter_constT<maxIteration_constT; ++iter_constT){
            //* Randomly choose two node to swap
            int index1 = nodeDistribution(randomEngine);
            int index2;
            do {
                index2 = nodeDistribution(randomEngine);
            } while(index2 == index1 || index2 == nodes[index1].m_nextIndex || index2 == nodes[index1].m_previousIndex);

            //* Get energy difference
            const double newEnergy = getEnergy(index1, index2);
            const double deltaE = newEnergy - energy;

            //* If dE > 0, don't change state with probability 1-exp(-dE/T). Otherwise, swap to nodes
            if (probDistribution(randomEngine) < std::exp(-1.0 * deltaE / temperature)){
                int temp1, temp2;
                temp1 = nodes[index1].m_previousIndex; temp2 = nodes[index2].m_previousIndex;
                std::swap(nodes[temp1].m_nextIndex, nodes[temp2].m_nextIndex);
                temp1 = nodes[index1].m_nextIndex; temp2 = nodes[index2].m_nextIndex;
                std::swap(nodes[temp1].m_previousIndex, nodes[temp2].m_previousIndex);
                std::swap(nodes[index1].m_nextIndex, nodes[index2].m_nextIndex);
                std::swap(nodes[index1].m_previousIndex, nodes[index2].m_previousIndex);
                energy = newEnergy;
            }
            saveData(t_file);
        }
    }//* End of funtion SA::run_T

    //* Run Simulated Annealing with exponential cooling strategy
    void run(){
        std::ofstream dataFile(directory + "data.txt", std::ios_base::app);

        for(unsigned iter=0; iter<maxIteration; ++iter){
            run_constT(dataFile);
            temperature /= cooling_power;
        }
    }//* End of function SA::run
}//* End of namespace SA