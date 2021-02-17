#include <vector>
#include <filesystem>

#include "TSP.hpp"


int main(int argc, char *argv[]){
    //* Get input
    SA::networkSize = std::stoi(argv[1]);
    SA::temperature = std::stod(argv[2]);
    const int randomEngineSeed = std::stoi(argv[3]);

    //* Set simulation parameters
    SA::cooling_power = 2.0;
    SA::maxIteration_constT = 100;
    SA::maxIteration = 7;

    //* Directory and File Name
    SA::directory = "../data/TSP/N" + std::to_string(SA::networkSize) + "-" + std::to_string(randomEngineSeed) + "/";

    //* Generate Directory or get initial condition
    namespace fs = std::filesystem;
    std::vector<EuclideanNode> initial;
    if (!fs::exists(SA::directory)){
        fs::create_directories(SA::directory);
    }
    else{
        initial = SA::setInitial();
    }

    //* Initialize and run
    SA::initialize(randomEngineSeed, initial);
    SA::run();

    return 0;
}