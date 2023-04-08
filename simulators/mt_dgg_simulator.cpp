#include <iostream>
#include <string> 

#include "DggFactory.hpp"
#include "PlantModel.hpp"

#include "simdjson.h"

int main(int argc, char* argv[])
{   

    std::cout << "Running Microtubule Dynamic Graph Grammar Simulator\n";
    
    std::string settings_filename;
    if(argc == 1)
        settings_filename = "settings.json"; // by default
    if(argc == 2)
    {
        settings_filename = argv[1];
    }

    simdjson::ondemand::parser parser;
    simdjson::padded_string json = simdjson::padded_string::load(settings_filename);
    simdjson::ondemand::document settings_file = parser.iterate(json);

    Cajete::DggFactory<Cajete::PlantModel, simdjson::ondemand::document> plant_factory;
    
    auto num_simulations = 1;
   
    //TODO: add some sort of ensemble simulation mode
    for(auto i = 0; i < num_simulations; i++)
        plant_factory.execute(settings_file);
}
