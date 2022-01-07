#ifndef CAJETE_PLANT_SSA_HPP
#define CAJETE_PLANT_SSA_HPP 

#include "PlantTypes.hpp"
#include "PlantGrammar.hpp"

#include "ExpandedComplex2D.hpp"
#include "YAGL_Graph.hpp"

#include <random>
#include <map>

namespace Cajete 
{
namespace Plant 
{
//template <typename BucketType>
//void plant_model_ssa(
//        BucketType& bucket,
//        ExpandedComplex2D<>& geoplex2D, 
//        YAGL::Graph<mt_key_type, MT_NodeData>& system_graph) 
template <typename BucketType, typename GeoplexType, typename GraphType>
void plant_model_ssa(BucketType& bucket, GeoplexType& geoplex2D, GraphType& system_graph)
{
    
    double DELTA = 0.1;
    int NUM_INTERNAL_STEPS = 10;
    double DELTA_DELTA_T = DELTA / NUM_INTERNAL_STEPS;

    double delta_t, exp_sample, tau, geocell_propensity;
    std::size_t events;

    auto all_matches = microtubule_growing_end_matcher(system_graph, bucket.second); 
    
    delta_t = 0.0; events = 0;
    while(delta_t < DELTA) 
    {
        //reset tau
        tau = 0.0;

        //sample the exponential variable
        exp_sample = 0.1; //-log(1-uniform_sample)
        
        
        while(delta_t < DELTA && tau < exp_sample)
        {
            // STEP(0) : find all the matches 
            auto k = bucket.first;
            //auto all_matches = microtubule_growing_end_matcher(system_graph, bucket.second); 
            
            // STEP(1) : solve the system of ODES
            microtubule_ode_solver(all_matches, system_graph, k); //be careful with solving, it could lead
            //to a segmentation fault in the sorting phase if a parameter is solved out to out of bounds

            // STEP(2) : sum all of the rule propensities
             //zero the geocell_propensity
            geocell_propensity = 0.0;

            // STEP(3) : use forward euler to solve the TAU ODE
            tau += geocell_propensity*DELTA_DELTA_T;
            
            // STEP(4) : advance the loop timer
            delta_t += DELTA_DELTA_T;
        }
    

        // if we get over our threshold an event can occur
        if(tau > exp_sample) 
        {
            //determine which rule to file and fire it
            auto k = bucket.first;
            auto all_matches = microtubule_growing_end_matcher(system_graph, bucket.second); 
            microtubule_rule_firing(all_matches, system_graph, k);
        }
    }
}

template <typename MatchType, typename GraphType, typename KeyType>
void microtubule_ode_solver(MatchType& all_matches, GraphType& system_graph, KeyType& k)
{
    for(auto match : all_matches)
    {
       bool bad_match = false;
       //check match integrity
       for(auto key : match)
       {
            auto dtag = system_graph.findNode(key)->second.getData().tagND[0];
            if(dtag != k)
            {
                bad_match = true;
                //std::cout << "Bad match found, it'll be skipped!\n";
                break;
            }
       }

       if(!bad_match)
       {
            microtubule_growing_end_polymerize_solve(system_graph, match);
       }
    }
}

template <typename MatchType, typename GraphType, typename KeyType>
void microtubule_rule_firing(MatchType& all_matches, GraphType& system_graph, KeyType& k)
{
for(auto match : all_matches)
    {
       bool bad_match = false;
       //check match integrity
       for(auto key : match)
       {
            auto dtag = system_graph.findNode(key)->second.getData().tagND[0];
            if(dtag != k)
            {
                bad_match = true;
                //std::cout << "Bad match found, it'll be skipped!\n";
                break;
            }
       }

       if(!bad_match)
       {
            microtubule_growing_end_polymerize_rewrite(system_graph, match); 
            break;
       }
    }

}


} //end namespace plant 
} //end namespace cajete

#endif