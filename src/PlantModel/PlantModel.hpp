#ifndef PLANT_MODEL_HPP
#define PLANT_MODEL_HPP

#include <iostream>

#include "PlantTypes.hpp"
#include "PlantUtils.hpp"
#include "PlantGrammar.hpp"
#include "PlantSSA.hpp"

#include "DggModel.hpp"
#include "YAGL_Graph.hpp"
#include "YAGL_Node.hpp"
#include "YAGL_Algorithms.hpp"
#include "VtkWriter.hpp"

#include "ExpandedComplex2D.hpp"

#include "CartesianComplex2D.hpp"

#include "CartesianHashFunctions.hpp"

#include <map>
#include <random>
#include <chrono>
#include <string>
#include <filesystem>
#include <numeric>

namespace Cajete
{
    double compute_correlation(const std::vector<double>& x, const std::vector<double>& y) {
        // Compute means
        double mean_x = std::accumulate(x.begin(), x.end(), 0.0) / x.size();
        double mean_y = std::accumulate(y.begin(), y.end(), 0.0) / y.size();

        // Compute variances
        double var_x = std::accumulate(x.begin(), x.end(), 0.0,
            [mean_x](double acc, double xi) { return acc + (xi - mean_x) * (xi - mean_x); }) / (x.size() - 1);
        double var_y = std::accumulate(y.begin(), y.end(), 0.0,
            [mean_y](double acc, double yi) { return acc + (yi - mean_y) * (yi - mean_y); }) / (y.size() - 1);

        // Check for zero variance
        if (var_x == 0 || var_y == 0) {
            return 0; // or any other value you choose to return in this case
        }

        // Compute covariance
        double cov = 0;
        for (size_t i = 0; i < x.size(); ++i) {
            cov += (x[i] - mean_x) * (y[i] - mean_y);
        }
        cov /= x.size() - 1;

        // Compute correlation coefficient
        double corr = cov / std::sqrt(var_x * var_y);

        return corr;
}


    template<typename GraphType>
    double compute_correlation(GraphType& system_graph)
    {
        std::vector<double> x;
        std::vector<double> y;
        for(auto& node : system_graph.getNodeSetRef())
        {
            auto& u = node.second.getData().unit_vec;
            x.push_back(u[0]);
            y.push_back(u[1]);
            //y.push_back(abs(u[1]));
        }
        
        return compute_correlation(x, y);
        /*std::cout << "C: " << correlation;
        return correlation;
        auto avg_x = std::reduce(x.begin(), x.end(), 0.0) / x.size();
        auto avg_y = std::reduce(y.begin(), y.end(), 0.0) / y.size();
        //std::cout << "avg_x: " << avg_x << ", avg_y: " << avg_y << "\n";

        double e_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
        double n = x.size();
        double e_x_y = n * avg_x * avg_y;
        std::cout << e_xy << "\n";
        std::cout << e_x_y << "\n";
        std::vector<double> r1(x.size());
        std::vector<double> r2(y.size());
        
        std::transform(x.begin(), x.end(), r1.begin(),
                [&](double x) { return (x - avg_x)*(x - avg_x); });
        double s1 = std::sqrt(std::accumulate(r1.begin(), r1.end(), 0.0));
        std::cout << s1 << "\n"; 
        std::transform(y.begin(), y.end(), r2.begin(), 
                [&](double y){ return (y - avg_y)*(y - avg_y); });
        double s2 = std::sqrt(std::accumulate(r2.begin(), r2.end(), 0.0));
        std::cout << s2 << "\n";
        /double correlation = (e_xy - e_x_y) / (s1*s2);
        std::cout << "C: " << correlation << "\n";
        return correlation;*/
    }

    struct Parameters
    {
        std::string EXPERIMENT_NAME;
        double DELTA;
        double DELTA_DELTA_T;
        int NUM_INTERNAL_STEPS;
        std::size_t CELL_NX;
        std::size_t CELL_NY;
        double CELL_DX;
        double CELL_DY;
        bool GHOSTED;
        std::size_t NUM_MT;
        double MT_MIN_SEGMENT_INIT;
        double MT_MAX_SEGMENT_INIT;
        std::size_t NUM_STEPS;
        double LENGTH_DIV_FACTOR;
        double DIV_LENGTH;
        double DIV_LENGTH_RETRACT;
        double V_PLUS;
        double V_MINUS;
        double SIGMOID_K;
        double TOTAL_TIME;
        double MAXIMAL_REACTION_RADIUS;
        double DELTA_T_MIN;
        std::size_t CHECKPOINT_FREQUENCY;
        double CATASTROPHE_RATE_FACTOR; 
        double ZIPPER_RATE_FACTOR;
        double CROSSOVER_RATE_FACTOR;
        double RESCUE_RATE_FACTOR;
        double INSTABILITY_RATE_FACTOR;
        double WOBBLE_ANGLE_FACTOR;
        double CRITICAL_ANGLE;
    };
  
    template <typename DataType>
    void print_numpy_array_stats(DataType& data, std::string var_name)
    {
        std::cout << var_name << " = np.asarray([";
            for(auto i = 0; i < data.size(); i++)
            {
                if(i != 0 && i % 20 == 0) std::cout << "\n";
                if(i != data.size() - 1)
                    std::cout << data[i] << ", ";
                else
                    std::cout << data[i];
            }
            std::cout << "]);\n";
    }

    template <typename ParamType, typename InterfaceType>
    void set_parameters(ParamType& settings, InterfaceType& interface)
    {
        std::string_view temp = interface["META"]["EXPERIMENT"];
        settings.EXPERIMENT_NAME = static_cast<std::string>(temp);

        std::cout << settings.EXPERIMENT_NAME << "+++\n";
                settings.CELL_NX = std::size_t(interface["SETTINGS"]["CELL_NX"]); 
        settings.CELL_NY = std::size_t(interface["SETTINGS"]["CELL_NY"]);
        
        settings.CELL_DX = double(interface["SETTINGS"]["CELL_DX"]);
        settings.CELL_DY = double(interface["SETTINGS"]["CELL_DY"]);
        settings.GHOSTED = bool(interface["SETTINGS"]["GHOSTED"]);

        settings.NUM_MT = std::size_t(interface["SETTINGS"]["NUM_MT"]);
        settings.MT_MIN_SEGMENT_INIT = double(interface["SETTINGS"]["MT_MIN_SEGMENT_INIT"]);
        settings.MT_MAX_SEGMENT_INIT = double(interface["SETTINGS"]["MT_MAX_SEGMENT_INIT"]);

        settings.LENGTH_DIV_FACTOR = double(interface["SETTINGS"]["LENGTH_DIV_FACTOR"]);
        settings.DIV_LENGTH = double(interface["SETTINGS"]["DIV_LENGTH"]);
        settings.DIV_LENGTH_RETRACT = double(interface["SETTINGS"]["DIV_LENGTH_RETRACT"]);
        settings.V_PLUS = double(interface["SETTINGS"]["V_PLUS"]);
        settings.V_MINUS = double(interface["SETTINGS"]["V_MINUS"]);

        settings.SIGMOID_K = double(interface["SETTINGS"]["SIGMOID_K"]);
        
        settings.MAXIMAL_REACTION_RADIUS = settings.DIV_LENGTH*2.0;
        
        //Simulate until the specified unit time
        settings.TOTAL_TIME = double(interface["SETTINGS"]["TOTAL_TIME"]);
        settings.NUM_INTERNAL_STEPS = double(interface["SETTINGS"]["NUM_INTERNAL_STEPS"]);
        //Delta should be big, but not to big. In this case, the maximum amount of time it would
        //take one MT to grow a single unit of MT
        settings.DELTA = 
            0.25*settings.MAXIMAL_REACTION_RADIUS / std::max(settings.V_PLUS, settings.V_MINUS);
        //The internal step of the solver should be at least this small
        settings.DELTA_DELTA_T = settings.DELTA / settings.NUM_INTERNAL_STEPS; 
        settings.DELTA_T_MIN = settings.DELTA_DELTA_T;
        settings.NUM_STEPS = settings.TOTAL_TIME / settings.DELTA;
        settings.CHECKPOINT_FREQUENCY = std::size_t(interface["SETTINGS"]["CHECKPOINT_FREQUENCY"]);
        settings.CATASTROPHE_RATE_FACTOR = double(interface["RATE_FACTORS"]["CATASTROPHE_RATE_FACTOR"]);
        settings.ZIPPER_RATE_FACTOR = double(interface["RATE_FACTORS"]["ZIPPER_RATE_FACTOR"]);
        settings.CROSSOVER_RATE_FACTOR = double(interface["RATE_FACTORS"]["CROSSOVER_RATE_FACTOR"]);
        settings.RESCUE_RATE_FACTOR = double(interface["RATE_FACTORS"]["RESCUE_RATE_FACTOR"]);
        settings.INSTABILITY_RATE_FACTOR = double(interface["RATE_FACTORS"]["INSTABILITY_RATE_FACTOR"]);
        settings.WOBBLE_ANGLE_FACTOR = double(interface["RATE_FACTORS"]["WOBBLE_ANGLE_FACTOR"]);
        settings.CRITICAL_ANGLE = double(interface["RATE_FACTORS"]["CRITICAL_ANGLE"]);
    }

    //Models are inteded to be designed based on the 
    //DggModel specification. Right now it's very loose
    //and capable of handling almost anything
    template <typename InterfaceType>
    class PlantModel : public DggModel<InterfaceType> {
    public:
        using key_type = Plant::mt_key_type;
        using gplex_key_type = typename Cajete::ExpandedComplex2D<>::graph_type::key_type;
        using data_type = Plant::MT_NodeData;
        using graph_type = YAGL::Graph<key_type, data_type>;
        using node_type = typename graph_type::node_type;

        void init(InterfaceType& interface) override {

            std::cout << "\n\n-----------------------------------------------------------------------\n";
            //TODO: implement timers to monitor start up phase
            std::cout << "Initializing the plant model simulation\n";
            
            std::cout << "Parsing the input interface and setting configuration settings\n";
            //TODO: handle the interface input
            set_parameters(settings, interface); 
            
            std::cout << "Cleaning up old results folder if it exists and creating a new one\n";
            results_dir_name = settings.EXPERIMENT_NAME + "_results";
            std::filesystem::remove_all(results_dir_name);
            std::filesystem::create_directory(results_dir_name);
            
            std::cout << "Generating the expanded cell complex\n";
            geoplex2D.init(settings.CELL_NX, 
                    settings.CELL_NY, 
                    settings.CELL_DX, 
                    settings.CELL_DY, 
                    settings.GHOSTED); //ghosted
            std::cout << geoplex2D;
            
            //Save expanded cell complex graph
            Cajete::VtkFileWriter<typename Cajete::ExpandedComplex2D<>::types::graph_type> writer;
            writer.save(geoplex2D.getGraph(), results_dir_name+"/factory_geoplex");
            
            std::cout << "Initializing the system graph\n";
            Plant::microtubule_uniform_scatter(system_graph, geoplex2D, settings); 
            
            //std::cout << "Generating the grammar\n";
            //TODO: implement a grammar setup phase
    
        }

        void run() override {
            std::cout << "Running the plant model simulation\n";

            Cajete::VtkFileWriter<graph_type> vtk_writer;
            std::vector<std::size_t> con_com;
            con_com.push_back(YAGL::connected_components(system_graph));
            std::vector<std::size_t> total_nodes;
            std::vector<std::size_t> type_counts[5];
            std::vector<double> time_count; 
            std::vector<double> correlation;
            std::size_t junction_count = 0;
            std::size_t positive_count = 0;
            std::size_t negative_count = 0;
            std::size_t zipper_count = 0;
            std::size_t intermediate_count = 0;
            for(auto iter = system_graph.node_list_begin(); 
                    iter != system_graph.node_list_end(); iter++) {
                    auto itype = iter->second.getData().type;
                    if(itype == Plant::negative)
                        negative_count++;
                    if(itype == Plant::positive)
                        positive_count++;
                    if(itype == Plant::intermediate)
                        intermediate_count++;
                    if(itype == Plant::junction)
                        junction_count++;
                    if(itype == Plant::zipper)
                        zipper_count++;
            }
            type_counts[Plant::negative].push_back(negative_count);
            type_counts[Plant::positive].push_back(positive_count);
            type_counts[Plant::intermediate].push_back(intermediate_count);
            type_counts[Plant::junction].push_back(junction_count);
            type_counts[Plant::zipper].push_back(zipper_count);
            
            total_nodes.push_back(negative_count +
                    positive_count + intermediate_count + junction_count + zipper_count);
            //compute the intial correlation
            correlation.push_back(compute_correlation(system_graph));  
            std::string title = results_dir_name+"/simulation_step_";
            std::cout << "Saving the initial state of the system graph\n";
            vtk_writer.save(system_graph, title+std::to_string(0));

            //TODO: move the simulation algorithm to its own class
            //is this where we run the simulation?
            double total_run_time = 0.0;
            for(auto i = 1; i <= settings.NUM_STEPS; i++)
            {
                std::cout << "Running step " << i << std::endl;
                
                std::map<gplex_key_type, std::vector<key_type>> bucketsND[3];
                std::size_t complementND[3] = {0, 0, 0};
                
                double dim_time = 0.0;
                double tot_time = 0.0;
                std::cout << "Binning the graph into 2D partitions\n";
                
                auto start = std::chrono::high_resolution_clock::now();
                //TODO: optimize this to work only for 2D
                Cajete::expanded_cartesian_complex_sort_stl(bucketsND, complementND, geoplex2D, system_graph);

                //TODO: hoist the initial system pattern matcher code outside of the ssa phase 
                //      since we only want to smartly recomputed rule updates 
                auto stop = std::chrono::high_resolution_clock::now();
                auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
                std::cout << "Sorting took " << duration.count() << " milliseconds\n";
                dim_time += duration.count();
                 
                std::cout << "Running the Hybrid ODES/SSA inner loop 2D phase\n";
                for(auto& bucket : bucketsND[0])
                {
                   auto k = bucket.first; //check to see if this is a domain to simulate
                   if(geoplex2D.getGraph().findNode(k)->second.getData().interior)
                   {
                        auto start = std::chrono::high_resolution_clock::now();
                        plant_model_ssa(bucket, geoplex2D, system_graph, settings);
                        auto stop = std::chrono::high_resolution_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
                        //std::cout << "Cell " << k << " took " << duration.count() << " milliseconds\n";
                        dim_time += duration.count();
                   }
                }
                
                tot_time += dim_time;
                std::cout << "2D took " << dim_time << " milliseconds\n";

                //TODO: remove, right now connected_components should remain constant with only 
                //growth rules
                //std::cout << "----------------\n";
                //std::cout << "CC: " << YAGL::connected_components(system_graph); 
                //std::cout << "\n---------------\n";
                                
                std::cout << "Synchronizing work\n";
                //TODO: this is where a barrier would be for a parallel code
                for(auto& item : bucketsND) item.clear();
               
                dim_time = 0.0;
                
                std::cout << "Binning the graph into 1D partitions\n";
                start = std::chrono::high_resolution_clock::now();
                //TODO: optimize this to work only for 1D
                Cajete::expanded_cartesian_complex_sort_stl(bucketsND, complementND, geoplex2D, system_graph);
                stop = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
                std::cout << "Sorting took " << duration.count() << " milliseconds\n";
                dim_time += duration.count();

                std::cout << "Running the Hybrid ODES/SSA inner loop 1D phase\n";
                for(auto& bucket : bucketsND[1])
                {
                   auto k = bucket.first; //check to see if this is a domain to simulate
                   if(geoplex2D.getGraph().findNode(k)->second.getData().interior)
                   {
                        auto start = std::chrono::high_resolution_clock::now();
                        plant_model_ssa(bucket, geoplex2D, system_graph, settings);
                        auto stop = std::chrono::high_resolution_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
                        //std::cout << "Cell " << k << " took " << duration.count() << " milliseconds\n";
                        dim_time += duration.count();
                   }
                }
                
                tot_time += dim_time;
                std::cout << "1D took " << dim_time << " milliseconds\n";

                std::cout << "Synchronizing work\n";
                //TODO: this is where a barrier would be for a parallel code
                for(auto& item : bucketsND) item.clear();
            
                dim_time = 0.0;

                std::cout << "Binning the graph into 0D partitions\n";
                start = std::chrono::high_resolution_clock::now();
                //TODO: optimize this to work only for 0D
                Cajete::expanded_cartesian_complex_sort_stl(bucketsND, complementND, geoplex2D, system_graph);
                stop = std::chrono::high_resolution_clock::now();
                duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
                std::cout << "Sorting took " << duration.count() << " milliseconds\n";
                dim_time += duration.count();

                std::cout << "Running the Hybrid ODES/SSA inner loop 0D phase\n";
                for(auto& bucket : bucketsND[2])
                {
                   auto k = bucket.first; //check to see if this is a domain to simulate
                   if(geoplex2D.getGraph().findNode(k)->second.getData().interior)
                   {
                        auto start = std::chrono::high_resolution_clock::now();
                        plant_model_ssa(bucket, geoplex2D, system_graph, settings);
                        auto stop = std::chrono::high_resolution_clock::now();
                        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
                        //std::cout << "Cell " << k << " took " << duration.count() << " milliseconds\n";
                        dim_time += duration.count();
                   }
                }
                
                tot_time += dim_time;
                std::cout << "0D took " << dim_time << " milliseconds\n";
        
                std::cout << "Synchronizing work\n";
                //TODO: this is where a barrier would be for a parallel code
                for(auto item : bucketsND) item.clear();
        
                std::cout << "Running the checkpointer\n";
                if(i % settings.CHECKPOINT_FREQUENCY == 0)
                {
                    vtk_writer.save(system_graph, title+std::to_string(i));
                    //compute the correlation 
                    correlation.push_back(compute_correlation(system_graph));
                }
                con_com.push_back(YAGL::connected_components(system_graph)); 
                
                std::size_t junction_count = 0;
                std::size_t positive_count = 0;
                std::size_t negative_count = 0;
                std::size_t zipper_count = 0;
                std::size_t intermediate_count = 0;
                for(auto iter = system_graph.node_list_begin(); 
                        iter != system_graph.node_list_end(); iter++) {
                        auto itype = iter->second.getData().type;
                        if(itype == Plant::negative)
                            negative_count++;
                        if(itype == Plant::positive)
                            positive_count++;
                        if(itype == Plant::intermediate)
                            intermediate_count++;
                        if(itype == Plant::junction)
                            junction_count++;
                        if(itype == Plant::zipper)
                            zipper_count++;
                }
                type_counts[Plant::negative].push_back(negative_count);
                type_counts[Plant::positive].push_back(positive_count);
                type_counts[Plant::intermediate].push_back(intermediate_count);
                type_counts[Plant::junction].push_back(junction_count);
                type_counts[Plant::zipper].push_back(zipper_count);
                
                total_nodes.push_back(negative_count +
                        positive_count + intermediate_count + junction_count + zipper_count);
                
                std::cout << "Total dimensional time is " << tot_time << " milliseconds\n";
                time_count.push_back(tot_time);
                total_run_time += tot_time;
            }
            std::cout << "-----------------------------------------------------------------------\n\n";
            std::cout << "total run time of simulation including io: " << (total_run_time/1000)/60 << " minutes\n"; 
            print_numpy_array_stats(con_com, "con_com");
            print_numpy_array_stats(type_counts[Plant::negative], "negative");
            print_numpy_array_stats(type_counts[Plant::positive], "positive");
            print_numpy_array_stats(type_counts[Plant::intermediate], "intermediate");
            print_numpy_array_stats(type_counts[Plant::junction], "junction");
            print_numpy_array_stats(type_counts[Plant::zipper], "zipper");
            print_numpy_array_stats(total_nodes, "total_nodes");
            print_numpy_array_stats(time_count, "time_count");
            print_numpy_array_stats(correlation, "correlation");
        }


    private:
        Parameters settings;
        CartesianComplex2D<> cplex2D;
        ExpandedComplex2D<> geoplex2D;
        YAGL::Graph<Plant::mt_key_type, Plant::MT_NodeData> system_graph; 
        std::string results_dir_name;
};

} //end namespace Cajete

#endif 
