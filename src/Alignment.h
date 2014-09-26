/*
 * Alignment.h
 *
 *  Created on: Jul 20, 2014
 *      Author: kgori
 */

#ifndef DISTANCE_H_
#define DISTANCE_H_

#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Seq/DistanceMatrix.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Simulation/HomogeneousSequenceSimulator.h>

#include <iostream>
#include <map>
#include <memory>
#include <vector>

using namespace std;
using namespace bpp;

class Alignment {
    public :
        Alignment(vector<pair<string, string>>, string datatype);
        Alignment(string filename, string file_format, string datatype, bool interleaved=true);
        Alignment(string filename, string file_format, string datatype, string model_name, bool interleaved=true);
        void read_alignment(string filename, string file_format, string datatype, bool interleaved=true);
        void set_model(string model_name);
        void set_alpha(int ncat=4, double alpha=1.0);
        void set_rates(vector<double>, string order="acgt");
        void set_frequencies(vector<double>);
        double get_alpha();
        vector<double> get_rates(string order);
        vector<double> get_frequencies();
        vector<string> get_names();
        bool is_dna();
        bool is_protein();
        string get_model();
        
        // Distance
        void compute_distances();
        void fast_compute_distances();
        void set_distance_matrix(vector<vector<double>> matrix);
        string get_nj_tree();
        string get_nj_tree(vector<vector<double>> matrix);
        vector<vector<double>> get_distances();
        vector<vector<double>> get_variances();
        vector<vector<double>> get_distance_variance_matrix();
        
        // Likelihood
        void initialise_likelihood(string tree);
        void optimise_parameters(bool fix_branch_lengths);
        double get_likelihood();
        string get_tree();

        // Simulator
        void write_alignment(string filename, string file_format, bool interleaved);
        void set_simulator(string tree);
        vector<pair<string, string>> simulate(unsigned int nsites, string tree);
        vector<pair<string, string>> simulate(unsigned int nsites);
        vector<pair<string, string>> get_simulated_sequences();

    private :
        void _set_dna();
        void _set_protein();
        void _write_phylip(string filename, bool interleaved);
        void _write_fasta(string filename);
        map<int, double> _vector_to_map(vector<double>);
        void _check_distances_exist();
        void _check_compatible_model(string datatype, string model);
        void _clear_distances();
        void _clear_likelihood();
        bool _is_file(string filename);
        bool _is_tree_string(string tree_string);
        double _jcdist(double d, double g, double s);
        double _jcvar(double d, double g, double s);
        shared_ptr<DistanceMatrix> _create_distance_matrix(vector<vector<double>> matrix);
        shared_ptr<VectorSiteContainer> sequences = nullptr;
        shared_ptr<VectorSiteContainer> simulated_sequences = nullptr;
        shared_ptr<SubstitutionModel> model = nullptr;
        shared_ptr<GammaDiscreteDistribution> rates = nullptr;
        shared_ptr<DistanceMatrix> distances = nullptr;
        shared_ptr<DistanceMatrix> variances = nullptr;
        shared_ptr<RHomogeneousTreeLikelihood> likelihood = nullptr;
        shared_ptr<HomogeneousSequenceSimulator> simulator = nullptr;
        bool dna{false};
        bool protein{false};
        string _model;
};

#endif /* DISTANCE_H_ */
