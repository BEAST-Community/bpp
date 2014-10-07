/*
 * Alignment.h
 *
 *  Created on: Jul 20, 2014
 *      Author: kgori
 */

#ifndef _ALIGNMENT_H_
#define _ALIGNMENT_H_

#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Numeric/Prob/AbstractDiscreteDistribution.h>
#include <Bpp/Seq/DistanceMatrix.h>
#include <Bpp/Phyl/Likelihood/NNIHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Simulation/HomogeneousSequenceSimulator.h>

#include <iostream>
#include <map>
#include <memory>
#include <vector>

using namespace std;
using namespace bpp;

class Alignment {
    public :
        Alignment();
        Alignment(vector<Alignment> alignments);
        Alignment(vector<pair<string, string>> headers_sequences, string datatype);
        Alignment(string filename, string file_format, bool interleaved=true);
        Alignment(string filename, string file_format, string datatype, bool interleaved=true);
        Alignment(string filename, string file_format, string datatype, string model_name, bool interleaved=true);
        void read_alignment(string filename, string file_format, bool interleaved);
        void read_alignment(string filename, string file_format, string datatype, bool interleaved=true);
        void sort_alignment(bool ascending=true);
        void write_alignment(string filename, string file_format, bool interleaved=true);
        void set_substitution_model(string model_name);
        void set_gamma_rate_model(size_t ncat=4, double alpha=1.0);
        void set_constant_rate_model();
        void set_alpha(double alpha);
        void set_number_of_gamma_categories(size_t ncat);
        void set_rates(vector<double>, string order="acgt");
        void set_frequencies(vector<double>);
        void set_namespace(string name);
        vector<pair<string, string>> get_sequences();
        double get_alpha();
        size_t get_number_of_gamma_categories();
        vector<double> get_rates(string order);
        vector<double> get_rate_model_categories();
        vector<double> get_frequencies();
        vector<string> get_names();
        size_t get_number_of_sequences();
        vector<string> get_sites();
        size_t get_number_of_sites();
        size_t get_number_of_free_parameters();
        string get_substitution_model();
        vector<vector<double>> get_exchangeabilities();
        string get_namespace();
        vector<string> get_informative_sites(bool exclude_gaps);
        size_t get_number_of_informative_sites(bool exclude_gaps);
        size_t get_number_of_distinct_sites();
        bool is_dna();
        bool is_protein();
        void _print_params();

        // Distance
        void compute_distances();
        void fast_compute_distances();
        void set_distance_matrix(vector<vector<double>> matrix);
        string get_bionj_tree();
        string get_bionj_tree(vector<vector<double>> matrix);
        vector<vector<double>> get_distances();
        vector<vector<double>> get_variances();
        vector<vector<double>> get_distance_variance_matrix();

        // Likelihood
        void initialise_likelihood();
        void initialise_likelihood(string tree);
        void optimise_parameters(bool fix_branch_lengths);
        void optimise_topology(bool fix_model_params);
        double get_likelihood();
        string get_tree();

        // Simulator
        void write_simulation(size_t nsites, string filename, string file_format, bool interleaved=true);
        void set_simulator(string tree);
        vector<pair<string, string>> simulate(size_t nsites, string tree);
        vector<pair<string, string>> simulate(size_t nsites);
        vector<pair<string, string>> get_simulated_sequences();

        // Bootstrap
        vector<pair<string, string>> get_bootstrapped_sequences();
        void chkdst();

        // Misc
        static string get_mrp_supertree(vector<string> trees);

    private :
        string _get_datatype();
        vector<pair<string, string>> _get_sequences(VectorSiteContainer *seqs);
        void _write_fasta(shared_ptr<VectorSiteContainer> seqs, string filename);
        void _write_phylip(shared_ptr<VectorSiteContainer> seqs, string filename, bool interleaved=true);
        map<int, double> _vector_to_map(vector<double>);
        void _check_compatible_model(string model);
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
        shared_ptr<AbstractDiscreteDistribution> rates = nullptr;
        shared_ptr<DistanceMatrix> distances = nullptr;
        shared_ptr<DistanceMatrix> variances = nullptr;
        shared_ptr<NNIHomogeneousTreeLikelihood> likelihood = nullptr;
        shared_ptr<HomogeneousSequenceSimulator> simulator = nullptr;
        string _name;
        string _computeTree(DistanceMatrix dists, DistanceMatrix vars) throw (Exception);
};

#endif /* _ALIGNMENT_H_ */
