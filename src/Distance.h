/*
 * Distance.h
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

#include <iostream>
#include <map>
#include <memory>
#include <vector>

using namespace std;
using namespace bpp;

class Distance {
    public :
        Distance(string filename, string file_format, string datatype, string model_name, bool interleaved=true);
        void read_alignment(string filename, string file_format, string datatype, bool interleaved=true);
        void set_model(string model_name);
        void set_alpha(int ncat=4, double alpha=1.0);
        void set_rates(vector<double>, string order="acgt");
        void set_frequencies(vector<double>);
        double get_alpha();
        vector<double> get_rates(string order);
        vector<double> get_frequencies();
        vector<vector<double>> get_distances();
        vector<string> get_names();
        void initialise_likelihood(string tree);
        void optimise_parameters(bool fix_branch_lengths);
        double get_likelihood();
        void compute_distances();
        bool is_dna();
        bool is_protein();
        string get_tree();

    private :
        void _set_dna();
        void _set_protein();
        map<int, double> _vector_to_map(vector<double>);
        void _check_distances_exist();
        void _check_compatible_model(string datatype, string model);
        void _clear_distances();
        void _clear_likelihood();
        bool _is_file(string filename);
        bool _is_tree_string(string tree_string);

        shared_ptr<VectorSiteContainer> sequences = nullptr;
        shared_ptr<SubstitutionModel> model = nullptr;
        shared_ptr<GammaDiscreteDistribution> rates = nullptr;
        shared_ptr<DistanceMatrix> distances = nullptr;
        shared_ptr<RHomogeneousTreeLikelihood> likelihood = nullptr;
        bool dna{false};
        bool protein{false};
        bool _model = false;
        bool _rates = false;
        bool _distances{false};
        bool _likelihood{false};
};

#endif /* DISTANCE_H_ */
