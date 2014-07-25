/*
 * Simulate.h
 *
 *  Created on: Jul 23, 2014
 *      Author: kgori
 */

#ifndef SIMULATOR_H_
#define SIMULATOR_H_

#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Phyl/Model/SubstitutionModel.h>
#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Phyl/Simulation/HomogeneousSequenceSimulator.h>

#include <iostream>
#include <map>
#include <memory>
#include <vector>

using namespace std;
using namespace bpp;

class Simulator {
    public :
        Simulator(string model_name);
        void set_model(string model_name);
        void set_alpha(int ncat=4, double alpha=1.0);
        void set_rates(vector<double>, string order="acgt");
        void set_frequencies(vector<double>);
        double get_alpha();
        vector<double> get_rates(string order);
        vector<double> get_frequencies();
        bool is_dna();
        bool is_protein();
        void write_alignment(string filename, string file_format, bool interleaved);
        void set_simulator(string tree);
        vector<pair<string, string>> simulate(unsigned int nsites, string tree);
        vector<pair<string, string>> simulate(unsigned int nsites);
        vector<pair<string, string>> get_sequences();

    private :
        void _set_datatype(string model);
        map<int, double> _vector_to_map(vector<double>);
        void _write_phylip(string filename, bool interleaved);
        void _write_fasta(string filename);
        bool _is_file(string filename);
        bool _is_tree_string(string tree_string);

        shared_ptr<VectorSiteContainer> sequences;
        shared_ptr<SubstitutionModel> model;
        shared_ptr<GammaDiscreteDistribution> rates;
        shared_ptr<HomogeneousSequenceSimulator> simulator = nullptr;
        bool _dna{false};
        bool _protein{false};
        bool _model = false;
        bool _rates = false;
};

#endif /* SIMULATOR_H_ */
