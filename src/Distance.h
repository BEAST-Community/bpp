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
//        double get_alpha();
//        vector<double> get_rates();
//        vector<double> get_frequencies();
        vector<vector<double>> get_distances();
        vector<string> get_names();
        void compute_distances();
        bool is_dna();
        bool is_protein();

    private :
        void _set_dna();
        void _set_protein();
        map<int, double> _vector_to_map(vector<double>);
        void _check_distances_exist();
        void _check_compatible_model(string datatype, string model);
	void _clear_distances();

        shared_ptr<VectorSiteContainer> sequences;
        shared_ptr<SubstitutionModel> model;
        shared_ptr<GammaDiscreteDistribution> rates;
        shared_ptr<DistanceMatrix> distances;
        bool dna{false};
        bool protein{false};
        bool _distances{false};
};

#endif /* DISTANCE_H_ */
