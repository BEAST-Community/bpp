/*
 * Simulator.cpp
 *
 *  Created on: Jul 20, 2014
 *      Author: kgori
 */

#include "Simulator.h"
#include "SiteContainerBuilder.h"
#include "ModelFactory.h"

#include <Bpp/Numeric/Prob/GammaDiscreteDistribution.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Seq/SymbolListTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Simulation/HomogeneousSequenceSimulator.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Io/Phylip.h>

#include <fstream>
#include <memory>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <utility>

using namespace bpp;
using namespace std;

Simulator::Simulator(string model_name) {
    set_model(model_name);
}

void Simulator::set_model(string model_name) {
    unique_ptr<ModelFactory> factory(new ModelFactory());
    model = factory->create(model_name);
    _set_datatype(model_name);
    _model = true;
}

bool Simulator::is_dna() {
    return (_dna && !_protein);
}

bool Simulator::is_protein() {
    return (_protein && !_dna);
}

void Simulator::set_alpha(int ncat, double alpha) {
    rates = make_shared<GammaDiscreteDistribution>(ncat, alpha, alpha);
    rates->aliasParameters("alpha", "beta");
    _rates = true;
}

void Simulator::set_rates(vector<double> rates, string order) {
    if (is_dna()) {
        if (order == "acgt" || order == "ACGT") {
            double normaliser = rates[1];
            model->setParameterValue("a", rates[4] / normaliser);
            model->setParameterValue("b", rates[2] / normaliser);
            model->setParameterValue("c", rates[5] / normaliser);
            model->setParameterValue("d", rates[0] / normaliser);
            model->setParameterValue("e", rates[3] / normaliser);
        }
        else if (order == "tcag" || order == "TCAG") {
            model->setParameterValue("a", rates[0]);
            model->setParameterValue("b", rates[1]);
            model->setParameterValue("c", rates[2]);
            model->setParameterValue("d", rates[3]);
            model->setParameterValue("e", rates[4]);
        }
        else {
            throw Exception("Unrecognised order for rates: " + order);
        }
    }
}

void Simulator::set_frequencies(vector<double> freqs) {
    size_t reqd = is_dna() ? 4 : 20;
    if (freqs.size() != reqd) {
        throw Exception("Frequencies vector is the wrong length (dna: 4; aa: 20)");
    }
    map<int, double> m = _vector_to_map(freqs);
    model->setFreq(m);
}

map<int, double> Simulator::_vector_to_map(vector<double> vec) {
    map<int, double> m;
    size_t l = vec.size();
    for (size_t i = 0; i < l; ++i) {
        m[i] = vec[i];
    }
    return m;
}

double Simulator::get_alpha() {
    if (rates) {
        return rates->getParameterValue("alpha");
    }
    return -1;
}

vector<double> Simulator::get_rates(string order) {
    if (is_dna()) {
        vector<double> rates_vec;
        if (order == "acgt" || order == "ACGT") { //{a-c, a-g, a-t, c-g, c-t, g-t=1}
            double normaliser = model->getParameterValue("c");
            rates_vec.push_back(model->getParameterValue("d") / normaliser);
            rates_vec.push_back(1.0 / normaliser);
            rates_vec.push_back(model->getParameterValue("b") / normaliser);
            rates_vec.push_back(model->getParameterValue("e") / normaliser);
            rates_vec.push_back(model->getParameterValue("a") / normaliser);
            rates_vec.push_back(1.0);
        }
        else if (order == "tcag" || order == "TCAG") { //{a=t-c, b=t-a, c=t-g, d=c-a, e=c-g, f=a-g=1}
            rates_vec.push_back(model->getParameterValue("a"));
            rates_vec.push_back(model->getParameterValue("b"));
            rates_vec.push_back(model->getParameterValue("c"));
            rates_vec.push_back(model->getParameterValue("d"));
            rates_vec.push_back(model->getParameterValue("e"));
            rates_vec.push_back(1.0);
        }
        else {
            cerr << "Unknown order: " << order << ". Accepted orders are {tcag, acgt}" << endl;
            throw exception();
        }
        return rates_vec;
    }
    else {
        cerr << "Getting and setting rates is not implemented for protein models" << endl;
        throw exception();
    }
}

vector<double> Simulator::get_frequencies() {
    return model->getFrequencies();
}


void Simulator::_set_datatype(string model) {
    if (model == "JTT92" || model == "JCprot" || model == "DSO78" || model == "WAG01" || model == "LG08") {
        _dna = false;
        _protein = true;
    }
    else if (model == "GTR") {
        _dna = true;
        _protein = false;
    }
    else {
        cerr << "Not sure how to set datatype for unknown model: " << model << endl;
        throw exception();
    }
}


void Simulator::write_alignment(string filename, string file_format, bool interleaved) {
    if (file_format == "fas" || file_format == "fasta") {
        _write_fasta(filename);
    }
    else if (file_format == "phy" || file_format == "phylip") {
        _write_phylip(filename, interleaved);
    }
    else {
        cerr << "Unrecognised file format: " << file_format << endl;
        throw exception();
    }
}


void Simulator::_write_fasta(string filename) {
    Fasta writer;
    writer.writeAlignment(filename, *sequences);
}


void Simulator::_write_phylip(string filename, bool interleaved) {
    Phylip writer{true, !interleaved, 100, true, "  "};
    writer.writeAlignment(filename, *sequences, true);
}


void Simulator::set_simulator(string tree) {
    if (!_model) {
        cerr << "Model not set" << endl;
        throw exception();
    }
    if (!_rates) {
        cerr << "Rates not set" << endl;
        throw exception();
    }
    Tree * simtree;
    Newick * reader = new Newick(false);
    if (_is_file(tree)) {
        simtree = reader->read(tree);
        delete reader;
    }
    else if (_is_tree_string(tree)) {
        stringstream ss{tree};
        simtree = reader->read(ss);
        delete reader;
    }
    else {
        cerr << "Couldn\'t understand this tree: " << tree << endl;
        delete reader;
        throw exception();
    }
    simulator = make_shared<HomogeneousSequenceSimulator>(model.get(), rates.get(), simtree);
    delete simtree;
}

vector<pair<string, string>> Simulator::simulate(unsigned int nsites, string tree) {
    set_simulator(tree);
    return simulate(nsites);
}

vector<pair<string, string>> Simulator::simulate(unsigned int nsites) {
    if (!simulator) {
        cout << "Tried to simulate without a simulator" << endl;
        throw exception();
    }
    size_t nsites_{nsites};
    SiteContainer * tmp = simulator->simulate(nsites_);
    sequences = make_shared<VectorSiteContainer>(*tmp);
    delete tmp;
    return get_sequences();
}

bool Simulator::_is_file(string filename) {
    ifstream fl(filename.c_str());
    bool result = true;
    if (!fl) {
        result = false;
    }
    fl.close();
    return result;
}

bool Simulator::_is_tree_string(string tree_string) {
    size_t l = tree_string.length();
    return (tree_string[0]=='(' && tree_string[l-1]==';');
}

vector<pair<string, string>> Simulator::get_sequences() {
    vector<pair<string, string>> ret;
    if (!sequences) {
        cerr << "No sequences exist yet" << endl;
        throw exception();
    }
    for (size_t i = 0; i < sequences->getNumberOfSequences(); ++i) {
        BasicSequence seq = sequences->getSequence(i);
        ret.push_back(make_pair(seq.getName(), seq.toString()));
    }
    return ret;
}
