/*
 * Distance.cpp
 *
 *  Created on: Jul 20, 2014
 *      Author: kgori
 */

#include "Distance.h"
#include "SiteContainerBuilder.h"
#include "ModelFactory.h"

#include "Bpp/Numeric/Prob/GammaDiscreteDistribution.h"
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include "Bpp/Seq/SymbolListTools.h"
#include "Bpp/Phyl/Distance/DistanceEstimation.h"
#include <Bpp/Phyl/Io/Newick.h>
#include "Bpp/Phyl/OptimizationTools.h"

#include <string>
#include <sstream>
#include <vector>
#include <map>

using namespace bpp;
using namespace std;

Distance::Distance(string filename, string file_format, string datatype, string model_name, bool interleaved) {
    _check_compatible_model(datatype, model_name);
    read_alignment(filename, file_format, datatype, interleaved);
    set_model(model_name);
    set_alpha();
}

void Distance::read_alignment(string filename, string file_format, string datatype, bool interleaved) {
    sequences = SiteContainerBuilder::read_alignment(filename, file_format, datatype, interleaved);
    string type = sequences->getAlphabet()->getAlphabetType();
    if (type == "DNA alphabet") {
        _set_dna();
    }
    else if (type == "Proteic alphabet") {
        _set_protein();
    }
    else {
        cout << "Type = " << type << endl;
    }
    _clear_distances();
    _clear_likelihood();
}

void Distance::set_model(string model_name) {
    unique_ptr<ModelFactory> factory(new ModelFactory());
    model = factory->create(model_name);
    _clear_distances();
    _clear_likelihood();
}

bool Distance::is_dna() {
    return dna && !protein;
}

bool Distance::is_protein() {
    return protein && !dna;
}

void Distance::set_alpha(int ncat, double alpha) {
    rates = make_shared<GammaDiscreteDistribution>(ncat, alpha, alpha);
    rates->aliasParameters("alpha", "beta");
    _clear_distances();
    _clear_likelihood();
}

void Distance::set_rates(vector<double> rates, string order) {
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
        _clear_distances();
        _clear_likelihood();
    }
}

void Distance::_clear_distances() {
    if (distances) {
        distances.reset();
    }
}

void Distance::_clear_likelihood() {
    if (likelihood) {
        likelihood.reset();
    }
}

void Distance::set_frequencies(vector<double> freqs) {
    size_t reqd = is_dna() ? 4 : 20;
    if (freqs.size() != reqd) {
        throw Exception("Frequencies vector is the wrong length (dna: 4; aa: 20)");
    }
    map<int, double> m = _vector_to_map(freqs);
    model->setFreq(m);
    _clear_distances();
    _clear_likelihood();
}

void Distance::compute_distances() {
    VectorSiteContainer* sites_ = sequences->clone();
    SiteContainerTools::changeGapsToUnknownCharacters(*sites_);
    size_t n = sites_->getNumberOfSequences();
    vector<string> names = get_names();

    _clear_distances();
    distances = make_shared<DistanceMatrix>(names);
    for (size_t i = 0; i < n; i++) {
        (*distances)(i, i) = 0;
        for (size_t j = i + 1; j < n; j++) {
            TwoTreeLikelihood* lik =
                new TwoTreeLikelihood(names[i], names[j], *sites_, model.get(), rates.get(), false);
            lik->initialize();
            lik->enableDerivatives(true);
            size_t d = SymbolListTools::getNumberOfDistinctPositions(sites_->getSequence(i), sites_->getSequence(j));
            size_t g = SymbolListTools::getNumberOfPositionsWithoutGap(sites_->getSequence(i), sites_->getSequence(j));
            lik->setParameterValue("BrLen", g == 0 ? lik->getMinimumBranchLength() : std::max(lik->getMinimumBranchLength(), static_cast<double>(d) / static_cast<double>(g)));
            // Optimization:
            ParameterList params = lik->getBranchLengthsParameters();
            OptimizationTools::optimizeNumericalParameters(lik, params, 0, 1, 0.000001, 1000000, NULL, NULL, false, 0, OptimizationTools::OPTIMIZATION_NEWTON, OptimizationTools::OPTIMIZATION_BRENT);
            // Store results:
            (*distances)(i, j) = lik->getParameterValue("BrLen");
            (*distances)(j, i) = 1.0 / lik->d2f("BrLen", params);
            //cout << "d g i j dist var : " << d << " " << g << " " << i << " " << j << " " << (*distances)(i, j) << " " << (*distances)(j, i) << endl;
            delete lik;
        }
    }
    delete sites_;
}

void Distance::_check_distances_exist() {
    if (!distances) {
        compute_distances();
    }
}

vector<vector<double>> Distance::get_distances() {
    _check_distances_exist();
    vector<vector<double>> vec;
    size_t nrow = distances->getNumberOfRows();
    for (size_t i = 0; i < nrow; ++i) {
        vec.push_back(distances->row(i));
    }
    return vec;
}

vector<string> Distance::get_names() {
    return sequences->getSequencesNames();
}

map<int, double> Distance::_vector_to_map(vector<double> vec) {
    map<int, double> m;
    size_t l = vec.size();
    for (size_t i = 0; i < l; ++i) {
        m[i] = vec[i];
    }
    return m;
}

void Distance::_set_dna() {
    dna = true;
    protein = false;
}

void Distance::_set_protein() {
    dna = false;
    protein = true;
}

double Distance::get_alpha() {
    if (rates) {
        return rates->getParameterValue("alpha");
    }
    return -1;
}

vector<double> Distance::get_rates(string order) {
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

vector<double> Distance::get_frequencies() {
    return model->getFrequencies();
}

void Distance::_check_compatible_model(string datatype, string model) {
    bool incompat = false;
    if (datatype == "dna" and (model == "JTT92" || model == "JCprot" || model == "DSO78" || model == "WAG01" || model == "LG08")) {
        incompat = true;
    }
    else if (datatype == "aa" and (model == "JCnuc" || model == "JC69" || model == "K80" || model == "HKY85" || model == "TN93"  || model == "GTR" || model == "T92" || model == "F84")) {
        incompat = true;
    }
    if (incompat) {
        cerr << "Incompatible model (" << model << ") and datatype (" << datatype << ")" << endl;
        throw exception();
    }
}

double Distance::get_likelihood() {
    if (!likelihood) {
        cerr << "Likelihood calculator not set - call initialise_likelihood" << endl;
        throw exception();
    }
    return likelihood->getLogLikelihood();
}

bool Distance::_is_file(string filename) {
    ifstream fl(filename.c_str());
    bool result = true;
    if (!fl) {
        result = false;
    }
    fl.close();
    return result;
}

bool Distance::_is_tree_string(string tree_string) {
    size_t l = tree_string.length();
    return (tree_string[0]=='(' && tree_string[l-1]==';');
}

void Distance::initialise_likelihood(string tree) {
    if (!model) {
        cerr << "Model not set" << endl;
        throw exception();
    }
    if (!rates) {
        cerr << "Rates not set" << endl;
        throw exception();
    }
    Tree * liktree;
    auto reader = make_shared<Newick>(false);
    if (_is_file(tree)) {
        liktree = reader->read(tree);
    }
    else if (_is_tree_string(tree)) {
        stringstream ss{tree};
        liktree = reader->read(ss);
    }
    else {
        cerr << "Couldn\'t understand this tree: " << tree << endl;
        throw exception();
    }
    VectorSiteContainer* sites_ = sequences->clone();
    SiteContainerTools::changeGapsToUnknownCharacters(*sites_);
    likelihood = make_shared<RHomogeneousTreeLikelihood>(*liktree, *sites_, model.get(), rates.get(), true, false, true);
    likelihood->initialize();
    delete liktree;
}

void Distance::optimise_parameters(bool fix_branch_lengths) {
    if (!likelihood) {
        cerr << "Likelihood calculator not set - call initialise_likelihood" << endl;
        throw exception();
    }
    ParameterList pl;
    if (fix_branch_lengths) {
        pl = likelihood->getSubstitutionModelParameters();
        pl.addParameters(likelihood->getRateDistributionParameters());
    }
    else {
        pl = likelihood->getParameters();
    }
    OptimizationTools::optimizeNumericalParameters2(likelihood.get(), pl, 0, 0.0001, 1000000, NULL, NULL, false, false, 0);
}

string Distance::get_tree() {
    if (!likelihood) {
        cerr << "Likelihood calculator not set - call initialise_likelihood" << endl;
        throw exception();
    }
    auto *tree = likelihood->getTree().clone();
    stringstream ss;
    Newick treeWriter;
    treeWriter.write(*tree, ss);
    delete tree;
    return ss.str();
}
