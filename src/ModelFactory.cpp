/*
 * ModelFactory.cpp
 *
 *  Created on: Jul 20, 2014
 *      Author: kgori
 */

#include "ModelFactory.h"
#include <map>
#include <string>

map<string, Model> ModelMap{
    {"JCnuc", Model::JCnuc},
    {"JC69", Model::JCnuc},
    {"K80", Model::K80},
    {"HKY85", Model::HKY85},
    {"TN93", Model::TN93},
    {"GTR", Model::GTR},
    {"T92", Model::T92},
    {"F84", Model::F84},
    {"JTT92", Model::JTT92},
    {"JTT92+F", Model::JTT92},
    {"JCprot", Model::JCprot},
    {"DSO78",  Model::DSO78},
    {"DSO78+F",  Model::DSO78},
    {"WAG01", Model::WAG01},
    {"WAG01+F", Model::WAG01},
    {"LG08", Model::LG08},
    {"LG08+F", Model::LG08},
};

bool hasEnding (std::string const &fullString, std::string const &ending)
{
    if (fullString.length() >= ending.length()) {
        return (0 == fullString.compare (fullString.length() - ending.length(), ending.length(), ending));
    } else {
        return false;
    }
}

Model string_to_model(string name) throw (Exception) {
    if ( ModelMap.find(name) == ModelMap.end() ) {
      throw Exception("ModelFactory::string_to_model() - unknown model \"" + name + "\"");
    } else {
      return ModelMap[name];
    }
}

ModelFactory::ModelFactory() {}

ModelFactory::~ModelFactory() {}

shared_ptr<SubstitutionModel> ModelFactory::create(string model_name) throw (Exception) {
    try {
        Model model = string_to_model(model_name);
        bool pf = hasEnding(model_name, "+F");
        return create(model, pf);
    }
    catch (exception &e){
        cout << e.what()
             << " ---> using JCnuc model due to error finding model \"" << model_name << "\"" << endl;
        Model model = Model::JCnuc;
        return create(model, false);
    }
}

shared_ptr<SubstitutionModel> ModelFactory::create(Model model, bool parameterise_freqs) throw (Exception) {
    switch (model) {
    case Model::JCnuc:
        return make_shared<JCnuc>(&AlphabetTools::DNA_ALPHABET);
        break;
    case Model::K80:
        return make_shared<K80>(&AlphabetTools::DNA_ALPHABET);
        break;
    case Model::HKY85:
        return make_shared<HKY85>(&AlphabetTools::DNA_ALPHABET);
        break;
    case Model::TN93:
        return make_shared<TN93>(&AlphabetTools::DNA_ALPHABET);
        break;
    case Model::GTR:
        return make_shared<GTR>(&AlphabetTools::DNA_ALPHABET, 1.0, 0.1, 0.02, 0.1, 0.02, 0.25, 0.25, 0.25, 0.25);
        break;
    case Model::T92:
        return make_shared<T92>(&AlphabetTools::DNA_ALPHABET);
        break;
    case Model::F84:
        return make_shared<F84>(&AlphabetTools::DNA_ALPHABET);
        break;
    case Model::JTT92:
        if (parameterise_freqs) {
            auto freqs_set = new FullProteinFrequenciesSet(&AlphabetTools::PROTEIN_ALPHABET);
            return make_shared<JTT92>(&AlphabetTools::PROTEIN_ALPHABET, freqs_set, true);
        }
        return make_shared<JTT92>(&AlphabetTools::PROTEIN_ALPHABET);
        break;
    case Model::JCprot:
        if (parameterise_freqs) {
            auto freqs_set = new FullProteinFrequenciesSet(&AlphabetTools::PROTEIN_ALPHABET);
            return make_shared<JCprot>(&AlphabetTools::PROTEIN_ALPHABET, freqs_set, true);
        }
        return make_shared<JCprot>(&AlphabetTools::PROTEIN_ALPHABET);
        break;
    case Model::DSO78:
        if (parameterise_freqs) {
            auto freqs_set = new FullProteinFrequenciesSet(&AlphabetTools::PROTEIN_ALPHABET);
            return make_shared<DSO78>(&AlphabetTools::PROTEIN_ALPHABET, freqs_set, true);
        }
        return make_shared<DSO78>(&AlphabetTools::PROTEIN_ALPHABET);
        break;
    case Model::WAG01:
        if (parameterise_freqs) {
            auto freqs_set = new FullProteinFrequenciesSet(&AlphabetTools::PROTEIN_ALPHABET);
            return make_shared<WAG01>(&AlphabetTools::PROTEIN_ALPHABET, freqs_set, true);
        }
        return make_shared<WAG01>(&AlphabetTools::PROTEIN_ALPHABET);
        break;
    case Model::LG08:
        if (parameterise_freqs) {
            auto freqs_set = new FullProteinFrequenciesSet(&AlphabetTools::PROTEIN_ALPHABET);
            return make_shared<LG08>(&AlphabetTools::PROTEIN_ALPHABET, freqs_set, true);
        }
        return make_shared<LG08>(&AlphabetTools::PROTEIN_ALPHABET);
        break;
    default:
        throw Exception("ModelFactory::create() - unknown model. ");
    }
}

shared_ptr<SubstitutionModel> ModelFactory::create(string model_name, vector<double> freqs) throw (Exception) {
    Model model = string_to_model(model_name);
    FullProteinFrequenciesSet *freqs_set = nullptr;
    switch (model) {
    case Model::JTT92:
        freqs_set = new FullProteinFrequenciesSet(&AlphabetTools::PROTEIN_ALPHABET, freqs);
        return make_shared<JTT92>(&AlphabetTools::PROTEIN_ALPHABET, freqs_set, false);
    case Model::JCprot:
        freqs_set = new FullProteinFrequenciesSet(&AlphabetTools::PROTEIN_ALPHABET, freqs);
        return make_shared<JCprot>(&AlphabetTools::PROTEIN_ALPHABET, freqs_set, false);
    case Model::DSO78:
        freqs_set = new FullProteinFrequenciesSet(&AlphabetTools::PROTEIN_ALPHABET, freqs);
        return make_shared<DSO78>(&AlphabetTools::PROTEIN_ALPHABET, freqs_set, false);
    case Model::WAG01:
        freqs_set = new FullProteinFrequenciesSet(&AlphabetTools::PROTEIN_ALPHABET, freqs);
        return make_shared<WAG01>(&AlphabetTools::PROTEIN_ALPHABET, freqs_set, false);
    case Model::LG08:
        freqs_set = new FullProteinFrequenciesSet(&AlphabetTools::PROTEIN_ALPHABET, freqs);
        return make_shared<LG08>(&AlphabetTools::PROTEIN_ALPHABET, freqs_set, false);
    default:
        throw Exception("ModelFactory::create() - unknown model. ");
    }
}
