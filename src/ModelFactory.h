/*
 * ModelFactory.h
 *
 *  Created on: Jul 20, 2014
 *      Author: kgori
 */

#ifndef MODELFACTORY_H_
#define MODELFACTORY_H_

#include "Bpp/Seq/Alphabet/Alphabet.h"
#include "Bpp/Seq/Alphabet/AlphabetTools.h"
#include "Bpp/Seq/Alphabet/AlphabetExceptions.h"
#include "Bpp/Seq/GeneticCode/GeneticCode.h"
#include "Bpp/Phyl/Model/SubstitutionModel.h"
#include <Bpp/Phyl/Model/Nucleotide/JCnuc.h>
#include <Bpp/Phyl/Model/Nucleotide/K80.h>
#include <Bpp/Phyl/Model/Nucleotide/HKY85.h>
#include <Bpp/Phyl/Model/Nucleotide/TN93.h>
#include <Bpp/Phyl/Model/Nucleotide/GTR.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Phyl/Model/Nucleotide/F84.h>
#include <Bpp/Phyl/Model/Protein/JTT92.h>
#include <Bpp/Phyl/Model/Protein/JCprot.h>
#include <Bpp/Phyl/Model/Protein/DSO78.h>
#include <Bpp/Phyl/Model/Protein/WAG01.h>
#include <Bpp/Phyl/Model/Protein/LG08.h>
#include <Bpp/Phyl/Model/FrequenciesSet/ProteinFrequenciesSet.h>
#include <memory>
#include <map>
#include <string>

using namespace bpp;
using namespace std;


enum class Model {
    JCnuc,
    K80,
    HKY85,
    TN93,
    GTR,
    T92,
    F84,
    JTT92,
    JCprot,
    DSO78,
    WAG01,
    LG08,
};

Model string_to_model(string name) throw (Exception);

class ModelFactory {

public:
    ModelFactory();
    virtual ~ModelFactory();
    static shared_ptr<SubstitutionModel> create(string model_name) throw (Exception);
    static shared_ptr<SubstitutionModel> create(Model model, bool parameterise_freqs) throw (Exception);
    static shared_ptr<SubstitutionModel> create(string model_name, vector<double> freqs) throw (Exception);
};

#endif /* MODELFACTORY_H_ */
