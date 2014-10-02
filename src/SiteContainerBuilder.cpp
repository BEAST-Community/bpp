/*
 * SiteContainerBuilder.cpp
 *
 *  Created on: Jul 20, 2014
 *      Author: kgori
 */

#include "SiteContainerBuilder.h"
#include <Bpp/Exceptions.h>
#include <Bpp/Seq/Alphabet/AlphabetExceptions.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/LetterAlphabet.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Io/Phylip.h>
#include <Bpp/Seq/Sequence.h>
#include <algorithm>
#include <functional>
#include <numeric>
#include <stdexcept>

/*
Takes two vectors and returns the union of their elements, in sorted order
*/
template<typename T>
vector<T> merge_vectors(vector<T> &v1, vector<T> &v2) {
    sort(v1.begin(), v1.end());
    sort(v2.begin(), v2.end());
    vector<T> result;
    set_union(v1.begin(), v1.end(), v2.begin(), v2.end(),
              back_inserter(result));
    return result;
}

/*
Merges two VectorSiteContainers into a single VectorSiteContainer
*/
shared_ptr<VectorSiteContainer> merge_vscs(shared_ptr<VectorSiteContainer> s1, shared_ptr<VectorSiteContainer> s2)
throw (AlphabetMismatchException, Exception)
{
    if (s1->getAlphabet()->getAlphabetType() != s2->getAlphabet()->getAlphabetType())
        throw AlphabetMismatchException("SiteContainerTools::merge.", s1->getAlphabet(), s2->getAlphabet());

    auto alphabet = s1->getAlphabet();
    const int unknown_char = alphabet->getUnknownCharacterCode();
    size_t l1 = s1->getNumberOfSites();
    size_t l2 = s2->getNumberOfSites();
    auto n1 = s1->getSequencesNames();
    auto n2 = s2->getSequencesNames();
    vector<string> allnames = merge_vectors(n1, n2);

    unique_ptr<VectorSequenceContainer> tmp_container(new VectorSequenceContainer(alphabet));
    for (auto &n : allnames) {
        vector<int> seqvec1;
        vector<int> seqvec2;
        seqvec1.reserve(l1+l2);
        seqvec2.reserve(l2);

        if (s1->hasSequence(n)) {
            auto tmpvec = s1->getSequence(n).getContent();
            seqvec1.insert( seqvec1.end(), tmpvec.begin(), tmpvec.end());
        }
        else {
            vector<int> tmpvec(l1, unknown_char);
            seqvec1.insert( seqvec1.end(), tmpvec.begin(), tmpvec.end());
        }

        if (s2->hasSequence(n)) {
            auto tmpvec = s2->getSequence(n).getContent();
            seqvec2.insert( seqvec2.end(), tmpvec.begin(), tmpvec.end());
        }
        else {
            vector<int> tmpvec(l2, unknown_char);
            seqvec2.insert( seqvec2.end(), tmpvec.begin(), tmpvec.end());
        }
        seqvec1.insert(seqvec1.end(), seqvec2.begin(), seqvec2.end());
        unique_ptr<BasicSequence> seq(new BasicSequence(n, seqvec1, alphabet));

        tmp_container->addSequence(*seq, false);
    }

    auto output = make_shared<VectorSiteContainer>(*tmp_container);
    return output;
}

/*
Folds vector of VSC into single VSC using merge()
*/
//shared_ptr<VectorSiteContainer> fold(vector<shared_ptr<VectorSiteContainer>> vec_of_vsc) {
//    auto first_vsc = make_shared<VectorSiteContainer>(*vec_of_vsc[0]);
//    //return accumulate(vec_of_vsc.begin()+1, vec_of_vsc.end(), first_vsc, merge_vscs);
//    return merge_vscs(first_vsc, vec_of_vsc[1]);
//}

shared_ptr<VectorSiteContainer> SiteContainerBuilder::construct_alignment_from_strings(vector<pair<string, string>> headers_sequences, string datatype)
        throw (Exception) {
    if (asking_for_dna(datatype)) {
        return make_shared<VectorSiteContainer>(*_convert_vector_to_dna_sequence_container(headers_sequences));
    }
    else if (asking_for_protein(datatype)) {
        return make_shared<VectorSiteContainer>(*_convert_vector_to_protein_sequence_container(headers_sequences));
    }
    else {
        throw Exception(datatype);
    }
}

shared_ptr<VectorSiteContainer> SiteContainerBuilder::read_alignment(string filename,
        string file_format, bool interleaved)
                throw (Exception) {
    try {
        return read_alignment(filename, file_format, "nt", interleaved);
    }
    catch (AlphabetException &e) {
        return read_alignment(filename, file_format, "aa", interleaved);
    }
}


shared_ptr<VectorSiteContainer> SiteContainerBuilder::read_alignment(string filename,
        string file_format, string datatype, bool interleaved)
                throw (Exception) {
    if (asking_for_fasta(file_format)) {
        if (asking_for_dna(datatype)) {
            return read_fasta_dna_file(filename);
        }
        else if (asking_for_protein(datatype)) {
            return read_fasta_protein_file(filename);
        }
        else {
            throw Exception(datatype);
        }
    }

    else if (asking_for_phylip(file_format)) {
        if (asking_for_dna(datatype)) {
            return read_phylip_dna_file(filename, interleaved);

        }
        else if (asking_for_protein(datatype)) {
            return read_phylip_protein_file(filename, interleaved);
        }
        else {
            throw Exception(datatype);
        }
    }
    else {
        throw Exception(file_format);
    }
}

shared_ptr<VectorSiteContainer> SiteContainerBuilder::construct_sorted_alignment(VectorSiteContainer *sites,
        bool ascending) {
    VectorSequenceContainer *tmp = new VectorSequenceContainer(sites->getAlphabet());
    vector<string> names = sites->getSequencesNames();
    if (ascending) sort(names.begin(), names.end());
    else sort(names.begin(), names.end(), greater<string>());
    for (string &name : names) {
        tmp->addSequence(sites->getSequence(name));
    }
    auto ret = make_shared<VectorSiteContainer>(*tmp);
    delete tmp;
    return ret;
}

/*
Folds vector of VSC into single VSC using merge_vscs()
*/
shared_ptr<VectorSiteContainer> SiteContainerBuilder::concatenate_alignments(vector<shared_ptr<VectorSiteContainer>> vec_of_vsc) {
    auto first_vsc = make_shared<VectorSiteContainer>(*vec_of_vsc[0]);
    return accumulate(vec_of_vsc.begin()+1, vec_of_vsc.end(), first_vsc, merge_vscs);
//    return merge_vscs(first_vsc, vec_of_vsc[1]);
}

bool SiteContainerBuilder::asking_for_fasta(string file_format) {
    return (file_format == "fasta" || file_format == "fas" || file_format == ".fasta" || file_format == ".fas");
}

bool SiteContainerBuilder::asking_for_phylip(string file_format) {
    return (file_format == "phylip" || file_format == "phy" || file_format == ".phylip" || file_format == ".phy");
}

bool SiteContainerBuilder::asking_for_dna(string datatype) {
    return (datatype == "dna" || datatype == "nucleotide" || datatype == "nt");
}

bool SiteContainerBuilder::asking_for_protein(string datatype) {
    return (datatype == "protein" || datatype == "aminoacid" || datatype == "aa");
}

shared_ptr<VectorSiteContainer> SiteContainerBuilder::read_fasta_dna_file(
        string filename) {
    Fasta reader;
    SequenceContainer* alignment = reader.readSequences(filename, &AlphabetTools::DNA_ALPHABET);
    shared_ptr<VectorSiteContainer> sequences(new VectorSiteContainer(*alignment));
    delete alignment;
    if (sequences->getNumberOfSequences() == 0) {
        sequences.reset();
        throw Exception("The alignment is empty - did you specify the right file format?");
    }
    return sequences;
}

shared_ptr<VectorSiteContainer> SiteContainerBuilder::read_fasta_protein_file(
        string filename) {
    Fasta reader;
    SequenceContainer* alignment = reader.readSequences(filename, &AlphabetTools::PROTEIN_ALPHABET);
    shared_ptr<VectorSiteContainer> sequences(new VectorSiteContainer(*alignment));
    delete alignment;
    if (sequences->getNumberOfSequences() == 0) {
        sequences.reset();
        throw Exception("The alignment is empty - did you specify the right file format?");
    }
    return sequences;
}

shared_ptr<VectorSiteContainer> SiteContainerBuilder::read_phylip_dna_file(
        string filename, bool interleaved) {
    Phylip reader{true, !interleaved, 100, true, "  "};
    SequenceContainer* alignment = reader.readSequences(filename, &AlphabetTools::DNA_ALPHABET);
    shared_ptr<VectorSiteContainer> sequences(new VectorSiteContainer(*alignment));
    delete alignment;
    if (sequences->getNumberOfSequences() == 0) {
        sequences.reset();
        throw Exception("The alignment is empty - did you specify the right file format?");
    }
    return sequences;
}

shared_ptr<VectorSiteContainer> SiteContainerBuilder::read_phylip_protein_file(
        string filename, bool interleaved) {
    //const ProteicAlphabet *alpha = new ProteicAlphabet();
    Phylip reader{true, !interleaved, 100, true, "  "};
    SequenceContainer* alignment = reader.readSequences(filename, &AlphabetTools::PROTEIN_ALPHABET);
    shared_ptr<VectorSiteContainer> sequences(new VectorSiteContainer(*alignment));
    delete alignment;
    if (sequences->getNumberOfSequences() == 0) {
        sequences.reset();
        throw Exception("The alignment is empty - did you specify the right file format?");
    }
    return sequences;
}

shared_ptr<BasicSequence> SiteContainerBuilder::_convert_pair_to_dna_sequence(pair<string, string> h_s) {
    return make_shared<BasicSequence>(h_s.first, h_s.second, &AlphabetTools::DNA_ALPHABET);
}

shared_ptr<BasicSequence> SiteContainerBuilder::_convert_pair_to_protein_sequence(pair<string, string> h_s) {
    return make_shared<BasicSequence>(h_s.first, h_s.second, &AlphabetTools::PROTEIN_ALPHABET);
}

shared_ptr<VectorSequenceContainer> SiteContainerBuilder::_convert_vector_to_dna_sequence_container(vector<pair<string, string>> v) {
    auto container = make_shared<VectorSequenceContainer>(&AlphabetTools::DNA_ALPHABET);
    for (auto item : v) {
        container->addSequence(*_convert_pair_to_dna_sequence(item));
    }
    return container;
}

shared_ptr<VectorSequenceContainer> SiteContainerBuilder::_convert_vector_to_protein_sequence_container(vector<pair<string, string>> v) {
    auto container = make_shared<VectorSequenceContainer>(&AlphabetTools::PROTEIN_ALPHABET);
    for (auto item : v) {
        container->addSequence(*_convert_pair_to_protein_sequence(item));
    }
    return container;
}


