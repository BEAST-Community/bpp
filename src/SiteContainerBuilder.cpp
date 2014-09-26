/*
 * SiteContainerBuilder.cpp
 *
 *  Created on: Jul 20, 2014
 *      Author: kgori
 */

#include "SiteContainerBuilder.h"
#include <Bpp/Exceptions.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Alphabet/LetterAlphabet.h>
#include <Bpp/Seq/Io/Fasta.h>
#include <Bpp/Seq/Io/Phylip.h>
#include <Bpp/Seq/Sequence.h>
#include <stdexcept>

shared_ptr<VectorSiteContainer> SiteContainerBuilder::construct_alignment_from_strings(vector<pair<string, string> headers_sequences, 
    string datatype) 
        throw (Exception) {
    vector<Sequence*> v;
    if asking_for_dna(datatype) {
        for (pair &h_s : headers_sequences) {
            v.push_back(_convert_pair_to_dna_sequence(h_s).get());
        }
    }
    else if asking_for_protein(datatype) {
        for (pair &h_s : headers_sequences) {
            v.push_back(_convert_pair_to_protein_sequence(h_s).get());
        }
    }
    else {
        throw Exception(datatype);
    }
    auto tmp_seq_container = make_shared<VectorSequenceContainer>(v);
    auto site_container = make_shared<VectorSiteContainer>(*tmp_seq_container);
    return site_container;
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
    return sequences;
}

shared_ptr<VectorSiteContainer> SiteContainerBuilder::read_fasta_protein_file(
        string filename) {
    Fasta reader;
    SequenceContainer* alignment = reader.readSequences(filename, &AlphabetTools::PROTEIN_ALPHABET);
    shared_ptr<VectorSiteContainer> sequences(new VectorSiteContainer(*alignment));
    delete alignment;
    return sequences;
}

shared_ptr<VectorSiteContainer> SiteContainerBuilder::read_phylip_dna_file(
        string filename, bool interleaved) {
    Phylip reader{true, !interleaved, 100, true, "  "};
    SequenceContainer* alignment = reader.readSequences(filename, &AlphabetTools::DNA_ALPHABET);
    shared_ptr<VectorSiteContainer> sequences(new VectorSiteContainer(*alignment));
    delete alignment;
    return sequences;
}

shared_ptr<VectorSiteContainer> SiteContainerBuilder::read_phylip_protein_file(
        string filename, bool interleaved) {
    //const ProteicAlphabet *alpha = new ProteicAlphabet();
    Phylip reader{true, !interleaved, 100, true, "  "};
    SequenceContainer* alignment = reader.readSequences(filename, &AlphabetTools::PROTEIN_ALPHABET);
    shared_ptr<VectorSiteContainer> sequences(new VectorSiteContainer(*alignment));
    delete alignment;
    return sequences;
}

shared_ptr<BasicSequence> SiteContainerBuilder::_convert_pair_to_dna_sequence(pair<string, string>) {
    return make_shared<BasicSequence>(pair.first, pair.second, &AlphabetTools::DNA_ALPHABET);
}
shared_ptr<BasicSequence> SiteContainerBuilder::_convert_pair_to_protein_sequence() {
    return make_shared<BasicSequence>(pair.first, pair.second, &AlphabetTools::DNA_ALPHABET);
}
