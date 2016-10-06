#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/find.h>
#include <seqan/index.h>

struct FindSeedOptions
{
    bool canonical;
    bool bruteforce;
    bool suffixarray;
    bool horspool;

    seqan::CharString mirna_fa;
    seqan::CharString mrna_fa;
    seqan::CharString ofile;

    FindSeedOptions() : canonical(false), bruteforce(false), suffixarray(false), horspool(false) {}
};

seqan::ArgumentParser::ParseResult parseCommandLine(
        FindSeedOptions &options,
        int const argc,
        char const **argv);

int read_fasta(
        seqan::CharString const &fasta,
        seqan::StringSet<seqan::CharString> &ids,
        seqan::StringSet<seqan::RnaString> &seqs);

int get_seed(seqan::RnaString const &mirna, seqan::RnaString &seed);

int search_by_bruteforce(
        std::ofstream &ofile,
        seqan::StringSet<seqan::CharString> const &mirna_ids,
        seqan::StringSet<seqan::RnaString> const &mirnas,
        seqan::StringSet<seqan::CharString> const &mrna_ids,
        seqan::StringSet<seqan::RnaString> &mrnas);

int search_by_horspool(
        std::ofstream &ofile,
        seqan::StringSet<seqan::CharString> const &mirna_ids,
        seqan::StringSet<seqan::RnaString> const &mirnas,
        seqan::StringSet<seqan::CharString> const &mrna_ids,
        seqan::StringSet<seqan::RnaString> &mrnas);

int search_by_suffixarray(
        std::ofstream &ofile,
        seqan::StringSet<seqan::CharString> const &mirna_ids,
        seqan::StringSet<seqan::RnaString> const &mirnas,
        seqan::StringSet<seqan::CharString> const &mrna_ids,
        seqan::StringSet<seqan::RnaString> &mrnas);

int search_by_ngram(
        std::ofstream &ofile,
        seqan::StringSet<seqan::CharString> const &mirna_ids,
        seqan::StringSet<seqan::RnaString> const &mirnas,
        seqan::StringSet<seqan::CharString> const &mrna_ids,
        seqan::StringSet<seqan::RnaString> &mrnas);

