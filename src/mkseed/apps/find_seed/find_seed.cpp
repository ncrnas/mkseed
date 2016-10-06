#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/find.h>
#include <seqan/index.h>
#include "mkseed_core.hpp"

using namespace seqan;

int main(int argc, char const ** argv)
{
    // Declare variables
    StringSet<CharString> mirna_ids;
    StringSet<CharString> mrna_ids;
    StringSet<RnaString> mirna_seqs;
    StringSet<RnaString> mrna_seqs;

    // Parse the command line.
    FindSeedOptions options;
    ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
    if (res != seqan::ArgumentParser::PARSE_OK)
    {
        return res == seqan::ArgumentParser::PARSE_ERROR;
    }

    // Read fasta files
    int fread_res;
    fread_res = read_fasta(options.mirna_fa, mirna_ids, mirna_seqs);
    if (fread_res != 0)
    {
        return 1;
    }
    fread_res = read_fasta(options.mrna_fa, mrna_ids, mrna_seqs);
    if (fread_res != 0)
    {
        return 1;
    }

    // Open output file
    std::ofstream ofile(toCString(options.ofile));
    if (!ofile.good())
    {
        std::cerr << "ERROR: Could not open output file " << toCString(options.ofile) << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    // Search seed sites for each miRNA
    int search_res;
    if (options.bruteforce)
    {
        search_res = search_by_bruteforce(ofile, mirna_ids, mirna_seqs, mrna_ids, mrna_seqs);
    }
    else if (options.horspool)
    {
        search_res = search_by_horspool(ofile, mirna_ids, mirna_seqs, mrna_ids, mrna_seqs);
    }
    else if (options.suffixarray)
    {
        search_res = search_by_suffixarray(ofile, mirna_ids, mirna_seqs, mrna_ids, mrna_seqs);
    }
    else
    {
        search_res = search_by_ngram(ofile, mirna_ids, mirna_seqs, mrna_ids, mrna_seqs);
    }

    if (search_res != 0)
    {
        return 1;
    }

    return 0;
}
