#include <iostream>
#include <seqan/arg_parse.h>
#include <seqan/seq_io.h>
#include <seqan/find.h>
#include <seqan/index.h>

using namespace seqan;

struct FindSeedOptions
{
    bool canonical;
    bool bruteforce;
    bool suffixarray;
    bool horspool;

    CharString mirna_fa;
    CharString mrna_fa;
    CharString ofile;

    FindSeedOptions() : canonical(false), bruteforce(false), suffixarray(false), horspool(false) {}
};


ArgumentParser::ParseResult parseCommandLine(
        FindSeedOptions &options,
        int const argc,
        char const **argv)
{
    // Setup ArgumentParser
    ArgumentParser parser("fine_seed");

    // Set short description, version, and date
    setShortDescription(parser, "Find miRNA seed sites in genomic regions");
    setVersion(parser, "1.0");
    setDate(parser, "January 2014");

    // Define usage line and long description
    addUsageLine(parser,
            "[\\fIOPTIONS\\fP] \"\\fIMIRNA FILE\\fP\" \"\\fIMRNA FILE\\fP\" \"\\fIOUTPUT FILE\\fP\"");
    addDescription(parser,
            "This program finds canonical and non-canonical seed sites of miRNAs "
            "in given genomic regions.");

    // Define Arguments
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::INPUTFILE));
    addArgument(parser, ArgParseArgument(ArgParseArgument::OUTPUTFILE));

    // Define Options
    addSection(parser, "Find Seed Options");
    addOption(parser, ArgParseOption("c", "canonical", "Select canonical seeds only."));
    addOption(parser, ArgParseOption("b", "bruteforce", "Use brute-force."));
    addOption(parser, ArgParseOption("s", "suffixarray", "Use suffix-array."));
    addOption(parser, ArgParseOption("l", "horspool", "Use Horspool."));

    // Add Examples Section
    addTextSection(parser, "Examples");
    addListItem(parser,
            "\\fBfind_seed\\fP \\fImirna_fasta\\fP \\fImrna_fasta\\fP \\fIoutput_file\\fP",
            "find both canonical and non-canonical seed sites of \\fImiRNAs\\fP  "
            "in \\fImRNA\\fP regions and write the results to \\fIoutput\\fP.");
    addListItem(parser,
            "\\fBfind_seed\\fP \\fB-c\\fP \\fImirna_fasta\\fP \\fImrna_fasta\\fP \\fIoutput\\fP",
            "find only canonical seed sites of miRNAs.");

    // Parse command line
    ArgumentParser::ParseResult res = parse(parser, argc, argv);
    if (res != ArgumentParser::PARSE_OK)
    {
        return res;
    }

    // Extract arguments
    getArgumentValue(options.mirna_fa, parser, 0);
    getArgumentValue(options.mrna_fa, parser, 1);
    getArgumentValue(options.ofile, parser, 2);
    std::fstream mirna_fa(toCString(options.mirna_fa), std::ios::binary | std::ios::in);
    if (!mirna_fa.good())
    {
        std::cerr << "ERROR: Could not open input file " << toCString(options.mirna_fa) << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }
    std::fstream mrna_fa(toCString(options.mrna_fa), std::ios::binary | std::ios::in);
    if (!mrna_fa.good())
    {
        std::cerr << "ERROR: Could not open input file " << toCString(options.mrna_fa) << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }
    std::fstream ofile(toCString(options.ofile), std::ios::binary | std::ios::out);
    if (!ofile.good())
    {
        std::cerr << "ERROR: Could not open output file " << toCString(options.ofile) << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    // Extract options
    options.canonical = isSet(parser, "canonical");
    options.bruteforce = isSet(parser, "bruteforce");
    options.suffixarray = isSet(parser, "suffixarray");
    options.horspool = isSet(parser, "horspool");
    if (options.bruteforce && options.suffixarray)
    {
        std::cerr << "ERROR: You cannot specify both brute-force and suffix-array!" << std::endl;
        return ArgumentParser::PARSE_ERROR;
    }

    return ArgumentParser::PARSE_OK;
}

int read_fasta(
        CharString const &fasta,
        StringSet<CharString> &ids,
        StringSet<RnaString> &seqs)
{
    CharString id;
    CharString seq;

    SequenceStream seqStream(toCString(fasta));
    if (!isGood(seqStream))
    {
        std::cerr << "ERROR: Could not open the file!" << std::endl;
        return 1;
    }

    while (!atEnd(seqStream))
    {
        if (readRecord(id, seq, seqStream) != 0)
        {
            std::cerr << "ERROR: Could not read from " << toCString(fasta) << "!" << std::endl;
            return 1;
        }

        toUpper(seq);
        for (unsigned i = 0; i < length(seq); ++i)
        {
            if (seq[i] == 'T')
            {
                seq[i] = 'U';
            }
        }
        appendValue(ids, id);
        appendValue(seqs, seq);
    }

    return 0;
}

int get_seed(RnaString const &mirna, RnaString &seed)
{
    for (unsigned i = 0; i < length(seed); ++i)
    {
        seed[i] = mirna[i+1];
    }
    reverseComplement(seed);
    return 0;

}

int search_by_bruteforce(
        std::ofstream &ofile,
        StringSet<CharString> const &mirna_ids,
        StringSet<RnaString> const &mirnas,
        StringSet<CharString> const &mrna_ids,
        StringSet<RnaString> &mrnas)
{
    RnaString seed = "AAAAAA";
    int get_seed_res;

    for (unsigned i = 0; i < length(mirnas); ++i)
    {

        get_seed_res = get_seed(mirnas[i], seed);
        if (get_seed_res != 0)
        {
            std::cerr << "ERROR: Could not get the seed sequence of " << toCString((CharString)mirnas[i]) << std::endl;
            return 1;
        }

        Pattern<RnaString, Simple> pattern(seed);

        for (unsigned j = 0; j < length(mrnas); ++j)
        {
            Finder<RnaString> finder(mrnas[j]);
            while (find(finder, pattern))
            {
                ofile << toCString((CharString)mirna_ids[i]) << "\t" << toCString((CharString)mrna_ids[j]);
                ofile << "\t" << beginPosition(finder) << std::endl;
            }
        }

    }

    return 0;
}

int search_by_horspool(
        std::ofstream &ofile,
        StringSet<CharString> const &mirna_ids,
        StringSet<RnaString> const &mirnas,
        StringSet<CharString> const &mrna_ids,
        StringSet<RnaString> &mrnas)
{
    RnaString seed = "AAAAAA";
    int get_seed_res;

    for (unsigned i = 0; i < length(mirnas); ++i)
    {

        get_seed_res = get_seed(mirnas[i], seed);
        if (get_seed_res != 0)
        {
            std::cerr << "ERROR: Could not get the seed sequence of " << toCString((CharString)mirnas[i]) << std::endl;
            return 1;
        }

        Pattern<RnaString, Horspool> pattern(seed);

        for (unsigned j = 0; j < length(mrnas); ++j)
        {
            Finder<RnaString> finder(mrnas[j]);
            while (find(finder, pattern))
            {
                ofile << toCString((CharString)mirna_ids[i]) << "\t" << toCString((CharString)mrna_ids[j]);
                ofile << "\t" << beginPosition(finder) << std::endl;
            }
        }

    }

    return 0;
}
int search_by_suffixarray(
        std::ofstream &ofile,
        StringSet<CharString> const &mirna_ids,
        StringSet<RnaString> const &mirnas,
        StringSet<CharString> const &mrna_ids,
        StringSet<RnaString> &mrnas)
{

    typedef Index<StringSet<RnaString>, IndexEsa<> > TIndexEsa;
    TIndexEsa index(mrnas);
    Finder<TIndexEsa> finder(index);

    RnaString seed = "AAAAAA";
    int get_seed_res;
    int mrna_id;
    int site_pos;

    for (unsigned i = 0; i < length(mirnas); ++i)
    {

        get_seed_res = get_seed(mirnas[i], seed);
        if (get_seed_res != 0)
        {
            std::cerr << "ERROR: Could not get the seed sequence of " << toCString((CharString)mirnas[i]) << std::endl;
            return 1;
        }

        goBegin(finder);
        clear(finder);
        while (find(finder, seed))
        {
            mrna_id = position(finder).i1;
            site_pos= position(finder).i2;

            ofile << toCString((CharString)mirna_ids[i]) << "\t" << toCString((seqan::CharString)mrna_ids[mrna_id]);
            ofile << "\t" << site_pos << std::endl;
        }

    }

    return 0;
}

int search_by_ngram(
        std::ofstream &ofile,
        StringSet<CharString> const &mirna_ids,
        StringSet<RnaString> const &mirnas,
        StringSet<CharString> const &mrna_ids,
        StringSet<RnaString> &mrnas)
{
    typedef Index<StringSet<RnaString>, IndexQGram<UngappedShape<6> > > TIndexQGram;
    TIndexQGram index(mrnas);
    Finder<TIndexQGram> finder(index);

    RnaString seed = "AAAAAA";
    int get_seed_res;
    int mrna_id;
    int site_pos;

    for (unsigned i = 0; i < length(mirnas); ++i)
    {

        get_seed_res = get_seed(mirnas[i], seed);
        if (get_seed_res != 0)
        {
            std::cerr << "ERROR: Could not get the seed sequence of " << toCString((CharString)mirnas[i]) << std::endl;
            return 1;
        }

        goBegin(finder);
        clear(finder);
        while (find(finder, seed))
        {
            mrna_id = position(finder).i1;
            site_pos= position(finder).i2;

            ofile << toCString((CharString)mirna_ids[i]) << "\t" << toCString((seqan::CharString)mrna_ids[mrna_id]);
            ofile << "\t" << site_pos << std::endl;
        }
    }

    return 0;
}

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
