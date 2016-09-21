#include<string>
#include "gtest/gtest.h"
#include "mkseed_core.hpp"
#include "get_data_path.hpp"
#include "basic_check.hpp"

class TestDS : public ::testing::Test
{
protected:
    char *IFNAME1;
    char *IFNAME2;
    char *OFNAME1;
    char *OFNAME2;

    virtual void SetUp() {
        dfile = STRINGIZE(TEST_DATA_PATH);

        ifile1 = dfile;
        ifile2 = dfile;
        ifile1 += IFNAME1;
        ifile2 += IFNAME2;

        ofile1 = dfile;
        ofile2 = dfile;
        ofile1 += OFNAME1;
        ofile2 += OFNAME2;

        seqan::clear(mirna_ids);
        seqan::clear(mrna_ids);
        seqan::clear(mirna_seqs);
        seqan::clear(mrna_seqs);
    }

    virtual void TearDown() {
    }

    void read_files() {
        fread_res = read_fasta(ifile1, mirna_ids, mirna_seqs);
        EXPECT_EQ(0, fread_res);

        fread_res = read_fasta(ifile2, mrna_ids, mrna_seqs);
        EXPECT_EQ(0, fread_res);
    }

    int fread_res;
    seqan::CharString dfile;
    seqan::CharString ifile1;
    seqan::CharString ifile2;
    seqan::CharString ofile1;
    seqan::CharString ofile2;
    seqan::StringSet<seqan::CharString> mirna_ids;
    seqan::StringSet<seqan::CharString> mrna_ids;
    seqan::StringSet<seqan::RnaString> mirna_seqs;
    seqan::StringSet<seqan::RnaString> mrna_seqs;
};
