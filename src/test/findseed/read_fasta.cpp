#include<string>
#include "gtest/gtest.h"
#include "mkseed_core.hpp"
#include "get_data_path.hpp"
#include "basic_check.hpp"

namespace {

    class fasta_test : public ::testing::Test
    {
    protected:
        virtual void SetUp() {
            dfile = STRINGIZE(TEST_DATA_PATH);
            seqan::clear(mirna_ids);
            seqan::clear(mrna_ids);
            seqan::clear(mirna_seqs);
            seqan::clear(mrna_seqs);
        }

        virtual void TearDown() {
        }

        int fread_res;
        seqan::CharString dfile;
        seqan::StringSet<seqan::CharString> mirna_ids;
        seqan::StringSet<seqan::CharString> mrna_ids;
        seqan::StringSet<seqan::RnaString> mirna_seqs;
        seqan::StringSet<seqan::RnaString> mrna_seqs;
    };

    TEST_F(fasta_test, mirna_fasta) {
        dfile += "mature_single_1.fa";

        fread_res = read_fasta(dfile, mirna_ids, mirna_seqs);
        EXPECT_EQ(0, fread_res);

        EXPECT_EQ(1u, length(mirna_ids));
        EXPECT_STREQ("hsa-miR-0001*MIMAT0000001", seqan::toCString(mirna_ids[0]));

        EXPECT_EQ(1u, length(mirna_seqs));
        const char *seq1 = "CGUGCCACCCUUUUCCCCAG";
        //seqan::CharString str3 = (seqan::CharString)mirna_seqs[0];
        //const char *seq2 = seqan::toCString(str3);
        EXPECT_STREQ(seq1, seqan::toCString((seqan::CharString)mirna_seqs[0]));
    }

    TEST_F(fasta_test, mrna_fasta) {
        dfile += "refseq_single_1.fa";

        fread_res = read_fasta(dfile, mrna_ids, mrna_seqs);
        EXPECT_EQ(0, fread_res);

        EXPECT_EQ(1u, length(mrna_ids));
        EXPECT_STREQ("NM_000001*chr1*-*0000001*9000001", seqan::toCString(mrna_ids[0]));

        EXPECT_EQ(1u, length(mrna_seqs));
        const char *text = "UUUGAGUUCUGGCCACCAACAAUUUAGUCAUAUCUGAUAGGUACAAAAGA"
                "AAACCAAGAUUUUGAUAUGACCACCUUUCAACACUUUUUACUGCAACUAG";
        EXPECT_STREQ(text, seqan::toCString((seqan::CharString)mrna_seqs[0]));
    }

    TEST_F(fasta_test, comp_oput) {
        seqan::CharString dfile1;
        seqan::CharString dfile2;

        dfile1 = dfile;
        dfile1 += "mature_single_1.fa";
        dfile2 = dfile;
        dfile2 += "mature_single_1.fa";
        int fres = gtest_compare_two_files(dfile1, dfile2);
        if (fres != 0)
        {
            std::cerr << "ERROR: Could not open the file!" << std::endl;
        }
    }

}