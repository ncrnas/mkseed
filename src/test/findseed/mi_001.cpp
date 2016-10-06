#include <iostream>
#include "test_ds.hpp"

namespace {

    class MI001 : public TestDS
    {
    protected:
        MI001() {
            IFNAME1 = (char *)"mir_001.txt";
            IFNAME2 = (char *)"utr3_001.txt";
            OFNAME1 = (char *)"test_output_1.txt";
            OFNAME2 = (char *)"test_ds1.txt";
        }
    };

    TEST_F(MI001, mirna_id) {
        read_files();

        EXPECT_EQ(2u, length(mirna_ids));

        EXPECT_STREQ("hsa-miR-124 MIMAT0000422", seqan::toCString(mirna_ids[0]));
        EXPECT_STREQ("hsa-miR-1 MIMAT0000416", seqan::toCString(mirna_ids[1]));

    }

    TEST_F(MI001, mirna_seq) {
        read_files();

        EXPECT_EQ(2u, length(mirna_seqs));

        const char *seq1 = "UAAGGCACGCGGUGAAUGCC";
        EXPECT_STREQ(seq1, seqan::toCString((seqan::CharString)mirna_seqs[0]));

        const char *seq2 = "UGGAAUGUAAAGAAGUAUGUAU";
        EXPECT_STREQ(seq2, seqan::toCString((seqan::CharString)mirna_seqs[1]));
    }

    TEST_F(MI001, get_seed) {
        read_files();

        int get_seed_res;
        seqan::RnaString seed = "AAAAAA";

        get_seed_res = get_seed(mirna_seqs[0], seed);
        EXPECT_EQ(0, get_seed_res);
        EXPECT_STREQ("UGCCUU", seqan::toCString((seqan::CharString)seed));

        get_seed_res = get_seed(mirna_seqs[1], seed);
        EXPECT_EQ(0, get_seed_res);
        EXPECT_STREQ("CAUUCC", seqan::toCString((seqan::CharString)seed));
    }
}