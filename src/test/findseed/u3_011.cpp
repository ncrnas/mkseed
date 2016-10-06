#include <iostream>
#include "test_ds.hpp"

namespace {

    class U3011 : public TestDS
    {
    protected:
        U3011() {
            IFNAME1 = (char *)"mir_002.fa";
            IFNAME2 = (char *)"utr3_011.fa";
            OFNAME1 = (char *)"test_output_11.txt";
            OFNAME2 = (char *)"test_ds11.txt";
        }
    };

    TEST_F(U3011, mrna_fasta) {
        read_files();

        EXPECT_EQ(1u, length(mrna_ids));

        const char *id1 = "NM_000001*chr1*-*0000001*9000001";
        EXPECT_STREQ(id1, seqan::toCString(mrna_ids[0]));

        EXPECT_EQ(1u, length(mrna_seqs));
        const char *seq1 = "UUUGAGUUCUGGCCACCAACAAUUUAGUCAUAUCUGAUAGGUACAAAAGA"
                "AAACCAAGAUUUUGAUAUGACCACCUUUCAACACUUUUUACUGCAACUAG";
        EXPECT_STREQ(seq1, seqan::toCString((seqan::CharString)mrna_seqs[0]));
    }

    TEST_F(U3011, comp_ngram) {
        comp_ngram();
    }

    TEST_F(U3011, comp_bruteforce) {
        comp_bruteforce();
    }

    TEST_F(U3011, comp_horspool) {
        comp_horspool();
    }

    TEST_F(U3011, comp_suffixarray) {
        comp_suffixarray();
    }
}